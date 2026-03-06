#pragma once
#include "ec.hpp"
#include <vector>
#include <cstdint>
#include <cstddef>
#include <thread>
#include <atomic>

// Pippenger bucket method MSM for G1.
// Uses unsigned digit decomposition (no Booth encoding) for V1 simplicity.
//
// NOT constant-time for bucket scatter (same as BLST).

// Window size selection matching BLST's heuristic.
inline size_t choose_wbits(size_t npoints) {
    size_t wbits = 0;
    for (size_t n = npoints; n >>= 1; wbits++);
    if (wbits > 12) return wbits - 3;
    if (wbits > 8)  return wbits - 2;
    if (wbits > 4)  return wbits - 1;
    return wbits ? 2 : 1;
}

// Extract wbits bits from a scalar byte array (little-endian) at start_bit.
// Zero-pads beyond scalar boundary.
inline uint32_t extract_bits(const uint8_t *scalar, size_t scalar_bytes,
                              size_t start_bit, size_t wbits) {
    uint32_t val = 0;
    for (size_t b = 0; b < wbits; b++) {
        size_t abs_bit = start_bit + b;
        size_t byte_idx = abs_bit / 8;
        if (byte_idx >= scalar_bytes) break;
        if ((scalar[byte_idx] >> (abs_bit % 8)) & 1) {
            val |= (1u << b);
        }
    }
    return val;
}

// Process a chunk of points for one window, scattering into the provided buckets.
template<class Curve, class Ring>
void scatter_chunk(
    const Curve &curve,
    const Ring &ring,
    const typename Curve::AffPoint *points,
    const uint8_t *const *scalars,
    size_t start, size_t end,
    size_t scalar_bytes, size_t bit_pos, size_t effective_wbits,
    typename Curve::ProjPoint *buckets,
    uint8_t *occupied,
    size_t /*nbuckets*/
) {
    using ProjPoint = typename Curve::ProjPoint;
    using AffPoint = typename Curve::AffPoint;

    for (size_t i = start; i < end; i++) {
        uint32_t digit = extract_bits(scalars[i], scalar_bytes,
                                       bit_pos, effective_wbits);
        if (digit == 0) continue;

        size_t bidx = digit - 1;
        AffPoint pt = points[i];
        if (!occupied[bidx]) {
            buckets[bidx] = ProjPoint(pt.x, pt.y, ring.one());
            occupied[bidx] = 1;
        } else {
            buckets[bidx] = curve.add_mixed_point(buckets[bidx], pt, ring);
        }
    }
}

// Integrate buckets using running sum (summation by parts).
template<class Curve, class Ring>
typename Curve::ProjPoint integrate_buckets(
    const Curve &curve,
    const Ring &ring,
    typename Curve::ProjPoint *buckets,
    const uint8_t *occupied,
    size_t nbuckets
) {
    using ProjPoint = typename Curve::ProjPoint;

    ProjPoint running_sum = curve.zero(ring);
    ProjPoint window_sum = curve.zero(ring);
    bool running_started = false;
    bool window_started = false;

    for (size_t j = nbuckets; j-- > 0; ) {
        if (occupied[j]) {
            if (!running_started) {
                running_sum = buckets[j];
                running_started = true;
            } else {
                running_sum = curve.add_point(running_sum, buckets[j], ring);
            }
        }
        if (running_started) {
            if (!window_started) {
                window_sum = running_sum;
                window_started = true;
            } else {
                window_sum = curve.add_point(window_sum, running_sum, ring);
            }
        }
    }

    return window_started ? window_sum : curve.zero(ring);
}

// Single-threaded Pippenger MSM.
template<class Curve, class Ring>
typename Curve::ProjPoint pippenger_msm(
    const Curve &curve,
    const Ring &ring,
    const typename Curve::AffPoint *points,
    const uint8_t *const *scalars,
    size_t npoints,
    size_t scalar_bits
) {
    using ProjPoint = typename Curve::ProjPoint;

    if (npoints == 0) return curve.zero(ring);

    size_t wbits = choose_wbits(npoints);
    size_t nbuckets = (1u << wbits) - 1;
    size_t scalar_bytes = (scalar_bits + 7) / 8;

    // Reusable bucket storage
    std::vector<ProjPoint> buckets(nbuckets, curve.zero(ring));
    std::vector<uint8_t> occupied(nbuckets, 0);

    ProjPoint result = curve.zero(ring);
    bool result_initialized = false;
    size_t num_windows = (scalar_bits + wbits - 1) / wbits;

    for (size_t w = num_windows; w-- > 0; ) {
        size_t bit_pos = w * wbits;

        if (result_initialized) {
            for (size_t d = 0; d < wbits; d++) {
                result = curve.double_point(result, ring);
            }
        }

        size_t effective_wbits = wbits;
        if (bit_pos + wbits > scalar_bits) {
            effective_wbits = scalar_bits - bit_pos;
        }

        std::fill(occupied.begin(), occupied.end(), 0);
        scatter_chunk(curve, ring, points, scalars, 0, npoints,
                      scalar_bytes, bit_pos, effective_wbits,
                      buckets.data(), occupied.data(), nbuckets);

        ProjPoint window_sum = integrate_buckets(curve, ring,
            buckets.data(), occupied.data(), nbuckets);

        // Check if window produced a non-trivial result
        bool has_points = false;
        for (size_t j = 0; j < nbuckets && !has_points; j++) {
            if (occupied[j]) has_points = true;
        }

        if (has_points) {
            if (!result_initialized) {
                result = window_sum;
                result_initialized = true;
            } else {
                result = curve.add_point(result, window_sum, ring);
            }
        }
    }

    return result;
}

// Multi-threaded Pippenger MSM.
// Splits points across threads within each window. Each thread scatters into
// its own bucket array. After scatter, partial buckets are merged and integrated.
// Only beneficial for large npoints where scatter dominates.
template<class Curve, class Ring>
typename Curve::ProjPoint pippenger_msm_parallel(
    const Curve &curve,
    const Ring &ring,
    const typename Curve::AffPoint *points,
    const uint8_t *const *scalars,
    size_t npoints,
    size_t scalar_bits,
    size_t num_threads = 0
) {
    using ProjPoint = typename Curve::ProjPoint;

    if (num_threads == 0) {
        num_threads = std::thread::hardware_concurrency();
        if (num_threads == 0) num_threads = 1;
    }

    if (num_threads == 1 || npoints < 1024) {
        return pippenger_msm(curve, ring, points, scalars, npoints, scalar_bits);
    }

    if (npoints == 0) return curve.zero(ring);

    size_t wbits = choose_wbits(npoints);
    size_t nbuckets = (1u << wbits) - 1;
    size_t scalar_bytes = (scalar_bits + 7) / 8;
    size_t num_windows = (scalar_bits + wbits - 1) / wbits;

    if (num_threads > npoints / 64) num_threads = std::max(npoints / 64, (size_t)1);

    // Pre-allocate per-thread bucket arrays (reused across windows)
    std::vector<std::vector<ProjPoint>> all_buckets(num_threads);
    std::vector<std::vector<uint8_t>> all_occupied(num_threads);
    for (size_t t = 0; t < num_threads; t++) {
        all_buckets[t].resize(nbuckets, curve.zero(ring));
        all_occupied[t].resize(nbuckets, 0);
    }

    // Compute chunk boundaries
    size_t chunk = npoints / num_threads;
    size_t remainder = npoints % num_threads;
    std::vector<size_t> starts(num_threads), ends(num_threads);
    for (size_t t = 0; t < num_threads; t++) {
        starts[t] = t * chunk + std::min(t, remainder);
        ends[t] = starts[t] + chunk + (t < remainder ? 1 : 0);
    }

    ProjPoint result = curve.zero(ring);
    bool result_initialized = false;

    for (size_t w = num_windows; w-- > 0; ) {
        size_t bit_pos = w * wbits;

        if (result_initialized) {
            for (size_t d = 0; d < wbits; d++) {
                result = curve.double_point(result, ring);
            }
        }

        size_t effective_wbits = wbits;
        if (bit_pos + wbits > scalar_bits) {
            effective_wbits = scalar_bits - bit_pos;
        }

        // Launch worker threads for scatter (threads 1..num_threads-1)
        std::vector<std::thread> threads;
        threads.reserve(num_threads - 1);

        for (size_t t = 1; t < num_threads; t++) {
            // Clear this thread's occupied flags
            std::fill(all_occupied[t].begin(), all_occupied[t].end(), 0);

            threads.emplace_back([&, t, bit_pos, effective_wbits]() {
                scatter_chunk(curve, ring, points, scalars,
                              starts[t], ends[t],
                              scalar_bytes, bit_pos, effective_wbits,
                              all_buckets[t].data(), all_occupied[t].data(),
                              nbuckets);
            });
        }

        // Main thread does chunk 0
        std::fill(all_occupied[0].begin(), all_occupied[0].end(), 0);
        scatter_chunk(curve, ring, points, scalars,
                      starts[0], ends[0],
                      scalar_bytes, bit_pos, effective_wbits,
                      all_buckets[0].data(), all_occupied[0].data(),
                      nbuckets);

        // Wait for all workers
        for (auto &th : threads) {
            th.join();
        }

        // Merge thread buckets into thread 0
        auto &merged = all_buckets[0];
        auto &merged_occ = all_occupied[0];

        for (size_t t = 1; t < num_threads; t++) {
            for (size_t b = 0; b < nbuckets; b++) {
                if (!all_occupied[t][b]) continue;
                if (!merged_occ[b]) {
                    merged[b] = all_buckets[t][b];
                    merged_occ[b] = 1;
                } else {
                    merged[b] = curve.add_point(merged[b], all_buckets[t][b], ring);
                }
            }
        }

        // Integrate merged buckets
        bool has_points = false;
        for (size_t j = 0; j < nbuckets && !has_points; j++) {
            if (merged_occ[j]) has_points = true;
        }

        if (has_points) {
            ProjPoint window_sum = integrate_buckets(curve, ring,
                merged.data(), merged_occ.data(), nbuckets);

            if (!result_initialized) {
                result = window_sum;
                result_initialized = true;
            } else {
                result = curve.add_point(result, window_sum, ring);
            }
        }
    }

    return result;
}
