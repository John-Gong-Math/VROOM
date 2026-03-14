#pragma once
#include "ec.hpp"
#include <vector>
#include <cstdint>
#include <cstddef>
#include <thread>
#include <atomic>

// Pippenger bucket method MSM for G1.
// Uses Booth (signed digit) encoding to halve bucket count, improving cache
// behavior during scatter. Prefetches future buckets to hide memory latency.
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

// Pre-encode all scalar digits using Booth (signed) encoding.
// Outputs: digits[w * npoints + i] = signed digit for point i in window w.
// An extra window is included (num_windows = ceil(scalar_bits/wbits) + 1) to
// absorb any carry from the topmost real window.
inline void booth_encode_scalars(
    const uint8_t *const *scalars,
    size_t npoints,
    size_t scalar_bytes,
    size_t wbits,
    size_t num_windows,
    int32_t *digits
) {
    uint32_t half = 1u << (wbits - 1);

    for (size_t i = 0; i < npoints; i++) {
        uint32_t carry = 0;
        for (size_t w = 0; w < num_windows; w++) {
            size_t bit_pos = w * wbits;
            uint32_t raw = extract_bits(scalars[i], scalar_bytes, bit_pos, wbits) + carry;
            if (raw > half) {
                digits[w * npoints + i] = static_cast<int32_t>(raw) - static_cast<int32_t>(1u << wbits);
                carry = 1;
            } else {
                digits[w * npoints + i] = static_cast<int32_t>(raw);
                carry = 0;
            }
        }
        // Carry should be 0 for valid scalars < 2^scalar_bits with the extra window.
    }
}

// Scatter with Booth-encoded (signed) digits and prefetching.
// Uses XYZZ buckets for cheaper mixed addition (no mul_3b).
// digits points to the start of this window's digit array (size >= end).
template<class Curve, class Ring>
void scatter_chunk(
    const Curve &curve,
    const Ring &ring,
    const typename Curve::AffPoint *points,
    const int32_t *digits,
    size_t start, size_t end,
    typename Curve::XYZZPt *buckets,
    uint8_t *occupied,
    size_t /*nbuckets*/
) {
    using AffPoint = typename Curve::AffPoint;

    static constexpr size_t PREFETCH_AHEAD = 4;

    for (size_t i = start; i < end; i++) {
        // Prefetch a future bucket to hide memory latency
        if (i + PREFETCH_AHEAD < end) {
            int32_t future_digit = digits[i + PREFETCH_AHEAD];
            if (future_digit != 0) {
                size_t future_bidx = static_cast<size_t>(
                    future_digit < 0 ? -future_digit : future_digit) - 1;
                __builtin_prefetch(&buckets[future_bidx], 1, 1);
                __builtin_prefetch(
                    reinterpret_cast<const char*>(&buckets[future_bidx]) + 64, 1, 1);
                __builtin_prefetch(
                    reinterpret_cast<const char*>(&buckets[future_bidx]) + 128, 1, 1);
                __builtin_prefetch(&occupied[future_bidx], 1, 1);
            }
        }

        int32_t digit = digits[i];
        if (digit == 0) continue;

        bool neg = digit < 0;
        size_t bidx = static_cast<size_t>(neg ? -digit : digit) - 1;

        AffPoint pt = points[i];
        if (neg) {
            pt = curve.negate_affine(pt, ring);
        }

        if (!occupied[bidx]) {
            buckets[bidx] = curve.xyzz_from_affine(pt, ring);
            occupied[bidx] = 1;
        } else {
            // Check if bucket is identity (ZZ=0) from prior point cancellation.
            // xyzz_add_affine degenerates when ZZ=0: ZZ3=ZZ*HH=0 traps at identity.
            if (Ring::is_zero(buckets[bidx].ZZ)) {
                buckets[bidx] = curve.xyzz_from_affine(pt, ring);
            } else {
                buckets[bidx] = curve.xyzz_add_affine(buckets[bidx], pt, ring);
            }
        }
    }
}

// Integrate XYZZ buckets using running sum (summation by parts).
// Converts each XYZZ bucket to ProjPoint and uses the complete PointAdd
// (Algorithm 7) for accumulation. This avoids the XYZZ addition formula's
// degeneration when P == Q (which happens regularly during integration when
// running_sum is added to window_sum without being updated between steps).
template<class Curve, class Ring>
typename Curve::ProjPoint integrate_buckets(
    const Curve &curve,
    const Ring &ring,
    typename Curve::XYZZPt *buckets,
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
            // Skip identity buckets (ZZ=0 from point cancellation, e.g. P + (-P)).
            // xyzz_to_proj maps identity XYZZ to (0,0,0) which is invalid for
            // Algorithm 7's complete addition.
            if (!Ring::is_zero(buckets[j].ZZ)) {
                if (!running_started) {
                    running_sum = curve.xyzz_to_proj(buckets[j], ring);
                    running_started = true;
                } else {
                    auto bucket_proj = curve.xyzz_to_proj(buckets[j], ring);
                    running_sum = curve.add_point(running_sum, bucket_proj, ring);
                }
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

// Single-threaded Pippenger MSM with Booth encoding and prefetching.
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
    size_t nbuckets = 1u << (wbits - 1);   // Halved by Booth encoding
    size_t scalar_bytes = (scalar_bits + 7) / 8;

    // +1 extra window to absorb carry from Booth encoding of the top window
    size_t num_windows = (scalar_bits + wbits - 1) / wbits + 1;

    // Pre-encode all scalar digits
    std::vector<int32_t> digits(num_windows * npoints);
    booth_encode_scalars(scalars, npoints, scalar_bytes, wbits,
                         num_windows, digits.data());

    // Reusable XYZZ bucket storage
    using XYZZPt = typename Curve::XYZZPt;
    std::vector<XYZZPt> buckets(nbuckets, curve.xyzz_zero(ring));
    std::vector<uint8_t> occupied(nbuckets, 0);

    ProjPoint result = curve.zero(ring);
    bool result_initialized = false;

    for (size_t w = num_windows; w-- > 0; ) {
        if (result_initialized) {
            for (size_t d = 0; d < wbits; d++) {
                result = curve.double_point(result, ring);
            }
        }

        std::fill(occupied.begin(), occupied.end(), 0);
        scatter_chunk(curve, ring, points,
                      digits.data() + w * npoints,
                      0, npoints,
                      buckets.data(), occupied.data(), nbuckets);

        ProjPoint window_sum = integrate_buckets(curve, ring,
            buckets.data(), occupied.data(), nbuckets);

        // Check if window produced a non-trivial result
        bool has_points = false;
        for (size_t j = 0; j < nbuckets; j++) {
            if (occupied[j]) { has_points = true; break; }
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

// Multi-threaded Pippenger MSM with per-window parallelism.
// Each thread processes complete windows (scatter all points + integrate),
// with zero inter-thread synchronization during compute. An atomic counter
// provides dynamic load balancing — fast-finishing threads grab the next window.
// Final Horner reduction is sequential (~255 doublings).
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
    using XYZZPt = typename Curve::XYZZPt;

    if (num_threads == 0) {
        num_threads = std::thread::hardware_concurrency();
        if (num_threads == 0) num_threads = 1;
    }

    if (num_threads == 1 || npoints < 1024) {
        return pippenger_msm(curve, ring, points, scalars, npoints, scalar_bits);
    }

    if (npoints == 0) return curve.zero(ring);

    size_t wbits = choose_wbits(npoints);
    size_t nbuckets = 1u << (wbits - 1);   // Halved by Booth encoding
    size_t scalar_bytes = (scalar_bits + 7) / 8;
    size_t num_windows = (scalar_bits + wbits - 1) / wbits + 1;

    // Cap threads at number of windows — no point having more
    if (num_threads > num_windows) num_threads = num_windows;

    // Pre-encode all scalar digits (shared read-only across threads)
    std::vector<int32_t> digits(num_windows * npoints);
    booth_encode_scalars(scalars, npoints, scalar_bytes, wbits,
                         num_windows, digits.data());

    // Per-window results (written by exactly one thread each)
    std::vector<ProjPoint> window_results(num_windows, curve.zero(ring));
    std::vector<uint8_t> window_has_result(num_windows, 0);

    // Atomic work counter for dynamic load balancing
    std::atomic<size_t> next_window{0};

    // Worker function: each thread processes complete windows
    auto worker = [&]() {
        // Thread-local bucket storage — reused across windows for cache warmth
        std::vector<XYZZPt> buckets(nbuckets, curve.xyzz_zero(ring));
        std::vector<uint8_t> occupied(nbuckets, 0);

        while (true) {
            size_t w = next_window.fetch_add(1, std::memory_order_relaxed);
            if (w >= num_windows) break;

            // Clear occupied flags for this window
            std::fill(occupied.begin(), occupied.end(), 0);

            // Scatter all points into buckets for this window
            scatter_chunk(curve, ring, points,
                          digits.data() + w * npoints,
                          0, npoints,
                          buckets.data(), occupied.data(), nbuckets);

            // Integrate buckets into window result
            window_results[w] = integrate_buckets(curve, ring,
                buckets.data(), occupied.data(), nbuckets);

            // Check if window produced a non-trivial result
            bool has_points = false;
            for (size_t j = 0; j < nbuckets; j++) {
                if (occupied[j]) { has_points = true; break; }
            }
            window_has_result[w] = has_points ? 1 : 0;
        }
    };

    // Spawn num_threads-1 workers, main thread runs worker() too
    std::vector<std::thread> threads;
    threads.reserve(num_threads - 1);
    for (size_t t = 1; t < num_threads; t++) {
        threads.emplace_back(worker);
    }
    worker();  // Main thread participates

    for (auto &th : threads) {
        th.join();
    }

    // Sequential Horner reduction: accumulate window results from high to low
    ProjPoint result = curve.zero(ring);
    bool result_initialized = false;

    for (size_t w = num_windows; w-- > 0; ) {
        if (result_initialized) {
            for (size_t d = 0; d < wbits; d++) {
                result = curve.double_point(result, ring);
            }
        }

        if (window_has_result[w]) {
            if (!result_initialized) {
                result = window_results[w];
                result_initialized = true;
            } else {
                result = curve.add_point(result, window_results[w], ring);
            }
        }
    }

    return result;
}
