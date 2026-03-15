#pragma once
#include "ec.hpp"
#include "batch_inversion.hpp"
#include <vector>
#include <cstdint>
#include <cstddef>
#include <thread>
#include <atomic>
#include <algorithm>

// ============================================================
// Pippenger MSM with batch affine schedule.
//
// Self-contained re-implementation of the Rust msm_best algorithm
// using VROOM's native types: AffinePoint, ProjectivePoint,
// BoundedRing, and batch_invert.
//
// All buckets are affine (2 field elements = 256 bytes each).
// Conflict points (same bucket hit twice in one batch) go to
// projective accumulators. At integration, both are combined.
//
// Optimizations over v1:
// - field_mul: multiply arbitrary-bounds elements in 1 mul
//   (replaces field_normalize + modmul = 2 muls)
// - Optimized affine addition: 6M per pair (was 9M)
// - PointMixedAdd in integration: 13M per bucket (was 15M)
// - Point-parallel threading with Horner evaluation
// ============================================================

// ---- Scalar encoding utilities ----

// Window size heuristic (BLST-style).
inline size_t msm_wbits(size_t npoints) {
    size_t wbits = 0;
    for (size_t n = npoints; n >>= 1; wbits++);
    if (wbits > 12) return wbits - 3;
    if (wbits > 8)  return wbits - 2;
    if (wbits > 4)  return wbits - 1;
    return wbits ? 2 : 1;
}

// Extract wbits bits from a little-endian scalar byte array at start_bit.
inline uint32_t msm_extract_bits(const uint8_t *scalar, size_t scalar_bytes,
                                  size_t start_bit, size_t wbits) {
    uint32_t val = 0;
    for (size_t b = 0; b < wbits; b++) {
        size_t abs_bit = start_bit + b;
        size_t byte_idx = abs_bit / 8;
        if (byte_idx >= scalar_bytes) break;
        if ((scalar[byte_idx] >> (abs_bit % 8)) & 1)
            val |= (1u << b);
    }
    return val;
}

// Booth (signed-digit) encoding of all scalars.
// Output: digits[w * npoints + i] = signed digit for point i in window w.
// An extra window absorbs any carry from the topmost real window.
inline void msm_booth_encode(
    const uint8_t *const *scalars,
    size_t npoints,
    size_t scalar_bytes,
    size_t wbits,
    size_t num_windows,
    int32_t *digits)
{
    uint32_t half = 1u << (wbits - 1);
    for (size_t i = 0; i < npoints; i++) {
        uint32_t carry = 0;
        for (size_t w = 0; w < num_windows; w++) {
            uint32_t raw = msm_extract_bits(scalars[i], scalar_bytes,
                                             w * wbits, wbits) + carry;
            if (raw > half) {
                digits[w * npoints + i] = static_cast<int32_t>(raw)
                    - static_cast<int32_t>(1u << wbits);
                carry = 1;
            } else {
                digits[w * npoints + i] = static_cast<int32_t>(raw);
                carry = 0;
            }
        }
    }
}

// ---- Field helpers ----

// Normalize any ExpandedElement to StandardElement via multiply-by-one.
template<class Ring, class ElemType, class RNSType>
typename Ring::StandardElement field_normalize(
    const Ring &ring,
    const ExpandedElement<ElemType, RNSType> &element)
{
    auto p = ring.prep_left(element);
    auto one_p = ring.prep(ring.one());
    auto product = p * one_p;
    auto r = ring.template ready<Ring::MAX_ADD>(product);
    using RE = typename Ring::template ReadyElement<Ring::MAX_ADD>;
    std::array<RE, 1> arr = {r};
    return ring.template batch_reduce_expand<1>(arr)[0];
}

// Multiply two arbitrary-bounds elements, returning StandardElement.
// prep_left/prep handle negative and wide bounds via reduce_auto.
// Replaces the pattern field_normalize(x) + modmul(result, y) (2 muls)
// with a single multiplication path (1 mul).
template<class Ring, class Elem1, class RNS1, class Elem2, class RNS2>
typename Ring::StandardElement field_mul(
    const Ring &ring,
    const ExpandedElement<Elem1, RNS1> &a,
    const ExpandedElement<Elem2, RNS2> &b)
{
    auto product = ring.prep_left(a) * ring.prep(b);
    auto r = ring.template ready<Ring::MAX_ADD>(product);
    using RE = typename Ring::template ReadyElement<Ring::MAX_ADD>;
    std::array<RE, 1> arr = {r};
    return ring.template batch_reduce_expand<1>(arr)[0];
}

// Field subtraction returning StandardElement.
template<class Ring>
typename Ring::StandardElement field_sub(
    const Ring &ring,
    const typename Ring::StandardElement &a,
    const typename Ring::StandardElement &b)
{
    return field_normalize(ring, a + ring.standard_negate(b));
}

// ---- Affine point helpers ----

template<class Ring>
AffinePoint<typename Ring::StandardElement> negate_affine(
    const AffinePoint<typename Ring::StandardElement> &pt,
    const Ring &ring)
{
    return AffinePoint<typename Ring::StandardElement>(
        pt.x, ring.standard_negate(pt.y));
}

// ---- Schedule: batch affine additions with projective conflict fallback ----

template<class Ring>
class AffineSchedule {
    using StdElem = typename Ring::StandardElement;
    using AffPt   = AffinePoint<StdElem>;
    using ProjPt  = ProjectivePoint<StdElem>;

    const Ring &ring;

    size_t capacity;
    size_t count;
    size_t nbuckets_;
    bool window_has_points_;

    // Per-slot metadata
    std::vector<size_t>  slot_bucket;
    std::vector<size_t>  slot_base;
    std::vector<uint8_t> slot_neg;

    // Membership: is this bucket in the current schedule batch?
    std::vector<uint8_t> in_sched;

    // Flush temporaries (reused across flushes)
    std::vector<size_t>  add_slots;
    std::vector<AffPt>   add_points;
    std::vector<StdElem> dx_vec;
    std::vector<uint8_t> fallback_flags;

    // Buckets touched in the previous window; used for O(touched) reset.
    std::vector<size_t> touched_affine;
    std::vector<size_t> touched_proj;

public:
    // Affine buckets: primary accumulation target
    std::vector<AffPt>   buckets;
    std::vector<uint8_t> occupied;

    // Projective conflict accumulators: rare fallback
    std::vector<ProjPt>  proj_buckets;
    std::vector<uint8_t> proj_occupied;

    AffineSchedule(const Ring &r, size_t nb, size_t cap)
        : ring(r), capacity(cap), count(0), nbuckets_(nb)
        , window_has_points_(false)
        , slot_bucket(cap), slot_base(cap), slot_neg(cap)
        , in_sched(nb, 0)
        , buckets(nb), occupied(nb, 0)
        , proj_buckets(nb, ProjPt(r.zero(), r.one(), r.zero()))
        , proj_occupied(nb, 0)
    {
        add_slots.reserve(capacity);
        add_points.reserve(capacity);
        dx_vec.reserve(capacity);
        fallback_flags.reserve(capacity);
        touched_affine.reserve(nbuckets_);
        touched_proj.reserve(nbuckets_);
    }

    bool contains(size_t bidx) const { return in_sched[bidx]; }

    // Add a point to the schedule. Auto-flushes when full.
    void add(const AffPt *bases, size_t base_idx, size_t bidx, bool neg) {
        if (count >= capacity) flush(bases);
        slot_bucket[count] = bidx;
        slot_base[count]   = base_idx;
        slot_neg[count]    = neg ? 1 : 0;
        in_sched[bidx]     = 1;
        window_has_points_ = true;
        count++;
    }

    // Greedy conflict: add point directly to projective accumulator.
    void conflict_add(const AffPt *bases, size_t base_idx, size_t bidx, bool neg) {
        AffPt pt = bases[base_idx];
        if (neg) pt = negate_affine(pt, ring);

        if (!proj_occupied[bidx]) {
            proj_buckets[bidx]  = ProjPt(pt.x, pt.y, ring.one());
            proj_occupied[bidx] = 1;
            touched_proj.push_back(bidx);
        } else {
            proj_buckets[bidx] = PointMixedAdd(proj_buckets[bidx], pt, ring);
        }
        window_has_points_ = true;
    }

    // Final flush at end of window.
    void execute(const AffPt *bases) {
        if (count > 0) flush(bases);
    }

    // Reset for a new window.
    void reset() {
        count = 0;
        window_has_points_ = false;

        for (size_t bidx : touched_affine)
            occupied[bidx] = 0;
        touched_affine.clear();

        for (size_t bidx : touched_proj)
            proj_occupied[bidx] = 0;
        touched_proj.clear();
    }

    bool has_points() const {
        if (!window_has_points_) return false;

        for (size_t bidx : touched_affine) {
            if (occupied[bidx]) return true;
        }
        for (size_t bidx : touched_proj) {
            if (proj_occupied[bidx]) return true;
        }
        return false;
    }

private:
    AffPt resolve(const AffPt *bases, size_t slot) const {
        AffPt pt = bases[slot_base[slot]];
        if (slot_neg[slot]) pt = negate_affine(pt, ring);
        return pt;
    }

    void flush(const AffPt *bases) {
        add_slots.clear();

        // Phase 1: empty affine buckets get first point directly.
        // Non-empty buckets are queued for batch affine addition.
        for (size_t s = 0; s < count; s++) {
            size_t bidx = slot_bucket[s];
            if (!occupied[bidx]) {
                buckets[bidx]  = resolve(bases, s);
                occupied[bidx] = 1;
                touched_affine.push_back(bidx);
            } else {
                add_slots.push_back(s);
            }
        }

        // Phase 2: batch affine-affine addition via batch inversion.
        if (!add_slots.empty()) {
            size_t n = add_slots.size();
            if (dx_vec.size() < n) dx_vec.resize(n);
            if (fallback_flags.size() < n) fallback_flags.resize(n);
            add_points.clear();

            // Resolve conflicted points once; reused across both loops below.
            for (size_t i = 0; i < n; i++) {
                add_points.push_back(resolve(bases, add_slots[i]));
                fallback_flags[i] = 0;
            }

            // 2a. Compute dx = Q.x - P.x for each pair.
            // Normalized to StandardElement for batch_invert. (1 mul each)
            for (size_t i = 0; i < n; i++) {
                size_t bidx = slot_bucket[add_slots[i]];
                const AffPt &pt = add_points[i];
                dx_vec[i] = field_normalize(ring, pt.x - buckets[bidx].x);
                if (Ring::is_zero(dx_vec[i])) {
                    fallback_flags[i] = 1;
                    dx_vec[i] = ring.one();  // placeholder
                }
            }

            // 2b. Single batch inversion.
            auto inv_dx = batch_invert(ring, dx_vec);

            // 2c. Optimized affine addition: 6M per pair (was 9M).
            //
            // Cost breakdown:
            //   dx normalize:          1M (step 2a above)
            //   field_mul(dy, inv_dx): 1M (lambda)
            //   modmul(lambda, lambda):1M (lam_sq)
            //   field_mul(diff, lam):  1M (y3_mul)
            //   batch normalize x3,y3: 2M (batched into 1 call)
            //   Total:                 6M
            for (size_t i = 0; i < n; i++) {
                size_t bidx = slot_bucket[add_slots[i]];
                const AffPt &pt = add_points[i];

                if (fallback_flags[i]) {
                    // dx = 0: either doubling (P == Q) or cancellation (P == -Q).
                    auto dy = field_sub(ring, pt.y, buckets[bidx].y);
                    if (Ring::is_zero(dy)) {
                        // P + (-P) = O. Mark bucket empty.
                        if (occupied[bidx]) {
                            occupied[bidx] = 0;
                        }
                    } else {
                        // P == Q -> doubling. Use projective formula, convert back.
                        ProjPt proj(buckets[bidx].x, buckets[bidx].y, ring.one());
                        ProjPt dbl = PointDouble(proj, ring);
                        auto z_inv = batch_invert(ring, std::vector<StdElem>{dbl.z});
                        buckets[bidx] = AffPt(
                            ring.modmul(dbl.x, z_inv[0]),
                            ring.modmul(dbl.y, z_inv[0]));
                    }
                    continue;
                }

                auto &P = buckets[bidx];

                // dy_wide: direct subtraction, no normalize (FREE)
                // prep_left/prep handle the wider bounds via reduce_auto.
                auto dy_wide = pt.y - P.y;

                // lambda = dy_wide * inv_dx: 1 mul
                // (replaces field_normalize(dy) + modmul = 2 muls)
                auto lambda = field_mul(ring, dy_wide, inv_dx[i]);

                // lam_sq = lambda^2: 1 modmul
                auto lam_sq = ring.modmul(lambda, lambda);

                // x3_wide = lam_sq - P.x - pt.x (FREE, wider bounds)
                auto x3_wide = lam_sq - P.x - pt.x;

                // diff = P.x - x3_wide (FREE, prep handles wider bounds)
                auto diff = P.x - x3_wide;

                // y3_mul = diff * lambda: 1 mul
                // (replaces field_normalize(P.x - x3) + modmul = 2 muls)
                auto y3_mul = field_mul(ring, diff, lambda);

                // y3_wide = y3_mul - P.y (FREE)
                auto y3_wide = y3_mul - P.y;

                // Batch normalize x3 and y3 together (2 normalizes in 1 batch)
                auto one = ring.one();
                auto x3_product = ring.prep_left(x3_wide) * ring.prep(one);
                auto y3_product = ring.prep_left(y3_wide) * ring.prep(one);
                auto [x3, y3] = ring.batch_reduce_expand(x3_product, y3_product);

                P = AffPt(x3, y3);
            }
        }

        // Phase 3: clear schedule membership.
        for (size_t s = 0; s < count; s++)
            in_sched[slot_bucket[s]] = 0;
        count = 0;
    }
};

// ---- Integration: summation-by-parts using PointMixedAdd for affine buckets ----
// Uses PointMixedAdd (13M) instead of converting to ProjPt + PointAdd (15M).

template<class Ring>
ProjectivePoint<typename Ring::StandardElement> integrate_buckets_v2(
    const Ring &ring,
    const AffinePoint<typename Ring::StandardElement>  *aff_buckets,
    const uint8_t                                      *aff_occupied,
    const ProjectivePoint<typename Ring::StandardElement> *proj_buckets,
    const uint8_t                                        *proj_occupied,
    size_t nbuckets)
{
    using StdElem  = typename Ring::StandardElement;
    using ProjPt   = ProjectivePoint<StdElem>;

    ProjPt zero_pt(ring.zero(), ring.one(), ring.zero());
    ProjPt running_sum = zero_pt;
    ProjPt window_sum  = zero_pt;
    bool running_started = false;
    bool window_started  = false;

    for (size_t j = nbuckets; j-- > 0; ) {
        bool has_a = aff_occupied[j];
        bool has_p = proj_occupied[j];

        if (has_a || has_p) {
            if (!running_started) {
                // First bucket: initialize running_sum
                if (has_a) {
                    running_sum = ProjPt(aff_buckets[j].x, aff_buckets[j].y, ring.one());
                    if (has_p) {
                        running_sum = PointAdd(running_sum, proj_buckets[j], ring);
                    }
                } else {
                    running_sum = proj_buckets[j];
                }
                running_started = true;
            } else {
                // PointMixedAdd for affine buckets: 13M vs PointAdd's 15M
                if (has_a) {
                    running_sum = PointMixedAdd(running_sum, aff_buckets[j], ring);
                }
                if (has_p) {
                    running_sum = PointAdd(running_sum, proj_buckets[j], ring);
                }
            }
        }

        if (running_started) {
            if (!window_started) {
                window_sum    = running_sum;
                window_started = true;
            } else {
                window_sum = PointAdd(window_sum, running_sum, ring);
            }
        }
    }

    return window_started ? window_sum : zero_pt;
}

// ============================================================
// Public API
// ============================================================

// Single-threaded Pippenger MSM.
// Horner evaluation (high-to-low windows) minimizes doublings.
template<class Ring>
ProjectivePoint<typename Ring::StandardElement> msm_v2(
    const Ring &ring,
    const AffinePoint<typename Ring::StandardElement> *points,
    const uint8_t *const *scalars,
    size_t npoints,
    size_t scalar_bits = 255)
{
    using StdElem = typename Ring::StandardElement;
    using ProjPt  = ProjectivePoint<StdElem>;

    ProjPt identity(ring.zero(), ring.one(), ring.zero());
    if (npoints == 0) return identity;

    size_t wbits        = msm_wbits(npoints);
    size_t nbuckets     = 1u << (wbits - 1);
    size_t scalar_bytes = (scalar_bits + 7) / 8;
    size_t num_windows  = (scalar_bits + wbits - 1) / wbits + 1;

    std::vector<int32_t> digits(num_windows * npoints);
    msm_booth_encode(scalars, npoints, scalar_bytes, wbits,
                      num_windows, digits.data());

    AffineSchedule<Ring> sched(ring, nbuckets, nbuckets);

    ProjPt result = identity;
    bool result_init = false;

    for (size_t w = num_windows; w-- > 0; ) {
        // Horner shift
        if (result_init) {
            for (size_t d = 0; d < wbits; d++)
                result = PointDouble(result, ring);
        }

        sched.reset();
        const int32_t *wd = digits.data() + w * npoints;

        // Scatter with schedule + prefetch
        for (size_t i = 0; i < npoints; i++) {
            static constexpr size_t PREFETCH_AHEAD = 4;
            if (i + PREFETCH_AHEAD < npoints) {
                int32_t fd = wd[i + PREFETCH_AHEAD];
                if (fd != 0) {
                    size_t fb = static_cast<size_t>(fd < 0 ? -fd : fd) - 1;
                    __builtin_prefetch(&sched.buckets[fb], 1, 1);
                    __builtin_prefetch(&sched.occupied[fb], 0, 1);
                }
            }

            int32_t digit = wd[i];
            if (digit == 0) continue;

            bool neg    = digit < 0;
            size_t bidx = static_cast<size_t>(neg ? -digit : digit) - 1;

            if (sched.contains(bidx)) {
                sched.conflict_add(points, i, bidx, neg);
            } else {
                sched.add(points, i, bidx, neg);
            }
        }

        sched.execute(points);

        // Integrate
        ProjPt window_sum = integrate_buckets_v2(
            ring,
            sched.buckets.data(), sched.occupied.data(),
            sched.proj_buckets.data(), sched.proj_occupied.data(),
            nbuckets);

        if (sched.has_points()) {
            if (!result_init) {
                result      = window_sum;
                result_init = true;
            } else {
                result = PointAdd(result, window_sum, ring);
            }
        }
    }

    return result;
}

// Multi-threaded: per-window work-stealing with batch affine arithmetic.
// Each thread owns a thread-local AffineSchedule, reused across windows
// for cache warmth. A single atomic counter provides zero-synchronization
// dynamic load balancing. Sequential Horner reduction after all windows
// are processed.
template<class Ring>
ProjectivePoint<typename Ring::StandardElement> msm_v2_parallel(
    const Ring &ring,
    const AffinePoint<typename Ring::StandardElement> *points,
    const uint8_t *const *scalars,
    size_t npoints,
    size_t scalar_bits = 255,
    size_t num_threads = 0)
{
    using StdElem = typename Ring::StandardElement;
    using ProjPt  = ProjectivePoint<StdElem>;

    ProjPt identity(ring.zero(), ring.one(), ring.zero());

    if (num_threads == 0) {
        num_threads = std::thread::hardware_concurrency();
        if (num_threads == 0) num_threads = 1;
    }

    if (num_threads == 1 || npoints < 1024)
        return msm_v2(ring, points, scalars, npoints, scalar_bits);

    if (npoints == 0) return identity;

    size_t wbits        = msm_wbits(npoints);
    size_t nbuckets     = 1u << (wbits - 1);
    size_t scalar_bytes = (scalar_bits + 7) / 8;
    size_t num_windows  = (scalar_bits + wbits - 1) / wbits + 1;

    // Cap threads at number of windows — no point having more
    if (num_threads > num_windows) num_threads = num_windows;

    std::vector<int32_t> digits(num_windows * npoints);
    msm_booth_encode(scalars, npoints, scalar_bytes, wbits,
                      num_windows, digits.data());

    // Per-window results (written by exactly one thread each)
    std::vector<ProjPt>  window_results(num_windows, identity);
    std::vector<uint8_t> window_has_result(num_windows, 0);

    // Atomic work counter for dynamic load balancing
    std::atomic<size_t> next_window{0};

    // Worker function: each thread processes complete windows
    auto worker = [&]() {
        // Thread-local AffineSchedule — reused across windows for cache warmth
        AffineSchedule<Ring> sched(ring, nbuckets, nbuckets);

        while (true) {
            size_t w = next_window.fetch_add(1, std::memory_order_relaxed);
            if (w >= num_windows) break;

            sched.reset();
            const int32_t *wd = digits.data() + w * npoints;

            // Scatter ALL npoints for this window (with prefetch)
            for (size_t i = 0; i < npoints; i++) {
                static constexpr size_t PREFETCH_AHEAD = 4;
                if (i + PREFETCH_AHEAD < npoints) {
                    int32_t fd = wd[i + PREFETCH_AHEAD];
                    if (fd != 0) {
                        size_t fb = static_cast<size_t>(fd < 0 ? -fd : fd) - 1;
                        __builtin_prefetch(&sched.buckets[fb], 1, 1);
                        __builtin_prefetch(&sched.occupied[fb], 0, 1);
                    }
                }

                int32_t digit = wd[i];
                if (digit == 0) continue;

                bool neg    = digit < 0;
                size_t bidx = static_cast<size_t>(neg ? -digit : digit) - 1;

                if (sched.contains(bidx)) {
                    sched.conflict_add(points, i, bidx, neg);
                } else {
                    sched.add(points, i, bidx, neg);
                }
            }

            sched.execute(points);

            // Integrate buckets into window result
            window_results[w] = integrate_buckets_v2(
                ring,
                sched.buckets.data(), sched.occupied.data(),
                sched.proj_buckets.data(), sched.proj_occupied.data(),
                nbuckets);

            window_has_result[w] = sched.has_points() ? 1 : 0;
        }
    };

    // Spawn num_threads-1 workers, main thread runs worker() too
    std::vector<std::thread> threads;
    threads.reserve(num_threads - 1);
    for (size_t t = 1; t < num_threads; t++)
        threads.emplace_back(worker);
    worker();  // Main thread participates

    for (auto &th : threads)
        th.join();

    // Sequential Horner reduction: accumulate window results from high to low
    ProjPt result = identity;
    bool result_init = false;

    for (size_t w = num_windows; w-- > 0; ) {
        if (result_init) {
            for (size_t d = 0; d < wbits; d++)
                result = PointDouble(result, ring);
        }

        if (window_has_result[w]) {
            if (!result_init) {
                result      = window_results[w];
                result_init = true;
            } else {
                result = PointAdd(result, window_results[w], ring);
            }
        }
    }

    return result;
}
