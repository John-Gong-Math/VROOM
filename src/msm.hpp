#pragma once
#include "pippenger_v2.hpp"
#include "scalar_mult.hpp"
#include <vector>

// Public MSM API with dispatch.
// Uses batch affine Pippenger (v2) for large point counts,
// windowed scalar mult for small counts.

// Single-threaded MSM.
template<class Curve, class Ring>
typename Curve::ProjPoint msm(
    const Curve &curve,
    const Ring &ring,
    const typename Curve::AffPoint *points,
    const uint8_t *const *scalars,
    size_t npoints,
    size_t scalar_bits = 255
) {
    using ProjPoint = typename Curve::ProjPoint;

    if (npoints == 0) {
        return curve.zero(ring);
    }

    if (npoints == 1) {
        // Single point: use existing windowed scalar mult
        BigInt scalar_bigint(0);
        for (int i = static_cast<int>((scalar_bits + 7) / 8) - 1; i >= 0; i--) {
            scalar_bigint = scalar_bigint * BigInt(256) + BigInt(static_cast<unsigned long>(scalars[0][i]));
        }
        ProjPoint P(points[0].x, points[0].y, ring.one());
        return scalar_mult_table<Curve, Ring, 255, 4>(scalar_bigint, P, curve, ring);
    }

    if (npoints < 32) {
        // Small count: naive accumulation of individual scalar mults
        ProjPoint result = curve.zero(ring);
        for (size_t i = 0; i < npoints; i++) {
            BigInt scalar_bigint(0);
            for (int j = static_cast<int>((scalar_bits + 7) / 8) - 1; j >= 0; j--) {
                scalar_bigint = scalar_bigint * BigInt(256) + BigInt(static_cast<unsigned long>(scalars[i][j]));
            }
            ProjPoint P(points[i].x, points[i].y, ring.one());
            ProjPoint sP = scalar_mult_table<Curve, Ring, 255, 4>(scalar_bigint, P, curve, ring);
            result = curve.add_point(result, sP, ring);
        }
        return result;
    }

    // Large count: batch affine Pippenger (v2)
    return msm_v2(ring, points, scalars, npoints, scalar_bits);
}

// Multi-threaded MSM.
// num_threads=0 means auto-detect via hardware_concurrency().
template<class Curve, class Ring>
typename Curve::ProjPoint msm_parallel(
    const Curve &curve,
    const Ring &ring,
    const typename Curve::AffPoint *points,
    const uint8_t *const *scalars,
    size_t npoints,
    size_t scalar_bits = 255,
    size_t num_threads = 0
) {
    if (npoints < 32) {
        return msm(curve, ring, points, scalars, npoints, scalar_bits);
    }
    return msm_v2_parallel(ring, points, scalars, npoints,
                            scalar_bits, num_threads);
}

// Point-parallel multi-threaded MSM.
// Partitions points across threads for better cache locality.
template<class Curve, class Ring>
typename Curve::ProjPoint msm_point_parallel(
    const Curve &curve,
    const Ring &ring,
    const typename Curve::AffPoint *points,
    const uint8_t *const *scalars,
    size_t npoints,
    size_t scalar_bits = 255,
    size_t num_threads = 0
) {
    if (npoints < 32) {
        return msm(curve, ring, points, scalars, npoints, scalar_bits);
    }
    return msm_v2_point_parallel(ring, points, scalars, npoints,
                                  scalar_bits, num_threads);
}

// Auto-dispatching parallel MSM.
// Selects point-parallel for large N, per-window parallel for small N.
template<class Curve, class Ring>
typename Curve::ProjPoint msm_auto_parallel(
    const Curve &curve,
    const Ring &ring,
    const typename Curve::AffPoint *points,
    const uint8_t *const *scalars,
    size_t npoints,
    size_t scalar_bits = 255,
    size_t num_threads = 0
) {
    if (npoints < 32) {
        return msm(curve, ring, points, scalars, npoints, scalar_bits);
    }
    return msm_v2_auto_parallel(ring, points, scalars, npoints,
                                 scalar_bits, num_threads);
}
