#pragma once
#include "pippenger.hpp"
#include "scalar_mult.hpp"
#include <vector>

// Public MSM API with dispatch.
// Selects the best algorithm based on the number of points.
//
// G1 only for V1. G2 Pippenger deferred to V2 (needs Fp2 batch inversion).

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

    // Large count: Pippenger
    return pippenger_msm(curve, ring, points, scalars, npoints, scalar_bits);
}

// Multi-threaded MSM.
// Uses parallel point scatter within Pippenger windows.
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
    return pippenger_msm_parallel(curve, ring, points, scalars, npoints,
                                   scalar_bits, num_threads);
}
