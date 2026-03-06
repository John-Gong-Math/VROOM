#pragma once
#include "batch_inversion.hpp"
#include "ec.hpp"
#include <vector>

// Batch projective-to-affine conversion.
// Uses batch inversion to convert N projective points to affine
// with only one field inversion + O(n) multiplications.
//
// For projective point (X, Y, Z), affine = (X/Z, Y/Z).
// BLS12-381 uses Jacobian coordinates where affine = (X/Z^2, Y/Z^3).
// VROOM's PointAdd/PointDouble use the "complete" formulas from
// eprint.iacr.org/2015/1060 which use standard projective coordinates
// (not Jacobian), so affine = (X/Z, Y/Z).
//
// Wait - looking at ec.hpp more carefully: the zero point is (0, 1, 0)
// and Algorithm 7/9 from that paper use projective (not Jacobian) coordinates.
// Affine = (X/Z, Y/Z).
template<class Curve, class Ring>
std::vector<typename Curve::AffPoint> batch_to_affine(
    const Curve &/*curve*/,
    const Ring &ring,
    const std::vector<typename Curve::ProjPoint> &points
) {
    using StdElem = typename Ring::StandardElement;
    using AffPoint = typename Curve::AffPoint;

    size_t n = points.size();
    if (n == 0) return {};

    // Extract Z-coordinates
    std::vector<StdElem> z_coords(n);
    for (size_t i = 0; i < n; i++) {
        z_coords[i] = points[i].z;
    }

    // Batch invert all Z-coordinates
    std::vector<StdElem> z_invs = batch_invert(ring, z_coords);

    // Convert each point: x_affine = X * Z_inv, y_affine = Y * Z_inv
    std::vector<AffPoint> result(n);
    for (size_t i = 0; i < n; i++) {
        // Check if point is at infinity (Z == 0)
        auto z_bi = ring.to_bigint(points[i].z);
        if (z_bi == BigInt(0)) {
            // Return (0, 0) for infinity
            result[i] = AffPoint{ring.zero(), ring.zero()};
        } else {
            result[i].x = ring.modmul(points[i].x, z_invs[i]);
            result[i].y = ring.modmul(points[i].y, z_invs[i]);
        }
    }

    return result;
}
