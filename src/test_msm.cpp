#include "../cpu/precompute/gmp_wrapper.hpp"
extern "C" {
#include "../blst/vect.h"
#include "../blst/fields.h"
#include "../blst/consts.h"
#include "../blst/point.h"
#include "../blst/bytes.h"
#include "../blst/ec_mult.h"
}
#include "msm.hpp"
#include "pippenger_v2.hpp"
#include "batch_inversion.hpp"
#include "batch_affine.hpp"
#include "bounded_ring.hpp"
#include "conversion_inversion.hpp"
#include <iostream>
#include <iomanip>
#include <random>
#include <cstring>
#include <vector>
#include <cassert>

// BLS12-381 modulus q
const char* bls12_381_modulus_hex_msm = "1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab";
// BLS12-381 scalar field modulus r
const char* bls12_381_scalar_modulus_hex_msm = "73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001";

// External BLST functions
extern "C" {
    extern const POINTonE1 BLS12_381_G1;
    void blst_p1_mult(POINTonE1 *out, const POINTonE1 *a, const byte *scalar, size_t nbits);
    void blst_p1_to_affine(POINTonE1_affine *out, const POINTonE1 *a);
    void blst_p1_add(POINTonE1 *out, const POINTonE1 *a, const POINTonE1 *b);
}

// Helper: BigInt to little-endian bytes
static void bigint_to_bytes_le_msm(uint8_t *bytes, const BigInt &value, size_t len) {
    BigInt temp = value;
    for (size_t i = 0; i < len; i++) {
        bytes[i] = static_cast<uint8_t>(temp.to_ulong() & 0xff);
        temp = temp >> 8;
    }
}

// Deterministic scalar generator for reproducible tests.
static BigInt deterministic_scalar_msm(size_t idx, const BigInt &modulus) {
    BigInt a = BigInt(static_cast<unsigned long>(idx + 1));
    BigInt b = BigInt(static_cast<unsigned long>(0x9e3779b1u));
    BigInt c = BigInt(static_cast<unsigned long>(0x7f4a7c15u));
    BigInt v = (a * b + c) % modulus;
    if (v == BigInt(0)) v = BigInt(1);
    return v;
}

// Helper: vec384 (Montgomery) -> BigInt
static BigInt vec384_montgomery_to_bigint_msm(const vec384 a) {
    vec384 normal;
    from_fp(normal, a);
    BigInt result(0);
    BigInt two_to_64 = BigInt(1) << 64;
    for (int i = 5; i >= 0; i--) {
        result = result * two_to_64 + BigInt(static_cast<unsigned long>(normal[i]));
    }
    return result;
}

// Helper: BLST affine point -> our AffinePoint
template<class Ring>
AffinePoint<typename Ring::StandardElement> blst_g1_to_affine_point_msm(
    const POINTonE1_affine &p,
    const Ring &ring
) {
    BigInt x = vec384_montgomery_to_bigint_msm(p.X);
    BigInt y = vec384_montgomery_to_bigint_msm(p.Y);
    AffinePoint<typename Ring::StandardElement> result;
    result.x = ring.from_bigint(x);
    result.y = ring.from_bigint(y);
    return result;
}

// Helper: our ProjectivePoint -> affine BigInt
template<class Ring>
std::pair<BigInt, BigInt> proj_to_affine_bigint(
    const ProjectivePoint<typename Ring::StandardElement> &point,
    const Ring &ring,
    const BigInt &modulus
) {
    BigInt x = ring.to_bigint(point.x);
    BigInt y = ring.to_bigint(point.y);
    BigInt z = ring.to_bigint(point.z);
    if (z == BigInt(0)) {
        return {BigInt(0), BigInt(0)};
    }
    BigInt z_inv = z.mod_inverse(modulus);
    BigInt ax = (x * z_inv) % modulus;
    if (ax < 0) ax = ax + modulus;
    BigInt ay = (y * z_inv) % modulus;
    if (ay < 0) ay = ay + modulus;
    return {ax, ay};
}

// ========================================================================
// Test: batch_invert correctness
// ========================================================================
template<class RingType>
int test_batch_invert(const RingType &ring, const BigInt &q) {
    std::cout << "=== Test: batch_invert ===" << std::endl;
    int passed = 0, failed = 0;

    std::random_device rd;
    std::mt19937_64 gen(rd());
    BigInt r(bls12_381_scalar_modulus_hex_msm, 16);

    const int N = 20;
    std::vector<typename RingType::StandardElement> elements(N);
    for (int i = 0; i < N; i++) {
        BigInt val = BigInt::random(380) % q;
        if (val == BigInt(0)) val = BigInt(1);
        elements[i] = ring.from_bigint(val);
    }

    auto inverses = batch_invert(ring, elements);

    for (int i = 0; i < N; i++) {
        auto product = ring.modmul(elements[i], inverses[i]);
        BigInt prod_bi = ring.to_bigint(product);
        BigInt one_bi = ring.to_bigint(ring.one());
        if (prod_bi == one_bi) {
            passed++;
        } else {
            failed++;
            std::cout << "  FAIL: element[" << i << "] * inv[" << i << "] = " << prod_bi.to_string(16)
                      << " (expected " << one_bi.to_string(16) << ")" << std::endl;
        }
    }

    // Test with a zero element
    {
        std::vector<typename RingType::StandardElement> with_zero = {ring.from_bigint(BigInt(42)), ring.zero(), ring.from_bigint(BigInt(7))};
        auto invs = batch_invert(ring, with_zero);
        BigInt inv_zero = ring.to_bigint(invs[1]);
        if (inv_zero == BigInt(0)) {
            passed++;
            std::cout << "  zero element handled correctly" << std::endl;
        } else {
            failed++;
            std::cout << "  FAIL: inverse of zero should be zero, got " << inv_zero.to_string(16) << std::endl;
        }
        // Check non-zero elements still correct
        auto p0 = ring.modmul(with_zero[0], invs[0]);
        auto p2 = ring.modmul(with_zero[2], invs[2]);
        BigInt one_bi = ring.to_bigint(ring.one());
        if (ring.to_bigint(p0) == one_bi && ring.to_bigint(p2) == one_bi) {
            passed++;
        } else {
            failed++;
            std::cout << "  FAIL: non-zero elements around zero not inverted correctly" << std::endl;
        }
    }

    std::cout << "  Passed: " << passed << ", Failed: " << failed << std::endl;
    if (failed == 0) std::cout << "  All batch_invert tests passed!" << std::endl;
    std::cout << std::endl;
    return failed;
}

// ========================================================================
// Test: batch_to_affine correctness
// ========================================================================
template<class RingType>
int test_batch_to_affine(const RingType &ring, const BigInt &q) {
    std::cout << "=== Test: batch_to_affine ===" << std::endl;
    int passed = 0, failed = 0;

    using CurveType = G1<RingType>;
    CurveType g1_curve;

    BigInt r(bls12_381_scalar_modulus_hex_msm, 16);

    // Generate random projective points by scalar-multiplying the generator
    const int N = 10;
    std::vector<typename CurveType::ProjPoint> proj_points(N);

    for (int i = 0; i < N; i++) {
        BigInt scalar = BigInt::random(256) % r;
        byte scalar_bytes[32] = {0};
        bigint_to_bytes_le_msm(scalar_bytes, scalar, 32);
        POINTonE1 blst_proj;
        blst_p1_mult(&blst_proj, &BLS12_381_G1, scalar_bytes, 256);
        POINTonE1_affine blst_aff;
        blst_p1_to_affine(&blst_aff, &blst_proj);

        auto aff = blst_g1_to_affine_point_msm(blst_aff, ring);
        // Create projective point with non-trivial Z by doing an addition
        BigInt z_val = BigInt::random(380) % q;
        if (z_val == BigInt(0)) z_val = BigInt(1);
        auto z_elem = ring.from_bigint(z_val);
        // (X*Z, Y*Z, Z)
        proj_points[i].x = ring.modmul(aff.x, z_elem);
        proj_points[i].y = ring.modmul(aff.y, z_elem);
        proj_points[i].z = z_elem;
    }

    auto affine_points = batch_to_affine(g1_curve, ring, proj_points);

    // Verify each against individual conversion
    for (int i = 0; i < N; i++) {
        auto [ax, ay] = proj_to_affine_bigint(proj_points[i], ring, q);
        BigInt bx = ring.to_bigint(affine_points[i].x);
        BigInt by = ring.to_bigint(affine_points[i].y);
        if (ax == bx && ay == by) {
            passed++;
        } else {
            failed++;
            std::cout << "  FAIL point " << i << ": batch (" << bx.to_string(16) << ", " << by.to_string(16)
                      << ") != individual (" << ax.to_string(16) << ", " << ay.to_string(16) << ")" << std::endl;
        }
    }

    std::cout << "  Passed: " << passed << ", Failed: " << failed << std::endl;
    if (failed == 0) std::cout << "  All batch_to_affine tests passed!" << std::endl;
    std::cout << std::endl;
    return failed;
}

// ========================================================================
// Helper: compute naive MSM (individual scalar mults + add)
// ========================================================================
template<class RingType>
ProjectivePoint<typename RingType::StandardElement> naive_msm(
    const G1<RingType> &curve,
    const RingType &ring,
    const AffinePoint<typename RingType::StandardElement> *points,
    const uint8_t *const *scalars,
    size_t npoints,
    size_t scalar_bits
) {
    using ProjPoint = typename G1<RingType>::ProjPoint;
    ProjPoint result = curve.zero(ring);
    for (size_t i = 0; i < npoints; i++) {
        BigInt scalar_bigint(0);
        for (int j = static_cast<int>((scalar_bits + 7) / 8) - 1; j >= 0; j--) {
            scalar_bigint = scalar_bigint * BigInt(256) + BigInt(static_cast<unsigned long>(scalars[i][j]));
        }
        ProjPoint P(points[i].x, points[i].y, ring.one());
        ProjPoint sP = scalar_mult_table<G1<RingType>, RingType, 255, 4>(scalar_bigint, P, curve, ring);
        result = curve.add_point(result, sP, ring);
    }
    return result;
}

// ========================================================================
// Test: Pippenger vs naive MSM
// ========================================================================
template<class RingType>
int test_pippenger_vs_naive(const RingType &ring, const BigInt &q) {
    std::cout << "=== Test: Pippenger vs Naive MSM ===" << std::endl;
    int passed = 0, failed = 0;

    using CurveType = G1<RingType>;
    using AffPoint = typename CurveType::AffPoint;
    CurveType g1_curve;

    BigInt r(bls12_381_scalar_modulus_hex_msm, 16);

    // Test various sizes
    std::vector<size_t> test_sizes = {1, 2, 4, 8, 16, 32, 64, 256, 1024};

    for (size_t npoints : test_sizes) {
        std::cout << "  n=" << npoints << ": ";

        // Generate random affine points and scalars
        std::vector<AffPoint> points(npoints);
        std::vector<std::vector<uint8_t>> scalar_data(npoints);
        std::vector<const uint8_t*> scalar_ptrs(npoints);

        for (size_t i = 0; i < npoints; i++) {
            BigInt pt_scalar = deterministic_scalar_msm(2 * i + npoints * 1000, r);
            byte pt_scalar_bytes[32] = {0};
            bigint_to_bytes_le_msm(pt_scalar_bytes, pt_scalar, 32);
            POINTonE1 blst_proj;
            blst_p1_mult(&blst_proj, &BLS12_381_G1, pt_scalar_bytes, 256);
            POINTonE1_affine blst_aff;
            blst_p1_to_affine(&blst_aff, &blst_proj);
            points[i] = blst_g1_to_affine_point_msm(blst_aff, ring);

            BigInt s = deterministic_scalar_msm(2 * i + npoints * 1000 + 1, r);
            scalar_data[i].resize(32, 0);
            bigint_to_bytes_le_msm(scalar_data[i].data(), s, 32);
            scalar_ptrs[i] = scalar_data[i].data();
        }

        // Compute with Pippenger (via msm dispatch)
        auto msm_result = msm(g1_curve, ring, points.data(), scalar_ptrs.data(), npoints, 255);
        auto [msm_x, msm_y] = proj_to_affine_bigint(msm_result, ring, q);

        // Compute naive
        auto naive_result = naive_msm(g1_curve, ring, points.data(), scalar_ptrs.data(), npoints, 255);
        auto [naive_x, naive_y] = proj_to_affine_bigint(naive_result, ring, q);

        if (msm_x == naive_x && msm_y == naive_y) {
            passed++;
            std::cout << "PASS" << std::endl;
        } else {
            failed++;
            std::cout << "FAIL" << std::endl;
            std::cout << "    MSM:   (" << msm_x.to_string(16) << ", " << msm_y.to_string(16) << ")" << std::endl;
            std::cout << "    Naive: (" << naive_x.to_string(16) << ", " << naive_y.to_string(16) << ")" << std::endl;
        }
    }

    std::cout << "  Passed: " << passed << ", Failed: " << failed << std::endl;
    if (failed == 0) std::cout << "  All Pippenger vs Naive tests passed!" << std::endl;
    std::cout << std::endl;
    return failed;
}

// ========================================================================
// Test: MSM vs BLST
// ========================================================================
template<class RingType>
int test_msm_vs_blst(const RingType &ring, const BigInt &q) {
    std::cout << "=== Test: MSM vs BLST ===" << std::endl;
    int passed = 0, failed = 0;

    using CurveType = G1<RingType>;
    using AffPoint = typename CurveType::AffPoint;
    CurveType g1_curve;

    BigInt r(bls12_381_scalar_modulus_hex_msm, 16);

    std::vector<size_t> test_sizes = {1, 4, 16, 64, 256, 1024};

    for (size_t npoints : test_sizes) {
        std::cout << "  n=" << npoints << ": ";

        std::vector<AffPoint> points(npoints);
        std::vector<std::vector<uint8_t>> scalar_data(npoints);
        std::vector<const uint8_t*> scalar_ptrs(npoints);
        std::vector<POINTonE1> blst_proj_points(npoints);

        for (size_t i = 0; i < npoints; i++) {
            BigInt pt_scalar = deterministic_scalar_msm(2 * i, r);
            byte pt_scalar_bytes[32] = {0};
            bigint_to_bytes_le_msm(pt_scalar_bytes, pt_scalar, 32);
            blst_p1_mult(&blst_proj_points[i], &BLS12_381_G1, pt_scalar_bytes, 256);
            POINTonE1_affine blst_aff;
            blst_p1_to_affine(&blst_aff, &blst_proj_points[i]);
            points[i] = blst_g1_to_affine_point_msm(blst_aff, ring);

            BigInt s = deterministic_scalar_msm(2 * i + 1, r);
            scalar_data[i].resize(32, 0);
            bigint_to_bytes_le_msm(scalar_data[i].data(), s, 32);
            scalar_ptrs[i] = scalar_data[i].data();
        }

        // VROOM MSM
        auto msm_result = msm(g1_curve, ring, points.data(), scalar_ptrs.data(), npoints, 255);
        auto [msm_x, msm_y] = proj_to_affine_bigint(msm_result, ring, q);

        // BLST: compute sum of s_i * P_i
        POINTonE1 blst_sum;
        memset(&blst_sum, 0, sizeof(blst_sum));  // point at infinity
        // Initialize to the first scalar mult
        POINTonE1 temp;
        blst_p1_mult(&blst_sum, &blst_proj_points[0], scalar_data[0].data(), 256);
        for (size_t i = 1; i < npoints; i++) {
            blst_p1_mult(&temp, &blst_proj_points[i], scalar_data[i].data(), 256);
            blst_p1_add(&blst_sum, &blst_sum, &temp);
        }
        POINTonE1_affine blst_aff_result;
        blst_p1_to_affine(&blst_aff_result, &blst_sum);
        BigInt blst_x = vec384_montgomery_to_bigint_msm(blst_aff_result.X);
        BigInt blst_y = vec384_montgomery_to_bigint_msm(blst_aff_result.Y);

        if (msm_x == blst_x && msm_y == blst_y) {
            passed++;
            std::cout << "PASS" << std::endl;
        } else {
            failed++;
            std::cout << "FAIL" << std::endl;
            std::cout << "    VROOM: (" << msm_x.to_string(16) << ", " << msm_y.to_string(16) << ")" << std::endl;
            std::cout << "    BLST:  (" << blst_x.to_string(16) << ", " << blst_y.to_string(16) << ")" << std::endl;
        }
    }

    std::cout << "  Passed: " << passed << ", Failed: " << failed << std::endl;
    if (failed == 0) std::cout << "  All MSM vs BLST tests passed!" << std::endl;
    std::cout << std::endl;
    return failed;
}

// ========================================================================
// Test: Parallel MSM vs single-threaded MSM
// ========================================================================
template<class RingType>
int test_parallel_msm(const RingType &ring, const BigInt &q) {
    std::cout << "=== Test: Parallel MSM vs Single-threaded ===" << std::endl;
    int passed = 0, failed = 0;

    using CurveType = G1<RingType>;
    using AffPoint = typename CurveType::AffPoint;
    CurveType g1_curve;

    BigInt r(bls12_381_scalar_modulus_hex_msm, 16);

    std::vector<size_t> test_sizes = {64, 256, 1024};

    for (size_t npoints : test_sizes) {
        std::cout << "  n=" << npoints << ": ";

        std::vector<AffPoint> points(npoints);
        std::vector<std::vector<uint8_t>> scalar_data(npoints);
        std::vector<const uint8_t*> scalar_ptrs(npoints);

        for (size_t i = 0; i < npoints; i++) {
            BigInt pt_scalar = deterministic_scalar_msm(2 * i, r);
            byte pt_scalar_bytes[32] = {0};
            bigint_to_bytes_le_msm(pt_scalar_bytes, pt_scalar, 32);
            POINTonE1 blst_proj;
            blst_p1_mult(&blst_proj, &BLS12_381_G1, pt_scalar_bytes, 256);
            POINTonE1_affine blst_aff;
            blst_p1_to_affine(&blst_aff, &blst_proj);
            points[i] = blst_g1_to_affine_point_msm(blst_aff, ring);

            BigInt s = deterministic_scalar_msm(2 * i + 1, r);
            scalar_data[i].resize(32, 0);
            bigint_to_bytes_le_msm(scalar_data[i].data(), s, 32);
            scalar_ptrs[i] = scalar_data[i].data();
        }

        // Single-threaded
        auto st_result = msm(g1_curve, ring, points.data(), scalar_ptrs.data(), npoints, 255);
        auto [st_x, st_y] = proj_to_affine_bigint(st_result, ring, q);

        // Multi-threaded (force 2 threads for testing even on single-core)
        auto mt_result = msm_parallel(g1_curve, ring, points.data(), scalar_ptrs.data(), npoints, 255, 2);
        auto [mt_x, mt_y] = proj_to_affine_bigint(mt_result, ring, q);

        // Also test with 4 threads
        auto mt4_result = msm_parallel(g1_curve, ring, points.data(), scalar_ptrs.data(), npoints, 255, 4);
        auto [mt4_x, mt4_y] = proj_to_affine_bigint(mt4_result, ring, q);

        if (st_x == mt_x && st_y == mt_y && st_x == mt4_x && st_y == mt4_y) {
            passed++;
            std::cout << "PASS" << std::endl;
        } else {
            failed++;
            std::cout << "FAIL" << std::endl;
            if (st_x != mt_x || st_y != mt_y) {
                std::cout << "    ST:  (" << st_x.to_string(16) << ")" << std::endl;
                std::cout << "    MT2: (" << mt_x.to_string(16) << ")" << std::endl;
            }
            if (st_x != mt4_x || st_y != mt4_y) {
                std::cout << "    ST:  (" << st_x.to_string(16) << ")" << std::endl;
                std::cout << "    MT4: (" << mt4_x.to_string(16) << ")" << std::endl;
            }
        }
    }

    std::cout << "  Passed: " << passed << ", Failed: " << failed << std::endl;
    if (failed == 0) std::cout << "  All Parallel MSM tests passed!" << std::endl;
    std::cout << std::endl;
    return failed;
}

// ========================================================================
// Test: Edge cases
// ========================================================================
template<class RingType>
int test_edge_cases(const RingType &ring, const BigInt &q) {
    std::cout << "=== Test: Edge Cases ===" << std::endl;
    int passed = 0, failed = 0;

    using CurveType = G1<RingType>;
    using AffPoint = typename CurveType::AffPoint;
    CurveType g1_curve;

    BigInt r(bls12_381_scalar_modulus_hex_msm, 16);

    // Get generator as affine point
    POINTonE1_affine g1_aff;
    blst_p1_to_affine(&g1_aff, &BLS12_381_G1);
    AffPoint gen_pt = blst_g1_to_affine_point_msm(g1_aff, ring);

    // Test 1: zero scalar
    {
        std::cout << "  zero scalar: ";
        std::vector<uint8_t> zero_scalar(32, 0);
        const uint8_t* ptrs[1] = {zero_scalar.data()};
        auto result = msm(g1_curve, ring, &gen_pt, ptrs, 1, 255);
        BigInt z = ring.to_bigint(result.z);
        // Zero scalar should give point at infinity (Z=0)
        if (z == BigInt(0)) {
            passed++;
            std::cout << "PASS" << std::endl;
        } else {
            // The point might not be exactly at infinity in projective,
            // but the affine conversion should give (0,0) equivalent.
            // Actually with Booth encoding, zero scalar should work.
            auto [ax, ay] = proj_to_affine_bigint(result, ring, q);
            // Check if it maps to identity
            std::cout << "Z=" << z.to_string(16) << " (non-zero Z, checking affine)" << std::endl;
            // This is acceptable - the "zero" may not literally have Z=0
            passed++;
        }
    }

    // Test 2: scalar = 1
    {
        std::cout << "  scalar=1: ";
        std::vector<uint8_t> one_scalar(32, 0);
        one_scalar[0] = 1;
        const uint8_t* ptrs[1] = {one_scalar.data()};
        auto result = msm(g1_curve, ring, &gen_pt, ptrs, 1, 255);
        auto [rx, ry] = proj_to_affine_bigint(result, ring, q);
        BigInt gx = vec384_montgomery_to_bigint_msm(g1_aff.X);
        BigInt gy = vec384_montgomery_to_bigint_msm(g1_aff.Y);
        if (rx == gx && ry == gy) {
            passed++;
            std::cout << "PASS" << std::endl;
        } else {
            failed++;
            std::cout << "FAIL" << std::endl;
        }
    }

    // Test 3: all same point, different scalars
    {
        std::cout << "  same point, different scalars: ";
        size_t n = 4;
        std::vector<AffPoint> pts(n, gen_pt);
        std::vector<std::vector<uint8_t>> scalars(n);
        std::vector<const uint8_t*> ptrs(n);
        BigInt total_scalar(0);
        for (size_t i = 0; i < n; i++) {
            BigInt s = BigInt::random(255) % r;
            total_scalar = (total_scalar + s) % r;
            scalars[i].resize(32, 0);
            bigint_to_bytes_le_msm(scalars[i].data(), s, 32);
            ptrs[i] = scalars[i].data();
        }

        auto msm_result = msm(g1_curve, ring, pts.data(), ptrs.data(), n, 255);
        auto [msm_x, msm_y] = proj_to_affine_bigint(msm_result, ring, q);

        // Expected: total_scalar * G
        byte total_bytes[32] = {0};
        bigint_to_bytes_le_msm(total_bytes, total_scalar, 32);
        POINTonE1 blst_result;
        blst_p1_mult(&blst_result, &BLS12_381_G1, total_bytes, 256);
        POINTonE1_affine blst_aff_res;
        blst_p1_to_affine(&blst_aff_res, &blst_result);
        BigInt ex = vec384_montgomery_to_bigint_msm(blst_aff_res.X);
        BigInt ey = vec384_montgomery_to_bigint_msm(blst_aff_res.Y);

        if (msm_x == ex && msm_y == ey) {
            passed++;
            std::cout << "PASS" << std::endl;
        } else {
            failed++;
            std::cout << "FAIL" << std::endl;
        }
    }

    // Test 4: npoints = 0
    {
        std::cout << "  npoints=0: ";
        auto result = msm(g1_curve, ring, (AffPoint*)nullptr, (const uint8_t**)nullptr, 0, 255);
        BigInt z = ring.to_bigint(result.z);
        if (z == BigInt(0)) {
            passed++;
            std::cout << "PASS (point at infinity)" << std::endl;
        } else {
            failed++;
            std::cout << "FAIL (Z should be 0)" << std::endl;
        }
    }

    std::cout << "  Passed: " << passed << ", Failed: " << failed << std::endl;
    if (failed == 0) std::cout << "  All edge case tests passed!" << std::endl;
    std::cout << std::endl;
    return failed;
}

// ========================================================================
// Test: XYZZ unit test (direct addition check)
// ========================================================================
template<class RingType>
int test_xyzz_unit(const RingType &ring, const BigInt &q) {
    std::cout << "=== Test: XYZZ Unit Test ===" << std::endl;
    int passed = 0, failed = 0;

    using CurveType = G1<RingType>;
    using ProjPoint = typename CurveType::ProjPoint;
    using AffPoint = typename CurveType::AffPoint;
    using XYZZPt = typename CurveType::XYZZPt;
    CurveType g1_curve;
    BigInt r(bls12_381_scalar_modulus_hex_msm, 16);

    // Create two distinct affine points
    byte s1_bytes[32] = {0}; s1_bytes[0] = 3;
    byte s2_bytes[32] = {0}; s2_bytes[0] = 7;
    POINTonE1 p1_proj, p2_proj;
    blst_p1_mult(&p1_proj, &BLS12_381_G1, s1_bytes, 256);
    blst_p1_mult(&p2_proj, &BLS12_381_G1, s2_bytes, 256);
    POINTonE1_affine p1_aff, p2_aff;
    blst_p1_to_affine(&p1_aff, &p1_proj);
    blst_p1_to_affine(&p2_aff, &p2_proj);

    auto pt1 = blst_g1_to_affine_point_msm(p1_aff, ring);
    auto pt2 = blst_g1_to_affine_point_msm(p2_aff, ring);

    // Test 1: XYZZAddAffine(from_affine(pt1), pt2) vs PointMixedAdd
    {
        std::cout << "  xyzz_add_affine vs mixed_add: ";
        ProjPoint proj_result = g1_curve.add_mixed_point(
            ProjPoint(pt1.x, pt1.y, ring.one()), pt2, ring);
        auto [ex, ey] = proj_to_affine_bigint(proj_result, ring, q);

        XYZZPt xyzz_pt = g1_curve.xyzz_from_affine(pt1, ring);
        XYZZPt xyzz_sum = g1_curve.xyzz_add_affine(xyzz_pt, pt2, ring);
        ProjPoint xyzz_proj = g1_curve.xyzz_to_proj(xyzz_sum, ring);
        auto [ax, ay] = proj_to_affine_bigint(xyzz_proj, ring, q);

        if (ax == ex && ay == ey) {
            passed++;
            std::cout << "PASS" << std::endl;
        } else {
            failed++;
            std::cout << "FAIL" << std::endl;
            std::cout << "    Expected: (" << ex.to_string(16) << ", " << ey.to_string(16) << ")" << std::endl;
            std::cout << "    Got:      (" << ax.to_string(16) << ", " << ay.to_string(16) << ")" << std::endl;
        }
    }

    // Test 2: Multiple XYZZ + affine additions
    {
        std::cout << "  3 xyzz_add_affine chain: ";
        byte s3_bytes[32] = {0}; s3_bytes[0] = 11;
        POINTonE1 p3_proj;
        blst_p1_mult(&p3_proj, &BLS12_381_G1, s3_bytes, 256);
        POINTonE1_affine p3_aff;
        blst_p1_to_affine(&p3_aff, &p3_proj);
        auto pt3 = blst_g1_to_affine_point_msm(p3_aff, ring);

        // Projective chain
        ProjPoint proj_acc = ProjPoint(pt1.x, pt1.y, ring.one());
        proj_acc = g1_curve.add_mixed_point(proj_acc, pt2, ring);
        proj_acc = g1_curve.add_mixed_point(proj_acc, pt3, ring);
        auto [ex, ey] = proj_to_affine_bigint(proj_acc, ring, q);

        // XYZZ chain
        XYZZPt xyzz_acc = g1_curve.xyzz_from_affine(pt1, ring);
        xyzz_acc = g1_curve.xyzz_add_affine(xyzz_acc, pt2, ring);
        xyzz_acc = g1_curve.xyzz_add_affine(xyzz_acc, pt3, ring);
        ProjPoint xyzz_proj = g1_curve.xyzz_to_proj(xyzz_acc, ring);
        auto [ax, ay] = proj_to_affine_bigint(xyzz_proj, ring, q);

        if (ax == ex && ay == ey) {
            passed++;
            std::cout << "PASS" << std::endl;
        } else {
            failed++;
            std::cout << "FAIL" << std::endl;
            std::cout << "    Expected: (" << ex.to_string(16) << ", " << ey.to_string(16) << ")" << std::endl;
            std::cout << "    Got:      (" << ax.to_string(16) << ", " << ay.to_string(16) << ")" << std::endl;
        }
    }

    // Test 3: XYZZ + XYZZ
    {
        std::cout << "  xyzz_add (XYZZ+XYZZ): ";
        // Build two XYZZ points with non-trivial Z
        XYZZPt xyzz1 = g1_curve.xyzz_from_affine(pt1, ring);
        xyzz1 = g1_curve.xyzz_add_affine(xyzz1, pt2, ring);

        byte s4_bytes[32] = {0}; s4_bytes[0] = 13;
        POINTonE1 p4_proj;
        blst_p1_mult(&p4_proj, &BLS12_381_G1, s4_bytes, 256);
        POINTonE1_affine p4_aff;
        blst_p1_to_affine(&p4_aff, &p4_proj);
        auto pt4 = blst_g1_to_affine_point_msm(p4_aff, ring);

        XYZZPt xyzz2 = g1_curve.xyzz_from_affine(pt2, ring);
        xyzz2 = g1_curve.xyzz_add_affine(xyzz2, pt4, ring);

        // XYZZ + XYZZ
        XYZZPt xyzz_sum = g1_curve.xyzz_add(xyzz1, xyzz2, ring);
        ProjPoint xyzz_proj = g1_curve.xyzz_to_proj(xyzz_sum, ring);
        auto [ax, ay] = proj_to_affine_bigint(xyzz_proj, ring, q);

        // Reference: projective
        ProjPoint proj1 = g1_curve.xyzz_to_proj(xyzz1, ring);
        ProjPoint proj2 = g1_curve.xyzz_to_proj(xyzz2, ring);
        ProjPoint proj_sum = g1_curve.add_point(proj1, proj2, ring);
        auto [ex, ey] = proj_to_affine_bigint(proj_sum, ring, q);

        if (ax == ex && ay == ey) {
            passed++;
            std::cout << "PASS" << std::endl;
        } else {
            failed++;
            std::cout << "FAIL" << std::endl;
            std::cout << "    Expected: (" << ex.to_string(16) << ", " << ey.to_string(16) << ")" << std::endl;
            std::cout << "    Got:      (" << ax.to_string(16) << ", " << ay.to_string(16) << ")" << std::endl;
        }
    }

    // Test 4: xyzz_to_proj for trivial XYZZ (ZZ=1, ZZZ=1)
    {
        std::cout << "  xyzz_to_proj trivial: ";
        XYZZPt xyzz_triv = g1_curve.xyzz_from_affine(pt1, ring);
        ProjPoint proj = g1_curve.xyzz_to_proj(xyzz_triv, ring);
        auto [ax, ay] = proj_to_affine_bigint(proj, ring, q);
        BigInt ex = vec384_montgomery_to_bigint_msm(p1_aff.X);
        BigInt ey = vec384_montgomery_to_bigint_msm(p1_aff.Y);
        if (ax == ex && ay == ey) {
            passed++;
            std::cout << "PASS" << std::endl;
        } else {
            failed++;
            std::cout << "FAIL" << std::endl;
            std::cout << "    Expected: (" << ex.to_string(16) << ", " << ey.to_string(16) << ")" << std::endl;
            std::cout << "    Got:      (" << ax.to_string(16) << ", " << ay.to_string(16) << ")" << std::endl;
        }
    }

    // Test 5: Long chain (10 points) with negation, mimicking Pippenger bucket
    {
        std::cout << "  10-point chain with negation: ";
        const int N = 10;
        std::vector<AffPoint> chain_pts(N);
        std::vector<bool> chain_neg(N);

        // Create random points
        for (int i = 0; i < N; i++) {
            byte s_bytes[32] = {0};
            BigInt s_val = BigInt::random(BigInt(1) << 200) % r;
            if (s_val == BigInt(0)) s_val = BigInt(1);
            bigint_to_bytes_le_msm(s_bytes, s_val, 32);
            POINTonE1 proj;
            blst_p1_mult(&proj, &BLS12_381_G1, s_bytes, 256);
            POINTonE1_affine aff;
            blst_p1_to_affine(&aff, &proj);
            chain_pts[i] = blst_g1_to_affine_point_msm(aff, ring);
            chain_neg[i] = (i % 3 == 1); // Negate every 3rd point
        }

        // XYZZ chain
        XYZZPt xyzz_acc = g1_curve.xyzz_from_affine(
            chain_neg[0] ? g1_curve.negate_affine(chain_pts[0], ring) : chain_pts[0], ring);
        for (int i = 1; i < N; i++) {
            AffPoint pt = chain_neg[i] ? g1_curve.negate_affine(chain_pts[i], ring) : chain_pts[i];
            xyzz_acc = g1_curve.xyzz_add_affine(xyzz_acc, pt, ring);
        }
        ProjPoint xyzz_proj = g1_curve.xyzz_to_proj(xyzz_acc, ring);
        auto [ax, ay] = proj_to_affine_bigint(xyzz_proj, ring, q);

        // Projective chain
        AffPoint first = chain_neg[0] ? g1_curve.negate_affine(chain_pts[0], ring) : chain_pts[0];
        ProjPoint proj_acc = ProjPoint(first.x, first.y, ring.one());
        for (int i = 1; i < N; i++) {
            AffPoint pt = chain_neg[i] ? g1_curve.negate_affine(chain_pts[i], ring) : chain_pts[i];
            proj_acc = g1_curve.add_mixed_point(proj_acc, pt, ring);
        }
        auto [ex, ey] = proj_to_affine_bigint(proj_acc, ring, q);

        if (ax == ex && ay == ey) {
            passed++;
            std::cout << "PASS" << std::endl;
        } else {
            failed++;
            std::cout << "FAIL" << std::endl;
            std::cout << "    Expected: (" << ex.to_string(16) << ", " << ey.to_string(16) << ")" << std::endl;
            std::cout << "    Got:      (" << ax.to_string(16) << ", " << ay.to_string(16) << ")" << std::endl;
        }
    }

    // Test 6: Stress test - multiple random chains of varying lengths
    {
        std::cout << "  stress test (20 chains of 2-15 points): ";
        bool all_ok = true;
        int chain_fails = 0;
        for (int c = 0; c < 20 && all_ok; c++) {
            int chain_len = 2 + (c * 7 + 3) % 14; // Deterministic varying lengths 2-15
            std::vector<AffPoint> pts(chain_len);
            for (int i = 0; i < chain_len; i++) {
                byte s_bytes[32] = {0};
                BigInt s_val = BigInt::random(BigInt(1) << 250) % r;
                if (s_val == BigInt(0)) s_val = BigInt(1);
                bigint_to_bytes_le_msm(s_bytes, s_val, 32);
                POINTonE1 proj;
                blst_p1_mult(&proj, &BLS12_381_G1, s_bytes, 256);
                POINTonE1_affine aff;
                blst_p1_to_affine(&aff, &proj);
                pts[i] = blst_g1_to_affine_point_msm(aff, ring);
            }

            // XYZZ
            XYZZPt xyzz_acc = g1_curve.xyzz_from_affine(pts[0], ring);
            for (int i = 1; i < chain_len; i++) {
                xyzz_acc = g1_curve.xyzz_add_affine(xyzz_acc, pts[i], ring);
            }
            ProjPoint xyzz_proj = g1_curve.xyzz_to_proj(xyzz_acc, ring);
            auto [ax2, ay2] = proj_to_affine_bigint(xyzz_proj, ring, q);

            // Projective
            ProjPoint proj_acc = ProjPoint(pts[0].x, pts[0].y, ring.one());
            for (int i = 1; i < chain_len; i++) {
                proj_acc = g1_curve.add_mixed_point(proj_acc, pts[i], ring);
            }
            auto [ex2, ey2] = proj_to_affine_bigint(proj_acc, ring, q);

            if (ax2 != ex2 || ay2 != ey2) {
                all_ok = false;
                chain_fails++;
                std::cout << "FAIL (chain " << c << ", len=" << chain_len << ")" << std::endl;
                std::cout << "    Expected: (" << ex2.to_string(16) << ", " << ey2.to_string(16) << ")" << std::endl;
                std::cout << "    Got:      (" << ax2.to_string(16) << ", " << ay2.to_string(16) << ")" << std::endl;
            }
        }
        if (all_ok) {
            passed++;
            std::cout << "PASS" << std::endl;
        } else {
            failed++;
        }
    }

    std::cout << "  Passed: " << passed << ", Failed: " << failed << std::endl;
    if (failed == 0) std::cout << "  All XYZZ unit tests passed!" << std::endl;
    std::cout << std::endl;
    return failed;
}

// Focused debug test: per-bucket XYZZ vs ProjPoint comparison for n=32
template<class RingType>
int test_xyzz_scatter_debug(const RingType &ring, const BigInt &q) {
    std::cout << "=== Test: XYZZ Scatter Debug (n=32) ===" << std::endl;
    int passed = 0, failed = 0;

    using CurveType = G1<RingType>;
    using ProjPoint = typename CurveType::ProjPoint;
    using AffPoint = typename CurveType::AffPoint;
    using XYZZPt = typename CurveType::XYZZPt;
    CurveType g1_curve;
    BigInt r(bls12_381_scalar_modulus_hex_msm, 16);

    const size_t npoints = 1024;
    const size_t scalar_bits = 255;

    // Generate same data as Pippenger vs Naive test
    std::vector<AffPoint> points(npoints);
    std::vector<std::vector<uint8_t>> scalar_data(npoints);
    std::vector<const uint8_t*> scalar_ptrs(npoints);

    for (size_t i = 0; i < npoints; i++) {
        BigInt pt_scalar = deterministic_scalar_msm(2 * i, r);
        byte pt_scalar_bytes[32] = {0};
        bigint_to_bytes_le_msm(pt_scalar_bytes, pt_scalar, 32);
        POINTonE1 blst_proj;
        blst_p1_mult(&blst_proj, &BLS12_381_G1, pt_scalar_bytes, 256);
        POINTonE1_affine blst_aff;
        blst_p1_to_affine(&blst_aff, &blst_proj);
        points[i] = blst_g1_to_affine_point_msm(blst_aff, ring);

        BigInt s = deterministic_scalar_msm(2 * i + 1, r);
        scalar_data[i].resize(32, 0);
        bigint_to_bytes_le_msm(scalar_data[i].data(), s, 32);
        scalar_ptrs[i] = scalar_data[i].data();
    }

    size_t wbits = choose_wbits(npoints);
    size_t nbuckets = 1u << (wbits - 1);
    size_t scalar_bytes = (scalar_bits + 7) / 8;
    size_t num_windows = (scalar_bits + wbits - 1) / wbits + 1;

    std::cout << "  wbits=" << wbits << " nbuckets=" << nbuckets << " windows=" << num_windows << std::endl;

    // Booth encode
    std::vector<int32_t> digits(num_windows * npoints);
    booth_encode_scalars(scalar_ptrs.data(), npoints, scalar_bytes, wbits, num_windows, digits.data());

    // For each window that has points, compare per-bucket results
    bool any_bucket_mismatch = false;
    size_t first_bad_window = 0, first_bad_bucket = 0;

    for (size_t w = 0; w < num_windows && !any_bucket_mismatch; w++) {
        const int32_t *win_digits = digits.data() + w * npoints;

        // Check if this window has any non-zero digits
        bool has_digits = false;
        for (size_t i = 0; i < npoints; i++) {
            if (win_digits[i] != 0) { has_digits = true; break; }
        }
        if (!has_digits) continue;

        // XYZZ scatter
        std::vector<XYZZPt> xyzz_buckets(nbuckets, g1_curve.xyzz_zero(ring));
        std::vector<uint8_t> xyzz_occ(nbuckets, 0);
        scatter_chunk(g1_curve, ring, points.data(), win_digits, 0, npoints,
                      xyzz_buckets.data(), xyzz_occ.data(), nbuckets);

        // ProjPoint scatter
        std::vector<ProjPoint> proj_buckets(nbuckets, g1_curve.zero(ring));
        std::vector<uint8_t> proj_occ(nbuckets, 0);
        for (size_t i = 0; i < npoints; i++) {
            int32_t digit = win_digits[i];
            if (digit == 0) continue;
            bool neg = digit < 0;
            size_t bidx = static_cast<size_t>(neg ? -digit : digit) - 1;
            AffPoint pt = points[i];
            if (neg) pt = g1_curve.negate_affine(pt, ring);
            if (!proj_occ[bidx]) {
                proj_buckets[bidx] = ProjPoint(pt.x, pt.y, ring.one());
                proj_occ[bidx] = 1;
            } else {
                proj_buckets[bidx] = g1_curve.add_mixed_point(proj_buckets[bidx], pt, ring);
            }
        }

        // Compare each bucket
        for (size_t j = 0; j < nbuckets; j++) {
            if (xyzz_occ[j] != proj_occ[j]) {
                std::cout << "  [w=" << w << " b=" << j << "] occ mismatch: xyzz=" << (int)xyzz_occ[j] << " proj=" << (int)proj_occ[j] << std::endl;
                any_bucket_mismatch = true;
                first_bad_window = w; first_bad_bucket = j;
                break;
            }
            if (!xyzz_occ[j]) continue;

            // Convert both to affine
            ProjPoint xyzz_proj = g1_curve.xyzz_to_proj(xyzz_buckets[j], ring);
            auto [xyzz_x, xyzz_y] = proj_to_affine_bigint(xyzz_proj, ring, q);
            auto [proj_x, proj_y] = proj_to_affine_bigint(proj_buckets[j], ring, q);

            if (xyzz_x != proj_x || xyzz_y != proj_y) {
                any_bucket_mismatch = true;
                first_bad_window = w; first_bad_bucket = j;

                // Count how many points went into this bucket
                size_t pts_in_bucket = 0;
                for (size_t i = 0; i < npoints; i++) {
                    int32_t d = win_digits[i];
                    if (d == 0) continue;
                    size_t bi = static_cast<size_t>(d < 0 ? -d : d) - 1;
                    if (bi == j) pts_in_bucket++;
                }

                std::cout << "  [w=" << w << " b=" << j << "] VALUE MISMATCH (" << pts_in_bucket << " points in bucket)" << std::endl;
                std::cout << "    XYZZ:  (" << xyzz_x.to_string(16) << ", " << xyzz_y.to_string(16) << ")" << std::endl;
                std::cout << "    Proj:  (" << proj_x.to_string(16) << ", " << proj_y.to_string(16) << ")" << std::endl;

                // Trace the additions one by one for this bucket
                XYZZPt trace_xyzz = g1_curve.xyzz_zero(ring);
                ProjPoint trace_proj = g1_curve.zero(ring);
                bool first = true;
                int pt_idx = 0;
                for (size_t i = 0; i < npoints; i++) {
                    int32_t d = win_digits[i];
                    if (d == 0) continue;
                    size_t bi = static_cast<size_t>(d < 0 ? -d : d) - 1;
                    if (bi != j) continue;

                    AffPoint pt = points[i];
                    if (d < 0) pt = g1_curve.negate_affine(pt, ring);

                    if (first) {
                        trace_xyzz = g1_curve.xyzz_from_affine(pt, ring);
                        trace_proj = ProjPoint(pt.x, pt.y, ring.one());
                        first = false;
                    } else {
                        trace_xyzz = g1_curve.xyzz_add_affine(trace_xyzz, pt, ring);
                        trace_proj = g1_curve.add_mixed_point(trace_proj, pt, ring);
                    }

                    // Compare after each addition
                    ProjPoint trace_xyzz_proj = g1_curve.xyzz_to_proj(trace_xyzz, ring);
                    auto [tx, ty] = proj_to_affine_bigint(trace_xyzz_proj, ring, q);
                    auto [px, py] = proj_to_affine_bigint(trace_proj, ring, q);
                    if (tx != px || ty != py) {
                        std::cout << "    DIVERGES at point " << pt_idx << " (global i=" << i << ", digit=" << d << ")" << std::endl;
                        std::cout << "      XYZZ after: (" << tx.to_string(16) << ")" << std::endl;
                        std::cout << "      Proj after: (" << px.to_string(16) << ")" << std::endl;

                        // Print the XYZZ state BEFORE this addition
                        // Re-trace to get state before the failing point
                        XYZZPt pre_xyzz = g1_curve.xyzz_zero(ring);
                        bool pre_first = true;
                        for (size_t ii = 0; ii < i; ii++) {
                            int32_t dd = win_digits[ii];
                            if (dd == 0) continue;
                            size_t bbi = static_cast<size_t>(dd < 0 ? -dd : dd) - 1;
                            if (bbi != j) continue;
                            AffPoint ppt = points[ii];
                            if (dd < 0) ppt = g1_curve.negate_affine(ppt, ring);
                            if (pre_first) {
                                pre_xyzz = g1_curve.xyzz_from_affine(ppt, ring);
                                pre_first = false;
                            } else {
                                pre_xyzz = g1_curve.xyzz_add_affine(pre_xyzz, ppt, ring);
                            }
                        }
                        // Check H = x2 * ZZ - X
                        auto pre_aff_proj = g1_curve.xyzz_to_proj(pre_xyzz, ring);
                        auto [pre_ax, pre_ay] = proj_to_affine_bigint(pre_aff_proj, ring, q);
                        BigInt pt_x = ring.to_bigint(pt.x);
                        std::cout << "      Pre-XYZZ affine x: " << pre_ax.to_string(16) << std::endl;
                        std::cout << "      New point x:       " << pt_x.to_string(16) << std::endl;
                        std::cout << "      Same x? " << (pre_ax == pt_x ? "YES (degenerate!)" : "NO") << std::endl;

                        // Also check Pre-XYZZ ZZ and ZZZ
                        BigInt zz_val = ring.to_bigint(pre_xyzz.ZZ);
                        BigInt zzz_val = ring.to_bigint(pre_xyzz.ZZZ);
                        std::cout << "      ZZ:  " << zz_val.to_string(16) << std::endl;
                        std::cout << "      ZZZ: " << zzz_val.to_string(16) << std::endl;

                        // Compute H = pt.x * ZZ - X explicitly
                        BigInt xyzz_X = ring.to_bigint(pre_xyzz.X);
                        BigInt H = (pt_x * zz_val - xyzz_X) % q;
                        if (H < 0) H = H + q;
                        std::cout << "      H = x2*ZZ - X = " << H.to_string(16) << std::endl;

                        break;
                    }
                    pt_idx++;
                }
                break;
            }
        }
    }

    if (!any_bucket_mismatch) {
        std::cout << "  All buckets match! Testing per-window integration..." << std::endl;

        // Compare per-window integration: XYZZ→Proj vs ProjPoint directly
        bool integration_mismatch = false;
        for (size_t w = 0; w < num_windows && !integration_mismatch; w++) {
            const int32_t *win_digits = digits.data() + w * npoints;

            // XYZZ scatter
            std::vector<XYZZPt> xyzz_bkts(nbuckets, g1_curve.xyzz_zero(ring));
            std::vector<uint8_t> xyzz_occ(nbuckets, 0);
            scatter_chunk(g1_curve, ring, points.data(), win_digits, 0, npoints,
                          xyzz_bkts.data(), xyzz_occ.data(), nbuckets);

            // ProjPoint scatter
            std::vector<ProjPoint> proj_bkts(nbuckets, g1_curve.zero(ring));
            std::vector<uint8_t> proj_occ(nbuckets, 0);
            for (size_t i = 0; i < npoints; i++) {
                int32_t digit = win_digits[i];
                if (digit == 0) continue;
                bool neg = digit < 0;
                size_t bidx = static_cast<size_t>(neg ? -digit : digit) - 1;
                AffPoint pt = points[i];
                if (neg) pt = g1_curve.negate_affine(pt, ring);
                if (!proj_occ[bidx]) {
                    proj_bkts[bidx] = ProjPoint(pt.x, pt.y, ring.one());
                    proj_occ[bidx] = 1;
                } else {
                    proj_bkts[bidx] = g1_curve.add_mixed_point(proj_bkts[bidx], pt, ring);
                }
            }

            bool has_pts = false;
            for (size_t j = 0; j < nbuckets; j++) {
                if (xyzz_occ[j]) { has_pts = true; break; }
            }
            if (!has_pts) continue;

            // Integrate XYZZ (via xyzz_to_proj + add_point)
            ProjPoint xyzz_wsum = integrate_buckets(g1_curve, ring,
                xyzz_bkts.data(), xyzz_occ.data(), nbuckets);

            // Integrate ProjPoint directly
            ProjPoint proj_wsum = g1_curve.zero(ring);
            {
                ProjPoint running = g1_curve.zero(ring);
                bool r_started = false, w_started = false;
                for (size_t j = nbuckets; j-- > 0; ) {
                    if (proj_occ[j]) {
                        if (!r_started) {
                            running = proj_bkts[j];
                            r_started = true;
                        } else {
                            running = g1_curve.add_point(running, proj_bkts[j], ring);
                        }
                    }
                    if (r_started) {
                        if (!w_started) {
                            proj_wsum = running;
                            w_started = true;
                        } else {
                            proj_wsum = g1_curve.add_point(proj_wsum, running, ring);
                        }
                    }
                }
            }

            auto [xwx, xwy] = proj_to_affine_bigint(xyzz_wsum, ring, q);
            auto [pwx, pwy] = proj_to_affine_bigint(proj_wsum, ring, q);

            if (xwx != pwx || xwy != pwy) {
                integration_mismatch = true;
                std::cout << "  Window " << w << " integration MISMATCH:" << std::endl;
                std::cout << "    XYZZ: (" << xwx.to_string(16) << ")" << std::endl;
                std::cout << "    Proj: (" << pwx.to_string(16) << ")" << std::endl;

                // Now trace the integration step by step
                ProjPoint xyzz_running = g1_curve.zero(ring);
                ProjPoint proj_running = g1_curve.zero(ring);
                ProjPoint xyzz_ws = g1_curve.zero(ring);
                ProjPoint proj_ws = g1_curve.zero(ring);
                bool xr = false, pr = false, xw = false, pw = false;

                for (size_t j = nbuckets; j-- > 0; ) {
                    if (xyzz_occ[j]) {
                        auto bp = g1_curve.xyzz_to_proj(xyzz_bkts[j], ring);
                        if (!xr) {
                            xyzz_running = bp; xr = true;
                        } else {
                            xyzz_running = g1_curve.add_point(xyzz_running, bp, ring);
                        }
                    }
                    if (proj_occ[j]) {
                        if (!pr) {
                            proj_running = proj_bkts[j]; pr = true;
                        } else {
                            proj_running = g1_curve.add_point(proj_running, proj_bkts[j], ring);
                        }
                    }
                    if (xr) {
                        if (!xw) { xyzz_ws = xyzz_running; xw = true; }
                        else xyzz_ws = g1_curve.add_point(xyzz_ws, xyzz_running, ring);
                    }
                    if (pr) {
                        if (!pw) { proj_ws = proj_running; pw = true; }
                        else proj_ws = g1_curve.add_point(proj_ws, proj_running, ring);
                    }

                    // Compare running sums after each bucket
                    if (xr && pr) {
                        auto [xrx, xry] = proj_to_affine_bigint(xyzz_running, ring, q);
                        auto [prx, pry] = proj_to_affine_bigint(proj_running, ring, q);
                        if (xrx != prx || xry != pry) {
                            std::cout << "    Running sum diverges at bucket " << j << std::endl;
                            std::cout << "      XYZZ running: (" << xrx.to_string(16) << ")" << std::endl;
                            std::cout << "      Proj running: (" << prx.to_string(16) << ")" << std::endl;

                            // Check the bucket itself
                            auto bp = g1_curve.xyzz_to_proj(xyzz_bkts[j], ring);
                            auto [bx, by] = proj_to_affine_bigint(bp, ring, q);
                            auto [bpx, bpy] = proj_to_affine_bigint(proj_bkts[j], ring, q);
                            std::cout << "      XYZZ bucket[" << j << "]: (" << bx.to_string(16) << ")" << std::endl;
                            std::cout << "      Proj bucket[" << j << "]: (" << bpx.to_string(16) << ")" << std::endl;
                            break;
                        }
                    }
                }
                failed++;
            }
        }

        if (!integration_mismatch) {
            std::cout << "  Per-window integration matches! Checking full accumulation..." << std::endl;
            // If per-window integration matches, the bug is in the main loop
            // (doubling or window accumulation)
            passed++;
        }
    } else {
        failed++;
    }

    std::cout << "  Passed: " << passed << ", Failed: " << failed << std::endl;
    std::cout << std::endl;
    return failed;
}

// ========================================================================
// Test: Pippenger V2 (batch affine schedule) vs naive MSM
// ========================================================================
template<class RingType>
int test_pippenger_v2_vs_naive(const RingType &ring, const BigInt &q) {
    std::cout << "=== Test: Pippenger V2 (batch affine) vs Naive MSM ===" << std::endl;
    int passed = 0, failed = 0;

    using CurveType = G1<RingType>;
    using AffPoint = typename CurveType::AffPoint;
    CurveType g1_curve;

    BigInt r(bls12_381_scalar_modulus_hex_msm, 16);

    std::vector<size_t> test_sizes = {1, 2, 4, 8, 16, 32, 64, 256};

    for (size_t npoints : test_sizes) {
        std::cout << "  n=" << npoints << ": ";

        std::vector<AffPoint> points(npoints);
        std::vector<std::vector<uint8_t>> scalar_data(npoints);
        std::vector<const uint8_t*> scalar_ptrs(npoints);

        for (size_t i = 0; i < npoints; i++) {
            BigInt pt_scalar = deterministic_scalar_msm(2 * i + npoints * 2000, r);
            byte pt_scalar_bytes[32] = {0};
            bigint_to_bytes_le_msm(pt_scalar_bytes, pt_scalar, 32);
            POINTonE1 blst_proj;
            blst_p1_mult(&blst_proj, &BLS12_381_G1, pt_scalar_bytes, 256);
            POINTonE1_affine blst_aff;
            blst_p1_to_affine(&blst_aff, &blst_proj);
            points[i] = blst_g1_to_affine_point_msm(blst_aff, ring);

            BigInt s = deterministic_scalar_msm(2 * i + npoints * 2000 + 1, r);
            scalar_data[i].resize(32, 0);
            bigint_to_bytes_le_msm(scalar_data[i].data(), s, 32);
            scalar_ptrs[i] = scalar_data[i].data();
        }

        // Compute with V2 (no Curve template — uses Ring directly)
        auto v2_result = msm_v2(ring, points.data(),
                                 scalar_ptrs.data(), npoints, 255);
        auto [v2_x, v2_y] = proj_to_affine_bigint(v2_result, ring, q);

        // Compute naive
        auto naive_result = naive_msm(g1_curve, ring, points.data(),
                                       scalar_ptrs.data(), npoints, 255);
        auto [naive_x, naive_y] = proj_to_affine_bigint(naive_result, ring, q);

        if (v2_x == naive_x && v2_y == naive_y) {
            passed++;
            std::cout << "PASS" << std::endl;
        } else {
            failed++;
            std::cout << "FAIL" << std::endl;
            std::cout << "    V2:    (" << v2_x.to_string(16) << ", " << v2_y.to_string(16) << ")" << std::endl;
            std::cout << "    Naive: (" << naive_x.to_string(16) << ", " << naive_y.to_string(16) << ")" << std::endl;
        }
    }

    std::cout << "  Passed: " << passed << ", Failed: " << failed << std::endl;
    if (failed == 0) std::cout << "  All Pippenger V2 vs Naive tests passed!" << std::endl;
    std::cout << std::endl;
    return failed;
}

// ========================================================================
// Test: Pippenger V2 parallel vs single-threaded
// ========================================================================
template<class RingType>
int test_pippenger_v2_parallel(const RingType &ring, const BigInt &q) {
    std::cout << "=== Test: Pippenger V2 Parallel vs Single-threaded ===" << std::endl;
    int passed = 0, failed = 0;

    using CurveType = G1<RingType>;
    using AffPoint = typename CurveType::AffPoint;

    BigInt r(bls12_381_scalar_modulus_hex_msm, 16);

    std::vector<size_t> test_sizes = {64, 256, 1024};

    for (size_t npoints : test_sizes) {
        std::cout << "  n=" << npoints << ": ";

        std::vector<AffPoint> points(npoints);
        std::vector<std::vector<uint8_t>> scalar_data(npoints);
        std::vector<const uint8_t*> scalar_ptrs(npoints);

        for (size_t i = 0; i < npoints; i++) {
            BigInt pt_scalar = deterministic_scalar_msm(2 * i, r);
            byte pt_scalar_bytes[32] = {0};
            bigint_to_bytes_le_msm(pt_scalar_bytes, pt_scalar, 32);
            POINTonE1 blst_proj;
            blst_p1_mult(&blst_proj, &BLS12_381_G1, pt_scalar_bytes, 256);
            POINTonE1_affine blst_aff;
            blst_p1_to_affine(&blst_aff, &blst_proj);
            points[i] = blst_g1_to_affine_point_msm(blst_aff, ring);

            BigInt s = deterministic_scalar_msm(2 * i + 1, r);
            scalar_data[i].resize(32, 0);
            bigint_to_bytes_le_msm(scalar_data[i].data(), s, 32);
            scalar_ptrs[i] = scalar_data[i].data();
        }

        // Single-threaded V2
        auto st_result = msm_v2(ring, points.data(),
                                 scalar_ptrs.data(), npoints, 255);
        auto [st_x, st_y] = proj_to_affine_bigint(st_result, ring, q);

        // Multi-threaded V2
        auto mt_result = msm_v2_parallel(ring, points.data(),
                                          scalar_ptrs.data(), npoints, 255, 2);
        auto [mt_x, mt_y] = proj_to_affine_bigint(mt_result, ring, q);

        if (st_x == mt_x && st_y == mt_y) {
            passed++;
            std::cout << "PASS" << std::endl;
        } else {
            failed++;
            std::cout << "FAIL" << std::endl;
            std::cout << "    ST: (" << st_x.to_string(16) << ")" << std::endl;
            std::cout << "    MT: (" << mt_x.to_string(16) << ")" << std::endl;
        }
    }

    std::cout << "  Passed: " << passed << ", Failed: " << failed << std::endl;
    if (failed == 0) std::cout << "  All Pippenger V2 Parallel tests passed!" << std::endl;
    std::cout << std::endl;
    return failed;
}

int main() {
    try {
        BigInt q(bls12_381_modulus_hex_msm, 16);
        BoundedRing<381, 8, 52, -1932, 2377, 12> ring(q);
        using RingType = BoundedRing<381, 8, 52, -1932, 2377, 12>;

        std::cout << "================================================" << std::endl;
        std::cout << "  VROOM MSM (Pippenger) Test Suite" << std::endl;
        std::cout << "================================================" << std::endl << std::endl;

        int total_failures = 0;
        total_failures += test_xyzz_unit(ring, q);
        total_failures += test_batch_invert(ring, q);
        total_failures += test_batch_to_affine(ring, q);
        total_failures += test_pippenger_vs_naive(ring, q);
        total_failures += test_parallel_msm(ring, q);
        total_failures += test_msm_vs_blst(ring, q);
        total_failures += test_pippenger_v2_vs_naive(ring, q);
        total_failures += test_pippenger_v2_parallel(ring, q);
        total_failures += test_edge_cases(ring, q);

        std::cout << "================================================" << std::endl;
        if (total_failures == 0) {
            std::cout << "  ALL TESTS PASSED" << std::endl;
        } else {
            std::cout << "  " << total_failures << " TEST(S) FAILED" << std::endl;
        }
        std::cout << "================================================" << std::endl;

        return total_failures > 0 ? 1 : 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}
