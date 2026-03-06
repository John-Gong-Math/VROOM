#include <benchmark/benchmark.h>
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
#include "bounded_ring.hpp"
#include "conversion_inversion.hpp"
#include <vector>
#include <cstring>

// BLS12-381 modulus q
const char* bench_bls12_381_modulus_hex = "1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab";
// BLS12-381 scalar field modulus r
const char* bench_bls12_381_scalar_hex = "73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001";

// BLST external functions
extern "C" {
    extern const POINTonE1 BLS12_381_G1;
    void blst_p1_mult(POINTonE1 *out, const POINTonE1 *a, const byte *scalar, size_t nbits);
    void blst_p1_to_affine(POINTonE1_affine *out, const POINTonE1 *a);
    void blst_p1_add(POINTonE1 *out, const POINTonE1 *a, const POINTonE1 *b);

    // BLST Pippenger MSM
    size_t blst_p1s_mult_pippenger_scratch_sizeof(size_t npoints);
    void blst_p1s_mult_pippenger(POINTonE1 *ret,
                                  const POINTonE1_affine *const points[],
                                  size_t npoints,
                                  const byte *const scalars[], size_t nbits,
                                  void *scratch);
}

static void bigint_to_bytes_le_bench(uint8_t *bytes, const BigInt &value, size_t len) {
    BigInt temp = value;
    for (size_t i = 0; i < len; i++) {
        bytes[i] = static_cast<uint8_t>(temp.to_ulong() & 0xff);
        temp = temp >> 8;
    }
}

static BigInt vec384_montgomery_to_bigint_bench(const vec384 a) {
    vec384 normal;
    from_fp(normal, a);
    BigInt result(0);
    BigInt two_to_64 = BigInt(1) << 64;
    for (int i = 5; i >= 0; i--) {
        result = result * two_to_64 + BigInt(static_cast<unsigned long>(normal[i]));
    }
    return result;
}

using RingType = BoundedRing<381, 8, 52, -1932, 2377, 12>;

// Shared test data: BLST affine points + scalars + VROOM affine points
struct MSMTestData {
    // BLST format
    std::vector<POINTonE1_affine> blst_points;
    std::vector<const POINTonE1_affine*> blst_point_ptrs;
    // VROOM format
    std::vector<AffinePoint<typename RingType::StandardElement>> vroom_points;
    // Shared scalars
    std::vector<std::vector<uint8_t>> scalar_data;
    std::vector<const uint8_t*> scalar_ptrs;
};

static MSMTestData generate_test_data(const RingType &ring, size_t npoints) {
    MSMTestData data;
    BigInt r(bench_bls12_381_scalar_hex, 16);

    data.blst_points.resize(npoints);
    data.blst_point_ptrs.resize(npoints);
    data.vroom_points.resize(npoints);
    data.scalar_data.resize(npoints);
    data.scalar_ptrs.resize(npoints);

    for (size_t i = 0; i < npoints; i++) {
        // Generate point: random_scalar * G
        BigInt pt_scalar = BigInt::random(256) % r;
        if (pt_scalar == BigInt(0)) pt_scalar = BigInt(1);
        byte pt_scalar_bytes[32] = {0};
        bigint_to_bytes_le_bench(pt_scalar_bytes, pt_scalar, 32);
        POINTonE1 blst_proj;
        blst_p1_mult(&blst_proj, &BLS12_381_G1, pt_scalar_bytes, 256);
        blst_p1_to_affine(&data.blst_points[i], &blst_proj);
        data.blst_point_ptrs[i] = &data.blst_points[i];

        // Convert to VROOM format
        BigInt x = vec384_montgomery_to_bigint_bench(data.blst_points[i].X);
        BigInt y = vec384_montgomery_to_bigint_bench(data.blst_points[i].Y);
        data.vroom_points[i] = AffinePoint<typename RingType::StandardElement>(
            ring.from_bigint(x), ring.from_bigint(y));

        // Generate random scalar
        BigInt s = BigInt::random(255) % r;
        data.scalar_data[i].resize(32, 0);
        bigint_to_bytes_le_bench(data.scalar_data[i].data(), s, 32);
        data.scalar_ptrs[i] = data.scalar_data[i].data();
    }
    return data;
}

// ---- VROOM MSM Benchmark ----

static void BM_VROOM_MSM(benchmark::State& state) {
    size_t npoints = static_cast<size_t>(state.range(0));
    BigInt q(bench_bls12_381_modulus_hex, 16);
    RingType ring(q);
    G1<RingType> g1_curve;

    auto data = generate_test_data(ring, npoints);
    benchmark::DoNotOptimize(data);

    for (auto _ : state) {
        auto result = msm(g1_curve, ring, data.vroom_points.data(),
                          data.scalar_ptrs.data(), npoints, 255);
        benchmark::DoNotOptimize(result);
    }
}

// ---- VROOM Parallel MSM Benchmark ----

static void BM_VROOM_MSM_Parallel(benchmark::State& state) {
    size_t npoints = static_cast<size_t>(state.range(0));
    BigInt q(bench_bls12_381_modulus_hex, 16);
    RingType ring(q);
    G1<RingType> g1_curve;

    auto data = generate_test_data(ring, npoints);
    benchmark::DoNotOptimize(data);

    for (auto _ : state) {
        auto result = msm_parallel(g1_curve, ring, data.vroom_points.data(),
                                    data.scalar_ptrs.data(), npoints, 255, 0);
        benchmark::DoNotOptimize(result);
    }
}

// ---- BLST Pippenger MSM Benchmark ----

static void BM_BLST_Pippenger(benchmark::State& state) {
    size_t npoints = static_cast<size_t>(state.range(0));
    BigInt q(bench_bls12_381_modulus_hex, 16);
    RingType ring(q);

    auto data = generate_test_data(ring, npoints);

    // Allocate scratch buffer for BLST Pippenger
    size_t scratch_sz = blst_p1s_mult_pippenger_scratch_sizeof(npoints);
    std::vector<uint8_t> scratch(scratch_sz);
    benchmark::DoNotOptimize(data);

    for (auto _ : state) {
        POINTonE1 result;
        blst_p1s_mult_pippenger(&result,
                                 data.blst_point_ptrs.data(),
                                 npoints,
                                 data.scalar_ptrs.data(), 255,
                                 scratch.data());
        benchmark::DoNotOptimize(result);
    }
}

// Register benchmarks for 2^20 points
BENCHMARK(BM_VROOM_MSM)->Arg(1048576)->Unit(benchmark::kMillisecond)->MinWarmUpTime(0.5);
BENCHMARK(BM_VROOM_MSM_Parallel)->Arg(1048576)->Unit(benchmark::kMillisecond)->MinWarmUpTime(0.5);
BENCHMARK(BM_BLST_Pippenger)->Arg(1048576)->Unit(benchmark::kMillisecond)->MinWarmUpTime(0.5);

BENCHMARK_MAIN();
