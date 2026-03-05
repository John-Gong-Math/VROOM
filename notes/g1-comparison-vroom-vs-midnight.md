# G1 Benchmark Comparison: VROOM vs midnight-zk

## VM Setup

- **Cloud**: Google Cloud Platform
- **Instance type**: `c3-standard-4`
- **Zone**: us-central1-a
- **CPU**: Intel Sapphire Rapids, 4 vCPUs @ 2700 MHz
- **CPU Caches**: L1d 48 KiB (x2), L1i 32 KiB (x2), L2 2048 KiB (x2), L3 107520 KiB (x1)
- **ISA extensions**: AVX512-IFMA confirmed (`/proc/cpuinfo`)
- **OS**: Debian 12 (debian-cloud/debian-12)
- **Date**: 2026-03-05

### VROOM (C++)
- **Compiler**: clang-14 (Debian package), flags: `-std=c++20 -O3 -march=native -mavx512ifma -D__ADX__`
- **Configuration**: 50-bit RNS Montgomery (`BoundedRing<381, 8, 50, -1932, 2377, 12>`), GLV endomorphism
- **Benchmark framework**: Google Benchmark (`libbenchmark-dev`)

### midnight-zk (Rust)
- **Toolchain**: rustc 1.90.0 (1159e78c4 2025-09-14), x86_64-unknown-linux-gnu
- **Backend**: BLST 0.3.15 (`blst_p1_mult` — C/x86-64 assembly)
- **Build profile**: `opt-level=3`, `lto=true`, `codegen-units=1`, `overflow-checks=false`
- **Benchmark framework**: Criterion.rs 0.7 (1000 samples, 3s warm-up)

## Results

| Operation | VROOM (C++) | midnight-zk (Rust/BLST) | Speedup |
|---|---|---|---|
| G1 Addition | 225 ns | 797 ns | **3.5x** |
| G1 Mixed Addition | 216 ns | 616 ns | **2.9x** |
| G1 Doubling | 203 ns | 349 ns | **1.7x** |
| G1 Scalar Multiplication | 45,000 ns | 97,454 ns | **2.2x** |

## Raw Output

### VROOM
```
BM_G1_Add                         225 ns          225 ns      3115479
BM_G1_MixedAdd                    216 ns          216 ns      3238115
BM_G1_Double                      203 ns          203 ns      3443975
BM_G1_ScalarMult_GLV_50bit      45001 ns        45000 ns        15768
```

### midnight-zk
```
Bn256-G1 addition          time: [797.29 ns 797.43 ns 797.59 ns]
Bn256-G1 mixed addition    time: [615.80 ns 615.94 ns 616.08 ns]
Bn256-G1 doubling          time: [349.20 ns 349.32 ns 349.46 ns]
Bn256-G1 scalar mult       time: [97.435 µs 97.454 µs 97.476 µs]
```

## Analysis

VROOM is faster across all operations. The speedup varies by operation:

- **G1 Addition (3.5x)**: Largest gap. VROOM's 8-wide AVX512-IFMA vectorization parallelizes the field multiplications within a single point addition far more effectively than BLST's scalar Montgomery approach.

- **G1 Mixed Addition (2.9x)**: Similar advantage. Mixed addition saves one field multiplication (since the affine point has z=1), benefiting both implementations proportionally.

- **G1 Doubling (1.7x)**: Smallest gap. Doubling has fewer field multiplications than addition, so there is less opportunity for AVX512 parallelism to shine.

- **G1 Scalar Multiplication (2.2x)**: VROOM uses the GLV endomorphism to decompose a 256-bit scalar into two 128-bit half-scalars, halving the number of point doublings. Combined with AVX512-IFMA field arithmetic, this gives a 2.2x advantage over BLST's `blst_p1_mult` (which uses optimized x86-64 assembly but without GLV or AVX512-IFMA for field ops).

## Implementation Differences

| Aspect | VROOM | midnight-zk |
|---|---|---|
| Language | C++20 | Rust (wrapping C/asm) |
| Field arithmetic | RNS Montgomery, 50-bit limbs, AVX512-IFMA | Montgomery, 64-bit limbs, x86-64 assembly |
| Scalar mult algorithm | GLV (2 × 128-bit, window=4) | Standard (1 × 255-bit) |
| Vectorization | 8-wide AVX512 (`__m512i`) | Scalar (no SIMD for field ops) |

## Reproducing

```bash
# From the VROOM repo root:
./bench_g1_comparison_gcp.sh
```

Requires: gcloud CLI authenticated with billing-enabled project.
