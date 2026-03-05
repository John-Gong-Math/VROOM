# G1 Scalar Multiplication Benchmark Results

## Environment

- **Platform**: Google Cloud `c3-standard-4` (Intel Sapphire Rapids, AVX512-IFMA)
- **Region**: us-central1-a
- **Compiler**: clang-14, `-O3 -march=native -mavx512ifma`
- **Configuration**: 50-bit RNS Montgomery (BLS12-381)
- **Date**: 2026-03-05

## Results

```
---------------------------------------------------------------------
Benchmark                           Time             CPU   Iterations
---------------------------------------------------------------------
BM_G1_Add                         224 ns          224 ns      3116729
BM_G1_MixedAdd                    216 ns          216 ns      3237521
BM_G1_Double                      203 ns          203 ns      3448418
BM_G1_ScalarMult_GLV_50bit      44303 ns        44302 ns        15786
```

## Column Explanations

- **Time**: Wall-clock time per call (average). Includes everything: CPU work, cache misses, context switches.
- **CPU**: CPU time per call (average). Excludes time when the process was not scheduled on a CPU. In practice, for short compute-bound benchmarks like these, Time and CPU are nearly identical.
- **Iterations**: The number of times Google Benchmark invoked the function. This is NOT a user-controlled parameter. Google Benchmark automatically increases iterations until the total runtime is long enough (~0.5s) to produce a statistically stable average. Higher iterations = the function is faster. For example:
  - `BM_G1_Double` at 203 ns/call ran 3,448,418 times (~0.7s total)
  - `BM_G1_ScalarMult_GLV_50bit` at 44,303 ns/call ran 15,786 times (~0.7s total)

## Benchmark Descriptions

- **BM_G1_Add**: Projective point addition on G1 (P + Q where P != Q)
- **BM_G1_MixedAdd**: Mixed addition on G1 (projective P + affine Q, saves one multiplication since Q.z = 1)
- **BM_G1_Double**: Point doubling on G1 (2P)
- **BM_G1_ScalarMult_GLV_50bit**: Full scalar multiplication on G1 using the GLV endomorphism. Decomposes a 256-bit scalar into two 128-bit scalars via the BLS12-381 endomorphism, then runs a 2-table windowed multi-scalar multiplication (window size = 4 bits). This is the end-to-end cost of computing `[k]P` for a random 256-bit scalar `k`.

## Comparison with AWS (from README)

| Benchmark | GCP c3-standard-4 (clang-14) | AWS c7i.metal-24xl (clang-21) |
|---|---|---|
| G1_Add | 224 ns | 174 ns |
| G1_MixedAdd | 216 ns | 168 ns |
| G1_Double | 203 ns | 169 ns |
| G1_ScalarMult_GLV | 44,303 ns | 35,893 ns |

The AWS numbers are faster due to: (1) bare-metal instance (no virtualization overhead), (2) newer clang-21 compiler generating better code, (3) higher base clock speed.

## Running with Statistics

To get mean/median/stddev across multiple runs, use:

```bash
./bench_pairing_50bit --benchmark_filter='BM_G1_' --benchmark_repetitions=10
```

This reports per-run timings plus summary lines:
- `_mean` -- average across all repetitions
- `_median` -- middle value
- `_stddev` -- standard deviation
- `_cv` -- coefficient of variation (stddev / mean), lower = more stable
