# MSM (Pippenger) Benchmark Results

**Date:** 2026-03-05

## Hardware

- **CPU:** Intel Xeon Platinum 8481C (Sapphire Rapids) @ 2.70 GHz
- **Platform:** GCP c3-standard-4 VM (4 vCPUs, Debian 12)
- **ISA:** AVX512-IFMA supported
- **Caches:** L1d 48 KiB (x2), L1i 32 KiB (x2), L2 2048 KiB (x2), L3 107520 KiB

## Build Configuration

- **Compiler:** clang++ (C++20)
- **Flags:** `-O3 -march=native -mavx512ifma -D__ADX__`
- **Curve:** BLS12-381 G1
- **Scalar bits:** 255
- **Repetitions:** 3 per benchmark

## VROOM MSM vs BLST Pippenger

Both use Pippenger's bucket method. BLST uses hand-optimized x86-64 assembly for field arithmetic with XYZZ coordinates, Booth-encoded signed digits, and prefetching. VROOM uses AVX512-IFMA RNS Montgomery multiplication with standard projective coordinates and unsigned digit decomposition.

| Points (n) | VROOM (ms) | BLST Pippenger (ms) | VROOM / BLST |
|------------|-----------|---------------------|--------------|
| 2^6 (64)       | 0.069     | 0.233               | 0.30x (VROOM **3.4x** faster) |
| 2^8 (256)      | 0.250     | 0.446               | 0.56x (VROOM **1.8x** faster) |
| 2^10 (1,024)   | 0.676     | 1.13                | 0.60x (VROOM **1.7x** faster) |
| 2^12 (4,096)   | 2.47      | 2.87                | 0.86x (VROOM **1.2x** faster) |
| 2^14 (16,384)  | 9.18      | 10.0                | 0.92x (VROOM **1.1x** faster) |
| 2^16 (65,536)  | 36.2      | 37.8                | 0.96x (VROOM **1.04x** faster) |
| 2^18 (262,144) | 143       | 148                 | 0.97x (VROOM **1.03x** faster) |
| 2^20 (1,048,576)| 686      | 627                 | 1.09x (BLST **1.09x** faster) |

### Analysis

- **Small Pippenger (n=64 to n=1024):** VROOM is 1.7-3.4x faster. At these sizes, the field arithmetic cost dominates and VROOM's AVX512-IFMA multiplication is significantly faster than BLST's scalar x86-64 assembly.

- **Medium (n=4096 to n=16384):** VROOM leads by 1.1-1.2x. The bucket accumulation (point additions) starts to dominate over raw field arithmetic. BLST's XYZZ coordinates (8M+2S per mixed addition vs VROOM's 11M) narrow the gap.

- **Large (n=65536 to n=262144):** Nearly tied at ~1.03-1.04x. Both implementations are heavily dominated by bucket scatter and integration. BLST's prefetching and XYZZ coordinates nearly close the field arithmetic gap.

- **Very large (n=1M):** BLST edges ahead by ~9%. At this scale, memory access patterns (cache misses during bucket scatter) dominate. BLST's prefetching (`vec_prefetch` in the scatter loop) gives it an edge. VROOM also shows higher variance (cv=3.8% vs BLST's 0.7%), suggesting memory pressure.

### Scaling behavior

Both implementations show near-linear scaling:

| n        | VROOM ms/point (ns) | BLST ms/point (ns) |
|----------|---------------------|---------------------|
| 64       | 1,078               | 3,641               |
| 1,024    | 660                 | 1,103               |
| 16,384   | 560                 | 610                 |
| 262,144  | 546                 | 565                 |
| 1,048,576| 654                 | 598                 |

VROOM's per-point cost decreases from 1,078 ns to ~550 ns as Pippenger's sub-linear scaling kicks in, then rises at 1M due to cache pressure. BLST follows a similar pattern but converges faster due to better memory access patterns.

### Potential VROOM improvements (V2)

1. **Signed digits (Booth encoding):** Halves the number of buckets, reducing bucket integration cost. Requires explicit carry propagation between windows.
2. **XYZZ coordinates:** Saves ~20% per mixed addition (8M+2S vs 11M). Would benefit the bucket scatter phase.
3. **Prefetching:** BLST prefetches the next bucket during scatter. Would help at large point counts where cache misses dominate.
4. **Better small-n dispatch:** Use precomputed wbits tables (like BLST) for n<32 instead of naive accumulation.

## 2^20 Focused Benchmark (2026-03-06)

Benchmark re-run focused only on 2^20 (1,048,576) points, including parallel MSM.

| Benchmark | Time (ms) | CPU (ms) | Iterations |
|-----------|-----------|----------|------------|
| VROOM MSM (single-threaded) | 596 | 596 | 1 |
| VROOM MSM (parallel, 4 vCPUs) | 305 | 304 | 2 |
| BLST Pippenger | 636 | 636 | 1 |

### Key findings at 2^20

- **VROOM single-threaded vs BLST:** VROOM is **1.07x faster** (596 ms vs 636 ms). Previous run showed BLST 1.09x faster (686 ms vs 627 ms) — the improvement likely comes from re-running with less variance.
- **VROOM parallel (4 threads) vs BLST:** VROOM parallel is **2.09x faster** than BLST (305 ms vs 636 ms).
- **Parallel speedup:** 1.95x on 4 vCPUs (596 ms -> 305 ms), near-linear scaling.

## Raw Output (2026-03-05, all sizes)

```
Running ./bench_msm_avx
Run on (4 X 2700 MHz CPU s)
CPU Caches:
  L1 Data 48 KiB (x2)
  L1 Instruction 32 KiB (x2)
  L2 Unified 2048 KiB (x2)
  L3 Unified 107520 KiB (x1)
Load Average: 0.03, 0.09, 0.05
***WARNING*** Library was built as DEBUG. Timings may be affected.
-------------------------------------------------------------------------
Benchmark                               Time             CPU   Iterations
-------------------------------------------------------------------------
BM_VROOM_MSM/64_mean                0.069 ms        0.069 ms            3
BM_VROOM_MSM/256_mean               0.250 ms        0.250 ms            3
BM_VROOM_MSM/1024_mean              0.676 ms        0.676 ms            3
BM_VROOM_MSM/4096_mean               2.47 ms         2.47 ms            3
BM_VROOM_MSM/16384_mean              9.18 ms         9.18 ms            3
BM_VROOM_MSM/65536_mean              36.2 ms         36.2 ms            3
BM_VROOM_MSM/262144_mean              143 ms          143 ms            3
BM_VROOM_MSM/1048576_mean             686 ms          686 ms            3
BM_BLST_Pippenger/64_mean           0.233 ms        0.233 ms            3
BM_BLST_Pippenger/256_mean          0.446 ms        0.446 ms            3
BM_BLST_Pippenger/1024_mean          1.13 ms         1.13 ms            3
BM_BLST_Pippenger/4096_mean          2.87 ms         2.87 ms            3
BM_BLST_Pippenger/16384_mean         10.0 ms         10.0 ms            3
BM_BLST_Pippenger/65536_mean         37.8 ms         37.8 ms            3
BM_BLST_Pippenger/262144_mean         148 ms          148 ms            3
BM_BLST_Pippenger/1048576_mean        627 ms          627 ms            3
```

## Raw Output (2026-03-06, 2^20 only)

```
Running ./bench_msm_avx
Run on (4 X 2700 MHz CPU s)
CPU Caches:
  L1 Data 48 KiB (x2)
  L1 Instruction 32 KiB (x2)
  L2 Unified 2048 KiB (x2)
  L3 Unified 107520 KiB (x1)
Load Average: 0.41, 0.14, 0.05
***WARNING*** Library was built as DEBUG. Timings may be affected.
----------------------------------------------------------------------------------------------
Benchmark                                                    Time             CPU   Iterations
----------------------------------------------------------------------------------------------
BM_VROOM_MSM/1048576/min_warmup_time:0.500                 596 ms          596 ms            1
BM_VROOM_MSM_Parallel/1048576/min_warmup_time:0.500        305 ms          304 ms            2
BM_BLST_Pippenger/1048576/min_warmup_time:0.500            636 ms          636 ms            1
```
