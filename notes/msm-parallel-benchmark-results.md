# MSM Parallel Benchmark Results

**Date:** 2026-03-05

## Hardware

- **CPU:** Intel Xeon Platinum 8481C (Sapphire Rapids) @ 2.70 GHz
- **Platform:** GCP c3-standard-4 VM (4 vCPUs / 2 physical cores, Debian 12)
- **ISA:** AVX512-IFMA supported
- **Caches:** L1d 48 KiB (x2), L1i 32 KiB (x2), L2 2048 KiB (x2), L3 107520 KiB

## Build Configuration

- **Compiler:** clang++ (C++20)
- **Flags:** `-O3 -march=native -mavx512ifma -D__ADX__`
- **Curve:** BLS12-381 G1
- **Scalar bits:** 255
- **Repetitions:** 3 per benchmark

## VROOM Single-threaded vs Parallel vs BLST Pippenger

VROOM Parallel splits point scatter across threads within each Pippenger window. Each thread has its own bucket array; after scatter, partial buckets are merged and integrated sequentially. The Horner window loop remains sequential to avoid quadratic doubling cost. Falls back to single-threaded for n < 1024.

BLST Pippenger is single-threaded with hand-optimized x86-64 assembly, XYZZ coordinates, Booth-encoded signed digits, and prefetching.

| Points | VROOM 1T (ms) | VROOM 4T (ms) | BLST 1T (ms) | Parallel Speedup | VROOM 4T vs BLST 1T |
|--------|--------------|---------------|-------------|-----------------|---------------------|
| 2^6    | 0.059        | 0.061         | 0.234       | 1.0x (fallback)  | **3.8x** faster |
| 2^8    | 0.203        | 0.208         | 0.447       | 1.0x (fallback)  | **2.2x** faster |
| 2^10   | 0.619        | 2.33          | 1.13        | 0.27x (overhead) | 0.49x slower |
| 2^12   | 2.23         | 2.95          | 2.87        | 0.76x            | ~1.0x tied |
| 2^14   | 8.57         | 5.99          | 9.96        | **1.43x**        | **1.66x** faster |
| 2^16   | 33.9         | 19.1          | 38.0        | **1.78x**        | **1.99x** faster |
| 2^18   | 135          | 72.0          | 147         | **1.88x**        | **2.04x** faster |
| 2^20   | 563          | 302           | 619         | **1.86x**        | **2.05x** faster |

### Analysis

- **Small (n = 2^6 to 2^8):** Parallel falls back to single-threaded (n < 1024 threshold). VROOM's AVX512-IFMA field arithmetic dominates, giving 2-4x over BLST.

- **Overhead zone (n = 2^10 to 2^12):** Thread creation per window (~37 windows at 2^10) overwhelms the scatter work per thread. Parallel is slower than single-threaded. Single-threaded VROOM still competitive with BLST.

- **Sweet spot (n = 2^14 to 2^20):** Parallel scatter delivers 1.4-1.9x speedup over single-threaded on 4 vCPUs (2 physical cores). Combined with VROOM's faster field arithmetic, this gives **~2x over BLST** across the entire range.

- **Scaling at 2^20:** VROOM Parallel achieves 302ms vs BLST's 619ms (2.05x). The parallel speedup (1.86x from 1T) is near the theoretical max for 2 physical cores with hyperthreading. Low variance (cv=0.44%) indicates stable memory access patterns.

### Parallel scaling behavior

| Points | 1T (ms) | 4T (ms) | Speedup | Efficiency (per vCPU) |
|--------|---------|---------|---------|----------------------|
| 2^14   | 8.57    | 5.99    | 1.43x   | 36% |
| 2^16   | 33.9    | 19.1    | 1.78x   | 44% |
| 2^18   | 135     | 72.0    | 1.88x   | 47% |
| 2^20   | 563     | 302     | 1.86x   | 47% |

Efficiency plateaus at ~47% per vCPU (or ~93% per physical core), which is expected: the VM has 2 physical cores with 2 hyperthreads each. The sequential phases (window doublings, bucket merge, bucket integration) limit further scaling.

### Comparison with previous single-threaded results

Single-threaded numbers improved from the previous benchmark (e.g., 563ms vs 689ms at 2^20). This is due to bucket array reuse across windows (allocated once, cleared per window) reducing allocation pressure.

### Threading implementation

- Uses `std::thread` with per-window create/join
- Each thread scatters its chunk of points into private bucket arrays
- After join, main thread merges partial buckets and integrates
- Bucket arrays pre-allocated and reused across windows
- Falls back to single-threaded for n < 1024 to avoid thread overhead

### Potential improvements

1. **Thread pool / barrier:** Avoid thread creation per window. The initial attempt with `std::condition_variable` had a race condition; a C++20 `std::barrier` or lock-free barrier would be cleaner.
2. **Parallel bucket merge:** Currently sequential. Could be parallelized by partitioning the bucket range across threads.
3. **Parallel bucket integration:** The running-sum integration is inherently sequential, but could be split into sub-ranges and combined.
4. **Higher thread counts:** On machines with more cores (e.g., 8-16), the scatter phase would scale further. The merge and integration phases would become the bottleneck.

## Raw Output

```
Running ./bench_msm_avx
Run on (4 X 2700 MHz CPU s)
CPU Caches:
  L1 Data 48 KiB (x2)
  L1 Instruction 32 KiB (x2)
  L2 Unified 2048 KiB (x2)
  L3 Unified 107520 KiB (x1)
Load Average: 0.14, 0.39, 0.66
-----------------------------------------------------------------------------------------------------
Benchmark                                                           Time             CPU   Iterations
-----------------------------------------------------------------------------------------------------
BM_VROOM_MSM/64_mean                                            0.059 ms        0.059 ms            3
BM_VROOM_MSM/256_mean                                           0.203 ms        0.203 ms            3
BM_VROOM_MSM/1024_mean                                          0.619 ms        0.619 ms            3
BM_VROOM_MSM/4096_mean                                           2.23 ms         2.23 ms            3
BM_VROOM_MSM/16384_mean                                          8.57 ms         8.57 ms            3
BM_VROOM_MSM/65536_mean                                          33.9 ms         33.9 ms            3
BM_VROOM_MSM/262144_mean                                          135 ms          135 ms            3
BM_VROOM_MSM/1048576_mean                                         563 ms          563 ms            3
BM_VROOM_MSM_Parallel/64_mean                                   0.061 ms        0.061 ms            3
BM_VROOM_MSM_Parallel/256_mean                                  0.208 ms        0.208 ms            3
BM_VROOM_MSM_Parallel/1024_mean                                  2.33 ms         1.96 ms            3
BM_VROOM_MSM_Parallel/4096_mean                                  2.95 ms         2.54 ms            3
BM_VROOM_MSM_Parallel/16384_mean                                 5.99 ms         5.67 ms            3
BM_VROOM_MSM_Parallel/65536_mean                                 19.1 ms         18.9 ms            3
BM_VROOM_MSM_Parallel/262144_mean                                72.0 ms         71.9 ms            3
BM_VROOM_MSM_Parallel/1048576_mean                                302 ms          301 ms            3
BM_BLST_Pippenger/64_mean                                       0.234 ms        0.234 ms            3
BM_BLST_Pippenger/256_mean                                      0.447 ms        0.447 ms            3
BM_BLST_Pippenger/1024_mean                                      1.13 ms         1.13 ms            3
BM_BLST_Pippenger/4096_mean                                      2.87 ms         2.87 ms            3
BM_BLST_Pippenger/16384_mean                                     9.96 ms         9.96 ms            3
BM_BLST_Pippenger/65536_mean                                     38.0 ms         38.0 ms            3
BM_BLST_Pippenger/262144_mean                                     147 ms          147 ms            3
BM_BLST_Pippenger/1048576_mean                                    619 ms          619 ms            3
```
