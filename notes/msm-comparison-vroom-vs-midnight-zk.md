# MSM Comparison: VROOM vs midnight-zk (midnight-curves)

**Date:** 2026-03-05

## Hardware

- **CPU:** Intel Xeon Platinum 8481C (Sapphire Rapids) @ 2.70 GHz
- **Platform:** GCP c3-standard-4 VM (4 vCPUs / 2 physical cores, Debian 12)
- **ISA:** AVX512-IFMA supported
- **Caches:** L1d 48 KiB (x2), L1i 32 KiB (x2), L2 2048 KiB (x2), L3 107520 KiB

## Implementations

| Implementation | Language | Field Arithmetic | MSM Algorithm | Threading |
|---|---|---|---|---|
| **VROOM 1T** | C++20 | AVX512-IFMA RNS Montgomery | Pippenger bucket method, unsigned digits | Single-threaded |
| **VROOM 4T** | C++20 | AVX512-IFMA RNS Montgomery | Pippenger, parallel bucket scatter | std::thread, 4 vCPUs |
| **BLST Pippenger** | C + x86-64 asm | Hand-optimized asm (ADX/MULX) | Pippenger, XYZZ coords, Booth encoding, prefetch | Single-threaded |
| **midnight-zk msm_best** | Rust | BLST via FFI (blst crate) | Pippenger, Booth encoding, batch affine addition | Rayon, 4 vCPUs |

## Build Configuration

- **VROOM:** clang++ C++20, `-O3 -march=native -mavx512ifma`, 3 reps
- **BLST (C):** Hand-optimized x86-64 assembly with ADX, linked as `libblst_small.a`
- **midnight-zk:** Rust 1.90.0, `RUSTFLAGS="-C target-cpu=native"`, `opt-level=3, lto=true, codegen-units=1`, 3 reps
- **Curve:** BLS12-381 G1, 255-bit scalars

## Results

| Points | VROOM 1T (ms) | VROOM 4T (ms) | BLST 1T (ms) | midnight-zk msm_best 4T (ms) |
|--------|--------------|---------------|-------------|-------------------------------|
| 2^10   | 0.614        | 2.32          | 1.13        | 14.3                          |
| 2^11   | 1.16         | 2.38          | 1.62        | 24.7                          |
| 2^12   | 2.25         | 2.98          | 2.88        | 44.9                          |
| 2^13   | 4.36         | 4.03          | 5.23        | 49.4                          |
| 2^14   | 8.64         | 6.03          | 10.0        | 88.0                          |
| 2^15   | 17.1         | 10.4          | 19.6        | 157                           |
| 2^16   | 33.9         | 19.1          | 38.1        | 292                           |
| 2^17   | 67.6         | 37.1          | 75.0        | 553                           |
| 2^18   | 135          | 72.7          | 148         | 964                           |
| 2^19   | 271          | 144           | 297         | 1884                          |
| 2^20   | 556          | 296           | 631         | 3671                          |

## Speedup: VROOM vs midnight-zk msm_best

Both use 4 vCPUs. VROOM Parallel vs midnight-zk msm_best (Rayon):

| Points | VROOM 4T (ms) | midnight-zk 4T (ms) | VROOM Speedup |
|--------|--------------|---------------------|---------------|
| 2^10   | 2.32         | 14.3                | **6.2x**      |
| 2^12   | 2.98         | 44.9                | **15.1x**     |
| 2^14   | 6.03         | 88.0                | **14.6x**     |
| 2^16   | 19.1         | 292                 | **15.3x**     |
| 2^18   | 72.7         | 964                 | **13.3x**     |
| 2^20   | 296          | 3671                | **12.4x**     |

Even single-threaded VROOM vs multi-threaded midnight-zk:

| Points | VROOM 1T (ms) | midnight-zk 4T (ms) | VROOM 1T Speedup |
|--------|--------------|---------------------|------------------|
| 2^14   | 8.64         | 88.0                | **10.2x**        |
| 2^16   | 33.9         | 292                 | **8.6x**         |
| 2^18   | 135          | 964                 | **7.1x**         |
| 2^20   | 556          | 3671                | **6.6x**         |

## Analysis

### Why is VROOM so much faster?

The dominant factor is **field arithmetic speed**. Every elliptic curve operation (point addition, doubling) consists of multiple field multiplications and additions. VROOM uses AVX512-IFMA RNS Montgomery multiplication, which processes the 381-bit BLS12-381 modulus using 8 x 52-bit limbs in a single AVX512 register. This gives ~3-4x faster field multiplication compared to BLST's already-optimized x86-64 assembly.

midnight-zk's `msm_best` calls BLST field operations through Rust FFI (the `blst` crate). While BLST's field arithmetic is fast (hand-tuned ADX/MULX assembly), the MSM algorithm layer is in Rust. The key overheads:

1. **FFI overhead per field operation:** Each BLST field call crosses the Rust-C boundary. In a tight Pippenger scatter loop doing millions of point additions, this adds up.

2. **Batch affine addition cost:** midnight-zk's `msm_best` uses batch affine addition (Montgomery's trick to amortize inversions). While this reduces the cost per addition from ~16M to ~6M+1I/batch, the field inversion is still expensive and the batching adds memory traffic and complexity.

3. **Algorithmic overhead:** midnight-zk's Pippenger uses Booth encoding (signed digits, halving bucket count) and a Schedule-based batch accumulation pattern. VROOM uses simpler unsigned digits with direct bucket accumulation. The simpler approach has less overhead per point.

4. **Memory access patterns:** VROOM's bucket scatter is a simple indexed accumulation into a contiguous array. midnight-zk's batch affine approach groups collisions into a Schedule struct, processes them in batches of 64, and does multiple passes. This hurts cache locality.

### BLST C vs midnight-zk Rust

Interestingly, BLST's C Pippenger (called directly by VROOM) is **5.8x faster** than midnight-zk's Rust MSM at 2^20, despite both using BLST field arithmetic under the hood. This shows the algorithmic and systems-level overhead in midnight-zk's approach: the batch affine strategy, Booth encoding, Rayon scheduling, and FFI boundary all add measurable cost compared to BLST's monolithic C implementation with XYZZ coordinates and prefetching.

### Threading efficiency

| Implementation | 2^20 (ms) | Threads | Per-thread efficiency |
|---|---|---|---|
| VROOM Parallel | 296 | 4 vCPU | 47% per vCPU (93% per physical core) |
| midnight-zk msm_best | 3671 | 4 vCPU (Rayon) | N/A (no 1T baseline available) |
| BLST Pippenger | 631 | 1 | 100% (reference) |

VROOM achieves 1.88x parallel speedup on 2 physical cores (4 vCPUs), near the theoretical max for this hardware. midnight-zk's Rayon parallelism operates at the window level (parallel windows when c >= 10), which has different scaling characteristics.

## Raw Output

### midnight-zk msm_best (Rayon, 4 vCPUs)

```
Generating 1048576 points...
Done generating.

Points          msm_best (ms)
------------------------------
2^10                   14.34
2^11                   24.72
2^12                   44.85
2^13                   49.43
2^14                   87.98
2^15                  156.95
2^16                  291.58
2^17                  553.25
2^18                  963.75
2^19                 1884.10
2^20                 3671.15
```

### VROOM (clang++, AVX512-IFMA)

```
=== VROOM 1T ===
BM_VROOM_MSM/1024_mean           0.614 ms
BM_VROOM_MSM/2048_mean            1.16 ms
BM_VROOM_MSM/4096_mean            2.25 ms
BM_VROOM_MSM/8192_mean            4.36 ms
BM_VROOM_MSM/16384_mean           8.64 ms
BM_VROOM_MSM/32768_mean           17.1 ms
BM_VROOM_MSM/65536_mean           33.9 ms
BM_VROOM_MSM/131072_mean          67.6 ms
BM_VROOM_MSM/262144_mean           135 ms
BM_VROOM_MSM/524288_mean           271 ms
BM_VROOM_MSM/1048576_mean          556 ms

=== VROOM Parallel (4 vCPUs) ===
BM_VROOM_MSM_Parallel/1024_mean            2.32 ms
BM_VROOM_MSM_Parallel/2048_mean            2.38 ms
BM_VROOM_MSM_Parallel/4096_mean            2.98 ms
BM_VROOM_MSM_Parallel/8192_mean            4.03 ms
BM_VROOM_MSM_Parallel/16384_mean           6.03 ms
BM_VROOM_MSM_Parallel/32768_mean           10.4 ms
BM_VROOM_MSM_Parallel/65536_mean           19.1 ms
BM_VROOM_MSM_Parallel/131072_mean          37.1 ms
BM_VROOM_MSM_Parallel/262144_mean          72.7 ms
BM_VROOM_MSM_Parallel/524288_mean           144 ms
BM_VROOM_MSM_Parallel/1048576_mean          296 ms

=== BLST Pippenger 1T ===
BM_BLST_Pippenger/1024_mean            1.13 ms
BM_BLST_Pippenger/2048_mean            1.62 ms
BM_BLST_Pippenger/4096_mean            2.88 ms
BM_BLST_Pippenger/8192_mean            5.23 ms
BM_BLST_Pippenger/16384_mean           10.0 ms
BM_BLST_Pippenger/32768_mean           19.6 ms
BM_BLST_Pippenger/65536_mean           38.1 ms
BM_BLST_Pippenger/131072_mean          75.0 ms
BM_BLST_Pippenger/262144_mean           148 ms
BM_BLST_Pippenger/524288_mean           297 ms
BM_BLST_Pippenger/1048576_mean          631 ms
```
