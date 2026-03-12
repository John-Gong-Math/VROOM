# MSM Implementation Comparison: VROOM Parallel vs BLST Pippenger

**Date:** 2026-03-09

## Overview

Both implementations use **Pippenger's bucket method** for multi-scalar multiplication on BLS12-381 G1. This document compares the architectural and algorithmic differences.

## Algorithm Structure

Both share the same high-level flow: decompose scalars into windows, scatter points into buckets per window, integrate buckets via running sum, then combine windows with doublings.

| Aspect | VROOM `msm_parallel` | BLST `blst_p1s_mult_pippenger` |
|---|---|---|
| **Window size heuristic** | `choose_wbits` (`pippenger.hpp:15`) — identical to BLST | `pippenger_window_size` (`multi_scalar.c:287`) |
| **Window iteration** | High-to-low, doubling result between windows | Same |
| **Bucket integration** | Running sum + weighted sum (`integrate_buckets`) | Same pattern (`ptype_integrate_buckets`) |

The window size functions are functionally identical — both use the same log2-based breakpoints (>12→-3, >8→-2, >4→-1).

## Digit Encoding

| | VROOM | BLST |
|---|---|---|
| **Encoding** | **Unsigned digits** — extracts raw `wbits` bits, digit range [0, 2^w - 1], needs 2^w - 1 buckets | **Booth (signed) digits** — digit range [-(2^(w-1)), 2^(w-1)], needs only 2^(w-1) buckets |
| **Impact** | ~2× more buckets for same window size, more memory, more integration work | Half the buckets; negative digits handled by subtracting from bucket instead of adding |

Booth encoding in BLST (`ec_mult.h:46`) halves bucket count, reducing both memory and the O(2^w) integration cost. VROOM's comment at `pippenger.hpp:10` acknowledges this is a "V1 simplicity" choice.

## Coordinate Systems & Point Addition Cost

| | VROOM | BLST |
|---|---|---|
| **Bucket coordinates** | **Projective (X, Y, Z)** — 3 field elements per bucket | **XYZZ (X, Y, ZZ, ZZZ)** — 4 field elements, but faster mixed addition |
| **Mixed add (affine into bucket)** | `PointMixedAdd` (`ec.hpp:172`): Algorithm 7 adapted with Z_Q=1. **12M** (6 batch-reduced muls + 6 post-expand muls) + adds/subs | `XYZZ_DADD_AFFINE`: **8M + 2S** with inline doubling detection |
| **Proj+Proj add (bucket merge)** | `PointAdd` (`ec.hpp:114`): **12M** via Algorithm 7 from eprint 2015/1060 | Same cost class but in XYZZ coords |
| **Doubling** | `PointDouble` (`ec.hpp:220`): Algorithm 9, **4M + 4S** (batch-reduced) | Handled inline within DADD_AFFINE when P==Q |

BLST's XYZZ coordinates are specifically optimized for Pippenger: mixed addition with affine input is cheaper (8M+2S vs VROOM's ~12M). VROOM compensates by batching Montgomery reductions via AVX512.

## Field Arithmetic

| | VROOM | BLST |
|---|---|---|
| **Representation** | **RNS Montgomery** — 381-bit modulus split across 8×52-bit limbs in two RNS channels | **Standard Montgomery** — 6×64-bit limbs |
| **Multiplication** | AVX512-IFMA `vpmadd52` instructions; multiplications are batched (`batch_reduce`, `batch_expand`) to amortize reduction cost | x86-64 assembly (`mulx_mont_384` in `assembly_small.S`), using MULX/ADCX/ADOX (ADX extension) |
| **Batch amortization** | Core advantage — 6 multiplications reduced together in `batch_reduce` saves ~30% per-mul cost by sharing base-change overhead | No batching; each multiplication is independent |

VROOM's point addition formulas in `ec.hpp` are carefully structured to group multiplications into batches of 3–6 for `batch_reduce` / `batch_reduce_expand` calls. This is where the AVX512-IFMA speedup materializes.

## Parallelism

| | VROOM | BLST |
|---|---|---|
| **Threading** | **Multi-threaded** — splits points across threads within each window (`pippenger.hpp:182`) | **Single-threaded** |
| **Strategy** | Each thread gets its own bucket array, scatters its point chunk independently, then buckets are merged sequentially | N/A |
| **Thread count** | `hardware_concurrency()` by default, capped at `npoints/64` | N/A |
| **Overhead** | Per-window: spawn threads → scatter → join → merge buckets → integrate. Threads are re-created each window (no thread pool) | None |

The parallel merge phase (`pippenger.hpp:279-289`) is sequential: thread 0's buckets are updated by iterating through each other thread's buckets. For `num_threads × nbuckets` bucket entries, this is O(T × 2^w) point additions — a non-trivial cost that partially offsets the parallelism gain.

## Memory & Cache Behavior

| | VROOM | BLST |
|---|---|---|
| **Bucket memory** | `num_threads × (2^w - 1)` projective points (3 RNS elements each, 8×52-bit × 2 channels = large) | `2^(w-1)` XYZZ points (4 × 384-bit = 192 bytes each) |
| **Prefetching** | None | `ptype_prefetch(buckets, wnxt, cbits)` — prefetches next bucket during scatter |
| **Occupied tracking** | Explicit `uint8_t[]` flags to skip empty buckets | Buckets zeroed after integration; implicit zero = identity |
| **Scratch reuse** | Per-thread bucket arrays reused across windows | Buckets cleared in-place during `integrate_buckets` |

BLST's Booth encoding means ~half the buckets, better cache utilization. BLST also prefetches the next bucket during scatter — important since bucket access is essentially random. VROOM's per-thread replication multiplies the memory footprint by the thread count.

## Inversion (for batch affine conversion)

Both use Montgomery's trick for batch inversion, but the contexts differ:

| | VROOM | BLST |
|---|---|---|
| **Implementation** | `batch_inversion.hpp` — forward/backward pass with single `invert_via_blst()` call | `POINTS_TO_AFFINE_IMPL` in `multi_scalar.c:16` — same algorithm, chunked for cache |
| **Used in MSM?** | `batch_affine.hpp` exists but is **not used** in the Pippenger path — buckets stay projective | Used for precomputed table normalization |

## Small-case Dispatch

| | VROOM | BLST |
|---|---|---|
| **n = 0** | Returns identity | Returns identity |
| **n = 1** | Windowed scalar mult (4-bit window) | `POINTonE1_mult_w5` (5-bit window) |
| **n < 32** | Naive: individual scalar mults + accumulate | Precomputed wbits tables if memory allows, else Pippenger |
| **n < 1024** | Falls back to single-threaded Pippenger | N/A (always single-threaded) |

## Summary of Key Tradeoffs

| Advantage | VROOM | BLST |
|---|---|---|
| **Field mul throughput** | AVX512-IFMA + batched reductions amortize cost across 6 muls | Hand-tuned x86-64 assembly (ADX) but no batching |
| **Parallelism** | Multi-threaded scatter | Single-threaded |
| **Bucket efficiency** | 2× more buckets (unsigned digits) | Half the buckets (Booth encoding) |
| **Mixed addition cost** | ~12M (complete formula) | 8M+2S (XYZZ specialized) |
| **Cache optimization** | No prefetching | Prefetching next bucket |
| **Thread overhead** | Thread creation per window + sequential merge | None |

VROOM bets on **raw arithmetic throughput** (AVX512-IFMA batching + threading) to overcome its algorithmic disadvantages (more buckets, costlier point additions, no prefetching). BLST bets on **algorithmic efficiency** (Booth encoding, XYZZ coords, prefetching, assembly) within a single thread.
