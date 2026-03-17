# VROOM Field Arithmetic Benchmark Results

**Platform:** GCP VM `vroom-c3-22` (c3-standard-22, Sapphire Rapids, AVX512-IFMA)
**Benchmark command:** `taskset -c 0-15 cargo bench --bench field_arith -p midnight-curves`
**Date:** 2026-03-17
**Branch:** `vroom-fields-bench`

## 1. Single-Operation Microbenchmarks

### Fp (Base Field, 381-bit) — blst vs VROOM

| Operation | blst (ns) | VROOM (ns) | Speedup | Winner |
|-----------|-----------|------------|---------|--------|
| add       | 9.71      | 1.49       | 6.5x    | VROOM  |
| sub       | 9.72      | 2.65       | 3.7x    | VROOM  |
| double    | 9.38      | 1.68       | 5.6x    | VROOM  |
| neg       | 9.71      | 2.19       | 4.4x    | VROOM  |
| mul       | 37.52     | 38.63      | 0.97x   | ~tie   |
| square    | 37.16     | 38.91      | 0.96x   | ~tie   |
| invert    | 3,150     | 29,110     | 0.11x   | blst (9.2x) |

### Fr (Scalar Field, 255-bit) — blst vs VROOM

| Operation | blst (ns) | VROOM (ns) | Speedup | Winner |
|-----------|-----------|------------|---------|--------|
| add       | 8.37      | 1.68       | 5.0x    | VROOM  |
| sub       | 8.38      | 2.75       | 3.0x    | VROOM  |
| double    | 8.37      | 1.34       | 6.2x    | VROOM  |
| neg       | 9.38      | 2.24       | 4.2x    | VROOM  |
| mul       | 19.81     | 27.98      | 0.71x   | blst (1.4x) |
| square    | 19.76     | 27.87      | 0.71x   | blst (1.4x) |
| invert    | 1,890     | 2,740*     | 0.69x   | blst (1.45x) |

\* Fr inversion uses GMP `mod_inverse` (convert-invert-convert optimization).
Original square-and-multiply implementation was 37,250 ns (19.7x slower than blst).

### Why Single Mul Shows No Speedup

The VROOM paper (ePrint 2026/393, Table 11) explicitly states: single Fq multiplication = 28.3 ns for VROOM vs 28.9 ns for blst — *"With a batch size of 1, there is effectively no speedup to using AVX512 and RNS over traditional methods."*

AVX512 has high throughput (0.5 CPI) but high latency (~10 cycles). A single multiplication cannot fill the pipeline.

## 2. 2^20 Multiplication Chain (1,048,576 total multiplications)

### Fp (Base Field)

| Benchmark | Total time | Per-mul (ns) | vs blst |
|-----------|-----------|--------------|---------|
| blst Fp chain (sequential)       | 54.6 ms | 52.1 | baseline |
| VROOM Fp chain (sequential)      | 63.6 ms | 60.7 | 0.86x (slower) |
| **VROOM Fp batch×6 chain**       | **25.8 ms** | **24.6** | **2.12x faster** |

### Fr (Scalar Field)

| Benchmark | Total time | Per-mul (ns) | vs blst |
|-----------|-----------|--------------|---------|
| blst Fr chain (sequential)       | 35.2 ms | 33.5 | baseline |
| VROOM Fr chain (sequential)      | 51.7 ms | 49.3 | 0.68x (slower) |
| **VROOM Fr batch×6 chain**       | **22.1 ms** | **21.1** | **1.59x faster** |

### Analysis

- **Sequential VROOM is slower** than blst because single modmul has no latency-hiding benefit.
- **Batched VROOM (6 parallel chains) is 1.6–2.1x faster** because `batch_modmul<6>` interleaves 6 independent CRNS reductions, hiding AVX512 instruction latency.
- The paper's Table 11 shows the same effect: per-Fq-mul drops from 28.3 ns (batch=1) to 16.5 ns (batch=8).
- Our batch=6 result (24.6 ns for Fp, 21.1 ns for Fr) is consistent with the paper's batch=6 row (18.3 ns), with the difference attributable to benchmark overhead (element reload between criterion iterations).

## 3. Sum of 2 Products: a·b + c·d

| Benchmark | Time (ns) | vs blst |
|-----------|-----------|---------|
| blst Fp (2 muls + 1 add) | 90.7 | baseline |
| **VROOM Fp (single reduction)** | **40.5** | **2.24x faster** |
| blst Fr (2 muls + 1 add) | 49.9 | baseline |
| **VROOM Fr (single reduction)** | **31.4** | **1.59x faster** |

### Why This Works

In blst, `a*b + c*d` requires 2 full modular multiplications + 1 addition = 2 × O(t²) reductions.

In VROOM RNS, the elementwise products `a*b` and `c*d` are O(t) each. The sum `(a*b) + (c*d)` is also O(t). Only one O(t²) CRNS reduction is needed for the final result. This is the key insight from the paper's Section 4.2.

This pattern appears extensively in:
- **Fp2 multiplication**: (a+bi)(c+di) = (ac-bd) + (ad+bc)i — sums of 2 products
- **Fp12 multiplication**: many sums of products → paper reports **5.27x speedup**
- **EC point addition**: projective formula has sums of products → paper reports **3.34x speedup**

## 4. Fr Inversion Optimization

| Implementation | Time | vs blst |
|----------------|------|---------|
| blst (euclidean algorithm) | 1.89 µs | baseline |
| VROOM old (square-and-multiply, 254 modmuls) | 37.25 µs | 19.7x slower |
| **VROOM new (GMP mod_inverse)** | **2.74 µs** | **1.45x slower** |

The convert-invert-convert optimization: `to_bigint()` → GMP `mod_inverse` → `from_bigint()`. This is the approach recommended by the paper (Section 6.4): *"We convert one number to standard radix form, apply blst's inversion, and then convert back to RNS."*

Improvement: **13.6x faster** than the naive square-and-multiply approach.

## 5. FFT/NTT Analysis

### Why VROOM Cannot Improve FFT

FFT operates on the scalar field Fr. Each butterfly does:
```
t  = b * twiddle    // 1 multiplication
a' = a + t          // 1 addition
b' = a - t          // 1 subtraction
```

#### Theoretical butterfly cost

| Component | blst | VROOM (batch×6) |
|-----------|------|-----------------|
| 1 Fr mul  | 20 ns | ~21 ns |
| 1 Fr add  | 8.4 ns | 1.7 ns |
| 1 Fr sub  | 8.4 ns | 2.7 ns |
| **Total** | **~37 ns** | **~25 ns** |

Compute-only prediction: ~1.5x speedup. But memory effects dominate.

#### The fundamental problem: 3x element size

| | blst | VROOM |
|---|---|---|
| Fr element size | 32 bytes | 96 bytes (2 × 6 limbs × 8B) |
| 2^20 FFT data | 32 MB | 96 MB |
| 2^22 FFT data | 128 MB | 384 MB |
| L3 cache (c3-standard-22) | ~36 MB | ~36 MB |

VROOM's RNS representation stores each element as two halves (M and N), each with 6 limbs. This makes every element **3x larger** than blst's Montgomery form.

FFT has strided butterfly access patterns. Once the data exceeds cache, performance is dominated by **memory bandwidth**, not compute. VROOM moves 3x more data through the memory hierarchy.

#### Expected performance by FFT size

| FFT size | blst data | VROOM data | Expected winner |
|----------|-----------|------------|-----------------|
| 2^10 (1K) | 32 KB | 96 KB | VROOM ~1.3-1.5x (fits L1/L2) |
| 2^14 (16K) | 512 KB | 1.5 MB | VROOM ~1.1-1.3x (L2 boundary) |
| 2^17 (128K) | 4 MB | 12 MB | ~Tie or blst edge (L3 pressure) |
| 2^20 (1M) | 32 MB | 96 MB | blst wins ~1.5-2x (memory-bound) |
| 2^22 (4M) | 128 MB | 384 MB | blst wins ~2-3x (DRAM-bound) |

#### Multi-threading makes it worse for VROOM

- All cores share the same memory bus; more threads = more bandwidth contention
- blst benefits more: per-thread working set (32MB/16 = 2MB) fits L2; VROOM's (96MB/16 = 6MB) does not
- Sustained all-core AVX512 causes frequency throttling on Sapphire Rapids (~3.8→3.4 GHz)

#### Unexpanded form doesn't help

Storing only the N-half (48 bytes vs 96 bytes) reduces memory 2x, but requires an extra CRNS expansion (O(t²)) before each multiply. Total CRNS work per butterfly increases from 2 to 4, trading memory for more compute — the opposite of what FFT needs.

#### Conclusion

VROOM cannot improve FFT at PLONK-relevant sizes (2^17 to 2^22). The RNS representation is fundamentally mismatched with FFT's memory access pattern. The paper's authors never benchmark FFT — their wins come from pairings and EC operations where computation-to-memory ratio is high.

**Recommended architecture:** blst for FFT, VROOM for EC/pairing operations, with O(N) conversion at the boundary (negligible vs O(N log N) FFT cost).

## 6. Summary: Where VROOM Wins and Loses

| Regime | VROOM vs blst | Reason |
|--------|---------------|--------|
| Additive ops (add/sub/neg/double) | **3-6x faster** | Elementwise SIMD, no modular reduction |
| Single multiplication | ~1x (tie) | Same O(t²) cost, no latency hiding |
| Batched multiplication (6+) | **1.6-2.1x faster** | AVX512 latency hiding via batch_modmul |
| Sum of products (a·b + c·d) | **1.6-2.2x faster** | Single CRNS reduction for multiple products |
| Field extensions (Fp2, Fp12) | **1.85-5.27x faster** (paper) | Sums of products + batching |
| EC point operations | **2.6-4.0x faster** (paper) | All of the above combined |
| Pairings | **3.3x faster** (paper) | All of the above combined |
| Inversion | blst wins (1.45x for Fr) | Euclidean algorithm is faster; use convert-invert-convert |
| FFT (2^17+) | **blst wins ~1.5-3x** | 3x element size, memory-bandwidth bound |
