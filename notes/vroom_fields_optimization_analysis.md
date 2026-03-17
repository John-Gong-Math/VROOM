# VROOM Field Arithmetic: Performance Analysis & Optimization Opportunities

## Current Benchmark Results (GCP VM vroom-c3-22)

### Fp (Base Field, 381-bit)

| Operation | blst (native) | VROOM | Ratio |
|-----------|--------------|-------|-------|
| add       | 9.71 ns      | 1.49 ns | **6.5x VROOM** |
| sub       | 9.72 ns      | 2.65 ns | **3.7x VROOM** |
| double    | 9.38 ns      | 1.68 ns | **5.6x VROOM** |
| neg       | 9.71 ns      | 2.19 ns | **4.4x VROOM** |
| mul       | 37.52 ns     | 38.63 ns | ~1x tie |
| square    | 37.16 ns     | 38.91 ns | ~1x tie |
| invert    | 3.15 µs      | 29.11 µs | **9.2x blst** |

### Fr (Scalar Field, 255-bit)

| Operation | blst (native) | VROOM | Ratio |
|-----------|--------------|-------|-------|
| add       | 8.37 ns      | 1.68 ns | **5.0x VROOM** |
| sub       | 8.38 ns      | 2.75 ns | **3.0x VROOM** |
| double    | 8.37 ns      | 1.34 ns | **6.2x VROOM** |
| neg       | 9.38 ns      | 2.24 ns | **4.2x VROOM** |
| mul       | 19.81 ns     | 27.98 ns | **1.4x blst** |
| square    | 19.76 ns     | 27.87 ns | **1.4x blst** |
| invert    | 1.89 µs      | 37.25 µs | **19.7x blst** |

## Analysis: Why These Results Match the Paper

The VROOM paper (ePrint 2026/393) explicitly states in Table 11:

> **Batch size 1: Fq multiplication = 28.3 ns** (vs blst's 28.9 ns — no speedup)

This is exactly what we observe. The paper explains: *"With a batch size of 1, there is effectively no speedup to using AVX512 and RNS over traditional methods."*

VROOM's advantage comes from **three sources that our microbenchmark doesn't exercise**:

## Optimization 1: Batched Multiplication (Latency Hiding)

**Paper Table 11:** Amortized time per Fq multiplication by batch size:

| Batch | Time/mul |
|-------|----------|
| 1     | 28.3 ns  |
| 2     | 23.3 ns  |
| 4     | 19.3 ns  |
| 6     | 18.3 ns  |
| 8     | 16.5 ns  |
| 10    | 15.8 ns  |

AVX512 instructions have **high throughput (0.5 CPI) but high latency (~10 cycles)**. Batching multiple independent multiplications lets the CPU pipeline fill the latency gaps. With batch=8, VROOM achieves **1.75x speedup over blst** on Fq multiplication.

**Implementation:** VROOM already has `batch_modmul<N>()`. We should benchmark it:

```cpp
// Batch of 6 independent multiplications
ring.batch_modmul<6>(a_array, b_array);
```

**Expected gain:** ~1.75x over blst for Fp mul at batch=8.

## Optimization 2: Sum-of-Products (Shared Reduction)

The paper's key insight: a single modular multiplication costs `2t² + 13t`, where `2t²` is the unreduced product and `13t` is the reduction. For a **sum of k products**, schoolbook costs `(k+1)t²` (k products + 1 reduction), but VROOM costs only `2t² + (12.5+k)t`.

For `k=2` (a·b + c·d): schoolbook = 3t², VROOM = 2t² — **1.5x speedup**.

This is because in RNS:
- Elementwise products are O(t) — you can accumulate them
- Only one CRNS reduction (O(t²)) is needed for the final sum

**This is the dominant win for field extensions:**
- Fp2 multiply: 1.85x over blst
- Fp12 multiply: **5.27x over blst**

**Implementation:** VROOM's `batch_reduce_expand` already supports accumulating multiple products before reducing:

```cpp
auto prod1 = a * b;  // elementwise, O(t)
auto prod2 = c * d;  // elementwise, O(t)
auto sum = prod1 + prod2;  // elementwise add
auto result = ring.batch_reduce_expand<1>({ring.ready<MAX_ADD>(sum)})[0];  // one reduction
```

## Optimization 3: Unexpanded Form Computation (Section 4.3)

For expressions like `(a+b) · (c+d)`, instead of expanding all 4 values to full RNS-MN (4 CRNS calls), we can:
1. Add `a+b` and `c+d` in **unexpanded form** (RNS-N only, half the work)
2. Expand only 2 values before multiplication (2 CRNS calls)

This saves 2 expensive CRNS calls per such expression, which is significant since CRNS is `t² + 4t` operations.

## Optimization 4: 50-bit Residues Instead of 52-bit

The paper specifically recommends 50-bit residues for BLS12-381 (Section 5.2):

> "It is more efficient to use 50-bit residues, which works since 8 × 50 = 400 > 395. This allows us to utilize the redundant form outputs from Section 5.1.1 and Section 5.1.2, eliminating many operations."

With 52-bit moduli, Montgomery reduction output has bounds `Bounds<0, (1+⌈k/4⌉)m>` which may require an extra elementwise reduction. With 50-bit moduli, the output `Bounds<0, 4m>` fits in 52 bits (2 bits of headroom), eliminating the reduction step entirely for k ≤ 12.

**Current code:** `BoundedRing<381, 8, 52, ...>` — uses 52-bit elements.
**Paper recommendation:** Use 50-bit elements for BLS12-381.

This would require creating a new ring type and adjusting the LOG_MULTIPLES parameter.

## Optimization 5: Fr Inversion via blst (Convert-Invert-Convert)

Our Fr inversion uses naive square-and-multiply: 254 squarings + ~127 multiplications = ~381 modmuls at ~28ns each ≈ 10.7µs theoretical minimum (we measure 37.25µs due to BigInt overhead in the bit-scanning loop).

The paper's approach for inversion (Section 6.4):
> "We convert one number to standard radix form, apply blst's inversion, and then convert back to RNS. This results in a total 1.55× speedup versus blst's inversion algorithm."

For Fr, blst's euclidean inversion is 1.89µs. Adding two RNS↔radix conversions (~28ns each), the total would be ~1.95µs — a **19x improvement** over our current 37.25µs.

For Fp, the addition-chain inversion (425 steps ≈ 425 × 28ns ≈ 11.9µs) vs blst's 3.15µs + 2×28ns ≈ 3.2µs. The convert-invert-convert approach would be **3.6x faster** than our current 29.1µs.

**Implementation:** Use `ring.to_bigint()` / `ring.from_bigint()` to convert, call blst's inversion, convert back.

## Optimization 6: The Real Benchmark Should Be Higher-Level

The paper's headline results aren't from field microbenchmarks — they're from **protocol-level operations**:

| Operation | blst | VROOM | Speedup |
|-----------|------|-------|---------|
| G1 point add | 586 ns | 175 ns | **3.34x** |
| G2 point add | 1420 ns | 356 ns | **3.99x** |
| Fp2 multiply | 87 ns | 46.9 ns | **1.85x** |
| Fp12 multiply | 1580 ns | 300 ns | **5.27x** |
| Pairing | 504 µs | 152 µs | **3.29x** |
| BLS verify | 723 µs | 211 µs | **3.43x** |

These speedups come from the combination of:
- Lazy reduction (sums of products share reductions)
- Unexpanded form (defer expensive CRNS expansion)
- Batched reductions (hide AVX512 latency)
- 50-bit residues (eliminate elementwise reductions)

## Recommended Benchmarks to Add

1. **Batch Fp/Fr multiply** (batch sizes 1, 2, 4, 6, 8) — shows latency hiding
2. **Sum of 2 products** (a·b + c·d) — shows shared reduction
3. **Fp2 multiply** (complex multiplication using 3 Fp muls) — shows field extension advantage
4. **EC point addition** — shows the full protocol-level speedup
5. **Inversion via convert-blst-convert** — dramatic improvement over current approach

## Summary: What Matters

For **isolated field operations**, VROOM wins on additive ops (3-6x) and ties/loses slightly on multiplicative ops. This is expected and explained by the paper.

The real question is: **where in the midnight-zk stack would VROOM deliver wins?**

The answer is: anywhere that involves **batched multiplications**, **sums of products**, or **field extension arithmetic** — which is essentially all of elliptic curve cryptography and pairing computation. The paper demonstrates 3-5x speedups at the protocol level, even though individual Fq multiplication shows no improvement.
