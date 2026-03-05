# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

VROOM implements RNS (Residue Number System) Montgomery multiplication using AVX512-IFMA instructions for BLS12-381 elliptic curve pairing operations. It is a C++20 header-only library with benchmarks using Google Benchmark and correctness tests validated against GMP and BLST.

This is research-grade code, not audited for production use.

## Build Commands

Build BLST dependency first, then the main source:

```bash
cd blst && make && cd ../src && make
```

### Key make targets (run from `src/`):

| Target | Description |
|---|---|
| `make test_pairing_avx` | Pairing test (AVX512 IFMA) |
| `make test_pairing` | Pairing test (integer fallback, no AVX512 needed) |
| `make bench_pairing_50bit` | Primary benchmark (50-bit params, AVX512 IFMA) |
| `make test_elementwise_reduce` | Elementwise reduce test |
| `make test_inversion` | Inversion test (requires GMP) |
| `make test_ec` | Elliptic curve test |
| `make test_scalar_mult` | Scalar multiplication test |

Run a single test: `cd src && make test_pairing && ./test_pairing`

### Compiler

- Preferred: `clang-21` (or any modern clang). Set via `CXX=clang++ make <target>`.
- GCC works but is significantly slower at runtime.
- C++20 (`-std=c++20` with clang, `-std=c++2a` with GCC).

### Dependencies

- **GMP/GMPXX**: Required for tests (`-lgmp -lgmpxx`)
- **Google Benchmark**: Required for benchmarks (`-lbenchmark -lpthread`)
- **BLST**: Built from `blst/` directory, produces `libblst_small.a`

## Architecture

### Compile-Time Bounds Tracking

The central design pattern is **compile-time bounds propagation**. Every arithmetic operation tracks value ranges through template parameters so overflow is caught at compile time:

- `Bounds<lower, upper>` (`bounds.hpp`) — interval arithmetic on `int64_t` template params
- `BoundedElement<BoundsType, limbs, element_bits>` (`bounded_element.hpp`) — a vector of limbs with static bounds, wrapping `AVXVector<limbs>`
- `SmallElement` / `ExpandedElement` (`ring_element.hpp`) — RNS-level wrappers with `RNSBounds` tracking bounds across both RNS channels
- `BoundedRing` (`bounded_ring.hpp`) — the ring of integers mod p, parameterized by modulus bits, limb count, element bits, and max negative/positive RNS bounds

### AVX512 / Fallback Dispatch

`cpu/vector/vector_impl.hpp` conditionally includes either `avx512ifma.hpp` or `fallback.hpp` based on `USE_AVX512_FALLBACK`. Both provide the same `AVXVector<limbs>` interface. The Makefile uses `-mavx512ifma` for AVX targets and `-DUSE_AVX512_FALLBACK` for fallback targets.

### Pairing Pipeline (BLS12-381)

Built as a tower of algebraic types, each templated on the ring:

1. **FP** — base field elements via `BoundedRing` (RNS Montgomery representation)
2. **FP2** (`fp2.hpp`) — quadratic extension, pair of FP elements
3. **FP12** (`fp12.hpp`) — degree-12 extension for pairing target group
4. **EC** (`ec.hpp`) — projective/affine elliptic curve points (G1 over FP, G2 over FP2)
5. **Miller loop** (`miller.hpp`) — computes the Miller function
6. **Final exponentiation** (`final_exponentiation.hpp`) — raises Miller output to the pairing power
7. **Pairing** (`pairing.hpp`) — combines Miller loop + final exponentiation

### CPU Layer (`cpu/`)

- `cpu/vector/` — AVX512-IFMA and fallback vector implementations, base change, multiplication
- `cpu/reduction/` — Montgomery reduction
- `cpu/precompute/` — GMP-based precomputation of ring constants

### BLST Integration (`blst/`)

A trimmed subset of the [BLST](https://github.com/supranational/blst) library providing FP inversion (used as `FP12BLSTInverter` / `FP12RNS_BLSTInverter`). Built as `libblst_small.a`. Contains x86-64 assembly (`assembly_small.S`) and C source.

### Test Data

`test_data/` contains precomputed test vectors (`.hpp` files) for pairing values, test points, and FP12 ring test values.

## Taking notes

When required to take notes with regard to the conversation, write the conversation contents into a markdown file with the specified name under the folder `./notes`.
