# Compiler Issue: Pippenger MSM Performance

**Date:** 2026-03-13

## Issue

The Booth-encoded Pippenger scatter loop in `src/pippenger.hpp` runs **2.3x slower** when compiled with g++ 11.4 compared to clang 14.

| Compiler | VROOM MSM (ST) | BLST |
|----------|---------------|------|
| g++ 11.4 | 1398 ms | 628 ms |
| clang 14 | 614 ms | 631 ms |

## Root Cause

g++ produces significantly worse code for the scatter loop, which involves:
- Booth-encoded signed digit reads (`int32_t` array)
- Conditional point negation (`AffPoint` copy + `negate_affine`)
- Random bucket access with `__builtin_prefetch`
- Mixed-point addition (`add_mixed_point`)

The CLAUDE.md already documents this: "GCC works but is significantly slower at runtime."

## Resolution

Always compile with clang for benchmarks and production:

```bash
CXX=clang++ make bench_msm_avx
```

## VM Details

- VM: msm-bench-large (n2-standard-16, Sapphire Rapids)
- g++: Ubuntu 11.4.0
- clang++: Ubuntu clang 14.0.0
- n = 2^20 points, BLS12-381 G1
