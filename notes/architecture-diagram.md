# VROOM Architecture Diagram

## High-Level System Overview

```
+===========================================================================+
|                           VROOM Pairing Pipeline                          |
|                                                                           |
|  pairing.hpp                                                              |
|  +---------------------------------------------------------------------+  |
|  |                      pairing(P, Q, inverter, ...)                   |  |
|  |                                                                     |  |
|  |   +---------------------------+   +-------------------------------+ |  |
|  |   |      Miller Loop          |   |    Final Exponentiation       | |  |
|  |   |      miller.hpp           |   |    final_exponentiation.hpp   | |  |
|  |   |                           |   |                               | |  |
|  |   |  miller_loop(ring,        |   |  final_exp(f, inverter,       | |  |
|  |   |    fp2ring, fp12ring)     |   |    fp12ring, ring)            | |  |
|  |   +-------------+-------------+   +-------+-----------+-----------+ |  |
|  |                 |                         |           |             |  |
|  +---------------------------------------------------------------------+  |
|                    |                         |           |                 |
+--------------------|-------------------------|-----------|----------------+
                     |                         |           |
    +----------------v-------------------------v-----------v-----------+
    |                    Algebraic Type Tower                          |
    |                                                                  |
    |  +------------------------------------------------------------+  |
    |  |  FP12        fp12.hpp                                      |  |
    |  |  Degree-12 extension field (target group of pairing)       |  |
    |  +------------------------------+-----------------------------+  |
    |                                 |                                |
    |  +------------------------------v-----------------------------+  |
    |  |  FP2         fp2.hpp                                       |  |
    |  |  Quadratic extension field  (x + y*u)                     |  |
    |  +------------------------------+-----------------------------+  |
    |                                 |                                |
    |  +------------------------------v-----------------------------+  |
    |  |  FP (Base Field)   bounded_ring.hpp                        |  |
    |  |  BoundedRing<381, 8, 52, ...>                              |  |
    |  |  RNS Montgomery representation of F_p                      |  |
    |  +------------------------------------------------------------+  |
    +------------------------------------------------------------------+
                     |
    +----------------v-------------------------------------------------+
    |               Elliptic Curve Layer     ec.hpp                    |
    |                                                                  |
    |  ProjectivePoint<FP>   -- G1 points (over base field)           |
    |  ProjectivePoint<FP2>  -- G2 points (over quadratic extension)  |
    |  AffinePoint<FP/FP2>   -- Affine representations                |
    |                                                                  |
    |  scalar_mult.hpp  -- GLV/GLS scalar multiplication              |
    +------------------------------------------------------------------+


## Compile-Time Bounds Tracking (Core Design Pattern)

```
  Bounds<lower, upper>                     bounds.hpp
  Compile-time interval arithmetic         [lo, hi] + [lo, hi] = ...
         |
         v
  BoundedElement<Bounds, limbs, bits>      bounded_element.hpp
  AVXVector with static bounds             Overflow caught at compile time
         |
         v
  SmallElement / ExpandedElement           ring_element.hpp
  RNS-level wrappers                       Tracks bounds across both RNS channels
         |
         v
  BoundedRing                              bounded_ring.hpp
  Full ring mod p                          Montgomery mult, reduce, change base
```


## CPU Backend (Hardware Abstraction)

```
  cpu/
  +------------------------------------------------------------------+
  |                                                                  |
  |  vector/                                                         |
  |  +------------------------------------------------------------+  |
  |  |  vector_impl.hpp  ---- dispatches via USE_AVX512_FALLBACK  |  |
  |  |       |                        |                            |  |
  |  |       v                        v                            |  |
  |  |  avx512ifma.hpp          fallback.hpp                       |  |
  |  |  __m512i intrinsics      uint64_t arrays                    |  |
  |  |  _mm512_madd52lo_epu64   portable emulation                 |  |
  |  |                                                             |  |
  |  |  multiplication.hpp  -- modular multiplication              |  |
  |  |  changebase.hpp      -- RNS base conversion                 |  |
  |  |  conversion.hpp      -- format conversions                  |  |
  |  +------------------------------------------------------------+  |
  |                                                                  |
  |  reduction/                                                      |
  |  +------------------------------------------------------------+  |
  |  |  montgomery.hpp  -- Montgomery reduction (2-channel RNS)   |  |
  |  +------------------------------------------------------------+  |
  |                                                                  |
  |  precompute/                                                     |
  |  +------------------------------------------------------------+  |
  |  |  precompute.hpp   -- RNS moduli & constant generation      |  |
  |  |  gmp_wrapper.hpp  -- GMP BigInt utilities                   |  |
  |  +------------------------------------------------------------+  |
  |                                                                  |
  +------------------------------------------------------------------+
```


## BLST Integration

```
  blst/
  +------------------------------------------------------------------+
  |  libblst_small.a    (trimmed BLST subset)                        |
  |                                                                  |
  |  server_small.c     -- C source (pairing, field ops, recip)     |
  |  assembly_small.S   -- x86-64 assembly (field arithmetic)       |
  |                                                                  |
  |  Used by src/ for:                                               |
  |    FP12BLSTInverter      -- FP12 inversion via BLST             |
  |    FP12RNS_BLSTInverter  -- RNS FP12 inversion (hybrid)         |
  |    conversion_inversion.hpp  -- RNS <-> BLST format conversion  |
  +------------------------------------------------------------------+
```


## Data Flow: Full Pairing Computation

```
  G1 AffinePoint (P)     G2 AffinePoint (Q)
        |                       |
        +----------+------------+
                   |
                   v
           Miller Loop (miller.hpp)
           Iterates over BLS12-381 parameter bits
           Each step: double_step / add_step
           Uses FP, FP2, FP12 arithmetic
                   |
                   v
           FP12 element (f)
                   |
                   v
           Final Exponentiation (final_exponentiation.hpp)
           f^((p^12 - 1) / r)
           = easy_part * hard_part
           Uses FP12 squaring, Frobenius, inversion
                   |
                   v
           FP12 pairing result
```


## File Dependency Graph (src/)

```
  pairing.hpp
    +-- miller.hpp
    +-- final_exponentiation.hpp
    +-- conversion_inversion.hpp
          +-- inversion.hpp

  miller.hpp, final_exponentiation.hpp
    +-- fp12.hpp
          +-- fp2.hpp
                +-- ring_element.hpp
                      +-- bounded_element.hpp
                      |     +-- bounds.hpp
                      |     +-- cpu/vector/vector_impl.hpp
                      +-- preprocess.hpp

  bounded_ring.hpp
    +-- ring_element.hpp
    +-- elementwise_reduce.hpp
    +-- change_base.hpp
    +-- cpu/precompute/gmp_wrapper.hpp
    +-- cpu/vector/multiplication.hpp

  ec.hpp  -- ProjectivePoint, AffinePoint (used by miller.hpp)
  scalar_mult.hpp  -- GLV/GLS scalar multiplication
```
