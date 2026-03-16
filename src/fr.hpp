#pragma once
#include "bounded_ring.hpp"

// BLS12-381 scalar field modulus r
// r = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
inline const char* bls12_381_scalar_modulus_hex = "73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001";

// Fr uses 6 limbs of 52 bits each (RNS_BITS = 311, needs 271 for default bounds).
// 5 limbs is impossible: 255 + ceil_log2(9*4309) = 255 + 16 = 271 > 259.
using FrRing = BoundedRing<255, 6, 52, -1932, 2377>;

inline FrRing make_fr_ring() {
    BigInt r(bls12_381_scalar_modulus_hex, 16);
    return FrRing(r);
}
