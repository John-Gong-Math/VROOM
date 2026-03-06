#pragma once
#include "conversion_inversion.hpp"
#include <vector>

// Montgomery batch inversion using Montgomery's trick.
// Computes inverses of all elements in a single batch using only one
// field inversion (via BLST) plus O(n) multiplications.
//
// NOTE: G1 only for V1. G2 would require Fp2 batch inversion (norm trick).
template<class Ring>
std::vector<typename Ring::StandardElement> batch_invert(
    const Ring &ring,
    const std::vector<typename Ring::StandardElement> &elements
) {
    using StdElem = typename Ring::StandardElement;
    size_t n = elements.size();
    if (n == 0) return {};

    std::vector<StdElem> partials(n);
    std::vector<StdElem> result(n);

    // Track which elements are zero (point at infinity).
    // Substitute 1 for zero elements to keep the product chain valid.
    auto one = ring.one();

    // Forward pass: partials[i] = product of elements[0..i]
    // For zero elements, substitute 1 so the product stays non-zero.
    auto is_zero = [&](const StdElem &e) -> bool {
        auto bi = ring.to_bigint(e);
        return bi == BigInt(0);
    };

    std::vector<bool> zero_flags(n);
    zero_flags[0] = is_zero(elements[0]);
    partials[0] = zero_flags[0] ? one : elements[0];

    for (size_t i = 1; i < n; i++) {
        zero_flags[i] = is_zero(elements[i]);
        StdElem ei = zero_flags[i] ? one : elements[i];
        partials[i] = ring.modmul(partials[i - 1], ei);
    }

    // Single inversion of the accumulated product
    StdElem inv = invert_via_blst(ring, partials[n - 1]);

    // Backward pass: recover individual inverses
    for (size_t i = n - 1; i > 0; i--) {
        if (zero_flags[i]) {
            result[i] = ring.zero();
        } else {
            result[i] = ring.modmul(inv, partials[i - 1]);
        }
        StdElem ei = zero_flags[i] ? one : elements[i];
        inv = ring.modmul(inv, ei);
    }
    // i == 0
    if (zero_flags[0]) {
        result[0] = ring.zero();
    } else {
        result[0] = inv;
    }

    return result;
}
