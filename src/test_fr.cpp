#include "../cpu/precompute/gmp_wrapper.hpp"
#include "fr.hpp"
#include <iostream>
#include <cassert>

int main() {
    std::cout << "=== Testing Fr (BLS12-381 Scalar Field) ===" << std::endl;

    BigInt r(bls12_381_scalar_modulus_hex, 16);
    auto ring = make_fr_ring();

    // Test 1: basic ring construction sanity check
    // Note: check_all() is skipped because check_matrix_bound assumes limbs == indigits,
    // which only holds for the 381-bit base field (limbs=8, indigits=8), not Fr (limbs=6, indigits=5).
    std::cout << "\nTest 1: Ring construction + one() sanity" << std::endl;
    assert(ring.to_bigint(ring.one()) == 1);
    std::cout << "  PASSED" << std::endl;

    // Test 2: Round-trip from_bigint -> to_bigint
    std::cout << "\nTest 2: Round-trip from_bigint -> to_bigint" << std::endl;
    {
        // Test 0
        auto zero_elem = ring.from_bigint(BigInt(0));
        assert(ring.to_bigint(zero_elem) == 0);

        // Test 1
        auto one_elem = ring.from_bigint(BigInt(1));
        assert(ring.to_bigint(one_elem) == 1);

        // Test r-1
        BigInt r_minus_1 = r - 1;
        auto rmin1_elem = ring.from_bigint(r_minus_1);
        assert(ring.to_bigint(rmin1_elem) == r_minus_1);

        // Random values
        for (int i = 0; i < 100; i++) {
            BigInt val = BigInt::random(r);
            auto elem = ring.from_bigint(val);
            BigInt recovered = ring.to_bigint(elem);
            if (recovered != val) {
                std::cerr << "Round-trip failed for value: " << val.to_string(16) << std::endl;
                std::cerr << "Recovered: " << recovered.to_string(16) << std::endl;
                assert(false);
            }
        }
        std::cout << "  PASSED" << std::endl;
    }

    // Test 3: Addition (use check_bounds to extract value from non-StandardElement)
    std::cout << "\nTest 3: Addition (a + b mod r)" << std::endl;
    {
        for (int i = 0; i < 100; i++) {
            BigInt a_val = BigInt::random(r);
            BigInt b_val = BigInt::random(r);

            auto a_elem = ring.from_bigint(a_val);
            auto b_elem = ring.from_bigint(b_val);

            auto sum = a_elem + b_elem;
            auto [result, ok] = ring.check_bounds(sum, "add");
            assert(ok);
            BigInt expected = (a_val + b_val) % r;

            if (result != expected) {
                std::cerr << "Addition failed:" << std::endl;
                std::cerr << "  a = " << a_val.to_string(16) << std::endl;
                std::cerr << "  b = " << b_val.to_string(16) << std::endl;
                std::cerr << "  got = " << result.to_string(16) << std::endl;
                std::cerr << "  expected = " << expected.to_string(16) << std::endl;
                assert(false);
            }
        }
        std::cout << "  PASSED" << std::endl;
    }

    // Test 4: Multiplication
    std::cout << "\nTest 4: Multiplication (a * b mod r)" << std::endl;
    {
        for (int i = 0; i < 100; i++) {
            BigInt a_val = BigInt::random(r);
            BigInt b_val = BigInt::random(r);

            auto a_elem = ring.from_bigint(a_val);
            auto b_elem = ring.from_bigint(b_val);

            auto product = ring.modmul(a_elem, b_elem);
            BigInt result = ring.to_bigint(product);
            BigInt expected = (a_val * b_val) % r;

            if (result != expected) {
                std::cerr << "Multiplication failed:" << std::endl;
                std::cerr << "  a = " << a_val.to_string(16) << std::endl;
                std::cerr << "  b = " << b_val.to_string(16) << std::endl;
                std::cerr << "  got = " << result.to_string(16) << std::endl;
                std::cerr << "  expected = " << expected.to_string(16) << std::endl;
                assert(false);
            }
        }
        std::cout << "  PASSED" << std::endl;
    }

    // Test 5: Negation
    std::cout << "\nTest 5: Negation (r - a mod r)" << std::endl;
    {
        for (int i = 0; i < 100; i++) {
            BigInt a_val = BigInt::random(r);
            auto a_elem = ring.from_bigint(a_val);

            auto neg_elem = ring.standard_negate(a_elem);
            BigInt result = ring.to_bigint(neg_elem);
            BigInt expected = (a_val == 0) ? BigInt(0) : (r - a_val);

            if (result != expected) {
                std::cerr << "Negation failed:" << std::endl;
                std::cerr << "  a = " << a_val.to_string(16) << std::endl;
                std::cerr << "  got = " << result.to_string(16) << std::endl;
                std::cerr << "  expected = " << expected.to_string(16) << std::endl;
                assert(false);
            }
        }
        std::cout << "  PASSED" << std::endl;
    }

    // Test 6: one() identity for multiplication
    std::cout << "\nTest 6: Multiplicative identity (a * 1 = a)" << std::endl;
    {
        auto one = ring.one();
        for (int i = 0; i < 20; i++) {
            BigInt a_val = BigInt::random(r);
            auto a_elem = ring.from_bigint(a_val);
            auto product = ring.modmul(a_elem, one);
            assert(ring.to_bigint(product) == a_val);
        }
        std::cout << "  PASSED" << std::endl;
    }

    // Test 7: Edge cases near modulus boundary
    std::cout << "\nTest 7: Edge cases near modulus boundary" << std::endl;
    {
        // r-1 + 1 = 0 mod r (via check_bounds)
        auto rmin1 = ring.from_bigint(r - 1);
        auto one = ring.from_bigint(BigInt(1));
        auto sum = rmin1 + one;
        auto [sum_val, sum_ok] = ring.check_bounds(sum, "edge_add");
        assert(sum_ok);
        assert(sum_val == 0);

        // (r-1) * (r-1) = 1 mod r
        auto product = ring.modmul(rmin1, rmin1);
        BigInt result = ring.to_bigint(product);
        assert(result == 1);

        // 0 * x = 0
        auto zero = ring.from_bigint(BigInt(0));
        BigInt x_val = BigInt::random(r);
        auto x = ring.from_bigint(x_val);
        auto zero_product = ring.modmul(zero, x);
        assert(ring.to_bigint(zero_product) == 0);

        std::cout << "  PASSED" << std::endl;
    }

    std::cout << "\n=== All Fr tests passed! ===" << std::endl;
    return 0;
}
