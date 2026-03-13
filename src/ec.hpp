#pragma once
#include "preprocess.hpp"
#include <type_traits>

template<class RingElement>
class ProjectivePoint {
    public:
    RingElement x;
    RingElement y;
    RingElement z;
    
    ProjectivePoint() = default;
    ProjectivePoint(const RingElement &x_val, const RingElement &y_val, const RingElement &z_val)
        : x(x_val), y(y_val), z(z_val) {}

    INLINE ProjectivePoint ctime_select(bool selection, const ProjectivePoint &other) const {
        return ProjectivePoint(x.ctime_select(selection, other.x), y.ctime_select(selection, other.y), z.ctime_select(selection, other.z));
    }

    // Clearer API: replace self with replacement if condition is true
    // ctime_select(selection, b) returns *this when selection is true, b when false
    // So we need to negate: if condition is true, we want replacement (not *this)
    INLINE ProjectivePoint ctime_replace_if(bool condition, const ProjectivePoint &replacement) const {
        return ctime_select(!condition, replacement);
    }

    // Clearer API: replace self with replacement if condition is false
    INLINE ProjectivePoint ctime_replace_unless(bool condition, const ProjectivePoint &replacement) const {
        return ctime_select(condition, replacement);
    }
};

template<class RingElement>
class AffinePoint {
    public:
    RingElement x;
    RingElement y;

    AffinePoint() = default;
    AffinePoint(const RingElement &x_, const RingElement &y_) : x(x_), y(y_) {}

    INLINE AffinePoint ctime_select(bool selection, const AffinePoint &other) const {
        return AffinePoint(x.ctime_select(selection, other.x), y.ctime_select(selection, other.y));
    }

    // Clearer API: replace self with replacement if condition is true
    // ctime_select(selection, b) returns *this when selection is true, b when false
    // So we need to negate: if condition is true, we want replacement (not *this)
    INLINE AffinePoint ctime_replace_if(bool condition, const AffinePoint &replacement) const {
        return ctime_select(!condition, replacement);
    }

    // Clearer API: replace self with replacement if condition is false
    INLINE AffinePoint ctime_replace_unless(bool condition, const AffinePoint &replacement) const {
        return ctime_select(condition, replacement);
    }
};

// XYZZ coordinates: (X, Y, ZZ, ZZZ) where ZZ = Z², ZZZ = Z³.
// Represents affine point (X/ZZ, Y/ZZZ). Identity when ZZ == 0 && ZZZ == 0.
template<class RingElement>
class XYZZPoint {
    public:
    RingElement X;
    RingElement Y;
    RingElement ZZ;
    RingElement ZZZ;

    XYZZPoint() = default;
    XYZZPoint(const RingElement &x, const RingElement &y, const RingElement &zz, const RingElement &zzz)
        : X(x), Y(y), ZZ(zz), ZZZ(zzz) {}
};

// XYZZ + Affine addition.
// P = (X1, Y1, ZZ1, ZZZ1), Q = (x2, y2) affine.
// 11 multiplications (including 1 identity for X3 normalization), 4 batch rounds.
// No mul_3b needed (unlike projective Algorithm 7).
template<class Ring>
XYZZPoint<typename Ring::StandardElement> XYZZAddAffine(
    const XYZZPoint<typename Ring::StandardElement> &P,
    const AffinePoint<typename Ring::StandardElement> &Q,
    const Ring &ring)
{
    // Round 1 (2M): U2 = x2*ZZ1, S2 = y2*ZZZ1
    auto U2_wide = ring.prep_left(Q.x) * P.ZZ;
    auto S2_wide = ring.prep_left(Q.y) * P.ZZZ;
    auto [U2, S2] = ring.batch_reduce_expand(U2_wide, S2_wide);

    auto H = U2 - P.X;
    auto R = S2 - P.Y;

    // Round 2 (2M): H², R²
    auto H_p = ring.prep(H);
    auto R_p = ring.prep(R);
    auto HH_wide = H_p * H_p;
    auto Rsq_wide = R_p * R_p;
    auto [HH, Rsq] = ring.batch_reduce_expand(HH_wide, Rsq_wide);

    // Round 3 (3M): H*HH, X1*HH, ZZ1*HH
    auto HH_p = ring.prep(HH);
    auto HHH_wide = H_p * HH_p;
    auto V_wide = ring.prep_left(P.X) * HH_p;
    auto ZZ3_wide = ring.prep_left(P.ZZ) * HH_p;
    auto [HHH, V, ZZ3] = ring.batch_reduce_expand(HHH_wide, V_wide, ZZ3_wide);

    // X3 = Rsq - HHH - 2V (non-standard bounds, normalized via identity mul)
    auto X3_sub = Rsq - HHH - V - V;
    auto V_minus_X3 = V - X3_sub;

    // Round 4 (4M): X3 normalization, Y3 (accumulated), ZZZ3
    auto one = ring.one();
    auto X3_wide = ring.prep_left(X3_sub) * ring.prep(one);
    auto HHH_p2 = ring.prep(HHH);
    auto R_vmx_wide = ring.prep_left(R) * ring.prep(V_minus_X3);
    auto Y1_negHHH_wide = ring.prep_left(P.Y) * ring.negate(HHH_p2);
    auto Y3_wide = R_vmx_wide + Y1_negHHH_wide;
    auto ZZZ3_wide = ring.prep_left(P.ZZZ) * HHH_p2;

    auto [X3, Y3, ZZZ3] = ring.batch_reduce_expand(X3_wide, Y3_wide, ZZZ3_wide);
    return XYZZPoint<typename Ring::StandardElement>(X3, Y3, ZZ3, ZZZ3);
}

// XYZZ + XYZZ addition for integration phase.
// P = (X1,Y1,ZZ1,ZZZ1), Q = (X2,Y2,ZZ2,ZZZ2).
// 15 multiplications (including 1 identity for X3), 4 batch rounds.
template<class Ring>
XYZZPoint<typename Ring::StandardElement> XYZZAdd(
    const XYZZPoint<typename Ring::StandardElement> &P,
    const XYZZPoint<typename Ring::StandardElement> &Q,
    const Ring &ring)
{
    // Round 1 (6M): cross products
    auto U1_wide = ring.prep_left(P.X) * Q.ZZ;
    auto S1_wide = ring.prep_left(P.Y) * Q.ZZZ;
    auto U2_wide = ring.prep_left(Q.X) * P.ZZ;
    auto S2_wide = ring.prep_left(Q.Y) * P.ZZZ;
    auto ZZcross_wide = ring.prep_left(P.ZZ) * Q.ZZ;
    auto ZZZcross_wide = ring.prep_left(P.ZZZ) * Q.ZZZ;
    auto [U1, S1, U2, S2, ZZcross, ZZZcross] = ring.batch_reduce_expand(
        U1_wide, S1_wide, U2_wide, S2_wide, ZZcross_wide, ZZZcross_wide);

    auto H = U2 - U1;
    auto R = S2 - S1;

    // Round 2 (2M): H², R²
    auto H_p = ring.prep(H);
    auto R_p = ring.prep(R);
    auto HH_wide = H_p * H_p;
    auto Rsq_wide = R_p * R_p;
    auto [HH, Rsq] = ring.batch_reduce_expand(HH_wide, Rsq_wide);

    // Round 3 (4M): H*HH, U1*HH, ZZcross*HH, ZZZcross*ZZZ deferred
    auto HH_p = ring.prep(HH);
    auto HHH_wide = H_p * HH_p;
    auto V_wide = ring.prep_left(U1) * HH_p;
    auto ZZ3_wide = ring.prep_left(ZZcross) * HH_p;
    auto [HHH, V, ZZ3] = ring.batch_reduce_expand(HHH_wide, V_wide, ZZ3_wide);

    // X3 = Rsq - HHH - 2V
    auto X3_sub = Rsq - HHH - V - V;
    auto V_minus_X3 = V - X3_sub;

    // Round 4 (5M): X3 normalization, Y3 (accumulated), ZZZ3
    auto one = ring.one();
    auto X3_wide = ring.prep_left(X3_sub) * ring.prep(one);
    auto HHH_p2 = ring.prep(HHH);
    auto R_vmx_wide = ring.prep_left(R) * ring.prep(V_minus_X3);
    auto S1_negHHH_wide = ring.prep_left(S1) * ring.negate(HHH_p2);
    auto Y3_wide = R_vmx_wide + S1_negHHH_wide;
    auto ZZZ3_wide = ring.prep_left(ZZZcross) * HHH_p2;

    auto [X3, Y3, ZZZ3] = ring.batch_reduce_expand(X3_wide, Y3_wide, ZZZ3_wide);
    return XYZZPoint<typename Ring::StandardElement>(X3, Y3, ZZ3, ZZZ3);
}

// Convert XYZZ → Projective.
// (X, Y, ZZ, ZZZ) → (X*ZZZ, Y*ZZ, ZZ*ZZZ)
// since affine = (X/ZZ, Y/ZZZ) and projective = (X_p/Z_p, Y_p/Z_p).
template<class Ring>
ProjectivePoint<typename Ring::StandardElement> XYZZToProj(
    const XYZZPoint<typename Ring::StandardElement> &P,
    const Ring &ring)
{
    auto Xp_wide = ring.prep_left(P.X) * P.ZZZ;
    auto Yp_wide = ring.prep_left(P.Y) * P.ZZ;
    auto Zp_wide = ring.prep_left(P.ZZ) * P.ZZZ;
    auto [Xp, Yp, Zp] = ring.batch_reduce_expand(Xp_wide, Yp_wide, Zp_wide);
    return ProjectivePoint<typename Ring::StandardElement>(Xp, Yp, Zp);
}

// ~5% faster, but doesn't work with new type system yet because it needs to preserve rns bounds across expansion.
/*
template<class Ring, class OffsetElement>
ProjectivePoint<Ring::StandardElement> PointAdd(const ProjectivePoint<Ring::StandardElement> &P, const ProjectivePoint<Ring::StandardElement> &Q, const Ring &ring, const OffsetElement &offset) {
    auto t3_4 = P.x + P.y;
    auto t3_14 = P.x + P.z;
    auto t4_9 = P.y + P.z;
    auto t4_5 = Q.x + Q.y;
    auto x3_10 = Q.y + Q.z;
    auto y3_15 = Q.x + Q.z;

    auto t0_1 = ring.prep_left(P.x) * Q.x;
    auto t1_2 = ring.prep_left(P.y) * Q.y;
    auto t2_3 = ring.prep_left(P.z) * Q.z;
    auto t3_6 = ring.prep_left(t3_4) * ring.prep(t4_5);
    auto t5_11 = ring.prep_left(t4_9) * ring.prep(x3_10);
    auto x3_16 = ring.prep_left(t3_14) * ring.prep(y3_15);

    auto [t0_1_reduced, t1_2_reduced, t2_3_reduced, t3_6_reduced, t5_11_reduced, x3_16_reduced] = ring.batch_reduce(t0_1, t1_2, t2_3, t3_6, t5_11, x3_16);
    auto t4_7 = t0_1_reduced + t1_2_reduced;
    auto x3_12 = t1_2_reduced + t2_3_reduced;
    auto y3_17 = t0_1_reduced + t2_3_reduced;
    auto y3_18 = x3_16_reduced - y3_17;
    auto x3_19 = t0_1_reduced + t0_1_reduced;
    auto t2_21 = ring.mul_3b(t2_3_reduced);

    auto t3_8_s = t3_6_reduced - t4_7;
    auto t4_13_s = t5_11_reduced - x3_12;
    auto t0_20_s = t0_1_reduced + x3_19;
    auto z3_22_s = t1_2_reduced + t2_21;
    auto t1_23_s = t1_2_reduced - t2_21;
    auto y3_24_s = ring.mul_3b(y3_18);

    auto [t4_13, t3_8, y3_24, t1_23, t0_20, z3_22] = ring.batch_expand(t4_13_s, t3_8_s, y3_24_s, t1_23_s, t0_20_s, z3_22_s);
    
    auto t4_13l = ring.prep_left(t4_13);
    auto t1_23l = ring.prep_left(t1_23);
    auto t0_20l = ring.prep_left(t0_20);
    auto x3_25 = t4_13l * ring.negate(y3_24);
    auto t2_26 = t1_23l * t3_8;
    auto y3_28 = t0_20l * y3_24;
    auto t1_29 = t1_23l * z3_22;
    auto t0_31 = t0_20l * t3_8;
    auto z3_32 = t4_13l * z3_22;
    auto x3_27_wide = offset + t2_26 + x3_25;
    auto y3_30_wide = offset + t1_29 + y3_28;
    auto z3_33_wide = offset + z3_32 + t0_31;
    
    auto [x3_27, y3_30, z3_33] = ring.batch_reduce_expand(x3_27_wide, y3_30_wide, z3_33_wide);
    return ProjectivePoint<Ring::StandardElement>(x3_27, y3_30, z3_33);
}
*/

// https://eprint.iacr.org/2015/1060.pdf Algorithm 7
template<class Ring>
ProjectivePoint<typename Ring::StandardElement> PointAdd(const ProjectivePoint<typename Ring::StandardElement> &P, const ProjectivePoint<typename Ring::StandardElement> &Q, const Ring &ring) {
    auto t3_4 = P.x + P.y;
    auto t3_14 = P.x + P.z;
    auto t4_9 = P.y + P.z;
    auto t4_5 = Q.x + Q.y;
    auto x3_10 = Q.y + Q.z;
    auto y3_15 = Q.x + Q.z;

    auto t0_1_wide = ring.prep_left(P.x) * Q.x;
    auto t1_2_wide = ring.prep_left(P.y) * Q.y;
    auto t2_3_wide = ring.prep_left(P.z) * Q.z;
    auto t3_6_wide = ring.prep_left(t3_4) * ring.prep(t4_5);
    auto t5_11_wide = ring.prep_left(t4_9) * ring.prep(x3_10);
    auto x3_16_wide = ring.prep_left(t3_14) * ring.prep(y3_15);

    // The simplest solution is to use tab completion or scripting to write out the cases for 1..12 inputs
    // 
    // expand(SmallElement<BoundedElement<Bounds<l, u>>, )
    // ring.prep() -> mont 2 auto
    // batch<> with avx limbs as inputs
    // reattach rns bounds.
    auto [t0_1, t1_2, t2_3, t3_6, t5_11, x3_16] = ring.batch_reduce(t0_1_wide, t1_2_wide, t2_3_wide, t3_6_wide, t5_11_wide, x3_16_wide);
    auto t4_7 = t0_1 + t1_2;
    auto x3_12 = t1_2 + t2_3;
    auto y3_17 = t0_1 + t2_3;
    auto y3_18 = x3_16 - y3_17;
    auto x3_19 = t0_1 + t0_1;
    auto t2_21 = ring.mul_3b(t2_3);

    auto t3_8_s = ring.prep(t3_6 - t4_7);
    auto y3_24_s = ring.prep(ring.mul_3b(y3_18));
    auto z3_22_s = ring.prep(t1_2 + t2_21);
    auto t1_23_s = ring.prep(t1_2 - t2_21);
    auto t4_13_s = ring.prep(t5_11 - x3_12);
    auto t0_20_s = ring.prep(t0_1 + x3_19);

    auto [t3_8, y3_24, z3_22, t1_23, t4_13, t0_20] = ring.batch_expand(t3_8_s, y3_24_s, z3_22_s, t1_23_s, t4_13_s, t0_20_s);

    auto t4_13l = ring.prep_left(t4_13);
    auto t1_23l = ring.prep_left(t1_23);
    auto t0_20l = ring.prep_left(t0_20);
    
    auto x3_25 = t4_13l * ring.negate(y3_24);
    auto t2_26 = t1_23l * t3_8;
    auto y3_28 = t0_20l * y3_24;
    auto t1_29 = t1_23l * z3_22;
    auto t0_31 = t0_20l * t3_8;
    auto z3_32 = t4_13l * z3_22;
    auto x3_27_wide =  t2_26 + x3_25;
    auto y3_30_wide =  t1_29 + y3_28;
    auto z3_33_wide =  z3_32 + t0_31;
    
    auto [x3_27, y3_30, z3_33] = ring.batch_reduce_expand(x3_27_wide, y3_30_wide, z3_33_wide);
    return ProjectivePoint<typename Ring::StandardElement>(x3_27, y3_30, z3_33);
}

// Not enough time to figure out Algorithm 8; this saves a mult and you could save more additions (by also using small form)
template<class Ring>
ProjectivePoint<typename Ring::StandardElement> PointMixedAdd(const ProjectivePoint<typename Ring::StandardElement> &P, const AffinePoint<typename Ring::StandardElement> &Q, const Ring &ring) {
    auto t3_4 = P.x + P.y;
    auto t3_14 = P.x + P.z;
    auto t4_9 = P.y + P.z;
    auto t4_5 = Q.x + Q.y;
    auto one = ring.one();
    auto x3_10 = Q.y + one;
    auto y3_15 = Q.x + one;

    auto t0_1_wide = ring.prep_left(P.x) * Q.x;
    auto t1_2_wide = ring.prep_left(P.y) * Q.y;
    // auto t2_3_wide = ring.prep_left(P.z) * Q.z;
    auto t3_6_wide = ring.prep_left(t3_4) * ring.prep(t4_5);
    auto t5_11_wide = ring.prep_left(t4_9) * ring.prep(x3_10);
    auto x3_16_wide = ring.prep_left(t3_14) * ring.prep(y3_15);

    auto [t0_1, t1_2, t3_6, t5_11, x3_16] = ring.batch_reduce_expand(t0_1_wide, t1_2_wide, t3_6_wide, t5_11_wide, x3_16_wide);
    auto t2_3 = P.z;
    auto t4_7 = t0_1 + t1_2;
    auto x3_12 = t1_2 + t2_3;
    auto y3_17 = t0_1 + t2_3;
    auto y3_18 = x3_16 - y3_17;
    auto x3_19 = t0_1 + t0_1;
    auto t2_21 = ring.prep(ring.mul_3b(t2_3));

    auto t3_8 = ring.prep(t3_6 - t4_7);
    auto y3_24 = ring.prep(ring.mul_3b(y3_18));
    auto z3_22 = ring.prep(t1_2 + t2_21);
    auto t1_23 = ring.prep_left(t1_2 - t2_21);
    auto t4_13 = ring.prep_left(t5_11 - x3_12);
    auto t0_20 = ring.prep_left(t0_1 + x3_19);
    
    auto x3_25 = t4_13 * ring.negate(y3_24);
    auto t2_26 = t1_23 * t3_8;
    auto y3_28 = t0_20 * y3_24;
    auto t1_29 = t1_23 * z3_22;
    auto t0_31 = t0_20 * t3_8;
    auto z3_32 = t4_13 * z3_22;
    auto x3_27_wide =  t2_26 + x3_25;
    auto y3_30_wide =  t1_29 + y3_28;
    auto z3_33_wide =  z3_32 + t0_31;
    
    auto [x3_27, y3_30, z3_33] = ring.batch_reduce_expand(x3_27_wide, y3_30_wide, z3_33_wide);
    return ProjectivePoint<typename Ring::StandardElement>(x3_27, y3_30, z3_33);
}

// https://eprint.iacr.org/2015/1060.pdf Algorithm 9
template<class Ring>
ProjectivePoint<typename Ring::StandardElement> PointDouble(const ProjectivePoint<typename Ring::StandardElement> &P,const Ring &ring) {
    auto py = ring.prep_left(P.y);
    auto pz = ring.prep_left(P.z);
    constexpr auto _2 = std::integral_constant<int, 2>{};
    constexpr auto _3 = std::integral_constant<int, 3>{};
    constexpr auto _8 = std::integral_constant<int, 8>{};

    auto yy = py * P.y;
    auto zz = pz * P.z;
    auto xy = py * P.x;
    auto yz = py * P.z;

    auto [yy_reduced, zz_reduced, xy_reduced, yz_reduced] = ring.batch_reduce_expand(yy, zz, xy, yz);

    // Add/sub operations
    auto _8yy = ring.prep_left(yy_reduced * _8);
    auto b3_zz = ring.prep(ring.mul_3b(zz_reduced));
    auto yy_plus_b3_zz = ring.prep(yy_reduced + b3_zz);
    auto b9_zz = b3_zz * _3;
    auto yy_minus_b9_zz = ring.prep_left(yy_reduced - b9_zz);
    auto two_xy = ring.prep(xy_reduced * _2);

    // Multiplication round 2
    auto y3_14 = yy_minus_b9_zz * yy_plus_b3_zz;
    auto x3_18 = yy_minus_b9_zz * two_xy;
    auto z3_10 = _8yy * yz_reduced;
    auto y3_15 = y3_14 + (_8yy * b3_zz);

    auto [x3, y3, z3] = ring.batch_reduce_expand(x3_18, y3_15, z3_10);
    return ProjectivePoint<typename Ring::StandardElement>(x3, y3, z3);
}


template<class Ring>
class G1 {
    public:
    using ProjPoint = ProjectivePoint<typename Ring::StandardElement>;
    using AffPoint = AffinePoint<typename Ring::StandardElement>;
    using XYZZPt = XYZZPoint<typename Ring::StandardElement>;

    ProjPoint add_point(const ProjPoint &P, const ProjPoint &Q, const Ring &ring) const {
        return PointAdd<Ring>(P, Q, ring);
    }

    ProjPoint add_mixed_point(const ProjPoint &P, const AffPoint &Q, const Ring &ring) const {
        return PointMixedAdd<Ring>(P, Q, ring);
    }

    ProjPoint double_point(const ProjPoint &P, const Ring &ring) const {
        return PointDouble<Ring>(P, ring);
    }

    INLINE ProjPoint negate(const ProjPoint &P, const Ring &ring) const {
        return ProjPoint(P.x, ring.standard_negate(P.y), P.z);
    }

    INLINE AffPoint negate_affine(const AffPoint &P, const Ring &ring) const {
        return AffPoint(P.x, ring.standard_negate(P.y));
    }

    INLINE ProjPoint zero(const Ring &ring) const {
        return ProjPoint(ring.zero(), ring.one(), ring.zero());
    }

    // XYZZ operations for Pippenger MSM
    INLINE XYZZPt xyzz_from_affine(const AffPoint &P, const Ring &ring) const {
        return XYZZPt(P.x, P.y, ring.one(), ring.one());
    }

    INLINE XYZZPt xyzz_zero(const Ring &ring) const {
        return XYZZPt(ring.zero(), ring.one(), ring.zero(), ring.zero());
    }

    XYZZPt xyzz_add_affine(const XYZZPt &P, const AffPoint &Q, const Ring &ring) const {
        return XYZZAddAffine<Ring>(P, Q, ring);
    }

    XYZZPt xyzz_add(const XYZZPt &P, const XYZZPt &Q, const Ring &ring) const {
        return XYZZAdd<Ring>(P, Q, ring);
    }

    ProjPoint xyzz_to_proj(const XYZZPt &P, const Ring &ring) const {
        return XYZZToProj<Ring>(P, ring);
    }

};

// Forward declaration
template<class Ring> class FP2Ring;

template<class Ring>
class G2 {
    public:
    using FP2RingType = FP2Ring<Ring>;
    using ProjPoint = ProjectivePoint<typename FP2RingType::StandardElement>;
    using AffPoint = AffinePoint<typename FP2RingType::StandardElement>;

    ProjPoint add_point(const ProjPoint &P, const ProjPoint &Q, const FP2RingType &ring) const {
        return PointAdd<FP2RingType>(P, Q, ring);
    }

    ProjPoint add_mixed_point(const ProjPoint &P, const AffPoint &Q, const FP2RingType &ring) const {
        return PointMixedAdd<FP2RingType>(P, Q, ring);
    }

    ProjPoint double_point(const ProjPoint &P, const FP2RingType &ring) const {
        return PointDouble<FP2RingType>(P, ring);
    }

    INLINE ProjPoint negate(const ProjPoint &P, const FP2RingType &ring) const {
        return ProjPoint(P.x, ring.standard_negate(P.y), P.z);
    }

    INLINE ProjPoint zero(const FP2RingType &ring) const {
        return ProjPoint(ring.zero(), ring.one(), ring.zero());
    }
};