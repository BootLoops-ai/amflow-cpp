// SPDX-License-Identifier: MIT
// Tests for amflow::numeric::AcbValue / AcbVector / AcbMatrix / AcbPoly.
//
// Domain contract: docs/REFACTOR_DESIGN.md §3.1, §5.1.

#include <gtest/gtest.h>

#include <flint/acb.h>
#include <flint/arb.h>
#include <flint/fmpq.h>
#include <flint/fmpz.h>

#include "amflow/numeric/acb_wrap.hpp"
#include "amflow/numeric/options.hpp"

#include <string>
#include <utility>

namespace dsn = amflow::numeric;

// ---------------------------------------------------------------------------
//  AcbValue
// ---------------------------------------------------------------------------

TEST(AcbValueLifecycle, DefaultIsZeroAndFinite) {
    dsn::AcbValue v;
    EXPECT_TRUE(v.is_zero());
    EXPECT_TRUE(v.is_finite());
    EXPECT_FALSE(v.is_one());
    EXPECT_TRUE(v.is_real());
}

TEST(AcbValueLifecycle, SetOneRoundTrips) {
    dsn::AcbValue v;
    v.set_one();
    EXPECT_TRUE(v.is_one());
    EXPECT_FALSE(v.is_zero());
    EXPECT_TRUE(v.is_real());
}

TEST(AcbValueLifecycle, SetSiRoundTrips) {
    dsn::AcbValue v;
    v.set_si(42);
    EXPECT_FALSE(v.is_zero());
    EXPECT_TRUE(v.is_real());
    // Compare via FLINT to avoid string fragility.
    fmpq_t expected;
    fmpq_init(expected);
    fmpq_set_si(expected, 42, 1);
    dsn::AcbValue tmp;
    tmp.set_fmpq(expected);
    EXPECT_NE(acb_eq(v.raw(), tmp.raw()), 0);
    fmpq_clear(expected);
}

TEST(AcbValueLifecycle, SetSiSiHasImaginaryPart) {
    dsn::AcbValue v;
    v.set_si_si(3, 5);
    EXPECT_FALSE(v.is_real());
    EXPECT_FALSE(v.is_zero());
}

TEST(AcbValueLifecycle, SetFmpqRespectsRationals) {
    fmpq_t q;
    fmpq_init(q);
    fmpq_set_si(q, 7, 3);

    dsn::AcbValue v;
    v.set_fmpq(q);

    // Multiply by 3 to clear the fraction; should equal 7.
    acb_t three, prod, seven;
    acb_init(three); acb_init(prod); acb_init(seven);
    acb_set_si(three, 3);
    acb_mul(prod, v.raw(), three, dsn::working_prec_bits());
    acb_set_si(seven, 7);
    EXPECT_NE(acb_overlaps(prod, seven), 0);   // overlap (both should be exact 7)
    acb_clear(three); acb_clear(prod); acb_clear(seven);
    fmpq_clear(q);
}

TEST(AcbValueLifecycle, MoveSemantics) {
    dsn::AcbValue a;
    a.set_si(11);

    dsn::AcbValue b = std::move(a);
    EXPECT_TRUE(b.is_finite());
    // Verify b == 11 by re-setting another and comparing.
    dsn::AcbValue eleven;
    eleven.set_si(11);
    EXPECT_NE(acb_eq(b.raw(), eleven.raw()), 0);

    dsn::AcbValue c;
    c = std::move(b);
    EXPECT_NE(acb_eq(c.raw(), eleven.raw()), 0);
}

TEST(AcbValueLifecycle, CloneIsDeepCopy) {
    dsn::AcbValue a;
    a.set_si(7);
    dsn::AcbValue b = a.clone();
    // Mutate a; b stays at 7.
    a.set_si(999);
    dsn::AcbValue seven;
    seven.set_si(7);
    EXPECT_NE(acb_eq(b.raw(), seven.raw()), 0);
    dsn::AcbValue nine_nine_nine;
    nine_nine_nine.set_si(999);
    EXPECT_NE(acb_eq(a.raw(), nine_nine_nine.raw()), 0);
}

TEST(AcbValueChop, ZeroIsChoppedRegardlessOfDigits) {
    dsn::AcbValue v;     // default 0
    EXPECT_TRUE(v.is_chop_zero(20));
    EXPECT_TRUE(v.is_chop_zero(50));
}

TEST(AcbValueChop, OneIsNotChopped) {
    dsn::AcbValue v;
    v.set_one();
    EXPECT_FALSE(v.is_chop_zero(20));
}

TEST(AcbValueChop, SmallNonZeroChoppedAtCoarseDigits) {
    dsn::AcbValue v;
    // 10^-30 is chop-zero at 20 digits, not at 50 digits.
    arb_set_si(acb_realref(v.raw()), 1);
    arb_t denom; arb_init(denom);
    arb_set_si(denom, 10);
    arb_pow_ui(denom, denom, 30u, dsn::working_prec_bits());
    arb_div(acb_realref(v.raw()), acb_realref(v.raw()), denom, dsn::working_prec_bits());
    arb_clear(denom);

    EXPECT_TRUE(v.is_chop_zero(20));
    EXPECT_FALSE(v.is_chop_zero(50));
}

TEST(AcbValueChop, NonPositiveDigitsAlwaysFalse) {
    dsn::AcbValue v;       // zero
    EXPECT_FALSE(v.is_chop_zero(0));
    EXPECT_FALSE(v.is_chop_zero(-3));
}

TEST(AcbValueString, ToStringContainsImaginaryUnit) {
    dsn::AcbValue v;
    v.set_si_si(2, 3);
    auto s = v.to_string(10);
    EXPECT_NE(s.find("*I"), std::string::npos);
}

// ---------------------------------------------------------------------------
//  AcbVector
// ---------------------------------------------------------------------------

TEST(AcbVectorLifecycle, EmptyConstruction) {
    dsn::AcbVector v;
    EXPECT_TRUE(v.empty());
    EXPECT_EQ(v.size(), 0u);
    EXPECT_EQ(v.raw(), nullptr);
}

TEST(AcbVectorLifecycle, SizedConstructionAndDefaultZero) {
    dsn::AcbVector v(5);
    EXPECT_EQ(v.size(), 5u);
    for (std::size_t i = 0; i < 5; ++i) {
        EXPECT_NE(acb_is_zero(v.at(i)), 0);
    }
}

TEST(AcbVectorLifecycle, MoveTransfersOwnership) {
    dsn::AcbVector a(3);
    acb_set_si(a.at(1), 17);
    dsn::AcbVector b = std::move(a);
    EXPECT_EQ(b.size(), 3u);
    EXPECT_EQ(a.size(), 0u);

    acb_t seventeen;
    acb_init(seventeen);
    acb_set_si(seventeen, 17);
    EXPECT_NE(acb_eq(b.at(1), seventeen), 0);
    acb_clear(seventeen);
}

TEST(AcbVectorLifecycle, CloneIsDeepCopy) {
    dsn::AcbVector a(2);
    acb_set_si(a.at(0), 4);
    acb_set_si(a.at(1), 5);

    dsn::AcbVector b = a.clone();
    acb_set_si(a.at(0), 999);

    acb_t four;
    acb_init(four);
    acb_set_si(four, 4);
    EXPECT_NE(acb_eq(b.at(0), four), 0);
    acb_clear(four);
}

// ---------------------------------------------------------------------------
//  AcbMatrix
// ---------------------------------------------------------------------------

TEST(AcbMatrixLifecycle, EmptyConstruction) {
    dsn::AcbMatrix m;
    EXPECT_TRUE(m.empty());
    EXPECT_EQ(m.rows(), 0);
    EXPECT_EQ(m.cols(), 0);
}

TEST(AcbMatrixLifecycle, SizedConstructionDefaultsToZero) {
    dsn::AcbMatrix m(3, 4);
    EXPECT_EQ(m.rows(), 3);
    EXPECT_EQ(m.cols(), 4);
    EXPECT_FALSE(m.empty());
    for (long i = 0; i < 3; ++i) {
        for (long j = 0; j < 4; ++j) {
            EXPECT_NE(acb_is_zero(m.entry(i, j)), 0);
        }
    }
}

TEST(AcbMatrixLifecycle, SetIdentity) {
    dsn::AcbMatrix m(3, 3);
    m.set_identity();
    for (long i = 0; i < 3; ++i) {
        for (long j = 0; j < 3; ++j) {
            if (i == j) EXPECT_NE(acb_is_one(m.entry(i, j)), 0);
            else        EXPECT_NE(acb_is_zero(m.entry(i, j)), 0);
        }
    }
}

TEST(AcbMatrixLifecycle, MoveTransfersOwnership) {
    dsn::AcbMatrix a(2, 2);
    acb_set_si(a.entry(0, 1), 9);
    dsn::AcbMatrix b = std::move(a);
    EXPECT_EQ(b.rows(), 2);
    EXPECT_EQ(b.cols(), 2);
    EXPECT_TRUE(a.empty());

    acb_t nine;
    acb_init(nine);
    acb_set_si(nine, 9);
    EXPECT_NE(acb_eq(b.entry(0, 1), nine), 0);
    acb_clear(nine);
}

TEST(AcbMatrixLifecycle, CloneIsDeepCopy) {
    dsn::AcbMatrix a(2, 2);
    acb_set_si(a.entry(0, 0), 7);
    dsn::AcbMatrix b = a.clone();
    acb_set_si(a.entry(0, 0), 999);

    acb_t seven;
    acb_init(seven);
    acb_set_si(seven, 7);
    EXPECT_NE(acb_eq(b.entry(0, 0), seven), 0);
    acb_clear(seven);
}

// ---------------------------------------------------------------------------
//  AcbPoly
// ---------------------------------------------------------------------------

TEST(AcbPolyLifecycle, DefaultIsZero) {
    dsn::AcbPoly p;
    EXPECT_TRUE(p.is_zero());
    EXPECT_EQ(p.length(), 0);
    EXPECT_EQ(p.degree(), -1);
}

TEST(AcbPolyLifecycle, MoveSemantics) {
    dsn::AcbPoly a;
    acb_t coeff;
    acb_init(coeff);
    acb_set_si(coeff, 3);
    acb_poly_set_coeff_acb(a.raw(), 2, coeff);
    EXPECT_EQ(a.degree(), 2);

    dsn::AcbPoly b = std::move(a);
    EXPECT_EQ(b.degree(), 2);

    acb_clear(coeff);
}

TEST(AcbPolyLifecycle, CloneIsDeepCopy) {
    dsn::AcbPoly a;
    acb_t coeff;
    acb_init(coeff);
    acb_set_si(coeff, 5);
    acb_poly_set_coeff_acb(a.raw(), 1, coeff);

    dsn::AcbPoly b = a.clone();
    a.set_zero();

    EXPECT_TRUE(a.is_zero());
    EXPECT_FALSE(b.is_zero());
    EXPECT_EQ(b.degree(), 1);

    acb_clear(coeff);
}

// ---------------------------------------------------------------------------
//  Free helper acb_is_chop_zero
// ---------------------------------------------------------------------------

TEST(FreeHelperChop, DelegatesToMember) {
    dsn::AcbValue v;
    v.set_one();
    EXPECT_FALSE(dsn::acb_is_chop_zero(v.raw(), 20));
    v.set_zero();
    EXPECT_TRUE(dsn::acb_is_chop_zero(v.raw(), 20));
}

// ---------------------------------------------------------------------------
//  arb_to_fmpq_simplest / acb_real_to_fmpq_simplest
// ---------------------------------------------------------------------------

namespace {

struct FmpqHolder {
    fmpq_t q;
    FmpqHolder()  { fmpq_init(q); }
    ~FmpqHolder() { fmpq_clear(q); }
};

// A ball around the exact rational p/q with the given radius 10^-rad_digits.
void ball_around(arb_t out, long p, long q, int rad_digits, long prec = 400) {
    FmpqHolder v;
    fmpq_set_si(v.q, p, q);
    arb_set_fmpq(out, v.q, prec);
    if (rad_digits > 0) {
        arb_t r;
        arb_init(r);
        arb_set_si(r, 10);
        arb_pow_ui(r, r, static_cast<unsigned long>(rad_digits), prec);
        arb_inv(r, r, prec);
        arb_add_error(out, r);
        arb_clear(r);
    }
}

std::string fmpq_str(const fmpq_t q) {
    char* s = fmpq_get_str(nullptr, 10, q);
    std::string out(s);
    flint_free(s);
    return out;
}

}  // namespace

TEST(RationalizeSimplest, TightBallAroundOneOver113IsExact) {
    arb_t x;
    arb_init(x);
    ball_around(x, 1, 113, 90);
    FmpqHolder got, want;
    ASSERT_TRUE(dsn::arb_to_fmpq_simplest(got.q, x, 60));
    fmpq_set_si(want.q, 1, 113);
    EXPECT_TRUE(fmpq_equal(got.q, want.q)) << "got " << fmpq_str(got.q);
    arb_clear(x);
}

TEST(RationalizeSimplest, NineFifthsIsExact) {
    arb_t x;
    arb_init(x);
    ball_around(x, 9, 5, 0, 200);
    FmpqHolder got, want;
    ASSERT_TRUE(dsn::arb_to_fmpq_simplest(got.q, x, 30));
    fmpq_set_si(want.q, 9, 5);
    EXPECT_TRUE(fmpq_equal(got.q, want.q)) << "got " << fmpq_str(got.q);
    arb_clear(x);
}

TEST(RationalizeSimplest, PiWithinToleranceAndNotTallerThanTruncation) {
    arb_t pi;
    arb_init(pi);
    arb_const_pi(pi, 200);
    FmpqHolder got;
    ASSERT_TRUE(dsn::arb_to_fmpq_simplest(got.q, pi, 30));

    // |out - pi| < 10^-30, checked in a 200-bit ball.
    arb_t diff, bound;
    arb_init(diff);
    arb_init(bound);
    arb_set_fmpq(diff, got.q, 200);
    arb_sub(diff, diff, pi, 200);
    arb_abs(diff, diff);
    arb_set_si(bound, 10);
    arb_pow_ui(bound, bound, 30, 200);
    arb_inv(bound, bound, 200);
    EXPECT_TRUE(arb_lt(diff, bound)) << "got " << fmpq_str(got.q);

    // Denominator below 10^30 (the truncation's denominator).
    fmpz_t pow10;
    fmpz_init(pow10);
    fmpz_set_ui(pow10, 10);
    fmpz_pow_ui(pow10, pow10, 30);
    EXPECT_LT(fmpz_cmp(fmpq_denref(got.q), pow10), 0) << "got " << fmpq_str(got.q);
    fmpz_clear(pow10);

    arb_clear(diff);
    arb_clear(bound);
    arb_clear(pi);
}

TEST(RationalizeSimplest, WideBallRadiusIsIgnored) {
    // Radius 1e-10 around a midpoint that is 1/3 to 400 bits, at 60 digits:
    // the radius plays no part, the midpoint alone is rationalized, and the
    // simplest rational within 10^-60 of it is 1/3.
    arb_t x;
    arb_init(x);
    ball_around(x, 1, 3, 10);
    FmpqHolder got, want;
    ASSERT_TRUE(dsn::arb_to_fmpq_simplest(got.q, x, 60));
    fmpq_set_si(want.q, 1, 3);
    EXPECT_TRUE(fmpq_equal(got.q, want.q)) << "got " << fmpq_str(got.q);
    arb_clear(x);
}

TEST(RationalizeSimplest, RadiusLargerThanToleranceKeepsMidpointDigits) {
    // Midpoint = pi truncated to 60 decimal digits, radius 1e-25, digits 60:
    // the result is within 10^-60 of the midpoint (the radius, 35 orders of
    // magnitude wider, is ignored).
    fmpz_t pow10, numer;
    fmpz_init(pow10);
    fmpz_init(numer);
    fmpz_set_ui(pow10, 10);
    fmpz_pow_ui(pow10, pow10, 60);
    {
        arb_t pi;
        arb_init(pi);
        arb_const_pi(pi, 400);
        arb_mul_fmpz(pi, pi, pow10, 400);
        arf_get_fmpz(numer, arb_midref(pi), ARF_RND_FLOOR);
        arb_clear(pi);
    }
    FmpqHolder mid;
    fmpq_set_fmpz_frac(mid.q, numer, pow10);

    arb_t x;
    arb_init(x);
    arb_set_fmpq(x, mid.q, 400);
    {
        arb_t r;
        arb_init(r);
        arb_set_si(r, 10);
        arb_pow_ui(r, r, 25, 400);
        arb_inv(r, r, 400);
        arb_add_error(x, r);
        arb_clear(r);
    }
    FmpqHolder got;
    ASSERT_TRUE(dsn::arb_to_fmpq_simplest(got.q, x, 60));

    // |got - midpoint of x| < 10^-60, exactly in fmpq arithmetic.
    FmpqHolder actual_mid, diff, bound;
    arf_get_fmpq(actual_mid.q, arb_midref(x));
    fmpq_sub(diff.q, got.q, actual_mid.q);
    fmpq_abs(diff.q, diff.q);
    fmpq_one(bound.q);
    fmpq_div_fmpz(bound.q, bound.q, pow10);
    EXPECT_LT(fmpq_cmp(diff.q, bound.q), 0) << "got " << fmpq_str(got.q);
    // And far closer to the midpoint than the radius would allow: the
    // difference is below 10^-50.
    FmpqHolder bound50;
    {
        fmpz_t p50;
        fmpz_init(p50);
        fmpz_set_ui(p50, 10);
        fmpz_pow_ui(p50, p50, 50);
        fmpq_one(bound50.q);
        fmpq_div_fmpz(bound50.q, bound50.q, p50);
        fmpz_clear(p50);
    }
    EXPECT_LT(fmpq_cmp(diff.q, bound50.q), 0) << "got " << fmpq_str(got.q);

    arb_clear(x);
    fmpz_clear(pow10);
    fmpz_clear(numer);
}

TEST(RationalizeSimplest, ZeroBallIsZero) {
    arb_t x;
    arb_init(x);
    arb_zero(x);
    FmpqHolder got;
    fmpq_set_si(got.q, 7, 3);
    ASSERT_TRUE(dsn::arb_to_fmpq_simplest(got.q, x, 60));
    EXPECT_TRUE(fmpq_is_zero(got.q));

    // A tiny midpoint within 10^-digits of 0 also gives 0.
    ball_around(x, 1, 1, 0);
    arb_mul_2exp_si(x, x, -300);
    ASSERT_TRUE(dsn::arb_to_fmpq_simplest(got.q, x, 60));
    EXPECT_TRUE(fmpq_is_zero(got.q));

    // A midpoint of 2^-100000 at 60 digits is 0 as well, answered before
    // the midpoint is materialized as a fraction.
    arb_one(x);
    arb_mul_2exp_si(x, x, -100000);
    fmpq_set_si(got.q, 7, 3);
    ASSERT_TRUE(dsn::arb_to_fmpq_simplest(got.q, x, 60));
    EXPECT_TRUE(fmpq_is_zero(got.q));
    arb_clear(x);
}

TEST(RationalizeSimplest, NegativeValueIsExact) {
    arb_t x;
    arb_init(x);
    ball_around(x, -227, 113, 90);
    FmpqHolder got, want;
    ASSERT_TRUE(dsn::arb_to_fmpq_simplest(got.q, x, 60));
    fmpq_set_si(want.q, -227, 113);
    EXPECT_TRUE(fmpq_equal(got.q, want.q)) << "got " << fmpq_str(got.q);
    arb_clear(x);
}

TEST(RationalizeSimplest, NonFiniteIsRejected) {
    arb_t x;
    arb_init(x);
    arb_indeterminate(x);
    FmpqHolder got;
    fmpq_set_si(got.q, 7, 3);
    EXPECT_FALSE(dsn::arb_to_fmpq_simplest(got.q, x, 60));
    FmpqHolder untouched;
    fmpq_set_si(untouched.q, 7, 3);
    EXPECT_TRUE(fmpq_equal(got.q, untouched.q));

    arb_pos_inf(x);
    EXPECT_FALSE(dsn::arb_to_fmpq_simplest(got.q, x, 60));
    arb_clear(x);
}

TEST(RationalizeSimplest, AcbWrapperRejectsImaginaryPart) {
    acb_t z;
    acb_init(z);
    ball_around(acb_realref(z), 1, 113, 90);
    FmpqHolder got, want;
    fmpq_set_si(want.q, 1, 113);

    // Real input: accepted, exact.
    ASSERT_TRUE(dsn::acb_real_to_fmpq_simplest(got.q, z, 60));
    EXPECT_TRUE(fmpq_equal(got.q, want.q)) << "got " << fmpq_str(got.q);

    // Imaginary part 10^-20 at digits 60: rejected, output untouched.
    arb_set_si(acb_imagref(z), 10);
    arb_pow_ui(acb_imagref(z), acb_imagref(z), 20, 200);
    arb_inv(acb_imagref(z), acb_imagref(z), 200);
    fmpq_set_si(got.q, 7, 3);
    EXPECT_FALSE(dsn::acb_real_to_fmpq_simplest(got.q, z, 60));
    FmpqHolder untouched;
    fmpq_set_si(untouched.q, 7, 3);
    EXPECT_TRUE(fmpq_equal(got.q, untouched.q));

    // The same imaginary part at digits 10 is below threshold: accepted.
    EXPECT_TRUE(dsn::acb_real_to_fmpq_simplest(got.q, z, 10));
    EXPECT_TRUE(fmpq_equal(got.q, want.q)) << "got " << fmpq_str(got.q);

    // Negative real part through the wrapper.
    arb_zero(acb_imagref(z));
    ball_around(acb_realref(z), -9, 5, 0, 200);
    ASSERT_TRUE(dsn::acb_real_to_fmpq_simplest(got.q, z, 30));
    fmpq_set_si(want.q, -9, 5);
    EXPECT_TRUE(fmpq_equal(got.q, want.q)) << "got " << fmpq_str(got.q);
    acb_clear(z);
}
