// SPDX-License-Identifier: MIT
// Tests for amflow::numeric::AcbValue / AcbVector / AcbMatrix / AcbPoly.
//
// Domain contract: docs/REFACTOR_DESIGN.md §3.1, §5.1.

#include <gtest/gtest.h>

#include <flint/acb.h>
#include <flint/arb.h>
#include <flint/fmpq.h>

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
