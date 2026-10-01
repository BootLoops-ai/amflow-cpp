// SPDX-License-Identifier: MIT
// Tests for amflow::numeric FmpqPoly / RationalComplex / RationalFunction
// / RationalMatrix / get_poles.

#include <gtest/gtest.h>

#include <flint/fmpq.h>
#include <flint/fmpq_poly.h>

#include "amflow/numeric/options.hpp"
#include "amflow/numeric/rational.hpp"

#include <stdexcept>
#include <string>

namespace dsn = amflow::numeric;

// ---------------------------------------------------------------------------
//  FmpqPoly
// ---------------------------------------------------------------------------

TEST(FmpqPoly, DefaultIsZero) {
    dsn::FmpqPoly p;
    EXPECT_TRUE(p.is_zero());
    EXPECT_FALSE(p.is_one());
    EXPECT_EQ(p.length(), 0);
    EXPECT_EQ(p.degree(), -1);
    EXPECT_EQ(p.valuation(), -1);
}

TEST(FmpqPoly, SetSiAndOne) {
    dsn::FmpqPoly p;
    p.set_si(7);
    EXPECT_FALSE(p.is_zero());
    EXPECT_FALSE(p.is_one());
    EXPECT_EQ(p.degree(), 0);
    EXPECT_EQ(p.valuation(), 0);

    dsn::FmpqPoly q;
    q.set_one();
    EXPECT_TRUE(q.is_one());
    EXPECT_EQ(q.degree(), 0);
}

TEST(FmpqPoly, ValuationSkipsZeroCoeffs) {
    dsn::FmpqPoly p;
    p.set_coeff_si(0, 0);
    p.set_coeff_si(1, 0);
    p.set_coeff_si(3, 5);   // η^3 * 5
    EXPECT_EQ(p.degree(), 3);
    EXPECT_EQ(p.valuation(), 3);
}

TEST(FmpqPoly, MoveSemantics) {
    dsn::FmpqPoly a;
    a.set_coeff_si(2, 3);
    dsn::FmpqPoly b = std::move(a);
    EXPECT_EQ(b.degree(), 2);
}

TEST(FmpqPoly, CloneIsDeepCopy) {
    dsn::FmpqPoly a;
    a.set_si(11);
    dsn::FmpqPoly b = a.clone();
    a.set_zero();
    EXPECT_FALSE(b.is_zero());
    EXPECT_TRUE(a.is_zero());
}

TEST(FmpqPoly, ToStringUsesGivenVariable) {
    dsn::FmpqPoly p;
    p.set_coeff_si(0, 1);
    p.set_coeff_si(2, -3);
    auto s = p.to_string("x");
    EXPECT_NE(s.find("x"), std::string::npos);
}

// ---------------------------------------------------------------------------
//  RationalComplex
// ---------------------------------------------------------------------------

TEST(RationalComplex, DefaultIsZero) {
    dsn::RationalComplex z;
    EXPECT_TRUE(z.is_zero());
    EXPECT_TRUE(z.is_real());
}

TEST(RationalComplex, FromFractions) {
    auto z = dsn::RationalComplex::from_fractions(1, 2, 3, 4);   // 1/2 + 3/4 I
    EXPECT_FALSE(z.is_real());
    EXPECT_FALSE(z.is_zero());
}

TEST(RationalComplex, ZeroDenominatorThrows) {
    EXPECT_THROW(dsn::RationalComplex::from_fractions(1, 0, 0, 1),
                 std::invalid_argument);
}

TEST(RationalComplex, EqualityIsExact) {
    auto a = dsn::RationalComplex::from_fractions(2, 4, 1, 1);   // 1/2 + I
    auto b = dsn::RationalComplex::from_fractions(1, 2, 1, 1);   // 1/2 + I
    EXPECT_EQ(a, b);                                              // canonical-equal
}

TEST(RationalComplex, ToAcbRoundsRationalsExactly) {
    auto z = dsn::RationalComplex::from_fractions(3, 1, 0, 1);
    auto acb = z.to_acb();
    dsn::AcbValue three;
    three.set_si(3);
    EXPECT_NE(acb_eq(acb.raw(), three.raw()), 0);
}

// ---------------------------------------------------------------------------
//  RationalFunction
// ---------------------------------------------------------------------------

TEST(RationalFunction, DefaultIsZero) {
    dsn::RationalFunction r;
    EXPECT_TRUE(r.is_zero());
    EXPECT_FALSE(r.is_one());
    EXPECT_TRUE(r.is_polynomial());           // 0/1 has constant denominator
    EXPECT_TRUE(r.is_constant());
}

TEST(RationalFunction, ZeroDenominatorThrows) {
    dsn::FmpqPoly num, den;
    num.set_si(1);
    EXPECT_THROW(dsn::RationalFunction(std::move(num), std::move(den)),
                 std::invalid_argument);
}

TEST(RationalFunction, FromSiSi) {
    auto r = dsn::RationalFunction::from_si_si(2, 4);   // 2/4 = 1/2
    fmpq_t expected;
    fmpq_init(expected);
    fmpq_set_si(expected, 1, 2);
    auto eq_r = dsn::RationalFunction::from_fmpq(expected);
    EXPECT_EQ(r, eq_r);
    fmpq_clear(expected);
}

TEST(RationalFunction, MonomialPositiveIsPolynomial) {
    auto eta3 = dsn::RationalFunction::monomial(3);   // η^3
    EXPECT_TRUE(eta3.is_polynomial());
    EXPECT_EQ(eta3.degree_num(), 3);
    EXPECT_EQ(eta3.degree_den(), 0);
}

TEST(RationalFunction, MonomialNegativeIsRational) {
    auto inv_eta2 = dsn::RationalFunction::monomial(-2);   // 1/η^2
    EXPECT_FALSE(inv_eta2.is_polynomial());
    EXPECT_EQ(inv_eta2.degree_num(), 0);
    EXPECT_EQ(inv_eta2.degree_den(), 2);
    EXPECT_EQ(inv_eta2.net_valuation(), -2);
}

TEST(RationalFunction, AdditionReducesAutomatically) {
    auto half = dsn::RationalFunction::from_si_si(1, 2);
    auto third = dsn::RationalFunction::from_si_si(1, 3);
    auto sum = half + third;   // = 5/6
    auto five_six = dsn::RationalFunction::from_si_si(5, 6);
    EXPECT_EQ(sum, five_six);
}

TEST(RationalFunction, MultiplyDivideInverse) {
    auto two = dsn::RationalFunction::from_si(2);
    auto three = dsn::RationalFunction::from_si(3);
    auto product = two * three;
    auto six = dsn::RationalFunction::from_si(6);
    EXPECT_EQ(product, six);

    auto quotient = product / two;
    EXPECT_EQ(quotient, three);
}

TEST(RationalFunction, DivideByZeroThrows) {
    auto two = dsn::RationalFunction::from_si(2);
    dsn::RationalFunction zero;
    EXPECT_THROW(two / zero, std::domain_error);
}

TEST(RationalFunction, ReciprocalSwapsNumDen) {
    auto eta = dsn::RationalFunction::monomial(1);
    auto inv = eta.reciprocal();
    EXPECT_EQ(inv, dsn::RationalFunction::monomial(-1));
}

TEST(RationalFunction, ReciprocalOfZeroThrows) {
    dsn::RationalFunction zero;
    EXPECT_THROW(zero.reciprocal(), std::domain_error);
}

TEST(RationalFunction, DerivativeOfMonomial) {
    auto eta3 = dsn::RationalFunction::monomial(3);
    auto d = eta3.derivative();
    // d(η^3)/dη = 3η^2
    fmpq_t three;
    fmpq_init(three);
    fmpq_set_si(three, 3, 1);
    auto expected = dsn::RationalFunction::monomial(2, three);
    EXPECT_EQ(d, expected);
    fmpq_clear(three);
}

TEST(RationalFunction, SubstituteInverseEta) {
    // (η + 1) -> (1/η + 1) = (1 + η)/η, so substitute_inverse_eta should
    // yield (1 + η) / η  (after canonicalisation, num degree 1, den degree 1).
    dsn::FmpqPoly num, den;
    num.set_coeff_si(0, 1);
    num.set_coeff_si(1, 1);
    den.set_one();
    dsn::RationalFunction r(std::move(num), std::move(den));

    auto inv = r.substitute_inverse_eta();
    EXPECT_EQ(inv.degree_num(), 1);
    EXPECT_EQ(inv.degree_den(), 1);
    // Verify numerically: at η=2, original is 3, inverse-substituted at η=2 is
    // (1+2)/2 = 3/2.
    dsn::AcbValue x;
    x.set_si(2);
    auto val = inv.evaluate(x.raw());
    fmpq_t three_halves;
    fmpq_init(three_halves);
    fmpq_set_si(three_halves, 3, 2);
    dsn::AcbValue expected;
    expected.set_fmpq(three_halves);
    EXPECT_NE(acb_eq(val.raw(), expected.raw()), 0);
    fmpq_clear(three_halves);
}

TEST(RationalFunction, EvaluateAtRegularPoint) {
    // r = (η + 2) / (η + 3); at η=0, r = 2/3.
    dsn::FmpqPoly num, den;
    num.set_coeff_si(0, 2);
    num.set_coeff_si(1, 1);
    den.set_coeff_si(0, 3);
    den.set_coeff_si(1, 1);
    dsn::RationalFunction r(std::move(num), std::move(den));

    dsn::AcbValue x;
    x.set_zero();
    auto val = r.evaluate(x.raw());
    fmpq_t two_thirds;
    fmpq_init(two_thirds);
    fmpq_set_si(two_thirds, 2, 3);
    dsn::AcbValue expected;
    expected.set_fmpq(two_thirds);
    EXPECT_NE(acb_overlaps(val.raw(), expected.raw()), 0);
    fmpq_clear(two_thirds);
}

TEST(RationalFunction, EvaluateAtPoleThrows) {
    // r = 1 / η; at η=0, denominator vanishes.
    auto inv = dsn::RationalFunction::monomial(-1);
    dsn::AcbValue x;
    x.set_zero();
    EXPECT_THROW(inv.evaluate(x.raw()), std::domain_error);
}

// ---------------------------------------------------------------------------
//  RationalMatrix
// ---------------------------------------------------------------------------

TEST(RationalMatrix, IdentityHasOnesOnDiagonal) {
    auto m = dsn::RationalMatrix::identity(3);
    EXPECT_EQ(m.rows(), 3u);
    EXPECT_EQ(m.cols(), 3u);
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            if (i == j) EXPECT_TRUE(m(i, j).is_one());
            else        EXPECT_TRUE(m(i, j).is_zero());
        }
    }
}

TEST(RationalMatrix, DensityCounts) {
    dsn::RationalMatrix m(2, 2);
    EXPECT_DOUBLE_EQ(m.density(), 0.0);
    m(0, 0) = dsn::RationalFunction::from_si(1);
    EXPECT_DOUBLE_EQ(m.density(), 0.25);
    m(1, 1) = dsn::RationalFunction::from_si(2);
    EXPECT_DOUBLE_EQ(m.density(), 0.5);
}

TEST(RationalMatrix, AtBoundsCheck) {
    dsn::RationalMatrix m(2, 2);
    EXPECT_THROW(m.at(2, 0), std::out_of_range);
    EXPECT_THROW(m.at(0, 5), std::out_of_range);
}

TEST(RationalMatrix, DenominatorLcmIsMonic) {
    dsn::RationalMatrix m(1, 2);
    m(0, 0) = dsn::RationalFunction::monomial(-1);   // 1 / η
    m(0, 1) = dsn::RationalFunction::monomial(-2);   // 1 / η^2
    auto lcm = m.denominator_lcm();
    // LCM(η, η^2) = η^2.
    EXPECT_EQ(lcm.degree(), 2);
    fmpq_t lead;
    fmpq_init(lead);
    fmpq_poly_get_coeff_fmpq(lead, lcm.raw(), 2);
    EXPECT_NE(fmpq_is_one(lead), 0);
    fmpq_clear(lead);
}

// ---------------------------------------------------------------------------
//  get_poles
// ---------------------------------------------------------------------------

TEST(GetPoles, EmptyMatrixReturnsEmpty) {
    dsn::RationalMatrix m;
    auto poles = dsn::get_poles(m);
    EXPECT_TRUE(poles.empty());
}

TEST(GetPoles, ConstantMatrixHasNoPoles) {
    auto m = dsn::RationalMatrix::identity(2);
    auto poles = dsn::get_poles(m);
    EXPECT_TRUE(poles.empty());
}

TEST(GetPoles, SimplePoleAtZero) {
    dsn::RationalMatrix m(1, 1);
    m(0, 0) = dsn::RationalFunction::monomial(-1);   // 1/η
    auto poles = dsn::get_poles(m);
    ASSERT_EQ(poles.size(), 1u);
    EXPECT_TRUE(poles[0].is_zero());
}

TEST(GetPoles, IntegerPolesAreRationalised) {
    // r = 1 / ((η - 2)(η + 3)); poles at 2 and -3.
    dsn::FmpqPoly num, den;
    num.set_one();
    // (η - 2)(η + 3) = η^2 + η - 6
    den.set_coeff_si(0, -6);
    den.set_coeff_si(1, 1);
    den.set_coeff_si(2, 1);

    dsn::RationalMatrix m(1, 1);
    m(0, 0) = dsn::RationalFunction(std::move(num), std::move(den));

    auto poles = dsn::get_poles(m);
    EXPECT_EQ(poles.size(), 2u);

    bool found_two = false, found_neg_three = false;
    auto two = dsn::RationalComplex::from_fractions(2, 1, 0, 1);
    auto neg_three = dsn::RationalComplex::from_fractions(-3, 1, 0, 1);
    for (const auto& p : poles) {
        if (p == two)        found_two = true;
        if (p == neg_three)  found_neg_three = true;
    }
    EXPECT_TRUE(found_two);
    EXPECT_TRUE(found_neg_three);
}
