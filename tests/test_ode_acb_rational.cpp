// SPDX-License-Identifier: MIT
// Tests for amflow::ode::acb_rational (Layer 7 helper).

#include <gtest/gtest.h>

#include "amflow/ode/acb_rational.hpp"
#include "amflow/numeric/acb_wrap.hpp"
#include "amflow/numeric/options.hpp"
#include "amflow/numeric/rational.hpp"

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/acb_poly.h>
#include <flint/fmpq.h>

#include <cmath>
#include <stdexcept>

namespace nm  = amflow::numeric;
namespace ode = amflow::ode;

namespace {

class AcbRationalTest : public ::testing::Test {
protected:
    void SetUp() override    { nm::set_default_options(); }
    void TearDown() override { nm::set_default_options(); }
};

bool acb_close_to_double(const nm::AcbValue& v, double r, double i,
                         double tol = 1e-15) {
    arb_t target_re, target_im;
    arb_init(target_re); arb_init(target_im);
    arb_set_d(target_re, r);
    arb_set_d(target_im, i);
    arb_t diff_re, diff_im;
    arb_init(diff_re); arb_init(diff_im);
    arb_sub(diff_re, acb_realref(v.raw()), target_re, 200);
    arb_sub(diff_im, acb_imagref(v.raw()), target_im, 200);
    arf_t mag_re, mag_im;
    arf_init(mag_re); arf_init(mag_im);
    arf_abs(mag_re, arb_midref(diff_re));
    arf_abs(mag_im, arb_midref(diff_im));
    bool ok = (arf_cmp_d(mag_re, tol) < 0) && (arf_cmp_d(mag_im, tol) < 0);
    arf_clear(mag_re); arf_clear(mag_im);
    arb_clear(diff_re); arb_clear(diff_im);
    arb_clear(target_re); arb_clear(target_im);
    return ok;
}

nm::FmpqPoly poly_si(std::initializer_list<long> coeffs) {
    nm::FmpqPoly p;
    long n = 0;
    for (long c : coeffs) { p.set_coeff_si(n, c); ++n; }
    return p;
}

}  // namespace

// ============================================================================
//  Construction & predicates
// ============================================================================

TEST_F(AcbRationalTest, DefaultIsZero) {
    ode::AcbRationalFunction r;
    EXPECT_TRUE(r.is_zero());
    EXPECT_FALSE(r.is_one());
}

TEST_F(AcbRationalTest, FromSi_OneIsOne) {
    auto r = ode::AcbRationalFunction::from_si(1);
    EXPECT_FALSE(r.is_zero());
    EXPECT_TRUE(r.is_one());
}

TEST_F(AcbRationalTest, FromRationalLifts) {
    // 5 / (1 + eta) -> AcbRational with same num/den (coeffs as acb).
    auto r = nm::RationalFunction(poly_si({5}), poly_si({1, 1}));
    auto a = ode::AcbRationalFunction::from_rational(r);
    EXPECT_FALSE(a.is_zero());
    auto v = a.value_at_zero();
    EXPECT_TRUE(acb_close_to_double(v, 5.0, 0.0));
}

// ============================================================================
//  Monomial construction (regression: monomial(k<0) used to leak den[0] = 1)
// ============================================================================

TEST_F(AcbRationalTest, MonomialNegativeK_DenominatorIsPureEtaPower) {
    // monomial(-2, c) should produce c / eta^2.  At eta -> 0 the entry blows
    // up; check that strip_eta_prefix sees both num and den having 0 coeffs
    // up to k=2 -> after multiplying by eta^2 it should equal c.
    nm::AcbValue c; c.set_si(7);
    auto r = ode::AcbRationalFunction::monomial(-2, c.raw());
    r.multiply_by_eta_power(2);   // r * eta^2 -> 7
    auto v = r.value_at_zero();
    EXPECT_TRUE(acb_close_to_double(v, 7.0, 0.0));
}

TEST_F(AcbRationalTest, MonomialPositiveK) {
    nm::AcbValue c; c.set_si(3);
    auto r = ode::AcbRationalFunction::monomial(2, c.raw());
    // r = 3*eta^2;  multiply by eta^-2 -> 3.
    r.multiply_by_eta_power(-2);
    auto v = r.value_at_zero();
    EXPECT_TRUE(acb_close_to_double(v, 3.0, 0.0));
}

// ============================================================================
//  Arithmetic (no auto-reduction; check results via value_at_zero only)
// ============================================================================

TEST_F(AcbRationalTest, AdditionEvaluates) {
    auto a = ode::AcbRationalFunction::from_si(2);
    auto b = ode::AcbRationalFunction::from_si(3);
    a += b;
    EXPECT_TRUE(acb_close_to_double(a.value_at_zero(), 5.0, 0.0));
}

TEST_F(AcbRationalTest, MultiplicationEvaluates) {
    // (1 + eta) * (1 - eta) = 1 - eta^2  -> at eta=0 = 1.
    auto a = ode::AcbRationalFunction::from_rational(
        nm::RationalFunction::from_polynomial(poly_si({1, 1})));
    auto b = ode::AcbRationalFunction::from_rational(
        nm::RationalFunction::from_polynomial(poly_si({1, -1})));
    a *= b;
    EXPECT_TRUE(acb_close_to_double(a.value_at_zero(), 1.0, 0.0));
}

TEST_F(AcbRationalTest, DivisionByZeroThrows) {
    auto a = ode::AcbRationalFunction::from_si(1);
    ode::AcbRationalFunction zero;
    EXPECT_THROW(a /= zero, std::domain_error);
}

TEST_F(AcbRationalTest, NegationFlipsSign) {
    auto a = ode::AcbRationalFunction::from_si(5);
    auto b = -a;
    EXPECT_TRUE(acb_close_to_double(b.value_at_zero(), -5.0, 0.0));
}

// ============================================================================
//  Derivative
// ============================================================================

TEST_F(AcbRationalTest, DerivativeOfPolynomial) {
    // d/deta (eta^2 + 3*eta + 5) = 2*eta + 3 -> at eta=0 = 3.
    auto a = ode::AcbRationalFunction::from_rational(
        nm::RationalFunction::from_polynomial(poly_si({5, 3, 1})));
    auto d = a.derivative();
    EXPECT_TRUE(acb_close_to_double(d.value_at_zero(), 3.0, 0.0));
}

// ============================================================================
//  strip_eta_prefix / value_at_zero
// ============================================================================

TEST_F(AcbRationalTest, ValueAtZero_WithCancelablePrefix) {
    // (eta * 7) / (eta * 1) -> after strip = 7 / 1 = 7.
    auto a = ode::AcbRationalFunction::from_acb(([](){
        nm::AcbValue v; v.set_si(7); return v;
    })().raw());
    a.multiply_by_eta_power(1);
    a.multiply_by_eta_power(-1);   // num *= eta, den *= eta
    auto v = a.value_at_zero();
    EXPECT_TRUE(acb_close_to_double(v, 7.0, 0.0));
}

TEST_F(AcbRationalTest, ValueAtZero_RemainingPoleThrows) {
    // 1 / eta -> after strip still has den vanishing.
    auto a = ode::AcbRationalFunction::from_si(1);
    a.multiply_by_eta_power(-1);
    EXPECT_THROW(a.value_at_zero(), std::domain_error);
}

// ============================================================================
//  Matrix ops
// ============================================================================

TEST_F(AcbRationalTest, MatrixIdentityIsIdentity) {
    auto I = ode::AcbRationalMatrix::identity(2);
    ASSERT_EQ(I.rows(), 2u);
    ASSERT_EQ(I.cols(), 2u);
    EXPECT_TRUE(I(0, 0).is_one());
    EXPECT_TRUE(I(0, 1).is_zero());
    EXPECT_TRUE(I(1, 0).is_zero());
    EXPECT_TRUE(I(1, 1).is_one());
}

TEST_F(AcbRationalTest, MatrixMatmulMatchesIdentityProduct) {
    auto A = ode::AcbRationalMatrix(2, 2);
    A(0, 0) = ode::AcbRationalFunction::from_si(1);
    A(0, 1) = ode::AcbRationalFunction::from_si(2);
    A(1, 0) = ode::AcbRationalFunction::from_si(3);
    A(1, 1) = ode::AcbRationalFunction::from_si(4);
    auto I = ode::AcbRationalMatrix::identity(2);
    auto P = A.matmul(I);
    EXPECT_TRUE(acb_close_to_double(P(0, 0).value_at_zero(), 1.0, 0.0));
    EXPECT_TRUE(acb_close_to_double(P(0, 1).value_at_zero(), 2.0, 0.0));
    EXPECT_TRUE(acb_close_to_double(P(1, 0).value_at_zero(), 3.0, 0.0));
    EXPECT_TRUE(acb_close_to_double(P(1, 1).value_at_zero(), 4.0, 0.0));
}

TEST_F(AcbRationalTest, MatrixSubAndAdd) {
    auto A = ode::AcbRationalMatrix(1, 1);
    A(0, 0) = ode::AcbRationalFunction::from_si(7);
    auto B = ode::AcbRationalMatrix(1, 1);
    B(0, 0) = ode::AcbRationalFunction::from_si(3);
    auto S = A + B;
    auto D = A - B;
    EXPECT_TRUE(acb_close_to_double(S(0, 0).value_at_zero(), 10.0, 0.0));
    EXPECT_TRUE(acb_close_to_double(D(0, 0).value_at_zero(),  4.0, 0.0));
}

TEST_F(AcbRationalTest, ValueAtZero_Matrix) {
    nm::RationalMatrix r(1, 1);
    r(0, 0) = nm::RationalFunction::from_si(7);
    auto a = ode::AcbRationalMatrix::from_rational(r);
    acb_mat_t out;
    acb_mat_init(out, 1, 1);
    ode::acb_rational_matrix_value_at_zero(out, a);
    nm::AcbValue v;
    acb_set(v.raw(), acb_mat_entry(out, 0, 0));
    EXPECT_TRUE(acb_close_to_double(v, 7.0, 0.0));
    acb_mat_clear(out);
}

TEST_F(AcbRationalTest, EtaAtZero_MultipliesByEtaThenEvaluates) {
    // m(0,0) = 5/eta -> eta * m(0,0) = 5 -> at eta=0 = 5.
    nm::RationalMatrix r(1, 1);
    r(0, 0) = nm::RationalFunction(poly_si({5}), poly_si({0, 1}));
    auto a = ode::AcbRationalMatrix::from_rational(r);
    acb_mat_t out;
    acb_mat_init(out, 1, 1);
    ode::acb_rational_matrix_eta_at_zero(out, a);
    nm::AcbValue v;
    acb_set(v.raw(), acb_mat_entry(out, 0, 0));
    EXPECT_TRUE(acb_close_to_double(v, 5.0, 0.0));
    acb_mat_clear(out);
}
