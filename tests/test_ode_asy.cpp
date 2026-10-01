// SPDX-License-Identifier: MIT
// Tests for amflow::ode::asy (Layer 4).

#include <gtest/gtest.h>

#include "amflow/ode/asy.hpp"
#include "amflow/numeric/options.hpp"
#include "amflow/numeric/rational.hpp"

#include <flint/acb.h>
#include <flint/arb.h>
#include <flint/fmpq.h>

#include <cmath>
#include <vector>

namespace nm = amflow::numeric;
namespace ode = amflow::ode;

namespace {

class AsyTest : public ::testing::Test {
protected:
    void SetUp() override    { nm::set_default_options(); }
    void TearDown() override { nm::set_default_options(); }
};

nm::FmpqPoly make_poly(std::initializer_list<long> coeffs) {
    nm::FmpqPoly p;
    long n = 0;
    for (long c : coeffs) { p.set_coeff_si(n, c); ++n; }
    return p;
}

nm::AcbValue acb_si(long n) { nm::AcbValue v; v.set_si(n); return v; }

bool acb_close_to_double(const nm::AcbValue& v,
                         double r, double i, double tol = 1e-15) {
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

std::vector<std::vector<nm::AcbValue>> exp_from_doubles(
    std::initializer_list<std::initializer_list<double>> rows) {
    std::vector<std::vector<nm::AcbValue>> out;
    out.reserve(rows.size());
    for (const auto& r : rows) {
        std::vector<nm::AcbValue> row;
        row.reserve(r.size());
        for (double d : r) {
            nm::AcbValue v;
            v.set_d_d(d, 0.0);
            row.push_back(std::move(v));
        }
        out.push_back(std::move(row));
    }
    return out;
}

}  // namespace

// =========================================================================
//  to_power_series
// =========================================================================

TEST_F(AsyTest, ToPowerSeries_Polynomial) {
    auto r = nm::RationalFunction::from_polynomial(make_poly({1, 2, 3}));
    auto ps = ode::to_power_series(r);
    EXPECT_EQ(ps.leading_offset, 0L);
    ASSERT_EQ(ps.coeffs.size(), 3u);
    EXPECT_TRUE(acb_close_to_double(ps.coeffs[0], 1.0, 0.0));
    EXPECT_TRUE(acb_close_to_double(ps.coeffs[1], 2.0, 0.0));
    EXPECT_TRUE(acb_close_to_double(ps.coeffs[2], 3.0, 0.0));
}

TEST_F(AsyTest, ToPowerSeries_LeadingNegativeOffset) {
    // r(eta) = (1 + eta) / eta^2 -> offset = -2, coeffs = {1, 1}
    nm::FmpqPoly num = make_poly({1, 1});
    nm::FmpqPoly den = make_poly({0, 0, 1});
    auto r = nm::RationalFunction(std::move(num), std::move(den));
    auto ps = ode::to_power_series(r);
    EXPECT_EQ(ps.leading_offset, -2L);
    ASSERT_EQ(ps.coeffs.size(), 2u);
    EXPECT_TRUE(acb_close_to_double(ps.coeffs[0], 1.0, 0.0));
    EXPECT_TRUE(acb_close_to_double(ps.coeffs[1], 1.0, 0.0));
}

TEST_F(AsyTest, ToPowerSeries_LeadingPositiveOffset) {
    // r(eta) = eta^3 * (2 + 5 eta) -> offset = 3, coeffs = {2, 5}
    nm::FmpqPoly num = make_poly({0, 0, 0, 2, 5});
    auto r = nm::RationalFunction::from_polynomial(std::move(num));
    auto ps = ode::to_power_series(r);
    EXPECT_EQ(ps.leading_offset, 3L);
    ASSERT_EQ(ps.coeffs.size(), 2u);
    EXPECT_TRUE(acb_close_to_double(ps.coeffs[0], 2.0, 0.0));
    EXPECT_TRUE(acb_close_to_double(ps.coeffs[1], 5.0, 0.0));
}

TEST_F(AsyTest, ToPowerSeries_ZeroIsEmpty) {
    auto ps = ode::to_power_series(nm::RationalFunction());
    EXPECT_TRUE(ps.is_zero());
    EXPECT_TRUE(ps.coeffs.empty());
}

TEST_F(AsyTest, ToPowerSeries_RejectsNonZeroPole) {
    // r = 1/(eta - 1):  pole at eta = 1, not at eta = 0.
    nm::FmpqPoly num = make_poly({1});
    nm::FmpqPoly den = make_poly({-1, 1});
    auto r = nm::RationalFunction(std::move(num), std::move(den));
    EXPECT_THROW(ode::to_power_series(r), std::invalid_argument);
}

TEST_F(AsyTest, ToPowerSeriesMatrix_AppliesElementwise) {
    nm::RationalMatrix m(1, 2);
    m(0, 0) = nm::RationalFunction::from_polynomial(make_poly({1, 2}));
    m(0, 1) = nm::RationalFunction::from_polynomial(make_poly({0, 0, 4}));
    auto out = ode::to_power_series_matrix(m);
    ASSERT_EQ(out.size(), 1u);
    ASSERT_EQ(out[0].size(), 2u);
    EXPECT_EQ(out[0][0].leading_offset, 0L);
    ASSERT_EQ(out[0][0].coeffs.size(), 2u);
    EXPECT_TRUE(acb_close_to_double(out[0][0].coeffs[0], 1.0, 0.0));
    EXPECT_TRUE(acb_close_to_double(out[0][0].coeffs[1], 2.0, 0.0));
    EXPECT_EQ(out[0][1].leading_offset, 2L);
    ASSERT_EQ(out[0][1].coeffs.size(), 1u);
    EXPECT_TRUE(acb_close_to_double(out[0][1].coeffs[0], 4.0, 0.0));
}

// =========================================================================
//  times_poly / inverse_poly / rational_expansion
// =========================================================================

TEST_F(AsyTest, TimesPoly_BasicConvolution) {
    // f = (1, 2, 3, 4),  poly = (5, 6) -> length-4 truncation of f * poly.
    std::vector<nm::AcbValue> f, poly;
    for (long c : {1, 2, 3, 4}) f.push_back(acb_si(c));
    for (long c : {5, 6})       poly.push_back(acb_si(c));
    auto r = ode::times_poly(f, poly);
    ASSERT_EQ(r.size(), 4u);
    EXPECT_TRUE(acb_close_to_double(r[0],  5.0, 0.0));
    EXPECT_TRUE(acb_close_to_double(r[1], 16.0, 0.0));
    EXPECT_TRUE(acb_close_to_double(r[2], 27.0, 0.0));
    EXPECT_TRUE(acb_close_to_double(r[3], 38.0, 0.0));
}

TEST_F(AsyTest, InversePoly_BasicDivision) {
    // 1 / (1 - eta) = 1 + eta + eta^2 + eta^3
    std::vector<nm::AcbValue> f(4);
    f[0].set_si(1);
    std::vector<nm::AcbValue> poly;
    poly.push_back(acb_si(1));
    poly.push_back(acb_si(-1));
    auto r = ode::inverse_poly(f, poly);
    for (std::size_t k = 0; k < 4; ++k) {
        EXPECT_TRUE(acb_close_to_double(r[k], 1.0, 0.0)) << "k = " << k;
    }
}

TEST_F(AsyTest, InversePoly_ThrowsWhenLeadingZero) {
    std::vector<nm::AcbValue> f(2);
    f[0].set_si(1);
    std::vector<nm::AcbValue> poly;
    poly.push_back(acb_si(0));
    poly.push_back(acb_si(1));
    EXPECT_THROW(ode::inverse_poly(f, poly), std::domain_error);
}

TEST_F(AsyTest, RationalExpansion_DecayingGeometricSeries) {
    // (num/den) * exp where num=1, den=1-eta, exp = (1, 0, 0, 0)  -> 1/(1-eta)
    std::vector<nm::AcbValue> num; num.push_back(acb_si(1));
    std::vector<nm::AcbValue> den; den.push_back(acb_si(1)); den.push_back(acb_si(-1));
    std::vector<nm::AcbValue> exp(4);
    exp[0].set_si(1);
    auto r = ode::rational_expansion(num, den, exp);
    for (std::size_t k = 0; k < 4; ++k) {
        EXPECT_TRUE(acb_close_to_double(r[k], 1.0, 0.0));
    }
}

TEST_F(AsyTest, MapRationalExpansion_SumsPerRow) {
    // Single row, two columns:
    //   nums[0][0]=1, dens[0][0]=1-eta;  expansion[0]=(1,0,0,0)  -> 1/(1-eta)
    //   nums[0][1]=1, dens[0][1]=1-eta;  expansion[1]=(0,1,0,0)  -> eta/(1-eta)
    // Sum = (1 + eta) / (1 - eta) = 1 + 2 eta + 2 eta^2 + 2 eta^3
    std::vector<std::vector<std::vector<nm::AcbValue>>> nums(1);
    std::vector<std::vector<std::vector<nm::AcbValue>>> dens(1);
    nums[0].resize(2); dens[0].resize(2);
    nums[0][0].push_back(acb_si(1));
    dens[0][0].push_back(acb_si(1)); dens[0][0].push_back(acb_si(-1));
    nums[0][1].push_back(acb_si(1));
    dens[0][1].push_back(acb_si(1)); dens[0][1].push_back(acb_si(-1));

    std::vector<std::vector<nm::AcbValue>> expansions(2);
    expansions[0].resize(4); expansions[0][0].set_si(1);
    expansions[1].resize(4); expansions[1][1].set_si(1);

    auto out = ode::map_rational_expansion(nums, dens, expansions);
    ASSERT_EQ(out.size(), 1u);
    ASSERT_EQ(out[0].size(), 4u);
    EXPECT_TRUE(acb_close_to_double(out[0][0], 1.0, 0.0));
    EXPECT_TRUE(acb_close_to_double(out[0][1], 2.0, 0.0));
    EXPECT_TRUE(acb_close_to_double(out[0][2], 2.0, 0.0));
    EXPECT_TRUE(acb_close_to_double(out[0][3], 2.0, 0.0));
}

// =========================================================================
//  extend / plus / rescale
// =========================================================================

TEST_F(AsyTest, ExtendExpansion_PadsRows) {
    nm::GlobalScope s; s.expansion.x_order = 4; s.commit();    // row length = 5

    auto exp = exp_from_doubles({{1, 2, 0, 0, 0}});
    ode::extend_expansion(exp, 3);
    ASSERT_EQ(exp.size(), 3u);
    EXPECT_EQ(exp[1].size(), 5u);
    EXPECT_TRUE(exp[1][0].is_zero());
    EXPECT_TRUE(exp[2][0].is_zero());
}

TEST_F(AsyTest, ExtendExpansion_NoOpWhenAlreadyLargeEnough) {
    nm::GlobalScope s; s.expansion.x_order = 3; s.commit();
    auto exp = exp_from_doubles({{1, 0, 0, 0}, {2, 0, 0, 0}});
    ode::extend_expansion(exp, 1);
    EXPECT_EQ(exp.size(), 2u);
}

TEST_F(AsyTest, PlusExpansion_SumsAndPads) {
    nm::GlobalScope s; s.expansion.x_order = 3; s.commit();   // row length = 4

    auto a = exp_from_doubles({{1, 2, 0, 0}});
    auto b = exp_from_doubles({{0, 1, 0, 0}, {3, 0, 0, 0}});
    std::vector<std::vector<std::vector<nm::AcbValue>>> exps;
    exps.push_back(std::move(a));
    exps.push_back(std::move(b));

    auto sum = ode::plus_expansion(exps);
    ASSERT_EQ(sum.size(), 2u);
    EXPECT_TRUE(acb_close_to_double(sum[0][0], 1.0, 0.0));
    EXPECT_TRUE(acb_close_to_double(sum[0][1], 3.0, 0.0));
    EXPECT_TRUE(acb_close_to_double(sum[1][0], 3.0, 0.0));
}

TEST_F(AsyTest, RescaleExpansion_ShiftsToTheRight) {
    nm::GlobalScope s; s.expansion.x_order = 4; s.commit();    // row length 5

    auto exp = exp_from_doubles({{1, 2, 3, 0, 0}});
    auto scaled = ode::rescale_expansion(exp, 2);
    ASSERT_EQ(scaled.size(), 1u);
    ASSERT_EQ(scaled[0].size(), 5u);
    EXPECT_TRUE(scaled[0][0].is_zero());
    EXPECT_TRUE(scaled[0][1].is_zero());
    EXPECT_TRUE(acb_close_to_double(scaled[0][2], 1.0, 0.0));
    EXPECT_TRUE(acb_close_to_double(scaled[0][3], 2.0, 0.0));
    EXPECT_TRUE(acb_close_to_double(scaled[0][4], 3.0, 0.0));
}

// =========================================================================
//  evaluate
// =========================================================================

TEST_F(AsyTest, EvaluateExpansion_PolynomialNoLog) {
    auto exp = exp_from_doubles({{2, 3, 5}});
    nm::AcbValue x; x.set_si(4);
    auto v = ode::evaluate_expansion(exp, x.raw());
    EXPECT_TRUE(acb_close_to_double(v, 94.0, 0.0));
}

TEST_F(AsyTest, EvaluateExpansion_WithLogPower) {
    // exp = {{0, 1}, {2, 0}}.  At eta = e:
    //    row 0 contributes e
    //    row 1 contributes 2 * log(e) = 2
    //    total = e + 2
    nm::AcbValue x;
    arb_const_e(acb_realref(x.raw()), nm::working_prec_bits());
    arb_zero(acb_imagref(x.raw()));

    auto exp = exp_from_doubles({{0, 1}, {2, 0}});
    auto v = ode::evaluate_expansion(exp, x.raw());
    EXPECT_TRUE(acb_close_to_double(v, std::exp(1.0) + 2.0, 0.0));
}

TEST_F(AsyTest, EvaluateAsyTerm_ScalesByX0Mu) {
    // mu = 1/2, exp = {{1}}, x0 = 9 -> sqrt(9) = 3
    nm::AcbValue mu;
    fmpq_t half;
    fmpq_init(half);
    fmpq_set_si(half, 1, 2);
    arb_set_fmpq(acb_realref(mu.raw()), half, nm::working_prec_bits());
    fmpq_clear(half);

    nm::AcbValue x; x.set_si(9);
    auto exp = exp_from_doubles({{1}});
    auto v = ode::evaluate_asy_term(mu.raw(), exp, x.raw());
    EXPECT_TRUE(acb_close_to_double(v, 3.0, 0.0));
}

TEST_F(AsyTest, EvaluateAsyExpansion_SumsAllTerms) {
    nm::GlobalScope s; s.expansion.x_order = 3; s.commit();
    // term 0: mu = 0,  exp = {{2}}      -> 2
    // term 1: mu = 1,  exp = {{3}}      -> 3 * x0
    // x0 = 4 -> total = 2 + 3*4 = 14
    ode::AsyExpansion asy;
    nm::AcbValue mu0, mu1;
    mu0.set_si(0); mu1.set_si(1);
    asy.emplace_back(std::move(mu0), exp_from_doubles({{2}}));
    asy.emplace_back(std::move(mu1), exp_from_doubles({{3}}));
    nm::AcbValue x; x.set_si(4);
    auto v = ode::evaluate_asy_expansion(asy, x.raw());
    EXPECT_TRUE(acb_close_to_double(v, 14.0, 0.0));
}

// =========================================================================
//  mu_differ_by_integer
// =========================================================================

TEST_F(AsyTest, MuDifferByInteger_BasicRealCases) {
    nm::AcbValue m1, m2, m3;
    m1.set_d_d(0.5, 0.0);
    m2.set_d_d(2.5, 0.0);   // diff = 2 -> integer
    m3.set_d_d(0.7, 0.0);   // diff = 0.2 -> not integer
    EXPECT_TRUE(ode::mu_differ_by_integer(m1.raw(), m2.raw()));
    EXPECT_FALSE(ode::mu_differ_by_integer(m1.raw(), m3.raw()));
}

TEST_F(AsyTest, MuDifferByInteger_ImaginaryShiftIsRejected) {
    nm::AcbValue m1, m2;
    m1.set_d_d(0.5, 1.0);
    m2.set_d_d(0.5, 2.0);
    EXPECT_FALSE(ode::mu_differ_by_integer(m1.raw(), m2.raw()));
}

TEST_F(AsyTest, MuDifferByInteger_EqualMusAreInteger) {
    nm::AcbValue m1, m2;
    m1.set_d_d(0.7, 0.4);
    m2.set_d_d(0.7, 0.4);
    EXPECT_TRUE(ode::mu_differ_by_integer(m1.raw(), m2.raw()));
}

// =========================================================================
//  plus_rule_set / union_rule_set
// =========================================================================

TEST_F(AsyTest, PlusRuleSet_CombinesIntegerOffsets) {
    nm::GlobalScope s; s.expansion.x_order = 3; s.commit();   // row length = 4

    // Region: mu=0 -> exp = {{1, 0, 0, 0}};  mu=1 -> exp = {{2, 0, 0, 0}}
    // Combined: min_mu = 0, exp = {{1, 2, 0, 0}}
    nm::AcbValue m0, m1;
    m0.set_si(0); m1.set_si(1);
    ode::AsyExpansion region;
    region.emplace_back(std::move(m0), exp_from_doubles({{1, 0, 0, 0}}));
    region.emplace_back(std::move(m1), exp_from_doubles({{2, 0, 0, 0}}));

    auto term = ode::plus_rule_set(region);
    EXPECT_TRUE(acb_close_to_double(term.mu, 0.0, 0.0));
    ASSERT_EQ(term.exp.size(), 1u);
    ASSERT_EQ(term.exp[0].size(), 4u);
    EXPECT_TRUE(acb_close_to_double(term.exp[0][0], 1.0, 0.0));
    EXPECT_TRUE(acb_close_to_double(term.exp[0][1], 2.0, 0.0));
}

TEST_F(AsyTest, PlusRuleSet_NegativeShiftFromMinFirstKey) {
    nm::GlobalScope s; s.expansion.x_order = 3; s.commit();

    // mu=2 -> exp[0][0]=1; mu=-1 -> exp[0][0]=4.  min_mu = -1; the second
    // term enters at offset 0; first term shifts by 3.
    //  combined exp = {{4, 0, 0, 1}}.
    nm::AcbValue m2, mneg1;
    m2.set_si(2); mneg1.set_si(-1);
    ode::AsyExpansion region;
    region.emplace_back(std::move(m2),    exp_from_doubles({{1, 0, 0, 0}}));
    region.emplace_back(std::move(mneg1), exp_from_doubles({{4, 0, 0, 0}}));

    auto term = ode::plus_rule_set(region);
    EXPECT_TRUE(acb_close_to_double(term.mu, -1.0, 0.0));
    ASSERT_EQ(term.exp[0].size(), 4u);
    EXPECT_TRUE(acb_close_to_double(term.exp[0][0], 4.0, 0.0));
    EXPECT_TRUE(acb_close_to_double(term.exp[0][1], 0.0, 0.0));
    EXPECT_TRUE(acb_close_to_double(term.exp[0][2], 0.0, 0.0));
    EXPECT_TRUE(acb_close_to_double(term.exp[0][3], 1.0, 0.0));
}

TEST_F(AsyTest, UnionRuleSet_GroupsByIntegerEquivalence) {
    nm::GlobalScope s; s.expansion.x_order = 3; s.commit();

    ode::AsyExpansion asy;
    nm::AcbValue mu0, mu05, mum1;
    mu0.set_si(0);
    fmpq_t half;
    fmpq_init(half); fmpq_set_si(half, 1, 2);
    arb_set_fmpq(acb_realref(mu05.raw()), half, nm::working_prec_bits());
    fmpq_clear(half);
    mum1.set_si(-1);

    asy.emplace_back(std::move(mu0),  exp_from_doubles({{1, 0, 0, 0}}));
    asy.emplace_back(std::move(mu05), exp_from_doubles({{1, 0, 0, 0}}));
    asy.emplace_back(std::move(mum1), exp_from_doubles({{4, 0, 0, 0}}));

    auto out = ode::union_rule_set(asy);
    ASSERT_EQ(out.size(), 2u);

    bool found_A = false, found_B = false;
    for (const auto& t : out) {
        if (acb_close_to_double(t.mu, -1.0, 0.0)) {
            found_A = true;
            EXPECT_TRUE(acb_close_to_double(t.exp[0][0], 4.0, 0.0));
            EXPECT_TRUE(acb_close_to_double(t.exp[0][1], 1.0, 0.0));
        }
        if (acb_close_to_double(t.mu, 0.5, 0.0)) {
            found_B = true;
            EXPECT_TRUE(acb_close_to_double(t.exp[0][0], 1.0, 0.0));
        }
    }
    EXPECT_TRUE(found_A);
    EXPECT_TRUE(found_B);
}

// =========================================================================
//  ps_times_rule_set
// =========================================================================

TEST_F(AsyTest, PSTimesRuleSet_MultipliesEachTerm) {
    nm::GlobalScope s; s.expansion.x_order = 3; s.commit();

    ode::PowerSeries ps;
    ps.leading_offset = 1;
    ps.coeffs.push_back(acb_si(2));   // ps = 2 * eta^1

    ode::AsyExpansion asy;
    nm::AcbValue mu0; mu0.set_si(0);
    asy.emplace_back(std::move(mu0), exp_from_doubles({{3, 0, 0, 0}}));

    auto out = ode::ps_times_rule_set(ps, asy);
    ASSERT_EQ(out.size(), 1u);
    EXPECT_TRUE(acb_close_to_double(out[0].mu, 1.0, 0.0));
    EXPECT_TRUE(acb_close_to_double(out[0].exp[0][0], 6.0, 0.0));
}

TEST_F(AsyTest, PSTimesRuleSet_ZeroPowerSeries_GivesEmpty) {
    ode::PowerSeries ps;
    ode::AsyExpansion asy;
    nm::AcbValue mu0; mu0.set_si(0);
    asy.emplace_back(std::move(mu0), exp_from_doubles({{1}}));
    auto out = ode::ps_times_rule_set(ps, asy);
    EXPECT_TRUE(out.empty());
}

TEST_F(AsyTest, PSTimesRuleSet_MultipleCoeffsCrossProduct) {
    nm::GlobalScope s; s.expansion.x_order = 3; s.commit();

    // ps = 1 + 2*eta (offset 0, coeffs {1, 2})
    // asy = single term mu=0, exp[0][0]=1
    // expected: { (mu=0, c=1), (mu=1, c=2) }
    ode::PowerSeries ps;
    ps.leading_offset = 0;
    ps.coeffs.push_back(acb_si(1));
    ps.coeffs.push_back(acb_si(2));

    ode::AsyExpansion asy;
    nm::AcbValue mu0; mu0.set_si(0);
    asy.emplace_back(std::move(mu0), exp_from_doubles({{1, 0, 0, 0}}));

    auto out = ode::ps_times_rule_set(ps, asy);
    ASSERT_EQ(out.size(), 2u);
    EXPECT_TRUE(acb_close_to_double(out[0].mu, 0.0, 0.0));
    EXPECT_TRUE(acb_close_to_double(out[0].exp[0][0], 1.0, 0.0));
    EXPECT_TRUE(acb_close_to_double(out[1].mu, 1.0, 0.0));
    EXPECT_TRUE(acb_close_to_double(out[1].exp[0][0], 2.0, 0.0));
}

// =========================================================================
//  times_asy_exp / plus_asy_exp
// =========================================================================

TEST_F(AsyTest, TimesAsyExp_PolynomialFactor) {
    nm::GlobalScope s; s.expansion.x_order = 3; s.commit();   // row length = 4

    // rational = 1 + 2 eta;  asy = { mu = 0, exp = {{1, 0, 0, 0}} }.
    auto rat = nm::RationalFunction::from_polynomial(make_poly({1, 2}));
    ode::AsyExpansion asy;
    nm::AcbValue mu0; mu0.set_si(0);
    asy.emplace_back(std::move(mu0), exp_from_doubles({{1, 0, 0, 0}}));

    auto out = ode::times_asy_exp(rat, asy);
    ASSERT_EQ(out.size(), 1u);
    EXPECT_TRUE(acb_close_to_double(out[0].mu, 0.0, 0.0));
    ASSERT_EQ(out[0].exp.size(), 1u);
    EXPECT_TRUE(acb_close_to_double(out[0].exp[0][0], 1.0, 0.0));
    EXPECT_TRUE(acb_close_to_double(out[0].exp[0][1], 2.0, 0.0));
}

TEST_F(AsyTest, TimesAsyExp_EtaPowerShiftsMu) {
    nm::GlobalScope s; s.expansion.x_order = 3; s.commit();

    // rational = 1 / eta^2 -> shift mu by -2; expansion unchanged
    auto rat = nm::RationalFunction::monomial(-2);
    ode::AsyExpansion asy;
    nm::AcbValue mu0; mu0.set_si(0);
    asy.emplace_back(std::move(mu0), exp_from_doubles({{1, 0, 0, 0}}));
    auto out = ode::times_asy_exp(rat, asy);
    ASSERT_EQ(out.size(), 1u);
    EXPECT_TRUE(acb_close_to_double(out[0].mu, -2.0, 0.0));
    EXPECT_TRUE(acb_close_to_double(out[0].exp[0][0], 1.0, 0.0));
}

TEST_F(AsyTest, TimesAsyExp_RationalDenominatorEvaluatesGeometricSeries) {
    nm::GlobalScope s; s.expansion.x_order = 3; s.commit();

    // rational = 1/(1 - eta); asy = { mu = 0, exp = {{1, 0, 0, 0}} }
    // expected: out[0].exp = {{1, 1, 1, 1}}, mu unchanged
    auto rat = nm::RationalFunction(make_poly({1}), make_poly({1, -1}));
    ode::AsyExpansion asy;
    nm::AcbValue mu0; mu0.set_si(0);
    asy.emplace_back(std::move(mu0), exp_from_doubles({{1, 0, 0, 0}}));
    auto out = ode::times_asy_exp(rat, asy);
    ASSERT_EQ(out.size(), 1u);
    EXPECT_TRUE(acb_close_to_double(out[0].mu, 0.0, 0.0));
    for (std::size_t k = 0; k < 4; ++k) {
        EXPECT_TRUE(acb_close_to_double(out[0].exp[0][k], 1.0, 0.0));
    }
}

TEST_F(AsyTest, PlusAsyExp_ConcatenatesAndUnions) {
    nm::GlobalScope s; s.expansion.x_order = 3; s.commit();

    ode::AsyExpansion a, b;
    nm::AcbValue mu_a, mu_b;
    mu_a.set_si(0);
    mu_b.set_si(1);
    a.emplace_back(std::move(mu_a), exp_from_doubles({{1, 0, 0, 0}}));
    b.emplace_back(std::move(mu_b), exp_from_doubles({{5, 0, 0, 0}}));

    std::vector<ode::AsyExpansion> list;
    list.push_back(std::move(a));
    list.push_back(std::move(b));

    auto out = ode::plus_asy_exp(list);
    ASSERT_EQ(out.size(), 1u);
    EXPECT_TRUE(acb_close_to_double(out[0].mu, 0.0, 0.0));
    EXPECT_TRUE(acb_close_to_double(out[0].exp[0][0], 1.0, 0.0));
    EXPECT_TRUE(acb_close_to_double(out[0].exp[0][1], 5.0, 0.0));
}

// =========================================================================
//  find_log_power
// =========================================================================

TEST_F(AsyTest, FindLogPower_DropsTrailingZeroRows) {
    auto exp = exp_from_doubles({{1, 0}, {0, 0}, {0, 0}});
    EXPECT_EQ(ode::find_log_power(exp), 0L);

    auto exp2 = exp_from_doubles({{1, 0}, {0, 1}, {0, 0}});
    EXPECT_EQ(ode::find_log_power(exp2), 1L);

    auto exp3 = exp_from_doubles({{0, 0}, {0, 0}});
    EXPECT_EQ(ode::find_log_power(exp3), -1L);
}

// =========================================================================
//  ps_map_rule_set
// =========================================================================

TEST_F(AsyTest, PSMapRuleSet_IdentityIsAPassThrough) {
    nm::GlobalScope s; s.expansion.x_order = 3; s.commit();

    ode::PowerSeries ps_one;
    ps_one.leading_offset = 0;
    ps_one.coeffs.push_back(acb_si(1));
    std::vector<std::vector<ode::PowerSeries>> psmap(1);
    psmap[0].push_back(std::move(ps_one));

    ode::AsyExpansion asy;
    nm::AcbValue mu; mu.set_si(0);
    asy.emplace_back(std::move(mu), exp_from_doubles({{7, 0, 0, 0}}));

    std::vector<ode::AsyExpansion> asy_list;
    asy_list.push_back(std::move(asy));

    auto out = ode::ps_map_rule_set(psmap, asy_list);
    ASSERT_EQ(out.size(), 1u);
    ASSERT_EQ(out[0].size(), 1u);
    EXPECT_TRUE(acb_close_to_double(out[0][0].mu, 0.0, 0.0));
    EXPECT_TRUE(acb_close_to_double(out[0][0].exp[0][0], 7.0, 0.0));
}

TEST_F(AsyTest, PSMapRuleSet_RowDimensionMismatchThrows) {
    std::vector<std::vector<ode::PowerSeries>> psmap(1);
    psmap[0].resize(2);     // expects K=2
    std::vector<ode::AsyExpansion> asy_list;
    asy_list.resize(1);     // K=1 -> mismatch
    EXPECT_THROW(ode::ps_map_rule_set(psmap, asy_list), std::invalid_argument);
}

// =========================================================================
//  Cross-validation: times_asy_exp followed by evaluate
// =========================================================================

TEST_F(AsyTest, TimesAsyExp_EvaluationMatchesDirectComputation) {
    nm::GlobalScope s; s.expansion.x_order = 5; s.commit();

    // rational = (1 + 2 eta) / (1 - eta);  asy = { mu = 0, exp = {{1, 0, 0, 0, 0, 0}} }
    // expected truncated series = (1 + 2 eta) * (1 + eta + eta^2 + eta^3 + eta^4 + eta^5)
    //   coefficients: c0=1, c1=1+2=3, c2=1+2=3, c3=1+2=3, c4=1+2=3, c5=1+2=3
    auto rat = nm::RationalFunction(make_poly({1, 2}), make_poly({1, -1}));
    ode::AsyExpansion asy;
    nm::AcbValue mu0; mu0.set_si(0);
    asy.emplace_back(std::move(mu0), exp_from_doubles({{1, 0, 0, 0, 0, 0}}));
    auto out = ode::times_asy_exp(rat, asy);
    ASSERT_EQ(out.size(), 1u);
    ASSERT_EQ(out[0].exp[0].size(), 6u);
    EXPECT_TRUE(acb_close_to_double(out[0].exp[0][0], 1.0, 0.0));
    EXPECT_TRUE(acb_close_to_double(out[0].exp[0][1], 3.0, 0.0));
    EXPECT_TRUE(acb_close_to_double(out[0].exp[0][2], 3.0, 0.0));
    EXPECT_TRUE(acb_close_to_double(out[0].exp[0][3], 3.0, 0.0));
    EXPECT_TRUE(acb_close_to_double(out[0].exp[0][4], 3.0, 0.0));
    EXPECT_TRUE(acb_close_to_double(out[0].exp[0][5], 3.0, 0.0));
}

TEST_F(AsyTest, CloneAsy_DeepCopiesWithoutAliasing) {
    nm::GlobalScope s; s.expansion.x_order = 3; s.commit();

    ode::AsyExpansion src;
    nm::AcbValue mu; mu.set_si(2);
    src.emplace_back(std::move(mu), exp_from_doubles({{7, 0, 0, 0}}));
    auto cp = ode::clone_asy(src);
    ASSERT_EQ(cp.size(), 1u);
    EXPECT_TRUE(acb_close_to_double(cp[0].mu, 2.0, 0.0));
    // Mutate clone; original untouched.
    cp[0].mu.set_si(99);
    EXPECT_TRUE(acb_close_to_double(src[0].mu,  2.0, 0.0));
    EXPECT_TRUE(acb_close_to_double(cp[0].mu,  99.0, 0.0));
}
