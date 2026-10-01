// SPDX-License-Identifier: MIT
// Tests for amflow::ode::zero (Layer 7 calcx00 / calc_zero / pick_zero).

#include <gtest/gtest.h>

#include "amflow/ode/asy.hpp"
#include "amflow/ode/blocks.hpp"
#include "amflow/ode/normalize.hpp"
#include "amflow/ode/zero.hpp"
#include "amflow/numeric/options.hpp"
#include "amflow/numeric/rational.hpp"

#include <flint/acb.h>
#include <flint/arb.h>
#include <flint/fmpq.h>

#include <cmath>
#include <cstdlib>
#include <vector>

namespace nm  = amflow::numeric;
namespace ode = amflow::ode;

namespace {

class ZeroTest : public ::testing::Test {
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
nm::AcbValue acb_dd(double r, double i) { nm::AcbValue v; v.set_d_d(r, i); return v; }
std::vector<nm::AcbValue> acb_row_si(std::initializer_list<long> values) {
    std::vector<nm::AcbValue> out;
    out.reserve(values.size());
    for (long v : values) out.push_back(acb_si(v));
    return out;
}
nm::AcbValue acb_frac(long num, long den) {
    nm::AcbValue v;
    fmpq_t q;
    fmpq_init(q);
    fmpq_set_si(q, num, den);
    v.set_fmpq(q);
    fmpq_clear(q);
    return v;
}

void set_poly_const_fraction(nm::FmpqPoly& poly, long num, long den) {
    fmpq_t q;
    fmpq_init(q);
    fmpq_set_si(q, num, den);
    poly.set_coeff_fmpq(0, q);
    fmpq_clear(q);
}

ode::BlockEquation make_direct_block(std::size_t dim,
                                     std::vector<std::size_t> block,
                                     std::vector<std::size_t> sub = {}) {
    ode::BlockEquation eq;
    eq.dx = make_poly({1});
    eq.ax.resize(dim);
    for (auto& row : eq.ax) row.resize(dim);
    eq.bx.resize(dim);
    for (auto& row : eq.bx) row.resize(sub.size());
    eq.block = std::move(block);
    eq.sub = std::move(sub);
    return eq;
}

ode::BehaviorEntry make_behavior_si(long mu, long log_power = 0) {
    ode::BehaviorEntry entry;
    entry.mu.set_si(mu);
    entry.log_power = log_power;
    return entry;
}

ode::BehaviorEntry make_behavior_frac(long num, long den, long log_power = 0) {
    ode::BehaviorEntry entry;
    fmpq_t q;
    fmpq_init(q);
    fmpq_set_si(q, num, den);
    entry.mu.set_fmpq(q);
    fmpq_clear(q);
    entry.log_power = log_power;
    return entry;
}

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

const ode::AsyTerm* find_mu_term(const ode::AsyExpansion& asy, double mu) {
    for (const auto& term : asy) {
        if (acb_close_to_double(term.mu, mu, 0.0, 1e-12)) return &term;
    }
    return nullptr;
}

}  // namespace

// ============================================================================
//  CalcZero -- ODE  f' = 0,  f(x0) = 5 -> f(0) = 5
// ============================================================================

TEST_F(ZeroTest, CalcZero_TrivialConstant) {
    nm::GlobalScope s;
    s.expansion.x_order = 8;
    s.expansion.extra_x_order = 5;
    s.global.silent_mode = true;
    s.commit();

    nm::RationalMatrix de(1, 1);

    std::vector<nm::AcbValue> bc;
    bc.push_back(acb_si(5));
    nm::AcbValue x0; x0.set_si(2);

    auto out = ode::calc_zero(de, bc, x0.raw());
    ASSERT_EQ(out.size(), 1u);

    nm::AcbValue v = ode::pick_zero_rule_s(out[0]);
    EXPECT_TRUE(acb_close_to_double(v, 5.0, 0.0, 1e-12)) << "got " << v.to_string(20);
}

// ============================================================================
//  CalcZero -- f' = (1/eta) f, f(x0) = 7 -> f(eta) = (7/2) eta. mu = 1
// ============================================================================

TEST_F(ZeroTest, CalcZero_OneOverEta_PositiveLeadingExp) {
    nm::GlobalScope s;
    s.expansion.x_order = 10;
    s.expansion.extra_x_order = 5;
    s.global.silent_mode = true;
    s.commit();

    nm::RationalMatrix de(1, 1);
    {
        nm::FmpqPoly num = make_poly({1});
        nm::FmpqPoly den = make_poly({0, 1});
        de(0, 0) = nm::RationalFunction(std::move(num), std::move(den));
    }
    std::vector<nm::AcbValue> bc;
    bc.push_back(acb_si(7));

    nm::AcbValue x0; x0.set_si(2);

    auto out = ode::calc_zero(de, bc, x0.raw());
    ASSERT_EQ(out.size(), 1u);

    // pick at zero gives 0 (mu = 1 positive).
    nm::AcbValue v = ode::pick_zero_rule_s(out[0]);
    EXPECT_TRUE(acb_close_to_double(v, 0.0, 0.0, 1e-12)) << "got " << v.to_string(20);

    // f(0.5) = (7/2) * 0.5 = 1.75
    nm::AcbValue ev_pt; ev_pt.set_d_d(0.5, 0.0);
    nm::AcbValue val = ode::evaluate_asy_expansion(out[0], ev_pt.raw());
    EXPECT_TRUE(acb_close_to_double(val, 1.75, 0.0, 1e-10))
        << "got " << val.to_string(20);
}

// ============================================================================
//  CalcZero -- f' = -1/eta f -> f = a*x0/eta.  mu = -1, divergent.
// ============================================================================

TEST_F(ZeroTest, CalcZero_NegativeLeadingExp_DivergesAtZero) {
    nm::GlobalScope s;
    s.expansion.x_order = 10;
    s.expansion.extra_x_order = 5;
    s.global.silent_mode = true;
    s.commit();

    nm::RationalMatrix de(1, 1);
    {
        nm::FmpqPoly num = make_poly({-1});
        nm::FmpqPoly den = make_poly({0, 1});
        de(0, 0) = nm::RationalFunction(std::move(num), std::move(den));
    }
    std::vector<nm::AcbValue> bc;
    bc.push_back(acb_si(7));
    nm::AcbValue x0; x0.set_si(2);

    auto out = ode::calc_zero(de, bc, x0.raw());
    ASSERT_EQ(out.size(), 1u);
    ASSERT_EQ(out[0].size(), 1u);

    EXPECT_TRUE(acb_close_to_double(out[0][0].mu, -1.0, 0.0));
    EXPECT_TRUE(acb_close_to_double(out[0][0].exp[0][0], 14.0, 0.0, 1e-12))
        << "got " << out[0][0].exp[0][0].to_string(20);

    nm::AcbValue v = ode::pick_zero_rule_s(out[0]);
    EXPECT_TRUE(acb_close_to_double(v, 0.0, 0.0, 1e-12)) << "got " << v.to_string(20);

    // f(0.5) = 7 * 2 / 0.5 = 28
    nm::AcbValue ev_pt; ev_pt.set_d_d(0.5, 0.0);
    nm::AcbValue val = ode::evaluate_asy_expansion(out[0], ev_pt.raw());
    EXPECT_TRUE(acb_close_to_double(val, 28.0, 0.0, 1e-10))
        << "got " << val.to_string(20);
}

// ============================================================================
//  CalcZero -- Jordan log block.  f1' = (1/eta) f2, f2' = 0.
//  Exact: f2 = b, f1 = a + b log(eta/x0).
// ============================================================================

TEST_F(ZeroTest, CalcZero_JordanLogBlock_MatchesExactSolution) {
    nm::GlobalScope s;
    s.expansion.x_order = 12;
    s.expansion.extra_x_order = 6;
    s.global.silent_mode = true;
    s.commit();

    nm::RationalMatrix de(2, 2);
    {
        nm::FmpqPoly num = make_poly({1});
        nm::FmpqPoly den = make_poly({0, 1});
        de(0, 1) = nm::RationalFunction(std::move(num), std::move(den));
    }

    std::vector<nm::AcbValue> bc;
    bc.push_back(acb_si(3));
    bc.push_back(acb_si(5));

    nm::AcbValue x0; x0.set_si(2);

    auto out = ode::calc_zero(de, bc, x0.raw());
    ASSERT_EQ(out.size(), 2u);
    ASSERT_EQ(out[0].size(), 1u);
    ASSERT_EQ(out[1].size(), 1u);

    EXPECT_TRUE(acb_close_to_double(out[0][0].mu, 0.0, 0.0));
    EXPECT_GE(out[0][0].exp.size(), 2u);
    EXPECT_TRUE(acb_close_to_double(
        out[0][0].exp[0][0], 3.0 - 5.0 * std::log(2.0), 0.0, 1e-10))
        << "got " << out[0][0].exp[0][0].to_string(30);
    EXPECT_TRUE(acb_close_to_double(out[0][0].exp[1][0], 5.0, 0.0, 1e-10));

    EXPECT_TRUE(acb_close_to_double(out[1][0].mu, 0.0, 0.0));
    EXPECT_FALSE(out[1][0].exp.empty());
    EXPECT_TRUE(acb_close_to_double(out[1][0].exp[0][0], 5.0, 0.0, 1e-10));

    nm::AcbValue x_eval1; x_eval1.set_d_d(0.5, 0.0);
    nm::AcbValue v10 = ode::evaluate_asy_expansion(out[0], x_eval1.raw());
    nm::AcbValue v11 = ode::evaluate_asy_expansion(out[1], x_eval1.raw());
    EXPECT_TRUE(acb_close_to_double(v10, 3.0 + 5.0 * std::log(0.25), 0.0, 1e-10))
        << "got " << v10.to_string(30);
    EXPECT_TRUE(acb_close_to_double(v11, 5.0, 0.0, 1e-10));
}

// ============================================================================
//  Calcx00 -- direct resonant mixed block (mu = 0, mu = 1/2).
// ============================================================================

TEST_F(ZeroTest, Calcx00_ResonantMixedBlock_UsesPairConstraintsAndCompensate) {
    nm::GlobalScope s;
    s.expansion.x_order = 8;
    s.expansion.extra_x_order = 4;
    s.global.silent_mode = true;
    s.commit();

    ode::BlockEquation eq = make_direct_block(2, {0, 1});
    eq.ax[0][1] = make_poly({0, 1});
    set_poly_const_fraction(eq.ax[1][1], 1, 2);

    std::vector<ode::BlockEquation> nheq;
    nheq.push_back(std::move(eq));
    auto nheqn = ode::nh_equations_num(nheq);
    ASSERT_EQ(nheq.size(), 1u);
    ASSERT_EQ(nheqn.size(), 1u);

    ode::AsymptoticBehaviorList behavior(1);
    behavior[0].push_back(make_behavior_si(0));
    behavior[0].push_back(make_behavior_frac(1, 2));

    // Exact solution at x0 = 1/4:
    //   f2(eta) = 12 sqrt(eta), f1(eta) = 3 + 8 eta^(3/2)
    //   bc = {4, 6}.
    std::vector<nm::AcbValue> bc;
    bc.push_back(acb_si(4));
    bc.push_back(acb_si(6));
    nm::AcbValue x0; x0.set_d_d(0.25, 0.0);

    auto out = ode::calcx00(nheq, nheqn, bc, x0.raw(), behavior);
    ASSERT_EQ(out.size(), 2u);

    const ode::AsyTerm* f1_mu0 = find_mu_term(out[0], 0.0);
    const ode::AsyTerm* f1_mu_half = find_mu_term(out[0], 0.5);
    const ode::AsyTerm* f2_mu_half = find_mu_term(out[1], 0.5);

    ASSERT_NE(f1_mu0, nullptr);
    ASSERT_NE(f1_mu_half, nullptr);
    ASSERT_NE(f2_mu_half, nullptr);

    EXPECT_TRUE(acb_close_to_double(f1_mu0->exp[0][0], 3.0, 0.0, 1e-12));
    EXPECT_TRUE(acb_close_to_double(f1_mu_half->exp[0][1], 8.0, 0.0, 1e-12));
    EXPECT_TRUE(acb_close_to_double(f2_mu_half->exp[0][0], 12.0, 0.0, 1e-12));

    nm::AcbValue x_eval1; x_eval1.set_d_d(1.0 / 16.0, 0.0);
    nm::AcbValue f1_val1 = ode::evaluate_asy_expansion(out[0], x_eval1.raw());
    nm::AcbValue f2_val1 = ode::evaluate_asy_expansion(out[1], x_eval1.raw());
    EXPECT_TRUE(acb_close_to_double(f1_val1, 3.125, 0.0, 1e-10));
    EXPECT_TRUE(acb_close_to_double(f2_val1, 3.0, 0.0, 1e-10));
}

// ============================================================================
//  Calcx00 -- resonant mixed block + sub-integral input (bx coupling).
// ============================================================================

TEST_F(ZeroTest, Calcx00_ResonantMixedBlock_WithSubIntegralInputMatchesExactSolution) {
    nm::GlobalScope s;
    s.expansion.x_order = 8;
    s.expansion.extra_x_order = 4;
    s.global.silent_mode = true;
    s.commit();

    // Block 0: g' = 0 -> eta g' = 0.
    ode::BlockEquation g_eq = make_direct_block(1, {0});

    // Block 1:
    //   f2' = (1/2) f2 / eta
    //   f1' = f2 + g                              (sub from block 0)
    ode::BlockEquation f_eq = make_direct_block(2, {1, 2}, {0});
    f_eq.ax[0][1] = make_poly({0, 1});
    set_poly_const_fraction(f_eq.ax[1][1], 1, 2);
    f_eq.bx[0][0] = nm::RationalFunction::from_polynomial(make_poly({0, 1}));

    std::vector<ode::BlockEquation> nheq;
    nheq.reserve(2);
    nheq.push_back(std::move(g_eq));
    nheq.push_back(std::move(f_eq));
    auto nheqn = ode::nh_equations_num(nheq);

    ode::AsymptoticBehaviorList behavior(2);
    behavior[0].push_back(make_behavior_si(0));
    behavior[1].push_back(make_behavior_si(0));
    behavior[1].push_back(make_behavior_frac(1, 2));

    // Exact at x0 = 1/4:  g=2, f2 = 12*sqrt(eta), f1 = 3 + 2*eta + 8*eta^(3/2).
    // bc = {2, 9/2, 6}.
    std::vector<nm::AcbValue> bc;
    bc.push_back(acb_si(2));
    bc.push_back(acb_dd(4.5, 0.0));
    bc.push_back(acb_si(6));
    nm::AcbValue x0; x0.set_d_d(0.25, 0.0);

    auto out = ode::calcx00(nheq, nheqn, bc, x0.raw(), behavior);
    ASSERT_EQ(out.size(), 3u);

    const ode::AsyTerm* g_mu0 = find_mu_term(out[0], 0.0);
    const ode::AsyTerm* f1_mu0 = find_mu_term(out[1], 0.0);
    const ode::AsyTerm* f1_mu_half = find_mu_term(out[1], 0.5);
    const ode::AsyTerm* f2_mu_half = find_mu_term(out[2], 0.5);
    ASSERT_NE(g_mu0, nullptr);
    ASSERT_NE(f1_mu0, nullptr);
    ASSERT_NE(f1_mu_half, nullptr);
    ASSERT_NE(f2_mu_half, nullptr);

    EXPECT_TRUE(acb_close_to_double(g_mu0->exp[0][0], 2.0, 0.0, 1e-12));
    EXPECT_TRUE(acb_close_to_double(f1_mu0->exp[0][0], 3.0, 0.0, 1e-12));
    EXPECT_TRUE(acb_close_to_double(f1_mu0->exp[0][1], 2.0, 0.0, 1e-12));
    EXPECT_TRUE(acb_close_to_double(f1_mu_half->exp[0][1], 8.0, 0.0, 1e-12));
    EXPECT_TRUE(acb_close_to_double(f2_mu_half->exp[0][0], 12.0, 0.0, 1e-12));

    nm::AcbValue x_eval; x_eval.set_d_d(1.0 / 16.0, 0.0);
    EXPECT_TRUE(acb_close_to_double(
        ode::evaluate_asy_expansion(out[0], x_eval.raw()), 2.0, 0.0, 1e-10));
    EXPECT_TRUE(acb_close_to_double(
        ode::evaluate_asy_expansion(out[1], x_eval.raw()), 3.25, 0.0, 1e-10));
    EXPECT_TRUE(acb_close_to_double(
        ode::evaluate_asy_expansion(out[2], x_eval.raw()), 3.0, 0.0, 1e-10));
}

// ============================================================================
//  Calcx00 -- larger 3x3 resonant block (mu = 0 plus mu = 1/2).
// ============================================================================

TEST_F(ZeroTest, Calcx00_LargerResonantBlock_MatchesExactSolution) {
    nm::GlobalScope s;
    s.expansion.x_order = 8;
    s.expansion.extra_x_order = 4;
    s.global.silent_mode = true;
    s.commit();

    // Direct 3x3 block:
    //   w' = (1/2) w / eta
    //   v' = w
    //   u' = v
    ode::BlockEquation eq = make_direct_block(3, {0, 1, 2});
    eq.ax[0][1] = make_poly({0, 1});
    eq.ax[1][2] = make_poly({0, 1});
    set_poly_const_fraction(eq.ax[2][2], 1, 2);

    std::vector<ode::BlockEquation> nheq;
    nheq.push_back(std::move(eq));
    auto nheqn = ode::nh_equations_num(nheq);

    ode::AsymptoticBehaviorList behavior(1);
    behavior[0].push_back(make_behavior_si(0));
    behavior[0].push_back(make_behavior_frac(1, 2));

    // Exact at x0 = 1/4:  w = 15*sqrt(eta), v = 2 + 10*eta^(3/2),
    //                     u = 1 + 2*eta + 4*eta^(5/2).
    // bc = {13/8, 13/4, 15/2}.
    std::vector<nm::AcbValue> bc;
    bc.push_back(acb_dd(1.625, 0.0));
    bc.push_back(acb_dd(3.25, 0.0));
    bc.push_back(acb_dd(7.5, 0.0));
    nm::AcbValue x0; x0.set_d_d(0.25, 0.0);

    auto out = ode::calcx00(nheq, nheqn, bc, x0.raw(), behavior);
    ASSERT_EQ(out.size(), 3u);

    const ode::AsyTerm* u_mu0 = find_mu_term(out[0], 0.0);
    const ode::AsyTerm* u_mu_half = find_mu_term(out[0], 0.5);
    const ode::AsyTerm* v_mu0 = find_mu_term(out[1], 0.0);
    const ode::AsyTerm* v_mu_half = find_mu_term(out[1], 0.5);
    const ode::AsyTerm* w_mu_half = find_mu_term(out[2], 0.5);
    ASSERT_NE(u_mu0, nullptr);
    ASSERT_NE(u_mu_half, nullptr);
    ASSERT_NE(v_mu0, nullptr);
    ASSERT_NE(v_mu_half, nullptr);
    ASSERT_NE(w_mu_half, nullptr);

    EXPECT_TRUE(acb_close_to_double(u_mu0->exp[0][0], 1.0, 0.0, 1e-12));
    EXPECT_TRUE(acb_close_to_double(u_mu0->exp[0][1], 2.0, 0.0, 1e-12));
    EXPECT_TRUE(acb_close_to_double(u_mu_half->exp[0][2], 4.0, 0.0, 1e-12));
    EXPECT_TRUE(acb_close_to_double(v_mu0->exp[0][0], 2.0, 0.0, 1e-12));
    EXPECT_TRUE(acb_close_to_double(v_mu_half->exp[0][1], 10.0, 0.0, 1e-12));
    EXPECT_TRUE(acb_close_to_double(w_mu_half->exp[0][0], 15.0, 0.0, 1e-12));

    nm::AcbValue x_eval; x_eval.set_d_d(1.0 / 16.0, 0.0);
    EXPECT_TRUE(acb_close_to_double(
        ode::evaluate_asy_expansion(out[0], x_eval.raw()),
        1.0 + 2.0 / 16.0 + 4.0 / 1024.0, 0.0, 1e-10));
    EXPECT_TRUE(acb_close_to_double(
        ode::evaluate_asy_expansion(out[1], x_eval.raw()),
        2.0 + 10.0 / 64.0, 0.0, 1e-10));
    EXPECT_TRUE(acb_close_to_double(
        ode::evaluate_asy_expansion(out[2], x_eval.raw()),
        15.0 / 4.0, 0.0, 1e-10));
}

// ============================================================================
//  Calcx00 -- multi-resonance block (integer-spaced exponents).
//
//  This is a 5-point production top-block structure in miniature.  A_00 = diag(0,1,1),
//  so behavior mu=0 resonates at n=0 (posi={0}) AND n=1 (posi={1,2}), and
//  behavior mu=1 resonates at n=0 (posi={1,2}).  Pre-patch, the second
//  resonance of mu=0 overwrote the first's pairing data and reset
//  kk[mu=0][0][1] = I, leaving n_pair=2 < n_var-Nblk=3 and throwing
//  "boundary matching needs more match equations than block rows".
//
//  System:
//    eta u' =          eta v        ->  u' = v
//    eta v' =   v                   ->  v  = c2 eta
//    eta w' =       w               ->  w  = c3 eta
//  so  v = c2 eta,  w = c3 eta,  u = c1 + (c2/2) eta^2.
//
//  At x0 = 1/4 with (c1,c2,c3) = (3,8,5):  bc = (3+8/32, 8/4, 5/4)
//                                              = (3.25, 2.0, 1.25).
// ============================================================================

TEST_F(ZeroTest, Calcx00_MultiResonanceIntegerSpacedBlock_MatchesExactSolution) {
    nm::GlobalScope s;
    s.expansion.x_order = 8;
    s.expansion.extra_x_order = 4;
    s.global.silent_mode = true;
    s.commit();

    ode::BlockEquation eq = make_direct_block(3, {0, 1, 2});
    eq.ax[0][1] = make_poly({0, 1});   // eta * v feeds into u'
    eq.ax[1][1] = make_poly({1});      // diag entry 1
    eq.ax[2][2] = make_poly({1});      // diag entry 1

    std::vector<ode::BlockEquation> nheq;
    nheq.push_back(std::move(eq));
    auto nheqn = ode::nh_equations_num(nheq);

    ode::AsymptoticBehaviorList behavior(1);
    behavior[0].push_back(make_behavior_si(0));
    behavior[0].push_back(make_behavior_si(1));

    std::vector<nm::AcbValue> bc;
    bc.push_back(acb_dd(3.25, 0.0));
    bc.push_back(acb_dd(2.0,  0.0));
    bc.push_back(acb_dd(1.25, 0.0));
    nm::AcbValue x0; x0.set_d_d(0.25, 0.0);

    // Pre-patch this would throw "boundary matching needs more match
    // equations than block rows"; now it must complete and give the exact
    // closed form.
    ode::AsymptoticBehaviorList behavior_clone(1);
    behavior_clone[0].push_back(make_behavior_si(0));
    behavior_clone[0].push_back(make_behavior_si(1));
    auto out = ode::calcx00(nheq, nheqn, bc, x0.raw(), behavior_clone);
    ASSERT_EQ(out.size(), 3u);

    const ode::AsyTerm* u_mu0 = find_mu_term(out[0], 0.0);
    const ode::AsyTerm* u_mu1 = find_mu_term(out[0], 1.0);
    const ode::AsyTerm* v_mu1 = find_mu_term(out[1], 1.0);
    const ode::AsyTerm* w_mu1 = find_mu_term(out[2], 1.0);
    ASSERT_NE(u_mu0, nullptr);
    ASSERT_NE(v_mu1, nullptr);
    ASSERT_NE(w_mu1, nullptr);

    // u: c1 at mu=0 n=0; (c2/2) at total power 2 (either mu=0 n=2 or mu=1 n=1).
    EXPECT_TRUE(acb_close_to_double(u_mu0->exp[0][0], 3.0, 0.0, 1e-10));
    double u_eta2 =
        (u_mu0->exp[0].size() > 2
             ? arf_get_d(arb_midref(acb_realref(u_mu0->exp[0][2].raw())),
                         ARF_RND_NEAR)
             : 0.0)
        + (u_mu1 != nullptr && u_mu1->exp[0].size() > 1
             ? arf_get_d(arb_midref(acb_realref(u_mu1->exp[0][1].raw())),
                         ARF_RND_NEAR)
             : 0.0);
    EXPECT_NEAR(u_eta2, 4.0, 1e-10);  // c2/2 = 4
    // v: c2 at mu=1 n=0.
    EXPECT_TRUE(acb_close_to_double(v_mu1->exp[0][0], 8.0, 0.0, 1e-10));
    // w: c3 at mu=1 n=0.
    EXPECT_TRUE(acb_close_to_double(w_mu1->exp[0][0], 5.0, 0.0, 1e-10));

    // Round-trip evaluation at x0 must reproduce the BC.
    EXPECT_TRUE(acb_close_to_double(
        ode::evaluate_asy_expansion(out[0], x0.raw()), 3.25, 0.0, 1e-10));
    EXPECT_TRUE(acb_close_to_double(
        ode::evaluate_asy_expansion(out[1], x0.raw()), 2.0, 0.0, 1e-10));
    EXPECT_TRUE(acb_close_to_double(
        ode::evaluate_asy_expansion(out[2], x0.raw()), 1.25, 0.0, 1e-10));

    // And at a SECOND point (x=1/16) the closed form must match too --
    // this guards against the pre-patch chain corruption (kk reset at n=1)
    // which can accidentally satisfy the BC at x0 while diverging elsewhere.
    nm::AcbValue x_eval; x_eval.set_d_d(1.0 / 16.0, 0.0);
    EXPECT_TRUE(acb_close_to_double(
        ode::evaluate_asy_expansion(out[0], x_eval.raw()),
        3.0 + 4.0 / 256.0, 0.0, 1e-10));
    EXPECT_TRUE(acb_close_to_double(
        ode::evaluate_asy_expansion(out[1], x_eval.raw()),
        8.0 / 16.0, 0.0, 1e-10));
    EXPECT_TRUE(acb_close_to_double(
        ode::evaluate_asy_expansion(out[2], x_eval.raw()),
        5.0 / 16.0, 0.0, 1e-10));
}

// ============================================================================
//  Calcx00 -- Jordan/log auto-rerun.  Same 2x2 nilpotent system as
//  CalcZero_JordanLogBlock (eta u' = v, eta v' = 0; exact u=c1+c2*log eta,
//  v=c2) but the supplied behavior DELIBERATELY under-reports log_power=0.
//  At logk=0 the post-solve log^1 residue at the resonant row is
//  kk_logp1[0,:].sol = [0,1].[c1,c2] = c2 != 0, which must trigger the
//  auto-rerun bump to logk=1 and produce the correct log term instead of a
//  silently truncated (log-free) expansion.
// ============================================================================

TEST_F(ZeroTest, Calcx00_ForcedLogTerm_AutoRerunBumpsLogPower) {
    nm::GlobalScope s;
    s.expansion.x_order = 8;
    s.expansion.extra_x_order = 4;
    s.global.silent_mode = true;
    s.commit();

    ode::BlockEquation eq = make_direct_block(2, {0, 1});
    eq.ax[0][1] = make_poly({1});   // A_00 = [[0,1],[0,0]]

    std::vector<ode::BlockEquation> nheq;
    nheq.push_back(std::move(eq));
    auto nheqn = ode::nh_equations_num(nheq);

    // WRONG on purpose: true log_power is 1 (Jordan block), supply 0.
    ode::AsymptoticBehaviorList behavior(1);
    behavior[0].push_back(make_behavior_si(0, /*log_power=*/0));

    // Exact at x0=1/4: c1=2, c2=3 -> u = 2 + 3*log(1/4), v = 3.
    const double log_x0 = std::log(0.25);
    std::vector<nm::AcbValue> bc;
    bc.push_back(acb_dd(2.0 + 3.0 * log_x0, 0.0));
    bc.push_back(acb_si(3));
    nm::AcbValue x0; x0.set_d_d(0.25, 0.0);

    // Capture the auto-rerun diagnostic (gated on AMF_DBG_X00 under
    // silent_mode).
    ::setenv("AMF_DBG_X00", "1", 1);
    testing::internal::CaptureStderr();
    std::vector<ode::AsyExpansion> out;
    ASSERT_NO_THROW(out = ode::calcx00(nheq, nheqn, bc, x0.raw(), behavior));
    std::string elog = testing::internal::GetCapturedStderr();
    ::unsetenv("AMF_DBG_X00");

    // The auto-rerun must have fired.
    EXPECT_NE(elog.find("bumping log_power"), std::string::npos)
        << "auto-rerun diagnostic not seen; stderr was:\n" << elog;

    ASSERT_EQ(out.size(), 2u);
    const ode::AsyTerm* f1_mu0 = find_mu_term(out[0], 0.0);
    const ode::AsyTerm* f2_mu0 = find_mu_term(out[1], 0.0);
    ASSERT_NE(f1_mu0, nullptr);
    ASSERT_NE(f2_mu0, nullptr);

    // log_power must have been bumped to >=1: exp now has at least 2 log-rows.
    ASSERT_GE(f1_mu0->exp.size(), 2u)
        << "log tower not extended; auto-rerun did not bump log_power";
    EXPECT_TRUE(acb_close_to_double(f1_mu0->exp[0][0], 2.0, 0.0, 1e-10))
        << "got " << f1_mu0->exp[0][0].to_string(30);
    EXPECT_TRUE(acb_close_to_double(f1_mu0->exp[1][0], 3.0, 0.0, 1e-10))
        << "got " << f1_mu0->exp[1][0].to_string(30);
    EXPECT_TRUE(acb_close_to_double(f2_mu0->exp[0][0], 3.0, 0.0, 1e-10));

    // Round-trip at x0 reproduces the BC, and at a second point matches the
    // closed form (guards against a fit that satisfies x0 only).
    EXPECT_TRUE(acb_close_to_double(
        ode::evaluate_asy_expansion(out[0], x0.raw()),
        2.0 + 3.0 * log_x0, 0.0, 1e-10));
    EXPECT_TRUE(acb_close_to_double(
        ode::evaluate_asy_expansion(out[1], x0.raw()), 3.0, 0.0, 1e-10));

    nm::AcbValue x_eval; x_eval.set_d_d(1.0 / 16.0, 0.0);
    EXPECT_TRUE(acb_close_to_double(
        ode::evaluate_asy_expansion(out[0], x_eval.raw()),
        2.0 + 3.0 * std::log(1.0 / 16.0), 0.0, 1e-10));
    EXPECT_TRUE(acb_close_to_double(
        ode::evaluate_asy_expansion(out[1], x_eval.raw()), 3.0, 0.0, 1e-10));
}

// ============================================================================
//  Calcx00 -- log row + sub-integral coupling
// ============================================================================

TEST_F(ZeroTest, Calcx00_LogRowAndSubIntegralInput_MatchExactSolution) {
    nm::GlobalScope s;
    s.expansion.x_order = 8;
    s.expansion.extra_x_order = 4;
    s.global.silent_mode = true;
    s.commit();

    // Block 0: g' = 0.
    ode::BlockEquation g_eq = make_direct_block(1, {0});

    // Block 1:  w' = w/(2 eta);  v' = w;  u' = v/eta + g.
    // Singular form:
    //   eta u' = v + eta g
    //   eta v' = eta w
    //   eta w' = (1/2) w
    ode::BlockEquation f_eq = make_direct_block(3, {1, 2, 3}, {0});
    f_eq.ax[0][1] = make_poly({1});
    f_eq.ax[1][2] = make_poly({0, 1});
    set_poly_const_fraction(f_eq.ax[2][2], 1, 2);
    f_eq.bx[0][0] = nm::RationalFunction::from_polynomial(make_poly({0, 1}));

    std::vector<ode::BlockEquation> nheq;
    nheq.reserve(2);
    nheq.push_back(std::move(g_eq));
    nheq.push_back(std::move(f_eq));
    auto nheqn = ode::nh_equations_num(nheq);

    ode::AsymptoticBehaviorList behavior(2);
    behavior[0].push_back(make_behavior_si(0));
    behavior[1].push_back(make_behavior_si(0, 1));
    behavior[1].push_back(make_behavior_frac(1, 2));

    // Exact at x0 = 1/4:
    //   g = 2, w = 9*sqrt(eta), v = 5 + 6*eta^(3/2),
    //   u = 3 + 5*log(eta/(1/4)) + 4*eta^(3/2) + 2*eta
    // bc = {2, 4, 23/4, 9/2}.
    std::vector<nm::AcbValue> bc;
    bc.push_back(acb_si(2));
    bc.push_back(acb_si(4));
    bc.push_back(acb_frac(23, 4));
    bc.push_back(acb_frac(9, 2));
    nm::AcbValue x0; x0.set_d_d(0.25, 0.0);

    auto out = ode::calcx00(nheq, nheqn, bc, x0.raw(), behavior);
    ASSERT_EQ(out.size(), 4u);

    const ode::AsyTerm* g_mu0 = find_mu_term(out[0], 0.0);
    const ode::AsyTerm* u_mu0 = find_mu_term(out[1], 0.0);
    const ode::AsyTerm* u_mu_half = find_mu_term(out[1], 0.5);
    const ode::AsyTerm* v_mu0 = find_mu_term(out[2], 0.0);
    const ode::AsyTerm* v_mu_half = find_mu_term(out[2], 0.5);
    const ode::AsyTerm* w_mu_half = find_mu_term(out[3], 0.5);
    ASSERT_NE(g_mu0, nullptr);
    ASSERT_NE(u_mu0, nullptr);
    ASSERT_NE(u_mu_half, nullptr);
    ASSERT_NE(v_mu0, nullptr);
    ASSERT_NE(v_mu_half, nullptr);
    ASSERT_NE(w_mu_half, nullptr);

    EXPECT_TRUE(acb_close_to_double(g_mu0->exp[0][0], 2.0, 0.0, 1e-12));
    ASSERT_EQ(u_mu0->exp.size(), 2u);
    EXPECT_TRUE(acb_close_to_double(u_mu0->exp[0][0], 3.0 + 5.0 * std::log(4.0), 0.0, 1e-12));
    EXPECT_TRUE(acb_close_to_double(u_mu0->exp[1][0], 5.0, 0.0, 1e-12));
    EXPECT_TRUE(acb_close_to_double(u_mu0->exp[0][1], 2.0, 0.0, 1e-12));
    EXPECT_TRUE(acb_close_to_double(u_mu_half->exp[0][1], 4.0, 0.0, 1e-12));
    EXPECT_TRUE(acb_close_to_double(v_mu0->exp[0][0], 5.0, 0.0, 1e-12));
    EXPECT_TRUE(acb_close_to_double(v_mu_half->exp[0][1], 6.0, 0.0, 1e-12));
    EXPECT_TRUE(acb_close_to_double(w_mu_half->exp[0][0], 9.0, 0.0, 1e-12));

    nm::AcbValue x_eval; x_eval.set_d_d(1.0 / 16.0, 0.0);
    EXPECT_TRUE(acb_close_to_double(
        ode::evaluate_asy_expansion(out[1], x_eval.raw()),
        3.0 + 5.0 * std::log(0.25) + 4.0 / 64.0 + 2.0 / 16.0, 0.0, 1e-10));
    EXPECT_TRUE(acb_close_to_double(
        ode::evaluate_asy_expansion(out[2], x_eval.raw()),
        5.0 + 6.0 / 64.0, 0.0, 1e-10));
    EXPECT_TRUE(acb_close_to_double(
        ode::evaluate_asy_expansion(out[3], x_eval.raw()),
        9.0 / 4.0, 0.0, 1e-10));
}

// ============================================================================
//  Calcx00 -- second-order log row (3x3 nilpotent Jordan, log^2 term).
// ============================================================================

TEST_F(ZeroTest, Calcx00_SecondOrderLogRow_HasDirectNonzeroLogSquaredTerm) {
    nm::GlobalScope s;
    s.expansion.x_order = 8;
    s.expansion.extra_x_order = 4;
    s.global.silent_mode = true;
    s.commit();

    // Direct: f3' = 0;  f2' = f3/eta;  f1' = f2/eta.
    // Solution: f1 = a + b*log(eta/x0) + (c/2)*log^2(eta/x0),
    //           f2 = b + c*log(eta/x0), f3 = c.
    ode::BlockEquation eq = make_direct_block(3, {0, 1, 2});
    eq.ax[0][1] = make_poly({1});
    eq.ax[1][2] = make_poly({1});

    std::vector<ode::BlockEquation> nheq;
    nheq.push_back(std::move(eq));
    auto nheqn = ode::nh_equations_num(nheq);

    ode::AsymptoticBehaviorList behavior(1);
    behavior[0].push_back(make_behavior_si(0, 2));

    // bc = {3, 5, 7} at x0 = 2.
    std::vector<nm::AcbValue> bc;
    bc.push_back(acb_si(3));
    bc.push_back(acb_si(5));
    bc.push_back(acb_si(7));
    nm::AcbValue x0; x0.set_si(2);

    auto out = ode::calcx00(nheq, nheqn, bc, x0.raw(), behavior);
    ASSERT_EQ(out.size(), 3u);

    const ode::AsyTerm* f1_mu0 = find_mu_term(out[0], 0.0);
    const ode::AsyTerm* f2_mu0 = find_mu_term(out[1], 0.0);
    const ode::AsyTerm* f3_mu0 = find_mu_term(out[2], 0.0);
    ASSERT_NE(f1_mu0, nullptr);
    ASSERT_NE(f2_mu0, nullptr);
    ASSERT_NE(f3_mu0, nullptr);

    const double log2 = std::log(2.0);
    ASSERT_EQ(f1_mu0->exp.size(), 3u);
    EXPECT_TRUE(acb_close_to_double(
        f1_mu0->exp[0][0], 3.0 - 5.0 * log2 + 3.5 * log2 * log2, 0.0, 1e-12));
    EXPECT_TRUE(acb_close_to_double(f1_mu0->exp[1][0], 5.0 - 7.0 * log2, 0.0, 1e-12));
    EXPECT_TRUE(acb_close_to_double(f1_mu0->exp[2][0], 3.5, 0.0, 1e-12));
    EXPECT_TRUE(acb_close_to_double(f2_mu0->exp[0][0], 5.0 - 7.0 * log2, 0.0, 1e-12));
    EXPECT_TRUE(acb_close_to_double(f2_mu0->exp[1][0], 7.0, 0.0, 1e-12));
    EXPECT_TRUE(acb_close_to_double(f2_mu0->exp[2][0], 0.0, 0.0, 1e-12));
    EXPECT_TRUE(acb_close_to_double(f3_mu0->exp[0][0], 7.0, 0.0, 1e-12));

    nm::AcbValue x_eval; x_eval.set_d_d(0.5, 0.0);
    EXPECT_TRUE(acb_close_to_double(
        ode::evaluate_asy_expansion(out[0], x_eval.raw()),
        3.0 + 5.0 * std::log(0.25) + 3.5 * std::log(0.25) * std::log(0.25),
        0.0, 1e-10));
    EXPECT_TRUE(acb_close_to_double(
        ode::evaluate_asy_expansion(out[1], x_eval.raw()),
        5.0 + 7.0 * std::log(0.25), 0.0, 1e-10));
    EXPECT_TRUE(acb_close_to_double(
        ode::evaluate_asy_expansion(out[2], x_eval.raw()),
        7.0, 0.0, 1e-10));
}

// ============================================================================
//  Calcx00 -- Jordan log block (mu = 0, log_power = 1)
// ============================================================================

TEST_F(ZeroTest, Calcx00_JordanLogBlock_HasDirectNonzeroLogRow) {
    nm::GlobalScope s;
    s.expansion.x_order = 8;
    s.expansion.extra_x_order = 4;
    s.global.silent_mode = true;
    s.commit();

    ode::BlockEquation eq = make_direct_block(2, {0, 1});
    eq.ax[0][1] = make_poly({1});

    std::vector<ode::BlockEquation> nheq;
    nheq.push_back(std::move(eq));
    auto nheqn = ode::nh_equations_num(nheq);

    ode::AsymptoticBehaviorList behavior(1);
    behavior[0].push_back(make_behavior_si(0, 1));

    // Exact at x0=2: f2=5, f1=3+5*log(eta/2).  bc = {3,5}.
    std::vector<nm::AcbValue> bc;
    bc.push_back(acb_si(3));
    bc.push_back(acb_si(5));
    nm::AcbValue x0; x0.set_si(2);

    auto out = ode::calcx00(nheq, nheqn, bc, x0.raw(), behavior);
    ASSERT_EQ(out.size(), 2u);

    const ode::AsyTerm* f1_mu0 = find_mu_term(out[0], 0.0);
    const ode::AsyTerm* f2_mu0 = find_mu_term(out[1], 0.0);
    ASSERT_NE(f1_mu0, nullptr);
    ASSERT_NE(f2_mu0, nullptr);

    ASSERT_EQ(f1_mu0->exp.size(), 2u);
    EXPECT_TRUE(acb_close_to_double(f1_mu0->exp[0][0], 3.0 - 5.0 * std::log(2.0), 0.0, 1e-12));
    EXPECT_TRUE(acb_close_to_double(f1_mu0->exp[1][0], 5.0, 0.0, 1e-12));
    EXPECT_TRUE(acb_close_to_double(f2_mu0->exp[0][0], 5.0, 0.0, 1e-12));
    EXPECT_TRUE(acb_close_to_double(f2_mu0->exp[1][0], 0.0, 0.0, 1e-12));

    nm::AcbValue x_eval; x_eval.set_d_d(0.5, 0.0);
    nm::AcbValue f1_val = ode::evaluate_asy_expansion(out[0], x_eval.raw());
    nm::AcbValue f2_val = ode::evaluate_asy_expansion(out[1], x_eval.raw());
    EXPECT_TRUE(acb_close_to_double(f1_val, 3.0 + 5.0 * std::log(0.25), 0.0, 1e-10));
    EXPECT_TRUE(acb_close_to_double(f2_val, 5.0, 0.0, 1e-10));
}

// ============================================================================
//  LearnFromRuleS / LearnFromRuleSAll
// ============================================================================

TEST_F(ZeroTest, LearnFromRuleS_PrunesTrailingZeroLogRows) {
    nm::GlobalScope s;
    s.expansion.learn_x_order = 3;
    s.expansion.test_x_order = 2;
    s.global.silent_mode = true;
    s.commit();

    ode::AsyExpansion asy;

    ode::AsyTerm mu0;
    mu0.mu.set_si(0);
    mu0.exp.resize(3);
    mu0.exp[0] = acb_row_si({1, 0, 0, 0});
    mu0.exp[1] = acb_row_si({2, 0, 0, 0});
    mu0.exp[2] = acb_row_si({0, 0, 0, 0});
    asy.push_back(std::move(mu0));

    ode::AsyTerm muhalf;
    muhalf.mu = acb_frac(1, 2);
    muhalf.exp.resize(2);
    muhalf.exp[0] = acb_row_si({3, 0, 0, 0});
    muhalf.exp[1] = acb_row_si({0, 0, 0, 0});
    asy.push_back(std::move(muhalf));

    ode::BlockBehavior learned = ode::learn_from_rule_s(asy);
    ASSERT_EQ(learned.size(), 2u);
    EXPECT_TRUE(acb_close_to_double(learned[0].mu, 0.0, 0.0));
    EXPECT_EQ(learned[0].log_power, 1L);
    EXPECT_TRUE(acb_close_to_double(learned[1].mu, 0.5, 0.0));
    EXPECT_EQ(learned[1].log_power, 0L);
}

TEST_F(ZeroTest, LearnFromRuleSAll_MergesSameMuAcrossBlockAndKeepsMaxLogPower) {
    nm::GlobalScope s;
    s.expansion.learn_x_order = 3;
    s.expansion.test_x_order = 2;
    s.global.silent_mode = true;
    s.commit();

    std::vector<ode::AsyExpansion> asy_list(3);

    ode::AsyTerm a0;
    a0.mu.set_si(0);
    a0.exp.resize(1);
    a0.exp[0] = acb_row_si({1, 0, 0, 0});
    asy_list[0].push_back(std::move(a0));

    ode::AsyTerm a1_mu0;
    a1_mu0.mu.set_si(0);
    a1_mu0.exp.resize(3);
    a1_mu0.exp[0] = acb_row_si({2, 0, 0, 0});
    a1_mu0.exp[1] = acb_row_si({4, 0, 0, 0});
    a1_mu0.exp[2] = acb_row_si({0, 0, 0, 0});
    asy_list[1].push_back(std::move(a1_mu0));

    ode::AsyTerm a1_muhalf;
    a1_muhalf.mu = acb_frac(1, 2);
    a1_muhalf.exp.resize(2);
    a1_muhalf.exp[0] = acb_row_si({5, 0, 0, 0});
    a1_muhalf.exp[1] = acb_row_si({0, 0, 0, 0});
    asy_list[1].push_back(std::move(a1_muhalf));

    ode::AsyTerm a2_mu0;
    a2_mu0.mu.set_si(0);
    a2_mu0.exp.resize(2);
    a2_mu0.exp[0] = acb_row_si({3, 0, 0, 0});
    a2_mu0.exp[1] = acb_row_si({6, 0, 0, 0});
    asy_list[2].push_back(std::move(a2_mu0));

    std::vector<std::vector<std::size_t>> blocks = {{0, 1}, {2}};
    ode::AsymptoticBehaviorList learned = ode::learn_from_rule_s_all(asy_list, blocks);
    ASSERT_EQ(learned.size(), 2u);
    ASSERT_EQ(learned[0].size(), 2u);
    ASSERT_EQ(learned[1].size(), 1u);

    EXPECT_TRUE(acb_close_to_double(learned[0][0].mu, 0.0, 0.0));
    EXPECT_EQ(learned[0][0].log_power, 1L);
    EXPECT_TRUE(acb_close_to_double(learned[0][1].mu, 0.5, 0.0));
    EXPECT_EQ(learned[0][1].log_power, 0L);

    EXPECT_TRUE(acb_close_to_double(learned[1][0].mu, 0.0, 0.0));
    EXPECT_EQ(learned[1][0].log_power, 1L);
}

// ============================================================================
//  PickZeroRuleS
// ============================================================================

TEST_F(ZeroTest, PickZeroRuleS_NoIntegerMu_ReturnsZero) {
    ode::AsyExpansion asy;
    ode::AsyTerm t;
    t.mu = acb_dd(0.5, 0.0);
    t.exp.resize(1);
    t.exp[0].push_back(acb_si(99));
    asy.push_back(std::move(t));

    nm::AcbValue v = ode::pick_zero_rule_s(asy);
    EXPECT_TRUE(v.is_zero());
}

TEST_F(ZeroTest, PickZeroRuleS_MuZero_ReturnsConstantTerm) {
    ode::AsyExpansion asy;
    ode::AsyTerm t;
    t.mu = acb_si(0);
    t.exp.resize(1);
    t.exp[0].push_back(acb_si(42));
    asy.push_back(std::move(t));

    nm::AcbValue v = ode::pick_zero_rule_s(asy);
    EXPECT_TRUE(acb_close_to_double(v, 42.0, 0.0));
}

TEST_F(ZeroTest, PickZeroRuleS_PositiveIntegerMu_ReturnsZero) {
    ode::AsyExpansion asy;
    ode::AsyTerm t;
    t.mu = acb_si(2);
    t.exp.resize(1);
    t.exp[0].push_back(acb_si(99));
    asy.push_back(std::move(t));

    nm::AcbValue v = ode::pick_zero_rule_s(asy);
    EXPECT_TRUE(v.is_zero());
}

TEST_F(ZeroTest, PickZeroRuleS_NegativeMu_ReturnsLowerCoefficient) {
    ode::AsyExpansion asy;
    ode::AsyTerm t;
    t.mu = acb_si(-2);
    t.exp.resize(1);
    t.exp[0].push_back(acb_si(11));   // c0 (eta^-2 coeff)
    t.exp[0].push_back(acb_si(13));   // c1 (eta^-1 coeff)
    t.exp[0].push_back(acb_si(99));   // c2 (eta^0  coeff)
    asy.push_back(std::move(t));

    nm::AcbValue v = ode::pick_zero_rule_s(asy);
    EXPECT_TRUE(acb_close_to_double(v, 99.0, 0.0));
}

// ============================================================================
//  CalcZero with algebraic-Jordan fallback (irrational eigenvalues already
//  in principal strip).
// ============================================================================

TEST_F(ZeroTest, CalcZero_AlgebraicSemisimpleBlockUsesLocalJordanFallback) {
    nm::GlobalScope s;
    s.global.working_pre = 120;
    s.global.chop_pre = 20;
    s.global.rationalize_pre = 40;
    s.global.silent_mode = true;
    s.expansion.x_order = 6;
    s.expansion.extra_x_order = 4;
    s.expansion.learn_x_order = -1;
    s.expansion.test_x_order = 5;
    s.commit();

    auto eta_pole = [](long num, long den) {
        nm::FmpqPoly n;
        set_poly_const_fraction(n, num, den);
        nm::FmpqPoly d = make_poly({0, 1});
        return nm::RationalFunction(std::move(n), std::move(d));
    };

    // dI/deta = (1/eta) [[0,-1/8],[1,1]] I; eigenvalues (1 ± 1/sqrt(2)) / 2,
    // both already in [0,1).  Choose the larger-eigenvalue eigenvector.
    nm::RationalMatrix de(2, 2);
    de(0, 1) = eta_pole(-1, 8);
    de(1, 0) = eta_pole(1, 1);
    de(1, 1) = eta_pole(1, 1);

    EXPECT_THROW((void)ode::normalize_mat(de), std::runtime_error);

    const double lambda_plus = 0.5 * (1.0 + 1.0 / std::sqrt(2.0));
    std::vector<nm::AcbValue> bc;
    bc.push_back(acb_dd(1.0, 0.0));
    bc.push_back(acb_dd(-8.0 * lambda_plus, 0.0));

    nm::AcbValue x0; x0.set_one();

    auto out = ode::calc_zero(de, bc, x0.raw());
    ASSERT_EQ(out.size(), 2u);

    const ode::AsyTerm* term0 = find_mu_term(out[0], lambda_plus);
    const ode::AsyTerm* term1 = find_mu_term(out[1], lambda_plus);
    ASSERT_NE(term0, nullptr);
    ASSERT_NE(term1, nullptr);

    ASSERT_EQ(term0->exp.size(), 1u);
    ASSERT_EQ(term1->exp.size(), 1u);
    ASSERT_GE(term0->exp[0].size(), 2u);
    ASSERT_GE(term1->exp[0].size(), 2u);

    EXPECT_TRUE(acb_close_to_double(term0->exp[0][0], 1.0, 0.0, 1e-10));
    EXPECT_TRUE(acb_close_to_double(term1->exp[0][0], -8.0 * lambda_plus, 0.0, 1e-10));
    EXPECT_TRUE(term0->exp[0][1].is_chop_zero(18));
    EXPECT_TRUE(term1->exp[0][1].is_chop_zero(18));
}

// ============================================================================
//  CalcZero_BubbleDirectSystem_MatchesMathematica
// ============================================================================
//
//  Direct DESolver bench against Mathematica for the bubble-style 2x2 system
//  with eps = 1/100.  Boundary at eta = -I/8 was prepared from MMA by
//  InfToRegular + RegularRun.

TEST_F(ZeroTest, CalcZero_BubbleDirectSystem_MatchesMathematica) {
    nm::GlobalScope s;
    s.global.working_pre = 100;
    s.global.chop_pre = 20;
    s.global.rationalize_pre = 20;
    s.global.silent_mode = true;
    s.expansion.x_order = 100;
    s.expansion.extra_x_order = 20;
    s.expansion.learn_x_order = -1;
    s.expansion.test_x_order = 5;
    s.commit();

    const nm::RationalFunction one = nm::RationalFunction::from_si(1);
    const nm::RationalFunction eta = nm::RationalFunction::monomial(1);
    const nm::RationalFunction eps_1_100 = nm::RationalFunction::from_si_si(1, 100);

    nm::RationalMatrix de(2, 2);
    de(0, 0) = (one - eps_1_100) / (one + eta);
    de(0, 1) = nm::RationalFunction();
    de(1, 0) =
        nm::RationalFunction::from_si_si(-99, 200) /
        (eta * eta
         + nm::RationalFunction::from_si_si(3, 4) * eta
         - nm::RationalFunction::from_si_si(1, 4));
    de(1, 1) =
        nm::RationalFunction::from_si_si(49, 100) /
        (eta - nm::RationalFunction::from_si_si(1, 4));

    std::vector<nm::AcbValue> bc;
    bc.push_back(acb_dd(100.44470243046783, -12.428747763120662));
    bc.push_back(acb_dd(100.65097129361438, 1.3510415486571369));

    nm::AcbValue x0;
    x0.set_d_d(0.0, -0.125);

    auto out = ode::calc_zero(de, bc, x0.raw());
    ASSERT_EQ(out.size(), 2u);

    nm::AcbValue v0 = ode::pick_zero_rule_s(out[0]);
    nm::AcbValue v1 = ode::pick_zero_rule_s(out[1]);
    EXPECT_TRUE(acb_close_to_double(v0, 100.43695466580869, 0.0, 1e-6))
        << "got " << v0.to_string(30);
    EXPECT_TRUE(acb_close_to_double(v1, 100.98911701865423, 1.4250286148998885, 1e-6))
        << "got " << v1.to_string(30);
}
