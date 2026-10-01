// SPDX-License-Identifier: MIT
// Tests for amflow::ode::regular (Layer 5).

#include <gtest/gtest.h>

#include "amflow/ode/regular.hpp"
#include "amflow/ode/blocks.hpp"
#include "amflow/numeric/options.hpp"
#include "amflow/numeric/rational.hpp"

#include <flint/acb.h>
#include <flint/arb.h>
#include <flint/fmpq.h>

#include <algorithm>
#include <cmath>
#include <vector>

namespace nm  = amflow::numeric;
namespace ode = amflow::ode;

namespace {

class RegularTest : public ::testing::Test {
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

::testing::AssertionResult acb_close_relative(const nm::AcbValue& v,
                                              double r, double i,
                                              double rel_tol = 1e-12,
                                              double abs_floor = 1e-20) {
    const double got_r = arf_get_d(arb_midref(acb_realref(v.raw())), ARF_RND_NEAR);
    const double got_i = arf_get_d(arb_midref(acb_imagref(v.raw())), ARF_RND_NEAR);
    const double err_r = std::abs(got_r - r);
    const double err_i = std::abs(got_i - i);
    const double scale_r = std::max(std::abs(r), abs_floor);
    const double scale_i = std::max(std::abs(i), abs_floor);
    if (err_r <= rel_tol * scale_r && err_i <= rel_tol * scale_i) {
        return ::testing::AssertionSuccess();
    }
    return ::testing::AssertionFailure()
        << "got (" << got_r << ", " << got_i << ")"
        << " vs expected (" << r << ", " << i << ")"
        << " with relative tol " << rel_tol;
}

}  // namespace

// ============================================================================
//  expand_nh_equations_num
// ============================================================================

TEST_F(RegularTest, ExpandAtZero_IsClone) {
    nm::RationalMatrix de(1, 1);
    de(0, 0) = nm::RationalFunction::from_si(2);
    auto eqs  = ode::nh_equations(de, ode::EquationMode::Regular);
    auto neqs = ode::nh_equations_num(eqs);

    nm::AcbValue zero;
    auto out = ode::expand_nh_equations_num(neqs, zero.raw());
    ASSERT_EQ(out.size(), neqs.size());
    EXPECT_EQ(out[0].dxexp.size(), neqs[0].dxexp.size());
    EXPECT_EQ(out[0].axexp.size(), neqs[0].axexp.size());
    EXPECT_EQ(out[0].block,         neqs[0].block);
}

TEST_F(RegularTest, ExpandAtX0_ShiftsLinearTerm) {
    // de = 5/(1 + eta) -> Regular form dx = 1+eta, ax = 5.
    nm::RationalMatrix de(1, 1);
    {
        nm::FmpqPoly num = make_poly({5});
        nm::FmpqPoly den = make_poly({1, 1});
        de(0, 0) = nm::RationalFunction(std::move(num), std::move(den));
    }
    auto eqs  = ode::nh_equations(de, ode::EquationMode::Regular);
    auto neqs = ode::nh_equations_num(eqs);

    ASSERT_EQ(neqs[0].dxexp.size(), 2u);
    EXPECT_TRUE(acb_close_to_double(neqs[0].dxexp[0], 1.0, 0.0));
    EXPECT_TRUE(acb_close_to_double(neqs[0].dxexp[1], 1.0, 0.0));

    // Shift to x0 = 3:  (1 + eta) -> (4 + eta);  ax stays.
    nm::AcbValue x0; x0.set_si(3);
    auto out = ode::expand_nh_equations_num(neqs, x0.raw());
    ASSERT_EQ(out[0].dxexp.size(), 2u);
    EXPECT_TRUE(acb_close_to_double(out[0].dxexp[0], 4.0, 0.0));
    EXPECT_TRUE(acb_close_to_double(out[0].dxexp[1], 1.0, 0.0));

    ASSERT_EQ(out[0].axexp.size(), 1u);
    ASSERT_EQ(out[0].axexp[0].size(), 1u);
    ASSERT_EQ(out[0].axexp[0][0].size(), 1u);
    EXPECT_TRUE(acb_close_to_double(out[0].axexp[0][0][0], 5.0, 0.0));
}

// ============================================================================
//  evaluate_taylor
// ============================================================================

TEST_F(RegularTest, EvaluateTaylor_HornerScheme) {
    // f(eta) = 1 + 2 eta + 3 eta^2  -> at dh=4 yields 57.
    ode::TaylorCoefficients c(1);
    c[0].push_back(acb_si(1));
    c[0].push_back(acb_si(2));
    c[0].push_back(acb_si(3));
    nm::AcbValue dh; dh.set_si(4);
    auto v = ode::evaluate_taylor(c, dh.raw());
    ASSERT_EQ(v.size(), 1u);
    EXPECT_TRUE(acb_close_to_double(v[0], 1.0 + 2*4 + 3*16, 0.0));
}

// --- `evaluate_taylor` strips arb radii ----
//
// `evaluate_taylor` (`src/ode/regular.cpp`) deliberately calls
// `midpoint_only` on every intermediate Horner accumulator and on
// every input coefficient.  The source comment documents this as
// "midpoint-only Horner: mirrors Mathematica's point arithmetic
// and avoids catastrophic interval blow-up on cancellation-heavy
// regular contours."  This test locks the no-radius-propagation
// contract by feeding non-zero arb radii into a coefficient and
// asserting the output radius is exactly zero.
TEST_F(RegularTest, EvaluateTaylor_StripsArbRadii_LocksMidpointOnlyContract) {
    // Coefficient c = 1 ± 0.001 (non-zero radius).  Without the
    // midpoint-only stripping, evaluating this at dh=2 with two
    // additional zero-coefficient terms would propagate the radius
    // (Horner: ((0)*dh + 0)*dh + 1 → still has the radius from c0).
    ode::TaylorCoefficients c(1);
    {
        nm::AcbValue v;
        v.set_si(1);
        // Inflate the real-part radius via mag (1e-3).  This
        // simulates the kind of radius an upstream arb computation
        // would accumulate.
        mag_t r;
        mag_init(r);
        mag_set_d(r, 1e-3);
        arb_add_error_mag(acb_realref(v.raw()), r);
        mag_clear(r);
        c[0].push_back(std::move(v));
    }
    {
        nm::AcbValue zero;
        c[0].push_back(std::move(zero));
    }

    nm::AcbValue dh; dh.set_si(2);
    auto out = ode::evaluate_taylor(c, dh.raw());
    ASSERT_EQ(out.size(), 1u);

    // Value should be 1 (since the second coefficient is 0:
    // f(2) = c0 + c1*2 = 1 + 0).
    EXPECT_TRUE(acb_close_to_double(out[0], 1.0, 0.0, 1e-12));

    // The midpoint-only contract: radius is exactly zero on output.
    arb_t re_rad, im_rad;
    arb_init(re_rad); arb_init(im_rad);
    arb_get_rad_arb(re_rad, acb_realref(out[0].raw()));
    arb_get_rad_arb(im_rad, acb_imagref(out[0].raw()));
    EXPECT_TRUE(arb_is_zero(re_rad))
        << "midpoint-only Horner must strip the input radius";
    EXPECT_TRUE(arb_is_zero(im_rad));
    arb_clear(re_rad); arb_clear(im_rad);
}

// ============================================================================
//  calcx1x2 -- f' = 2 f, f(0) = 1  ->  f_k = 2^k / k!
// ============================================================================

TEST_F(RegularTest, Calcx1x2_ExponentialODE_Coefficients) {
    nm::GlobalScope s; s.expansion.x_order = 5; s.commit();

    nm::RationalMatrix de(1, 1);
    de(0, 0) = nm::RationalFunction::from_si(2);
    auto eqs  = ode::nh_equations(de, ode::EquationMode::Regular);
    auto neqs = ode::nh_equations_num(eqs);

    std::vector<nm::AcbValue> bc; bc.push_back(acb_si(1));
    nm::AcbValue x0;  // zero

    auto coeffs = ode::calcx1x2(neqs, bc, x0.raw());
    ASSERT_EQ(coeffs.size(), 1u);
    ASSERT_EQ(coeffs[0].size(), 6u);

    double fact = 1.0;
    for (long k = 0; k <= 5; ++k) {
        if (k > 0) fact *= k;
        double expected = std::pow(2.0, k) / fact;
        EXPECT_TRUE(acb_close_to_double(coeffs[0][k], expected, 0.0))
            << "k=" << k << " got " << coeffs[0][k].to_string(20);
    }
}

TEST_F(RegularTest, Calcx1x2_FromNonZeroCenter) {
    nm::GlobalScope s; s.expansion.x_order = 5; s.commit();

    nm::RationalMatrix de(1, 1);
    de(0, 0) = nm::RationalFunction::from_si(2);
    auto eqs  = ode::nh_equations(de, ode::EquationMode::Regular);
    auto neqs = ode::nh_equations_num(eqs);

    nm::AcbValue bc_val; bc_val.set_d_d(std::exp(2.0), 0.0);
    std::vector<nm::AcbValue> bc; bc.push_back(std::move(bc_val));
    nm::AcbValue x0; x0.set_si(1);

    auto coeffs = ode::calcx1x2(neqs, bc, x0.raw());
    ASSERT_EQ(coeffs[0].size(), 6u);

    double fact = 1.0;
    for (long k = 0; k <= 5; ++k) {
        if (k > 0) fact *= k;
        double expected = std::exp(2.0) * std::pow(2.0, k) / fact;
        EXPECT_TRUE(acb_close_to_double(coeffs[0][k], expected, 0.0, 1e-12))
            << "k=" << k << " got " << coeffs[0][k].to_string(20);
    }
}

TEST_F(RegularTest, Calcx1x2_RejectsVanishingDxAtZero) {
    // Hand-craft a BlockEquationNum with dxexp[0] = 0 to verify that
    // calcx1x2 catches that and throws (would mean the matrix is *not*
    // regular at the requested centre).
    ode::BlockEquationNum eq;
    eq.dxexp.resize(2);                  // dx = 0 + 1*eta
    eq.dxexp[0].set_si(0);
    eq.dxexp[1].set_si(1);
    eq.axexp.resize(1);
    eq.axexp[0].resize(1);
    eq.axexp[0][0].push_back(acb_si(2)); // ax(0,0) = 2 (constant)
    eq.block = {0};

    std::vector<ode::BlockEquationNum> neqs;
    neqs.push_back(std::move(eq));

    std::vector<nm::AcbValue> bc; bc.push_back(acb_si(1));
    nm::AcbValue x0;  // zero
    EXPECT_THROW(ode::calcx1x2(neqs, bc, x0.raw()), std::runtime_error);
}

// ============================================================================
//  calc_run -- f' = 2 f  on a chain of points
// ============================================================================

TEST_F(RegularTest, CalcRun_ExponentialAcrossThreePoints) {
    nm::GlobalScope s;
    s.expansion.x_order  = 60;
    s.global.silent_mode = true;
    s.commit();

    nm::RationalMatrix de(1, 1);
    de(0, 0) = nm::RationalFunction::from_si(2);

    std::vector<nm::AcbValue> bc; bc.push_back(acb_si(1));

    std::vector<nm::AcbValue> run;
    run.push_back(acb_dd(0.0, 0.0));
    run.push_back(acb_dd(0.3, 0.0));
    run.push_back(acb_dd(0.7, 0.0));
    run.push_back(acb_dd(1.0, 0.0));

    auto out = ode::calc_run(de, bc, run);
    ASSERT_EQ(out.size(), 1u);
    EXPECT_TRUE(acb_close_to_double(out[0], std::exp(2.0), 0.0, 1e-10))
        << "got " << out[0].to_string(20);
}

TEST_F(RegularTest, CalcRun_EmptyRun_ReturnsEmpty) {
    nm::RationalMatrix de(1, 1);
    de(0, 0) = nm::RationalFunction::from_si(2);
    std::vector<nm::AcbValue> bc; bc.push_back(acb_si(1));
    std::vector<nm::AcbValue> run;
    auto out = ode::calc_run(de, bc, run);
    EXPECT_TRUE(out.empty());
}

TEST_F(RegularTest, CalcRun_SinglePointRun_ReturnsBC) {
    nm::RationalMatrix de(1, 1);
    de(0, 0) = nm::RationalFunction::from_si(2);
    std::vector<nm::AcbValue> bc; bc.push_back(acb_si(7));
    std::vector<nm::AcbValue> run; run.push_back(acb_si(0));
    auto out = ode::calc_run(de, bc, run);
    ASSERT_EQ(out.size(), 1u);
    EXPECT_TRUE(acb_close_to_double(out[0], 7.0, 0.0));
}

// ============================================================================
//  Coupled 2x2: harmonic oscillator at pi/2.  f' = M f with M = [[0,-1],[1,0]].
// ============================================================================

TEST_F(RegularTest, CalcRun_HarmonicOscillator_AtPiOver2) {
    nm::GlobalScope s;
    s.expansion.x_order  = 80;
    s.global.silent_mode = true;
    s.commit();

    nm::RationalMatrix de(2, 2);
    de(0, 0) = nm::RationalFunction();
    de(0, 1) = nm::RationalFunction::from_si(-1);
    de(1, 0) = nm::RationalFunction::from_si(1);
    de(1, 1) = nm::RationalFunction();

    std::vector<nm::AcbValue> bc;
    bc.push_back(acb_si(1));   // f1(0) = 1
    bc.push_back(acb_si(0));   // f2(0) = 0

    nm::AcbValue piov2;
    arb_const_pi(acb_realref(piov2.raw()), nm::working_prec_bits());
    arb_mul_2exp_si(acb_realref(piov2.raw()), acb_realref(piov2.raw()), -1);
    arb_zero(acb_imagref(piov2.raw()));

    nm::AcbValue quarter;
    arb_set(acb_realref(quarter.raw()), acb_realref(piov2.raw()));
    arb_mul_2exp_si(acb_realref(quarter.raw()), acb_realref(quarter.raw()), -1);
    arb_zero(acb_imagref(quarter.raw()));

    std::vector<nm::AcbValue> run;
    run.push_back(acb_dd(0.0, 0.0));
    {
        nm::AcbValue x;
        arb_set(acb_realref(x.raw()), acb_realref(quarter.raw()));
        arb_zero(acb_imagref(x.raw()));
        run.push_back(std::move(x));
    }
    run.push_back(std::move(piov2));

    auto out = ode::calc_run(de, bc, run);
    ASSERT_EQ(out.size(), 2u);
    EXPECT_TRUE(acb_close_to_double(out[0], 0.0, 0.0, 1e-10))
        << "f1(pi/2) = " << out[0].to_string(20);
    EXPECT_TRUE(acb_close_to_double(out[1], 1.0, 0.0, 1e-10))
        << "f2(pi/2) = " << out[1].to_string(20);
}

// ============================================================================
//  Block-coupled system:
//   f1' = (1/(eta+1)) f1                    (block 0)
//   f2' = (1/(eta+1)) f2 + (1/(eta+1)) f1   (block 1, depends on f1)
//   f1(0)=1, f2(0)=0  ->  f1=1+eta, f2=(1+eta)*log(1+eta)
// ============================================================================

TEST_F(RegularTest, CalcRun_TwoBlocks_BxCoupling) {
    nm::GlobalScope s;
    s.expansion.x_order  = 120;
    s.global.silent_mode = true;
    s.commit();

    nm::RationalMatrix de(2, 2);
    {
        nm::FmpqPoly num = make_poly({1});
        nm::FmpqPoly den = make_poly({1, 1});
        de(0, 0) = nm::RationalFunction(std::move(num), std::move(den));
    }
    de(0, 1) = nm::RationalFunction();
    {
        nm::FmpqPoly num = make_poly({1});
        nm::FmpqPoly den = make_poly({1, 1});
        de(1, 0) = nm::RationalFunction(std::move(num), std::move(den));
    }
    {
        nm::FmpqPoly num = make_poly({1});
        nm::FmpqPoly den = make_poly({1, 1});
        de(1, 1) = nm::RationalFunction(std::move(num), std::move(den));
    }

    std::vector<nm::AcbValue> bc;
    bc.push_back(acb_si(1));
    bc.push_back(acb_si(0));

    std::vector<nm::AcbValue> run;
    run.push_back(acb_dd(0.0, 0.0));
    run.push_back(acb_dd(0.25, 0.0));
    run.push_back(acb_dd(0.5, 0.0));

    auto out = ode::calc_run(de, bc, run);
    ASSERT_EQ(out.size(), 2u);
    EXPECT_TRUE(acb_close_to_double(out[0], 1.5, 0.0, 1e-10))
        << "f1 = " << out[0].to_string(20);
    EXPECT_TRUE(acb_close_to_double(out[1], 1.5 * std::log(1.5), 0.0, 1e-10))
        << "f2 = " << out[1].to_string(20);
}

// ============================================================================
//  Sunrise direct regular run regression (Mathematica reference).
// ============================================================================

TEST_F(RegularTest, CalcRun_SunriseDirectSystem_MatchesMathematica) {
    nm::GlobalScope s;
    s.expansion.x_order       = 120;
    s.expansion.extra_x_order = 20;
    s.global.silent_mode      = true;
    s.commit();

    const nm::RationalFunction x   = nm::RationalFunction::monomial(1);
    const nm::RationalFunction one = nm::RationalFunction::from_si(1);
    const nm::RationalFunction denom = (x - one) * x;

    nm::RationalMatrix de(2, 2);
    de(0, 0) = nm::RationalFunction::from_si_si(4999, 5000) / x;
    de(0, 1) = nm::RationalFunction::from_si(-1) / x;
    de(1, 0) = nm::RationalFunction::from_si_si(-99965003, 50000000) / denom;
    de(1, 1) = (nm::RationalFunction::from_si_si(19997, 10000)
                - nm::RationalFunction::from_si_si(1, 10000) * x) / denom;

    std::vector<nm::AcbValue> bc;
    bc.push_back(acb_dd(5.0007982346841556e7, 0.0));
    bc.push_back(acb_dd(4.9999230841935509e7, 0.0));

    std::vector<nm::AcbValue> run;
    run.push_back(acb_dd(0.5, 0.0));
    run.push_back(acb_dd(0.1, 0.0));

    auto out = ode::calc_run(de, bc, run);
    ASSERT_EQ(out.size(), 2u);
    EXPECT_TRUE(acb_close_relative(out[0], 5.0008982403510377e7, 0.0, 1e-9))
        << "master[0] = " << out[0].to_string(30);
    EXPECT_TRUE(acb_close_relative(out[1], 4.9999230617393328e7, 0.0, 1e-9))
        << "master[1] = " << out[1].to_string(30);
}
