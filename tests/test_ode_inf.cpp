// SPDX-License-Identifier: MIT
// Tests for amflow::ode::inf (Layer 6).

#include <gtest/gtest.h>

#include "amflow/ode/asy.hpp"
#include "amflow/ode/inf.hpp"
#include "amflow/numeric/options.hpp"
#include "amflow/numeric/rational.hpp"

#include <flint/acb.h>
#include <flint/arb.h>
#include <flint/fmpq.h>

#include <cmath>
#include <vector>

namespace nm  = amflow::numeric;
namespace ode = amflow::ode;

namespace {

class InfTest : public ::testing::Test {
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

nm::FmpqPoly make_poly_vec(std::vector<long> coeffs) {
    nm::FmpqPoly p;
    for (std::size_t i = 0; i < coeffs.size(); ++i) {
        fmpq_t c; fmpq_init(c);
        fmpq_set_si(c, coeffs[i], 1);
        p.set_coeff_fmpq(static_cast<long>(i), c);
        fmpq_clear(c);
    }
    return p;
}

}  // namespace

// ============================================================================
//  ReverseBCS / UnionBCS / ReadBCS
// ============================================================================

TEST_F(InfTest, ReverseBCS_FlipsSignOfMu) {
    ode::BoundarySpec bc;
    {
        ode::BoundaryEntry e; e.mu = acb_si(2); e.value = acb_si(7);
        bc.push_back(std::move(e));
    }
    {
        ode::BoundaryEntry e; e.mu = acb_si(-1); e.value = acb_si(3);
        bc.push_back(std::move(e));
    }
    auto out = ode::reverse_bcs(bc);
    ASSERT_EQ(out.size(), 2u);
    EXPECT_TRUE(acb_close_to_double(out[0].mu, -2.0, 0.0));
    EXPECT_TRUE(acb_close_to_double(out[0].value, 7.0, 0.0));
    EXPECT_TRUE(acb_close_to_double(out[1].mu, 1.0, 0.0));
    EXPECT_TRUE(acb_close_to_double(out[1].value, 3.0, 0.0));
}

TEST_F(InfTest, UnionBCS_MergesEqualMu) {
    ode::BoundarySpec bc;
    {
        ode::BoundaryEntry e; e.mu = acb_si(1); e.value = acb_si(5);
        bc.push_back(std::move(e));
    }
    {
        ode::BoundaryEntry e; e.mu = acb_si(2); e.value = acb_si(7);
        bc.push_back(std::move(e));
    }
    {
        ode::BoundaryEntry e; e.mu = acb_si(1); e.value = acb_si(8);
        bc.push_back(std::move(e));
    }
    auto out = ode::union_bcs(bc);
    ASSERT_EQ(out.size(), 2u);

    bool found_one = false, found_two = false;
    for (const auto& e : out) {
        if (acb_close_to_double(e.mu, 1.0, 0.0)) {
            EXPECT_TRUE(acb_close_to_double(e.value, 13.0, 0.0));
            found_one = true;
        }
        if (acb_close_to_double(e.mu, 2.0, 0.0)) {
            EXPECT_TRUE(acb_close_to_double(e.value, 7.0, 0.0));
            found_two = true;
        }
    }
    EXPECT_TRUE(found_one);
    EXPECT_TRUE(found_two);
}

TEST_F(InfTest, ReadBCS_GroupsRegionAndReturnsMin) {
    // mu = 1, 2, 4, 1.5; region = 0 -> integer-equivalent set {1, 2, 4}, min = 1.
    ode::BoundarySpec bc;
    for (long m : {1, 2, 4}) {
        ode::BoundaryEntry e; e.mu = acb_si(m); e.value = acb_si(m * 10);
        bc.push_back(std::move(e));
    }
    {
        ode::BoundaryEntry e; e.mu = acb_dd(1.5, 0.0); e.value = acb_si(99);
        bc.push_back(std::move(e));
    }
    nm::AcbValue region = acb_si(0);
    auto r = ode::read_bcs(bc, region.raw());
    EXPECT_TRUE(acb_close_to_double(r.ini, 1.0, 0.0));
    ASSERT_EQ(r.shifted.size(), 3u);
    bool ok0 = false, ok1 = false, ok3 = false;
    for (const auto& [o, v] : r.shifted) {
        if (o == 0 && acb_close_to_double(v, 10.0, 0.0)) ok0 = true;
        if (o == 1 && acb_close_to_double(v, 20.0, 0.0)) ok1 = true;
        if (o == 3 && acb_close_to_double(v, 40.0, 0.0)) ok3 = true;
    }
    EXPECT_TRUE(ok0);
    EXPECT_TRUE(ok1);
    EXPECT_TRUE(ok3);
}

TEST_F(InfTest, ReadBCS_NoMatchYieldsRegionWithEmptyList) {
    ode::BoundarySpec bc;
    {
        ode::BoundaryEntry e; e.mu = acb_dd(0.5, 0.0); e.value = acb_si(7);
        bc.push_back(std::move(e));
    }
    nm::AcbValue region = acb_si(0);
    auto r = ode::read_bcs(bc, region.raw());
    EXPECT_TRUE(acb_close_to_double(r.ini, 0.0, 0.0));
    EXPECT_TRUE(r.shifted.empty());
}

TEST_F(InfTest, CloneBoundary_DeepCopiesIndependently) {
    ode::BoundarySpec bc;
    ode::BoundaryEntry e; e.mu = acb_si(3); e.value = acb_si(11);
    bc.push_back(std::move(e));

    auto cp = ode::clone_boundary(bc);
    cp[0].mu.set_si(99);
    EXPECT_TRUE(acb_close_to_double(bc[0].mu, 3.0, 0.0));
    EXPECT_TRUE(acb_close_to_double(cp[0].mu, 99.0, 0.0));
}

// ============================================================================
//  CalcInf -- 1D ODE  I'(eta) = -I(eta) / eta^2
//
//  Substitute x = 1/eta:  dI/dx = I, so I(x) = I(0) * exp(x).
//  Boundary at eta=infinity: I -> 1 -> mu = 0, value = 1.
//  Expected: single mu=0 term with coefficients[0][n] = 1/n!.
// ============================================================================

TEST_F(InfTest, CalcInf_OneDimensional_ExpInverseEta) {
    nm::GlobalScope s;
    s.expansion.x_order        = 30;
    s.expansion.extra_x_order  = 20;
    s.global.silent_mode       = true;
    s.commit();

    nm::RationalMatrix de(1, 1);
    {
        nm::FmpqPoly num = make_poly({-1});
        nm::FmpqPoly den = make_poly({0, 0, 1});  // eta^2
        de(0, 0) = nm::RationalFunction(std::move(num), std::move(den));
    }

    std::vector<ode::BoundarySpec> bcs(1);
    {
        ode::BoundaryEntry e; e.mu = acb_si(0); e.value = acb_si(1);
        bcs[0].push_back(std::move(e));
    }

    auto out = ode::calc_inf(de, bcs);
    ASSERT_EQ(out.size(), 1u);
    ASSERT_EQ(out[0].size(), 1u);
    const auto& term = out[0][0];
    EXPECT_TRUE(acb_close_to_double(term.mu, 0.0, 0.0));
    ASSERT_EQ(term.exp.size(), 1u);
    ASSERT_GE(term.exp[0].size(), 5u);

    double fact = 1.0;
    for (long n = 0; n <= 5; ++n) {
        if (n > 0) fact *= n;
        double expected = 1.0 / fact;
        EXPECT_TRUE(acb_close_to_double(term.exp[0][n], expected, 0.0, 1e-12))
            << "n=" << n << " got " << term.exp[0][n].to_string(20);
    }

    // Evaluate at x = 0.1 (eta = 10) -> exp(0.1).
    nm::AcbValue x_eval; x_eval.set_d_d(0.1, 0.0);
    nm::AcbValue v = ode::evaluate_asy_expansion(out[0], x_eval.raw());
    EXPECT_TRUE(acb_close_to_double(v, std::exp(0.1), 0.0, 1e-10))
        << "got " << v.to_string(20);
}

// ============================================================================
//  CalcInf -- ODE with nontrivial leading exponent.
//  I'(eta) = (a/eta) * I(eta), a = 2;  I = C * eta^a; boundary mu=a, value=1.
//  Expected: ini = -a, single term with leading coeff 1, rest zero.
// ============================================================================

TEST_F(InfTest, CalcInf_OneDimensional_PowerLaw) {
    nm::GlobalScope s;
    s.expansion.x_order        = 10;
    s.expansion.extra_x_order  = 10;
    s.global.silent_mode       = true;
    s.commit();

    const long a = 2;
    nm::RationalMatrix de(1, 1);
    {
        nm::FmpqPoly num = make_poly({a});
        nm::FmpqPoly den = make_poly({0, 1});
        de(0, 0) = nm::RationalFunction(std::move(num), std::move(den));
    }

    std::vector<ode::BoundarySpec> bcs(1);
    {
        ode::BoundaryEntry e; e.mu = acb_si(a); e.value = acb_si(1);
        bcs[0].push_back(std::move(e));
    }

    auto out = ode::calc_inf(de, bcs);
    ASSERT_EQ(out.size(), 1u);
    ASSERT_EQ(out[0].size(), 1u);
    const auto& term = out[0][0];
    EXPECT_TRUE(acb_close_to_double(term.mu, -static_cast<double>(a), 0.0));
    ASSERT_EQ(term.exp.size(), 1u);
    EXPECT_TRUE(acb_close_to_double(term.exp[0][0], 1.0, 0.0));
    for (std::size_t n = 1; n < term.exp[0].size(); ++n) {
        EXPECT_TRUE(acb_close_to_double(term.exp[0][n], 0.0, 0.0, 1e-12))
            << "n=" << n;
    }

    nm::AcbValue x_eval; x_eval.set_d_d(0.1, 0.0);
    nm::AcbValue v = ode::evaluate_asy_expansion(out[0], x_eval.raw());
    EXPECT_TRUE(acb_close_to_double(v, 100.0, 0.0, 1e-10))
        << "got " << v.to_string(20);
}

// ============================================================================
//  determine_boundary_order -- structural check
// ============================================================================

TEST_F(InfTest, DetermineBoundaryOrder_OneFreeCoefficient) {
    nm::GlobalScope s; s.expansion.extra_x_order = 10; s.commit();

    nm::RationalMatrix de_inf(1, 1);
    {
        nm::FmpqPoly num = make_poly({-2});
        nm::FmpqPoly den = make_poly({0, 1});
        de_inf(0, 0) = nm::RationalFunction(std::move(num), std::move(den));
    }
    std::vector<long> power = {-2};
    auto out = ode::determine_boundary_order(de_inf, power);
    ASSERT_EQ(out.size(), 1u);
    EXPECT_EQ(out[0], 0L);
}

TEST_F(InfTest, DetermineBoundaryOrder_ThreeMastersNegativePowers) {
    // Sanity: doesn't crash with negative powers like [0, -1, -1].
    nm::GlobalScope s; s.expansion.extra_x_order = 10; s.commit();

    nm::RationalMatrix m(3, 3);
    for (std::size_t i = 0; i < 3; ++i) {
        nm::FmpqPoly num = make_poly_vec({1});
        nm::FmpqPoly den = make_poly_vec({0, 1});
        m(i, i) = nm::RationalFunction(std::move(num), std::move(den));
    }
    m(0, 1) = nm::RationalFunction::from_si(1);
    m(0, 2) = nm::RationalFunction::from_si(1);
    m(1, 0) = nm::RationalFunction::from_si(0);
    m(1, 2) = nm::RationalFunction::from_si(0);
    m(2, 0) = nm::RationalFunction::from_si(0);
    m(2, 1) = nm::RationalFunction::from_si(0);

    std::vector<long> power = {0, -1, -1};
    auto out = ode::determine_boundary_order(m, power);
    ASSERT_EQ(out.size(), 3u);
}

TEST_F(InfTest, DetermineBoundaryOrder_RationalCoefficientsNegPower) {
    // Mimic the bubble deinf structure with negative-power request.
    // Sanity: must not crash on rational coefficients with mixed denominators
    // and a power vector containing negative entries.
    nm::GlobalScope s; s.expansion.extra_x_order = 10; s.commit();

    auto rf = [](std::vector<long> nums, long nden,
                 std::vector<long> dens, long dden) {
        nm::FmpqPoly num, den;
        for (std::size_t i = 0; i < nums.size(); ++i) {
            fmpq_t c; fmpq_init(c);
            fmpq_set_si(c, nums[i], nden);
            num.set_coeff_fmpq(static_cast<long>(i), c);
            fmpq_clear(c);
        }
        for (std::size_t i = 0; i < dens.size(); ++i) {
            fmpq_t c; fmpq_init(c);
            fmpq_set_si(c, dens[i], dden);
            den.set_coeff_fmpq(static_cast<long>(i), c);
            fmpq_clear(c);
        }
        return nm::RationalFunction(std::move(num), std::move(den));
    };

    nm::RationalMatrix m(3, 3);
    m(0, 0) = rf({0, 49, -490}, 500, {1, -10, 5}, 5);
    m(0, 1) = rf({0, 0, 99}, 250, {1, -10, 5}, 5);
    m(0, 2) = rf({0, -99, 297}, 500, {1, -9, -5, 1}, 5);
    m(1, 0) = nm::RationalFunction::from_si(0);
    m(1, 1) = nm::RationalFunction::from_si(0);
    m(1, 2) = nm::RationalFunction::from_si(0);
    m(2, 0) = nm::RationalFunction::from_si(0);
    m(2, 1) = nm::RationalFunction::from_si(0);
    m(2, 2) = rf({0, 99}, 100, {1, 1}, 1);

    std::vector<long> power = {0, -1, -1};
    auto out = ode::determine_boundary_order(m, power);
    ASSERT_EQ(out.size(), 3u);
}

TEST_F(InfTest, DetermineBoundaryOrder_ZeroDiagonalNegPower) {
    nm::GlobalScope s; s.expansion.extra_x_order = 10; s.commit();

    nm::RationalMatrix m(1, 1);
    m(0, 0) = nm::RationalFunction::from_si(0);
    std::vector<long> power = {-1};
    auto out = ode::determine_boundary_order(m, power);
    ASSERT_EQ(out.size(), 1u);
}

TEST_F(InfTest, DetermineBoundaryOrder_SunriseTopPatternMatchesMathBoundary) {
    nm::GlobalScope s;
    s.global.working_pre        = 232;
    s.global.chop_pre           = 20;
    s.global.rationalize_pre    = 20;
    s.global.silent_mode        = true;
    s.expansion.x_order         = 464;
    s.expansion.extra_x_order   = 464;
    s.expansion.learn_x_order   = -1;
    s.expansion.test_x_order    = 5;
    s.commit();

    const nm::RationalFunction eta = nm::RationalFunction::monomial(1);
    const nm::RationalFunction one = nm::RationalFunction::from_si(1);
    const nm::RationalFunction eps_1_100 = nm::RationalFunction::from_si_si(1, 100);

    nm::RationalMatrix de(3, 3);
    de(0, 0) = (nm::RationalFunction::from_si(2) * (one - eps_1_100)) / (one + eta);
    de(0, 1) = nm::RationalFunction();
    de(0, 2) = nm::RationalFunction();

    de(1, 0) =
        (nm::RationalFunction::from_si(3)
         * (nm::RationalFunction::from_si(30)
            - nm::RationalFunction::from_si(22) * eps_1_100
            - nm::RationalFunction::from_si(4) * eta
            + nm::RationalFunction::from_si(3) * eps_1_100 * eta))
        / (nm::RationalFunction::from_si(2)
           * (nm::RationalFunction::from_si(16)
              + nm::RationalFunction::from_si(32) * eta
              + nm::RationalFunction::from_si(19) * eta * eta
              + nm::RationalFunction::from_si(3) * eta * eta * eta));
    de(1, 1) =
        (nm::RationalFunction::from_si(-32)
         + nm::RationalFunction::from_si(72) * eta
         - nm::RationalFunction::from_si(84) * eps_1_100 * eta
         + nm::RationalFunction::from_si(2) * eta * eta
         - nm::RationalFunction::from_si(9) * eps_1_100 * eta * eta)
        / (nm::RationalFunction::from_si(2)
           * (nm::RationalFunction::from_si(16)
              + nm::RationalFunction::from_si(32) * eta
              + nm::RationalFunction::from_si(19) * eta * eta
              + nm::RationalFunction::from_si(3) * eta * eta * eta));
    de(1, 2) =
        (nm::RationalFunction::from_si(4)
         - nm::RationalFunction::from_si(3) * eps_1_100)
        / (nm::RationalFunction::from_si(2)
           * (nm::RationalFunction::from_si(16)
              + nm::RationalFunction::from_si(32) * eta
              + nm::RationalFunction::from_si(19) * eta * eta
              + nm::RationalFunction::from_si(3) * eta * eta * eta));

    de(2, 0) =
        (nm::RationalFunction::from_si(3)
         * (nm::RationalFunction::from_si(480)
            - nm::RationalFunction::from_si(192) * eps_1_100
            - nm::RationalFunction::from_si(152) * eta
            + nm::RationalFunction::from_si(280) * eps_1_100 * eta
            + nm::RationalFunction::from_si(168) * eta * eta
            - nm::RationalFunction::from_si(86) * eps_1_100 * eta * eta
            - nm::RationalFunction::from_si(10) * eta * eta * eta
            + nm::RationalFunction::from_si(9) * eps_1_100 * eta * eta * eta))
        / (nm::RationalFunction::from_si(2)
           * (nm::RationalFunction::from_si(4) + eta)
           * (nm::RationalFunction::from_si(4)
              + nm::RationalFunction::from_si(7) * eta
              + nm::RationalFunction::from_si(3) * eta * eta));
    de(2, 1) =
        (nm::RationalFunction::from_si(-1792)
         + nm::RationalFunction::from_si(768) * eps_1_100
         + nm::RationalFunction::from_si(256) * eta
         - nm::RationalFunction::from_si(768) * eps_1_100 * eta
         - nm::RationalFunction::from_si(464) * eta * eta
         + nm::RationalFunction::from_si(384) * eps_1_100 * eta * eta
         + nm::RationalFunction::from_si(232) * eta * eta * eta
         - nm::RationalFunction::from_si(96) * eps_1_100 * eta * eta * eta
         - nm::RationalFunction::from_si(10) * eta * eta * eta * eta
         + nm::RationalFunction::from_si(9) * eps_1_100 * eta * eta * eta * eta)
        / (nm::RationalFunction::from_si(2)
           * (nm::RationalFunction::from_si(4) + eta)
           * (nm::RationalFunction::from_si(4)
              + nm::RationalFunction::from_si(7) * eta
              + nm::RationalFunction::from_si(3) * eta * eta));
    de(2, 2) =
        (nm::RationalFunction::from_si(128)
         - nm::RationalFunction::from_si(96) * eps_1_100
         + nm::RationalFunction::from_si(48) * eta
         - nm::RationalFunction::from_si(36) * eps_1_100 * eta
         + nm::RationalFunction::from_si(28) * eta * eta
         - nm::RationalFunction::from_si(21) * eps_1_100 * eta * eta)
        / (nm::RationalFunction::from_si(2)
           * (nm::RationalFunction::from_si(4) + eta)
           * (nm::RationalFunction::from_si(4)
              + nm::RationalFunction::from_si(7) * eta
              + nm::RationalFunction::from_si(3) * eta * eta));

    nm::RationalMatrix de_inf(3, 3);
    nm::RationalFunction prefactor = nm::RationalFunction::from_si(-1);
    prefactor.multiply_by_eta_power(-2);
    for (std::size_t i = 0; i < de.rows(); ++i) {
        for (std::size_t j = 0; j < de.cols(); ++j) {
            de_inf(i, j) = de(i, j).substitute_inverse_eta() * prefactor;
        }
    }

    std::vector<nm::RationalFunction> power_q;
    power_q.push_back(nm::RationalFunction::from_si_si(-99, 50));
    power_q.push_back(nm::RationalFunction::from_si_si(-49, 50));
    power_q.push_back(nm::RationalFunction::from_si_si(-149, 50));

    const auto out = ode::determine_boundary_order(de_inf, power_q);
    ASSERT_EQ(out.size(), 3u);
    EXPECT_EQ(out[0], 0L);
    EXPECT_EQ(out[1], -1L);
    EXPECT_EQ(out[2], 0L);
}

// ============================================================================
//  Pentabox failure-pattern reproducer
//
//  Goal: replicate the pentabox-2L bug where a master with all-zero BC
//  is coupled (via DE) to other masters that have non-zero BC in
//  different fractional Frobenius regions.
//
//  Setup: 3-master Cauchy-Euler DE in eta frame
//      I_a' = (mu_a/eta) I_a                              (mu_a = -9/5)
//      I_b' = (mu_b/eta) I_b + (1/eta) I_a                (mu_b = -9/10)
//      I_c' = (mu_c/eta) I_c + (1/eta) I_b                (mu_c = 0)
//
//      A=B=1 BC for a, b; c has NO BC.
//
//  Analytical solution (homogeneous + particular):
//      I_a(eta) = eta^mu_a
//      I_b(eta) = 1/(mu_a - mu_b) * eta^mu_a + eta^mu_b
//      I_c(eta) = K_cb*K_ba/[(mu_a-mu_b)(mu_a-mu_c)] * eta^mu_a
//                 + K_cb/(mu_b-mu_c) * eta^mu_b
//                 + C * eta^mu_c   (C = 0 since c has no BC)
//
//  For mu_a = -9/5, mu_b = -9/10, mu_c = 0:
//      c_a := 1/[(mu_a-mu_b)(mu_a-mu_c)]
//           = 1/[(-9/10)(-9/5)]
//           = 1/(81/50) = 50/81
//      c_b := 1/(mu_b-mu_c) = 1/(-9/10) = -10/9
//
//  After reverse_bcs, in x = 1/eta frame:
//      a's BC at mu = 9/5 (frac 4/5)
//      b's BC at mu = 9/10 (frac 9/10)
//      c: empty
//
//  Two distinct fractional regions: 4/5 and 9/10.
//
//  Expected calc_inf output for integral c:
//      Region 4/5 (ini = 4/5):
//        J_c(x) = x^(-4/5) * (4/5-class part of I_c)
//               = x^(-4/5) * (c_a * x^(9/5))
//               = c_a * x   (order-1 coefficient)
//        => exp[0] = [0, 50/81, 0, 0, ...]
//
//      Region 9/10 (ini = 9/10):
//        J_c(x) = x^(-9/10) * (9/10-class part of I_c)
//               = x^(-9/10) * (c_b * x^(9/10))
//               = c_b   (order-0 coefficient)
//        => exp[0] = [-10/9, 0, 0, ...]
//
//  If C++ produces these correctly, calc_inf handles multi-fractional
//  + all-zero-BC + DE-coupling correctly.  Then the pentabox bug is
//  either in a higher-level pipeline (boundary computation, sub-system
//  setup) or specific to the system size.  If C++ doesn't match, this
//  is the fast reproducer for the pentabox failure.
// ============================================================================

TEST_F(InfTest, CalcInf_MultiFractional_CrossCoupling_AllZeroBC) {
    // The original pentabox-2L config (working_pre=120, rationalize_pre=100)
    // tripped a precision-mismatch bug: rationalize_pre exceeded what
    // working_pre could resolve, so the rationalization captured binary
    // representation noise as a "real" rational, leaking a 1e-25 residual
    // into m_pure's diagonal and propagating into 1e7-1e21 errors on
    // masters with all-zero BC.  Now defended in acb_real_to_fmpq_local
    // (caps rationalize_digits at working_prec * log10(2) - safety).
    nm::GlobalScope s;
    s.expansion.x_order        = 8;
    s.expansion.extra_x_order  = 8;
    s.global.silent_mode       = true;
    s.global.working_pre       = 120;   // pentabox default; defensive cap kicks in
    s.global.chop_pre          = 20;
    s.global.rationalize_pre   = 100;   // pentabox default; intentionally too high
    s.commit();

    // This hand-built 3×3 DE produces legitimate cross-coupling residuals
    // in the ~10⁻³⁰..10⁻⁴⁰ range after `sparse_gaussian` propagation; the
    // new sparse-layer default `sparse_chop_digits = max(chop_pre,
    // working_pre−40)` = 80 chops them away and the cross-coupling
    // collapses to 0.  Production matrices (oracle benches, bn3_4mass)
    // do not hit this pathology because their legitimate coefficients
    // are O(1) — the chop floor only catches Gauss-accumulated noise.
    // For the toy DE, pin the chop floor to the user's chop_pre = 20 via
    // the dedicated env var so the legitimate small intermediates survive.
    ::setenv("AMFLOW_SPARSE_CHOP_DIGITS", "20", /*overwrite=*/1);
    struct RestoreChopDigits {
        ~RestoreChopDigits() { ::unsetenv("AMFLOW_SPARSE_CHOP_DIGITS"); }
    } restore_chop;

    auto rf_si = [](long n) { return nm::RationalFunction::from_si(n); };
    auto rf_frac = [](long n, long d) { return nm::RationalFunction::from_si_si(n, d); };

    const nm::RationalFunction eta = nm::RationalFunction::monomial(1);

    nm::RationalMatrix de(3, 3);
    // de[0,0] = mu_a / eta = -9/5 / eta
    de(0, 0) = rf_frac(-9, 5) / eta;
    de(0, 1) = nm::RationalFunction();
    de(0, 2) = nm::RationalFunction();
    // de[1,0] = 1/eta
    de(1, 0) = rf_si(1) / eta;
    // de[1,1] = mu_b/eta = -9/10 / eta
    de(1, 1) = rf_frac(-9, 10) / eta;
    de(1, 2) = nm::RationalFunction();
    de(2, 0) = nm::RationalFunction();
    // de[2,1] = 1/eta
    de(2, 1) = rf_si(1) / eta;
    // de[2,2] = mu_c/eta = 0
    de(2, 2) = nm::RationalFunction();

    std::vector<ode::BoundarySpec> bcs(3);
    {
        // a: BC mu = -9/5, value = 1 (in eta frame before reverse)
        ode::BoundaryEntry e;
        e.mu.set_si(0);
        // build -9/5 via acb_set_fmpq
        fmpq_t q; fmpq_init(q); fmpq_set_si(q, -9, 5);
        acb_set_fmpq(e.mu.raw(), q, 80);
        fmpq_clear(q);
        e.value = acb_si(1);
        bcs[0].push_back(std::move(e));
    }
    {
        // b: BC mu = -9/10, value = 1
        ode::BoundaryEntry e;
        fmpq_t q; fmpq_init(q); fmpq_set_si(q, -9, 10);
        acb_set_fmpq(e.mu.raw(), q, 80);
        fmpq_clear(q);
        e.value = acb_si(1);
        bcs[1].push_back(std::move(e));
    }
    // c: empty

    auto asy = ode::calc_inf(de, bcs);
    ASSERT_EQ(asy.size(), 3u);

    // We expect each integral has 2 terms (2 fractional regions: 4/5 and 9/10).
    EXPECT_EQ(asy[0].size(), 2u) << "integral a should have 2 region terms";
    EXPECT_EQ(asy[1].size(), 2u) << "integral b should have 2 region terms";
    EXPECT_EQ(asy[2].size(), 2u) << "integral c should have 2 region terms";

    // For integral c (the all-zero-BC master), verify cross-coupling values.
    // c_a = 50/81 in Region 4/5 (order 1)
    // c_b = -10/9 in Region 9/10 (order 0)
    auto find_region_term = [&](const ode::AsyExpansion& asy_i, double frac_target) -> const ode::AsyTerm* {
        for (const auto& t : asy_i) {
            arb_t fpart;
            arb_init(fpart);
            // mu - floor(Re mu); for 0 < mu < 1, this is mu itself.
            arb_set(fpart, acb_realref(t.mu.raw()));
            fmpz_t z; fmpz_init(z);
            arf_get_fmpz(z, arb_midref(fpart), ARF_RND_FLOOR);
            arb_t zb; arb_init(zb); arb_set_fmpz(zb, z);
            arb_sub(fpart, fpart, zb, 200);
            double mid = arf_get_d(arb_midref(fpart), ARF_RND_NEAR);
            arb_clear(zb); fmpz_clear(z); arb_clear(fpart);
            if (std::abs(mid - frac_target) < 1e-10) return &t;
        }
        return nullptr;
    };

    const ode::AsyTerm* c_45 = find_region_term(asy[2], 0.8);
    const ode::AsyTerm* c_910 = find_region_term(asy[2], 0.9);
    ASSERT_NE(c_45, nullptr) << "no term in Region 4/5 for integral c";
    ASSERT_NE(c_910, nullptr) << "no term in Region 9/10 for integral c";

    // In Region 4/5: order-1 coeff of c = 50/81
    ASSERT_GE(c_45->exp.size(), 1u);
    ASSERT_GE(c_45->exp[0].size(), 2u);
    EXPECT_TRUE(acb_close_to_double(c_45->exp[0][0], 0.0, 0.0, 1e-12))
        << "c Region 4/5 order 0 should be 0; got " << c_45->exp[0][0].to_string(20);
    EXPECT_TRUE(acb_close_to_double(c_45->exp[0][1], 50.0 / 81.0, 0.0, 1e-12))
        << "c Region 4/5 order 1 should be 50/81 ≈ 0.617283...; got "
        << c_45->exp[0][1].to_string(20);

    // In Region 9/10: order-0 coeff of c = -10/9
    ASSERT_GE(c_910->exp.size(), 1u);
    ASSERT_GE(c_910->exp[0].size(), 1u);
    EXPECT_TRUE(acb_close_to_double(c_910->exp[0][0], -10.0 / 9.0, 0.0, 1e-12))
        << "c Region 9/10 order 0 should be -10/9 ≈ -1.1111...; got "
        << c_910->exp[0][0].to_string(20);
}
