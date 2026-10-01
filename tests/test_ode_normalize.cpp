// SPDX-License-Identifier: MIT
// Tests for amflow::ode::normalize (Layer 7).

#include <gtest/gtest.h>

#include "amflow/ode/blocks.hpp"
#include "amflow/ode/normalize.hpp"
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

class NormalizeTest : public ::testing::Test {
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

void set_poly_const_fraction(nm::FmpqPoly& poly, long num, long den) {
    fmpq_t q;
    fmpq_init(q);
    fmpq_set_si(q, num, den);
    poly.set_coeff_fmpq(0, q);
    fmpq_clear(q);
}

nm::RationalMatrix rat_mul(const nm::RationalMatrix& A,
                           const nm::RationalMatrix& B) {
    nm::RationalMatrix out(A.rows(), B.cols());
    for (std::size_t i = 0; i < A.rows(); ++i) {
        for (std::size_t j = 0; j < B.cols(); ++j) {
            nm::RationalFunction acc;
            for (std::size_t k = 0; k < A.cols(); ++k) {
                acc += A(i, k) * B(k, j);
            }
            out(i, j) = std::move(acc);
        }
    }
    return out;
}

nm::RationalMatrix rat_sub(const nm::RationalMatrix& A,
                           const nm::RationalMatrix& B) {
    nm::RationalMatrix out(A.rows(), A.cols());
    for (std::size_t i = 0; i < A.rows(); ++i) {
        for (std::size_t j = 0; j < A.cols(); ++j) {
            out(i, j) = A(i, j) - B(i, j);
        }
    }
    return out;
}

nm::RationalMatrix rat_derivative(const nm::RationalMatrix& A) {
    nm::RationalMatrix out(A.rows(), A.cols());
    for (std::size_t i = 0; i < A.rows(); ++i) {
        for (std::size_t j = 0; j < A.cols(); ++j) {
            out(i, j) = A(i, j).derivative();
        }
    }
    return out;
}

}  // namespace

// ============================================================================
//  AsymptoticBehavior
// ============================================================================

TEST_F(NormalizeTest, AsymptoticBehavior_ConstantDiagonal_GivesZeroBehavior) {
    nm::RationalMatrix de(1, 1);
    auto beh = ode::asymptotic_behavior(de);
    ASSERT_EQ(beh.size(), 1u);
    ASSERT_EQ(beh[0].size(), 1u);
    EXPECT_TRUE(acb_close_to_double(beh[0][0].mu, 0.0, 0.0));
    EXPECT_EQ(beh[0][0].log_power, 0L);
}

TEST_F(NormalizeTest, AsymptoticBehavior_OnePoleOneOverEta) {
    // mat = 2/eta -> L0 = 2.  Single behavior entry mu = 2, log_power = 0.
    nm::RationalMatrix de(1, 1);
    {
        nm::FmpqPoly num = make_poly({2});
        nm::FmpqPoly den = make_poly({0, 1});
        de(0, 0) = nm::RationalFunction(std::move(num), std::move(den));
    }
    auto beh = ode::asymptotic_behavior(de);
    ASSERT_EQ(beh.size(), 1u);
    ASSERT_EQ(beh[0].size(), 1u);
    EXPECT_TRUE(acb_close_to_double(beh[0][0].mu, 2.0, 0.0));
    EXPECT_EQ(beh[0][0].log_power, 0L);
}

// ============================================================================
//  NormalizeMat
// ============================================================================

TEST_F(NormalizeTest, NormalizeMat_AlreadyNormalised_IdentityTransform) {
    nm::RationalMatrix de(2, 2);
    auto nr = ode::normalize_mat(de);
    EXPECT_EQ(nr.T.rows(), 2u);
    EXPECT_EQ(nr.invT.rows(), 2u);
    EXPECT_EQ(nr.B.rows(), 2u);
    EXPECT_TRUE(nr.T(0, 0).is_one());
    EXPECT_TRUE(nr.T(1, 1).is_one());
    EXPECT_TRUE(nr.T(0, 1).is_zero());
    EXPECT_TRUE(nr.T(1, 0).is_zero());
}

TEST_F(NormalizeTest, NormalizeMat_RationalizesJordanBasisSimply) {
    nm::GlobalScope s;
    s.global.rationalize_pre = 100;
    s.global.working_pre = 140;
    s.global.silent_mode = true;
    s.commit();

    auto eta_pole = [](long num, long den) {
        nm::FmpqPoly n;
        set_poly_const_fraction(n, num, den);
        nm::FmpqPoly d = make_poly({0, 1});
        return nm::RationalFunction(std::move(n), std::move(d));
    };

    nm::RationalMatrix de(2, 2);
    de(0, 0) = eta_pole(3, 2);
    de(0, 1) = eta_pole(-1, 2);
    de(1, 0) = eta_pole(3, 2);
    de(1, 1) = eta_pole(-1, 2);

    auto nr = ode::normalize_mat(de);
    EXPECT_LE(ode::poincare_rank(nr.B), 0);

    nm::RationalMatrix lhs = rat_sub(rat_mul(rat_mul(nr.invT, de), nr.T),
                                     rat_mul(nr.invT, rat_derivative(nr.T)));
    EXPECT_EQ(lhs, nr.B);
    EXPECT_EQ(rat_mul(nr.invT, nr.T), nm::RationalMatrix::identity(2));
}

TEST_F(NormalizeTest, NormalizeMat_ExactIntegerEigenFloor) {
    nm::GlobalScope s;
    s.global.working_pre = 30;
    s.global.chop_pre = 8;
    s.global.rationalize_pre = 20;
    s.global.silent_mode = true;
    s.commit();

    auto eta_pole = [](long num, long den) {
        nm::FmpqPoly n;
        set_poly_const_fraction(n, num, den);
        nm::FmpqPoly d = make_poly({0, 1});
        return nm::RationalFunction(std::move(n), std::move(d));
    };

    nm::RationalMatrix de(1, 1);
    de(0, 0) = eta_pole(1, 1);

    auto nr = ode::normalize_mat(de);
    nm::RationalFunction lead = nr.B(0, 0);
    lead.multiply_by_eta_power(1);
    nm::AcbValue zero;
    nm::AcbValue lead0 = lead.evaluate(zero.raw(), nm::working_prec_bits());
    EXPECT_TRUE(nm::acb_is_chop_zero(lead0.raw(), nm::chop_pre()));
}

TEST_F(NormalizeTest, NormalizeMat_NearIntegerRationalUsesExactFloor) {
    nm::GlobalScope s;
    s.global.working_pre = 30;
    s.global.chop_pre = 8;
    s.global.rationalize_pre = 20;
    s.global.silent_mode = true;
    s.commit();

    auto eta_pole = [](long num, long den) {
        nm::FmpqPoly n;
        set_poly_const_fraction(n, num, den);
        nm::FmpqPoly d = make_poly({0, 1});
        return nm::RationalFunction(std::move(n), std::move(d));
    };

    nm::RationalMatrix de(1, 1);
    de(0, 0) = eta_pole(999999999999L, 1000000000000L);

    auto nr = ode::normalize_mat(de);
    nm::RationalFunction lead = nr.B(0, 0);
    lead.multiply_by_eta_power(1);
    nm::AcbValue zero;
    nm::AcbValue lead0 = lead.evaluate(zero.raw(), nm::working_prec_bits());
    EXPECT_FALSE(nm::acb_is_chop_zero(lead0.raw(), nm::chop_pre()));
}

TEST_F(NormalizeTest, NormalizeMat_AlgebraicProjectorShearingReachesJordanGap) {
    nm::GlobalScope s;
    s.global.working_pre = 120;
    s.global.chop_pre = 20;
    s.global.silent_mode = true;
    s.commit();

    auto eta_pole = [](long num, long den) {
        nm::FmpqPoly n;
        set_poly_const_fraction(n, num, den);
        nm::FmpqPoly d = make_poly({0, 1});
        return nm::RationalFunction(std::move(n), std::move(d));
    };

    nm::RationalMatrix de(3, 3);
    de(0, 0) = eta_pole(5, 2);
    de(0, 1) = eta_pole(-3, 2);
    de(0, 2) = eta_pole(-3, 2);
    de(1, 0) = eta_pole(7, 2);
    de(1, 1) = eta_pole(-3, 2);
    de(1, 2) = eta_pole(-7, 2);
    de(2, 0) = eta_pole(2, 1);
    de(2, 1) = eta_pole(-1, 1);
    de(2, 2) = eta_pole(-1, 1);

    EXPECT_THROW((void)ode::normalize_mat(de), std::runtime_error);
}

TEST_F(NormalizeTest, NormalizeMat_AlreadyJordanBlockAvoidsNumericBasis) {
    nm::GlobalScope s;
    s.global.silent_mode = true;
    s.commit();

    auto eta_pole = [](long num, long den) {
        nm::FmpqPoly n;
        set_poly_const_fraction(n, num, den);
        nm::FmpqPoly d = make_poly({0, 1});
        return nm::RationalFunction(std::move(n), std::move(d));
    };

    nm::RationalMatrix de(3, 3);
    de(0, 1) = eta_pole(1, 1);
    de(1, 2) = eta_pole(1, 1);

    auto nr = ode::normalize_mat(de);
    EXPECT_LE(ode::poincare_rank(nr.B), 0);
    EXPECT_TRUE(nr.T(0, 0).is_one());
    EXPECT_TRUE(nr.T(1, 1).is_one());
    EXPECT_TRUE(nr.T(2, 2).is_one());
    EXPECT_TRUE(nr.B(0, 1) == de(0, 1));
    EXPECT_TRUE(nr.B(1, 2) == de(1, 2));
}

TEST_F(NormalizeTest, NormalizeMat_OffDiagonalFuchsianCleanup) {
    nm::GlobalScope s;
    s.global.silent_mode = true;
    s.commit();

    nm::RationalMatrix de(2, 2);
    de(1, 0) = nm::RationalFunction::monomial(-2);

    auto nr = ode::normalize_mat(de);
    EXPECT_LE(ode::poincare_rank(nr.B), 0);
    EXPECT_TRUE(nr.B(1, 0).is_zero());
}

TEST_F(NormalizeTest, NormalizeMat_DiagonalToFuchsian_ReducesNilpotentDoublePole) {
    nm::GlobalScope s;
    s.global.silent_mode = true;
    s.commit();

    nm::RationalMatrix de(2, 2);
    de(0, 1) = nm::RationalFunction::monomial(-2);

    auto nr = ode::normalize_mat(de);
    EXPECT_LE(ode::poincare_rank(nr.B), 0);
}

TEST_F(NormalizeTest, NormalizeMat_DiagonalToFuchsian_RejectsIrreducibleDoublePole) {
    nm::GlobalScope s;
    s.global.silent_mode = true;
    s.commit();

    nm::RationalMatrix de(1, 1);
    de(0, 0) = nm::RationalFunction::monomial(-2);

    EXPECT_THROW((void)ode::normalize_mat(de), std::runtime_error);
}

// ============================================================================
//  normalize_mat_for_calc_zero
// ============================================================================

TEST_F(NormalizeTest, NormalizeMatForCalcZero_AlreadyJordanProducesIdentityRotations) {
    nm::GlobalScope s;
    s.global.silent_mode = true;
    s.commit();

    nm::RationalMatrix de(2, 2);
    de(0, 1) = nm::RationalFunction::monomial(-1);

    auto nr = ode::normalize_mat_for_calc_zero(de);
    EXPECT_LE(ode::poincare_rank(nr.base.B), 0);
    EXPECT_GE(nr.rotations.size(), 1u);
}
