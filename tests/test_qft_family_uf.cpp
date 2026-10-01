// SPDX-License-Identifier: MIT
// Tests for amflow::qft::family_uf.

#include <gtest/gtest.h>

#include "amflow/algebra/mpoly.hpp"
#include "amflow/qft/family_config.hpp"
#include "amflow/qft/family_uf.hpp"

#include <cmath>

namespace alg = amflow::algebra;
namespace qft = amflow::qft;

namespace {

alg::Mfrac sub_var(const alg::Mfrac& f, long var, long value) {
    return f.substitute(var, value);
}

}  // namespace

// ---------------------------------------------------------------------------
//  One-loop bubble:  D1 = l^2 - msq, D2 = (l - p)^2 - msq, p^2 -> s
//  Expected: u = x1+x2,  f = -x1*x2*s,  f0 = msq*(x1+x2).
// ---------------------------------------------------------------------------

TEST(FamilyUFTest, OneLoopBubble_UFCorrect) {
    auto fc = qft::FamilyConfig::build(
        "bubble",
        {"l"}, {"p"},
        {},
        {{"p^2", "s"}},
        {"l^2 - msq", "(l - p)^2 - msq"});

    auto uf = qft::evaluate_uf(fc, fc.propagators_after_conservation);
    ASSERT_FALSE(uf.degenerate);
    ASSERT_TRUE(uf.u.is_polynomial());

    {
        alg::Mfrac u_eval = sub_var(uf.u, uf.first_x_var, 1);
        u_eval = sub_var(u_eval, uf.first_x_var + 1, 1);
        EXPECT_EQ(u_eval.to_string(), "2");
    }
    EXPECT_EQ(sub_var(sub_var(uf.u, uf.first_x_var, 1),
                      uf.first_x_var + 1, 0).to_string(), "1");
    EXPECT_EQ(sub_var(sub_var(uf.u, uf.first_x_var, 0),
                      uf.first_x_var + 1, 1).to_string(), "1");
    EXPECT_EQ(sub_var(sub_var(uf.u, uf.first_x_var, 0),
                      uf.first_x_var + 1, 0).to_string(), "0");

    long inv_s   = fc.var_index("s");
    long inv_msq = fc.var_index("msq");
    {
        alg::Mfrac fe = uf.f.clone();
        fe = sub_var(fe, uf.first_x_var, 1);
        fe = sub_var(fe, uf.first_x_var + 1, 1);
        fe = sub_var(fe, inv_s, 5);
        fe = sub_var(fe, inv_msq, 7);
        EXPECT_EQ(fe.to_string(), "-5");
    }
    {
        alg::Mfrac fe = uf.f.clone();
        fe = sub_var(fe, uf.first_x_var, 2);
        fe = sub_var(fe, uf.first_x_var + 1, 3);
        fe = sub_var(fe, inv_s, 5);
        fe = sub_var(fe, inv_msq, 7);
        EXPECT_EQ(fe.to_string(), "-30");
    }
    {
        alg::Mfrac fe = uf.f.clone();
        fe = sub_var(fe, uf.first_x_var, 1);
        fe = sub_var(fe, uf.first_x_var + 1, 2);
        fe = sub_var(fe, inv_s, 3);
        fe = sub_var(fe, inv_msq, 4);
        EXPECT_EQ(fe.to_string(), "-6");
    }
}

// ---------------------------------------------------------------------------
//  Two-loop sunrise: U = x1 x2 + x1 x3 + x2 x3, F = -s x1 x2 x3
//  (mass terms separated into f0 = m1*x1 + m2*x2 + m3*x3).
// ---------------------------------------------------------------------------

TEST(FamilyUFTest, TwoLoopSunrise_UFCorrect) {
    auto fc = qft::FamilyConfig::build(
        "sunrise",
        {"l1", "l2"}, {"p"},
        {},
        {{"p^2", "s"}},
        {"l1^2 - m1sq", "l2^2 - m2sq", "(l1 + l2 - p)^2 - m3sq"});

    auto uf = qft::evaluate_uf(fc, fc.propagators_after_conservation);
    ASSERT_FALSE(uf.degenerate);

    long x1 = uf.first_x_var;
    long x2 = uf.first_x_var + 1;
    long x3 = uf.first_x_var + 2;

    long inv_s    = fc.var_index("s");
    long inv_m1sq = fc.var_index("m1sq");
    long inv_m2sq = fc.var_index("m2sq");
    long inv_m3sq = fc.var_index("m3sq");

    auto eval = [&](const alg::Mfrac& f, long X1, long X2, long X3,
                    long S, long M1, long M2, long M3) -> std::string {
        alg::Mfrac fe = f.clone();
        fe = sub_var(fe, x1, X1);
        fe = sub_var(fe, x2, X2);
        fe = sub_var(fe, x3, X3);
        fe = sub_var(fe, inv_s, S);
        fe = sub_var(fe, inv_m1sq, M1);
        fe = sub_var(fe, inv_m2sq, M2);
        fe = sub_var(fe, inv_m3sq, M3);
        return fe.to_string();
    };

    EXPECT_EQ(eval(uf.u, 1, 1, 1, 0, 0, 0, 0), "3");
    EXPECT_EQ(eval(uf.u, 2, 3, 5, 0, 0, 0, 0), "31");
    EXPECT_EQ(eval(uf.u, 1, 0, 0, 0, 0, 0, 0), "0");
    EXPECT_EQ(eval(uf.u, 1, 1, 0, 0, 0, 0, 0), "1");
    EXPECT_EQ(eval(uf.u, 0, 1, 1, 0, 0, 0, 0), "1");

    EXPECT_EQ(eval(uf.f, 1, 1, 1, 0, 1, 1, 1), "0");
    EXPECT_EQ(eval(uf.f, 1, 1, 1, 4, 0, 0, 0), "-4");
    EXPECT_EQ(eval(uf.f, 1, 1, 1, 5, 1, 2, 3), "-5");
    EXPECT_EQ(eval(uf.f, 2, 3, 5, 1, 0, 0, 0), "-30");
    EXPECT_EQ(eval(uf.f, 2, 3, 5, 0, 1, 0, 0), "0");
    EXPECT_EQ(eval(uf.f, 2, 3, 5, 2, 0, 1, 0), "-60");

    EXPECT_EQ(eval(uf.f0, 1, 1, 1, 0, 1, 1, 1), "3");
    EXPECT_EQ(eval(uf.f0, 1, 2, 3, 0, 4, 5, 6), "32");

    {
        alg::Mfrac F_phys = uf.f.clone();
        alg::Mfrac u_f0 = uf.u.clone();
        u_f0 *= uf.f0;
        F_phys += u_f0;
        alg::Mfrac fe = F_phys.clone();
        fe = sub_var(fe, x1, 1);
        fe = sub_var(fe, x2, 1);
        fe = sub_var(fe, x3, 1);
        fe = sub_var(fe, inv_s,    0);
        fe = sub_var(fe, inv_m1sq, 1);
        fe = sub_var(fe, inv_m2sq, 1);
        fe = sub_var(fe, inv_m3sq, 1);
        EXPECT_EQ(fe.to_string(), "9");
    }
}

// ---------------------------------------------------------------------------
//  Massless one-loop tadpole.
// ---------------------------------------------------------------------------

TEST(FamilyUFTest, OneLoopMasslessTadpole) {
    auto fc = qft::FamilyConfig::build(
        "tadpole",
        {"l"}, {},
        {},
        {},
        {"l^2"});

    auto uf = qft::evaluate_uf(fc, fc.propagators_after_conservation);
    ASSERT_FALSE(uf.degenerate);

    long x1 = uf.first_x_var;
    EXPECT_EQ(uf.u.to_string(), "__feyn_x1");
    EXPECT_TRUE(uf.f.is_zero());

    alg::Mfrac u_eval = sub_var(uf.u, x1, 7);
    EXPECT_EQ(u_eval.to_string(), "7");
}

// ---------------------------------------------------------------------------
//  Degenerate Ax (two-loop family but only l1 in propagators).
// ---------------------------------------------------------------------------

TEST(FamilyUFTest, DegenerateAxIsReported) {
    auto fc = qft::FamilyConfig::build(
        "deg",
        {"l1", "l2"}, {},
        {},
        {},
        {"l1^2", "l1^2 - 1"});

    auto uf = qft::evaluate_uf(fc, fc.propagators_after_conservation);
    EXPECT_TRUE(uf.degenerate);
}

// ---------------------------------------------------------------------------
//  Conservation-applied 4-point box.
// ---------------------------------------------------------------------------

TEST(FamilyUFTest, ConservationApplied) {
    auto fc = qft::FamilyConfig::build(
        "box",
        {"l"}, {"p1", "p2", "p3", "p4"},
        {{"p4", "-p1 - p2 - p3"}},
        {{"p1*p1", "0"},
         {"p2*p2", "0"},
         {"p3*p3", "0"},
         {"p1*p2", "s/2"},
         {"p2*p3", "t/2"},
         {"p1*p3", "(-s-t)/2"}},
        {"l^2", "(l + p1)^2", "(l + p1 + p2)^2", "(l - p4)^2"});

    auto uf = qft::evaluate_uf(fc, fc.propagators_after_conservation);
    ASSERT_FALSE(uf.degenerate);

    long x1 = uf.first_x_var;
    long x2 = uf.first_x_var + 1;
    long x3 = uf.first_x_var + 2;
    long x4 = uf.first_x_var + 3;

    {
        alg::Mfrac u_eval = uf.u.clone();
        u_eval = sub_var(u_eval, x1, 1);
        u_eval = sub_var(u_eval, x2, 1);
        u_eval = sub_var(u_eval, x3, 1);
        u_eval = sub_var(u_eval, x4, 1);
        EXPECT_EQ(u_eval.to_string(), "4");
    }
    {
        alg::Mfrac u_eval = uf.u.clone();
        u_eval = sub_var(u_eval, x1, 2);
        u_eval = sub_var(u_eval, x2, 3);
        u_eval = sub_var(u_eval, x3, 5);
        u_eval = sub_var(u_eval, x4, 7);
        EXPECT_EQ(u_eval.to_string(), "17");
    }

    long inv_s = fc.var_index("s");
    long inv_t = fc.var_index("t");
    {
        alg::Mfrac fe = uf.f.clone();
        fe = sub_var(fe, x1, 1);
        fe = sub_var(fe, x2, 1);
        fe = sub_var(fe, x3, 1);
        fe = sub_var(fe, x4, 1);
        fe = sub_var(fe, inv_s, 2);
        fe = sub_var(fe, inv_t, 3);
        EXPECT_EQ(fe.to_string(), "-5");
    }
    {
        alg::Mfrac fe = uf.f.clone();
        fe = sub_var(fe, x1, 2);
        fe = sub_var(fe, x2, 3);
        fe = sub_var(fe, x3, 5);
        fe = sub_var(fe, x4, 7);
        fe = sub_var(fe, inv_s, 1);
        fe = sub_var(fe, inv_t, 0);
        EXPECT_EQ(fe.to_string(), "-10");
    }
    {
        alg::Mfrac fe = uf.f.clone();
        fe = sub_var(fe, x1, 1);
        fe = sub_var(fe, x2, 1);
        fe = sub_var(fe, x3, 1);
        fe = sub_var(fe, x4, 1);
        fe = sub_var(fe, inv_s, 0);
        fe = sub_var(fe, inv_t, 4);
        EXPECT_EQ(fe.to_string(), "-4");
    }
}

TEST(FamilyUFTest, ABCResultShape) {
    auto fc = qft::FamilyConfig::build(
        "shape",
        {"l"}, {"p"},
        {},
        {{"p^2", "s"}},
        {"l^2", "(l-p)^2"});

    auto abc = qft::evaluate_abc(fc, fc.propagators_after_conservation);
    EXPECT_EQ(abc.Ax.rows(), 1u);
    EXPECT_EQ(abc.Ax.cols(), 1u);
    EXPECT_EQ(abc.Bx.size(), 1u);
    EXPECT_EQ(abc.first_x_var, fc.ctx->n_vars());
    EXPECT_EQ(abc.uf_ctx->n_vars(), fc.ctx->n_vars() + 2);
}
