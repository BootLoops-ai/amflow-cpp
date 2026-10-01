// SPDX-License-Identifier: MIT
// Tests for amflow::qft::topology.

#include <gtest/gtest.h>

#include "amflow/algebra/mpoly.hpp"
#include "amflow/qft/family_config.hpp"
#include "amflow/qft/family_uf.hpp"
#include "amflow/qft/topology.hpp"

#include <flint/fmpz_mpoly.h>

namespace alg = amflow::algebra;
namespace qft = amflow::qft;

// ---------------------------------------------------------------------------
//  factor_mpoly
// ---------------------------------------------------------------------------

TEST(TopologyTest, FactorMpoly_Trivial) {
    auto ctx = std::make_shared<alg::MpolyContext>(
        std::vector<std::string>{"x", "y"});
    alg::Mpoly x = alg::Mpoly::variable(ctx, 0);
    alg::Mpoly y = alg::Mpoly::variable(ctx, 1);

    alg::Mpoly p = (x + y) * (x - y);
    auto factors = qft::factor_mpoly(p);
    int nontrivial = 0;
    for (const auto& f : factors) {
        long len = fmpz_mpoly_length(f.raw(), ctx->raw());
        if (len > 1 || (len == 1 && (f - alg::Mpoly::one(ctx)).is_zero() == false
                        && (f + alg::Mpoly::one(ctx)).is_zero() == false)) {
            ++nontrivial;
        }
    }
    EXPECT_GE(nontrivial, 2);
}

TEST(TopologyTest, FactorMpoly_LinearTimesLinear) {
    auto ctx = std::make_shared<alg::MpolyContext>(
        std::vector<std::string>{"a", "b", "c"});
    alg::Mpoly a = alg::Mpoly::variable(ctx, 0);
    alg::Mpoly b = alg::Mpoly::variable(ctx, 1);
    alg::Mpoly c = alg::Mpoly::variable(ctx, 2);

    alg::Mpoly p = (a + b) * (a + c);
    auto factors = qft::factor_mpoly(p);
    EXPECT_GE(factors.size(), 2u);
}

TEST(TopologyTest, FactorMpoly_Irreducible) {
    auto ctx = std::make_shared<alg::MpolyContext>(
        std::vector<std::string>{"x", "y"});
    alg::Mpoly x = alg::Mpoly::variable(ctx, 0);
    alg::Mpoly y = alg::Mpoly::variable(ctx, 1);
    alg::Mpoly p = x * x + y * y + alg::Mpoly::one(ctx);
    auto factors = qft::factor_mpoly(p);
    EXPECT_EQ(factors.size(), 1u);
}

// ---------------------------------------------------------------------------
//  feynman_vars_in
// ---------------------------------------------------------------------------

TEST(TopologyTest, FeynmanVarsIn_OneLoopBubble) {
    auto fc = qft::FamilyConfig::build(
        "bubble",
        {"l"}, {"p"},
        {},
        {{"p^2", "s"}},
        {"l^2 - msq", "(l - p)^2 - msq"});
    auto uf = qft::evaluate_uf(fc, fc.propagators_after_conservation);

    auto u_vars = qft::feynman_vars_in(uf.u, uf.first_x_var, 2);
    EXPECT_EQ(u_vars.size(), 2u);
    EXPECT_EQ(u_vars[0], uf.first_x_var);
    EXPECT_EQ(u_vars[1], uf.first_x_var + 1);

    auto f_vars = qft::feynman_vars_in(uf.f, uf.first_x_var, 2);
    EXPECT_EQ(f_vars.size(), 2u);
}

// ---------------------------------------------------------------------------
//  ZeroSectorQ
// ---------------------------------------------------------------------------

TEST(TopologyTest, ZeroSectorQ_NonZero_OneLoopBubble) {
    auto fc = qft::FamilyConfig::build(
        "bubble",
        {"l"}, {"p"},
        {},
        {{"p^2", "s"}},
        {"l^2 - msq", "(l - p)^2 - msq"});
    EXPECT_FALSE(qft::zero_sector_q(fc, fc.propagators_after_conservation));
}

TEST(TopologyTest, ZeroSectorQ_NonZero_OneLoopMassiveTadpole) {
    auto fc = qft::FamilyConfig::build(
        "tadpole",
        {"l"}, {},
        {},
        {},
        {"l^2 - msq"});
    EXPECT_FALSE(qft::zero_sector_q(fc, fc.propagators_after_conservation));
}

TEST(TopologyTest, ZeroSectorQ_True_OneLoopMasslessTadpole) {
    auto fc = qft::FamilyConfig::build(
        "masslessTadpole",
        {"l"}, {},
        {},
        {},
        {"l^2"});
    EXPECT_TRUE(qft::zero_sector_q(fc, fc.propagators_after_conservation));
}

TEST(TopologyTest, ZeroSectorQ_True_TwoLoopMasslessSunsetSubsector) {
    auto fc = qft::FamilyConfig::build(
        "subsector",
        {"l1", "l2"}, {"p"},
        {},
        {{"p^2", "s"}},
        {"l1^2", "l2^2"});
    EXPECT_TRUE(qft::zero_sector_q(fc, fc.propagators_after_conservation));
}

TEST(TopologyTest, ZeroSectorQ_NonZero_TwoLoopSunrise) {
    auto fc = qft::FamilyConfig::build(
        "sunrise",
        {"l1", "l2"}, {"p"},
        {},
        {{"p^2", "s"}},
        {"l1^2 - m1sq", "l2^2 - m2sq", "(l1 + l2 - p)^2 - m3sq"});
    EXPECT_FALSE(qft::zero_sector_q(fc, fc.propagators_after_conservation));
}

// ---------------------------------------------------------------------------
//  AnalyzeTopology
// ---------------------------------------------------------------------------

TEST(TopologyTest, AnalyzeTopology_OneLoopBubble) {
    auto fc = qft::FamilyConfig::build(
        "bubble",
        {"l"}, {"p"},
        {},
        {{"p^2", "s"}},
        {"l^2 - msq", "(l - p)^2 - msq"});

    auto info = qft::analyze_topology(fc, fc.propagators_after_conservation);

    ASSERT_EQ(info.size(), 1u);
    EXPECT_EQ(info[0].loopnum, 1);
    EXPECT_EQ(info[0].var.size(), 2u);
}

TEST(TopologyTest, AnalyzeTopology_TwoLoopSunrise) {
    auto fc = qft::FamilyConfig::build(
        "sunrise",
        {"l1", "l2"}, {"p"},
        {},
        {{"p^2", "s"}},
        {"l1^2 - m1sq", "l2^2 - m2sq", "(l1 + l2 - p)^2 - m3sq"});

    auto info = qft::analyze_topology(fc, fc.propagators_after_conservation);

    ASSERT_EQ(info.size(), 1u);
    EXPECT_EQ(info[0].loopnum, 2);
    EXPECT_EQ(info[0].var.size(), 3u);

    ASSERT_EQ(info[0].mass.size(), 3u);
    for (const auto& m : info[0].mass) {
        EXPECT_FALSE(m.is_zero());
    }
}

TEST(TopologyTest, AnalyzeTopology_FactorisedU) {
    auto fc = qft::FamilyConfig::build(
        "doublebubble",
        {"l1", "l2"}, {"p"},
        {},
        {{"p^2", "s"}},
        {"l1^2 - m1sq", "(l1 - p)^2 - m1sq",
         "l2^2 - m2sq", "(l2 - p)^2 - m2sq"});

    auto uf = qft::evaluate_uf(fc, fc.propagators_after_conservation);
    ASSERT_FALSE(uf.degenerate);

    auto info = qft::analyze_topology(fc, fc.propagators_after_conservation);
    ASSERT_EQ(info.size(), 2u);
    EXPECT_EQ(info[0].loopnum, 1);
    EXPECT_EQ(info[1].loopnum, 1);
    EXPECT_EQ(info[0].var.size(), 2u);
    EXPECT_EQ(info[1].var.size(), 2u);
}
