// SPDX-License-Identifier: MIT
// Tests for amflow::qft::family_config.

#include <gtest/gtest.h>

#include "amflow/algebra/mpoly.hpp"
#include "amflow/numeric/options.hpp"
#include "amflow/qft/family_config.hpp"

#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace alg = amflow::algebra;
namespace nm  = amflow::numeric;
namespace qft = amflow::qft;

namespace {

class FamilyConfigTest : public ::testing::Test {
protected:
    void SetUp() override    { nm::set_default_options(); }
    void TearDown() override { nm::set_default_options(); }
};

qft::FamilyConfig make_box1() {
    return qft::FamilyConfig::build(
        "box1",
        {"l"},
        {"p1", "p2", "p3", "p4"},
        {{"p4", "-p1 - p2 - p3"}},
        {{"p1^2", "0"},
         {"p2^2", "0"},
         {"p3^2", "0"},
         {"p4^2", "0"},
         {"(p1 + p2)^2", "s"},
         {"(p1 + p3)^2", "t"}},
        {"l^2", "(l + p1)^2", "(l + p1 + p2)^2",
         "(l + p1 + p2 + p4)^2"});
}

}  // namespace

// =========================================================================
//  Construction basics
// =========================================================================

TEST_F(FamilyConfigTest, Build_Box1_BasicSizes) {
    auto fc = make_box1();
    EXPECT_EQ(fc.family, "box1");
    EXPECT_EQ(fc.n_loops(), 1u);
    EXPECT_EQ(fc.n_red_legs(), 3u);
    EXPECT_EQ(fc.n_invariants(), 2u);
    EXPECT_EQ(fc.n_propagators(), 4u);
    EXPECT_EQ(fc.n_sp(), 4u);
}

TEST_F(FamilyConfigTest, Build_FractionalPropagatorClearsCommonDenominator) {
    auto fc = qft::FamilyConfig::build(
        "fracprop",
        {"l"}, {"p"},
        {},
        {{"p^2", "0"}},
        {"1/2 + l*p"});

    ASSERT_EQ(fc.n_propagators(), 1u);
    ASSERT_EQ(fc.propagator_denominators, std::vector<std::string>({"2"}));
    EXPECT_EQ(fc.propagators_after_conservation[0],
              alg::Mpoly::from_string(fc.ctx, "1 + 2*l*p"));
}

TEST_F(FamilyConfigTest, Build_Box1_VariableOrdering) {
    auto fc = make_box1();
    EXPECT_EQ(fc.ctx->var_index("l"),  0);
    EXPECT_EQ(fc.ctx->var_index("p1"), 1);
    EXPECT_EQ(fc.ctx->var_index("p2"), 2);
    EXPECT_EQ(fc.ctx->var_index("p3"), 3);
    EXPECT_EQ(fc.ctx->var_index("s"),  4);
    EXPECT_EQ(fc.ctx->var_index("t"),  5);
    EXPECT_EQ(fc.ctx->var_index("p4"), -1);
}

TEST_F(FamilyConfigTest, Build_Box1_ReducedLegs) {
    auto fc = make_box1();
    EXPECT_EQ(fc.reduced_legs, (std::vector<std::string>{"p1", "p2", "p3"}));
}

TEST_F(FamilyConfigTest, Build_RejectsUnknownConservationLhs) {
    EXPECT_THROW(qft::FamilyConfig::build("f",
                                           {"l"}, {"p1", "p2"},
                                           {{"q1", "-p2"}},
                                           {}, {"l^2"}),
                 std::invalid_argument);
}

TEST_F(FamilyConfigTest, Build_DuplicateLoopName) {
    EXPECT_THROW(qft::FamilyConfig::build("f",
                                           {"l", "l"}, {"p1"},
                                           {}, {{"p1*p1", "0"}}, {"l^2"}),
                 std::invalid_argument);
}

// =========================================================================
//  Replacement algebra
// =========================================================================

TEST_F(FamilyConfigTest, Replacement_BasicLookup) {
    auto fc = make_box1();
    EXPECT_EQ(fc.reduced_replacement.count("p1*p2"), 1u);
    EXPECT_EQ(fc.reduced_replacement.count("p1*p1"), 1u);
    EXPECT_EQ(fc.reduced_replacement.count("p2*p3"), 1u);
    auto p1p2 = fc.reduced_replacement.at("p1*p2").clone();
    alg::Mfrac expected(alg::Mpoly::variable(fc.ctx, "s"),
                       alg::Mpoly::constant(fc.ctx, 2));
    EXPECT_EQ(p1p2, expected);
}

TEST_F(FamilyConfigTest, Replacement_InfersCrossTermsFromCompactQuadratics) {
    auto fc = make_box1();

    alg::Mfrac p1p2_expected(alg::Mpoly::variable(fc.ctx, "s"),
                            alg::Mpoly::constant(fc.ctx, 2));
    EXPECT_EQ(fc.reduced_replacement.at("p1*p2"), p1p2_expected);

    alg::Mfrac p1p3_expected(alg::Mpoly::variable(fc.ctx, "t"),
                            alg::Mpoly::constant(fc.ctx, 2));
    EXPECT_EQ(fc.reduced_replacement.at("p1*p3"), p1p3_expected);

    auto neg_st = alg::Mpoly::from_string(fc.ctx, "-s - t");
    alg::Mfrac p2p3_expected(std::move(neg_st),
                            alg::Mpoly::constant(fc.ctx, 2));
    EXPECT_EQ(fc.reduced_replacement.at("p2*p3"), p2p3_expected);
}

TEST_F(FamilyConfigTest, ApplyReplacement_OnSingleMonomial) {
    auto fc = make_box1();
    auto p1p2 = alg::Mpoly::from_string(fc.ctx, "p1*p2");
    auto r = fc.apply_replacement(p1p2);
    alg::Mfrac expected(alg::Mpoly::variable(fc.ctx, "s"),
                       alg::Mpoly::constant(fc.ctx, 2));
    EXPECT_EQ(r, expected);
}

TEST_F(FamilyConfigTest, ApplyReplacement_ExpandsAndReduces) {
    auto fc = make_box1();
    auto expr = alg::Mpoly::from_string(fc.ctx, "(p1 + p2)^2");
    auto r = fc.apply_replacement(expr);
    alg::Mfrac expected = alg::Mfrac::from_mpoly(alg::Mpoly::variable(fc.ctx, "s"));
    EXPECT_EQ(r, expected);
}

TEST_F(FamilyConfigTest, ApplyReplacement_LeavesLoopAlone) {
    auto fc = make_box1();
    auto expr = alg::Mpoly::from_string(fc.ctx, "l^2 + 3*l*p1");
    auto r = fc.apply_replacement(expr);
    alg::Mfrac expected = alg::Mfrac::from_mpoly(expr.clone());
    EXPECT_EQ(r, expected);
}

TEST_F(FamilyConfigTest, ApplyReplacement_OnPropagator_BoxOneLineThree) {
    auto fc = make_box1();
    auto third   = fc.propagators_after_conservation[2].clone();
    auto reduced = fc.apply_replacement(third);
    auto expected = alg::Mfrac::from_mpoly(
        alg::Mpoly::from_string(fc.ctx, "l^2 + 2*l*p1 + 2*l*p2 + s"));
    EXPECT_EQ(reduced, expected);
}

TEST_F(FamilyConfigTest, SPList_LoopLeg_Order) {
    auto fc = make_box1();
    ASSERT_EQ(fc.sp_list.size(), 4u);
    EXPECT_EQ(fc.sp_list[0], alg::Mpoly::from_string(fc.ctx, "l^2"));
    EXPECT_EQ(fc.sp_list[1], alg::Mpoly::from_string(fc.ctx, "l*p1"));
    EXPECT_EQ(fc.sp_list[2], alg::Mpoly::from_string(fc.ctx, "l*p2"));
    EXPECT_EQ(fc.sp_list[3], alg::Mpoly::from_string(fc.ctx, "l*p3"));
}

// =========================================================================
//  Conservation
// =========================================================================

TEST_F(FamilyConfigTest, Conservation_KeyP4_NotInContext) {
    auto fc = make_box1();
    EXPECT_EQ(fc.ctx->var_index("p4"), -1);
}

TEST_F(FamilyConfigTest, Conservation_FourthPropagator_IsLSquared) {
    auto fc = make_box1();
    auto& fourth = fc.propagators_after_conservation[3];
    auto expected = alg::Mpoly::from_string(fc.ctx, "l^2 - 2*l*p3 + p3^2");
    EXPECT_EQ(fourth, expected);
}

TEST_F(FamilyConfigTest, Conservation_StoredP4Rule) {
    auto fc = make_box1();
    ASSERT_EQ(fc.conservation.count("p4"), 1u);
    auto expected = alg::Mpoly::from_string(fc.ctx, "-p1 - p2 - p3");
    EXPECT_EQ(fc.conservation.at("p4"), expected);
}
