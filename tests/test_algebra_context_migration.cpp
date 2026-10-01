// SPDX-License-Identifier: MIT
// Tests for amflow::algebra::mfrac_to_ctx / project_mfrac_dropping
// (primitives that consolidate context-migration code).

#include <gtest/gtest.h>

#include "amflow/algebra/context_migration.hpp"
#include "amflow/algebra/mpoly.hpp"

#include <memory>
#include <stdexcept>

namespace alg = amflow::algebra;

namespace {

std::shared_ptr<alg::MpolyContext> make_ctx(std::vector<std::string> names) {
    return std::make_shared<alg::MpolyContext>(std::move(names));
}

}  // namespace

// ---------------------------------------------------------------------------
//  Strict lift (mpoly_to_ctx / mfrac_to_ctx)
// ---------------------------------------------------------------------------

TEST(MpolyToCtx, IdenticalCtxIsClone) {
    auto ctx = make_ctx({"x"});
    auto p = alg::Mpoly::from_string(ctx, "x^2 + 1");
    auto q = alg::mpoly_to_ctx(p, ctx);
    EXPECT_EQ(p, q);
    // Must be a clone, not the same object.
    EXPECT_NE(p.raw(), q.raw());
}

TEST(MpolyToCtx, NameLookupSurvivesCtxIdentityMismatch) {
    auto ctx_a = make_ctx({"x", "y"});
    auto ctx_b = make_ctx({"x", "y"});      // same names, different shared_ptr

    auto p_a = alg::Mpoly::from_string(ctx_a, "x^2 + y");
    auto p_b = alg::mpoly_to_ctx(p_a, ctx_b);

    auto expected_b = alg::Mpoly::from_string(ctx_b, "x^2 + y");
    EXPECT_EQ(p_b, expected_b);
}

TEST(MpolyToCtx, ReorderedVarsAreRespected) {
    auto ctx_a = make_ctx({"x", "y"});
    auto ctx_b = make_ctx({"y", "x"});      // reversed order

    auto p_a = alg::Mpoly::from_string(ctx_a, "x^2*y + 3*y");
    auto p_b = alg::mpoly_to_ctx(p_a, ctx_b);

    // In ctx_b, "x" is at index 1, "y" at 0.  The polynomial should still
    // evaluate the same.
    auto expected_b = alg::Mpoly::from_string(ctx_b, "x^2*y + 3*y");
    EXPECT_EQ(p_b, expected_b);
}

TEST(MpolyToCtx, MissingVariableThrows) {
    auto ctx_a = make_ctx({"x", "y"});
    auto ctx_b = make_ctx({"x"});           // no "y"

    auto p_a = alg::Mpoly::from_string(ctx_a, "x + y");
    EXPECT_THROW(alg::mpoly_to_ctx(p_a, ctx_b), std::invalid_argument);
}

TEST(MpolyToCtx, UnusedSourceVariableIsAllowed) {
    auto ctx_a = make_ctx({"x", "y"});
    auto ctx_b = make_ctx({"x"});

    // p has only x dependence (no y), so dst has all vars present in src's
    // monomials.  But our strict lift checks ALL src vars; so this should
    // throw because src's "y" is missing in dst.  This is the strict
    // semantics — match upstream MMA behaviour.
    auto p_a = alg::Mpoly::from_string(ctx_a, "x^2 + 1");
    EXPECT_THROW(alg::mpoly_to_ctx(p_a, ctx_b), std::invalid_argument);
}

TEST(MfracToCtx, ReorderedVarsRoundTrip) {
    auto ctx_a = make_ctx({"x", "y"});
    auto ctx_b = make_ctx({"y", "x"});

    auto r_a = alg::Mfrac(alg::Mpoly::from_string(ctx_a, "x + y"),
                            alg::Mpoly::from_string(ctx_a, "x - y"));
    auto r_b = alg::mfrac_to_ctx(r_a, ctx_b);
    auto expected = alg::Mfrac(alg::Mpoly::from_string(ctx_b, "x + y"),
                                 alg::Mpoly::from_string(ctx_b, "x - y"));
    EXPECT_EQ(r_b, expected);
}

// ---------------------------------------------------------------------------
//  Lossy projection (project_*_dropping)
// ---------------------------------------------------------------------------

TEST(ProjectDropping, SentinelPrefixIsDroppedFromMpoly) {
    auto ctx_a = make_ctx({"x", "__amf_eta", "y"});
    auto ctx_b = make_ctx({"x", "y"});

    // p = x^2 + __amf_eta * y + y^2 + 7
    auto p_a = alg::Mpoly::from_string(ctx_a,
        "x^2 + __amf_eta * y + y^2 + 7");
    auto p_b = alg::project_mpoly_dropping(p_a, ctx_b);

    // The __amf_eta * y monomial should drop.
    auto expected = alg::Mpoly::from_string(ctx_b, "x^2 + y^2 + 7");
    EXPECT_EQ(p_b, expected);
}

TEST(ProjectDropping, SentinelPrefixIsDroppedFromMfrac) {
    auto ctx_a = make_ctx({"x", "__feyn_0"});
    auto ctx_b = make_ctx({"x"});

    auto r_a = alg::Mfrac(
        alg::Mpoly::from_string(ctx_a, "x + __feyn_0"),
        alg::Mpoly::from_string(ctx_a, "1"));
    auto r_b = alg::project_mfrac_dropping(r_a, ctx_b);

    auto expected = alg::Mfrac::from_mpoly(
        alg::Mpoly::from_string(ctx_b, "x"));
    EXPECT_EQ(r_b, expected);
}

TEST(ProjectDropping, MultipleSentinelPrefixesAllDrop) {
    auto ctx_a = make_ctx({"x", "__amf_a", "__feyn_b", "__zsq_c"});
    auto ctx_b = make_ctx({"x"});

    auto p_a = alg::Mpoly::from_string(ctx_a,
        "x^2 + __amf_a + __feyn_b + __zsq_c + 5");
    auto p_b = alg::project_mpoly_dropping(p_a, ctx_b);

    auto expected = alg::Mpoly::from_string(ctx_b, "x^2 + 5");
    EXPECT_EQ(p_b, expected);
}

TEST(ProjectDropping, NonSentinelMissingVarStillThrows) {
    auto ctx_a = make_ctx({"x", "regular_other"});
    auto ctx_b = make_ctx({"x"});

    auto p_a = alg::Mpoly::from_string(ctx_a, "x + regular_other");
    EXPECT_THROW(alg::project_mpoly_dropping(p_a, ctx_b),
                 std::invalid_argument);
}

TEST(ProjectDropping, CustomDropPrefixesWork) {
    auto ctx_a = make_ctx({"x", "tmp_var"});
    auto ctx_b = make_ctx({"x"});

    auto p_a = alg::Mpoly::from_string(ctx_a, "x + tmp_var");
    auto p_b = alg::project_mpoly_dropping(p_a, ctx_b, {"tmp_"});
    auto expected = alg::Mpoly::from_string(ctx_b, "x");
    EXPECT_EQ(p_b, expected);
}

TEST(ProjectDropping, EmptyDropListBehavesLikeStrict) {
    auto ctx_a = make_ctx({"x", "__amf_eta"});
    auto ctx_b = make_ctx({"x"});
    auto p_a = alg::Mpoly::from_string(ctx_a, "x");
    EXPECT_THROW(alg::project_mpoly_dropping(p_a, ctx_b, {}),
                 std::invalid_argument);
}

TEST(ProjectDropping, IdentityCtxIsClone) {
    auto ctx = make_ctx({"x", "__amf_eta"});
    auto p = alg::Mpoly::from_string(ctx, "x + __amf_eta");
    auto q = alg::project_mpoly_dropping(p, ctx);
    // Same ctx → exact clone (no projection happens).
    EXPECT_EQ(p, q);
}

TEST(ProjectDropping, MfracDenominatorVanishingThrows) {
    auto ctx_a = make_ctx({"x", "__amf_eta"});
    auto ctx_b = make_ctx({"x"});

    // Numerator = x; denominator = __amf_eta (drops to zero on projection).
    auto r_a = alg::Mfrac(
        alg::Mpoly::variable(ctx_a, 0),
        alg::Mpoly::variable(ctx_a, 1));
    EXPECT_THROW(alg::project_mfrac_dropping(r_a, ctx_b), std::runtime_error);
}

TEST(ProjectDropping, MonomialMixingNormalAndSentinelDropsEntirely) {
    auto ctx_a = make_ctx({"x", "__amf_eta"});
    auto ctx_b = make_ctx({"x"});

    // p = x + 2*x*__amf_eta + 3
    // Project: the middle monomial drops because it touches __amf_eta.
    auto p_a = alg::Mpoly::from_string(ctx_a, "x + 2*x*__amf_eta + 3");
    auto p_b = alg::project_mpoly_dropping(p_a, ctx_b);
    auto expected = alg::Mpoly::from_string(ctx_b, "x + 3");
    EXPECT_EQ(p_b, expected);
}
