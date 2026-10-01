// SPDX-License-Identifier: MIT
// Tests for amflow::qft::region.

#include <gtest/gtest.h>

#include <flint/fmpz_mpoly.h>

#include "amflow/algebra/mpoly.hpp"
#include "amflow/qft/family_config.hpp"
#include "amflow/qft/jintegral.hpp"
#include "amflow/qft/region.hpp"

namespace alg = amflow::algebra;
namespace qft = amflow::qft;

namespace {

long max_degree_in_var_for_test(const alg::Mpoly& p, long var) {
    auto ctx = p.ctx();
    long len = fmpz_mpoly_length(p.raw(), ctx->raw());
    std::vector<unsigned long> exp((std::size_t)ctx->n_vars());
    long d = 0;
    for (long t = 0; t < len; ++t) {
        fmpz_mpoly_get_term_exp_ui(exp.data(), p.raw(), t, ctx->raw());
        if ((long)exp[(std::size_t)var] > d) d = (long)exp[(std::size_t)var];
    }
    return d;
}

}  // namespace

TEST(RegionTest, RegionContext_HasEtaAndHalfEta) {
    auto fc = qft::FamilyConfig::build(
        "bubble", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - msq", "(l - p)^2 - msq"});
    auto rctx = qft::make_region_context(fc);

    EXPECT_EQ(rctx.eta_var, fc.ctx->n_vars());
    EXPECT_EQ(rctx.half_eta_var, fc.ctx->n_vars() + 1);
    EXPECT_EQ(rctx.ctx->n_vars(), fc.ctx->n_vars() + 2);
    EXPECT_EQ(rctx.ctx->var_name(rctx.eta_var), "__amf_eta");
    EXPECT_EQ(rctx.ctx->var_name(rctx.half_eta_var), "__amf_half_eta");
}

TEST(RegionTest, BranchMomenta_OneLoopBubble) {
    auto fc = qft::FamilyConfig::build(
        "bubble", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - msq", "(l - p)^2 - msq"});
    auto br = qft::branch_momenta(fc, fc.propagators_after_conservation);

    ASSERT_EQ(br.size(), 2u);
    EXPECT_EQ(br[0].to_string(), "l");
    EXPECT_EQ(br[1].to_string(), "l");
}

TEST(RegionTest, BranchMomenta_Sunrise_Distinct) {
    auto fc = qft::FamilyConfig::build(
        "sunrise", {"l1", "l2"}, {"p"}, {}, {{"p^2", "s"}},
        {"l1^2 - m1sq", "l2^2 - m2sq", "(l1 + l2 - p)^2 - m3sq"});
    auto br = qft::branch_momenta(fc, fc.propagators_after_conservation);
    ASSERT_EQ(br.size(), 3u);
    EXPECT_EQ(br[0].to_string(), "l1");
    EXPECT_EQ(br[1].to_string(), "l2");
    EXPECT_EQ(br[2].to_string(), "l1+l2");
}

TEST(RegionTest, BranchToLoop_Invertible_Identity) {
    auto fc = qft::FamilyConfig::build(
        "two_loop", {"l1", "l2"}, {}, {}, {}, {"l1^2", "l2^2"});
    auto rctx = qft::make_region_context(fc);
    auto cands = qft::branch_momenta(fc, fc.propagators_after_conservation);
    ASSERT_EQ(cands.size(), 2u);
    auto tr = qft::branch_to_loop(fc, rctx, cands);
    ASSERT_TRUE(tr.ok);
    ASSERT_EQ(tr.map.size(), 2u);
    EXPECT_EQ(tr.map[0].numerator().to_string(), "l1");
    EXPECT_EQ(tr.map[0].denominator().to_string(), "1");
    EXPECT_EQ(tr.map[1].numerator().to_string(), "l2");
    EXPECT_EQ(tr.map[1].denominator().to_string(), "1");
}

TEST(RegionTest, BranchToLoop_Singular_FailsCleanly) {
    auto fc = qft::FamilyConfig::build(
        "two_loop", {"l1", "l2"}, {}, {}, {}, {"l1^2", "(l1 + l2)^2"});
    auto rctx = qft::make_region_context(fc);
    auto cands = qft::branch_momenta(fc, fc.propagators_after_conservation);
    std::vector<alg::Mpoly> dup;
    dup.push_back(cands[0].clone());
    dup.push_back(cands[0].clone());
    auto tr = qft::branch_to_loop(fc, rctx, dup);
    EXPECT_FALSE(tr.ok);
}

TEST(RegionTest, BranchToLoop_Invertible_MixedBasisMatchesMathematica) {
    auto fc = qft::FamilyConfig::build(
        "two_loop", {"l1", "l2"}, {}, {}, {}, {"l1^2", "l2^2"});
    auto rctx = qft::make_region_context(fc);

    std::vector<alg::Mpoly> cands;
    cands.push_back(alg::Mpoly::variable(fc.ctx, 1));
    cands.push_back(alg::Mpoly::variable(fc.ctx, 0) + alg::Mpoly::variable(fc.ctx, 1));

    auto tr = qft::branch_to_loop(fc, rctx, cands);
    ASSERT_TRUE(tr.ok);
    ASSERT_EQ(tr.map.size(), 2u);
    EXPECT_EQ(tr.map[0].to_string(), "-l1+l2");
    EXPECT_EQ(tr.map[1].to_string(), "l1");
}

TEST(RegionTest, ApplyRegionRule_OneLoopBubble_AllSmall) {
    auto fc = qft::FamilyConfig::build(
        "bubble", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - msq", "(l - p)^2 - msq"});
    auto rctx = qft::make_region_context(fc);
    auto cands = qft::branch_momenta(fc, fc.propagators_after_conservation);
    std::vector<alg::Mpoly> single_cand;
    single_cand.push_back(cands[0].clone());
    auto tr = qft::branch_to_loop(fc, rctx, single_cand);
    ASSERT_TRUE(tr.ok);

    auto rule = qft::region_rule(rctx, tr, {0});
    alg::Mpoly D2 = fc.propagators_after_conservation[1].clone();
    auto t = qft::apply_region_rule(fc, rctx, D2, rule);
    auto num = t.numerator();
    EXPECT_FALSE(num.is_zero());
}

TEST(RegionTest, ApplyRegionRule_OneLoopBubble_AllLarge) {
    auto fc = qft::FamilyConfig::build(
        "bubble", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - msq", "(l - p)^2 - msq"});
    auto rctx = qft::make_region_context(fc);
    auto cands = qft::branch_momenta(fc, fc.propagators_after_conservation);
    std::vector<alg::Mpoly> single_cand;
    single_cand.push_back(cands[0].clone());
    auto tr = qft::branch_to_loop(fc, rctx, single_cand);
    auto rule = qft::region_rule(rctx, tr, {1});

    alg::Mpoly D1 = fc.propagators_after_conservation[0].clone();
    auto t = qft::apply_region_rule(fc, rctx, D1, rule);
    auto num = t.numerator();
    long len = fmpz_mpoly_length(num.raw(), rctx.ctx->raw());
    bool found_he = false;
    std::vector<unsigned long> exp((std::size_t)rctx.ctx->n_vars());
    for (long i = 0; i < len; ++i) {
        fmpz_mpoly_get_term_exp_ui(exp.data(), num.raw(), i, rctx.ctx->raw());
        if (exp[(std::size_t)rctx.half_eta_var] > 0) { found_he = true; break; }
    }
    EXPECT_TRUE(found_he);
}

TEST(RegionTest, BranchScale_AllZeroIfAllSmall) {
    auto fc = qft::FamilyConfig::build(
        "two_loop", {"l1", "l2"}, {}, {}, {}, {"l1^2", "l2^2"});
    auto rctx = qft::make_region_context(fc);
    auto cands = qft::branch_momenta(fc, fc.propagators_after_conservation);
    auto tr = qft::branch_to_loop(fc, rctx, cands);
    auto rule_id = qft::region_rule(rctx, tr, {0, 0});
    std::vector<alg::Mfrac> tb;
    for (auto& c : cands) tb.push_back(qft::apply_region_rule(fc, rctx, c, rule_id));
    auto bs = qft::branch_scale(fc, tb, {0, 0});
    ASSERT_EQ(bs.size(), 2u);
    EXPECT_EQ(bs[0], 0);
    EXPECT_EQ(bs[1], 0);
}

TEST(RegionTest, BranchScale_LargeWhenLoopMatches) {
    auto fc = qft::FamilyConfig::build(
        "two_loop", {"l1", "l2"}, {}, {}, {}, {"l1^2", "l2^2"});
    auto rctx = qft::make_region_context(fc);
    auto cands = qft::branch_momenta(fc, fc.propagators_after_conservation);
    auto tr = qft::branch_to_loop(fc, rctx, cands);
    auto rule_id = qft::region_rule(rctx, tr, {0, 0});
    std::vector<alg::Mfrac> tb;
    for (auto& c : cands) tb.push_back(qft::apply_region_rule(fc, rctx, c, rule_id));
    auto bs = qft::branch_scale(fc, tb, {1, 0});
    ASSERT_EQ(bs.size(), 2u);
    EXPECT_EQ(bs[0], 1);
    EXPECT_EQ(bs[1], 0);
}

TEST(RegionTest, FindAllRegion_OneLoopBubble) {
    auto fc = qft::FamilyConfig::build(
        "bubble", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - msq", "(l - p)^2 - msq"});
    auto rctx = qft::make_region_context(fc);
    auto regs = qft::find_all_region(fc, rctx, {0, 1});
    EXPECT_GE(regs.size(), 1u);
}

TEST(RegionTest, FindAllRegion_TwoLoopSeparated) {
    auto fc = qft::FamilyConfig::build(
        "two_loop", {"l1", "l2"}, {}, {}, {}, {"l1^2", "l2^2"});
    auto rctx = qft::make_region_context(fc);
    auto regs = qft::find_all_region(fc, rctx, {0, 1});
    EXPECT_GE(regs.size(), 1u);
}

TEST(RegionTest, FindAllRegion_SunriseInjectedSubfamilyMatchesMathematica) {
    auto fc = qft::FamilyConfig::build(
        "sunrise_sub", {"l1", "l2"}, {"p"}, {}, {{"p^2", "s"}},
        {"l2^2 - 1 - eta", "(l1 + l2)^2 - 1", "l1^2",
         "(l1 + p)^2", "(l2 + p)^2"});
    auto rctx = qft::make_region_context(fc);

    auto regs = qft::find_all_region(fc, rctx, {0, 1, 2});
    ASSERT_EQ(regs.size(), 5u);

    EXPECT_EQ(regs[0].transform.map[0].to_string(), "-l1+l2");
    EXPECT_EQ(regs[0].transform.map[1].to_string(), "l1");
    EXPECT_EQ(regs[0].scale, (std::vector<int>{0, 0}));

    EXPECT_EQ(regs[1].transform.map[0].to_string(), "-l1+l2");
    EXPECT_EQ(regs[1].transform.map[1].to_string(), "l1");
    EXPECT_EQ(regs[1].scale, (std::vector<int>{0, 1}));

    EXPECT_EQ(regs[2].transform.map[0].to_string(), "-l1+l2");
    EXPECT_EQ(regs[2].transform.map[1].to_string(), "l1");
    EXPECT_EQ(regs[2].scale, (std::vector<int>{1, 0}));

    EXPECT_EQ(regs[3].transform.map[0].to_string(), "-l1+l2");
    EXPECT_EQ(regs[3].transform.map[1].to_string(), "l1");
    EXPECT_EQ(regs[3].scale, (std::vector<int>{1, 1}));

    EXPECT_EQ(regs[4].transform.map[0].to_string(), "l2");
    EXPECT_EQ(regs[4].transform.map[1].to_string(), "l1");
    EXPECT_EQ(regs[4].scale, (std::vector<int>{1, 0}));

    EXPECT_TRUE(qft::zero_region_q(fc, rctx, regs[0], {0, 1, 2}));
    EXPECT_TRUE(qft::zero_region_q(fc, rctx, regs[1], {0, 1, 2}));
    EXPECT_FALSE(qft::zero_region_q(fc, rctx, regs[2], {0, 1, 2}));
    EXPECT_FALSE(qft::zero_region_q(fc, rctx, regs[3], {0, 1, 2}));
    EXPECT_TRUE(qft::zero_region_q(fc, rctx, regs[4], {0, 1, 2}));
}

TEST(RegionTest, ZeroRegionQ_TrivialPasses) {
    auto fc = qft::FamilyConfig::build(
        "bubble", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - msq", "(l - p)^2 - msq"});
    auto rctx = qft::make_region_context(fc);
    auto regs = qft::find_all_region(fc, rctx, {0, 1});
    ASSERT_FALSE(regs.empty());
    bool any_nonzero = false;
    for (auto& r : regs) {
        if (!qft::zero_region_q(fc, rctx, r, {0, 1})) {
            any_nonzero = true; break;
        }
    }
    EXPECT_TRUE(any_nonzero);
}

TEST(RegionTest, ZeroRegionQ_OneLoopBubble_MatchesMathematicaRegions) {
    auto fc = qft::FamilyConfig::build(
        "bubble", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - msq - eta", "(l - p)^2 - msq - eta"});
    auto rctx = qft::make_region_context(fc);
    auto regs = qft::find_all_region(fc, rctx, {0, 1});
    ASSERT_EQ(regs.size(), 2u);

    bool saw_small = false;
    bool saw_large = false;
    for (const auto& r : regs) {
        ASSERT_EQ(r.scale.size(), 1u);
        bool is_zero = qft::zero_region_q(fc, rctx, r, {0, 1});
        if (r.scale[0] == 0) {
            saw_small = true;
            EXPECT_TRUE(is_zero);
        } else if (r.scale[0] == 1) {
            saw_large = true;
            EXPECT_FALSE(is_zero);
        } else {
            FAIL() << "unexpected bubble region scale";
        }
    }
    EXPECT_TRUE(saw_small);
    EXPECT_TRUE(saw_large);
}

TEST(RegionTest, RegionPower_OneLoopBubble) {
    auto fc = qft::FamilyConfig::build(
        "bubble", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - msq", "(l - p)^2 - msq"});
    auto rctx = qft::make_region_context(fc);
    auto pctx = qft::make_powers_context(fc);
    auto regs = qft::find_all_region(fc, rctx, {0, 1});
    ASSERT_FALSE(regs.empty());

    qft::JIntegral integ("bubble", {1, 1});

    auto powers = qft::region_power(fc, rctx, pctx, regs.front(), {integ});
    ASSERT_EQ(powers.size(), 1u);

    alg::Mpoly num = powers[0].numerator();
    long deg_eps = max_degree_in_var_for_test(num, pctx.eps_var);
    EXPECT_LE(deg_eps, 1);
}

// --- `region_power` skips `/.Numeric` -----
//
// `region_power` (`src/qft/findregion.cpp:463-529`) constructs the
// per-integral exponent value as
//     val = sum_scale * (2 - eps) - sum_p
// where `sum_scale` and `sum_p` are integer accumulators (the
// scale-vector sum and the eta-bearing-propagator index sum).
// `eps` is the only symbol that survives — no kinematic
// invariants (s, t, msq, ...) ever enter the expression.
// Upstream applies `/. Numeric` after constructing this expression
// (`AMFlow.m`); since the expression contains nothing that
// `Numeric` could substitute, the upstream operation is a no-op
// for any well-formed Numeric.  The C++ omission is therefore
// provably equivalent.
//
// This test locks the precondition: every term in the output
// polynomial has all non-eps variables at exponent 0.

TEST(RegionTest, RegionPower_OutputIsEpsOnlyPlusIntegers_AuditRow200Equivalence) {
    // Use a multi-invariant family so the `/. Numeric` no-op claim
    // is testable with non-trivial Numeric values.  The output
    // expression should still contain only `eps` and integers.
    auto fc = qft::FamilyConfig::build(
        "box_for_audit200", {"l"}, {"p1", "p2", "p3", "p4"},
        {{"p4", "-p1 - p2 - p3"}},
        {{"p1^2", "0"}, {"p2^2", "0"}, {"p3^2", "0"},
         {"(p1 + p2)^2", "s"}, {"(p2 + p3)^2", "t"}},
        {"l^2", "(l + p1)^2", "(l + p1 + p2)^2",
         "(l + p1 + p2 + p3)^2"});
    auto rctx = qft::make_region_context(fc);
    auto pctx = qft::make_powers_context(fc);
    auto regs = qft::find_all_region(fc, rctx, {0, 1, 2, 3});
    ASSERT_FALSE(regs.empty());

    qft::JIntegral integ("box_for_audit200", {1, 1, 1, 1});
    auto powers = qft::region_power(fc, rctx, pctx, regs.front(), {integ});
    ASSERT_EQ(powers.size(), 1u);

    // Walk every term of numerator and denominator; every variable
    // other than `eps` must have exponent 0.
    auto only_eps = [&](const alg::Mpoly& p) -> bool {
        auto ctx = p.ctx();
        long len = fmpz_mpoly_length(p.raw(), ctx->raw());
        std::vector<unsigned long> exp((std::size_t)ctx->n_vars());
        for (long t = 0; t < len; ++t) {
            fmpz_mpoly_get_term_exp_ui(exp.data(), p.raw(), t, ctx->raw());
            for (long v = 0; v < ctx->n_vars(); ++v) {
                if (v == pctx.eps_var) continue;
                if (exp[(std::size_t)v] != 0) return false;
            }
        }
        return true;
    };

    EXPECT_TRUE(only_eps(powers[0].numerator()))
        << "region_power output numerator should be eps-only — applying"
        << " `/.Numeric` would be a no-op (equivalence).";
    EXPECT_TRUE(only_eps(powers[0].denominator()))
        << "region_power output denominator should be eps-only.";
}

TEST(RegionTest, ToCompleteExplicit_OneLoopBubble_AlreadyComplete) {
    auto fc = qft::FamilyConfig::build(
        "bubble", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - msq", "(l - p)^2 - msq"});
    auto complete = qft::to_complete_explicit(fc, fc.propagators_after_conservation);
    EXPECT_EQ(complete.size(), 2u);
    EXPECT_EQ(fc.sp_list.size(), 2u);
}

TEST(RegionTest, ToCompleteExplicit_PartialAugmented) {
    auto fc = qft::FamilyConfig::build(
        "tad_with_leg", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - msq"});
    EXPECT_EQ(fc.sp_list.size(), 2u);
    auto complete = qft::to_complete_explicit(fc, fc.propagators_after_conservation);
    EXPECT_EQ(complete.size(), 2u);
}

TEST(RegionTest, SPListToDListSymbol_OneLoopBubble) {
    auto fc = qft::FamilyConfig::build(
        "bubble", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - msq", "(l - p)^2 - msq"});
    auto complete = qft::to_complete_explicit(fc, fc.propagators_after_conservation);
    auto dctx = qft::make_dlist_context(fc, complete.size());
    auto sp_to_d = qft::sp_list_to_dlist_symbol(fc, dctx, complete);
    ASSERT_EQ(sp_to_d.size(), 2u);
    for (auto& f : sp_to_d) {
        EXPECT_FALSE(f.is_zero());
    }
}

TEST(RegionTest, SquaredDenominators_RoundTrip) {
    auto fc = qft::FamilyConfig::build(
        "bubble", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - msq", "(l - p)^2 - msq"});
    auto sq = qft::squared_denominators(fc, fc.propagators_after_conservation);
    ASSERT_EQ(sq.size(), 2u);
    EXPECT_TRUE((sq[0] - fc.propagators_after_conservation[0]).is_zero());
    EXPECT_TRUE((sq[1] - fc.propagators_after_conservation[1]).is_zero());
}
