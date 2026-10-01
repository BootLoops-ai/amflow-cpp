// SPDX-License-Identifier: MIT
// Tests for amflow::qft::boundary.

#include <gtest/gtest.h>

#include <flint/fmpz_mpoly.h>

#include "amflow/algebra/mpoly.hpp"
#include "amflow/qft/boundary.hpp"
#include "amflow/qft/family_config.hpp"
#include "amflow/qft/jintegral.hpp"
#include "amflow/qft/region.hpp"

namespace alg = amflow::algebra;
namespace qft = amflow::qft;

namespace {

std::shared_ptr<alg::MpolyContext> single_var_ctx() {
    return std::make_shared<alg::MpolyContext>(std::vector<std::string>{"eps"});
}

alg::Mfrac mk_aff(const std::shared_ptr<alg::MpolyContext>& ctx, long a, long b) {
    alg::Mpoly p = alg::Mpoly::constant(ctx, a)
                 + alg::Mpoly::variable(ctx, 0) * alg::Mpoly::constant(ctx, b);
    return alg::Mfrac::from_mpoly(std::move(p));
}

}  // namespace

// ============================================================
// boundary_pattern
// ============================================================

TEST(BoundaryTest, BoundaryPattern_SingleRegion) {
    auto ctx = single_var_ctx();
    std::vector<std::vector<alg::Mfrac>> powers;
    {
        std::vector<alg::Mfrac> v;
        v.push_back(mk_aff(ctx, 3, 2));
        v.push_back(mk_aff(ctx, 1, -1));
        powers.push_back(std::move(v));
    }
    auto res = qft::boundary_pattern(powers);
    ASSERT_EQ(res.size(), 1u);
    ASSERT_EQ(res[0].size(), 2u);
    EXPECT_EQ(res[0][0].to_string(), mk_aff(ctx, 3, 2).to_string());
    EXPECT_EQ(res[0][1].to_string(), mk_aff(ctx, 1, -1).to_string());
}

TEST(BoundaryTest, BoundaryPattern_TwoRegionsIntegerShift) {
    auto ctx = single_var_ctx();
    std::vector<std::vector<alg::Mfrac>> powers;
    {
        std::vector<alg::Mfrac> v;
        v.push_back(mk_aff(ctx, 3, 2));
        powers.push_back(std::move(v));
    }
    {
        std::vector<alg::Mfrac> v;
        v.push_back(mk_aff(ctx, 5, 2));
        powers.push_back(std::move(v));
    }
    auto res = qft::boundary_pattern(powers);
    ASSERT_EQ(res.size(), 1u);
    ASSERT_EQ(res[0].size(), 1u);
    EXPECT_EQ(res[0][0].to_string(), mk_aff(ctx, 5, 2).to_string());
}

TEST(BoundaryTest, BoundaryPattern_TwoRegionsDistinctClasses) {
    auto ctx = single_var_ctx();
    std::vector<std::vector<alg::Mfrac>> powers;
    {
        std::vector<alg::Mfrac> v;
        v.push_back(mk_aff(ctx, 3, 2));
        powers.push_back(std::move(v));
    }
    {
        std::vector<alg::Mfrac> v;
        v.push_back(mk_aff(ctx, 5, 1));
        powers.push_back(std::move(v));
    }
    auto res = qft::boundary_pattern(powers);
    ASSERT_EQ(res.size(), 2u);
    EXPECT_EQ(res[0].size(), 1u);
    EXPECT_EQ(res[1].size(), 1u);
}

// ============================================================
// apart_one_var
// ============================================================

TEST(BoundaryTest, ApartOneVar_SingleLinear) {
    auto ctx = std::make_shared<alg::MpolyContext>(
        std::vector<std::string>{"D"});
    alg::Mpoly D = alg::Mpoly::variable(ctx, 0);
    alg::Mpoly poly = D - alg::Mpoly::one(ctx);

    auto pieces = qft::apart_one_var(poly, 0);
    ASSERT_EQ(pieces.size(), 1u);
    alg::Mfrac expected = alg::Mfrac::from_mpoly(alg::Mpoly::one(ctx)) /
                          alg::Mfrac::from_mpoly(poly.clone());
    EXPECT_EQ(pieces[0].to_string(), expected.to_string());
}

TEST(BoundaryTest, ApartOneVar_TwoLinearFactors) {
    auto ctx = std::make_shared<alg::MpolyContext>(
        std::vector<std::string>{"D"});
    alg::Mpoly D = alg::Mpoly::variable(ctx, 0);
    alg::Mpoly Dm1 = D - alg::Mpoly::one(ctx);
    alg::Mpoly Dm2 = D - alg::Mpoly::constant(ctx, 2);
    alg::Mpoly poly = Dm1 * Dm2;

    auto pieces = qft::apart_one_var(poly, 0);
    ASSERT_EQ(pieces.size(), 2u);

    alg::Mfrac sum = alg::Mfrac::zero(ctx);
    for (auto& p : pieces) sum += p;
    alg::Mfrac expected =
        alg::Mfrac::one(ctx) / alg::Mfrac::from_mpoly(poly.clone());
    EXPECT_EQ(sum.to_string(), expected.to_string());
}

TEST(BoundaryTest, ApartOneVar_ThreeLinearFactors) {
    auto ctx = std::make_shared<alg::MpolyContext>(
        std::vector<std::string>{"D"});
    alg::Mpoly D = alg::Mpoly::variable(ctx, 0);
    alg::Mpoly poly = (D - alg::Mpoly::one(ctx))
                    * (D - alg::Mpoly::constant(ctx, 2))
                    * (D - alg::Mpoly::constant(ctx, 5));

    auto pieces = qft::apart_one_var(poly, 0);
    ASSERT_EQ(pieces.size(), 3u);

    alg::Mfrac sum = alg::Mfrac::zero(ctx);
    for (auto& p : pieces) sum += p;
    alg::Mfrac expected =
        alg::Mfrac::one(ctx) / alg::Mfrac::from_mpoly(poly.clone());
    EXPECT_EQ(sum.to_string(), expected.to_string());
}

TEST(BoundaryTest, ApartOneVar_NegativeRoot) {
    auto ctx = std::make_shared<alg::MpolyContext>(
        std::vector<std::string>{"D"});
    alg::Mpoly D = alg::Mpoly::variable(ctx, 0);
    alg::Mpoly poly = (D - alg::Mpoly::one(ctx))
                    * (D + alg::Mpoly::constant(ctx, 3));

    auto pieces = qft::apart_one_var(poly, 0);
    ASSERT_EQ(pieces.size(), 2u);
    alg::Mfrac sum = alg::Mfrac::zero(ctx);
    for (auto& p : pieces) sum += p;
    alg::Mfrac expected =
        alg::Mfrac::one(ctx) / alg::Mfrac::from_mpoly(poly.clone());
    EXPECT_EQ(sum.to_string(), expected.to_string());
}

TEST(BoundaryTest, ApartOneVar_PolyParameterAlpha) {
    auto ctx = std::make_shared<alg::MpolyContext>(
        std::vector<std::string>{"D", "a", "b"});
    alg::Mpoly D = alg::Mpoly::variable(ctx, 0);
    alg::Mpoly a = alg::Mpoly::variable(ctx, 1);
    alg::Mpoly b = alg::Mpoly::variable(ctx, 2);
    alg::Mpoly poly = (D - a) * (D - b);

    auto pieces = qft::apart_one_var(poly, 0);
    ASSERT_EQ(pieces.size(), 2u);
    alg::Mfrac sum = alg::Mfrac::zero(ctx);
    for (auto& p : pieces) sum += p;
    alg::Mfrac expected =
        alg::Mfrac::one(ctx) / alg::Mfrac::from_mpoly(poly.clone());
    EXPECT_EQ(sum.to_string(), expected.to_string());
}

TEST(BoundaryTest, ApartOneVar_NoFactorMentionsVar) {
    auto ctx = std::make_shared<alg::MpolyContext>(
        std::vector<std::string>{"D", "a"});
    alg::Mpoly a = alg::Mpoly::variable(ctx, 1);
    alg::Mpoly poly = a + alg::Mpoly::constant(ctx, 2);

    auto pieces = qft::apart_one_var(poly, 0);
    ASSERT_EQ(pieces.size(), 1u);
    alg::Mfrac expected =
        alg::Mfrac::one(ctx) / alg::Mfrac::from_mpoly(poly.clone());
    EXPECT_EQ(pieces[0].to_string(), expected.to_string());
}

TEST(BoundaryTest, ApartOneVar_SingleRepeatedFactorReturnsAsIs) {
    auto ctx = std::make_shared<alg::MpolyContext>(
        std::vector<std::string>{"D"});
    alg::Mpoly D = alg::Mpoly::variable(ctx, 0);
    alg::Mpoly poly = (D - alg::Mpoly::one(ctx))
                    * (D - alg::Mpoly::one(ctx));
    auto pieces = qft::apart_one_var(poly, 0);
    ASSERT_EQ(pieces.size(), 1u);
    alg::Mfrac expected =
        alg::Mfrac::one(ctx) / alg::Mfrac::from_mpoly(poly.clone());
    EXPECT_EQ(pieces[0].to_string(), expected.to_string());
}

TEST(BoundaryTest, ApartOneVar_RepeatedFactorWithDistinctRootMatchesSum) {
    auto ctx = std::make_shared<alg::MpolyContext>(
        std::vector<std::string>{"D"});
    alg::Mpoly D = alg::Mpoly::variable(ctx, 0);
    alg::Mpoly poly = (D - alg::Mpoly::one(ctx))
                    * (D - alg::Mpoly::one(ctx))
                    * (D + alg::Mpoly::one(ctx));
    auto pieces = qft::apart_one_var(poly, 0);
    ASSERT_EQ(pieces.size(), 3u);

    alg::Mfrac sum = alg::Mfrac::zero(ctx);
    for (auto& p : pieces) sum += p;
    alg::Mfrac expected =
        alg::Mfrac::one(ctx) / alg::Mfrac::from_mpoly(poly.clone());
    EXPECT_EQ(sum.to_string(), expected.to_string());
}

TEST(BoundaryTest, ApartOneVar_RepeatedFactorWithPolynomialAlphaMatchesSum) {
    auto ctx = std::make_shared<alg::MpolyContext>(
        std::vector<std::string>{"D", "a", "b"});
    alg::Mpoly D = alg::Mpoly::variable(ctx, 0);
    alg::Mpoly a = alg::Mpoly::variable(ctx, 1);
    alg::Mpoly b = alg::Mpoly::variable(ctx, 2);
    alg::Mpoly poly = (D - a) * (D - a) * (D - b);

    auto pieces = qft::apart_one_var(poly, 0);
    ASSERT_EQ(pieces.size(), 3u);

    alg::Mfrac sum = alg::Mfrac::zero(ctx);
    for (auto& p : pieces) sum += p;
    alg::Mfrac expected =
        alg::Mfrac::one(ctx) / alg::Mfrac::from_mpoly(poly.clone());
    EXPECT_EQ(sum.to_string(), expected.to_string());
}

// ============================================================
// apart_rationals
// ============================================================

TEST(BoundaryTest, ApartRationals_SingleRationalSingleD) {
    auto fc = qft::FamilyConfig::build(
        "tad", {"l"}, {}, {}, {}, {"l^2 - 1"});
    auto dctx = qft::make_dlist_context(fc, 1);

    alg::Mpoly D0 = alg::Mpoly::variable(dctx.ctx, dctx.first_d_var);
    alg::Mfrac r =
        alg::Mfrac::one(dctx.ctx)
        / alg::Mfrac::from_mpoly(D0 - alg::Mpoly::one(dctx.ctx));

    std::vector<alg::Mfrac> in;
    in.push_back(r.clone());
    auto out = qft::apart_rationals(in, dctx);
    ASSERT_EQ(out.size(), 1u);
    ASSERT_EQ(out[0].size(), 1u);

    alg::Mfrac sum = alg::Mfrac::zero(dctx.ctx);
    for (auto& p : out[0]) sum += p;
    EXPECT_EQ(sum.to_string(), r.to_string());
}

TEST(BoundaryTest, ApartRationals_TwoDVariables) {
    auto fc = qft::FamilyConfig::build(
        "two_loop", {"l1", "l2"}, {}, {}, {}, {"l1^2", "l2^2"});
    auto dctx = qft::make_dlist_context(fc, 2);

    alg::Mpoly D0 = alg::Mpoly::variable(dctx.ctx, dctx.first_d_var);
    alg::Mpoly D1 = alg::Mpoly::variable(dctx.ctx, dctx.first_d_var + 1);
    alg::Mfrac r =
        alg::Mfrac::one(dctx.ctx) /
        alg::Mfrac::from_mpoly((D0 - alg::Mpoly::one(dctx.ctx))
                              * (D1 - alg::Mpoly::constant(dctx.ctx, 2)));

    std::vector<alg::Mfrac> in;
    in.push_back(r.clone());
    auto out = qft::apart_rationals(in, dctx);

    alg::Mfrac sum = alg::Mfrac::zero(dctx.ctx);
    for (auto& p : out[0]) sum += p;
    EXPECT_EQ(sum.to_string(), r.to_string());
}

TEST(BoundaryTest, ApartRationals_CrossTermInOneD) {
    auto fc = qft::FamilyConfig::build(
        "two_loop", {"l1", "l2"}, {}, {}, {}, {"l1^2", "l2^2"});
    auto dctx = qft::make_dlist_context(fc, 2);

    alg::Mpoly D0 = alg::Mpoly::variable(dctx.ctx, dctx.first_d_var);
    alg::Mpoly D1 = alg::Mpoly::variable(dctx.ctx, dctx.first_d_var + 1);
    alg::Mpoly poly = (D0 - alg::Mpoly::one(dctx.ctx))
                    * (D0 - alg::Mpoly::constant(dctx.ctx, 2))
                    * D1;
    alg::Mfrac r =
        alg::Mfrac::one(dctx.ctx) / alg::Mfrac::from_mpoly(poly.clone());

    std::vector<alg::Mfrac> in;
    in.push_back(r.clone());
    auto out = qft::apart_rationals(in, dctx);

    alg::Mfrac sum = alg::Mfrac::zero(dctx.ctx);
    for (auto& p : out[0]) sum += p;
    EXPECT_EQ(sum.to_string(), r.to_string());
}

// ============================================================
// laporta_integrals
// ============================================================

TEST(BoundaryTest, LaportaIntegrals_PureD) {
    auto fc = qft::FamilyConfig::build(
        "two_loop", {"l1", "l2"}, {}, {}, {}, {"l1^2", "l2^2"});
    auto dctx = qft::make_dlist_context(fc, 2);

    alg::Mpoly D0 = alg::Mpoly::variable(dctx.ctx, dctx.first_d_var);
    alg::Mpoly D1 = alg::Mpoly::variable(dctx.ctx, dctx.first_d_var + 1);
    alg::Mfrac term =
        alg::Mfrac::one(dctx.ctx) / alg::Mfrac::from_mpoly(D0 * D1 * D1);

    auto out = qft::laporta_integrals(term, dctx);
    ASSERT_EQ(out.size(), 1u);
    EXPECT_EQ(out[0].indices, (std::vector<long>{1, 2}));
    EXPECT_EQ(out[0].coeff.to_string(),
              alg::Mfrac::one(dctx.ctx).to_string());
}

TEST(BoundaryTest, LaportaIntegrals_NegativeIndex) {
    auto fc = qft::FamilyConfig::build(
        "two_loop", {"l1", "l2"}, {}, {}, {}, {"l1^2", "l2^2"});
    auto dctx = qft::make_dlist_context(fc, 2);

    alg::Mpoly D0 = alg::Mpoly::variable(dctx.ctx, dctx.first_d_var);
    alg::Mpoly D1 = alg::Mpoly::variable(dctx.ctx, dctx.first_d_var + 1);
    alg::Mfrac term =
        alg::Mfrac::from_mpoly(D0 * D0)
        / alg::Mfrac::from_mpoly(D1.clone());

    auto out = qft::laporta_integrals(term, dctx);
    ASSERT_EQ(out.size(), 1u);
    EXPECT_EQ(out[0].indices, (std::vector<long>{-2, 1}));
    EXPECT_EQ(out[0].coeff.to_string(),
              alg::Mfrac::one(dctx.ctx).to_string());
}

TEST(BoundaryTest, LaportaIntegrals_NumeratorWithCoefficient) {
    auto fc = qft::FamilyConfig::build(
        "two_loop", {"l1", "l2"}, {}, {}, {}, {"l1^2", "l2^2"});
    auto dctx = qft::make_dlist_context(fc, 2);

    alg::Mpoly D0 = alg::Mpoly::variable(dctx.ctx, dctx.first_d_var);
    alg::Mpoly D1 = alg::Mpoly::variable(dctx.ctx, dctx.first_d_var + 1);
    alg::Mfrac term =
        alg::Mfrac::from_mpoly(alg::Mpoly::constant(dctx.ctx, 5))
        / alg::Mfrac::from_mpoly(D0 * D1);

    auto out = qft::laporta_integrals(term, dctx);
    ASSERT_EQ(out.size(), 1u);
    EXPECT_EQ(out[0].indices, (std::vector<long>{1, 1}));
    EXPECT_EQ(out[0].coeff.to_string(), "5");
}

// ============================================================
// boundary_integrands  (smoke test)
// ============================================================

TEST(BoundaryTest, BoundaryIntegrands_OneLoopMassiveBubble_RunsToCompletion) {
    auto fc = qft::FamilyConfig::build(
        "bubble", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - msq", "(l - p)^2 - msq"});
    auto rctx = qft::make_region_context(fc);

    auto regs = qft::find_all_region(fc, rctx, {0, 1});
    ASSERT_FALSE(regs.empty());

    qft::JIntegral integ("bubble", {1, 1});
    std::vector<qft::JIntegral> ints;
    ints.push_back(integ);
    std::vector<long> border = {0};

    auto bi = qft::boundary_integrands(fc, rctx, ints, border, regs.front());
    ASSERT_EQ(bi.integrands.size(), 1u);
    ASSERT_FALSE(bi.integrands[0].empty());
    EXPECT_FALSE(bi.completede.empty());
}

TEST(BoundaryTest, BoundaryIntegrands_RespectsBorder) {
    auto fc = qft::FamilyConfig::build(
        "bubble", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - msq", "(l - p)^2 - msq"});
    auto rctx = qft::make_region_context(fc);
    auto regs = qft::find_all_region(fc, rctx, {0, 1});
    ASSERT_FALSE(regs.empty());

    qft::JIntegral integ("bubble", {1, 1});
    std::vector<qft::JIntegral> ints;
    ints.push_back(integ);
    std::vector<long> border = {2};

    auto bi = qft::boundary_integrands(fc, rctx, ints, border, regs.front());
    ASSERT_EQ(bi.integrands.size(), 1u);
    EXPECT_LE(bi.integrands[0].size(), 3u);
    EXPECT_GE(bi.integrands[0].size(), 1u);
}

// ============================================================
// boundary_integrals (smoke test + sunrise parity)
// ============================================================

TEST(BoundaryTest, BoundaryIntegrals_RunsForBubble) {
    auto fc = qft::FamilyConfig::build(
        "bubble", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - msq", "(l - p)^2 - msq"});
    auto rctx = qft::make_region_context(fc);

    auto regs = qft::find_all_region(fc, rctx, {0, 1});
    ASSERT_FALSE(regs.empty());

    qft::JIntegral integ("bubble", {1, 1});
    std::vector<qft::JIntegral> ints;
    ints.push_back(integ);

    std::vector<long> border = {0};
    auto bi = qft::boundary_integrands(fc, rctx, ints, border, regs.front());
    auto fams = qft::boundary_integrals(fc, rctx, bi);

    EXPECT_FALSE(fams.empty());
}

TEST(BoundaryTest, SunriseTopBoundaryIntegrals_MatchMathBoundaryLayout) {
    auto fc_eta = qft::FamilyConfig::build(
        "sunrise", {"l1", "l2"}, {"p"}, {}, {{"p^2", "s"}},
        {"(l1^2 - msq) - eta",
         "(l2^2 - msq) - eta",
         "((l1 + l2 - p)^2 - msq) - eta",
         "(l1 + p)^2",
         "(l2 + p)^2"});

    auto rctx = qft::make_region_context(fc_eta);
    const std::vector<std::size_t> top_posi = {0, 1, 2};
    auto regs = qft::find_all_region(fc_eta, rctx, top_posi);
    ASSERT_FALSE(regs.empty());

    const qft::Region* non_zero = nullptr;
    for (const auto& reg : regs) {
        if (!qft::zero_region_q(fc_eta, rctx, reg, top_posi)) {
            ASSERT_EQ(non_zero, nullptr);
            non_zero = &reg;
        }
    }
    ASSERT_NE(non_zero, nullptr);

    std::vector<qft::JIntegral> integrals = {
        qft::JIntegral("sunrise", {1, 1, 0, 0, 0}),
        qft::JIntegral("sunrise", {1, 1, 1, 0, 0}),
        qft::JIntegral("sunrise", {1, 1, 1, -2, 0}),
    };
    const std::vector<long> border = {0, -1, 0};

    const auto bi = qft::boundary_integrands(fc_eta, rctx, integrals, border, *non_zero);
    ASSERT_EQ(bi.integrands.size(), 3u);
    ASSERT_EQ(bi.integrands[0].size(), 1u);
    ASSERT_EQ(bi.integrands[1].size(), 0u);
    ASSERT_EQ(bi.integrands[2].size(), 1u);

    const auto fams = qft::boundary_integrals(fc_eta, rctx, bi);
    ASSERT_EQ(fams.size(), 1u);
    ASSERT_EQ(fams[0].terms.size(), 3u);
    ASSERT_EQ(fams[0].terms[0].size(), 1u);
    ASSERT_EQ(fams[0].terms[1].size(), 0u);
    ASSERT_EQ(fams[0].terms[2].size(), 1u);

    ASSERT_EQ(fams[0].terms[0][0].size(), 1u);
    EXPECT_EQ(fams[0].terms[0][0][0].indices, (std::vector<long>{1, 1, 0, 0, 0}));
    EXPECT_EQ(fams[0].terms[0][0][0].coeff.to_string(), "1");

    ASSERT_EQ(fams[0].terms[2][0].size(), 3u);
    EXPECT_EQ(fams[0].terms[2][0][0].indices, (std::vector<long>{-1, 1, 1, 0, 0}));
    EXPECT_EQ(fams[0].terms[2][0][0].coeff.to_string(), "1");
    EXPECT_EQ(fams[0].terms[2][0][1].indices, (std::vector<long>{0, 1, 1, 0, 0}));
    EXPECT_EQ(fams[0].terms[2][0][1].coeff.to_string(), "2");
    EXPECT_EQ(fams[0].terms[2][0][2].indices, (std::vector<long>{1, 1, 1, 0, 0}));
    EXPECT_EQ(fams[0].terms[2][0][2].coeff.to_string(), "1");
}
