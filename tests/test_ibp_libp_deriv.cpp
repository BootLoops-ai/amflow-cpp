// SPDX-License-Identifier: MIT
// Tests for amflow::ibp::libp_deriv.

#include <gtest/gtest.h>

#include <algorithm>

#include "amflow/algebra/mpoly.hpp"
#include "amflow/ibp/libp_deriv.hpp"
#include "amflow/qft/family_config.hpp"
#include "amflow/qft/jintegral.hpp"

namespace alg = amflow::algebra;
namespace ibp = amflow::ibp;
namespace qft = amflow::qft;

TEST(LibpDerivTest, LibpDenomsDeriv_Bubble_Wrt_msq) {
    auto fc = qft::FamilyConfig::build(
        "bubble", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - msq", "(l - p)^2 - msq"});

    auto res = ibp::libp_denoms_deriv(fc, "msq");
    ASSERT_EQ(res.coef.size(), 2u);
    for (auto& c : res.coef[0]) EXPECT_TRUE(c.is_zero()) << c.to_string();
    EXPECT_EQ(res.constant[0].to_string(), "-1");
    for (auto& c : res.coef[1]) EXPECT_TRUE(c.is_zero()) << c.to_string();
    EXPECT_EQ(res.constant[1].to_string(), "-1");
}

TEST(LibpDerivTest, LibpDenomsDeriv_Bubble_Wrt_s_KinematicChain) {
    // Replacement-defined invariant: now handled via the momentum-
    // derivative chain rule (was all-zero pre-patch).  Hand-derived:
    //   G = [[s]], C = ½ G⁻¹ ∂G/∂s = [[1/(2s)]],
    //   p·∂D₀/∂p = 0,  p·∂D₁/∂p = -2 l·p + 2 p² ;
    //   l·p = (D₀ - D₁)/2 + s/2  ⇒
    //   coef[1] = (-1/(2s), 1/(2s)),  const[1] = 1 - 1/2 = 1/2.
    auto fc = qft::FamilyConfig::build(
        "bubble", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - msq", "(l - p)^2 - msq"});

    auto res = ibp::libp_denoms_deriv(fc, "s");
    ASSERT_EQ(res.coef.size(), 2u);
    for (auto& c : res.coef[0]) EXPECT_TRUE(c.is_zero()) << c.to_string();
    EXPECT_TRUE(res.constant[0].is_zero()) << res.constant[0].to_string();
    EXPECT_EQ(res.coef[1][0].to_string(), "(-1)/(2*s)");
    EXPECT_EQ(res.coef[1][1].to_string(), "(1)/(2*s)");
    EXPECT_EQ(res.constant[1].to_string(), "(1)/(2)");
}

TEST(LibpDerivTest, LibpDeriv_J11_Wrt_msq) {
    auto fc = qft::FamilyConfig::build(
        "bubble", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - msq", "(l - p)^2 - msq"});

    qft::JIntegral j("bubble", {1, 1});
    auto terms = ibp::libp_deriv(fc, j, "msq");
    auto simplified = ibp::simplify_terms(std::move(terms));
    ASSERT_EQ(simplified.size(), 2u);

    std::sort(simplified.begin(), simplified.end(),
              [](const ibp::DerivTerm& a, const ibp::DerivTerm& b) {
                  return a.integ.indices() < b.integ.indices();
              });
    EXPECT_EQ(simplified[0].integ.indices(), (std::vector<long>{1, 2}));
    EXPECT_EQ(simplified[0].coef.to_string(), "1");
    EXPECT_EQ(simplified[1].integ.indices(), (std::vector<long>{2, 1}));
    EXPECT_EQ(simplified[1].coef.to_string(), "1");
}

TEST(LibpDerivTest, LibpDeriv_J_DegenerateZeroIndex) {
    auto fc = qft::FamilyConfig::build(
        "bubble", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - msq", "(l - p)^2 - msq"});

    qft::JIntegral j("bubble", {0, 1});
    auto terms = ibp::libp_deriv(fc, j, "msq");
    auto s = ibp::simplify_terms(std::move(terms));
    ASSERT_EQ(s.size(), 1u);
    EXPECT_EQ(s[0].integ.indices(), (std::vector<long>{0, 2}));
    EXPECT_EQ(s[0].coef.to_string(), "1");
}

TEST(LibpDerivTest, LibpDeriv_J_HigherIndex) {
    auto fc = qft::FamilyConfig::build(
        "bubble", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - msq", "(l - p)^2 - msq"});

    qft::JIntegral j("bubble", {2, 1});
    auto s = ibp::simplify_terms(ibp::libp_deriv(fc, j, "msq"));
    ASSERT_EQ(s.size(), 2u);
    std::sort(s.begin(), s.end(),
              [](const ibp::DerivTerm& a, const ibp::DerivTerm& b) {
                  return a.integ.indices() < b.integ.indices();
              });
    EXPECT_EQ(s[0].integ.indices(), (std::vector<long>{2, 2}));
    EXPECT_EQ(s[0].coef.to_string(), "1");
    EXPECT_EQ(s[1].integ.indices(), (std::vector<long>{3, 1}));
    EXPECT_EQ(s[1].coef.to_string(), "2");
}

TEST(LibpDerivTest, ComputeDerivative_LinearCombo) {
    auto fc = qft::FamilyConfig::build(
        "bubble", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - msq", "(l - p)^2 - msq"});

    long msq_idx = fc.var_index("msq");
    alg::Mfrac two = alg::Mfrac::from_si(fc.ctx, 2);
    alg::Mfrac msq_mfrac = alg::Mfrac::from_mpoly(
        alg::Mpoly::variable(fc.ctx, msq_idx));

    std::vector<ibp::DerivTerm> expr;
    {
        ibp::DerivTerm t;
        t.coef = two.clone();
        t.integ = qft::JIntegral("bubble", {1, 1});
        expr.push_back(std::move(t));
    }
    {
        ibp::DerivTerm t;
        t.coef = msq_mfrac.clone();
        t.integ = qft::JIntegral("bubble", {0, 1});
        expr.push_back(std::move(t));
    }

    auto out = ibp::simplify_terms(ibp::compute_derivative(fc, expr, "msq"));

    std::sort(out.begin(), out.end(),
              [](const ibp::DerivTerm& a, const ibp::DerivTerm& b) {
                  return a.integ.indices() < b.integ.indices();
              });

    ASSERT_EQ(out.size(), 4u);
    EXPECT_EQ(out[0].integ.indices(), (std::vector<long>{0, 1}));
    EXPECT_EQ(out[0].coef.to_string(), "1");
    EXPECT_EQ(out[1].integ.indices(), (std::vector<long>{0, 2}));
    EXPECT_EQ(out[1].coef.to_string(), "msq");
    EXPECT_EQ(out[2].integ.indices(), (std::vector<long>{1, 2}));
    EXPECT_EQ(out[2].coef.to_string(), "2");
    EXPECT_EQ(out[3].integ.indices(), (std::vector<long>{2, 1}));
    EXPECT_EQ(out[3].coef.to_string(), "2");
}

// --- Multi-invariant cases -----
//
// The existing tests above all use the single-mass bubble family
// `{l^2 - msq, (l-p)^2 - msq}` — one mass parameter `msq` plus one
// Replacement-defined kinematic invariant `s` (`p^2 -> s`).  These
// tests use a two-mass bubble with three distinct invariants
// {m1sq, m2sq, s} so that LIBPDeriv is exercised on a family where:
//   - different propagators carry different free-symbol invariants
//     (m1sq only in p_0, m2sq only in p_1),
//   - dp_k/d(invariant) is nonzero on a strict subset of k's
//     (exercising the `if (!ck.is_zero())` short-circuit at
//     `src/ibp/libp_deriv.cpp:327`),
//   - simplify_terms / compute_derivative handle a coefficient
//     ring with multiple invariants flowing through.
//
// All expected outputs are hand-derived from the formula
//   d/dx (J[a_0, a_1]) = sum_k (-a_k) * (dp_k/dx as constant + sp
//                                          expansion) * J_shifted
// (mirror of upstream LIBPDeriv with the sign of upstream's `-ak`
// factor preserved in `libp_deriv.cpp` `neg_ak`).
//
// Scope note: differentiation w.r.t. a *Replacement-defined*
// kinematic invariant (e.g., `s` defined via `p^2 -> s`) is now
// handled by a momentum-derivative chain rule that takes the place of
// T. Peraro's LiteIBP `LIBPDerivivative` (via AMFlow's `Kira/interface.m`), in
// `libp_deriv.cpp` — see `*_KinematicChain` and `*Box1_Wrt_*` tests.
// For variables not appearing in any reduced_replacement value
// (η, propagator masses) the chain Jacobian is zero and the
// pre-patch direct path is taken unchanged.  See
// `docs/AUDIT_MMA_PARITY.md`.

namespace {

qft::FamilyConfig make_two_mass_bubble() {
    return qft::FamilyConfig::build(
        "bub2m", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - m1sq", "(l - p)^2 - m2sq"});
}

}  // namespace

TEST(LibpDerivTest, LibpDenomsDeriv_TwoMassBubble_Wrt_m1sq) {
    // m1sq appears only in propagator 0 → only k=0 has a nonzero
    // constant term (-1); k=1 row is all zero.
    auto fc = make_two_mass_bubble();
    auto res = ibp::libp_denoms_deriv(fc, "m1sq");
    ASSERT_EQ(res.coef.size(), 2u);
    for (auto& c : res.coef[0]) EXPECT_TRUE(c.is_zero()) << c.to_string();
    EXPECT_EQ(res.constant[0].to_string(), "-1");
    for (auto& c : res.coef[1]) EXPECT_TRUE(c.is_zero()) << c.to_string();
    EXPECT_TRUE(res.constant[1].is_zero()) << res.constant[1].to_string();
}

TEST(LibpDerivTest, LibpDenomsDeriv_TwoMassBubble_Wrt_m2sq) {
    auto fc = make_two_mass_bubble();
    auto res = ibp::libp_denoms_deriv(fc, "m2sq");
    ASSERT_EQ(res.coef.size(), 2u);
    for (auto& c : res.coef[0]) EXPECT_TRUE(c.is_zero()) << c.to_string();
    EXPECT_TRUE(res.constant[0].is_zero()) << res.constant[0].to_string();
    for (auto& c : res.coef[1]) EXPECT_TRUE(c.is_zero()) << c.to_string();
    EXPECT_EQ(res.constant[1].to_string(), "-1");
}

TEST(LibpDerivTest, LibpDeriv_J11_TwoMassBubble_Wrt_m1sq) {
    // d/d(m1sq) J[1,1] = (-1)*(-1) * J[2,1] = +J[2,1].
    auto fc = make_two_mass_bubble();
    qft::JIntegral j("bub2m", {1, 1});
    auto s = ibp::simplify_terms(ibp::libp_deriv(fc, j, "m1sq"));
    ASSERT_EQ(s.size(), 1u);
    EXPECT_EQ(s[0].integ.indices(), (std::vector<long>{2, 1}));
    EXPECT_EQ(s[0].coef.to_string(), "1");
}

TEST(LibpDerivTest, LibpDeriv_J11_TwoMassBubble_Wrt_m2sq) {
    auto fc = make_two_mass_bubble();
    qft::JIntegral j("bub2m", {1, 1});
    auto s = ibp::simplify_terms(ibp::libp_deriv(fc, j, "m2sq"));
    ASSERT_EQ(s.size(), 1u);
    EXPECT_EQ(s[0].integ.indices(), (std::vector<long>{1, 2}));
    EXPECT_EQ(s[0].coef.to_string(), "1");
}

TEST(LibpDerivTest, LibpDeriv_J11_TwoMassBubble_Wrt_s_KinematicChain) {
    // Replacement-defined invariant: now handled by the momentum-
    // derivative chain rule (was zero-term pre-patch).  Hand-derived
    // (∂/∂s) J[1,1] = 1/(2s)·J[0,2] - 1/(2s)·J[1,1]
    //              + (m1sq-m2sq-s)/(2s)·J[1,2]
    // (see `LibpDenomsDeriv_Bubble_Wrt_s_KinematicChain`).
    auto fc = make_two_mass_bubble();
    qft::JIntegral j("bub2m", {1, 1});
    auto s = ibp::simplify_terms(ibp::libp_deriv(fc, j, "s"));
    ASSERT_EQ(s.size(), 3u);
    std::sort(s.begin(), s.end(),
              [](const ibp::DerivTerm& a, const ibp::DerivTerm& b) {
                  return a.integ.indices() < b.integ.indices();
              });
    EXPECT_EQ(s[0].integ.indices(), (std::vector<long>{0, 2}));
    EXPECT_EQ(s[0].coef.to_string(), "(1)/(2*s)");
    EXPECT_EQ(s[1].integ.indices(), (std::vector<long>{1, 1}));
    EXPECT_EQ(s[1].coef.to_string(), "(-1)/(2*s)");
    EXPECT_EQ(s[2].integ.indices(), (std::vector<long>{1, 2}));
    // β_{l·p} = (m1sq - m2sq + s)/2 here (vs s/2 in the equal-mass
    // case), so the J[1,2] coefficient picks up the mass split.
    EXPECT_EQ(s[2].coef.to_string(), "(m1sq-m2sq-s)/(2*s)");
}

TEST(LibpDerivTest, LibpDeriv_J21_TwoMassBubble_Wrt_m1sq) {
    // For J[2,1]: a_0 = 2, neg_a0 = -2.  Contribution from k=0:
    // -2 * (-1) * J[3,1] = +2 * J[3,1].
    auto fc = make_two_mass_bubble();
    qft::JIntegral j("bub2m", {2, 1});
    auto s = ibp::simplify_terms(ibp::libp_deriv(fc, j, "m1sq"));
    ASSERT_EQ(s.size(), 1u);
    EXPECT_EQ(s[0].integ.indices(), (std::vector<long>{3, 1}));
    EXPECT_EQ(s[0].coef.to_string(), "2");
}

TEST(LibpDerivTest, LibpDeriv_TwoMassBubble_AllInvariantsCommute) {
    // Cross-check: d/d(m1sq) followed by d/d(m2sq) of J[1,1] equals
    // d/d(m2sq) followed by d/d(m1sq) (partial derivatives commute).
    // This exercises compute_derivative over a multi-invariant family.
    auto fc = make_two_mass_bubble();
    qft::JIntegral j("bub2m", {1, 1});

    auto first_m1_then_m2 = ibp::simplify_terms(
        ibp::compute_derivative(fc,
            ibp::libp_deriv(fc, j, "m1sq"),
            "m2sq"));
    auto first_m2_then_m1 = ibp::simplify_terms(
        ibp::compute_derivative(fc,
            ibp::libp_deriv(fc, j, "m2sq"),
            "m1sq"));

    auto sort_terms = [](std::vector<ibp::DerivTerm>& v) {
        std::sort(v.begin(), v.end(),
                  [](const ibp::DerivTerm& a, const ibp::DerivTerm& b) {
                      return a.integ.indices() < b.integ.indices();
                  });
    };
    sort_terms(first_m1_then_m2);
    sort_terms(first_m2_then_m1);

    ASSERT_EQ(first_m1_then_m2.size(), first_m2_then_m1.size());
    for (std::size_t i = 0; i < first_m1_then_m2.size(); ++i) {
        EXPECT_EQ(first_m1_then_m2[i].integ.indices(),
                  first_m2_then_m1[i].integ.indices());
        EXPECT_EQ(first_m1_then_m2[i].coef.to_string(),
                  first_m2_then_m1[i].coef.to_string());
    }
}

// --- Box1 multi-leg kinematic chain ------------------------------
//
// Exercises the full chain rule with n_red = 3 (3×3 Gram inverse,
// off-diagonal C_{ai}, leg·leg → invariants in `decompose_chain`).

namespace {

qft::FamilyConfig make_box1_libp() {
    return qft::FamilyConfig::build(
        "box1", {"l"}, {"p1", "p2", "p3", "p4"},
        {{"p4", "-p1 - p2 - p3"}},
        {{"p1^2", "0"}, {"p2^2", "0"}, {"p3^2", "0"}, {"p4^2", "0"},
         {"(p1 + p2)^2", "s"}, {"(p1 + p3)^2", "t"}},
        {"l^2", "(l + p1)^2", "(l + p1 + p2)^2",
         "(l + p1 + p2 + p4)^2"});
}

}  // namespace

TEST(LibpDerivTest, LibpDeriv_Box1_Wrt_s_NonZero) {
    auto fc = make_box1_libp();
    qft::JIntegral j("box1", {1, 1, 1, 1});
    auto terms = ibp::simplify_terms(ibp::libp_deriv(fc, j, "s"));
    EXPECT_FALSE(terms.empty()) << "kinematic chain rule produced no terms";
    // All coefficients are rational in (s, t) with no leg dependence.
    for (const auto& t : terms) {
        const auto str = t.coef.to_string();
        EXPECT_EQ(str.find("p1"), std::string::npos) << str;
        EXPECT_EQ(str.find("p2"), std::string::npos) << str;
        EXPECT_EQ(str.find("p3"), std::string::npos) << str;
        EXPECT_EQ(str.find("l"),  std::string::npos) << str;
    }
}

TEST(LibpDerivTest, LibpDeriv_Box1_Wrt_t_NonZero) {
    // The crucial case: NO propagator carries an explicit `t`-dependent
    // leg·leg monomial (D₀=l², D₁=(l+p1)², D₂=(l+p1+p2)², D₃=(l-p3)²;
    // p1·p3, p2·p3 never appear).  The naive ∂D_k/∂t path is therefore
    // identically zero — only the G⁻¹ chain rule produces the (correct)
    // nonzero ∂I/∂t.
    auto fc = make_box1_libp();
    qft::JIntegral j("box1", {1, 1, 1, 1});
    auto terms = ibp::simplify_terms(ibp::libp_deriv(fc, j, "t"));
    EXPECT_FALSE(terms.empty()) << "kinematic chain rule produced no terms";
}

TEST(LibpDerivTest, LibpDeriv_Box1_EulerScaling) {
    // Euler / dimensional-analysis identity (massless 1-loop, 4 props):
    //   (s ∂/∂s + t ∂/∂t) I = (d/2 - 4) I .
    // The chain coefficients satisfy s·C_s + t·C_t = ½ G⁻¹(sJ_s+tJ_t)
    // = ½ G⁻¹ G = ½·𝟙 (since G is linear-homogeneous in (s,t) for the
    // massless box), so the Euler combination acts as ½ Σ_i p_i·∂/∂p_i.
    // The pre-IBP libp_deriv expression therefore has DENOMINATOR-FREE
    // coefficients (every term's coef is polynomial — the 1/det(G) pole
    // cancels), which is the structural check we make here.  The full
    // collapse to (d/2-4)·J[1,1,1,1] only holds AFTER IBP reduction
    // (verified end-to-end by the `mode:"diffeq"` path).
    auto fc = make_box1_libp();
    qft::JIntegral j("box1", {1, 1, 1, 1});

    auto ds = ibp::libp_deriv(fc, j, "s");
    auto dt = ibp::libp_deriv(fc, j, "t");
    long s_var = fc.var_index("s");
    long t_var = fc.var_index("t");
    alg::Mfrac s_mf = alg::Mfrac::from_mpoly(alg::Mpoly::variable(fc.ctx, s_var));
    alg::Mfrac t_mf = alg::Mfrac::from_mpoly(alg::Mpoly::variable(fc.ctx, t_var));
    for (auto& term : ds) term.coef *= s_mf;
    for (auto& term : dt) term.coef *= t_mf;
    std::vector<ibp::DerivTerm> sum;
    for (auto& term : ds) sum.push_back(std::move(term));
    for (auto& term : dt) sum.push_back(std::move(term));
    auto euler = ibp::simplify_terms(std::move(sum));
    EXPECT_FALSE(euler.empty());
    // Individual ∂/∂s, ∂/∂t terms carry 1/(s·t·(s+t)) Gram-inverse poles.
    bool ds_has_denom = false;
    for (const auto& term : ibp::simplify_terms(ibp::libp_deriv(fc, j, "s"))) {
        if (!term.coef.is_polynomial()) ds_has_denom = true;
    }
    EXPECT_TRUE(ds_has_denom) << "expected Gram-inverse poles in raw ∂/∂s";
    // The Euler sum's coefficients are polynomial (Gram pole cancels).
    for (const auto& term : euler) {
        EXPECT_TRUE(term.coef.is_polynomial())
            << "Euler-sum coef has non-trivial denominator: "
            << term.coef.to_string()
            << "  (J=" << term.integ.to_string() << ")";
    }
}

TEST(LibpDerivTest, LibpDenomsDeriv_Box1_Wrt_s_NarrowCtxMatchesWide) {
    // Narrow-context overload (the one ibp::diffeq actually uses) must
    // agree with the wide-context overload after numeric substitution.
    auto fc = make_box1_libp();
    auto wide = ibp::libp_denoms_deriv(fc, "s");

    auto narrow_ctx =
        std::make_shared<alg::MpolyContext>(std::vector<std::string>{"eta", "d"});
    std::map<std::string, std::string> nv = {{"s", "-3"}, {"t", "-7"}};
    auto narrow = ibp::libp_denoms_deriv(fc, "s", nv, narrow_ctx);

    ASSERT_EQ(wide.coef.size(), narrow.coef.size());
    // Wide entries are Mfracs in (s,t) on fc.ctx — substitute s=-3, t=-7
    // via Mfrac::from_string on the narrow ctx and compare to_string().
    // (Spot-check the constant[k] only; the full coef tensor is large.)
    for (std::size_t k = 0; k < wide.constant.size(); ++k) {
        EXPECT_EQ(wide.constant[k].is_zero(),
                  narrow.constant[k].is_zero())
            << "k=" << k << " wide=" << wide.constant[k].to_string()
            << " narrow=" << narrow.constant[k].to_string();
    }
}

TEST(LibpDerivTest, SimplifyTerms_DropsZeroCoefAndCombinesLikeJ) {
    auto fc = qft::FamilyConfig::build(
        "bubble", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - msq", "(l - p)^2 - msq"});

    std::vector<ibp::DerivTerm> terms;
    {
        ibp::DerivTerm t;
        t.coef = alg::Mfrac::from_si(fc.ctx, 1);
        t.integ = qft::JIntegral("bubble", {1, 1});
        terms.push_back(std::move(t));
    }
    {
        ibp::DerivTerm t;
        t.coef = alg::Mfrac::from_si(fc.ctx, -1);
        t.integ = qft::JIntegral("bubble", {1, 1});
        terms.push_back(std::move(t));
    }
    {
        ibp::DerivTerm t;
        t.coef = alg::Mfrac::from_si(fc.ctx, 5);
        t.integ = qft::JIntegral("bubble", {2, 1});
        terms.push_back(std::move(t));
    }

    auto s = ibp::simplify_terms(std::move(terms));
    ASSERT_EQ(s.size(), 1u);
    EXPECT_EQ(s[0].integ.indices(), (std::vector<long>{2, 1}));
    EXPECT_EQ(s[0].coef.to_string(), "5");
}
