// SPDX-License-Identifier: MIT
// Tests for amflow::pipeline::amfsystem.

#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include <flint/acb.h>
#include <flint/arb.h>
#include <flint/arf.h>
#include <flint/fmpq.h>

#include "amflow/pipeline/amfsystem.hpp"
#include "amflow/numeric/acb_wrap.hpp"
#include "amflow/qft/family_config.hpp"
#include "amflow/qft/jintegral.hpp"
#include "amflow/qft/topology.hpp"

namespace pipeline = amflow::pipeline;
namespace qft     = amflow::qft;
namespace numeric = amflow::numeric;
namespace fs      = std::filesystem;

namespace {

bool kira_available() {
    return fs::exists("/usr/local/bin/kira") &&
           fs::exists("/usr/share/Ferl7/fer64");
}

}  // namespace

TEST(AmfsystemTest, EndingQ_SingleMassScheme_GeneralVacuumStopsImmediately) {
    auto fc = qft::FamilyConfig::build(
        "vacuum2", {"l1", "l2"}, {}, {}, {},
        {"l1^2 - 1", "l2^2 - 2"});

    pipeline::AMFSystemOptions opts;
    std::vector<qft::JIntegral> preferred{qft::JIntegral("vacuum2", {1, 1})};

    EXPECT_TRUE(pipeline::ending_q(fc, preferred,
                                   pipeline::EndingScheme::SingleMass, opts));
}

// --- SingleMass enhancement: numeric-substituted single-mass detection -----
//
// Pinned by SingleMassEnhancement_* tests.
//
// Upstream's SingleMassQ (`AMFlow.m`) tests `info.mass[k]` literally
// (`Count[mass, 1] === 1 && Count[mass, 0] === Length-1`).  C++
// `pipeline::single_mass_q_numeric` (`src/pipeline/amfsystem.cpp:803`)
// is an enhancement: it applies `numeric_values` to the mass list
// before the literal test.  This means the C++ port recognises a
// family like `{l^2, (l-p)^2 - msq}` with `numeric_values={msq: 1}`
// as single-mass (post-substitution mass list `[0, 1]`), while
// upstream rejects it because the literal mass list still contains
// the symbol `msq` rather than `1`.
//
// Consequence: C++ flows this family via SingleMass scheme, while
// upstream falls through to the next scheme (Cutkosky / Tradition /
// Trivial).  The final integral value is the same — SingleMass
// scheme is mathematically valid for the substituted single-mass
// case — but the dispatcher path differs.  No oracle bench exposes
// this dispatcher difference because all 12 oracles use families
// whose component analysis already produces literal `0` / `1`
// masses.
//
// These tests lock the C++ enhancement contract:
//   1. Without numeric: pipeline matches upstream (literal-only).
//   2. With numeric resolving the symbol to `1`: C++ enhancement
//      fires, scheme reported as applicable.

TEST(AmfsystemTest, SingleMassEnhancement_NoNumeric_MatchesUpstreamLiteral) {
    // Tadpole with a single free-symbol mass `msq`, NO numeric_values
    // supplied.  Upstream's literal SingleMassQ rejects (mass list is
    // [msq], neither 1 nor 0), so the system is "already ending under
    // SingleMass" (= scheme NOT applicable, ending_q returns true).
    // C++ pipeline with empty numeric_q reduces single_mass_q_numeric
    // to the same literal test → also reports already-ending.
    auto fc = qft::FamilyConfig::build(
        "tad_sym_a", {"l"}, {}, {}, {}, {"l^2 - msq"});

    pipeline::AMFSystemOptions opts;   // empty numeric_values
    std::vector<qft::JIntegral> preferred{qft::JIntegral("tad_sym_a", {1})};

    EXPECT_TRUE(pipeline::ending_q(fc, preferred,
                                   pipeline::EndingScheme::SingleMass, opts))
        << "without numeric_values, both literal-upstream and C++ pipeline"
        << " should report this tadpole as already-ending under SingleMass"
        << " (mass list is symbolic [msq], neither literally 1 nor 0).";
}

TEST(AmfsystemTest, SingleMassEnhancement_NumericResolvesSymbolToOne_FlowsViaSingleMass) {
    // Same tadpole + numeric_values = {msq: 1}.  Post-substitution
    // mass list is [1] — a literal single-mass shape (n_one = 1,
    // n_zero = 0 = Length - 1 ✓).  C++ enhancement detects this and
    // reports SingleMass scheme as applicable (ending_q returns
    // false).  Upstream's literal SingleMassQ on the un-substituted
    // mass list [msq] would still reject this case — this is the
    // load-bearing C++ enhancement contract.
    auto fc = qft::FamilyConfig::build(
        "tad_sym_b", {"l"}, {}, {}, {}, {"l^2 - msq"});

    pipeline::AMFSystemOptions opts;
    opts.bb.numeric_values["msq"] = "1";
    std::vector<qft::JIntegral> preferred{qft::JIntegral("tad_sym_b", {1})};

    EXPECT_FALSE(pipeline::ending_q(fc, preferred,
                                    pipeline::EndingScheme::SingleMass, opts))
        << "with numeric_values={msq:1}, the C++ enhancement should detect"
        << " the substituted [1] mass list as single-mass and report"
        << " SingleMass scheme as applicable (the literal upstream check"
        << " on [msq] would reject this case).";
}

// --- factorize_family mass=-1 detection enhancement
//
// Same enhancement pattern as SingleMassQ above: C++
// `find_mass_minus_one` (`src/pipeline/amfsystem.cpp:2394`) applies
// `numeric_values` to each mass entry before testing for literal
// `-1`, while upstream's `Position[..., -1]` (AMFlow.m
// SingleMassSetupMaster) tests literally.  This allows a family
// like `{l^2 - msq}` with `Numeric = {msq -> 1}` to flow through
// SingleMass scheme: after factorize_family the loop-promotion
// produces a propagator whose substituted mass is the literal
// `-1` the algorithm expects.  Without the enhancement, only
// families written with literal `1` would work; the enhancement
// transparently lifts symbolic mass parameters to the same
// behavior when Numeric resolves them.
//
// End-to-end test: run `pipeline::amf_system_setup_master` on
// the symbolic-mass tadpole with `numeric_values={msq:1}` and
// verify the SingleMass scheme produces a strictly-smaller
// ending sub-system (the test pattern matches
// `SingleMassSetup_TadpoleStrictlyShrinksToEndingSystem` above
// but with a symbolic-mass family + Numeric).

TEST(AmfsystemTest, FactorizeFamilyMassMinusOne_NumericResolvesSymbolic) {
    // Tadpole `{l^2 - msq}` with Numeric={msq: 1}.  The enhancement
    // path:
    //   1. ending_q(SingleMass) returns false (single_mass_q_numeric
    //      detects [1] as single-mass after substitution).
    //   2. single_mass_setup_master runs factorize_family on the
    //      legs-zeroed propagators.
    //   3. find_mass_minus_one substitutes msq -> 1 into the
    //      factorize_family-output mass list and identifies the
    //      drop propagator (whose post-loop-promotion mass is
    //      literal -1 after substitution).
    //   4. setup completes, producing an ending sub-system.
    auto fc = qft::FamilyConfig::build(
        "tad_sym_2b4", {"l"}, {}, {}, {}, {"l^2 - msq"});

    pipeline::AMFSystemOptions opts;
    opts.ending_schemes = {pipeline::EndingScheme::SingleMass};
    opts.bb.numeric_values["msq"] = "1";
    std::vector<qft::JIntegral> preferred{
        qft::JIntegral("tad_sym_2b4", {2})};

    auto systems = pipeline::amf_system_setup_master(fc, preferred, opts);
    ASSERT_EQ(systems.size(), 1u);

    const auto& child = *systems[0];
    EXPECT_TRUE(child.is_ending())
        << "with Numeric={msq:1}, SingleMass scheme should produce an"
        << " ending sub-system on this symbolic-mass tadpole.  The C++"
        << " find_mass_minus_one enhancement is what allows this case"
        << " to flow (literal-upstream Position[..., -1] would miss).";
    EXPECT_LT(child.family().n_loops(), fc.n_loops());
}

TEST(AmfsystemTest, CutkoskyEndingQ_PhaseVolumeComponentSelectsScheme) {
    auto fc = qft::FamilyConfig::build(
        "cutbubble", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2", "(l+p)^2"},
        /*cut=*/{1, 1});

    pipeline::AMFSystemOptions opts;
    std::vector<qft::JIntegral> preferred{qft::JIntegral("cutbubble", {1, 1})};

    EXPECT_TRUE(pipeline::ending_q(fc, preferred,
                                   pipeline::EndingScheme::Tradition, opts));
    EXPECT_FALSE(pipeline::ending_q(fc, preferred,
                                    pipeline::EndingScheme::Cutkosky, opts));
}

TEST(AmfsystemTest, CutkoskySetupMaster_ClearsCutAndPrescriptionForSystem) {
    auto fc = qft::FamilyConfig::build(
        "cutbubble", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2", "(l+p)^2"},
        /*cut=*/{1, 1});

    pipeline::AMFSystemOptions opts;
    opts.ending_schemes = {
        pipeline::EndingScheme::Tradition,
        pipeline::EndingScheme::Cutkosky,
        pipeline::EndingScheme::SingleMass
    };
    std::vector<qft::JIntegral> preferred{qft::JIntegral("cutbubble", {1, 1})};

    auto systems = pipeline::amf_system_setup_master(fc, preferred, opts);
    ASSERT_EQ(systems.size(), 1u);

    const auto& root = *systems[0];
    EXPECT_FALSE(root.is_ending());
    EXPECT_EQ(root.ending_scheme(), pipeline::EndingScheme::Cutkosky);
    ASSERT_EQ(root.family().cut.size(), 2u);
    EXPECT_EQ(root.family().cut[0], 0);
    EXPECT_EQ(root.family().cut[1], 0);
    ASSERT_EQ(root.family().prescription.size(), 1u);
    EXPECT_EQ(root.family().prescription[0], 1);
    EXPECT_EQ(std::count(root.etac().begin(), root.etac().end(), -1), 1);
}

// --- Cutkosky physical-mass safety check ---------
//
// Phase 1A added the upstream `AMFlow.m:1050` guard at
// `src/pipeline/amfsystem.cpp:2911-2953`: abort if a phase-volume
// component's mass is negative after `Numeric` substitution
// (negative squared masses are unphysical for a Cutkosky cut and
// would otherwise produce silently meaningless numbers).  No
// production-scale oracle bench triggers this branch because all 4
// committed Cutkosky benches use positive masses.
//
// This test constructs a Cutkosky family with `Numeric` set to a
// negative phase-volume-component mass and verifies the guard
// fires with the expected error message.

TEST(AmfsystemTest, CutkoskySetupMaster_NegativePhaseMassAbortsWithUpstreamMessage) {
    // 1-loop cutbubble with msq as the phase-volume mass.  Cutkosky
    // is the relevant scheme (Tradition would trigger the D3
    // tradition-with-cut projection, which is fine but tests a
    // different code path).
    auto fc = qft::FamilyConfig::build(
        "cutbubble_negm", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - msq", "(l + p)^2 - msq"},
        /*cut=*/{1, 1});

    pipeline::AMFSystemOptions opts;
    opts.ending_schemes = {pipeline::EndingScheme::Cutkosky};
    // Negative msq → unphysical squared mass on the phase-volume
    // component.  The guard at amfsystem.cpp:2944-2952 should fire.
    opts.bb.numeric_values["msq"] = "-1";
    opts.bb.numeric_values["s"]   = "100";
    std::vector<qft::JIntegral> preferred{
        qft::JIntegral("cutbubble_negm", {1, 1})};

    try {
        pipeline::amf_system_setup_master(fc, preferred, opts);
        FAIL() << "expected Cutkosky setup to abort on negative msq";
    } catch (const std::runtime_error& e) {
        const std::string what(e.what());
        EXPECT_NE(what.find("negative"), std::string::npos)
            << "error should explain the negative-mass cause; got: " << what;
        EXPECT_NE(what.find("AMFlow.m:1050"), std::string::npos)
            << "error should cite the upstream reference; got: " << what;
    }
}

TEST(AmfsystemTest, SingleMassSetup_TadpoleStrictlyShrinksToEndingSystem) {
    auto fc = qft::FamilyConfig::build(
        "tad", {"l"}, {}, {}, {}, {"l^2 - 1"});

    pipeline::AMFSystemOptions opts;
    opts.ending_schemes = {pipeline::EndingScheme::SingleMass};
    std::vector<qft::JIntegral> preferred{qft::JIntegral("tad", {2})};

    auto systems = pipeline::amf_system_setup_master(fc, preferred, opts);
    ASSERT_EQ(systems.size(), 1u);

    const auto& child = *systems[0];
    EXPECT_TRUE(child.is_ending());
    EXPECT_LT(child.family().n_loops(), fc.n_loops());
    EXPECT_LT(qft::get_top_position(child.preferred()).size(),
              qft::get_top_position(preferred).size());
}

TEST(AmfsystemTest, SingleMassSetup_TranslatesComponentPropagatorsByOwnContext) {
    auto fc = qft::FamilyConfig::build(
        "smctx", {"l1", "l2", "l3"}, {"n"}, {}, {{"n^2", "-1"}},
        {"l1^2 - 1", "l2^2", "l3^2", "(l1 + l2 + l3)^2"});

    pipeline::AMFSystemOptions opts;
    opts.ending_schemes = {pipeline::EndingScheme::SingleMass};
    std::vector<qft::JIntegral> preferred{
        qft::JIntegral("smctx", {1, 1, 1, 1})};

    auto systems = pipeline::amf_system_setup_master(fc, preferred, opts);
    ASSERT_EQ(systems.size(), 1u);

    const auto& child_fc = systems[0]->family();
    EXPECT_EQ(child_fc.n_loops(), 2u);
    EXPECT_EQ(child_fc.n_red_legs(), 1u);
    EXPECT_EQ(child_fc.ctx->var_index("n"), -1);

    bool saw_promoted_leg = false;
    for (const auto& prop : child_fc.propagators_after_conservation) {
        const std::string s = prop.to_string();
        EXPECT_EQ(s.find("n"), std::string::npos);
        if (s.find("l1") != std::string::npos) saw_promoted_leg = true;
    }
    EXPECT_TRUE(saw_promoted_leg);
}

TEST(AmfsystemTest, EndingTadpole_AtEps0p1) {
    auto fc = qft::FamilyConfig::build(
        "tad", {"l"}, {}, {}, {}, {"l^2 - msq"});

    pipeline::AMFSystemOptions opts;
    opts.bb.numeric_values["msq"] = "1";

    std::vector<qft::JIntegral> preferred;
    preferred.push_back(qft::JIntegral("tad", {1}));
    std::vector<int> zero_etac(1, 0);

    pipeline::AMFSystem sys(std::move(fc), preferred, zero_etac,
                            pipeline::EndingScheme::Tradition, opts);
    EXPECT_TRUE(sys.is_ending());
    sys.setup();

    numeric::AcbValue eps;
    fmpq_t q;
    fmpq_init(q);
    fmpq_set_si(q, 1, 10);
    acb_set_fmpq(eps.raw(), q, 200);
    fmpq_clear(q);

    std::vector<numeric::AcbValue> epslist;
    epslist.push_back(std::move(eps));

    sys.solve(epslist);

    ASSERT_EQ(sys.solutions().size(), 1u);
    ASSERT_EQ(sys.solutions()[0].master_values.size(), 1u);

    arf_t mid;
    arf_init(mid);
    arf_set(mid, arb_midref(acb_realref(
        sys.solutions()[0].master_values[0].raw())));
    double got = arf_get_d(mid, ARF_RND_NEAR);
    arf_clear(mid);

    EXPECT_NEAR(got, 10.5705641, 1e-6);
}

TEST(AmfsystemTest, BubblePipeline_RunsToCompletion) {
    if (!kira_available()) {
        GTEST_SKIP() << "Kira/Fermat not installed";
    }

    auto fc = qft::FamilyConfig::build(
        "bubblefam", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - msq", "(l - p)^2 - msq"});

    pipeline::AMFSystemOptions opts;
    opts.bb.kira_executable   = "/usr/local/bin/kira";
    opts.bb.fermat_executable = "/usr/share/Ferl7/fer64";
    opts.bb.numeric_values["s"]   = "5";
    opts.bb.numeric_values["msq"] = "1";
    opts.ending_schemes = {pipeline::EndingScheme::Tradition,
                           pipeline::EndingScheme::SingleMass};
    opts.cache_root = (fs::temp_directory_path() /
                         "amflow_v2_amfsys_bubble").string();
    opts.max_recursion_depth = 16;
    fs::remove_all(opts.cache_root);

    std::vector<qft::JIntegral> preferred;
    preferred.push_back(qft::JIntegral("bubblefam", {1, 1}));

    auto systems = pipeline::amf_systems_setup(fc, preferred, opts);
    ASSERT_FALSE(systems.empty());

    numeric::AcbValue eps;
    fmpq_t q;
    fmpq_init(q);
    fmpq_set_si(q, 1, 100);
    acb_set_fmpq(eps.raw(), q, 200);
    fmpq_clear(q);

    std::vector<numeric::AcbValue> epslist;
    epslist.push_back(std::move(eps));

    auto sols = pipeline::amf_systems_solution(systems, epslist);
    ASSERT_EQ(sols.size(), 1u);
    ASSERT_FALSE(sols[0].master_values.empty());
}

// --- GRAVITYFLOW_FIXEDPOINT guard: true fixed point refuses by name --------
//
// The self-similarity guard fires when a boundary child
// reproduces its parent's propagator multiset with NO top-sector
// descent — the genuinely self-similar b_{k+1} == b_k recursion
// (cut-eikonal PM chains) that would grind to max_recursion_depth
// without producing boundary data.  This test constructs exactly the
// recursion state build_boundary records for such a child: the parent's
// canonical propagator signature (sorted propagators-after-conservation
// joined with ';') equals the child's own, and the parent's top-sector
// width equals the child's (4 active slots — zero descent).  setup()
// must refuse by name with the structured GRAVITYFLOW_FIXEDPOINT
// message BEFORE any engine work (the guard sits ahead of eta injection
// and the Kira-backed diffeq/boundary stages, so no Kira/Fermat is
// needed here).
//
// The legitimate width-DESCENDING reproduction — same propagator
// multiset, strictly smaller child top sector, e.g. the massless box's
// automatic-loop boundary child whose single master is a 2-slot
// sub-sector bubble — must NOT fire; that side is pinned end-to-end by
// SolveIntegralsMmaRefTest.
// BlackBoxAMFlowSingle_Box1AutomaticLoopTargetsMatchMathReferenceJson,
// which runs WITHOUT AMFLOW_ALLOW_FIXEDPOINT_RECURSION.

TEST(AmfsystemTest, GravityflowFixedpointGuard_TrueFixedPointRefusesByName) {
    // Arm the guard regardless of ambient environment: it stands down
    // under the pmflow pipeline knobs and the debug override.
    unsetenv("AMFLOW_FORCE_ENDING_DEPTH");
    unsetenv("AMFLOW_FIXEDPOINT_PROBE");
    unsetenv("AMFLOW_ALLOW_FIXEDPOINT_RECURSION");

    // A boundary child at depth 1 whose family reproduces its parent —
    // same shape the misfiring box1 child had, but here the preferred
    // masters span the FULL 4-slot top sector, exactly the parent's
    // recorded width: no descent, a true fixed point.
    auto fc = qft::FamilyConfig::build(
        "truefp_b0_r0_f0",
        {"l"}, {"p1", "p2", "p3", "p4"},
        {{"p4", "-p1 - p2 - p3"}},
        {{"p1^2", "0"}, {"p2^2", "0"}, {"p3^2", "0"}, {"p4^2", "0"},
         {"(p1 + p2)^2", "s"}, {"(p1 + p3)^2", "t"}},
        {"l^2", "(l + p1)^2", "(l + p1 + p2)^2", "(l + p1 + p2 + p4)^2"});

    // The parent's recorded signature: the child's own canonical
    // propagator signature (sorted propagators-after-conservation,
    // ';'-joined with a trailing ';') — the self-similar reproduction.
    std::vector<std::string> ps;
    for (const auto& p : fc.propagators_after_conservation) {
        ps.push_back(p.to_string());
    }
    std::sort(ps.begin(), ps.end());
    std::string parent_sig;
    for (const auto& s : ps) {
        parent_sig += s;
        parent_sig += ';';
    }

    pipeline::AMFSystemOptions opts;
    opts.current_depth = 1;
    opts.parent_prop_signature = parent_sig;
    opts.parent_top_sector_count = 4;   // == child width: zero descent

    // Non-ending Tradition system (nonzero etac) so setup() reaches the
    // guard — the recursion continuation point endings never see.
    std::vector<qft::JIntegral> preferred{
        qft::JIntegral("truefp_b0_r0_f0", {1, 1, 1, 1})};
    std::vector<int> etac{-1, 0, 0, 0};

    pipeline::AMFSystem sys(std::move(fc), preferred, std::move(etac),
                            pipeline::EndingScheme::Tradition, opts);
    ASSERT_FALSE(sys.is_ending());

    try {
        sys.setup();
        FAIL() << "expected the GRAVITYFLOW_FIXEDPOINT guard to refuse "
                  "the zero-descent self-similar boundary child";
    } catch (const std::runtime_error& e) {
        const std::string what(e.what());
        EXPECT_NE(what.find("GRAVITYFLOW_FIXEDPOINT: boundary family "
                            "'truefp_b0_r0_f0' at depth 1"),
                  std::string::npos)
            << "refusal must carry the structured code + family + depth; "
               "got: " << what;
        EXPECT_NE(what.find("reproduces its parent's propagator set "
                            "(4 propagators identical after conservation)"),
                  std::string::npos)
            << "refusal must keep the existing message; got: " << what;
        EXPECT_NE(what.find("AMFLOW_FIXEDPOINT_PROBE"), std::string::npos)
            << "refusal must point at the pmflow closure route; got: "
            << what;
    }
}
