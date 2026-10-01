// SPDX-License-Identifier: MIT
// Tests for amflow::qft::amfmode.

#include <algorithm>
#include <gtest/gtest.h>

#include "amflow/qft/amfmode.hpp"
#include "amflow/qft/family_config.hpp"

namespace qft = amflow::qft;

TEST(AMFModeTest, ModeNames_RoundTrip) {
    for (auto m : {qft::AMFMode::Prescription, qft::AMFMode::Mass,
                   qft::AMFMode::Propagator,   qft::AMFMode::Branch,
                   qft::AMFMode::Loop,         qft::AMFMode::All}) {
        EXPECT_EQ(qft::parse_amf_mode(qft::amf_mode_name(m)), m);
    }
}

TEST(AMFModeTest, ParseMode_ThrowsOnUnknown) {
    EXPECT_THROW(qft::parse_amf_mode("Bogus"), std::invalid_argument);
}

TEST(AMFModeTest, FamilyConfig_DefaultsForCutAndPrescription) {
    auto fc = qft::FamilyConfig::build(
        "bubble", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - msq", "(l - p)^2 - msq"});
    ASSERT_EQ(fc.cut.size(), 2u);
    EXPECT_EQ(fc.cut[0], 0);
    EXPECT_EQ(fc.cut[1], 0);
    ASSERT_EQ(fc.prescription.size(), 1u);
    EXPECT_EQ(fc.prescription[0], 1);
}

TEST(AMFModeTest, FamilyConfig_CutLengthMismatchThrows) {
    EXPECT_THROW(qft::FamilyConfig::build(
        "bubble", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - msq", "(l - p)^2 - msq"}, {0}),
                 std::invalid_argument);
}

TEST(AMFModeTest, FamilyConfig_PrescriptionLengthMismatchThrows) {
    EXPECT_THROW(qft::FamilyConfig::build(
        "bubble", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - msq", "(l - p)^2 - msq"}, {}, {1, 1}),
                 std::invalid_argument);
}

TEST(AMFModeTest, PrescriptionOfProp_AllPositive) {
    auto fc = qft::FamilyConfig::build(
        "sunrise", {"l1", "l2"}, {"p"}, {}, {{"p^2", "s"}},
        {"l1^2 - m1sq", "l2^2 - m2sq", "(l1 + l2 - p)^2 - m3sq"});
    auto p = fc.prescription_of_prop(fc.propagators_after_conservation[0]);
    ASSERT_TRUE(p.has_value());
    EXPECT_EQ(p.value(), 1);
}

TEST(AMFModeTest, PrescriptionOfProp_MixedFails) {
    auto fc = qft::FamilyConfig::build(
        "sunrise", {"l1", "l2"}, {"p"}, {}, {{"p^2", "s"}},
        {"l1^2 - m1sq", "l2^2 - m2sq", "(l1 + l2 - p)^2 - m3sq"},
        {}, {1, -1});
    auto p = fc.prescription_of_prop(fc.propagators_after_conservation[2]);
    EXPECT_FALSE(p.has_value());
}

TEST(AMFModeTest, AnalyzeTopSector_Sunrise_FullSector) {
    auto fc = qft::FamilyConfig::build(
        "sunrise", {"l1", "l2"}, {"p"}, {}, {{"p^2", "s"}},
        {"l1^2 - m1sq", "l2^2 - m2sq", "(l1 + l2 - p)^2 - m3sq"});

    auto info = qft::analyze_top_sector(fc, {0, 1, 2});
    ASSERT_EQ(info.size(), 1u);
    EXPECT_EQ(info[0].loopnum, 2);
    ASSERT_EQ(info[0].prop_index.size(), 3u);
    std::vector<std::size_t> pp = info[0].prop_index;
    std::sort(pp.begin(), pp.end());
    EXPECT_EQ(pp, (std::vector<std::size_t>{0, 1, 2}));
    for (int c : info[0].cut) EXPECT_EQ(c, 0);
    for (int p : info[0].pres) EXPECT_EQ(p, 1);
}

TEST(AMFModeTest, AnalyzeTopSector_DoubleBubble_TwoComponents) {
    auto fc = qft::FamilyConfig::build(
        "doublebubble", {"l1", "l2"}, {"p"}, {}, {{"p^2", "s"}},
        {"l1^2 - m1sq", "(l1 - p)^2 - m1sq",
         "l2^2 - m2sq", "(l2 - p)^2 - m2sq"});

    auto info = qft::analyze_top_sector(fc, {0, 1, 2, 3});
    ASSERT_EQ(info.size(), 2u);
    EXPECT_EQ(info[0].loopnum, 1);
    EXPECT_EQ(info[1].loopnum, 1);
    EXPECT_EQ(info[0].var.size(), 2u);
    EXPECT_EQ(info[1].var.size(), 2u);
}

TEST(AMFModeTest, AnalyzeTopSector_AbortsOnMixedPrescription) {
    auto fc = qft::FamilyConfig::build(
        "mixed", {"l1", "l2"}, {"p"}, {}, {{"p^2", "s"}},
        {"(l1 + l2 - p)^2"}, {}, {1, -1});
    EXPECT_THROW(qft::analyze_top_sector(fc, {0}), std::runtime_error);
}

TEST(AMFModeTest, AllPossiblePosition_Propagator_OneLoopBubble) {
    auto fc = qft::FamilyConfig::build(
        "bubble", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - msq", "(l - p)^2 - msq"});
    auto info = qft::analyze_top_sector(fc, {0, 1});

    auto groups = qft::all_possible_position(fc, info[0], qft::AMFMode::Propagator);
    ASSERT_EQ(groups.size(), 2u);
    for (const auto& g : groups) ASSERT_EQ(g.size(), 1u);
    std::vector<std::size_t> all;
    for (const auto& g : groups) all.push_back(g[0]);
    std::sort(all.begin(), all.end());
    EXPECT_EQ(all, (std::vector<std::size_t>{0, 1}));
}

TEST(AMFModeTest, AllPossiblePosition_Propagator_RespectsCut) {
    auto fc = qft::FamilyConfig::build(
        "bubble_cut", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - msq", "(l - p)^2 - msq"}, {1, 0});
    auto info = qft::analyze_top_sector(fc, {0, 1});

    auto groups = qft::all_possible_position(fc, info[0], qft::AMFMode::Propagator);
    ASSERT_EQ(groups.size(), 1u);
    EXPECT_EQ(groups[0][0], 1u);
}

TEST(AMFModeTest, AllPossiblePosition_Branch_Sunrise) {
    auto fc = qft::FamilyConfig::build(
        "sunrise", {"l1", "l2"}, {"p"}, {}, {{"p^2", "s"}},
        {"l1^2 - m1sq", "l2^2 - m2sq", "(l1 + l2 - p)^2 - m3sq"});
    auto info = qft::analyze_top_sector(fc, {0, 1, 2});

    auto groups = qft::all_possible_position(fc, info[0], qft::AMFMode::Branch);
    EXPECT_EQ(groups.size(), 3u);
    for (auto& g : groups) EXPECT_EQ(g.size(), 1u);
}

TEST(AMFModeTest, AllPossiblePosition_Loop_Sunrise) {
    auto fc = qft::FamilyConfig::build(
        "sunrise", {"l1", "l2"}, {"p"}, {}, {{"p^2", "s"}},
        {"l1^2 - m1sq", "l2^2 - m2sq", "(l1 + l2 - p)^2 - m3sq"});
    auto info = qft::analyze_top_sector(fc, {0, 1, 2});

    auto groups = qft::all_possible_position(fc, info[0], qft::AMFMode::Loop);
    EXPECT_EQ(groups.size(), 3u);
    for (auto& g : groups) EXPECT_EQ(g.size(), 2u);
}

TEST(AMFModeTest, AllPossiblePosition_Mass_Sunrise_DistinctMasses) {
    auto fc = qft::FamilyConfig::build(
        "sunrise", {"l1", "l2"}, {"p"}, {}, {{"p^2", "s"}},
        {"l1^2 - m1sq", "l2^2 - m2sq", "(l1 + l2 - p)^2 - m3sq"});
    auto info = qft::analyze_top_sector(fc, {0, 1, 2});

    auto groups = qft::all_possible_position(fc, info[0], qft::AMFMode::Mass);
    EXPECT_EQ(groups.size(), 3u);
    for (auto& g : groups) EXPECT_EQ(g.size(), 1u);
}

TEST(AMFModeTest, AllPossiblePosition_Mass_Sunrise_SharedMass) {
    auto fc = qft::FamilyConfig::build(
        "sunrise_shared", {"l1", "l2"}, {"p"}, {}, {{"p^2", "s"}},
        {"l1^2 - msq", "l2^2 - msq", "(l1 + l2 - p)^2 - msq"});
    auto info = qft::analyze_top_sector(fc, {0, 1, 2});

    auto groups = qft::all_possible_position(fc, info[0], qft::AMFMode::Mass);
    EXPECT_EQ(groups.size(), 1u);
    EXPECT_EQ(groups[0].size(), 3u);
}

TEST(AMFModeTest, AllPossiblePosition_Mass_PreservesGatherByOrder) {
    auto fc = qft::FamilyConfig::build(
        "boundary_test", {"l1", "l2", "l3"}, {"n"}, {}, {{"n^2", "-1"}},
        {"l1^2", "l2^2", "l3^2",
         "(l1 + 20*n)^2 + 440",
         "(l1 + l2 + l3 + 20*n)^2 + 420",
         "(l1 + l2)^2", "(l1 + l3)^2", "(l2 + l3)^2", "(l2 + n)^2"});
    auto info = qft::analyze_top_sector(fc, {0, 1, 2, 3, 4});
    ASSERT_FALSE(info.empty());

    auto groups = qft::all_possible_position(fc, info[0], qft::AMFMode::Mass);
    ASSERT_GE(groups.size(), 2u);
    EXPECT_EQ(groups[0], (std::vector<std::size_t>{3}));
    EXPECT_EQ(groups[1], (std::vector<std::size_t>{4}));

    auto pos = qft::amf_position(fc, {0, 1, 2, 3, 4},
                                  {qft::AMFMode::Prescription,
                                   qft::AMFMode::Mass,
                                   qft::AMFMode::Propagator});
    EXPECT_EQ(pos, (std::vector<std::size_t>{3}));
}

TEST(AMFModeTest, AllPossiblePosition_Prescription_HasNegativeLoop) {
    auto fc = qft::FamilyConfig::build(
        "presc", {"l1", "l2"}, {"p"}, {}, {{"p^2", "s"}},
        {"l1^2 - m1sq", "l2^2 - m2sq", "l1*l2 + 1"},
        {}, {1, -1});
    auto info = qft::analyze_top_sector(fc, {0, 1});
    ASSERT_FALSE(info.empty());

    bool found_neg_one = false;
    for (auto& ci : info) {
        auto groups = qft::all_possible_position(fc, ci, qft::AMFMode::Prescription);
        for (auto& g : groups) {
            for (std::size_t v : g) {
                if (v == 1) found_neg_one = true;
            }
        }
    }
    EXPECT_TRUE(found_neg_one);
}

TEST(AMFModeTest, AllPossiblePosition_All_Sunrise) {
    auto fc = qft::FamilyConfig::build(
        "sunrise", {"l1", "l2"}, {"p"}, {}, {{"p^2", "s"}},
        {"l1^2 - m1sq", "l2^2 - m2sq", "(l1 + l2 - p)^2 - m3sq"});
    auto info = qft::analyze_top_sector(fc, {0, 1, 2});

    auto groups = qft::all_possible_position(fc, info[0], qft::AMFMode::All);
    ASSERT_EQ(groups.size(), 1u);
    EXPECT_EQ(groups[0].size(), 3u);
}

TEST(AMFModeTest, AMFPosition_Sunrise_PropagatorMode) {
    auto fc = qft::FamilyConfig::build(
        "sunrise", {"l1", "l2"}, {"p"}, {}, {{"p^2", "s"}},
        {"l1^2 - m1sq", "l2^2 - m2sq", "(l1 + l2 - p)^2 - m3sq"});
    auto pos = qft::amf_position(fc, {0, 1, 2}, qft::AMFMode::Propagator);
    ASSERT_FALSE(pos.empty());
    EXPECT_EQ(pos.size(), 1u);
}

TEST(AMFModeTest, AMFPosition_ModeListFallback) {
    auto fc = qft::FamilyConfig::build(
        "sunrise", {"l1", "l2"}, {"p"}, {}, {{"p^2", "s"}},
        {"l1^2 - m1sq", "l2^2 - m2sq", "(l1 + l2 - p)^2 - m3sq"});

    auto pos = qft::amf_position(fc, {0, 1, 2},
                                  {qft::AMFMode::Prescription, qft::AMFMode::Mass});
    ASSERT_FALSE(pos.empty());
}

TEST(AMFModeTest, AMFPosition_EmptyTopposi) {
    auto fc = qft::FamilyConfig::build(
        "bubble", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - msq", "(l - p)^2 - msq"});
    auto pos = qft::amf_position(fc, {}, qft::AMFMode::Propagator);
    EXPECT_TRUE(pos.empty());
}

TEST(AMFModeTest, AMFEtaC_BasicNegativeOne) {
    auto fc = qft::FamilyConfig::build(
        "sunrise", {"l1", "l2"}, {"p"}, {}, {{"p^2", "s"}},
        {"l1^2 - m1sq", "l2^2 - m2sq", "(l1 + l2 - p)^2 - m3sq"});

    auto eta = qft::amf_eta_c(fc, {0, 2});
    ASSERT_EQ(eta.size(), 3u);
    EXPECT_EQ(eta[0], -1);
    EXPECT_EQ(eta[1], 0);
    EXPECT_EQ(eta[2], -1);
}

TEST(AMFModeTest, AMFEtaC_OutOfRangeThrows) {
    auto fc = qft::FamilyConfig::build(
        "bubble", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - msq", "(l - p)^2 - msq"});
    EXPECT_THROW(qft::amf_eta_c(fc, {0, 5}), std::out_of_range);
}

TEST(AMFModeTest, VacuumQ_SunriseIsNot) {
    auto fc = qft::FamilyConfig::build(
        "sunrise", {"l1", "l2"}, {"p"}, {}, {{"p^2", "s"}},
        {"l1^2 - m1sq", "l2^2 - m2sq", "(l1 + l2 - p)^2 - m3sq"});
    auto info = qft::analyze_top_sector(fc, {0, 1, 2});
    ASSERT_FALSE(info.empty());
    EXPECT_FALSE(info[0].vacQ);
    EXPECT_FALSE(qft::vacuum_q(info[0]));
    EXPECT_FALSE(qft::ending_q(info[0]));
}

TEST(AMFModeTest, EndingQ_OnlyTadpoleSector) {
    auto fc = qft::FamilyConfig::build(
        "subsector", {"l"}, {}, {}, {}, {"l^2 - 1"});
    auto info = qft::analyze_top_sector(fc, {0});
    ASSERT_FALSE(info.empty());
    EXPECT_TRUE(info[0].vacQ);
    EXPECT_TRUE(qft::vacuum_q(info[0]));
    EXPECT_TRUE(qft::single_mass_q(info[0]));
    EXPECT_TRUE(qft::ending_q(info[0]));
}

// --- SingleMassQ literal-vs-symbolic divergence -
//
// Lock the contract that `qft::single_mass_q` is **upstream-literal**:
// it tests `info.mass[k]` for being literal `0` or literal `1`
// without applying any substitution (mirrors upstream's
// `Count[mass, 1] === 1 && Count[mass, 0] === Length-1`,
// AMFlow.m).  A free-symbol mass like `msq` in the propagator does
// NOT count as a "mass = 1" here even if the user later supplies
// `Numeric = {msq -> 1}` — the pipeline-level wrapper
// `single_mass_q_numeric` (`src/pipeline/amfsystem.cpp:803`) is the
// C++ enhancement that substitutes Numeric before testing, see
// matching pipeline test
// `AmfsystemTest.SingleMassEnhancement_*`.
TEST(AMFModeTest, SingleMassQ_FreeSymbolMassReturnsFalse_LiteralUpstreamSemantic) {
    // Tadpole with a free-symbol mass `msq`, NO Numeric supplied.
    // Mass list literally is `[msq]`, which is neither 1 nor 0.
    // Upstream-faithful literal SingleMassQ → false.
    auto fc = qft::FamilyConfig::build(
        "tadpole_sym", {"l"}, {}, {}, {}, {"l^2 - msq"});
    auto info = qft::analyze_top_sector(fc, {0});
    ASSERT_FALSE(info.empty());
    EXPECT_TRUE(qft::vacuum_q(info[0]));   // cut all zero
    EXPECT_FALSE(qft::single_mass_q(info[0]))
        << "qft::single_mass_q must test literally; symbolic mass `msq`"
        << " is not literal `1`, so the component is vacuum-but-not-single-mass.";
    EXPECT_FALSE(qft::ending_q(info[0]));
}

TEST(AMFModeTest, SingleMassQ_MixedLiteralAndSymbolicMassReturnsFalse) {
    // Bubble with one massless and one symbolic-mass propagator:
    //   p_0 = l^2          (mass = 0 literal)
    //   p_1 = (l-p)^2 - msq (mass = msq symbolic)
    // qft::single_mass_q tests literally: count(==1) = 0 ≠ 1 → false.
    auto fc = qft::FamilyConfig::build(
        "bub_mixed_sym", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2", "(l - p)^2 - msq"});
    auto info = qft::analyze_top_sector(fc, {0, 1});
    ASSERT_FALSE(info.empty());
    EXPECT_FALSE(qft::single_mass_q(info[0]))
        << "Symbolic msq in mass list must not be detected as the literal `1`";
}

TEST(AMFModeTest, AMFCandidateComponent_GeneralVacuumFallsBackToBranch) {
    auto fc = qft::FamilyConfig::build(
        "vacuum_branch", {"l"}, {}, {}, {},
        {"l^2 - 1", "l^2 - 2"});

    auto info = qft::analyze_top_sector(fc, {0, 1});
    ASSERT_EQ(info.size(), 1u);
    EXPECT_TRUE(qft::vacuum_q(info[0]));
    EXPECT_FALSE(qft::single_mass_q(info[0]));

    auto cand = qft::amf_candidate_component(fc, info[0], qft::AMFMode::Mass);
    EXPECT_EQ(cand, (std::vector<std::size_t>{0, 1}));

    auto pos = qft::amf_position(fc, {0, 1}, qft::AMFMode::Mass);
    EXPECT_EQ(pos, (std::vector<std::size_t>{0, 1}));
}
