// SPDX-License-Identifier: MIT
// Tests for amflow::qft::jintegral.

#include <gtest/gtest.h>

#include "amflow/qft/jintegral.hpp"
#include "amflow/numeric/options.hpp"

#include <stdexcept>
#include <vector>

namespace nm  = amflow::numeric;
namespace qft = amflow::qft;

namespace {

class JIntegralTest : public ::testing::Test {
protected:
    void SetUp() override    { nm::set_default_options(); }
    void TearDown() override { nm::set_default_options(); }
};

qft::JIntegral J(const std::string& fam, std::vector<long> idx) {
    return qft::JIntegral(fam, std::move(idx));
}

}  // namespace

// =========================================================================
//  JIntegral basics
// =========================================================================

TEST_F(JIntegralTest, Construction_RejectsEmptyFamily) {
    EXPECT_THROW(J("", {1, 0, 1}), std::invalid_argument);
}

TEST_F(JIntegralTest, Indices_NPropsDotsRank) {
    auto j = J("box", {2, 0, 1, -3, 1});
    EXPECT_EQ(j.n_props(), 3);
    EXPECT_EQ(j.n_dots(),  1);
    EXPECT_EQ(j.rank(),    3);
}

TEST_F(JIntegralTest, Sector_PatternAndIntegral) {
    auto j = J("box", {2, 0, 1, -3, 1});
    EXPECT_EQ(j.sector_pattern(), (std::vector<int>{1, 0, 1, 0, 1}));

    auto sec = j.sector();
    EXPECT_EQ(sec.family(), "box");
    EXPECT_EQ(sec.indices(), (std::vector<long>{1, 0, 1, 0, 1}));
}

TEST_F(JIntegralTest, ToString_HumanReadable) {
    auto j = J("box1", {1, 0, -2});
    EXPECT_EQ(j.to_string(), "j[box1, 1, 0, -2]");
}

TEST_F(JIntegralTest, Equality_FamilyAndIndices) {
    EXPECT_EQ(J("a", {1, 2}), J("a", {1, 2}));
    EXPECT_NE(J("a", {1, 2}), J("a", {1, 3}));
    EXPECT_NE(J("a", {1, 2}), J("b", {1, 2}));
}

// =========================================================================
//  Sort
// =========================================================================

TEST_F(JIntegralTest, Sort_ByDecreasingProps) {
    auto a = J("f", {1, 1, 1, 0});
    auto b = J("f", {1, 1, 0, 0});
    auto sorted = qft::sort_integrals({b, a});
    ASSERT_EQ(sorted.size(), 2u);
    EXPECT_EQ(sorted[0], a);
    EXPECT_EQ(sorted[1], b);
}

TEST_F(JIntegralTest, Sort_TiebreakByDecreasingDots) {
    auto a = J("f", {1, 1});
    auto b = J("f", {2, 1});
    auto sorted = qft::sort_integrals({a, b});
    EXPECT_EQ(sorted[0], b);
    EXPECT_EQ(sorted[1], a);
}

TEST_F(JIntegralTest, Sort_TiebreakByDecreasingRank) {
    auto a = J("f", {1, 1,  0});
    auto b = J("f", {1, 1, -3});
    auto sorted = qft::sort_integrals({a, b});
    EXPECT_EQ(sorted[0], b);
    EXPECT_EQ(sorted[1], a);
}

TEST_F(JIntegralTest, Sort_TiebreakByFamilyName) {
    auto a = J("alpha", {1, 1, 0});
    auto b = J("beta",  {1, 1, 0});
    auto sorted = qft::sort_integrals({b, a});
    EXPECT_EQ(sorted[0], a);
    EXPECT_EQ(sorted[1], b);
}

TEST_F(JIntegralTest, SortWeight_FieldOrder) {
    auto j = J("f", {2, 0, 1, -3, 1});
    auto w = j.sort_weight();
    EXPECT_EQ(std::get<0>(w), -3);
    EXPECT_EQ(std::get<1>(w), -1);
    EXPECT_EQ(std::get<2>(w), -3);
    EXPECT_EQ(std::get<3>(w), -3);
}

// =========================================================================
//  Sector containment
// =========================================================================

TEST_F(JIntegralTest, IsSubsector_BasicCases) {
    EXPECT_TRUE (qft::is_subsector({1, 0, 0}, {1, 1, 0}));
    EXPECT_TRUE (qft::is_subsector({0, 0, 0}, {1, 1, 1}));
    EXPECT_TRUE (qft::is_subsector({1, 1, 0}, {1, 1, 0}));
    EXPECT_FALSE(qft::is_subsector({1, 1, 1}, {1, 1, 0}));
    EXPECT_FALSE(qft::is_subsector({1, 0},    {1, 0, 0}));
}

TEST_F(JIntegralTest, IsTopSector_OnlyMaximal) {
    std::vector<std::vector<int>> sectors = {
        {1, 0, 0},
        {1, 1, 0},
        {1, 1, 1},
        {0, 1, 1}
    };
    EXPECT_TRUE (qft::is_top_sector({1, 1, 1}, sectors));
    EXPECT_FALSE(qft::is_top_sector({1, 1, 0}, sectors));
    EXPECT_FALSE(qft::is_top_sector({1, 0, 0}, sectors));
    EXPECT_FALSE(qft::is_top_sector({0, 1, 1}, sectors));
}

TEST_F(JIntegralTest, IsTopSector_TwoIncomparableTops) {
    std::vector<std::vector<int>> sectors = {
        {1, 1, 0},
        {0, 1, 1},
        {1, 0, 0}
    };
    EXPECT_TRUE (qft::is_top_sector({1, 1, 0}, sectors));
    EXPECT_TRUE (qft::is_top_sector({0, 1, 1}, sectors));
    EXPECT_FALSE(qft::is_top_sector({1, 0, 0}, sectors));
}

// =========================================================================
//  GetTopSector / GetTopPosition
// =========================================================================

TEST_F(JIntegralTest, GetTopSector_PerColumnMaxAboveZero) {
    auto a = J("f", { 1, 0, 1});
    auto b = J("f", {-1, 0, 1});
    auto top = qft::get_top_sector({a, b});
    EXPECT_EQ(top, (std::vector<int>{1, 0, 1}));

    auto pos = qft::get_top_position({a, b});
    EXPECT_EQ(pos, (std::vector<std::size_t>{0, 2}));
}

TEST_F(JIntegralTest, GetTopSector_EmptyInput) {
    auto top = qft::get_top_sector({});
    EXPECT_TRUE(top.empty());
}

TEST_F(JIntegralTest, GetTopSector_MismatchedSizesThrow) {
    auto a = J("f", {1, 0});
    auto b = J("f", {1, 0, 1});
    EXPECT_THROW(qft::get_top_sector({a, b}), std::invalid_argument);
}

// =========================================================================
//  GetTopSectorList / SplitTarget
// =========================================================================

TEST_F(JIntegralTest, GetTopSectorList_FiltersToMaximal) {
    auto j1 = J("f", {1, 1, 0});
    auto j2 = J("f", {1, 1, 1});
    auto j3 = J("f", {1, 0, 0});
    auto tops = qft::get_top_sector_list({j1, j2, j3});
    ASSERT_EQ(tops.size(), 1u);
    EXPECT_EQ(tops[0], (std::vector<int>{1, 1, 1}));
}

TEST_F(JIntegralTest, GetTopSectorList_TwoIncomparable) {
    auto j1 = J("f", {1, 1, 0});
    auto j2 = J("f", {0, 1, 1});
    auto tops = qft::get_top_sector_list({j1, j2});
    EXPECT_EQ(tops.size(), 2u);
}

TEST_F(JIntegralTest, SplitTarget_PartitionsByTopSector) {
    auto a = J("f", {1, 1, 0, 0});
    auto b = J("f", {1, 0, 0, 0});
    auto c = J("f", {0, 0, 1, 1});
    auto d = J("f", {0, 0, 1, 0});

    auto split = qft::split_target({a, b, c, d});
    ASSERT_EQ(split.size(), 2u);

    bool found_chunk_a = false, found_chunk_b = false;
    for (const auto& chunk : split) {
        bool has_a = false, has_b = false, has_c = false, has_d = false;
        for (const auto& j : chunk) {
            if (j == a) has_a = true;
            if (j == b) has_b = true;
            if (j == c) has_c = true;
            if (j == d) has_d = true;
        }
        if (has_a && has_b && !has_c && !has_d) found_chunk_a = true;
        if (!has_a && !has_b && has_c && has_d) found_chunk_b = true;
    }
    EXPECT_TRUE(found_chunk_a);
    EXPECT_TRUE(found_chunk_b);
}

TEST_F(JIntegralTest, SplitTarget_AmbiguousSubsectorGoesToFirst) {
    auto a = J("f", {1, 1, 0});
    auto b = J("f", {1, 0, 0});
    auto c = J("f", {0, 1, 1});
    auto d = J("f", {0, 1, 0});

    auto split = qft::split_target({a, b, c, d});
    ASSERT_EQ(split.size(), 2u);
    long total = 0;
    for (const auto& chunk : split) total += static_cast<long>(chunk.size());
    EXPECT_EQ(total, 4);

    long count_d = 0;
    for (const auto& chunk : split) {
        for (const auto& j : chunk) {
            if (j == d) ++count_d;
        }
    }
    EXPECT_EQ(count_d, 1);
}
