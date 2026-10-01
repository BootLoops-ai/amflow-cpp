// SPDX-License-Identifier: MIT
// Tests for amflow::ode sparse linear algebra (Layer 3).

#include <gtest/gtest.h>

#include <flint/acb.h>

#include "amflow/numeric/acb_wrap.hpp"
#include "amflow/numeric/options.hpp"
#include "amflow/ode/sparse.hpp"

#include <stdexcept>
#include <utility>
#include <vector>

namespace ode = amflow::ode;
namespace dsn = amflow::numeric;

namespace {

dsn::AcbValue acb_si(long v) {
    dsn::AcbValue out;
    out.set_si(v);
    return out;
}

bool acb_eq_si(acb_srcptr x, long v) {
    dsn::AcbValue tmp;
    tmp.set_si(v);
    return acb_eq(x, tmp.raw()) != 0;
}

// Build a coefficient list: [c0, c1, c2, ...]  for poly  c0 + c1*eta + c2*eta^2 + ...
std::vector<dsn::AcbValue> coeff_list(std::initializer_list<long> coeffs) {
    std::vector<dsn::AcbValue> out;
    out.reserve(coeffs.size());
    for (long c : coeffs) out.push_back(acb_si(c));
    return out;
}

}  // namespace

// ---------------------------------------------------------------------------
//  SparseEntry / Row / sparsify
// ---------------------------------------------------------------------------

TEST(SparsifyRow, DropsZeros) {
    std::vector<dsn::AcbValue> dense;
    dense.emplace_back();              // 0
    dense.push_back(acb_si(7));        // 7
    dense.emplace_back();              // 0
    dense.push_back(acb_si(3));        // 3

    auto row = ode::sparsify_row(dense);
    ASSERT_EQ(row.size(), 2u);
    EXPECT_EQ(row[0].col, 1);
    EXPECT_TRUE(acb_eq_si(row[0].val.raw(), 7));
    EXPECT_EQ(row[1].col, 3);
    EXPECT_TRUE(acb_eq_si(row[1].val.raw(), 3));
}

TEST(SparsifyMatrix, BuildsPerRow) {
    std::vector<std::vector<dsn::AcbValue>> dense(2);
    dense[0].push_back(acb_si(1));
    dense[0].emplace_back();
    dense[1].emplace_back();
    dense[1].push_back(acb_si(2));

    auto sys = ode::sparsify_matrix(dense);
    ASSERT_EQ(sys.size(), 2u);
    ASSERT_EQ(sys[0].size(), 1u);
    EXPECT_EQ(sys[0][0].col, 0);
    ASSERT_EQ(sys[1].size(), 1u);
    EXPECT_EQ(sys[1][0].col, 1);
}

TEST(SparseEntry, CopyCloneIsDeep) {
    ode::SparseEntry a(2, acb_si(5));
    ode::SparseEntry b = a;            // copy ctor → deep clone via val.clone()
    acb_set_si(a.val.raw(), 99);       // mutate a
    EXPECT_TRUE(acb_eq_si(b.val.raw(), 5));
    EXPECT_EQ(b.col, 2);
}

// ---------------------------------------------------------------------------
//  sparse_translate
// ---------------------------------------------------------------------------

TEST(SparseTranslate, ShiftsAllColumns) {
    ode::SparseRow row;
    row.emplace_back(0, acb_si(1));
    row.emplace_back(2, acb_si(3));

    ode::sparse_translate(row, 5);
    EXPECT_EQ(row[0].col, 5);
    EXPECT_EQ(row[1].col, 7);
}

TEST(SparseTranslate, OperatesOverWholeSystem) {
    ode::SparseSystem sys(2);
    sys[0].emplace_back(0, acb_si(1));
    sys[1].emplace_back(3, acb_si(2));
    ode::sparse_translate(sys, 10);
    EXPECT_EQ(sys[0][0].col, 10);
    EXPECT_EQ(sys[1][0].col, 13);
}

// ---------------------------------------------------------------------------
//  scalar_sparse
// ---------------------------------------------------------------------------

TEST(ScalarSparse, MultipliesEveryEntry) {
    ode::SparseRow row;
    row.emplace_back(0, acb_si(2));
    row.emplace_back(2, acb_si(3));

    auto two = acb_si(2);
    ode::scalar_sparse(row, two.raw());
    EXPECT_TRUE(acb_eq_si(row[0].val.raw(), 4));
    EXPECT_TRUE(acb_eq_si(row[1].val.raw(), 6));
}

// ---------------------------------------------------------------------------
//  plus_sparse
// ---------------------------------------------------------------------------

TEST(PlusSparse, MergesDistinctColumns) {
    ode::SparseRow a;
    a.emplace_back(0, acb_si(1));
    a.emplace_back(3, acb_si(4));
    ode::SparseRow b;
    b.emplace_back(1, acb_si(2));
    b.emplace_back(4, acb_si(5));

    auto sum = ode::plus_sparse(a, b);
    ASSERT_EQ(sum.size(), 4u);
    EXPECT_EQ(sum[0].col, 0);
    EXPECT_TRUE(acb_eq_si(sum[0].val.raw(), 1));
    EXPECT_EQ(sum[1].col, 1);
    EXPECT_EQ(sum[2].col, 3);
    EXPECT_EQ(sum[3].col, 4);
}

TEST(PlusSparse, AddsMatchingColumns) {
    ode::SparseRow a;
    a.emplace_back(2, acb_si(3));
    ode::SparseRow b;
    b.emplace_back(2, acb_si(4));
    auto sum = ode::plus_sparse(a, b);
    ASSERT_EQ(sum.size(), 1u);
    EXPECT_EQ(sum[0].col, 2);
    EXPECT_TRUE(acb_eq_si(sum[0].val.raw(), 7));
}

TEST(PlusSparse, CancelsToZeroDropsEntry) {
    ode::SparseRow a;
    a.emplace_back(2, acb_si(5));
    ode::SparseRow b;
    b.emplace_back(2, acb_si(-5));
    auto sum = ode::plus_sparse(a, b);
    EXPECT_TRUE(sum.empty());
}

// ---------------------------------------------------------------------------
//  split_system
// ---------------------------------------------------------------------------

TEST(SplitSystem, FindsRowsContainingColumn) {
    ode::SparseSystem sys(3);
    sys[0].emplace_back(1, acb_si(1));      // contains col 1
    sys[0].emplace_back(3, acb_si(2));
    sys[1].emplace_back(0, acb_si(1));      // doesn't contain col 1
    sys[1].emplace_back(2, acb_si(1));
    sys[2].emplace_back(1, acb_si(7));      // contains col 1

    auto split = ode::split_system(sys, 1);
    EXPECT_EQ(split.dep,   (std::vector<std::size_t>{0, 2}));
    EXPECT_EQ(split.indep, (std::vector<std::size_t>{1}));
}

TEST(SplitSystem, EarlyOutOnSortedRows) {
    // Verify the sorted-ascending early-exit doesn't miss a hit at the start.
    ode::SparseSystem sys(1);
    sys[0].emplace_back(0, acb_si(1));
    sys[0].emplace_back(5, acb_si(1));
    auto split = ode::split_system(sys, 0);
    EXPECT_EQ(split.dep, (std::vector<std::size_t>{0}));
}

// ---------------------------------------------------------------------------
//  construct_matrix
// ---------------------------------------------------------------------------

TEST(ConstructMatrix, EmptyAxThrows) {
    std::vector<std::vector<std::vector<dsn::AcbValue>>> ax;
    auto dx = coeff_list({1});
    EXPECT_THROW(ode::construct_matrix(dx, ax, 1), std::invalid_argument);
}

TEST(ConstructMatrix, NonSquareAxThrows) {
    std::vector<std::vector<std::vector<dsn::AcbValue>>> ax(1);
    ax[0].resize(2);
    ax[0][0] = coeff_list({1});
    ax[0][1] = coeff_list({1});
    auto dx = coeff_list({1});
    EXPECT_THROW(ode::construct_matrix(dx, ax, 1), std::invalid_argument);
}

TEST(ConstructMatrix, ZeroDxThrows) {
    std::vector<std::vector<std::vector<dsn::AcbValue>>> ax(1);
    ax[0].resize(1);
    ax[0][0] = coeff_list({1});
    std::vector<dsn::AcbValue> dx;       // empty == zero
    EXPECT_THROW(ode::construct_matrix(dx, ax, 1), std::invalid_argument);
}

TEST(ConstructMatrix, ShapeAndColumnsForTrivialSystem) {
    // dx = 1, ax = [[1]] (constant), total_order = 2.
    auto dx = coeff_list({1});
    std::vector<std::vector<std::vector<dsn::AcbValue>>> ax(1);
    ax[0].resize(1);
    ax[0][0] = coeff_list({1});

    auto cm = ode::construct_matrix(dx, ax, 2);
    EXPECT_EQ(cm.block_size, 1);
    EXPECT_EQ(cm.total_order, 2);
    EXPECT_EQ(cm.total_columns, 4);   // (2 + 2) * 1
    EXPECT_EQ(cm.block_systems.size(), 3u);  // total_order + 1
}

TEST(ConstructMatrix, KnownValuesForLinearAxOneOverDx) {
    // dx = 1, ax = [[2*eta]] = [[[0, 2]]], total_order = 1.
    // total_columns = (1+2)*1 = 3.
    auto dx = coeff_list({1});
    std::vector<std::vector<std::vector<dsn::AcbValue>>> ax(1);
    ax[0].resize(1);
    ax[0][0] = coeff_list({0, 2});

    auto cm = ode::construct_matrix(dx, ax, 1);
    EXPECT_EQ(cm.total_columns, 3);
    ASSERT_EQ(cm.block_systems.size(), 2u);

    // i=0 row block (no translation):
    //   j=0: k_dn=0, k_an=-1.  dn_c=1, val=(1+1-0)*1=2 at local col 0.
    //   j=1: k_dn=1, k_an=0.  dn_c=null, axexp[0][0][0]=0 → all-zero; SKIP.
    //   j=2: k_dn=2, k_an=1.  dn_c=null, axexp[0][0][1]=2 → val=-2 at local col 1.
    ASSERT_EQ(cm.block_systems[0].size(), 1u);
    ASSERT_EQ(cm.block_systems[0][0].size(), 2u);
    EXPECT_EQ(cm.block_systems[0][0][0].col, 0);
    EXPECT_TRUE(acb_eq_si(cm.block_systems[0][0][0].val.raw(), 2));
    EXPECT_EQ(cm.block_systems[0][0][1].col, 1);
    EXPECT_TRUE(acb_eq_si(cm.block_systems[0][0][1].val.raw(), -2));

    // i=1 row block (translation N*1=1):
    //   j=1: k_dn=0, k_an=0.  dn_c=1, axexp[0][0][0]=0; val=(1+1-1)*1=1 at local col 0.
    //   j=2: k_dn=1, k_an=0.  dn_c=null, axexp[0][0][0]=0 → all-zero; SKIP.
    // Single entry at col 0 → translated col 1.
    ASSERT_EQ(cm.block_systems[1].size(), 1u);
    ASSERT_EQ(cm.block_systems[1][0].size(), 1u);
    EXPECT_EQ(cm.block_systems[1][0][0].col, 1);
    EXPECT_TRUE(acb_eq_si(cm.block_systems[1][0][0].val.raw(), 1));
}

// ---------------------------------------------------------------------------
//  forward_sparse_gaussian
// ---------------------------------------------------------------------------

TEST(ForwardGaussian, NoDepRowsLeavesSystemUntouched) {
    // System rows have nothing in column 0, vset = {0} — should leave system as-is.
    ode::SparseSystem sys(1);
    sys[0].emplace_back(2, acb_si(7));
    auto step = ode::forward_sparse_gaussian(std::move(sys), {0});
    EXPECT_TRUE(step.reduced.empty());
    ASSERT_EQ(step.residue.size(), 1u);
    ASSERT_EQ(step.residue[0].size(), 1u);
    EXPECT_EQ(step.residue[0][0].col, 2);
}

TEST(ForwardGaussian, SinglePivotNormalisesAndSeparates) {
    // One row [1, 2] — eliminate column 0, normalise leading to 1.
    // After: pivot row should be [(0, 1), (1, 2)] (already normalised since
    // leading is already 1).
    ode::SparseSystem sys(1);
    sys[0].emplace_back(0, acb_si(1));
    sys[0].emplace_back(1, acb_si(2));
    auto step = ode::forward_sparse_gaussian(std::move(sys), {0});
    ASSERT_EQ(step.reduced.size(), 1u);
    ASSERT_EQ(step.reduced[0].size(), 2u);
    EXPECT_TRUE(acb_eq_si(step.reduced[0][0].val.raw(), 1));
    EXPECT_TRUE(acb_eq_si(step.reduced[0][1].val.raw(), 2));
}

TEST(ForwardGaussian, EliminatesColumnFromOtherRows) {
    // Row 0: [2, 4]  (so after pivot normalise: [1, 2])
    // Row 1: [3, 1]  (eliminate col 0: factor = -3, row1 += -3 * pivot
    //                   → [3-3*1, 1-3*2] = [0, -5])
    ode::SparseSystem sys(2);
    sys[0].emplace_back(0, acb_si(2));
    sys[0].emplace_back(1, acb_si(4));
    sys[1].emplace_back(0, acb_si(3));
    sys[1].emplace_back(1, acb_si(1));

    auto step = ode::forward_sparse_gaussian(std::move(sys), {0});
    // Reduced has the pivot.
    ASSERT_EQ(step.reduced.size(), 1u);
    EXPECT_EQ(step.reduced[0][0].col, 0);
    EXPECT_TRUE(acb_eq_si(step.reduced[0][0].val.raw(), 1));
    EXPECT_EQ(step.reduced[0][1].col, 1);
    EXPECT_TRUE(acb_eq_si(step.reduced[0][1].val.raw(), 2));

    // Residue has the eliminated row1 → [(1, -5)].
    ASSERT_EQ(step.residue.size(), 1u);
    ASSERT_EQ(step.residue[0].size(), 1u);
    EXPECT_EQ(step.residue[0][0].col, 1);
    EXPECT_TRUE(acb_eq_si(step.residue[0][0].val.raw(), -5));
}

// ---------------------------------------------------------------------------
//  sparse_gaussian end-to-end smoke test
// ---------------------------------------------------------------------------

TEST(SparseGaussian, EmptyConstructedMatrixThrows) {
    ode::ConstructedMatrix cm;
    cm.block_size = 1;
    cm.total_order = 0;
    cm.total_columns = 2;
    // block_systems empty → throw.
    std::vector<dsn::AcbValue> nh;
    EXPECT_THROW(ode::sparse_gaussian(cm, nh), std::invalid_argument);
}

TEST(SparseGaussian, RunsOnTrivialSystemWithoutThrowing) {
    // dx = 1, ax = [[1]] (constant), total_order = 2.  nh = all zero.
    auto dx = coeff_list({1});
    std::vector<std::vector<std::vector<dsn::AcbValue>>> ax(1);
    ax[0].resize(1);
    ax[0][0] = coeff_list({1});
    auto cm = ode::construct_matrix(dx, ax, 2);

    // total_order+1 = 3 row blocks, each with 1 row → nh length must cover
    // every block_row.
    std::vector<dsn::AcbValue> nh(3);
    auto result = ode::sparse_gaussian(cm, nh);
    // Just verify no throw and reduced rows are returned.
    EXPECT_FALSE(result.reduce.empty());
}

// ---------------------------------------------------------------------------
//  chop_sparse
// ---------------------------------------------------------------------------

TEST(ChopSparse, RemovesValuesBelowThreshold) {
    ode::SparseRow row;
    row.emplace_back(0, acb_si(1));            // |x| = 1 → keep
    {
        // 10^-90 — well below the new sparse-layer chop floor
        // `max(chop_pre, working_pre - 40)` (= max(20, 60) = 60 at the
        // default working_pre = 100).
        dsn::AcbValue tiny;
        arb_set_si(acb_realref(tiny.raw()), 1);
        arb_t denom; arb_init(denom);
        arb_set_si(denom, 10);
        arb_pow_ui(denom, denom, 90u, dsn::working_prec_bits());
        arb_div(acb_realref(tiny.raw()),
                acb_realref(tiny.raw()), denom, dsn::working_prec_bits());
        arb_clear(denom);
        row.emplace_back(1, std::move(tiny));
    }

    // chop_sparse() picks up `sparse_chop_digits()`; at the default
    // working_pre=100 / chop_pre=20 that's 60 digits, so 10^-90 must be
    // chopped while the unit entry survives.
    ode::chop_sparse(row);
    ASSERT_EQ(row.size(), 1u);
    EXPECT_EQ(row[0].col, 0);
}

// ---------------------------------------------------------------------------
//  operator<<
// ---------------------------------------------------------------------------

TEST(SparseRowStream, RendersRuleStyle) {
    ode::SparseRow row;
    row.emplace_back(2, acb_si(7));
    std::ostringstream oss;
    oss << row;
    auto s = oss.str();
    EXPECT_NE(s.find("2 -> "), std::string::npos);
}
