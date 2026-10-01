// SPDX-License-Identifier: MIT
// Tests for amflow::ode::qqbar (Layer 7 helper).

#include <gtest/gtest.h>

#include "amflow/ode/qqbar.hpp"

#include <algorithm>
#include <utility>
#include <vector>

#include <flint/fmpq_mat.h>
#include <flint/fmpq_poly.h>

namespace ode = amflow::ode;

namespace {

void poly_set_si(fmpq_poly_t p, slong n, long v) {
    fmpq_poly_set_coeff_si(p, n, v);
}

void mat_set_si(fmpq_mat_t A, slong i, slong j, long v) {
    fmpq_set_si(fmpq_mat_entry(A, i, j), v, 1);
}

}  // namespace

TEST(QqbarTest, RealIrrationalRootFloorsAreExact) {
    // x^2 - x - 1 has roots golden ratio (~1.618) and -0.618.
    fmpq_poly_t p;
    fmpq_poly_init(p);
    poly_set_si(p, 2, 1);
    poly_set_si(p, 1, -1);
    poly_set_si(p, 0, -1);

    auto roots = ode::algebraic_roots_from_fmpq_poly(p);
    ASSERT_EQ(roots.size(), 2u);

    std::vector<long> floors;
    for (const auto& root : roots) floors.push_back(root.floor_real_si());
    std::sort(floors.begin(), floors.end());

    EXPECT_EQ(floors[0], -1);
    EXPECT_EQ(floors[1], 1);

    fmpq_poly_clear(p);
}

TEST(QqbarTest, ComplexRootsUseExactRealPartFloor) {
    // x^2 - 2x + 2 has roots 1 +- i.  Real part = 1, floor = 1.
    fmpq_poly_t p;
    fmpq_poly_init(p);
    poly_set_si(p, 2, 1);
    poly_set_si(p, 1, -2);
    poly_set_si(p, 0, 2);

    auto roots = ode::algebraic_roots_from_fmpq_poly(p);
    ASSERT_EQ(roots.size(), 2u);

    for (const auto& root : roots) {
        EXPECT_EQ(root.floor_real_si(), 1);
        EXPECT_FALSE(root.to_string().empty());
    }
    fmpq_poly_clear(p);
}

TEST(QqbarTest, MatrixEigenvaluesPreserveExactMultiplicity) {
    // Block-diag with two copies of companion(x^2 - 2). Eigenvalues +-sqrt(2),
    // each multiplicity 2.  Floors: -2 and 1.
    fmpq_mat_t A;
    fmpq_mat_init(A, 4, 4);

    mat_set_si(A, 0, 1, 1);
    mat_set_si(A, 1, 0, 2);
    mat_set_si(A, 2, 3, 1);
    mat_set_si(A, 3, 2, 2);

    auto eigs = ode::algebraic_eigenvalues_from_fmpq_mat(A);
    ASSERT_EQ(eigs.size(), 2u);

    std::vector<std::pair<long, long>> floor_mult;
    for (const auto& ev : eigs) {
        floor_mult.emplace_back(ev.value.floor_real_si(), ev.multiplicity);
        amflow::numeric::AcbValue z = ev.value.to_acb();
        EXPECT_TRUE(z.is_finite());
    }
    std::sort(floor_mult.begin(), floor_mult.end());

    EXPECT_EQ(floor_mult[0].first, -2);
    EXPECT_EQ(floor_mult[0].second, 2);
    EXPECT_EQ(floor_mult[1].first, 1);
    EXPECT_EQ(floor_mult[1].second, 2);

    fmpq_mat_clear(A);
}

TEST(QqbarTest, RationalValueRoundTrip) {
    // QqbarValue from fmpq_t roundtrips and is_rational.
    fmpq_t q;
    fmpq_init(q);
    fmpq_set_si(q, 7, 3);
    ode::QqbarValue v;
    v.set_fmpq(q);
    EXPECT_TRUE(v.is_rational());
    EXPECT_TRUE(v.is_real());
    EXPECT_FALSE(v.is_zero());
    EXPECT_EQ(v.floor_real_si(), 2L);   // floor(7/3) = 2
    fmpq_clear(q);
}

TEST(QqbarTest, ZeroPolyYieldsEmpty) {
    fmpq_poly_t p;
    fmpq_poly_init(p);              // zero poly, degree -1
    auto roots = ode::algebraic_roots_from_fmpq_poly(p);
    EXPECT_TRUE(roots.empty());
    fmpq_poly_clear(p);
}

TEST(QqbarTest, EmptyMatrixYieldsEmptyEigs) {
    fmpq_mat_t A;
    fmpq_mat_init(A, 0, 0);
    auto eigs = ode::algebraic_eigenvalues_from_fmpq_mat(A);
    EXPECT_TRUE(eigs.empty());
    fmpq_mat_clear(A);
}

TEST(QqbarTest, NonSquareMatrixThrows) {
    fmpq_mat_t A;
    fmpq_mat_init(A, 2, 3);
    EXPECT_THROW(ode::algebraic_eigenvalues_from_fmpq_mat(A),
                 std::invalid_argument);
    fmpq_mat_clear(A);
}
