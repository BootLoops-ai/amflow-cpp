// SPDX-License-Identifier: MIT
// Tests for amflow::ode::jordan (Layer 7 helper).

#include <gtest/gtest.h>

#include "amflow/ode/jordan.hpp"
#include "amflow/numeric/rational.hpp"

#include <flint/fmpq.h>
#include <flint/fmpq_mat.h>

#include <algorithm>
#include <functional>
#include <stdexcept>
#include <vector>

namespace ode = amflow::ode;
namespace nm  = amflow::numeric;

namespace {

void set_si(fmpq_mat_t M, slong i, slong j, long v) {
    fmpq_set_si(fmpq_mat_entry(M, i, j), v, 1);
}

void expect_reconstructs(const fmpq_mat_t A,
                         const fmpq_mat_t S,
                         const fmpq_mat_t J,
                         const fmpq_mat_t Sinv) {
    slong n = fmpq_mat_nrows(A);
    fmpq_mat_t SJ, recon, I_check, I_ref;
    fmpq_mat_init(SJ, n, n);
    fmpq_mat_init(recon, n, n);
    fmpq_mat_init(I_check, n, n);
    fmpq_mat_init(I_ref, n, n);

    fmpq_mat_mul(SJ, S, J);
    fmpq_mat_mul(recon, SJ, Sinv);
    EXPECT_TRUE(fmpq_mat_equal(recon, A));

    fmpq_mat_mul(I_check, S, Sinv);
    fmpq_mat_one(I_ref);
    EXPECT_TRUE(fmpq_mat_equal(I_check, I_ref));

    fmpq_mat_clear(I_ref);
    fmpq_mat_clear(I_check);
    fmpq_mat_clear(recon);
    fmpq_mat_clear(SJ);
}

}  // namespace

TEST(JordanTest, SimilarityTransformedJordanBlock) {
    fmpq_mat_t Jref, P, Pinv, PJ, A, S, J, Sinv;
    fmpq_mat_init(Jref, 2, 2);
    fmpq_mat_init(P, 2, 2);
    fmpq_mat_init(Pinv, 2, 2);
    fmpq_mat_init(PJ, 2, 2);
    fmpq_mat_init(A, 2, 2);
    fmpq_mat_init(S, 2, 2);
    fmpq_mat_init(J, 2, 2);
    fmpq_mat_init(Sinv, 2, 2);

    set_si(Jref, 0, 0, 3);
    set_si(Jref, 0, 1, 1);
    set_si(Jref, 1, 1, 3);
    set_si(P, 0, 0, 1);
    set_si(P, 0, 1, 2);
    set_si(P, 1, 0, 3);
    set_si(P, 1, 1, 4);
    ASSERT_NE(fmpq_mat_inv(Pinv, P), 0);
    fmpq_mat_mul(PJ, P, Jref);
    fmpq_mat_mul(A, PJ, Pinv);

    std::vector<long> blocks;
    ode::jordan_decomposition_exact(S, J, Sinv, blocks, A);

    ASSERT_EQ(blocks.size(), 1u);
    EXPECT_EQ(blocks[0], 2);
    EXPECT_TRUE(fmpq_mat_equal(J, Jref));
    expect_reconstructs(A, S, J, Sinv);

    fmpq_mat_clear(Sinv);
    fmpq_mat_clear(J);
    fmpq_mat_clear(S);
    fmpq_mat_clear(A);
    fmpq_mat_clear(PJ);
    fmpq_mat_clear(Pinv);
    fmpq_mat_clear(P);
    fmpq_mat_clear(Jref);
}

TEST(JordanTest, MixedBlocksSameEigenvalue) {
    // A = diag(2, 2, 2) with one off-diagonal making a 2-block + 1-block.
    fmpq_mat_t A, S, J, Sinv;
    fmpq_mat_init(A, 3, 3);
    fmpq_mat_init(S, 3, 3);
    fmpq_mat_init(J, 3, 3);
    fmpq_mat_init(Sinv, 3, 3);

    set_si(A, 0, 0, 2);
    set_si(A, 0, 1, 1);
    set_si(A, 1, 1, 2);
    set_si(A, 2, 2, 2);

    std::vector<long> blocks;
    ode::jordan_decomposition_exact(S, J, Sinv, blocks, A);

    ASSERT_EQ(blocks.size(), 2u);
    EXPECT_EQ(blocks[0], 2);
    EXPECT_EQ(blocks[1], 1);
    expect_reconstructs(A, S, J, Sinv);

    fmpq_mat_clear(Sinv);
    fmpq_mat_clear(J);
    fmpq_mat_clear(S);
    fmpq_mat_clear(A);
}

// --- Multi-distinct-eigenvalue cases -
//
// Existing tests above all use a single distinct eigenvalue.  The
// audit's row-190 concern is: C++ `jordan_decomposition_exact`
// groups blocks by eigenvalue (outer loop in `jordan.cpp:265`,
// `for (auto& [lambda, mult] : eigs)`), then within each
// eigenvalue places blocks in *descending size* order.  Upstream
// MMA `JordanDecomposition` sorts blocks globally by descending
// size, mixing eigenvalues.  The two orderings differ when a
// matrix has ≥2 distinct eigenvalues with non-trivial size mix.
//
// The audit calls this "same final algebra, different column
// permutation."  These tests lock that in: for ≥2-distinct-
// eigenvalue inputs we verify the decomposition's algebraic
// invariants (reconstruction `S J S^{-1} = A`, block-size
// multiset, eigenvalue multiset) without depending on the
// specific column order.  The end-to-end no-impact-on-final-
// integral is independently asserted by the 12 oracle benchmarks
// matching MMA at rel ~ 10^{-30}.

TEST(JordanTest, TwoDistinctEigenvalues_SimpleDiagonal) {
    // A = diag(3, 5) — two distinct simple eigenvalues.
    fmpq_mat_t A, S, J, Sinv;
    fmpq_mat_init(A, 2, 2);
    fmpq_mat_init(S, 2, 2);
    fmpq_mat_init(J, 2, 2);
    fmpq_mat_init(Sinv, 2, 2);
    set_si(A, 0, 0, 3);
    set_si(A, 1, 1, 5);

    std::vector<long> blocks;
    ode::jordan_decomposition_exact(S, J, Sinv, blocks, A);

    ASSERT_EQ(blocks.size(), 2u);
    EXPECT_EQ(blocks[0], 1);
    EXPECT_EQ(blocks[1], 1);
    expect_reconstructs(A, S, J, Sinv);

    fmpq_mat_clear(Sinv);
    fmpq_mat_clear(J);
    fmpq_mat_clear(S);
    fmpq_mat_clear(A);
}

TEST(JordanTest, TwoDistinctEigenvalues_MixedBlockSizes) {
    // A is similar to diag(J_2(2), J_1(7)).  C++ groups by
    // eigenvalue: λ=2 first (block size 2), then λ=7 (block size 1).
    // MMA would sort globally: [2-block @ λ=2, 1-block @ λ=7] —
    // happens to be the same order here.  But the algebraic
    // invariants (reconstruction, set of eigenvalues, multiset of
    // block sizes) are basis-invariant either way.
    fmpq_mat_t A, S, J, Sinv;
    fmpq_mat_init(A, 3, 3);
    fmpq_mat_init(S, 3, 3);
    fmpq_mat_init(J, 3, 3);
    fmpq_mat_init(Sinv, 3, 3);

    // J_2(2)+J_1(7) in upper-triangular form.
    set_si(A, 0, 0, 2); set_si(A, 0, 1, 1);
    set_si(A, 1, 1, 2);
    set_si(A, 2, 2, 7);

    std::vector<long> blocks;
    ode::jordan_decomposition_exact(S, J, Sinv, blocks, A);

    // Multiset of block sizes.
    ASSERT_EQ(blocks.size(), 2u);
    std::vector<long> sizes_sorted = blocks;
    std::sort(sizes_sorted.begin(), sizes_sorted.end(), std::greater<long>());
    EXPECT_EQ(sizes_sorted[0], 2);
    EXPECT_EQ(sizes_sorted[1], 1);

    // Eigenvalue multiset on the diagonal of J.
    std::vector<long> eig_diag;
    for (slong i = 0; i < 3; ++i) {
        const fmpq* e = fmpq_mat_entry(J, i, i);
        ASSERT_TRUE(fmpz_is_one(fmpq_denref(e))) << "diagonal not integer";
        eig_diag.push_back(fmpz_get_si(fmpq_numref(e)));
    }
    std::sort(eig_diag.begin(), eig_diag.end());
    EXPECT_EQ(eig_diag, (std::vector<long>{2, 2, 7}));

    expect_reconstructs(A, S, J, Sinv);

    fmpq_mat_clear(Sinv);
    fmpq_mat_clear(J);
    fmpq_mat_clear(S);
    fmpq_mat_clear(A);
}

TEST(JordanTest, ThreeDistinctEigenvalues_SimilarityTransformed) {
    // Build A = P * diag(1, 2, 3) * P^{-1} where P is a fixed
    // unimodular integer matrix.  Eigenvalues {1, 2, 3} — three
    // distinct simple ones.  Tests the multi-distinct-eigenvalue
    // outer-loop path with a non-trivial similarity transform so
    // that the result is not already diagonal.
    fmpq_mat_t Jref, P, Pinv, PJ, A, S, J, Sinv;
    fmpq_mat_init(Jref, 3, 3);
    fmpq_mat_init(P,    3, 3);
    fmpq_mat_init(Pinv, 3, 3);
    fmpq_mat_init(PJ,   3, 3);
    fmpq_mat_init(A,    3, 3);
    fmpq_mat_init(S,    3, 3);
    fmpq_mat_init(J,    3, 3);
    fmpq_mat_init(Sinv, 3, 3);

    set_si(Jref, 0, 0, 1);
    set_si(Jref, 1, 1, 2);
    set_si(Jref, 2, 2, 3);

    // P = [[1,1,0],[0,1,1],[0,0,1]]  (upper unitriangular, det = 1).
    set_si(P, 0, 0, 1); set_si(P, 0, 1, 1);
    set_si(P, 1, 1, 1); set_si(P, 1, 2, 1);
    set_si(P, 2, 2, 1);
    // P^{-1} for this upper unitriangular: invert directly.
    set_si(Pinv, 0, 0, 1); set_si(Pinv, 0, 1, -1); set_si(Pinv, 0, 2, 1);
    set_si(Pinv, 1, 1, 1); set_si(Pinv, 1, 2, -1);
    set_si(Pinv, 2, 2, 1);

    fmpq_mat_mul(PJ, P, Jref);
    fmpq_mat_mul(A,  PJ, Pinv);

    std::vector<long> blocks;
    ode::jordan_decomposition_exact(S, J, Sinv, blocks, A);

    ASSERT_EQ(blocks.size(), 3u);
    for (long b : blocks) EXPECT_EQ(b, 1);

    // Eigenvalue multiset on diagonal of J.
    std::vector<long> eig_diag;
    for (slong i = 0; i < 3; ++i) {
        const fmpq* e = fmpq_mat_entry(J, i, i);
        ASSERT_TRUE(fmpz_is_one(fmpq_denref(e)));
        eig_diag.push_back(fmpz_get_si(fmpq_numref(e)));
    }
    std::sort(eig_diag.begin(), eig_diag.end());
    EXPECT_EQ(eig_diag, (std::vector<long>{1, 2, 3}));

    expect_reconstructs(A, S, J, Sinv);

    fmpq_mat_clear(Sinv);
    fmpq_mat_clear(J);
    fmpq_mat_clear(S);
    fmpq_mat_clear(A);
    fmpq_mat_clear(PJ);
    fmpq_mat_clear(Pinv);
    fmpq_mat_clear(P);
    fmpq_mat_clear(Jref);
}

TEST(JordanTest, TwoDistinctEigenvalues_BothWithJordanBlocks) {
    // Most demanding case: λ=2 with a 2-block, λ=−1 with a 2-block.
    // 4×4 matrix where each eigenvalue has algebraic multiplicity 2
    // and a single Jordan chain.  Locks down that the "by-eigenvalue
    // outer loop, descending block size inner" ordering does not
    // miss any chains when ≥2 distinct eigenvalues each carry a
    // non-trivial chain.
    fmpq_mat_t A, S, J, Sinv;
    fmpq_mat_init(A, 4, 4);
    fmpq_mat_init(S, 4, 4);
    fmpq_mat_init(J, 4, 4);
    fmpq_mat_init(Sinv, 4, 4);

    // Block-diagonal: J_2(2) ⊕ J_2(-1).
    set_si(A, 0, 0, 2); set_si(A, 0, 1, 1);
    set_si(A, 1, 1, 2);
    set_si(A, 2, 2, -1); set_si(A, 2, 3, 1);
    set_si(A, 3, 3, -1);

    std::vector<long> blocks;
    ode::jordan_decomposition_exact(S, J, Sinv, blocks, A);

    ASSERT_EQ(blocks.size(), 2u);
    EXPECT_EQ(blocks[0], 2);
    EXPECT_EQ(blocks[1], 2);

    std::vector<long> eig_diag;
    for (slong i = 0; i < 4; ++i) {
        const fmpq* e = fmpq_mat_entry(J, i, i);
        ASSERT_TRUE(fmpz_is_one(fmpq_denref(e)));
        eig_diag.push_back(fmpz_get_si(fmpq_numref(e)));
    }
    std::sort(eig_diag.begin(), eig_diag.end());
    EXPECT_EQ(eig_diag, (std::vector<long>{-1, -1, 2, 2}));

    expect_reconstructs(A, S, J, Sinv);

    fmpq_mat_clear(Sinv);
    fmpq_mat_clear(J);
    fmpq_mat_clear(S);
    fmpq_mat_clear(A);
}

TEST(JordanTest, IrrationalEigenvalueThrows) {
    // A = [[0,1],[2,0]] -> eigenvalues +-sqrt(2) are not in Q.
    fmpq_mat_t A, S, J, Sinv;
    fmpq_mat_init(A, 2, 2);
    fmpq_mat_init(S, 2, 2);
    fmpq_mat_init(J, 2, 2);
    fmpq_mat_init(Sinv, 2, 2);

    set_si(A, 0, 1, 1);
    set_si(A, 1, 0, 2);

    std::vector<long> blocks;
    EXPECT_THROW(ode::jordan_decomposition_exact(S, J, Sinv, blocks, A),
                 std::runtime_error);

    fmpq_mat_clear(Sinv);
    fmpq_mat_clear(J);
    fmpq_mat_clear(S);
    fmpq_mat_clear(A);
}

TEST(JordanTest, FmpqMatFromRationalResidue_PlainConstantMatrix) {
    // mat = [[2, 0], [0, 3]] (constant entries) -> shift = 0 -> just reads
    // num/den constant terms.
    nm::RationalMatrix m(2, 2);
    m(0, 0) = nm::RationalFunction::from_si(2);
    m(0, 1) = nm::RationalFunction();
    m(1, 0) = nm::RationalFunction();
    m(1, 1) = nm::RationalFunction::from_si(3);

    fmpq_mat_t out;
    fmpq_mat_init(out, 2, 2);
    ode::fmpq_mat_from_rational_residue(out, m, 0);

    fmpq_t two, three;
    fmpq_init(two);   fmpq_set_si(two, 2, 1);
    fmpq_init(three); fmpq_set_si(three, 3, 1);
    EXPECT_TRUE(fmpq_equal(fmpq_mat_entry(out, 0, 0), two));
    EXPECT_TRUE(fmpq_is_zero(fmpq_mat_entry(out, 0, 1)));
    EXPECT_TRUE(fmpq_is_zero(fmpq_mat_entry(out, 1, 0)));
    EXPECT_TRUE(fmpq_equal(fmpq_mat_entry(out, 1, 1), three));
    fmpq_clear(two);
    fmpq_clear(three);
    fmpq_mat_clear(out);
}

TEST(JordanTest, FmpqMatFromRationalResidue_ShiftRemovesPole) {
    // mat = [[5/eta]] -> shift 1 -> [[5]].
    nm::FmpqPoly num; num.set_coeff_si(0, 5);
    nm::FmpqPoly den; den.set_coeff_si(1, 1);
    nm::RationalMatrix m(1, 1);
    m(0, 0) = nm::RationalFunction(std::move(num), std::move(den));

    fmpq_mat_t out;
    fmpq_mat_init(out, 1, 1);
    ode::fmpq_mat_from_rational_residue(out, m, 1);

    fmpq_t five;
    fmpq_init(five); fmpq_set_si(five, 5, 1);
    EXPECT_TRUE(fmpq_equal(fmpq_mat_entry(out, 0, 0), five));
    fmpq_clear(five);
    fmpq_mat_clear(out);
}

TEST(JordanTest, FmpqMatFromRationalResidue_RemainingPoleThrows) {
    // mat = [[5/eta]] -> shift 0 leaves the pole.
    nm::FmpqPoly num; num.set_coeff_si(0, 5);
    nm::FmpqPoly den; den.set_coeff_si(1, 1);
    nm::RationalMatrix m(1, 1);
    m(0, 0) = nm::RationalFunction(std::move(num), std::move(den));

    fmpq_mat_t out;
    fmpq_mat_init(out, 1, 1);
    EXPECT_THROW(ode::fmpq_mat_from_rational_residue(out, m, 0),
                 std::runtime_error);
    fmpq_mat_clear(out);
}

TEST(JordanTest, RationalMatrixFromFmpqMat_RoundTrip) {
    fmpq_mat_t A;
    fmpq_mat_init(A, 2, 2);
    set_si(A, 0, 0, 7);
    fmpq_set_si(fmpq_mat_entry(A, 1, 1), 3, 2);   // 3/2
    auto rm = ode::rational_matrix_from_fmpq_mat(A);
    ASSERT_EQ(rm.rows(), 2u);
    ASSERT_EQ(rm.cols(), 2u);
    EXPECT_TRUE(rm(0, 1).is_zero());
    EXPECT_TRUE(rm(1, 0).is_zero());

    fmpq_mat_t out;
    fmpq_mat_init(out, 2, 2);
    ode::fmpq_mat_from_rational_residue(out, rm, 0);
    EXPECT_TRUE(fmpq_mat_equal(A, out));
    fmpq_mat_clear(A);
    fmpq_mat_clear(out);
}

// ============================================================================
//  Eigenvector normalization regression test (pentabox 76-master sub-system)
//
//  Bug: FLINT's fmpz_mat_nullspace (used internally by
//  jordan_decomposition_exact's eigenvector search) returns vectors with
//  denominator-cleared INTEGER entries.  For an eigenvector that Mathematica
//  would emit as (6993/998, 1) FLINT returns (6993, 998).  Without
//  rescaling, downstream shearing / leading-Jordan T blocks accumulate huge
//  integer scale factors, cascading through the off-diagonal Sylvester step
//  in to_fuchsian_global to produce T entries up to 10^300+ in hard cases
//  (the pentabox 2L 5-leg 76-master sub-system) -- overwhelming any
//  practical working precision in PSMapRuleS.
//
//  Fix: normalize each null-space basis vector so its LAST non-zero entry
//  is 1, matching Mathematica's Eigenvectors / JordanDecomposition
//  convention.
//
//  Test matrix:  A = [[999/1000, 0], [499/3500, 0]]
//                (the post-shearing leading residue of the 2-master block
//                 {7,8} that surfaced the pentabox failure.)
//  Eigenvalues:  0 and 999/1000.
//  Expected (Mathematica-convention) eigenvectors:
//                col_for_zero      = (0, 1)
//                col_for_999/1000  = (6993/998, 1)
//
//  Without the fix, col_for_999/1000 = (6993, 998), 998x larger.
// ============================================================================
TEST(JordanTest, JordanEigenvectorsNormalizedToLastEntryOne) {
    fmpq_mat_t A;
    fmpq_mat_init(A, 2, 2);
    fmpq_set_si(fmpq_mat_entry(A, 0, 0), 999, 1000);
    fmpq_set_si(fmpq_mat_entry(A, 1, 0), 499, 3500);

    fmpq_mat_t S, J, Sinv;
    fmpq_mat_init(S, 2, 2);
    fmpq_mat_init(J, 2, 2);
    fmpq_mat_init(Sinv, 2, 2);
    std::vector<long> blocks;
    ode::jordan_decomposition_exact(S, J, Sinv, blocks, A);
    expect_reconstructs(A, S, J, Sinv);

    // Each column of S is a (generalized) eigenvector.  The last non-zero
    // entry must equal 1 in every column.
    for (slong c = 0; c < 2; ++c) {
        slong last_nz = -1;
        for (slong r = 1; r >= 0; --r) {
            if (!fmpq_is_zero(fmpq_mat_entry(S, r, c))) {
                last_nz = r;
                break;
            }
        }
        ASSERT_GE(last_nz, 0) << "S column " << c << " is all zero";
        EXPECT_TRUE(fmpq_is_one(fmpq_mat_entry(S, last_nz, c)))
            << "S col " << c << " last non-zero entry must equal 1";
    }

    // The column whose corresponding eigenvalue is 999/1000 must have its
    // first entry equal to 6993/998.
    fmpq_t expected_first;
    fmpq_init(expected_first);
    fmpq_set_si(expected_first, 6993, 998);

    bool found_nonzero_eig = false;
    for (slong c = 0; c < 2; ++c) {
        if (!fmpq_is_zero(fmpq_mat_entry(J, c, c))) {
            EXPECT_TRUE(fmpq_equal(fmpq_mat_entry(S, 0, c), expected_first))
                << "S col " << c << " first entry must equal 6993/998";
            found_nonzero_eig = true;
        }
    }
    EXPECT_TRUE(found_nonzero_eig);

    fmpq_clear(expected_first);
    fmpq_mat_clear(Sinv);
    fmpq_mat_clear(J);
    fmpq_mat_clear(S);
    fmpq_mat_clear(A);
}
