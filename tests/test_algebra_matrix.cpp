// SPDX-License-Identifier: MIT
// Tests for amflow::algebra MpolyMatrix / MfracMatrix / mpoly_matrix_inverse.

#include <gtest/gtest.h>

#include "amflow/algebra/mpoly.hpp"
#include "amflow/algebra/mpoly_matrix.hpp"

#include <memory>
#include <stdexcept>

namespace alg = amflow::algebra;

namespace {

std::shared_ptr<alg::MpolyContext> make_ctx(std::vector<std::string> names) {
    return std::make_shared<alg::MpolyContext>(std::move(names));
}

}  // namespace

// ---------------------------------------------------------------------------
//  MpolyMatrix construction
// ---------------------------------------------------------------------------

TEST(MpolyMatrix, IdentityHasOnesOnDiagonal) {
    auto ctx = make_ctx({"x"});
    auto m = alg::MpolyMatrix::identity(ctx, 3);
    EXPECT_EQ(m.rows(), 3u);
    EXPECT_EQ(m.cols(), 3u);
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            if (i == j) EXPECT_TRUE(m(i, j).is_one());
            else        EXPECT_TRUE(m(i, j).is_zero());
        }
    }
}

TEST(MpolyMatrix, ZerosIsZero) {
    auto ctx = make_ctx({"x"});
    auto m = alg::MpolyMatrix::zeros(ctx, 2, 3);
    EXPECT_EQ(m.rows(), 2u);
    EXPECT_EQ(m.cols(), 3u);
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 3; ++j)
            EXPECT_TRUE(m(i, j).is_zero());
}

TEST(MpolyMatrix, AtBoundsCheck) {
    auto ctx = make_ctx({"x"});
    auto m = alg::MpolyMatrix(ctx, 2, 2);
    EXPECT_THROW(m.at(2, 0), std::out_of_range);
    EXPECT_THROW(m.at(0, 5), std::out_of_range);
}

TEST(MpolyMatrix, Submatrix) {
    auto ctx = make_ctx({"x"});
    alg::MpolyMatrix m(ctx, 3, 3);
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j)
            m(i, j) = alg::Mpoly::constant(ctx, static_cast<long>(i * 3 + j));
    auto sub = m.submatrix({0, 2}, {1, 2});
    EXPECT_EQ(sub.rows(), 2u);
    EXPECT_EQ(sub.cols(), 2u);
    EXPECT_EQ(sub(0, 0), alg::Mpoly::constant(ctx, 1));
    EXPECT_EQ(sub(0, 1), alg::Mpoly::constant(ctx, 2));
    EXPECT_EQ(sub(1, 0), alg::Mpoly::constant(ctx, 7));
    EXPECT_EQ(sub(1, 1), alg::Mpoly::constant(ctx, 8));
}

TEST(MpolyMatrix, MinorRemovesRowAndColumn) {
    auto ctx = make_ctx({"x"});
    alg::MpolyMatrix m(ctx, 3, 3);
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j)
            m(i, j) = alg::Mpoly::constant(ctx, static_cast<long>(i * 3 + j));
    auto mm = m.minor(1, 1);   // drop row 1, col 1 → 2x2 of {0,2,6,8}
    EXPECT_EQ(mm.rows(), 2u);
    EXPECT_EQ(mm.cols(), 2u);
    EXPECT_EQ(mm(0, 0), alg::Mpoly::constant(ctx, 0L));
    EXPECT_EQ(mm(0, 1), alg::Mpoly::constant(ctx, 2));
    EXPECT_EQ(mm(1, 0), alg::Mpoly::constant(ctx, 6));
    EXPECT_EQ(mm(1, 1), alg::Mpoly::constant(ctx, 8));
}

TEST(MpolyMatrix, MinorOutOfRangeThrows) {
    auto ctx = make_ctx({"x"});
    alg::MpolyMatrix m(ctx, 2, 2);
    EXPECT_THROW(m.minor(3, 0), std::out_of_range);
    EXPECT_THROW(m.minor(0, 3), std::out_of_range);
}

// ---------------------------------------------------------------------------
//  Determinant (Bareiss)
// ---------------------------------------------------------------------------

TEST(MpolyMatrixDet, EmptyIsOne) {
    auto ctx = make_ctx({"x"});
    auto m = alg::MpolyMatrix::zeros(ctx, 0, 0);
    auto d = m.det();
    EXPECT_TRUE(d.is_one());
}

TEST(MpolyMatrixDet, OneByOneIsEntry) {
    auto ctx = make_ctx({"x"});
    alg::MpolyMatrix m(ctx, 1, 1);
    m(0, 0) = alg::Mpoly::from_string(ctx, "x + 5");
    auto d = m.det();
    auto expected = alg::Mpoly::from_string(ctx, "x + 5");
    EXPECT_EQ(d, expected);
}

TEST(MpolyMatrixDet, TwoByTwoSymbolic) {
    auto ctx = make_ctx({"x"});
    // [[x, 1], [2, x]] -> x^2 - 2
    alg::MpolyMatrix m(ctx, 2, 2);
    m(0, 0) = alg::Mpoly::from_string(ctx, "x");
    m(0, 1) = alg::Mpoly::from_string(ctx, "1");
    m(1, 0) = alg::Mpoly::from_string(ctx, "2");
    m(1, 1) = alg::Mpoly::from_string(ctx, "x");
    auto d = m.det();
    auto expected = alg::Mpoly::from_string(ctx, "x^2 - 2");
    EXPECT_EQ(d, expected);
}

TEST(MpolyMatrixDet, IdentityIsOne) {
    auto ctx = make_ctx({"x"});
    auto m = alg::MpolyMatrix::identity(ctx, 4);
    auto d = m.det();
    EXPECT_TRUE(d.is_one());
}

TEST(MpolyMatrixDet, SingularMatrixIsZero) {
    auto ctx = make_ctx({"x"});
    // Two equal rows.
    alg::MpolyMatrix m(ctx, 3, 3);
    m(0, 0) = alg::Mpoly::from_string(ctx, "1");
    m(0, 1) = alg::Mpoly::from_string(ctx, "x");
    m(0, 2) = alg::Mpoly::from_string(ctx, "x^2");
    m(1, 0) = alg::Mpoly::from_string(ctx, "1");
    m(1, 1) = alg::Mpoly::from_string(ctx, "x");
    m(1, 2) = alg::Mpoly::from_string(ctx, "x^2");
    m(2, 0) = alg::Mpoly::from_string(ctx, "0");
    m(2, 1) = alg::Mpoly::from_string(ctx, "1");
    m(2, 2) = alg::Mpoly::from_string(ctx, "0");
    auto d = m.det();
    EXPECT_TRUE(d.is_zero());
}

TEST(MpolyMatrixDet, ThreeByThreeMultivar) {
    auto ctx = make_ctx({"x", "y"});
    // [[x, y, 0], [0, x, y], [y, 0, x]]
    // det = x*(x*x - y*0) - y*(0*x - y*y) + 0 = x^3 + y^3
    alg::MpolyMatrix m(ctx, 3, 3);
    m(0, 0) = alg::Mpoly::from_string(ctx, "x");
    m(0, 1) = alg::Mpoly::from_string(ctx, "y");
    m(1, 1) = alg::Mpoly::from_string(ctx, "x");
    m(1, 2) = alg::Mpoly::from_string(ctx, "y");
    m(2, 0) = alg::Mpoly::from_string(ctx, "y");
    m(2, 2) = alg::Mpoly::from_string(ctx, "x");
    auto d = m.det();
    auto expected = alg::Mpoly::from_string(ctx, "x^3 + y^3");
    EXPECT_EQ(d, expected);
}

TEST(MpolyMatrixDet, NonSquareThrows) {
    auto ctx = make_ctx({"x"});
    auto m = alg::MpolyMatrix::zeros(ctx, 2, 3);
    EXPECT_THROW(m.det(), std::invalid_argument);
}

// ---------------------------------------------------------------------------
//  Adjugate
// ---------------------------------------------------------------------------

TEST(MpolyMatrixAdj, IdentityIsItself) {
    auto ctx = make_ctx({"x"});
    auto m = alg::MpolyMatrix::identity(ctx, 3);
    auto adj = m.adjugate();
    auto id = alg::MpolyMatrix::identity(ctx, 3);
    EXPECT_EQ(adj, id);
}

TEST(MpolyMatrixAdj, AdjAtimesAisDetTimesIdentity) {
    auto ctx = make_ctx({"x"});
    // [[x, 1], [2, x]] -- det = x^2 - 2.  adj * self = (x^2 - 2) * I.
    alg::MpolyMatrix a(ctx, 2, 2);
    a(0, 0) = alg::Mpoly::from_string(ctx, "x");
    a(0, 1) = alg::Mpoly::from_string(ctx, "1");
    a(1, 0) = alg::Mpoly::from_string(ctx, "2");
    a(1, 1) = alg::Mpoly::from_string(ctx, "x");

    auto adj = a.adjugate();
    auto det = a.det();

    // Compute adj * a manually.
    long n = 2;
    alg::MpolyMatrix prod(ctx, n, n);
    for (long i = 0; i < n; ++i)
        for (long j = 0; j < n; ++j) {
            alg::Mpoly s(ctx);
            for (long k = 0; k < n; ++k) {
                s += adj(i, k) * a(k, j);
            }
            prod(i, j) = std::move(s);
        }

    // Expect prod == det * identity.
    for (long i = 0; i < n; ++i) {
        for (long j = 0; j < n; ++j) {
            if (i == j) EXPECT_EQ(prod(i, j), det);
            else        EXPECT_TRUE(prod(i, j).is_zero());
        }
    }
}

TEST(MpolyMatrixAdj, NonSquareThrows) {
    auto ctx = make_ctx({"x"});
    auto m = alg::MpolyMatrix::zeros(ctx, 2, 3);
    EXPECT_THROW(m.adjugate(), std::invalid_argument);
}

// ---------------------------------------------------------------------------
//  MfracMatrix + mpoly_matrix_inverse
// ---------------------------------------------------------------------------

TEST(MfracMatrix, FromMpolyMatrix) {
    auto ctx = make_ctx({"x"});
    auto m = alg::MpolyMatrix::identity(ctx, 2);
    auto mf = alg::MfracMatrix::from_mpoly_matrix(m);
    EXPECT_EQ(mf.rows(), 2u);
    EXPECT_EQ(mf.cols(), 2u);
    EXPECT_TRUE(mf(0, 0).is_one());
    EXPECT_TRUE(mf(1, 0).is_zero());
}

TEST(MpolyMatrixInverse, IdentityIsItself) {
    auto ctx = make_ctx({"x"});
    auto id = alg::MpolyMatrix::identity(ctx, 3);
    alg::MfracMatrix inv;
    ASSERT_TRUE(alg::mpoly_matrix_inverse(inv, id));
    EXPECT_EQ(inv.rows(), 3u);
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            if (i == j) EXPECT_TRUE(inv(i, j).is_one());
            else        EXPECT_TRUE(inv(i, j).is_zero());
        }
    }
}

TEST(MpolyMatrixInverse, SingularReturnsFalse) {
    auto ctx = make_ctx({"x"});
    // Two equal rows -> singular.
    alg::MpolyMatrix m(ctx, 2, 2);
    m(0, 0) = alg::Mpoly::from_string(ctx, "x");
    m(0, 1) = alg::Mpoly::from_string(ctx, "x");
    m(1, 0) = alg::Mpoly::from_string(ctx, "x");
    m(1, 1) = alg::Mpoly::from_string(ctx, "x");
    alg::MfracMatrix inv;
    EXPECT_FALSE(alg::mpoly_matrix_inverse(inv, m));
}

TEST(MpolyMatrixInverse, InverseTimesOriginalIsIdentity) {
    auto ctx = make_ctx({"x"});
    // [[x, 1], [2, x]] -- det = x^2 - 2; inverse = adj / det.
    alg::MpolyMatrix a(ctx, 2, 2);
    a(0, 0) = alg::Mpoly::from_string(ctx, "x");
    a(0, 1) = alg::Mpoly::from_string(ctx, "1");
    a(1, 0) = alg::Mpoly::from_string(ctx, "2");
    a(1, 1) = alg::Mpoly::from_string(ctx, "x");

    alg::MfracMatrix inv;
    ASSERT_TRUE(alg::mpoly_matrix_inverse(inv, a));

    // inv * a should equal identity.
    long n = 2;
    auto a_mfrac = alg::MfracMatrix::from_mpoly_matrix(a);
    for (long i = 0; i < n; ++i) {
        for (long j = 0; j < n; ++j) {
            alg::Mfrac s = alg::Mfrac::zero(ctx);
            for (long k = 0; k < n; ++k) {
                s += inv(i, k) * a_mfrac(k, j);
            }
            if (i == j) {
                EXPECT_TRUE(s.is_one()) << "i=" << i << " j=" << j << " s=" << s.to_string();
            } else {
                EXPECT_TRUE(s.is_zero()) << "i=" << i << " j=" << j << " s=" << s.to_string();
            }
        }
    }
}

TEST(MpolyMatrixInverse, NonSquareThrows) {
    auto ctx = make_ctx({"x"});
    auto m = alg::MpolyMatrix::zeros(ctx, 2, 3);
    alg::MfracMatrix inv;
    EXPECT_THROW(alg::mpoly_matrix_inverse(inv, m), std::invalid_argument);
}
