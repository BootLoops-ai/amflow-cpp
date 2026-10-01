// SPDX-License-Identifier: MIT
// algebra::mpoly_matrix — matrices of multivariate polynomials and their
// basic linear algebra (determinant, adjugate, inverse).
//
//
// Used by qft::evaluate_uf for graph polynomial computation:
//
//     U(x)  = Det[A(x)]                       (Symanzik-1)
//     F(x)  = U * (-C + B^T A^{-1} B)         (Symanzik-2)
//
// Algorithms:
//   * det:      Bareiss fraction-free elimination.  Stays in ℤ[x] land by
//               exploiting the integer-divisibility property of Bareiss.
//   * adjugate: cofactor expansion via minors.
//   * inverse:  cofactor / determinant.  Returns MfracMatrix (the inverse
//               of a polynomial matrix is rational in general).

#ifndef AMFLOW_ALGEBRA_MPOLY_MATRIX_HPP
#define AMFLOW_ALGEBRA_MPOLY_MATRIX_HPP

#include <cstddef>
#include <iosfwd>
#include <memory>
#include <vector>

#include "amflow/algebra/mpoly.hpp"

namespace amflow::algebra {

// ---------------------------------------------------------------------------
//  MpolyMatrix  — matrix of Mpoly entries.  Move-only with explicit clone().
// ---------------------------------------------------------------------------

class MpolyMatrix {
public:
    MpolyMatrix() = default;
    MpolyMatrix(std::shared_ptr<MpolyContext> ctx,
                std::size_t rows, std::size_t cols);

    static MpolyMatrix identity(std::shared_ptr<MpolyContext> ctx, std::size_t n);
    static MpolyMatrix zeros  (std::shared_ptr<MpolyContext> ctx,
                                std::size_t rows, std::size_t cols);

    MpolyMatrix(const MpolyMatrix&)            = delete;
    MpolyMatrix& operator=(const MpolyMatrix&) = delete;

    MpolyMatrix(MpolyMatrix&&) noexcept            = default;
    MpolyMatrix& operator=(MpolyMatrix&&) noexcept = default;

    MpolyMatrix clone() const;

    std::size_t rows() const noexcept { return rows_; }
    std::size_t cols() const noexcept { return cols_; }
    bool        empty() const noexcept { return rows_ == 0 || cols_ == 0; }

    const std::shared_ptr<MpolyContext>& ctx() const noexcept { return ctx_; }

    Mpoly&       operator()(std::size_t i, std::size_t j)       { return data_[i * cols_ + j]; }
    const Mpoly& operator()(std::size_t i, std::size_t j) const { return data_[i * cols_ + j]; }
    Mpoly&       at(std::size_t i, std::size_t j);
    const Mpoly& at(std::size_t i, std::size_t j) const;

    MpolyMatrix submatrix(const std::vector<std::size_t>& rows,
                          const std::vector<std::size_t>& cols) const;

    MpolyMatrix minor(std::size_t row_to_drop, std::size_t col_to_drop) const;

    Mpoly       det() const;          // Bareiss
    MpolyMatrix adjugate() const;

    std::string to_string() const;

    bool operator==(const MpolyMatrix& other) const;
    bool operator!=(const MpolyMatrix& other) const { return !(*this == other); }

private:
    std::shared_ptr<MpolyContext> ctx_;
    std::size_t                   rows_ = 0;
    std::size_t                   cols_ = 0;
    std::vector<Mpoly>            data_;
};

std::ostream& operator<<(std::ostream& os, const MpolyMatrix& m);

// ---------------------------------------------------------------------------
//  MfracMatrix  — matrix of Mfrac entries.  Used for the result of inverse().
// ---------------------------------------------------------------------------

class MfracMatrix {
public:
    MfracMatrix() = default;
    MfracMatrix(std::shared_ptr<MpolyContext> ctx,
                std::size_t rows, std::size_t cols);

    static MfracMatrix from_mpoly_matrix(const MpolyMatrix& m);

    MfracMatrix(const MfracMatrix&)            = delete;
    MfracMatrix& operator=(const MfracMatrix&) = delete;

    MfracMatrix(MfracMatrix&&) noexcept            = default;
    MfracMatrix& operator=(MfracMatrix&&) noexcept = default;

    MfracMatrix clone() const;

    std::size_t rows() const noexcept { return rows_; }
    std::size_t cols() const noexcept { return cols_; }
    const std::shared_ptr<MpolyContext>& ctx() const noexcept { return ctx_; }

    Mfrac&       operator()(std::size_t i, std::size_t j)       { return data_[i * cols_ + j]; }
    const Mfrac& operator()(std::size_t i, std::size_t j) const { return data_[i * cols_ + j]; }

    std::string to_string() const;

private:
    std::shared_ptr<MpolyContext> ctx_;
    std::size_t                   rows_ = 0;
    std::size_t                   cols_ = 0;
    std::vector<Mfrac>            data_;
};

// Inverse of a square MpolyMatrix.  Returns false (and leaves out_inverse in
// an unspecified state) if the matrix is singular (det == 0); true otherwise.
bool mpoly_matrix_inverse(MfracMatrix& out_inverse, const MpolyMatrix& m);

}  // namespace amflow::algebra

#endif  // AMFLOW_ALGEBRA_MPOLY_MATRIX_HPP
