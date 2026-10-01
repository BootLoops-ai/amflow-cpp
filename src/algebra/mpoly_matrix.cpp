// SPDX-License-Identifier: MIT
// algebra::mpoly_matrix — implementation.
//
// amflow::algebra.
//
// Bareiss algorithm reference:
//   Bareiss (1968), "Sylvester's identity and multistep integer-preserving
//   Gaussian elimination."  Math. Comp. 22, 565-578.

#include "amflow/algebra/mpoly_matrix.hpp"

#include <algorithm>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <utility>

#include <flint/fmpz.h>
#include <flint/fmpz_mpoly.h>

namespace amflow::algebra {

// ===========================================================================
//  Mpoly::exact_divide  (Layer 11a addition that lives here for locality)
// ===========================================================================

bool Mpoly::exact_divide(Mpoly& out, const Mpoly& a, const Mpoly& b) {
    if (a.ctx().get() != b.ctx().get() || a.ctx().get() != out.ctx().get()) {
        throw std::invalid_argument("Mpoly::exact_divide: context mismatch");
    }
    int divides = fmpz_mpoly_divides(out.raw(), a.raw(), b.raw(), a.ctx()->raw());
    return divides != 0;
}

Mpoly Mpoly::gcd(const Mpoly& a, const Mpoly& b) {
    if (a.ctx().get() != b.ctx().get()) {
        throw std::invalid_argument("Mpoly::gcd: context mismatch");
    }
    Mpoly out(a.ctx());
    int ok = fmpz_mpoly_gcd(out.raw(), a.raw(), b.raw(), a.ctx()->raw());
    if (!ok) {
        throw std::runtime_error("Mpoly::gcd: FLINT gcd failed");
    }
    return out;
}

// ===========================================================================
//  MpolyMatrix
// ===========================================================================

MpolyMatrix::MpolyMatrix(std::shared_ptr<MpolyContext> ctx,
                          std::size_t rows, std::size_t cols)
    : ctx_(std::move(ctx)), rows_(rows), cols_(cols) {
    data_.reserve(rows_ * cols_);
    for (std::size_t k = 0; k < rows_ * cols_; ++k) {
        data_.emplace_back(ctx_);
    }
}

MpolyMatrix MpolyMatrix::identity(std::shared_ptr<MpolyContext> ctx, std::size_t n) {
    MpolyMatrix m(ctx, n, n);
    for (std::size_t i = 0; i < n; ++i) {
        m(i, i) = Mpoly::one(ctx);
    }
    return m;
}

MpolyMatrix MpolyMatrix::zeros(std::shared_ptr<MpolyContext> ctx,
                                std::size_t rows, std::size_t cols) {
    return MpolyMatrix(std::move(ctx), rows, cols);
}

MpolyMatrix MpolyMatrix::clone() const {
    MpolyMatrix out(ctx_, rows_, cols_);
    for (std::size_t k = 0; k < data_.size(); ++k) {
        out.data_[k] = data_[k].clone();
    }
    return out;
}

Mpoly& MpolyMatrix::at(std::size_t i, std::size_t j) {
    if (i >= rows_ || j >= cols_) throw std::out_of_range("MpolyMatrix::at");
    return data_[i * cols_ + j];
}
const Mpoly& MpolyMatrix::at(std::size_t i, std::size_t j) const {
    if (i >= rows_ || j >= cols_) throw std::out_of_range("MpolyMatrix::at");
    return data_[i * cols_ + j];
}

MpolyMatrix MpolyMatrix::submatrix(const std::vector<std::size_t>& rows,
                                    const std::vector<std::size_t>& cols) const {
    MpolyMatrix out(ctx_, rows.size(), cols.size());
    for (std::size_t i = 0; i < rows.size(); ++i) {
        for (std::size_t j = 0; j < cols.size(); ++j) {
            out(i, j) = at(rows[i], cols[j]).clone();
        }
    }
    return out;
}

MpolyMatrix MpolyMatrix::minor(std::size_t row_drop, std::size_t col_drop) const {
    if (row_drop >= rows_ || col_drop >= cols_) {
        throw std::out_of_range("MpolyMatrix::minor: drop index out of range");
    }
    MpolyMatrix out(ctx_, rows_ - 1, cols_ - 1);
    std::size_t dst_i = 0;
    for (std::size_t i = 0; i < rows_; ++i) {
        if (i == row_drop) continue;
        std::size_t dst_j = 0;
        for (std::size_t j = 0; j < cols_; ++j) {
            if (j == col_drop) continue;
            out(dst_i, dst_j) = at(i, j).clone();
            ++dst_j;
        }
        ++dst_i;
    }
    return out;
}

// ---------------------------------------------------------------------------
//  Bareiss fraction-free determinant
// ---------------------------------------------------------------------------
//
//   For k = 0, ..., n-2:
//     Choose a non-zero pivot in column k below row k.
//     Swap row into position k (sign flip if so).
//     For i, j > k:  M[i,j] = (M[k,k]*M[i,j] - M[i,k]*M[k,j]) / pivot_prev
//     (Sylvester guarantees exact divisibility by pivot_prev.)
//     pivot_prev = M[k,k]
//   det = sign * M[n-1,n-1]

Mpoly MpolyMatrix::det() const {
    if (rows_ != cols_) {
        throw std::invalid_argument("MpolyMatrix::det: matrix must be square");
    }
    long n = static_cast<long>(rows_);

    if (n == 0) {
        // Empty product -> 1 (mirrors Mathematica's Det[ {} ] = 1).
        return Mpoly::one(ctx_);
    }
    if (n == 1) {
        return data_[0].clone();
    }

    MpolyMatrix M = clone();

    Mpoly pivot_prev = Mpoly::one(ctx_);
    int sign = 1;

    for (long k = 0; k + 1 < n; ++k) {
        long pivot_row = -1;
        for (long r = k; r < n; ++r) {
            if (!M(r, k).is_zero()) { pivot_row = r; break; }
        }
        if (pivot_row < 0) return Mpoly::zero(ctx_);
        if (pivot_row != k) {
            for (long c = 0; c < n; ++c) {
                std::swap(M(k, c), M(pivot_row, c));
            }
            sign = -sign;
        }

        for (long i = k + 1; i < n; ++i) {
            for (long j = k + 1; j < n; ++j) {
                Mpoly t1 = M(k, k) * M(i, j);
                Mpoly t2 = M(i, k) * M(k, j);
                Mpoly diff = t1 - t2;
                if (pivot_prev.is_one()) {
                    M(i, j) = std::move(diff);
                } else {
                    Mpoly q(ctx_);
                    if (!Mpoly::exact_divide(q, diff, pivot_prev)) {
                        throw std::runtime_error(
                            "MpolyMatrix::det: Bareiss division was not exact "
                            "(implementation bug or pathological input).");
                    }
                    M(i, j) = std::move(q);
                }
            }
            M(i, k) = Mpoly::zero(ctx_);
        }
        pivot_prev = M(k, k).clone();
    }

    Mpoly result = M(n - 1, n - 1).clone();
    if (sign < 0) result = -std::move(result);
    return result;
}

// ---------------------------------------------------------------------------
//  Adjugate  (cofactor expansion)
//
//    adj(A)[i, j] = (-1)^(i + j) * det(M_{j, i})
//
//   where M_{j, i} is the minor with row j and column i removed.
// ---------------------------------------------------------------------------

MpolyMatrix MpolyMatrix::adjugate() const {
    if (rows_ != cols_) {
        throw std::invalid_argument("MpolyMatrix::adjugate: matrix must be square");
    }
    long n = static_cast<long>(rows_);
    MpolyMatrix out(ctx_, n, n);
    for (long i = 0; i < n; ++i) {
        for (long j = 0; j < n; ++j) {
            MpolyMatrix mm = minor(j, i);
            Mpoly d = mm.det();
            if ((i + j) % 2 != 0) d = -std::move(d);
            out(i, j) = std::move(d);
        }
    }
    return out;
}

// ---------------------------------------------------------------------------

std::string MpolyMatrix::to_string() const {
    std::ostringstream oss;
    oss << "{\n";
    for (std::size_t i = 0; i < rows_; ++i) {
        oss << "  {";
        for (std::size_t j = 0; j < cols_; ++j) {
            if (j != 0) oss << ", ";
            oss << at(i, j).to_string();
        }
        oss << "}";
        if (i + 1 < rows_) oss << ",";
        oss << "\n";
    }
    oss << "}";
    return oss.str();
}

bool MpolyMatrix::operator==(const MpolyMatrix& other) const {
    if (ctx_.get() != other.ctx_.get()) return false;
    if (rows_ != other.rows_ || cols_ != other.cols_) return false;
    for (std::size_t k = 0; k < data_.size(); ++k) {
        if (!(data_[k] == other.data_[k])) return false;
    }
    return true;
}

std::ostream& operator<<(std::ostream& os, const MpolyMatrix& m) {
    return os << m.to_string();
}

// ===========================================================================
//  MfracMatrix
// ===========================================================================

MfracMatrix::MfracMatrix(std::shared_ptr<MpolyContext> ctx,
                         std::size_t rows, std::size_t cols)
    : ctx_(std::move(ctx)), rows_(rows), cols_(cols) {
    data_.reserve(rows_ * cols_);
    for (std::size_t k = 0; k < rows_ * cols_; ++k) {
        data_.emplace_back(ctx_);
    }
}

MfracMatrix MfracMatrix::from_mpoly_matrix(const MpolyMatrix& m) {
    MfracMatrix out(m.ctx(), m.rows(), m.cols());
    for (std::size_t i = 0; i < m.rows(); ++i) {
        for (std::size_t j = 0; j < m.cols(); ++j) {
            out(i, j) = Mfrac::from_mpoly(m(i, j).clone());
        }
    }
    return out;
}

MfracMatrix MfracMatrix::clone() const {
    MfracMatrix out(ctx_, rows_, cols_);
    for (std::size_t k = 0; k < data_.size(); ++k) {
        out.data_[k] = data_[k].clone();
    }
    return out;
}

std::string MfracMatrix::to_string() const {
    std::ostringstream oss;
    oss << "{\n";
    for (std::size_t i = 0; i < rows_; ++i) {
        oss << "  {";
        for (std::size_t j = 0; j < cols_; ++j) {
            if (j != 0) oss << ", ";
            oss << data_[i * cols_ + j].to_string();
        }
        oss << "}";
        if (i + 1 < rows_) oss << ",";
        oss << "\n";
    }
    oss << "}";
    return oss.str();
}

// ---------------------------------------------------------------------------
//  Inverse via cofactor / determinant.
// ---------------------------------------------------------------------------

bool mpoly_matrix_inverse(MfracMatrix& out_inverse, const MpolyMatrix& m) {
    if (m.rows() != m.cols()) {
        throw std::invalid_argument("mpoly_matrix_inverse: matrix must be square");
    }
    Mpoly d = m.det();
    if (d.is_zero()) return false;

    MpolyMatrix adj = m.adjugate();
    long n = static_cast<long>(m.rows());

    // Build (1 / det) factor as Mfrac.
    Mfrac inv_det(Mpoly::one(m.ctx()), d.clone());

    out_inverse = MfracMatrix(m.ctx(), n, n);
    for (long i = 0; i < n; ++i) {
        for (long j = 0; j < n; ++j) {
            Mfrac f = Mfrac::from_mpoly(adj(i, j).clone());
            f *= inv_det;
            out_inverse(i, j) = std::move(f);
        }
    }
    return true;
}

}  // namespace amflow::algebra
