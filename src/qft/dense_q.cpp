// SPDX-License-Identifier: MIT
// qft::dense_q — implementation.
//

#include "amflow/qft/dense_q.hpp"

#include <algorithm>
#include <stdexcept>
#include <utility>

#include <flint/flint.h>

namespace amflow::qft {

DenseQMatrix::DenseQMatrix(std::size_t rows, std::size_t cols)
    : rows_(rows), cols_(cols), data_(rows * cols, nullptr) {
    for (std::size_t k = 0; k < rows * cols; ++k) {
        fmpq* p = (fmpq*)flint_malloc(sizeof(fmpq));
        fmpq_init(p);
        data_[k] = p;
    }
}

DenseQMatrix::~DenseQMatrix() {
    for (auto* p : data_) {
        if (p) {
            fmpq_clear(p);
            flint_free(p);
        }
    }
}

DenseQMatrix::DenseQMatrix(DenseQMatrix&& other) noexcept
    : rows_(other.rows_), cols_(other.cols_),
      data_(std::move(other.data_)) {
    other.rows_ = 0;
    other.cols_ = 0;
    other.data_.clear();
}

DenseQMatrix& DenseQMatrix::operator=(DenseQMatrix&& other) noexcept {
    if (this != &other) {
        for (auto* p : data_) {
            if (p) { fmpq_clear(p); flint_free(p); }
        }
        rows_  = other.rows_;
        cols_  = other.cols_;
        data_  = std::move(other.data_);
        other.rows_ = 0;
        other.cols_ = 0;
        other.data_.clear();
    }
    return *this;
}

void DenseQMatrix::set_si(std::size_t r, std::size_t c, long v) {
    fmpq_set_si(at(r, c), v, 1);
}

void DenseQMatrix::set_fmpq(std::size_t r, std::size_t c, const fmpq_t v) {
    fmpq_set(at(r, c), v);
}

void DenseQMatrix::rref() {
    std::size_t pivot_row = 0;
    for (std::size_t col = 0; col < cols_ && pivot_row < rows_; ++col) {
        std::size_t found = rows_;
        for (std::size_t r = pivot_row; r < rows_; ++r) {
            if (!fmpq_is_zero(at(r, col))) { found = r; break; }
        }
        if (found == rows_) continue;
        if (found != pivot_row) {
            for (std::size_t c = 0; c < cols_; ++c) {
                std::swap(data_[pivot_row * cols_ + c],
                          data_[found * cols_ + c]);
            }
        }

        fmpq_t pivot;
        fmpq_init(pivot);
        fmpq_set(pivot, at(pivot_row, col));
        for (std::size_t c = col; c < cols_; ++c) {
            fmpq_div(at(pivot_row, c), at(pivot_row, c), pivot);
        }

        for (std::size_t r = 0; r < rows_; ++r) {
            if (r == pivot_row) continue;
            if (fmpq_is_zero(at(r, col))) continue;
            fmpq_t mult;
            fmpq_init(mult);
            fmpq_set(mult, at(r, col));
            for (std::size_t c = col; c < cols_; ++c) {
                fmpq_t prod;
                fmpq_init(prod);
                fmpq_mul(prod, at(pivot_row, c), mult);
                fmpq_sub(at(r, c), at(r, c), prod);
                fmpq_clear(prod);
            }
            fmpq_clear(mult);
        }
        fmpq_clear(pivot);
        ++pivot_row;
    }
}

bool DenseQMatrix::has_unit_vector_at_last_columns(std::size_t n_lhs) const {
    for (std::size_t r = 0; r < rows_; ++r) {
        bool zero_lhs = true;
        for (std::size_t c = 0; c < n_lhs; ++c) {
            if (!fmpq_is_zero(at(r, c))) { zero_lhs = false; break; }
        }
        if (!zero_lhs) continue;
        for (std::size_t c = n_lhs; c < cols_; ++c) {
            if (!fmpq_is_zero(at(r, c))) return true;
        }
    }
    return false;
}

std::vector<std::size_t>
maximal_group_rows(std::size_t rows, std::size_t cols,
                    const std::vector<fmpq*>& flat) {
    if (flat.size() != rows * cols) {
        throw std::invalid_argument(
            "maximal_group_rows: flat array size mismatch");
    }
    if (rows == 0) return {};

    DenseQMatrix T(cols, rows);
    for (std::size_t r = 0; r < rows; ++r) {
        for (std::size_t c = 0; c < cols; ++c) {
            fmpq_set(T.at(c, r), flat[r * cols + c]);
        }
    }
    T.rref();

    std::vector<std::size_t> picks;
    for (std::size_t r = 0; r < cols; ++r) {
        for (std::size_t c = 0; c < rows; ++c) {
            const fmpq* v = T.at(r, c);
            if (!fmpq_is_zero(v)) {
                if (fmpq_is_one(v)) picks.push_back(c);
                break;
            }
        }
    }
    std::sort(picks.begin(), picks.end());
    picks.erase(std::unique(picks.begin(), picks.end()), picks.end());
    return picks;
}

}  // namespace amflow::qft
