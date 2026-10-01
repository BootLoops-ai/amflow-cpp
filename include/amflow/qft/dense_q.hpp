// SPDX-License-Identifier: MIT
// qft::dense_q — dense fmpq matrix utilities used by region/topology.

#ifndef AMFLOW_QFT_DENSE_Q_HPP
#define AMFLOW_QFT_DENSE_Q_HPP

#include <cstddef>
#include <vector>

#include <flint/fmpq.h>

namespace amflow::qft {

class DenseQMatrix {
public:
    DenseQMatrix(std::size_t rows, std::size_t cols);
    ~DenseQMatrix();

    DenseQMatrix(const DenseQMatrix&) = delete;
    DenseQMatrix& operator=(const DenseQMatrix&) = delete;

    DenseQMatrix(DenseQMatrix&& other) noexcept;
    DenseQMatrix& operator=(DenseQMatrix&& other) noexcept;

    std::size_t rows() const noexcept { return rows_; }
    std::size_t cols() const noexcept { return cols_; }

    fmpq*       at(std::size_t r, std::size_t c)        { return data_[r * cols_ + c]; }
    const fmpq* at(std::size_t r, std::size_t c) const  { return data_[r * cols_ + c]; }

    void set_si(std::size_t r, std::size_t c, long v);
    void set_fmpq(std::size_t r, std::size_t c, const fmpq_t v);

    void rref();
    bool has_unit_vector_at_last_columns(std::size_t n_lhs) const;

private:
    std::size_t rows_;
    std::size_t cols_;
    std::vector<fmpq*> data_;
};

std::vector<std::size_t>
maximal_group_rows(std::size_t rows, std::size_t cols,
                    const std::vector<fmpq*>& flat);

}  // namespace amflow::qft

#endif  // AMFLOW_QFT_DENSE_Q_HPP
