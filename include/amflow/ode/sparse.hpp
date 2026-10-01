// SPDX-License-Identifier: MIT
// ode::sparse — sparse linear algebra for the recurrence-coefficient solver
// (Layer 3).
//
//
// Mirrors MMA DESolver.m lines 589-680.
//
// Sparse rows are stored as `vector<(col, value)>` kept in ascending column
// order (matching the original `col -> value` rule list semantics).
//
// Variable indexing convention (verbatim from .m):
//
//   For an input matrix produced by ConstructMatrix with block size N
//   and totalorder T, the column count is (T + 2) * N.  Columns are
//   organised in (T + 2) blocks of N columns each.  Block index `j`
//   (0..T+1) corresponds to the series coefficient f_{T+1-j} of the
//   unknown vector.

#ifndef AMFLOW_ODE_SPARSE_HPP
#define AMFLOW_ODE_SPARSE_HPP

#include <cstddef>
#include <iosfwd>
#include <utility>
#include <vector>

#include "amflow/numeric/acb_wrap.hpp"
#include "amflow/numeric/options.hpp"
#include "amflow/ode/blocks.hpp"     // BlockEquationNum used by callers

namespace amflow::ode {

// ---------------------------------------------------------------------------
//  SparseEntry / SparseRow / SparseSystem
// ---------------------------------------------------------------------------

struct SparseEntry {
    long                col;
    numeric::AcbValue   val;

    SparseEntry() : col(0) {}
    SparseEntry(long c, numeric::AcbValue v) : col(c), val(std::move(v)) {}

    SparseEntry(SparseEntry&&) noexcept            = default;
    SparseEntry& operator=(SparseEntry&&) noexcept = default;

    SparseEntry(const SparseEntry& other) : col(other.col), val(other.val.clone()) {}
    SparseEntry& operator=(const SparseEntry& other) {
        if (this != &other) { col = other.col; val = other.val.clone(); }
        return *this;
    }
};

using SparseRow    = std::vector<SparseEntry>;
using SparseSystem = std::vector<SparseRow>;

// Build a row from a dense vector.  Drops exact zeros.  Keeps ascending col.
SparseRow sparsify_row(const std::vector<numeric::AcbValue>& dense);

// Build a system from a dense matrix.
SparseSystem sparsify_matrix(const std::vector<std::vector<numeric::AcbValue>>& dense);

// Shift every column index by `n`  (== SparseTranslate).
void sparse_translate(SparseRow& row, long n);
void sparse_translate(SparseSystem& sys, long n);

// In-place scalar multiply.  `s` is given as raw acb_srcptr because it is
// usually a function of the calling state (e.g. -leading_coeff).
void scalar_sparse(SparseRow& row, acb_srcptr s,
                   long prec = numeric::working_prec_bits());

// out := a + b   (PlusSparse).  Both rows must be sorted by col; the result
// is sorted, with explicit zeros removed.
SparseRow plus_sparse(const SparseRow& a,
                      const SparseRow& b,
                      long             prec = numeric::working_prec_bits());

// Split a system into rows that contain column m and rows that do not.
struct SplitResult {
    std::vector<std::size_t> dep;
    std::vector<std::size_t> indep;
};
SplitResult split_system(const SparseSystem& sys, long m);

// ---------------------------------------------------------------------------
//  ConstructedMatrix
// ---------------------------------------------------------------------------
//
//  Output of ConstructMatrix.
//
//  total_columns = (total_order + 2) * block_size
//
//  block_systems[i] holds the i-th "block row": all rows that involve the
//  variables in coefficient block i.

struct ConstructedMatrix {
    long                       block_size    = 0;
    long                       total_order   = 0;
    long                       total_columns = 0;
    std::vector<SparseSystem>  block_systems;
};

// Encode  dx(eta) f'(eta) - ax(eta) f(eta) = g(eta)  in the truncated
// coefficient basis.  `dxexp` is dx as a flat coefficient list; `axexp[i][j][k]`
// is the η^k coefficient of ax(i,j).
ConstructedMatrix construct_matrix(const std::vector<numeric::AcbValue>& dxexp,
                                   const std::vector<std::vector<std::vector<numeric::AcbValue>>>& axexp,
                                   long total_order,
                                   long prec = numeric::working_prec_bits());

// ---------------------------------------------------------------------------
//  ForwardSparseGaussian / SparseGaussian
// ---------------------------------------------------------------------------
//
//  ForwardSparseGaussian (rows, vset)  ->  (reduced, residue)
//
//  For each variable index n in vset:
//    1. Find rows that contain column n  ("dep"); others go to "indep".
//    2. If no dep, skip.
//    3. Take the first dep row, normalise its leading-n entry to 1.
//    4. Use it to eliminate n from every other dep row.
//    5. Move the normalised row to `reduced`; residual rows go back into the
//       system.

struct GaussianRowSet {
    std::vector<SparseRow> reduced;
    SparseSystem           residue;
};

GaussianRowSet forward_sparse_gaussian(SparseSystem        sys,
                                       const std::vector<long>& vset,
                                       long prec = numeric::working_prec_bits());

// Top-level SparseGaussian over the output of ConstructMatrix and the
// non-homogeneous list (one entry per row of every block_system, flat).
//
//   nh layout: length = block_size * total_count_of_block_rows
//
// Each non-zero nh entry is appended as an extra column at total_columns.
struct SparseGaussianResult {
    std::vector<SparseRow> reduce;     // reduced rows in variable order
    SparseSystem           residue;    // residue from the last block
};

SparseGaussianResult sparse_gaussian(const ConstructedMatrix& cm,
                                     const std::vector<numeric::AcbValue>& nh,
                                     long prec = numeric::working_prec_bits());

// Chop values in a sparse row using AMFLOW_SPARSE_CHOP_DIGITS env var
// (defaults to numeric::chop_pre()).  Set AMFLOW_NO_SPARSE_CHOP=1 to
// suppress all chopping.  Modifies `row` in place.
void chop_sparse(SparseRow& row);
void chop_sparse(SparseSystem& sys);

// Diagnostic.
std::ostream& operator<<(std::ostream& os, const SparseRow& row);

}  // namespace amflow::ode

#endif  // AMFLOW_ODE_SPARSE_HPP
