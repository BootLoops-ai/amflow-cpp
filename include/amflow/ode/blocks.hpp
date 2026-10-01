// SPDX-License-Identifier: MIT
// ode::blocks — equation-structure analysis (Layer 2).
//
//
// Mirrors the "equation analysis" subsection of MMA DESolver.m (lines
// 189-304):
//
//   AnalyzeBlockOld / AnalyzeBlock / AnalyzeBlock0 / SubBlockID
//   PoincareRank
//   NHEquations / NHEquationsNum
//   ToNum / PickElement / PickList / PickMat
//
// Block indices in this header are 0-based (C++ convention) even though
// MMA is 1-based.  Translation: subtract 1 reading, add 1 emitting human
// diagnostics.

#ifndef AMFLOW_ODE_BLOCKS_HPP
#define AMFLOW_ODE_BLOCKS_HPP

#include <cstddef>
#include <functional>
#include <vector>

#include "amflow/numeric/acb_wrap.hpp"
#include "amflow/numeric/options.hpp"
#include "amflow/numeric/rational.hpp"

namespace amflow::ode {

// ---------------------------------------------------------------------------
//  BlockPartition  -- result of AnalyzeBlock + SubBlockID.
// ---------------------------------------------------------------------------
//
//  blocks() : ordered list of contiguous "logical blocks", each given by the
//             *original* row indices it contains, in ascending order.  The
//             block ordering itself is the topological order: block i only
//             depends on blocks of strictly smaller index.
//
//  lower_dependencies(i) : list of block indices j < i such that block i
//                          actually has a non-zero entry in mat[block_i, block_j].
//                          Equals SubBlockID[mat][[i+1]] in MMA.
//
//  sub_rows(i)           : flat list of original row indices that appear in
//                          *any* block referenced by lower_dependencies(i).
//                          Equivalent to .m line  sub = Join @@ blocks[[#]] & /@ sub
//                          inside NHEquations.

class BlockPartition {
public:
    BlockPartition() = default;

    std::size_t                        size()   const { return blocks_.size(); }
    bool                               empty()  const { return blocks_.empty(); }

    const std::vector<std::size_t>&    block(std::size_t i) const { return blocks_.at(i); }
    const std::vector<std::vector<std::size_t>>& blocks() const { return blocks_; }

    const std::vector<std::size_t>&    lower_dependencies(std::size_t i) const {
        return sub_block_ids_.at(i);
    }

    // Resolve lower_dependencies into a flattened, ordered row-index list.
    std::vector<std::size_t> sub_rows(std::size_t i) const;

    // Allow analyser functions to populate us directly.
    void set_blocks(std::vector<std::vector<std::size_t>> b)         { blocks_         = std::move(b); }
    void set_sub_block_ids(std::vector<std::vector<std::size_t>> s)  { sub_block_ids_  = std::move(s); }

private:
    std::vector<std::vector<std::size_t>> blocks_;
    std::vector<std::vector<std::size_t>> sub_block_ids_;
};

// ---------------------------------------------------------------------------
//  AnalyzeBlock variants
// ---------------------------------------------------------------------------

// "Old" algorithm (.m line 193-201): contiguous-row "right-extent" partition.
// Best when the input is already block-lower-triangular and you only want to
// recover the diagonal block sizes.
std::vector<std::vector<std::size_t>>
analyze_block_old(const numeric::RationalMatrix& mat);

// Operates on any 0/1 sparsity pattern (e.g. extracted from a numeric Jordan
// block in AsymptoticBehavior).
std::vector<std::vector<std::size_t>>
analyze_block_old_pattern(std::size_t n,
                          std::function<bool(std::size_t, std::size_t)> is_nonzero);

// "Modern" algorithm (.m line 204-237): closure-based partition.  Used
// everywhere in DESolver.m to discover the block structure of a rational
// matrix.
std::vector<std::vector<std::size_t>>
analyze_block(const numeric::RationalMatrix& mat);

// Build a full BlockPartition (blocks + lower-dependencies).
BlockPartition build_partition(const numeric::RationalMatrix& mat);

// SubBlockID[mat] - just the lower_dependencies vector.
std::vector<std::vector<std::size_t>>
sub_block_ids(const numeric::RationalMatrix& mat,
              const std::vector<std::vector<std::size_t>>& blocks);

// ---------------------------------------------------------------------------
//  PoincareRank
// ---------------------------------------------------------------------------
//
//   PoincareRank[mat] = Max@Exponent[mat//Denominator, eta, Min] - 1
//
//   The denominator's eta-valuation is the order of the pole at eta=0; we
//   take the max across all entries and subtract 1.  Result -1 means the
//   matrix is regular at eta=0; 0 means Fuchsian-ready.

long poincare_rank(const numeric::RationalMatrix& mat);

// ---------------------------------------------------------------------------
//  NHEquations / NHEquationsNum
// ---------------------------------------------------------------------------

enum class EquationMode {
    Singular,   // factor = eta
    Regular,    // factor = 1
    Taylor      // factor = eta^(PoincareRank+1)
};

// Symbolic non-homogeneous block equation:
//
//     dx(eta) * f'(eta) = ax(eta) * f(eta) + bx(eta) * g(eta)
//
//   * dx :  fmpq polynomial (LCM-of-denominators step)
//   * ax :  square fmpq polynomial matrix (size = block.size())
//   * bx :  rectangular rational-function matrix (rows = block.size(),
//           cols = sub.size())
//   * block, sub : index lists, identical to BlockPartition's contents.
struct BlockEquation {
    numeric::FmpqPoly                                      dx;
    std::vector<std::vector<numeric::FmpqPoly>>            ax;
    std::vector<std::vector<numeric::RationalFunction>>    bx;
    std::vector<std::size_t>                                block;
    std::vector<std::size_t>                                sub;
};

// Numeric variant: every polynomial replaced by its acb coefficient list.
//
//   dxexp[k]            = coefficient of eta^k in dx
//   axexp[i][j][k]      = coefficient of eta^k in ax(i,j)
//   bxexpn[i][j][k]     = coefficient of eta^k in numerator(bx(i,j))
//   bxexpd[i][j][k]     = coefficient of eta^k in denominator(bx(i,j))
struct BlockEquationNum {
    std::vector<numeric::AcbValue>                                              dxexp;
    std::vector<std::vector<std::vector<numeric::AcbValue>>>                    axexp;
    std::vector<std::vector<std::vector<numeric::AcbValue>>>                    bxexpn;
    std::vector<std::vector<std::vector<numeric::AcbValue>>>                    bxexpd;
    std::vector<std::size_t>                                                    block;
    std::vector<std::size_t>                                                    sub;
};

// Build the symbolic form.  Returns one BlockEquation per block, in the
// partition's topological order.
std::vector<BlockEquation>
nh_equations(const numeric::RationalMatrix& mat, EquationMode mode);

// Variant taking a pre-computed partition (so callers that already analysed
// the matrix do not pay for it twice).
std::vector<BlockEquation>
nh_equations(const numeric::RationalMatrix& mat,
             const BlockPartition& partition,
             EquationMode mode);

// Convert a list of BlockEquations to numeric form.
std::vector<BlockEquationNum>
nh_equations_num(const std::vector<BlockEquation>& eqs,
                 long prec = numeric::working_prec_bits());

// ---------------------------------------------------------------------------
//  Coefficient utilities  (PickElement / PickList / PickMat, ToNum)
// ---------------------------------------------------------------------------

// CoefficientList[poly, eta]: copies coefficients into a fresh AcbValue
// vector at the given precision.  Length = degree+1, or 0 for the zero
// polynomial.
std::vector<numeric::AcbValue>
to_num(const numeric::FmpqPoly& poly, long prec = numeric::working_prec_bits());

// PickElement[list, order] : list[order] if it exists, else 0.
numeric::AcbValue
pick_element(const std::vector<numeric::AcbValue>& list, std::size_t order);

// PickList: PickElement for each list inside `lists`.
std::vector<numeric::AcbValue>
pick_list(const std::vector<std::vector<numeric::AcbValue>>& lists,
          std::size_t order);

// PickMat: PickList per row.
std::vector<std::vector<numeric::AcbValue>>
pick_mat(const std::vector<std::vector<std::vector<numeric::AcbValue>>>& mat,
         std::size_t order);

}  // namespace amflow::ode

#endif  // AMFLOW_ODE_BLOCKS_HPP
