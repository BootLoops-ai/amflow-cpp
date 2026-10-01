// SPDX-License-Identifier: MIT
// ode::inf — solving from the asymptotic boundary at eta = infinity (Layer 6).
//
//
// Mirrors MMA DESolver.m lines 684-794:
//
//   BuildTaylor                   line 687
//   DetermineBlockBoundaryOrder   line 692
//   DetermineBoundaryOrder        line 720
//   ReverseBCS                    line 723
//   UnionBCS                      line 724
//   ReadBCS                       line 727
//   CalcTaylor                    line 738
//   CalcInf                       line 782
//
// What this layer does
// --------------------
//
//   At eta = infinity AMFlow expresses every integral as a list of rules
//
//        bc[i] = { mu_1 -> v_1,  mu_2 -> v_2,  ... }
//
//   meaning  I_i(eta) = sum_k v_k * eta^{mu_k}  near infinity (one number
//   per asymptotic order, no log terms, no higher-order coefficients).
//
//   To turn this into a numerical solution at a chosen regular point x0
//   AMFlow:
//     1. flips eta -> x = 1/eta, so the singularity moves to x = 0.  This
//        rewrites the matrix:
//             deinf(x) = - x^{-2} * de(1/x)
//        and flips the boundary keys: ReverseBCS / UnionBCS.
//     2. groups the boundary keys by "differ-by-integer", choosing one
//        leading exponent ini[i] per region per integral.
//     3. peels off the leading exponent via the substitution
//             J(x) = x^{-ini} * I(x)
//        so that J satisfies a *Taylor*-form ODE around x = 0.
//     4. Sparse-Gauss-eliminates the truncated coefficient system and
//        plugs the boundary values into the unsolved positions.  This
//        gives the truncated Taylor series of every integral near x = 0.
//     5. Tags the series with its mu and merges the regions back into a
//        single AsyExpansion per integral.
//
//   The output of calc_inf is `vector<AsyExpansion>` -- one expansion (a
//   list of (mu, [[coefficients]])  rules) per integral.  The expansion
//   variable is x = 1/eta.

#ifndef AMFLOW_ODE_INF_HPP
#define AMFLOW_ODE_INF_HPP

#include <cstddef>
#include <utility>
#include <vector>

#include "amflow/numeric/acb_wrap.hpp"
#include "amflow/numeric/options.hpp"
#include "amflow/numeric/rational.hpp"
#include "amflow/ode/asy.hpp"
#include "amflow/ode/blocks.hpp"
#include "amflow/ode/regular.hpp"     // for TaylorCoefficients

namespace amflow::ode {

// ---------------------------------------------------------------------------
//  BoundarySpec
// ---------------------------------------------------------------------------
//
//  One boundary at eta = infinity for a single integral, expressed as a
//  list of (mu, value) pairs.  C++ analogue of the MMA list
//  { mu_1 -> v_1, mu_2 -> v_2, ... }.

struct BoundaryEntry {
    numeric::AcbValue mu;
    numeric::AcbValue value;
};

using BoundarySpec = std::vector<BoundaryEntry>;

// Deep-copy helper.
BoundarySpec clone_boundary(const BoundarySpec& bc);
std::vector<BoundarySpec> clone_boundaries(const std::vector<BoundarySpec>& bcs);

// ReverseBCS:  (mu, v) -> (-mu, v).  Used to translate the boundary from the
// eta-frame to the x = 1/eta-frame.
BoundarySpec reverse_bcs(const BoundarySpec& bc);
std::vector<BoundarySpec> reverse_bcs(const std::vector<BoundarySpec>& bcs);

// UnionBCS:  combine entries that share the same mu (exact equality on the
// midpoint, with chop tolerance ChopPre).  Sums the values.
BoundarySpec union_bcs(const BoundarySpec& bc,
                       long prec = numeric::working_prec_bits());
std::vector<BoundarySpec> union_bcs(const std::vector<BoundarySpec>& bcs,
                                    long prec = numeric::working_prec_bits());

// ReadBCS:  in `bcs` find every entry whose (mu - region) is an integer,
// return ini = the smallest such mu and the entries with mu replaced by
// (mu - ini), packaged as (integer_offset, value).
//
// If no entry is in the region, returns (ini = region, empty list).
struct ReadBcsResult {
    numeric::AcbValue                              ini;
    std::vector<std::pair<long, numeric::AcbValue>> shifted;
};
ReadBcsResult read_bcs(const BoundarySpec& bc,
                       acb_srcptr region,
                       long prec = numeric::working_prec_bits());

// ---------------------------------------------------------------------------
//  Boundary-order determination
// ---------------------------------------------------------------------------
//
//  Per integral i, returns the maximum "order" (= position in the truncated
//  series at which the system leaves a degree of freedom) plus 1.  -1 means
//  the integral has no free coefficient that the user must supply.

std::vector<long>
determine_boundary_order(const numeric::RationalMatrix& mat,
                         const std::vector<long>& power_per_integral);

std::vector<long>
determine_block_boundary_order(const numeric::RationalMatrix& block_mat,
                               const std::vector<long>& power);

// Rational-power overloads (needed when the pattern's μ value is not an
// integer, e.g. `1 - eps` evaluated at eps = 1/100 -> 99/100).
std::vector<long>
determine_boundary_order(const numeric::RationalMatrix& mat,
                         const std::vector<numeric::RationalFunction>& power_q);

std::vector<long>
determine_block_boundary_order(const numeric::RationalMatrix& block_mat,
                               const std::vector<numeric::RationalFunction>& power_q);

// ---------------------------------------------------------------------------
//  CalcTaylor
// ---------------------------------------------------------------------------
//
//  Inputs:
//     mat : the *post-BuildTaylor* matrix for a single region.  System ODE
//           is dJ/dx = mat * J near x = 0 with J Taylor-analytic at x = 0.
//     bc  : per-integral list of (order, value) pairs supplying boundary
//           values f0[i, order] for the unsolved coefficients.
//
//  Returns: per-integral truncated Taylor series of length (XOrder + 1).

using TaylorRegionInput =
    std::vector<std::vector<std::pair<long, numeric::AcbValue>>>;

TaylorCoefficients calc_taylor(const numeric::RationalMatrix& mat,
                               const TaylorRegionInput& bc,
                               long prec = numeric::working_prec_bits());

// Variant that takes the deinf matrix together with `ini` per integral.
// Internally calls build_taylor_symbolic on `de_inf` with the given `ini`,
// then dispatches to the main calc_taylor.
TaylorCoefficients
calc_taylor_with_ini(const numeric::RationalMatrix& de_inf,
                     const std::vector<numeric::AcbValue>& ini_per_integral,
                     const TaylorRegionInput& bc,
                     long prec = numeric::working_prec_bits());

// ---------------------------------------------------------------------------
//  CalcInf  (top-level for this layer)
// ---------------------------------------------------------------------------
//
//  Inputs:
//     de  : original rational ODE matrix in eta.
//     bcs : per-integral list of {mu -> value} boundary rules at eta=infty.
//
//  Returns: per-integral asymptotic expansion in x = 1/eta near x = 0.
//
//  IMPORTANT: The returned AsyExpansion has its mu values *already* in the
//  x-frame: a term  (mu_x, coeffs)  represents  x^{mu_x} * sum_n coeffs[0][n] x^n
//  with no log entries (boundary input has no log structure).

// ---------------------------------------------------------------------------
//  calc_inf_multi  (fp-probe column batching)
// ---------------------------------------------------------------------------
//
//  Multi-column CalcInf for the AMFLOW_FP_PROBE_BATCH=1 fixed-point-probe
//  path.  Every column shares the same de and the same BoundarySpec mu
//  structure (only entry values differ), so the column-independent work —
//  make_deinf, region detection, the canonical permutation,
//  build_taylor_symbolic, nh_equations(+num) of the pure Taylor matrix and
//  the per-block ConstructMatrix — is done ONCE per region and reused for
//  every column.  The per-column remainder (nh assembly, SparseGaussian,
//  back-substitution) is the exact operation sequence of calc_inf, so
//  results are bit-identical to per-column calc_inf calls.
//
//  Returns result[col][integral].
std::vector<std::vector<AsyExpansion>>
calc_inf_multi(const numeric::RationalMatrix& de,
               const std::vector<std::vector<BoundarySpec>>& bcs_cols,
               long prec = numeric::working_prec_bits());

std::vector<AsyExpansion>
calc_inf(const numeric::RationalMatrix& de,
         const std::vector<BoundarySpec>& bcs,
         long prec = numeric::working_prec_bits());

}  // namespace amflow::ode

#endif  // AMFLOW_ODE_INF_HPP
