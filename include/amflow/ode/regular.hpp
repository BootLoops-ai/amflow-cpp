// SPDX-License-Identifier: MIT
// ode::regular — regular-point recurrence and contour propagation (Layer 5).
//
//
// Mirrors MMA DESolver.m lines 798-849:
//
//   ExpandNHEquationsNum
//   Calcx1x2  (Taylor expansion at one regular point, returns coefficients
//              up to XOrder)
//   CalcRun   (chain Calcx1x2 across a list of regular points)
//
// What "regular point" means here
// -------------------------------
//
//   The matrix mat(eta) has no pole at eta = x0 -- equivalently, every
//   denominator is non-zero there.  This is what nh_equations(.., Regular)
//   prepares: it never multiplies by an explicit eta factor, so dx(0) is
//   identically the polynomial that vanishes at the singularities, evaluated
//   at the working point.

#ifndef AMFLOW_ODE_REGULAR_HPP
#define AMFLOW_ODE_REGULAR_HPP

#include <cstddef>
#include <vector>

#include "amflow/numeric/acb_wrap.hpp"
#include "amflow/numeric/options.hpp"
#include "amflow/numeric/rational.hpp"
#include "amflow/ode/blocks.hpp"

namespace amflow::ode {

// ---------------------------------------------------------------------------
//  expand_nh_equations_num
// ---------------------------------------------------------------------------
//
//  Takes a numeric block-equation list assembled at eta = 0 and re-expands
//  every polynomial as a Taylor series in (eta - x0), truncated to
//  (XOrder + 1) coefficients.  The block / sub index lists are unchanged.
//
//  When x0 == 0 the result is identical to the input.

std::vector<BlockEquationNum>
expand_nh_equations_num(const std::vector<BlockEquationNum>& nheqn,
                        acb_srcptr x0,
                        long prec = numeric::working_prec_bits());

// Convenience overload: x0 = 0 -> just returns a clone-ish copy.
std::vector<BlockEquationNum>
expand_nh_equations_num(const std::vector<BlockEquationNum>& nheqn,
                        long prec = numeric::working_prec_bits());

// ---------------------------------------------------------------------------
//  calcx1x2  (Taylor solve at a single regular point)
// ---------------------------------------------------------------------------
//
//  Inputs:
//    nheqn : block equations expanded at the *original* eta=0 (the function
//            re-expands them at x0 internally).
//    bc    : current values of all integrals at eta = x0  (one AcbValue per
//            integral, length = total number of integrals).
//    x0    : Taylor-expansion centre.
//
//  Output: result[i] is a vector<AcbValue> of length (XOrder + 1).
//          result[i][k] is the coefficient of (eta - x0)^k in integral i.

using TaylorCoefficients = std::vector<std::vector<numeric::AcbValue>>;

TaylorCoefficients
calcx1x2(const std::vector<BlockEquationNum>& nheqn,
         const std::vector<numeric::AcbValue>& bc,
         acb_srcptr x0,
         long prec = numeric::working_prec_bits());

// ---------------------------------------------------------------------------
//  calc_run  (propagate boundary values through a chain of regular points)
// ---------------------------------------------------------------------------
//
//  Propagates the boundary `bc` through the list run = [x_0, x_1, ..., x_M]
//  of regular points.  Returns the value of every integral at the *last*
//  point x_M.
//
//  Pre-condition: bc are the values at run[0]; the matrix `de` is regular
//  on every chord [run[i], run[i+1]].  Path planner's job (Layer 8).

std::vector<numeric::AcbValue>
calc_run(const numeric::RationalMatrix& de,
         const std::vector<numeric::AcbValue>& bc,
         const std::vector<numeric::AcbValue>& run,
         long prec = numeric::working_prec_bits());

// Variant taking a pre-computed nheqn (avoids redoing the symbolic block
// analysis on every step).
std::vector<numeric::AcbValue>
calc_run(const std::vector<BlockEquationNum>& nheqn,
         const std::vector<numeric::AcbValue>& bc,
         const std::vector<numeric::AcbValue>& run,
         long prec = numeric::working_prec_bits());

// ---------------------------------------------------------------------------
//  calcx1x2_multi / calc_run_multi  (fp-probe column batching)
// ---------------------------------------------------------------------------
//
//  Multi-column variants used by the AMFLOW_FP_PROBE_BATCH=1 fixed-point
//  probe path: all columns share the SAME de/nheqn (hence the same
//  expand_nh_equations_num re-expansion and the same mata[m] pick-matrices
//  per step); only the propagated value vectors differ.  Per-column
//  arithmetic is the exact operation sequence of the single-column
//  functions, so results are bit-identical to running each column through
//  calcx1x2 / calc_run separately.

std::vector<TaylorCoefficients>
calcx1x2_multi(const std::vector<BlockEquationNum>& nheqn,
               const std::vector<std::vector<numeric::AcbValue>>& bc_cols,
               acb_srcptr x0,
               long prec = numeric::working_prec_bits());

// Propagates every column through the run chain.  Returns [col][integral]
// values at the last run point.
std::vector<std::vector<numeric::AcbValue>>
calc_run_multi(const std::vector<BlockEquationNum>& nheqn,
               const std::vector<std::vector<numeric::AcbValue>>& bc_cols,
               const std::vector<numeric::AcbValue>& run,
               long prec = numeric::working_prec_bits());

// ---------------------------------------------------------------------------
//  evaluate_taylor
// ---------------------------------------------------------------------------
//
//  Given coeffs[i][k] = (eta - x0)^k coefficient of integral i, return
//      sum_k coeffs[i][k] * dh^k    for each i,
//  using midpoint-only Horner to mirror Mathematica's point arithmetic
//  (avoids catastrophic interval blow-up on cancellation-heavy contours).

std::vector<numeric::AcbValue>
evaluate_taylor(const TaylorCoefficients& coeffs,
                acb_srcptr dh,
                long prec = numeric::working_prec_bits());

}  // namespace amflow::ode

#endif  // AMFLOW_ODE_REGULAR_HPP
