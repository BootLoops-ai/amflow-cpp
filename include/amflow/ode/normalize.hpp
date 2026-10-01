// SPDX-License-Identifier: MIT
// ode::normalize — NormalizeMat + AsymptoticBehavior + Sylvester solver
// (Layer 7 normalisation pipeline).
//
//
// Mirrors MMA DESolver.m lines 251-515.
//
// The pipeline produces (T, invT, B) such that
//
//     B = invT . mat . T - invT . D[T, eta]
//
// (the standard similarity-with-derivative for the substitution I = T J in
// dI/deta = mat I), where B is in normalised Fuchsian form near eta = 0:
//
//     * Poincare rank 0 ;
//     * eigenvalues of (eta * B)|_{eta=0} have fractional parts in their
//       canonical (region) representative ;
//     * leading matrix is in Jordan form.

#ifndef AMFLOW_ODE_NORMALIZE_HPP
#define AMFLOW_ODE_NORMALIZE_HPP

#include <cstddef>
#include <utility>
#include <vector>

#include <flint/acb_mat.h>
#include <flint/fmpq_mat.h>

#include "amflow/numeric/acb_wrap.hpp"
#include "amflow/numeric/options.hpp"
#include "amflow/numeric/rational.hpp"

namespace amflow::ode {

// ---------------------------------------------------------------------------
//  AsymptoticBehavior
// ---------------------------------------------------------------------------
//
//  Per .m line 251-263.  Returns one entry per block; each entry is a list
//  of (mu, log_power) pairs giving the maximal allowed log-power for each
//  leading exponent detected by the diagonal-block analysis on the *Jordan-form
//  leading matrix*.

struct BehaviorEntry {
    numeric::AcbValue mu;
    long              log_power = 0;
};

using BlockBehavior          = std::vector<BehaviorEntry>;
using AsymptoticBehaviorList = std::vector<BlockBehavior>;

AsymptoticBehaviorList
asymptotic_behavior(const numeric::RationalMatrix& mat,
                    long prec = numeric::working_prec_bits());

// ---------------------------------------------------------------------------
//  NormalizeMat
// ---------------------------------------------------------------------------

struct NormalizationResult {
    numeric::RationalMatrix T;
    numeric::RationalMatrix invT;
    numeric::RationalMatrix B;
};

NormalizationResult
normalize_mat(const numeric::RationalMatrix& mat,
              long prec = numeric::working_prec_bits());

// Internal Layer-7 helper for CalcZero's algebraic-Jordan fallback.
struct BlockJordanRotation {
    std::vector<std::size_t>       block;
    std::vector<numeric::AcbValue> u;      // row-major constant matrix
    std::vector<numeric::AcbValue> invu;
    std::vector<numeric::AcbValue> jor;
};

struct CalcZeroNormalization {
    NormalizationResult              base;
    std::vector<BlockJordanRotation> rotations;
};

CalcZeroNormalization
normalize_mat_for_calc_zero(const numeric::RationalMatrix& mat,
                            long prec = numeric::working_prec_bits());

// ---------------------------------------------------------------------------
//  Sylvester-style solve used by SolveOffDiagonal
// ---------------------------------------------------------------------------
//
//   Given:    a0 (m x m), c0 (n x n), b0 (n x m), p (scalar)
//   Solve:    c0 . X - X . a0 + p X = -b0      for X (n x m)
//
// Returns false if the linear system is singular at the working precision.

namespace internal {

bool solve_off_diagonal_acb(acb_mat_t       X,
                            const acb_mat_t a0,
                            const acb_mat_t b0,
                            const acb_mat_t c0,
                            long            p,
                            long prec = numeric::working_prec_bits());

bool solve_off_diagonal_fmpq(fmpq_mat_t        X,
                             const fmpq_mat_t  a0,
                             const fmpq_mat_t  b0,
                             const fmpq_mat_t  c0,
                             long              p);

}  // namespace internal

}  // namespace amflow::ode

#endif  // AMFLOW_ODE_NORMALIZE_HPP
