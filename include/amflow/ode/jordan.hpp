// SPDX-License-Identifier: MIT
// ode::jordan — exact Jordan decomposition over Q (Layer 7 helper).
//
//
// Mirrors Mathematica's JordanDecomposition[A] for rational matrices.  Used
// by NormalizeMat where the residue stays over Q and routing through acb
// eigenvectors would introduce fragile rationalization.

#ifndef AMFLOW_ODE_JORDAN_HPP
#define AMFLOW_ODE_JORDAN_HPP

#include <vector>

#include <flint/fmpq.h>
#include <flint/fmpq_mat.h>

#include "amflow/numeric/rational.hpp"

namespace amflow::ode {

// Output-parameter API because fmpq_mat_t is a C array typedef.  Caller
// pre-initializes out_S, out_J, out_Sinv to the dimensions of A.
//
// Throws if A is not square, the characteristic polynomial has a non-linear
// irreducible factor over Q, or an internal exact linear-algebra invariant
// fails.
void jordan_decomposition_exact(fmpq_mat_t         out_S,
                                fmpq_mat_t         out_J,
                                fmpq_mat_t         out_Sinv,
                                std::vector<long>& out_block_sizes,
                                const fmpq_mat_t   A);

// Extract the constant term of eta^shift * mat at eta = 0.  After the shift
// each entry must be regular at eta = 0.
void fmpq_mat_from_rational_residue(fmpq_mat_t                       out,
                                    const numeric::RationalMatrix&   mat,
                                    long                              shift);

// Convert an fmpq matrix into a RationalMatrix with constant entries.
numeric::RationalMatrix
rational_matrix_from_fmpq_mat(const fmpq_mat_t in);

}  // namespace amflow::ode

#endif  // AMFLOW_ODE_JORDAN_HPP
