// SPDX-License-Identifier: MIT
// ode::zero — solving near eta = 0 (Layer 7 calcx00 / calc_zero).
//
//
// Mirrors MMA DESolver.m:
//   Calcx00                             line 856
//   LearnFromRuleS / LearnFromRuleSAll  lines 920-922
//   CalcZero                            line 952
//   PickZeroRuleS                       line 1091
//
// Pre-conditions
// --------------
//
//   `mat` is a *rational* matrix.  Its only allowed singularity at eta=0
//   may be of any Poincare rank; normalize_mat brings it to normalised
//   Fuchsian form (rank 0 + Jordan-form leading matrix + integer-shifted
//   eigenvalues).  After that calc_zero recovers the asymptotic expansion
//   of every integral near eta=0 starting from a numerical boundary at a
//   chosen regular point x0.

#ifndef AMFLOW_ODE_ZERO_HPP
#define AMFLOW_ODE_ZERO_HPP

#include <cstddef>
#include <utility>
#include <vector>

#include "amflow/numeric/acb_wrap.hpp"
#include "amflow/numeric/options.hpp"
#include "amflow/numeric/rational.hpp"
#include "amflow/ode/asy.hpp"
#include "amflow/ode/blocks.hpp"
#include "amflow/ode/normalize.hpp"

namespace amflow::ode {

// ---------------------------------------------------------------------------
//  Calcx00
// ---------------------------------------------------------------------------
//
//  Inputs (mirrors the .m signature exactly):
//     nheq      : block equations from nh_equations(B, EquationMode::Singular)
//     nheqn     : numeric variant of nheq
//     bc        : per-integral numerical boundary value at x0  (length = N)
//     x0        : regular point at which `bc` was evaluated
//     behavior  : per-block list of (mu, log_power) -- output of asymptotic_behavior
//
//  Returns: per-integral asymptotic expansion near eta=0.

std::vector<AsyExpansion>
calcx00(const std::vector<BlockEquation>&    nheq,
        const std::vector<BlockEquationNum>& nheqn,
        const std::vector<numeric::AcbValue>& bc,
        acb_srcptr                            x0,
        const AsymptoticBehaviorList&         behavior,
        long                                  prec = numeric::working_prec_bits());

// ---------------------------------------------------------------------------
//  LearnFromRuleS / LearnFromRuleSAll
// ---------------------------------------------------------------------------

BlockBehavior learn_from_rule_s(const AsyExpansion& asy);

AsymptoticBehaviorList
learn_from_rule_s_all(const std::vector<AsyExpansion>& asy_list,
                      const std::vector<std::vector<std::size_t>>& blocks);

// ---------------------------------------------------------------------------
//  CalcZero  (top level)
// ---------------------------------------------------------------------------

std::vector<AsyExpansion>
calc_zero(const numeric::RationalMatrix& de,
          const std::vector<numeric::AcbValue>& bc,
          acb_srcptr x0,
          long prec = numeric::working_prec_bits());

// ---------------------------------------------------------------------------
//  calc_zero_multi  (fp-probe column batching)
// ---------------------------------------------------------------------------
//
//  Multi-column CalcZero for the AMFLOW_FP_PROBE_BATCH=1 fixed-point-probe
//  path.  The column-independent (bc-independent) work — normalize_mat
//  (incl. ToFuchsian / kernel-rescue), the invT(x0) evaluation, the
//  nh_equations(+num) of the normalized B, asymptotic_behavior and the
//  T power-series map — is done ONCE and shared; per column only the bcT
//  rotation, the (optional) learn pass, calcx00 and the final ps_map run.
//  Per-column arithmetic is the exact operation sequence of calc_zero, so
//  results are bit-identical to per-column calc_zero calls.
//
//  Returns result[col][integral].
std::vector<std::vector<AsyExpansion>>
calc_zero_multi(const numeric::RationalMatrix& de,
                const std::vector<std::vector<numeric::AcbValue>>& bc_cols,
                acb_srcptr x0,
                long prec = numeric::working_prec_bits());

// ---------------------------------------------------------------------------
//  PickZeroRuleS
// ---------------------------------------------------------------------------
//
//  From a single-integral asymptotic expansion (after calc_zero), pick the
//  numerical leading value at eta = 0:
//    * If no integer mu is present: return 0.
//    * If integer mu = 0 is present: return (log^0, eta^0) coefficient.
//    * If the only integer mu is k > 0: return 0.
//    * If the only integer mu is k < 0: return (log^0, eta^{-k}) coefficient.
//    * If multiple integer mus are present: error.
numeric::AcbValue
pick_zero_rule_s(const AsyExpansion& asy,
                 long prec = numeric::working_prec_bits());

}  // namespace amflow::ode

#endif  // AMFLOW_ODE_ZERO_HPP
