// SPDX-License-Identifier: MIT
// amflow::algebra::numeric_subst — substitution helpers that bake
// user-supplied `numeric_values` rationals into Mfracs symbolically.
//
// Why: C++ AMFlow holds the IBP path on a wide polynomial-ring
// context (all family variables) where MMA operates on a much
// narrower effective context.  Multivariate-GCD inside
// `fmpz_mpoly_q_canonicalise` scales super-linearly with variable
// count, so applying numeric substitution + reprojection to a narrow
// context early in the pipeline keeps polynomial sizes bounded.
//
// Mirrors MMA's `/. Numeric` pattern at `Kira/interface.m:54` (Kira
// yaml-writing), `Kira/interface.m:406` (`ComputeDerivative` Together
// per j-coefficient), and `AMFlow.m:817` (`ReduceBoundary` per
// inner reduction).

#ifndef AMFLOW_ALGEBRA_NUMERIC_SUBST_HPP
#define AMFLOW_ALGEBRA_NUMERIC_SUBST_HPP

#include <map>
#include <memory>
#include <set>
#include <string>

#include "amflow/algebra/mpoly.hpp"

namespace amflow::algebra {

// Test if variable index `v` has positive exponent in any monomial of `p`.
bool mpoly_has_var(const Mpoly& p, long v);

// Parse a numeric-string ("3", "-5", "1/2", "0.1") into an fmpq_t.
// Returns true on success.  fmpq_t must already be initialised by the caller.
bool parse_rational_string(const std::string& s, fmpq_t v);

// Substitute numeric values for every variable in `src.ctx()` that is
// NOT in `keep_names`.  The returned Mfrac is still on `src.ctx()`
// (same MpolyContext) but the substituted variables now have exponent
// 0 in every monomial — they can then be safely projected to a
// narrower context via `mfrac_to_ctx_lenient`.
//
// Throws if a non-keep variable has a positive exponent somewhere AND
// is missing from `numeric_values` (caller bug).
Mfrac substitute_fc_vars(
    const Mfrac& src,
    const std::map<std::string, std::string>& numeric_values,
    const std::set<std::string>& keep_names);

// Project `src` from its current context onto `dst_ctx` by variable
// name.  Variables that exist in src.ctx() but NOT in dst_ctx are
// silently dropped IF their exponent is 0 in every monomial; throws
// with a caller-bug message if any monomial has a non-zero exponent on
// such a variable.  (Caller should `substitute_fc_vars` first to
// collapse fc-only variables to numeric values before projecting.)
Mfrac mfrac_to_ctx_lenient(const Mfrac& src,
                            const std::shared_ptr<MpolyContext>& dst_ctx);

// Composition helper: substitute + project in one call.  Equivalent to
// `mfrac_to_ctx_lenient(substitute_fc_vars(src, numeric_values, keep_names), dst_ctx)`.
// `keep_names` defaults to the set of dst_ctx's variable names.
Mfrac substitute_and_narrow(
    const Mfrac& src,
    const std::map<std::string, std::string>& numeric_values,
    const std::shared_ptr<MpolyContext>& dst_ctx);

}  // namespace amflow::algebra

#endif  // AMFLOW_ALGEBRA_NUMERIC_SUBST_HPP
