// SPDX-License-Identifier: MIT
// algebra::context_migration — strict and lossy projection between
// MpolyContexts.
//
// These two primitives consolidate inline monomial-walking code that
// would otherwise be duplicated across the pipeline / ibp / qft source
// files.
//
// Strict mode (`mfrac_to_ctx` / `mpoly_to_ctx`):
//   Every variable in the source must be present in the destination by
//   name.  Throws std::invalid_argument if a source variable is missing.
//
// Lossy mode (`project_mfrac_dropping` / `project_mpoly_dropping`):
//   Variables whose names start with any of the supplied prefixes are
//   dropped — monomials with positive exponent on such a variable are
//   discarded (set to 0 in the result).  Other missing variables still
//   throw.
//
// Default drop prefixes match the documented sentinel-prefix convention:
//   * "__amf_"   internal substitution sentinels
//   * "__feyn_"  Feynman parameters
//   * "__zsq_"   squared-mass auxiliaries

#ifndef AMFLOW_ALGEBRA_CONTEXT_MIGRATION_HPP
#define AMFLOW_ALGEBRA_CONTEXT_MIGRATION_HPP

#include <memory>
#include <string>
#include <vector>

#include "amflow/algebra/mpoly.hpp"

namespace amflow::algebra {

// ---------------------------------------------------------------------------
//  Strict context lift.  Every variable name in `src` must be present in
//  `dst`.  No-op fast path when src.ctx() and dst share identity.
// ---------------------------------------------------------------------------

Mpoly mpoly_to_ctx(const Mpoly& src,
                    const std::shared_ptr<MpolyContext>& dst);

Mfrac mfrac_to_ctx(const Mfrac& src,
                    const std::shared_ptr<MpolyContext>& dst);

// ---------------------------------------------------------------------------
//  Lossy projection.  Variables whose names start with any of the supplied
//  prefixes are dropped (monomials touching them go to zero).  Other
//  missing variables still throw.
// ---------------------------------------------------------------------------

Mpoly project_mpoly_dropping(
    const Mpoly& src,
    const std::shared_ptr<MpolyContext>& dst,
    std::vector<std::string> drop_prefixes = {"__amf_", "__feyn_", "__zsq_"});

Mfrac project_mfrac_dropping(
    const Mfrac& src,
    const std::shared_ptr<MpolyContext>& dst,
    std::vector<std::string> drop_prefixes = {"__amf_", "__feyn_", "__zsq_"});

}  // namespace amflow::algebra

#endif  // AMFLOW_ALGEBRA_CONTEXT_MIGRATION_HPP
