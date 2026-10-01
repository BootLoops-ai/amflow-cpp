// SPDX-License-Identifier: MIT
// numeric::log — env-var-gated diagnostic trace gate.
//
//
// All compile-time-disabled diagnostic stderr in this project is wrapped in
//
//     AMFLOW_TRACE("AMFLOW_DEBUG_FOO") {
//         std::cerr << "[FOO] ..." << ... << "\n";
//     }
//
// which fires only when the named environment variable is set non-empty.
//
// Trace categories (set the env var non-empty to enable diagnostic stderr):
//
//   AMFLOW_DEBUG_BC       — boundary-condition assembly (qft regions, amflow build/solve)
//   AMFLOW_DEBUG_DBO      — DetermineBoundaryOrder internals (ode)
//   AMFLOW_DEBUG_CT_PERM  — CalcTaylor row permutation (ode)
//   AMFLOW_DEBUG_REGULAR  — regular-point recurrence (ode)
//   AMFLOW_DEBUG_STAGES   — top-level stage labels (ode)
//   AMFLOW_TRACE_LAYER6   — fine-grained layer 6 trace (ode)
//   AMFLOW_DEBUG_PRECISE  — full-precision arb_get_str (numeric)
//   AMFLOW_DEBUG_SCHEME   — EndingScheme dispatch (amflow Cutkosky/SingleMass firing)
//
// Behavior-modifier env vars (NOT diagnostic; these change semantics —
// listed for completeness, not routed through trace_enabled):
//
//   AMFLOW_NO_SPARSE_CHOP
//   AMFLOW_SPARSE_CHOP_DIGITS
//   AMFLOW_JORDAN_DUMP
//   AMFLOW_KIRA_DUMP

#ifndef AMFLOW_NUMERIC_LOG_HPP
#define AMFLOW_NUMERIC_LOG_HPP

namespace amflow::numeric::log {

// True iff env var `env_name` is set to a non-empty string.
bool trace_enabled(const char* env_name);

}  // namespace amflow::numeric::log

#define AMFLOW_TRACE(env_name) \
    if (::amflow::numeric::log::trace_enabled(env_name))

#endif  // AMFLOW_NUMERIC_LOG_HPP
