// SPDX-License-Identifier: MIT
// ibp::reduce — top-level Kira-driven IBP reduction + diffeq.
//
//
// Mirrors AMFlow.m lines 1218-1247 plus ibp_interface/Kira/interface.m
// AnalyticReduction / DifferentialEquation (lines 454-505).
//
// The "BlackBox" abstraction (which originally existed to allow swapping
// reducers) has been flattened: this project is Kira-only per user
// constraint, so these helpers call Kira directly.

#ifndef AMFLOW_IBP_REDUCE_HPP
#define AMFLOW_IBP_REDUCE_HPP

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "amflow/algebra/mpoly.hpp"
#include "amflow/ibp/kira.hpp"
#include "amflow/ibp/libp_deriv.hpp"
#include "amflow/qft/family_config.hpp"
#include "amflow/qft/jintegral.hpp"

namespace amflow::ibp {

struct ReduceOptions {
    long ibp_rank = 5;
    long ibp_dot  = 0;
    long n_thread = 1;
    long integral_order = 5;
    std::map<std::string, std::string> numeric_values;
    std::string kira_executable;
    std::string fermat_executable;

    // Route the Reduce-mode IBP solve through FireFly finite-field
    // reconstruction instead of the Fermat algebraic backend.  See
    // KiraConfig::run_firefly.  Env-overridable in `make_kira_config`
    // via AMFLOW_RUN_FIREFLY (0/1) and AMFLOW_FIREFLY_BUNCH.
    bool run_firefly = false;
    long firefly_bunch = 0;

    // AUTO-FIREFLY ON AUX-MASS (η) DEFORMATION (default ON).
    //
    // The Fermat algebraic backend (`run_triangular` +
    // `run_back_substitution`) aborts on the auxiliary-mass-flow η-deformed
    // Reduce (SIGABRT — an internal assertion at ~35-45 GB peak in our
    // runs, or std::bad_alloc under a memory cap; SIGSEGV on some families):
    // the back-substitution swells the per-coefficient η-polynomials of
    // the eta-injected master tower until an internal invariant trips.
    // This is structural, not a tuning issue, and it is the same fault
    // that blocked a 15-propagator 3-loop self-energy reduction and a
    // 3-loop crossed-box top sector in our runs.
    //
    // The memory-light *fix* is FireFly finite-field reconstruction
    // (KiraConfig::run_firefly): it samples the IBP system over F_p at
    // numeric (η,d) points and rationally reconstructs the master
    // coefficients, never forming the dense symbolic swell.  Because the
    // η-deformed reduce is ONLY produced by AMFlow's own aux-mass DE
    // construction (the family context then carries the `eta` invariant),
    // we can detect it precisely — `fc.var_index("eta") >= 0` — and
    // default the Reduce backend to FireFly for exactly those reductions.
    //
    // Regression-safety: a NON-η family (box1, the bench suite, any
    // ordinary physical reduction) has no `eta` variable, so this is a
    // strict no-op there and the Fermat path is byte-identical.  An
    // explicit choice always wins: JSON `run_firefly:true|false`, the env
    // AMFLOW_RUN_FIREFLY (0/1), or setting `auto_firefly=false` / env
    // AMFLOW_NO_AUTO_FIREFLY=1 all suppress the auto-route.
    bool auto_firefly = true;

    // Kira 3 (arXiv:2505.20197) seeding levers; see KiraConfig.
    // env-overridable via AMFLOW_KIRA3_TRUNCATE_SP / AMFLOW_KIRA3_SECTORWISE.
    //
    // BOTH OFF by default.  The proven memory-light fix is the PLAIN
    // `run_firefly: true` finite-field reconstruction (whole solve after
    // run_initiate); it completes box1 in seconds and the 3-loop crossed-box
    // aux-mass+eta top-sector reduce at ~60 GB peak (vs the 680 GB Fermat
    // OOM-kill).  The two Kira-3 seeding levers were each MEASURED to
    // REGRESS and are therefore opt-in only:
    //
    //  * `kira3_sectorwise` -> emits `iterative_reduction: sectorwise` +
    //    `run_triangular: sectorwise` + `run_firefly: back`.  HANGS the
    //    aux-mass DE reduce even on the trivial 1-loop box1 (endless Kira
    //    "Regenerate:" sectorwise loop; 300 s timeout, no output).
    //  * `kira3_truncate_sp` -> emits a `truncate_sp {l=#loops}` seed
    //    truncation block.  It changes the master ENUMERATION, which
    //    destabilises the reduce-side UNION closure loop in reduce.cpp:
    //    the per-pass "extra master" discovery never closes, driving the
    //    dot/rank headroom to its ceiling (box1 diffeq seed ran away from
    //    r:5 to r:48,s:41).  FireFly then chokes on that enormous seeded
    //    system and times out.  Harmless-looking but couples to the
    //    closure loop; leave OFF until reduce.cpp's closure is made
    //    truncate_sp-aware.
    //
    // Re-enable per-run via AMFLOW_KIRA3_SECTORWISE=1 / AMFLOW_KIRA3_TRUNCATE_SP=1
    // only for experimentation; the production fix does not need either.
    bool kira3_truncate_sp = false;
    bool kira3_sectorwise  = false;

    // Symbolic-kinematics IBP (Layer B of `ibp_cache.hpp`): keep
    // kinematic invariants symbolic in the Kira YAML so every numeric
    // point shares ONE cached reduction.  See KiraConfig::symbolic_ibp.
    // env-overridable via AMFLOW_SYMBOLIC_IBP (0/1).
    bool symbolic_ibp = false;

    // Content-addressed Kira-results cache root (Layer A of
    // `ibp_cache.hpp`).  Empty → fall back to env `AMFLOW_IBP_CACHE`;
    // both empty → caching disabled, byte-identical to the pre-cache
    // path.  Surfaced from JSON `amf_options.blackbox.ibp_cache_dir`.
    std::string ibp_cache_dir;

    std::string work_dir;
    std::string log_file;

    // BOUNDARY-RECURSION STAGING MARKER (caps the preheat explosion on
    // region-expanded boundary sub-families).  Set ONLY by
    // AMFSystem::build_boundary for the recursive sub-family
    // BlackBoxReduce.  When set, `ibp::reduce()` honors
    // AMFLOW_BOUNDARY_PREHEAT_DOT / AMFLOW_BOUNDARY_PREHEAT_RANK as
    // CAPS (min with the jdot/jrank auto-floor) on the Masters-mode
    // preheat seed of that sub-family reduce — the target_reduce seed
    // keeps the floored values so target reducibility is untouched.
    // Default false and both envs unset ⇒ byte-identical staging.
    // Rationale: the boundary integrands (region-expanded, dots/rank up
    // to 3) auto-floor the sub-family select_mandatory_recursively to
    // (r:n_top+3, s:3, d:3); on a production 22-position family
    // (12-prop top sector 4095) that enumerates a 4,571,222-integral
    // mandatory closure -> kira std::bad_alloc at <1 G RSS
    // (measured).  The blackbox ibp_dot /
    // AMFLOW_PREHEAT_DOT knobs never reach this path (PREHEAT_* lives in
    // ibp::diffeq only; the floor overrides opts.ibp_dot upward).
    bool boundary_stage = false;

    // Upstream's `MasterRank` / `MasterDot` post-Kira master filter
    // (`Kira/interface.m:425-426,453`) is intentionally not exposed
    // here.  Upstream documents these as a manual override "only use
    // if you believe that some pseudo master integrals have appeared
    // in the list"; we pin them to the upstream default `Infinity`
    // (no filter).  If a pseudo-master problem ever surfaces, add a
    // one-line `select` post `kira_read_masters` plus a JSON knob.
};

struct ReductionContext {
    std::shared_ptr<algebra::MpolyContext> ctx;
    long d_var;
    std::vector<algebra::Mpoly> lift_gens;
};

ReductionContext make_reduction_context(const qft::FamilyConfig& fc);

struct ReduceResult {
    std::vector<qft::JIntegral>             masters;
    std::vector<std::vector<DerivTerm>>     rules;
    ReductionContext                         red_ctx;
};

ReduceResult
reduce(const qft::FamilyConfig& fc,
       const std::vector<qft::JIntegral>& targets,
       const std::vector<qft::JIntegral>& preferred,
       const std::vector<int>& top_pattern,
       const ReduceOptions& opts);

struct DiffeqResult {
    std::vector<qft::JIntegral>                              sortedmasters;
    std::vector<std::string>                                 vars;
    std::vector<std::vector<std::vector<algebra::Mfrac>>>   diffeq;
    ReductionContext                                          red_ctx;
};

DiffeqResult
diffeq(const qft::FamilyConfig& fc,
       const std::vector<qft::JIntegral>& jpreferred,
       const std::vector<std::string>& vars,
       const std::vector<int>& top_pattern,
       const ReduceOptions& opts);

// === AMFLOW_KIRA_FFSAVE (embedded FireFly save/resume persistence) ===
//
// True iff env AMFLOW_KIRA_FFSAVE is set non-empty and not "0".  When ON,
// the embedded-Kira work-dir lifecycle becomes RESUMABLE across an
// amflow-level relaunch:
//
//   * ibp::reduce / ibp::diffeq PRESERVE an existing `target_reduce`
//     working dir whenever the Kira input files they are about to write
//     (config/*.yaml, jobs.yaml[, jobs_kira2math.yaml], preferred,
//     target) byte-match what is already on disk AND the job
//     fingerprint sidecar (`.amflow_ffsave.job`: numeric_values incl.
//     the argv-only `-sd` eps, kira executable) matches.  A preserved
//     dir keeps `ff_save/` (Kira's always-on FireFly save states,
//     reductionFF.cpp `bool save = true`), `results/kira.db`, `tmp/
//     <fam>/SYSTEM_*.gz` and `sectormappings/`, so the re-fired Kira
//     RESUMES: FireFly fast-forwards past already-reconstructed
//     functions instead of restarting from probe zero.
//   * On an input MISMATCH (a different job reusing the dir) the old
//     dir is moved to `<dir>.ffsave_stale` (one slot, previous stale
//     removed) instead of being deleted, so an interrupted run's partial
//     state stays diagnosable/salvageable (harvest_ffsave_eta class).
//   * pipeline::black_box_amflow_single skips its top-level
//     `remove_all(reduce/amf)` hygiene wipes, and
//     AMFSystem ending-reduce skips its pre-wipe, so the preserved
//     dirs actually survive to the compare.
//   * ibp::kira_run appends to (rather than truncates) the per-wd
//     `kira_stdout.log`, so the streamed FireFly `Done: N / M`
//     counters of an interrupted attempt remain on disk after the relaunch.
//
// Default (env unset/0/empty): every call site takes the stock
// wipe-and-recreate path, byte-identical to the pre-knob build.
bool kira_ffsave_enabled();

}  // namespace amflow::ibp

#endif  // AMFLOW_IBP_REDUCE_HPP
