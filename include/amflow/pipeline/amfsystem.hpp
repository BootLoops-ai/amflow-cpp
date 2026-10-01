// SPDX-License-Identifier: MIT
// Layer 16: AMF system pipeline.
//
// Mirrors upstream AMFlow.m, lines 833-1204
// (https://gitlab.com/multiloop-pku/amflow ; line numbers track the
// upstream master branch at the time of the v1.0.0 release).
//
// Conceptual model (re-derived from the Mathematica source -- v2)
// ------------------------------------------------------------------
//
//   An `AMFSystem` is a node in the AMF tree.  It holds:
//
//   * `fc`               — the family this node operates on; in
//                          Mathematica it's the AMFlowInfo state.
//   * `preferred`        — the J integrals we want solved.
//   * `etac`             — per-propagator eta-coefficient.  All zero
//                          means this is an *ending* system.
//
//   Pipeline:
//
//   * `setup()`:
//       - if ending: nothing to do; solve() looks up Vacuum / explicit BC.
//       - else:
//           a. inject eta into propagators (fc -> fc_with_eta);
//           b. compute differential equation dM/dη on fc_with_eta via
//              Layer 15d's BlackBoxDiffeq;
//           c. find_all_region(fc_with_eta) and filter zero regions;
//           d. for each non-zero region:
//                 i. compute powers (Mfrac in fc + ε per master);
//                ii. compute boundary_integrands then boundary_integrals
//                    -> a list of (boundary_family_prop, integrand_terms);
//               iii. for each boundary family (sub-FamilyConfig):
//                       * create a sub-FamilyConfig that inherits
//                         loops/legs/conservation/replacement from fc
//                         but uses the boundary family's prop list;
//                       * collect every J appearing in this family's
//                         integrand_terms -> `all_ints`;
//                       * call BlackBoxReduce(all_ints, {}) on the
//                         sub-fc (Kira picks its own masters);
//                       * register a child AMFSystem with sub_masters
//                         as its preferred list, RECURSIVELY (it will
//                         re-run amf_system_setup_master and decide
//                         whether to inject or end).
//                       * record per-(master, integrand_term) "table[j,k]"
//                         coefficients in the sub_master basis -- this
//                         is what we evaluate at solve-time.
//
//   * `solve(epslist)`:  bottom-up.  For each ε:
//       - if ending: look up Vacuum[L, n] (or explicit_boundary).
//       - else: build the BC by walking each region's table:
//             BC[master_i] = list of (mu(ε), value(ε))
//             value(ε) = Σ_term coef_at_eps(ε) * Σ_subm  table[term, subm]
//                                                       * sub_solution[subm](ε)
//         then call Layer 9 amflow(de_at_eps, BC).

#ifndef AMFLOW_PIPELINE_AMFSYSTEM_HPP
#define AMFLOW_PIPELINE_AMFSYSTEM_HPP

#include <map>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "amflow/numeric/acb_wrap.hpp"
#include "amflow/qft/amfmode.hpp"
#include "amflow/ibp/reduce.hpp"
#include "amflow/qft/family_config.hpp"
#include "amflow/qft/jintegral.hpp"
#include "amflow/ibp/kira.hpp"

namespace amflow::pipeline {

enum class EndingScheme {
    Tradition,
    Cutkosky,
    SingleMass,
    // Always-applicable fallback: when none of the named schemes
    // produces an η-injection / Cutkosky / single-mass setup, the
    // dispatcher falls back to Trivial — set up an AMFSystem with no
    // η injection.  Mirrors upstream `AMFSystemSetupMaster[..., "Trivial"]`
    // (AMFlow.m:1086-1095) and the auto-append at AMFlow.m:1034.
    Trivial,
};
const char* ending_scheme_name(EndingScheme s) noexcept;

struct AMFSystemOptions {
    // Mathematica default $AMFMode = {"Prescription", "Mass",
    // "Propagator"}.  The fallback order matters: for the 1-loop bubble,
    // Prescription yields nothing, Mass selects both propagators, and
    // only then should Propagator be considered.
    std::vector<qft::AMFMode> amf_modes = {
        qft::AMFMode::Prescription,
        qft::AMFMode::Mass,
        qft::AMFMode::Propagator
    };

    // Mathematica default $EndingScheme = {"Tradition", "Cutkosky",
    // "SingleMass"}.  Recursion typically terminates at Tradition
    // (when sub-system preferred masters are vacuum integrals -- top
    // sector is empty).  Cutkosky/SingleMass only kick in for
    // phase-space integrals or single-mass families respectively.
    std::vector<EndingScheme> ending_schemes = {EndingScheme::Tradition};

    ibp::ReduceOptions bb;
    std::string cache_root;
    std::string direction;

    // ε-grid parallelism (Lever 1).  After setup() has built the diffeq /
    // boundary tables (all Kira/IBP work is done or cache-hit), every
    // ε-sample's DE solve is independent.  When eps_parallel > 1 the
    // per-ε `solve_one_eps` loop at each non-ending node is dispatched
    // to a thread pool of this size.  Default 1 = serial (bit-identical
    // to the pre-parallel code path).  Also overridable via env
    // `AMFLOW_EPS_PARALLEL`.
    long eps_parallel = 1;

    // Node-tree parallelism (Lever 2).  A non-ending node's sub-system
    // solves are mutually independent — each writes only its own
    // `solutions_` (and its own per-node checkpoint key).  When
    // node_parallel > 1 the `sub->solve(epslist)` loop at each
    // non-ending node dispatches the SIBLING sub-trees to a thread pool
    // of this size, with the same thread-local options-snapshot
    // inheritance as Lever 1.  Default 1 = serial (bit-identical to the
    // pre-parallel code path).  Also overridable via env
    // `AMFLOW_NODE_PARALLEL`.  See `tree_stats()` for the sibling-width
    // diagnostic that quantifies the headroom.
    //
    // Thread budget: Lever 1's "only one node runs its ε-loop at any
    // instant" invariant is relaxed when this knob is on — up to
    // node_parallel sub-trees run concurrently at EACH non-ending depth,
    // and each non-ending node still spawns its own ε pool, so size the
    // product (roughly node_parallel^depth × eps_parallel at the deepest
    // concurrent frontier) against the core count.  Ending sub-systems
    // may invoke Kira at solve time; that is safe here because all
    // reduce/cache coordination is file-based and already multi-process
    // hardened (per-key flock stampede lock + write-then-rename in the
    // ibp_cache), each ending reduce runs in a work dir keyed by
    // (system_id, integral), and the Kira launcher forks-and-execs
    // immediately with the chdir in the child — no process-global
    // mutation in the parent.
    long node_parallel = 1;

    // Safety net: maximum recursion depth in setup().  Default 8.
    // The Tradition scheme alone is not sufficient to terminate
    // the recursion for typical multi-loop families -- SingleMass
    // (FactorizeFamily) is required.  Until SingleMass is fully
    // implemented in Layer 16, deep families will throw a clear
    // error pointing to AUDIT.md instead of hanging.
    long max_recursion_depth = 8;
    long current_depth = 0;

    // GRAVITYFLOW_FIXEDPOINT guard: canonical propagator signature of the
    // immediate parent family (sorted propagators-after-conservation),
    // recorded by build_boundary when it creates a boundary sub-system,
    // together with the parent's top-sector width (next field).  setup()
    // compares the child's own signature against it and raises a
    // structured GRAVITYFLOW_FIXEDPOINT error when they are identical
    // AND the child shows no top-sector descent — the eta->inf boundary
    // recursion is then a fixed point (b_{k+1} == b_k) and would grind
    // to max_recursion_depth without ever producing boundary data.
    // Empty = top system / no parent tracking.
    std::string parent_prop_signature;

    // Companion to parent_prop_signature: the parent system's top-sector
    // width = number of active slots in qft::get_top_sector(preferred)
    // (propagator positions carrying a positive index in any preferred
    // master).  An identical propagator multiset alone does NOT prove
    // self-similarity: a boundary child may legitimately re-use the
    // parent's propagator basis while its masters live in a strictly
    // smaller top sector — the automatic-loop boundary of e.g. the
    // massless box is exactly this (the child's masters are sub-sector
    // bubble/tadpole integrals of the same 4-propagator family, and the
    // descent terminates because the width strictly decreases level by
    // level, a well-founded measure on the naturals).  The guard
    // therefore fires only when the signatures match AND the child's
    // width has NOT strictly decreased — no certified progress, i.e.
    // the recursion regenerates a system at least as large: the true
    // b_{k+1} == b_k fixed point (equal width at depth 1), or any
    // same-propagator chain the moment its width stops shrinking.
    // -1 = no parent tracking; a recorded signature arriving without a
    // width fires conservatively (pre-descent behavior).
    long parent_top_sector_count = -1;

    // Optional explicit boundary values for ending systems whose
    // masters are not in the Vacuum table.  Key = JIntegral key
    // ("<family>|<i1>|<i2>|...").  Value = per-ε vector of AcbValue
    // (wrapped in shared_ptr because AcbValue is move-only).
    std::map<std::string, std::shared_ptr<std::vector<numeric::AcbValue>>>
        explicit_boundary;
};

struct AMFSystemSolution {
    std::vector<numeric::AcbValue> master_values;   // one per preferred[i]

    // For SingleMass-style systems, the actual J integral that the
    // user originally requested may differ from `preferred[i]`.  We
    // record the "global preferred" mapping here so callers can
    // reconstruct user-requested values.
    //
    //   global_preferred[i] = the original J integral requested
    //   global_value[i]     = prefactor[i] * master_values[jpre_index[i]]
    //
    // When unused (Tradition scheme), these are empty and the caller
    // should use master_values directly.
    std::vector<qft::JIntegral>  global_preferred;
    std::vector<numeric::AcbValue>   global_values;
};

// One per-region BoundaryTerm coefficient.  Lives on the sub-system's
// reduction context (= sub_fc.ctx + 'd').  At solve-time we substitute
// invariants and ε to obtain a numeric value.
struct BoundaryTerm {
    algebra::Mfrac coef;            // Mfrac on sub-system's red_ctx
    long  sub_master_idx;  // index into sub-system's preferred[]
    std::string sub_master_key;  // stable lookup key for reordered outputs
};

// One reduced boundary entry, corresponding to one Mathematica
// {powers, config, masters, table} item from ReduceBoundary.  A single
// region may yield multiple such boundary entries/families; they are stored
// as separate RegionBoundary values so solve_one_eps can concatenate their
// contributions exactly like Mathematica's Join@@@Transpose[bc].
struct RegionBoundary {
    // mu(ε) for each preferred master (Mfrac on a "powers ctx" =
    // fc.ctx + ε, but actually only contains ε after substitution).
    // We keep them as Mfracs; substitute at solve-time.
    std::vector<algebra::Mfrac> powers;

    // Indices into AMFSystem::subsystems_ holding the root systems for
    // this region's boundary family.  Most regions have exactly one
    // root, but SingleMass factorization can legitimately return
    // multiple roots that must be combined multiplicatively (mirrors
    // Mathematica's AMFSystemsSolution / AMFSystemCombineSolution).
    std::vector<long> sub_system_indices;

    // For each preferred master i, a list of orders.
    // integrand_terms[i][k] = sum_terms list for the k-th order
    // boundary integrand of master i.
    std::vector<std::vector<std::vector<BoundaryTerm>>> integrand_terms;
};

class AMFSystem {
public:
    AMFSystem(qft::FamilyConfig fc,
                std::vector<qft::JIntegral> preferred,
                std::vector<int>       etac,
                EndingScheme           ending_used,
                const AMFSystemOptions& opts);

    // SingleMass-style multi-component system construction.  Each
    // component produces one AMFSystem rooted on a `(sub_fc, sub_preferred,
    // etac, gamma_params, original_preferred, original_jpre_idx)` tuple.
    //
    // gamma_params[i] = (n, m0, m_eps_coef, loop_num)
    //
    //   prefactor[i] = (-1)^(-n)
    //                * Gamma(2 - eps - m0 - m_eps_coef * eps)
    //                * Gamma(-2 + eps + m0 + m_eps_coef * eps + n)
    //                / (Gamma(2 - eps) * Gamma(n))
    //
    // (The `m = total - n - (loop_num - 1) * (2 - eps)` from .m
    //  expands to `m0 + m_eps_coef * eps` with
    //  m0 = total - n - 2*(loop_num-1), m_eps_coef = (loop_num-1).)
    struct GammaParams {
        long n;
        long m0;
        long m_eps_coef;
        long loop_num;
    };
    struct SingleMassConfig {
        std::vector<qft::JIntegral>   original_preferred;
        std::vector<GammaParams> gamma_params;     // size = original_preferred.size()
        std::vector<long>        jpre_idx;          // for each original master, index into preferred
    };
    void set_single_mass_config(SingleMassConfig cfg) { sm_cfg_ = std::move(cfg); }

    struct CutkoskyConfig {
        std::vector<qft::JIntegral> original_preferred;
        long phase_loop_num = 0;
    };
    void set_cutkosky_config(CutkoskyConfig cfg) { ck_cfg_ = std::move(cfg); }

    AMFSystem(const AMFSystem&) = delete;
    AMFSystem& operator=(const AMFSystem&) = delete;
    AMFSystem(AMFSystem&&) noexcept;
    AMFSystem& operator=(AMFSystem&&) noexcept;
    ~AMFSystem();

    // Recursively build sub-systems, diffeq, and boundary tables.
    void setup();

    // Bottom-up solve at every ε.
    void solve(const std::vector<numeric::AcbValue>& epslist);

    // ---- accessors ----
    const qft::FamilyConfig&            family() const noexcept { return *fc_; }
    const std::vector<qft::JIntegral>&  preferred() const noexcept { return preferred_; }
    const std::vector<int>&        etac() const noexcept { return etac_; }
    EndingScheme                   ending_scheme() const noexcept { return ending_; }
    bool                           is_ending() const noexcept { return is_ending_; }
    long                           system_id() const noexcept { return system_id_; }
    const AMFSystemOptions&        options()   const noexcept { return opts_; }

    const std::vector<AMFSystemSolution>& solutions() const noexcept {
        return solutions_;
    }

    long sub_system_count() const noexcept { return (long)subsystems_.size(); }
    const AMFSystem& sub_system(long i) const { return *subsystems_[(std::size_t)i]; }

    // ---- Fixed-point probe (η-flow self-similar boundary) ----
    //
    // For families whose boundary recursion is a FIXED POINT (the
    // η→∞ boundary family is the family itself — e.g. cut-eikonal PM
    // families), the recursion never
    // terminates and AMFlow has no independent boundary data.  But the
    // η=0 solution is LINEAR in the boundary sub-family value vector v:
    //     M(0) = A(ε) · v ,    and self-similarity gives  v ≡ M(0)
    // (the sub-family masters are the family's own masters), so M is the
    // eigenvector of A(ε) at eigenvalue 1, up to scale — pinned by one
    // analytic anchor.  When env AMFLOW_FIXEDPOINT_PROBE is set (to the
    // output JSON path) and this system is at depth 0, build_boundary
    // records the sub-family master basis (fp_sub_masters_) and SKIPS
    // the recursive sub-system setup; fixedpoint_probe() then assembles
    // A(ε) column-by-column via unit-vector injections into
    // solve_one_eps and writes the JSON table.
    const std::vector<qft::JIntegral>& fp_sub_masters() const noexcept {
        return fp_sub_masters_;
    }
    void fixedpoint_probe(const std::vector<numeric::AcbValue>& epslist,
                          const std::string& out_path);

    // Node-tree parallelism diagnostic: walk the sub-system tree.  `total` = node count
    // including this one; `leaves` = ending-system count; `max_siblings`
    // = the widest fan-out at any single node (the headroom for
    // node-level parallelism within one ε).
    void tree_stats(long& total, long& leaves, long& max_siblings) const;

private:
    std::unique_ptr<qft::FamilyConfig> fc_;
    std::vector<qft::JIntegral>       preferred_;
    std::vector<int>             etac_;
    EndingScheme                 ending_;
    bool                          is_ending_;
    AMFSystemOptions              opts_;
    long                          system_id_ = 0;

    // Per-system path direction, computed during setup() from the
    // prescriptions of η-touching loops (mirrors upstream
    // `AMFSystemDirection`, AMFlow.m:981-991).  Used in solve() to
    // override the global `numeric::run_direction()` when this
    // system's η contour orientation differs from the default.
    numeric::RunningOptions::Direction direction_
        = numeric::RunningOptions::Direction::NegIm;

    // Eta-injected version of fc_ (used for diffeq + boundary).  Only
    // built if !is_ending_.
    std::unique_ptr<qft::FamilyConfig> fc_with_eta_;

    ibp::DiffeqResult         diffeq_result_;
    std::vector<long>            pref_to_sorted_;     // preferred[i] -> sortedmasters index

    std::vector<RegionBoundary>  regions_;
    std::vector<std::unique_ptr<AMFSystem>> subsystems_;

    // Boundary "pattern groups": per-group, a list of mu_i (Mfrac on the
    // powers context) — one mu per preferred master.  Mirrors the MMA
    // variable `pattern = Select[pattern, AnyTrue[...]]` in
    // AMFSystemSolution.  At solve time each (pattern_group, master)
    // contributes a trivial  (mu_i -> 0)  entry to bc_sorted, supplying
    // the asymptotic orders that are NOT covered by a non-trivial
    // boundary integrand.  Without this, amflow sees an empty BC for any
    // sub-system whose regions all have border=-1 and returns 0.
    std::vector<std::vector<algebra::Mfrac>> bc_pattern_;

    std::vector<AMFSystemSolution> solutions_;

    // Optional scheme-specific global-preferred remapping.
    std::optional<SingleMassConfig> sm_cfg_;
    std::optional<CutkoskyConfig>   ck_cfg_;

    // Fixed-point probe state (see public fixedpoint_probe()).
    std::vector<qft::JIntegral>           fp_sub_masters_;
    const std::vector<numeric::AcbValue>* fp_injected_ = nullptr;
    bool                                  fp_probe_active_ = false;

    // Internal helpers.
    void build_diffeq();
    void build_boundary();
    AMFSystemSolution solve_one_eps(std::size_t eps_index, const numeric::AcbValue& eps);

    // Env-gated (AMFLOW_BOUNDARY_DUMP=<dir>) read-only JSON dump of the
    // eta=infinity boundary tables, written once per non-ending system at
    // the end of setup().  No-op when the env var is unset/empty; write
    // failures warn and continue (never abort a solve).
    void dump_boundary_json() const;

    // fp-probe column batching (AMFLOW_FP_PROBE_BATCH=1): solve
    // ALL unit-column injections of one eps point through ONE batched
    // fundamental-solution transport (ode::amflow_multi).  Returns
    // [col][preferred_i] — bit-identical to looping solve_one_eps over
    // fp_injected_ unit vectors.  Only valid on the probe path
    // (fp_probe_active_, all regions without sub-systems).
    std::vector<std::vector<numeric::AcbValue>>
    fp_probe_solve_eps_batched(std::size_t eps_index, const numeric::AcbValue& eps);
};

// Per-scheme ending check (mirrors AMFSystemEndingQ).
bool ending_q(const qft::FamilyConfig& fc,
              const std::vector<qft::JIntegral>& preferred,
              EndingScheme scheme,
              const AMFSystemOptions& opts);
bool ending_q_all(const qft::FamilyConfig& fc,
                  const std::vector<qft::JIntegral>& preferred,
                  const AMFSystemOptions& opts);

// Top-level: given a preferred-master list, build the root AMFSystem
// (or systems, for SingleMass which produces multiple roots).
std::vector<std::unique_ptr<AMFSystem>>
amf_system_setup_master(const qft::FamilyConfig& fc,
                          const std::vector<qft::JIntegral>& preferred,
                          const AMFSystemOptions& opts);

// Build the entire sub-tree (calls setup() on every node).
std::vector<std::unique_ptr<AMFSystem>>
amf_systems_setup(const qft::FamilyConfig& fc,
                  const std::vector<qft::JIntegral>& preferred,
                  const AMFSystemOptions& opts);

// Solve at every ε in epslist, return per-ε master values for the
// root system's preferred[].
std::vector<AMFSystemSolution>
amf_systems_solution(std::vector<std::unique_ptr<AMFSystem>>& systems,
                     const std::vector<numeric::AcbValue>& epslist);

}  // namespace amflow::pipeline

#endif  // AMFLOW_PIPELINE_AMFSYSTEM_HPP
