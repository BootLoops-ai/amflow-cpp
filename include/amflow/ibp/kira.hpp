// SPDX-License-Identifier: MIT
// ibp::kira — Kira IBP-reduction interface.
//
//
// Mirrors AMFlow.m / ibp_interface/Kira/interface.m.
//
// Workflow
// --------
//   1. Caller constructs a `KiraConfig` from a FamilyConfig + JIntegral
//      preferred masters + IBP rank/dot/topsector + numeric values for
//      kinematic invariants.
//   2. yaml writers populate <dir>/config/{integralfamilies,kinematics}.yaml,
//      <dir>/jobs.yaml, <dir>/preferred, <dir>/target.
//   3. kira_run invokes the kira binary.
//   4. parsers read masters and target table.

#ifndef AMFLOW_IBP_KIRA_HPP
#define AMFLOW_IBP_KIRA_HPP

#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "amflow/algebra/mpoly.hpp"
#include "amflow/qft/family_config.hpp"
#include "amflow/qft/jintegral.hpp"

namespace amflow::ibp {

struct KiraConfig {
    const qft::FamilyConfig* fc = nullptr;

    std::vector<int> top_pattern;

    long ibp_rank = 5;
    long ibp_dot  = 0;
    long n_thread = 1;
    long integral_order = 5;

    // When true the Reduce-mode jobs.yaml emits `run_firefly: true`
    // (FireFly finite-field reconstruction) in place of the Fermat
    // algebraic `run_triangular`/`run_back_substitution` pair.  FireFly
    // evaluates the IBP system over F_p at numeric points for every
    // symbolic variable (here `eta` + the dimension `d`), solves in
    // finite-field arithmetic (tiny RAM), and reconstructs the rational
    // master coefficients over Q from many evaluations.  This is the
    // memory-light fix for the auxiliary-mass+eta expression swell that
    // OOM-kills the symbolic Fermat path on the 3-loop crossed-box top
    // sector (the eta-power master tower).  Masters-mode enumeration is
    // unaffected (it never solves the dense system).
    bool run_firefly = false;
    // FireFly bunch size (`--bunch_size`); 0 = leave Kira default.
    long firefly_bunch = 0;

    // Per-subsector DOT-CAP on the Reduce-mode seed block.  When true the
    // reduce job emits `d: ibp_dot` so every k-prop sub-sector is seeded
    // at no more than `ibp_dot` dots (matching the Masters-mode
    // select_mandatory_recursively).  Without this the reduce seeds at
    // r = n_top + ibp_dot in EVERY sub-sector, exploding the mandatory
    // list (one 3-loop eta top sector: 499 -> 260) and the Fermat/FireFly
    // system size.  The dot-cap is the LOAD-BEARING half of the verified
    // eta-deformed SIG6 fix (dot-cap + lower seed + FireFly TOGETHER);
    // FireFly alone stalls on the uncapped 499-list (the documented
    // aux-mass-DE Regenerate stall).  Defaulted ON for eta families in
    // make_kira_config (see reduce.cpp); env AMFLOW_REDUCE_DCAP forces
    // on/off.  Reduce-mode only (Masters-mode is left byte-identical so
    // the canonical basis and ibp_cache key are unchanged).
    bool reduce_dcap = false;

    // --- Kira 3 (arXiv:2505.20197, Lange-Usovitsch-Wu) seeding levers ---
    //
    // The auxiliary-mass DE derivative targets carry dotted / negative-index
    // (numerator, i.e. ISP) integrals.  Under Kira-2-style flat (r,s) seeding
    // these promote a large spurious sub-sector master tower (the the 3-loop crossed-box
    // crossed-box top sector ran away to ~473 enumerated masters and OOM-killed
    // at 680 GB).  Kira 3 caps this at the source:
    //
    //  * `truncate_sp: {topologies, l}` -- the heuristic seed truncation /
    //    sector bound from the Kira 3 paper.  Restricts how far the numerator
    //    (scalar-product / ISP) seeding reaches; `l` = number of loops.
    //    Implements the s <= max(1, t-l+1, ...) sector bound.  When
    //    `kira3_truncate_sp` is true the Reduce/Masters jobs emit a
    //    `truncate_sp` block with `l = fc->n_loops()`.
    //  * `run_triangular: sectorwise` -- Kira 3's forward-elimination run
    //    sectorwise.  MEASURED REGRESSION on the aux-mass DE reduce: the
    //    `iterative_reduction: sectorwise` + `run_triangular: sectorwise`
    //    (+ `run_firefly: back`) combination stalls -- the symbolic eta+d
    //    triangular step grinds (observed ~250 s/iter, no progress) and on the
    //    high-seed closure loop it enters an endless "Regenerate:" sectorwise
    //    loop.  Plain `run_firefly: true` (whole finite-field solve after
    //    run_initiate) is the working memory-light fix and completes.  So
    //    `kira3_sectorwise` defaults OFF.
    //
    // `truncate_sp` is family-dependent: on the 3-loop crossed box it collapses
    // the 473-master / 680 GB tower to the genuine 39-master basis at a few GB,
    // and on two-mass 5-point-style families it prevents a Fermat SIG11 -- but on box1 it
    // *destabilises* the reduce.cpp UNION closure loop (the per-pass
    // "extra-master" discovery never closes; r:5 -> r:48,s:41).  We therefore
    // align the default to ReduceOptions::kira3_truncate_sp = false (the
    // measured-safe production value) and surface it as an explicit per-family
    // JSON / env knob.  See ReduceOptions for the full measured-regression
    // notes.
    bool kira3_truncate_sp = false;
    bool kira3_sectorwise  = false;

    // --- Kira 3.1 FireFly -> kira2math EXPORT split (the crossed-box export
    //     fault fix) ---
    //
    // In our runs with Kira 3.1, a combined jobs.yaml that runs
    // `reduce_sectors {... run_firefly: true}` followed by a `kira2math`
    // job in the SAME kira process returned rc=0 from the FireFly reduce
    // but left `results/<fam>/<target>.m` EMPTY (kira.db held only
    // metadata tables), and the in-process kira2math export then aborted
    // with SIG6 / std::bad_alloc despite free memory; splitting into two
    // invocations avoids it (observed in our runs with Kira 3.1 + FireFly;
    // we will report it to the Kira developers with a reproducer -- until
    // then treat it as an observation about our configuration, not a
    // confirmed Kira defect).
    //
    // The fix is to SPLIT the single combined jobs.yaml into TWO separate
    // kira invocations sharing the same kira.db / working directory:
    //   (1) a `run_firefly: true` reduce job (jobs.yaml) that commits the
    //       reconstructed coefficients to kira.db, then
    //   (2) a CLEAN standalone `kira2math` job (jobs_kira2math.yaml) run in
    //       a FRESH kira process that exports <target>.m from the committed
    //       DB.
    // A fresh process means the kira2math export reads the committed
    // EQUATION/FIREFLY tables from disk and never touches the corrupt
    // in-process FireFly regeneration state, so the bad_alloc cannot fire.
    //
    // Only the FireFly Reduce path is affected (the Fermat
    // run_triangular/run_back_substitution path exports kira2math
    // in-process without the fault).  Default: follow `run_firefly` (the
    // split is applied iff the FireFly path is taken), so the Fermat path
    // — and box1's 1-loop regression — is byte-identical.  Force on/off
    // via AMFLOW_KIRA_SPLIT_K2M (1/0) in make_kira_config.
    //
    // When false, `kira_write_jobs` emits the single combined jobs.yaml
    // (reduce + kira2math) exactly as before and the caller runs it once.
    // When true, `kira_write_jobs` emits jobs.yaml (reduce only) PLUS
    // jobs_kira2math.yaml (kira2math only); the caller runs jobs.yaml,
    // then kira_run(..., "jobs_kira2math.yaml").
    bool split_kira2math = false;

    // --- Symbolic-kinematics IBP (Layer B of the ibp_cache) ---
    //
    // When true, `kira_write_config` writes kinematics.yaml WITHOUT
    // baking `numeric_values` into the scalar-product rules /
    // propagator masses (i.e. `NumericSubs::build` returns an empty
    // substitution).  Kira then reduces over {d, η, s12, s23, …} fully
    // symbolically and the output `kira_target.m` carries rational
    // functions of the kinematic invariants.  The Layer-A content-hash
    // cache (`ibp_cache.hpp`) keys on the YAML, so with symbolic_ibp
    // every kinematic point HASHES TO THE SAME KEY → one Kira run is
    // shared across the whole oracle farm.  `reduce.cpp` parses the
    // symbolic output on a wide MpolyContext and `substitute_and_
    // narrow`s the per-point numerics back to {η, d}.
    //
    // Symbolic-kinematics Kira is slower per-reduction (more Fermat /
    // FireFly variables) and can blow up on multi-scale families; gate
    // it behind `AMFLOW_SYMBOLIC_IBP=1` / `"symbolic_ibp": true` and
    // fall back to per-point numerics when it doesn't close.
    bool symbolic_ibp = false;

    std::map<std::string, std::string> numeric_values;

    // --- AMFLOW_DIFFEQ_NUMD: numeric-d substitution inside the diffeq
    //     target_reduce (the etanumd_convert.sh sed INTERNALIZED) ---
    //
    // WHY (measured on a production eta-DE system): the aux-mass-DE
    // `target_reduce` FireFly reconstruction runs over the SYMBOLIC
    // variables of the Kira SYSTEM files — {d, eta} — even for a
    // single-eps black_box run (`kira -sd=...` only reaches the Fermat
    // path; FireFly reads SYSTEM_*.gz directly — a measured pitfall).
    // The 2-var reconstruction is the measured multi-day wall.
    // The equivalent external route
    // (tools/numkin/etarerun/etanumd_convert.sh) substitutes
    // d -> (P/Q) into the generated SYSTEM files and reruns FireFly
    // 1-var (eta).  This field internalizes that route: when non-empty
    // (a rational "P/Q", set from env AMFLOW_DIFFEQ_NUMD by
    // ibp::diffeq ONLY), the Reduce+FireFly leg is split into
    //   (1) kira `jobs_initiate.yaml` (run_initiate only) — writes
    //       tmp/<fam>/SYSTEM_*.gz natively,
    //   (2) in-place token-safe substitution d -> (P/Q) in every
    //       SYSTEM_*.gz (alphabet must be ⊆ {Eq, d, eta}; zero
    //       residual `d` tokens verified),
    //   (3) kira `jobs.yaml` as usual — resumes past initiate (the
    //       documented same-dir resume behavior) and FireFly reconstructs
    //       1-var (eta).
    // Because the substitution happens INSIDE the run, every consumer
    // (DE assembly, determine_boundary_order, boundary towers) is
    // computed natively at numeric d — nothing is injected, so the
    // e00-inject master=26/order=1 BC gap class cannot recur.  The
    // ibp-cache key is salted with the value (numeric-d entries never
    // collide with symbolic ones), and
    // black_box_amflow_single asserts d == 4-2*eps for the single eps
    // sample at startup (the silent-wrong-d killer).  Empty (default)
    // = byte-identical stock behaviour.
    std::string diffeq_numd;

    std::string kira_executable;
    std::string fermat_executable;
};

void kira_write_config(const KiraConfig& cfg, const std::string& dir);

enum class KiraReductionMode {
    Masters,
    Reduce,
};

void kira_write_jobs(const KiraConfig& cfg,
                     const std::string& dir,
                     KiraReductionMode mode);

void kira_write_preferred(const std::vector<qft::JIntegral>& preferred,
                            const qft::FamilyConfig& fc,
                            const std::string& dir);

void kira_write_targets(const std::vector<qft::JIntegral>& targets,
                         const qft::FamilyConfig& fc,
                         const std::string& dir);

// Runs the kira binary in `dir` on the job file `jobs_file` (default
// "jobs.yaml").  For the two-job FireFly export split the caller invokes
// this twice in the SAME `dir` (sharing kira.db): once with the default
// "jobs.yaml" (the run_firefly reduce) and once with
// "jobs_kira2math.yaml" (the standalone export in a fresh process).
double kira_run(const KiraConfig& cfg,
                const std::string& dir,
                const std::string& log_path = "",
                const std::string& jobs_file = "jobs.yaml");

std::vector<qft::JIntegral>
kira_read_masters(const qft::FamilyConfig& fc, const std::string& dir);

struct KiraReductionRule {
    qft::JIntegral lhs;
    std::vector<std::pair<std::string, qft::JIntegral>> rhs;
};

std::vector<KiraReductionRule>
kira_read_target_table(const qft::FamilyConfig& fc, const std::string& dir);

algebra::Mfrac
kira_parse_expression(const std::shared_ptr<algebra::MpolyContext>& ctx,
                       const std::string& expr);

// --- AMFLOW_DIFFEQ_NUMD support (see KiraConfig::diffeq_numd) ---
//
// Read + validate the env var.  Returns "" when unset/empty.  A set
// value must be an integer or reduced-form rational "P/Q" with each of
// |P|, Q < 9.2e18 (Kira parses coefficient tokens with std::stol —
// the etanumd_convert.sh stol guard) and Q > 0; anything else throws.
std::string diffeq_numd_env();

// The pre-FireFly staging: run `jobs_initiate.yaml` (derived from the
// already-written `dir`/jobs.yaml with the solve steps stripped), then
// substitute d -> (P/Q) in every `dir`/tmp/*/SYSTEM_*.gz with the
// etanumd token-safety checks.  Idempotent: a dir whose SYSTEM files
// already carry zero `d` tokens (ffsave resume) is left unchanged.
// Throws on any safety-check failure.
void diffeq_numd_pre_stage(const KiraConfig& cfg,
                           const std::string& dir,
                           const std::string& log_path);

}  // namespace amflow::ibp

#endif  // AMFLOW_IBP_KIRA_HPP
