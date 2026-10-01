// SPDX-License-Identifier: MIT
// ibp::reduce — implementation (flattened BlackBox; Kira-only).
//

#include "amflow/ibp/reduce.hpp"

#include "amflow/algebra/numeric_subst.hpp"
#include "amflow/ibp/ibp_cache.hpp"
#include "amflow/ibp/load_adaptive.hpp"

#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <functional>
#include <initializer_list>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>

#include <flint/fmpz.h>
#include <flint/fmpz_mpoly.h>
#include <flint/fmpz_mpoly_q.h>

namespace amflow::ibp {

using algebra::Mfrac;
using algebra::Mpoly;
using algebra::MpolyContext;

namespace fs = std::filesystem;

ReductionContext make_reduction_context(const qft::FamilyConfig& /*fc*/) {
    // Narrow the reduction MpolyContext to `{eta, d}` — the only
    // variables Kira's output rationals carry symbolically (numeric
    // kinematics are baked into the yaml by `kira_yaml.cpp`'s
    // `NumericSubs`).  Declaring the full `fc.ctx` here drags 9+
    // unused symbols through every FLINT multivariate-GCD inside
    // `fmpz_mpoly_q_canonicalise`, ballooning memory on multi-mass
    // families.
    //
    // This mirrors MMA: upstream `DifferentialEquation`
    // (`Kira/interface.m:505`) stores the diffeq tensor as symbolic
    // Mathematica expressions where unused symbols never enter the
    // algebra.  Variables in `fc.ctx` but not in `red_ctx` (loop
    // momenta, kinematic invariants, masses) must be substituted to
    // numeric values before any Mfrac is lifted — see
    // `substitute_fc_vars`.
    std::vector<std::string> names = {"eta", "d"};
    auto rctx = std::make_shared<MpolyContext>(std::move(names));

    ReductionContext out;
    out.ctx = rctx;
    out.d_var = 1;
    return out;
}

namespace {

using algebra::substitute_fc_vars;
using algebra::mfrac_to_ctx_lenient;

inline Mfrac mfrac_to_ctx(const Mfrac& src,
                           const std::shared_ptr<MpolyContext>& dst_ctx) {
    return algebra::mfrac_to_ctx_lenient(src, dst_ctx);
}

Mfrac lift_to_red_ctx(const Mfrac& src, const ReductionContext& red_ctx) {
    return mfrac_to_ctx(src, red_ctx.ctx);
}

KiraConfig make_kira_config(const qft::FamilyConfig& fc,
                              const std::vector<int>& top_pattern,
                              const ReduceOptions& opts) {
    KiraConfig cfg;
    cfg.fc = &fc;
    cfg.top_pattern = top_pattern;
    cfg.ibp_rank = opts.ibp_rank;
    cfg.ibp_dot  = opts.ibp_dot;
    cfg.n_thread = opts.n_thread;
    // === Upgrade F1: load-adaptive kira -p ===
    // When AMFLOW_KIRA_P_AUTO=1, pick -p from current spare capacity
    // (nproc − loadavg1)/FARM_N, clamped to [1, opts.n_thread].  See
    // amflow/ibp/load_adaptive.hpp for the rationale (14-job the 3-loop crossed-box farm
    // at load=252/96 cores ⇒ 2.6× oversubscribed ⇒ ~2× wall penalty).
    if (const char* v = std::getenv("AMFLOW_KIRA_P_AUTO");
        v && *v && !(v[0] == '0' && v[1] == '\0')) {
        long auto_p = load_adaptive_threads(opts.n_thread);
        if (auto_p != cfg.n_thread) {
            std::cerr << "[amflow][load-adaptive] kira -p"
                      << opts.n_thread << " -> -p" << auto_p
                      << " (nproc=" << load_adaptive_nproc()
                      << " load1=" << load_adaptive_loadavg1()
                      << " farm_n=" << load_adaptive_farm_n() << ")\n";
        }
        cfg.n_thread = auto_p;
    }
    cfg.integral_order = opts.integral_order;
    cfg.numeric_values = opts.numeric_values;
    cfg.kira_executable = opts.kira_executable;
    cfg.fermat_executable = opts.fermat_executable;

    // FireFly finite-field reconstruction for the Reduce-mode solve.
    // Default comes from the ReduceOptions, but is overridable per-run
    // via env without a rebuild (AMFLOW_RUN_FIREFLY=1 forces it on for
    // every Reduce pass; =0 forces it off).  This is the memory-light
    // fix for the aux-mass+eta symbolic coefficient swell.
    cfg.run_firefly   = opts.run_firefly;
    cfg.firefly_bunch = opts.firefly_bunch;

    // AUTO-FIREFLY ON AUX-MASS (η) DEFORMATION.  See
    // ReduceOptions::auto_firefly.  When the family context carries the
    // auxiliary-mass invariant `eta` (i.e. this is an AMFlow aux-mass DE
    // reduce, NOT an ordinary physical reduction), the Fermat back-
    // substitution backend aborts on the η-polynomial master tower
    // (SIGABRT — an internal assertion at ~35-45 GB peak in our runs, or
    // std::bad_alloc under a memory cap; SIGSEGV on some families).  FireFly finite-field reconstruction is the memory-
    // light fix.  Default the Reduce backend to FireFly for exactly
    // these η-deformed reductions.  A non-η family has var_index("eta")<0
    // and is a strict no-op (Fermat path byte-identical).  Explicit
    // choices override below: JSON/env AMFLOW_RUN_FIREFLY wins, and
    // auto_firefly=false / AMFLOW_NO_AUTO_FIREFLY=1 suppresses the route.
    {
        bool auto_on = opts.auto_firefly;
        if (const char* v = std::getenv("AMFLOW_NO_AUTO_FIREFLY"); v && *v) {
            auto_on = (v[0] == '0' && v[1] == '\0');
        }
        if (auto_on && !cfg.run_firefly && fc.var_index("eta") >= 0) {
            cfg.run_firefly = true;
        }
        // VERIFIED ETA-DEFORMED SIG6 FIX (the load-bearing companion to
        // auto-FireFly).  FireFly alone STALLS on the uncapped monolithic
        // eta top-sector reduce (the aux-mass-DE Regenerate loop, mandatory
        // list ~499 at r:17,s:3).  The verified fix is the COMBINATION:
        // per-subsector dot-cap (collapses the seed to d:ibp_dot, dropping
        // the eta top-sector mandatory list 499 -> 260 and the seed
        // r:17,s:3 -> r:11,s:1) PLUS FireFly.  Default the dot-cap ON for
        // exactly the eta families that take the FireFly route, so the fix
        // is automatic with no JSON/env change.  Non-eta families
        // (var_index("eta") < 0) are a strict no-op.  Env
        // AMFLOW_REDUCE_DCAP overrides per-run (kira_yaml.cpp).
        if (auto_on && fc.var_index("eta") >= 0) {
            cfg.reduce_dcap = true;
        }
    }

    if (const char* v = std::getenv("AMFLOW_RUN_FIREFLY"); v && *v) {
        cfg.run_firefly = !(v[0] == '0' && v[1] == '\0');
    }
    if (const char* v = std::getenv("AMFLOW_FIREFLY_BUNCH"); v && *v) {
        try { cfg.firefly_bunch = std::stol(v); } catch (...) {}
    }

    // Kira 3 seeding levers (arXiv:2505.20197): heuristic seed truncation
    // (truncate_sp, l = #loops) + sectorwise forward elimination.  BOTH
    // default OFF in ReduceOptions (truncate_sp regressed on box1, see
    // ReduceOptions comment); per-family overridable from the
    // amf_options.blackbox JSON or AMFLOW_KIRA3_* env without a rebuild.
    cfg.kira3_truncate_sp = opts.kira3_truncate_sp;
    cfg.kira3_sectorwise  = opts.kira3_sectorwise;
    if (const char* v = std::getenv("AMFLOW_KIRA3_TRUNCATE_SP"); v && *v) {
        cfg.kira3_truncate_sp = !(v[0] == '0' && v[1] == '\0');
    }
    if (const char* v = std::getenv("AMFLOW_KIRA3_SECTORWISE"); v && *v) {
        cfg.kira3_sectorwise = !(v[0] == '0' && v[1] == '\0');
    }

    // Kira 3.1 FireFly -> kira2math export split (see
    // KiraConfig::split_kira2math).  DEFAULT follows `run_firefly`: the
    // in-process kira2math-after-firefly SIG6/bad_alloc only afflicts the
    // FireFly Reduce path, so we split iff FireFly is active and leave the
    // Fermat path (box1 et al.) byte-identical.  AMFLOW_KIRA_SPLIT_K2M
    // (1/0) forces it on/off for experimentation.
    cfg.split_kira2math = cfg.run_firefly;
    if (const char* v = std::getenv("AMFLOW_KIRA_SPLIT_K2M"); v && *v) {
        cfg.split_kira2math = !(v[0] == '0' && v[1] == '\0');
    }

    // Symbolic-kinematics IBP (Layer B of the ibp_cache).  See
    // KiraConfig::symbolic_ibp.  Default from ReduceOptions; per-run
    // overridable via AMFLOW_SYMBOLIC_IBP=1/0 without a rebuild.
    cfg.symbolic_ibp = opts.symbolic_ibp;
    if (const char* v = std::getenv("AMFLOW_SYMBOLIC_IBP"); v && *v) {
        cfg.symbolic_ibp = !(v[0] == '0' && v[1] == '\0');
    }
    return cfg;
}

// === IBP results cache wrapper (Layer A of `ibp_cache.hpp`) ===
//
// Called AFTER all Kira input files have been written into `dir` and
// BEFORE `kira_run`.  Hashes the inputs, looks up the content-addressed
// cache, and on a HIT restores `results/<fam>/{masters,kira_target.m}`
// into `dir` so the subsequent `kira_read_*` calls see exactly what a
// fresh Kira run would have produced.  On a MISS, runs Kira (and the
// optional split kira2math export so the cached entry always carries a
// populated `kira_target.m`), then stores.  When the cache root is
// empty (no JSON `ibp_cache_dir` and no `AMFLOW_IBP_CACHE` env), this
// degenerates to a single `kira_run` — byte-identical to the pre-cache
// path bar the MISS log line.
// === SIG6/SIG11 -> FireFly auto-staging (the monolithic-reduce fix) ===
//
// The Fermat algebraic backend (`run_triangular`/`run_back_substitution`)
// builds the dense symbolic {eta,d} solution in one shot.  On an
// auxiliary-mass eta-deformed top-sector reduce with a large mandatory list
// (e.g. an 8-prop -> 15-prop blackbox basis, mandatory ~385-499) the
// symbolic eta+d expression swell makes the Kira/Fermat path abort (SIGABRT —
// an internal assertion at ~35-45 GB peak in our runs, or std::bad_alloc under
// a memory cap; SIGSEGV on some families), so raising memory caps alone does
// not clear it.
//
// The proven fix is to STAGE the same reduction through FireFly finite-field
// reconstruction: FireFly never materialises the monolithic symbolic system
// (it solves over F_p at numeric {eta,d} points and reconstructs the
// rationals), so peak RAM is tiny and the assert cannot fire.  This wrapper
// makes that staging AUTOMATIC: when the Fermat path aborts on a SIGNAL we
// re-emit jobs.yaml with `run_firefly: true` (+ the split kira2math export)
// into the SAME working dir and retry once.  A config that already requested
// FireFly is not retried (nothing to fall back to).
//
// Disable with AMFLOW_SIG6_FIREFLY_FALLBACK=0 (e.g. to surface the raw
// abort for diagnosis).  Byte-identical to the prior path on any run that
// does not abort.
namespace {
bool sig6_firefly_fallback_enabled() {
    const char* v = std::getenv("AMFLOW_SIG6_FIREFLY_FALLBACK");
    if (v == nullptr) return true;                 // default ON
    return !(v[0] == '0' && v[1] == '\0');
}
bool is_signal_abort(const std::exception& e) {
    const std::string m = e.what();
    return m.find("KILLED BY SIGNAL") != std::string::npos;
}

// AMFLOW_KIRA_FFSAVE resume-engagement guard support: strip a preserved
// dir back to its INPUTS (config/, jobs*.yaml, preferred, target,
// sidecar and the append-mode logs stay; all Kira/FireFly state goes),
// so a failed resume can be retried once from a clean slate.  A save
// killed inside FireFly's pre-shift window is UNRESUMABLE by design
// (FireFly hard-fails `Shift file not found!` -> kira rc=1, measured in
// a mid-save kill synthetic); without this guard that relaunch
// would be strictly WORSE than stock.
void ffsave_clear_state(const std::string& dir) {
    namespace fsx = std::filesystem;
    std::error_code ec;
    for (const char* rel : {"ff_save", "ff_save_tmp", "firefly_saves",
                            "results", "tmp", "sectormappings"}) {
        fsx::remove_all(fsx::path(dir) / rel, ec);
    }
}

// A monolithic-eta FireFly "Regenerate" deadlock surfaces NOT as a
// signal (the tree never dies on its own) but as the stall-watchdog
// error thrown by kira_run (see kira_run.cpp).  Treat it as a
// recoverable abort: re-stage through the dot-capped + FireFly path so
// the hang can never stand.  Unlike a signal abort, a stall can occur
// even when FireFly is ALREADY on -- FireFly alone stalls on the
// uncapped mandatory list -- so the recovery's load-bearing move is to
// add the dot-cap (reduce_dcap), not merely to enable FireFly.
bool is_stall_abort(const std::exception& e) {
    return std::string(e.what()).find("STALLED") != std::string::npos;
}
}  // namespace

void kira_run_cached(const KiraConfig& cfg,
                     const std::string& dir,
                     const ReduceOptions& opts,
                     bool run_split_export,
                     KiraReductionMode mode,
                     bool ffsave_resumed = false) {
    // AMFLOW_DIFFEQ_NUMD (KiraConfig::diffeq_numd, set by ibp::diffeq
    // only): active exactly when this is a Reduce-mode FireFly leg.
    // Masters-mode enumeration never solves the system (its result is
    // d-independent) and keeps the unsalted key so numd and symbolic
    // runs SHARE preheat cache entries; the Fermat path already gets
    // numeric d via `-sd` and is left untouched.
    const bool numd_active = !cfg.diffeq_numd.empty() &&
                             mode == KiraReductionMode::Reduce &&
                             cfg.run_firefly;
    // Run the Fermat (or already-FireFly) path, with a single automatic
    // re-stage if it (a) aborts on a SIGNAL (SIGABRT/SIGSEGV, see above) with
    // FireFly not yet in use, or (b) STALLS (0 CPU, observed with monolithic-
    // eta FireFly reduces;
    // caught by the kira_run watchdog) with the dot-cap not yet applied.
    // Lambda so both the cache-MISS and no-cache branches share it.
    auto run_with_staging = [&]() {
        try {
            // Numeric-d pre-stage: run_initiate only, then substitute
            // d -> (P/Q) in the generated SYSTEM files; the jobs.yaml
            // run below then resumes past initiate and FireFly
            // reconstructs 1-var (eta).  See diffeq_numd_pre_stage.
            if (numd_active) diffeq_numd_pre_stage(cfg, dir, opts.log_file);
            kira_run(cfg, dir, opts.log_file);
            if (run_split_export && cfg.split_kira2math) {
                kira_run(cfg, dir, opts.log_file, "jobs_kira2math.yaml");
            }
        } catch (const std::exception& e) {
            const bool sig   = is_signal_abort(e);
            const bool stall = is_stall_abort(e);
            // A signal abort is recoverable only if FireFly is not yet in
            // use (the OOM fix is to switch to finite-field FireFly).  A
            // stall is recoverable only if the dot-cap is not yet applied
            // (the deadlock fix is the dot-cap + FireFly combination;
            // FireFly alone stalls).  Otherwise we have already applied the
            // strongest staging we have -- re-throw.
            const bool recoverable_sig   = sig   && !cfg.run_firefly;
            const bool recoverable_stall = stall && !cfg.reduce_dcap;
            if ((!recoverable_sig && !recoverable_stall) ||
                !sig6_firefly_fallback_enabled()) {
                throw;
            }
            std::cerr << "[amflow][stall-staging] kira reduce "
                      << (stall ? "STALLED" : "aborted") << " ("
                      << e.what()
                      << "); RE-STAGING via dot-cap + FireFly finite-field "
                         "reconstruction in " << dir << "\n";
            // Re-emit jobs with the full verified eta fix: FireFly +
            // dot-cap + split export, then retry once.  Setting reduce_dcap
            // collapses the monolithic mandatory list (the load-bearing
            // half of the fix) so FireFly no longer enters the Regenerate
            // deadlock; run_firefly keeps peak RAM tiny.
            KiraConfig ff = cfg;
            ff.run_firefly    = true;
            ff.reduce_dcap    = true;
            ff.split_kira2math = true;
            kira_write_jobs(ff, dir, mode);
            kira_run(ff, dir, opts.log_file);
            if (run_split_export && ff.split_kira2math) {
                kira_run(ff, dir, opts.log_file, "jobs_kira2math.yaml");
            }
        }
    };

    // AMFLOW_KIRA_FFSAVE resume-engagement guard (the run_with_staging
    // wrapper actually invoked below): when THIS call runs Kira inside a
    // PRESERVED dir (ffsave_resumed) and Kira fails on a NON-signal
    // error — e.g. FireFly `Shift file not found!` on a save interrupted in
    // the pre-shift window, or a corrupt mid-initiate kira.db — the
    // preserved state was not resumable.  Mirror the standalone
    // engagement law (i168 discriminator: NOT_ENGAGED -> treat as
    // fresh): clear the Kira/FireFly state (inputs stay) and retry ONCE
    // from scratch.  Signal aborts rethrow unchanged: an external
    // external/OOM kill must leave the (still-valid) state for the next
    // relaunch, and the SIG6 staging above already owns that class.
    auto run_attempt = [&]() {
        if (!ffsave_resumed) {
            run_with_staging();
            return;
        }
        try {
            run_with_staging();
        } catch (const std::exception& e) {
            if (is_signal_abort(e)) throw;
            std::cerr << "[amflow][ffsave] RESUME NOT ENGAGED in " << dir
                      << " (" << e.what()
                      << "); clearing preserved state and retrying ONCE"
                         " fresh\n";
            ffsave_clear_state(dir);
            run_with_staging();
        }
    };

    const std::string root = ibp_cache_resolve_root(opts.ibp_cache_dir);
    const std::string fam  = cfg.fc->family;
    if (!root.empty()) {
        // Numeric-d legs salt the key with the substituted value: the
        // substitution rewrites tmp/SYSTEM files that are NOT hashed
        // inputs, so the unsalted key of a numd run is byte-identical
        // to the symbolic run's (the eps-independent-key
        // collision class) and the two must never share an entry.
        const std::string key = ibp_cache_key(
            dir, numd_active ? "AMFLOW_DIFFEQ_NUMD=" + cfg.diffeq_numd
                             : std::string());
        if (ibp_cache_lookup(root, key, fam, dir)) return;
        // === Stampede lock ===
        // First lookup MISSed.  Take the per-key advisory flock so that
        // when N farm jobs reach the same (s,t)-independent IBP stage at
        // once, only ONE runs Kira; the rest block here, then re-check
        // the cache → HIT.  When the lock is disabled (or the cache root
        // empty) `IbpCacheLock` is a no-op and the path is byte-identical
        // to the pre-upgrade flow.
        IbpCacheLock lock(root, key);
        if (lock.waited()) {
            // We blocked on another job; the entry is very likely now
            // populated.  Re-check before running Kira ourselves.
            if (ibp_cache_lookup(root, key, fam, dir)) return;
        }
        std::cerr << "[amflow][ibp-cache] MISS key=" << key
                  << " fam=" << fam
                  << " symbolic_ibp=" << (cfg.symbolic_ibp ? 1 : 0)
                  << " ... running Kira\n";
        run_attempt();
        ibp_cache_store(root, key, fam, dir);
        return;
    }
    run_attempt();
}

// === AMFLOW_KIRA_FFSAVE: embedded FireFly save/resume persistence ===
//
// WHY: the embedded Kira reduce runs FireFly with Kira's own save
// mechanism ALWAYS ON (the pinned fork's reductionFF.cpp hardcodes
// `bool save = true`; `ff_save/` accumulates per-function saved states
// during the run), but the amflow-level relaunch re-enters
// `run_reduce_pass` and hits the unconditional
// `fs::remove_all(reduce_dir)`, restarting a multi-hour reconstruction
// from probe zero.  A manual `kira -pN jobs.yaml` rerun in the SAME dir
// (the documented Kira resume behavior) works precisely because Kira
// resumes: `ff_save/states` present =>
// `reconst->resume_from_saved_state()` fast-forwards past every done
// function.  This block makes that route NATIVE, gated on
// AMFLOW_KIRA_FFSAVE=1 (default OFF = stock wipe, byte-identical).
//
// Job identity: the Kira input files (config/*.yaml, jobs.yaml,
// jobs_kira2math.yaml, preferred, target) are deterministic functions of
// the KiraConfig, so "same bytes" == "same reduction job" — the same
// invariant the Layer-A ibp_cache key relies on — EXCEPT for argv-only
// parameters that never reach the input files: the `-sd=4-2*<eps>`
// forwarding of a numeric eps (kira_run.cpp) and the Kira binary itself
// (a kira/kira128 weight-width swap changes the ff_save weight
// encoding).  Those go into a fingerprint sidecar `.amflow_ffsave.job`
// that must byte-match too.  The sidecar is NOT in ibp_cache's
// kInputFiles, so cache keys are unperturbed.
bool kira_ffsave_enabled_impl() {
    const char* v = std::getenv("AMFLOW_KIRA_FFSAVE");
    return v && *v && !(v[0] == '0' && v[1] == '\0');
}

namespace {

// Result-determining job identity that is NOT captured by the Kira
// input files: numeric_values (eps is argv-only via `-sd`; the others
// are baked into kinematics.yaml but hashing them again is harmless and
// defensive) and the Kira executable (ff_save state is weight-width
// specific).  Parallelism knobs (n_thread, firefly_bunch) are
// deliberately EXCLUDED: FireFly save states are valid across them.
std::string ffsave_fingerprint(const KiraConfig& cfg) {
    std::ostringstream fp;
    fp << "family=" << cfg.fc->family << "\n";
    fp << "kira=" << cfg.kira_executable << "\n";
    for (const auto& [name, value] : cfg.numeric_values) {
        fp << "nv:" << name << "=" << value << "\n";
    }
    // AMFLOW_DIFFEQ_NUMD: a numeric-d dir (substituted SYSTEM files +
    // 1-var ff_save state) must never be resumed by a symbolic run or
    // vice versa — the input files alone are byte-identical between the
    // two, so the fingerprint carries the substitution value.  Empty
    // (default) adds nothing: byte-identical legacy fingerprint.
    if (!cfg.diffeq_numd.empty()) {
        fp << "diffeq_numd=" << cfg.diffeq_numd << "\n";
    }
    return fp.str();
}

constexpr const char* kFfsaveSidecar = ".amflow_ffsave.job";

bool ffsave_read_file(const fs::path& p, std::string& out) {
    std::ifstream f(p, std::ios::binary);
    if (!f) return false;
    std::ostringstream ss;
    ss << f.rdbuf();
    out = ss.str();
    return true;
}

// Every regular file under `probe` must exist under `live` with equal
// bytes.  (Files present only in `live` — kira outputs, logs, ff_save —
// are the state we are preserving and do not count against the match.)
bool ffsave_inputs_match(const fs::path& probe, const fs::path& live) {
    std::error_code ec;
    for (auto it = fs::recursive_directory_iterator(probe, ec);
         !ec && it != fs::recursive_directory_iterator(); ++it) {
        if (!it->is_regular_file(ec)) continue;
        fs::path rel = fs::relative(it->path(), probe, ec);
        if (ec) return false;
        std::string a, b;
        if (!ffsave_read_file(it->path(), a)) return false;
        if (!ffsave_read_file(live / rel, b)) return false;
        if (a != b) return false;
    }
    return !ec;
}

std::size_t ffsave_count_states(const fs::path& dir) {
    std::error_code ec;
    std::size_t n = 0;
    for (auto it = fs::directory_iterator(dir / "ff_save" / "states", ec);
         !ec && it != fs::directory_iterator(); ++it) ++n;
    return n;
}

// Prepare `dir` for a Kira (re)run whose input files are produced by
// `write_inputs`.  Knob OFF: stock `remove_all` + recreate + write —
// byte-identical to the pre-knob path.  Knob ON: write the inputs into
// a sibling probe dir, byte-compare against the live dir; on a full
// match PRESERVE the live dir (Kira resumes from ff_save/kira.db/tmp),
// otherwise archive the old dir to ONE `.ffsave_stale` slot (partial
// state stays diagnosable) and promote the probe.  Returns true iff the
// live dir was preserved for resume.
bool ffsave_prepare_kira_dir(
        const std::string& dir,
        const std::string& fingerprint,
        const std::function<void(const std::string&)>& write_inputs) {
    if (!kira_ffsave_enabled_impl()) {
        fs::remove_all(dir);
        fs::create_directories(dir);
        write_inputs(dir);
        return false;
    }

    const fs::path live(dir);
    const fs::path probe(dir + ".ffsave_probe");
    std::error_code ec;
    fs::remove_all(probe, ec);
    fs::create_directories(probe);
    write_inputs(probe.string());
    {
        std::ofstream f(probe / kFfsaveSidecar, std::ios::binary);
        f << fingerprint;
    }

    if (fs::exists(live, ec) && ffsave_inputs_match(probe, live)) {
        std::size_t n_states = ffsave_count_states(live);
        std::cerr << "[amflow][ffsave] PRESERVE dir=" << dir
                  << " (inputs+fingerprint match; ff_save states="
                  << n_states << " — Kira resumes in place)\n";
        fs::remove_all(probe, ec);
        return true;
    }

    if (fs::exists(live, ec)) {
        std::size_t n_states = ffsave_count_states(live);
        if (n_states > 0) {
            const fs::path stale(dir + ".ffsave_stale");
            fs::remove_all(stale, ec);
            fs::rename(live, stale, ec);
            std::cerr << "[amflow][ffsave] FRESH dir=" << dir
                      << " (inputs differ — different job); old dir had "
                      << n_states << " ff_save states, archived to "
                      << stale.string() << " (one slot)\n";
            if (ec) fs::remove_all(live);  // rename failed: fall back
        } else {
            std::cerr << "[amflow][ffsave] FRESH dir=" << dir
                      << " (inputs differ or no prior state)\n";
            fs::remove_all(live);
        }
    }
    fs::rename(probe, live, ec);
    if (ec) {
        // Cross-device or races: fall back to write-in-place.
        fs::remove_all(live);
        fs::create_directories(live);
        write_inputs(dir);
        std::ofstream f(live / kFfsaveSidecar, std::ios::binary);
        f << fingerprint;
        fs::remove_all(probe, ec);
    }
    return false;
}

}  // namespace

// Parse a Kira coefficient string into the narrow {η, d} reduction
// context.  In the default (numeric-kinematics) path Kira's output
// already lives on {η, d} and this is just `kira_parse_expression`.
// Under `symbolic_ibp` the output additionally carries the kinematic
// invariants {s12, s23, …}; we parse on a wide context that names them
// and `substitute_and_narrow` the per-point `numeric_values` back to
// {η, d} so the rest of the pipeline (DE assembly, Mfrac arithmetic on
// the narrow ctx) is unchanged.
Mfrac parse_kira_coef(const KiraConfig& cfg,
                      const ReduceOptions& opts,
                      const std::shared_ptr<MpolyContext>& red_ctx,
                      const std::shared_ptr<MpolyContext>& parse_ctx,
                      const std::string& coef_str) {
    if (!cfg.symbolic_ibp) {
        return kira_parse_expression(red_ctx, coef_str);
    }
    Mfrac wide = kira_parse_expression(parse_ctx, coef_str);
    return algebra::substitute_and_narrow(wide, opts.numeric_values, red_ctx);
}

// Build the parse context for `parse_kira_coef`: {η, d} plus, under
// symbolic_ibp, every kinematic name the user supplied a numeric value
// for (those are exactly the symbols Kira can emit, since they are the
// `kinematic_invariants` left unsubstituted in the YAML).
std::shared_ptr<MpolyContext>
make_parse_ctx(const KiraConfig& cfg,
               const ReduceOptions& opts,
               const std::shared_ptr<MpolyContext>& red_ctx) {
    if (!cfg.symbolic_ibp) return red_ctx;
    std::vector<std::string> names;
    for (long i = 0; i < red_ctx->n_vars(); ++i) {
        names.push_back(red_ctx->var_name(i));
    }
    std::set<std::string> seen(names.begin(), names.end());
    // Kinematic invariants Kira may emit.  `numeric_values` keys are the
    // user-declared scales (s12, mtsq, …); also union in `fc->invariants`
    // so a symbol the user forgot to pin still parses (it will then throw
    // a clear "no numeric value for 'X'" in substitute_and_narrow rather
    // than the opaque "unknown identifier 'X'" from the parser).
    for (const auto& [k, _] : opts.numeric_values) {
        if (k == "eps") continue;
        if (seen.insert(k).second) names.push_back(k);
    }
    for (const auto& inv : cfg.fc->invariants) {
        if (seen.insert(inv).second) names.push_back(inv);
    }
    return std::make_shared<MpolyContext>(std::move(names));
}

ReduceOptions
apply_jdot_jrank_floor(const ReduceOptions& opts,
                       std::initializer_list<const std::vector<qft::JIntegral>*> lists,
                       bool dot_plus_one) {
    ReduceOptions out = opts;
    for (const auto* lst : lists) {
        for (const auto& j : *lst) {
            long d = j.n_dots() + (dot_plus_one ? 1 : 0);
            long r = j.rank();
            if (d > out.ibp_dot)  out.ibp_dot  = d;
            if (r > out.ibp_rank) out.ibp_rank = r;
        }
    }
    return out;
}

// Kira's raw master list is already in MMA-compatible ascending-by-
// sector order; MMA's `Preferred[preferred, dir]`
// (`Kira/interface.m:164-170`) passes it on unchanged.  This function
// is intentionally a no-op — kept as a hook so call sites stay
// visible if a future Kira release produces a non-MMA-compatible raw
// order.
void sort_integrals_like_amflow(std::vector<qft::JIntegral>& /*masters*/) {
}

}  // namespace

// Public knob probe (reduce.hpp): lets the pipeline call sites
// (solve_integrals/amfsystem hygiene wipes, kira_run log-append) share
// the ONE env definition.
bool kira_ffsave_enabled() {
    return kira_ffsave_enabled_impl();
}

ReduceResult
reduce(const qft::FamilyConfig& fc,
       const std::vector<qft::JIntegral>& targets,
       const std::vector<qft::JIntegral>& preferred,
       const std::vector<int>& top_pattern,
       const ReduceOptions& opts) {
    ReduceResult out;
    out.red_ctx = make_reduction_context(fc);

    if (targets.empty()) return out;

    if (top_pattern.size() != fc.propagators_after_conservation.size()) {
        throw std::invalid_argument(
            "ibp::reduce: top_pattern length != #propagators");
    }

    ReduceOptions opts_eff =
        apply_jdot_jrank_floor(opts, {&targets, &preferred}, false);

    // === Boundary-recursion Masters-preheat staging cap ===
    //
    // WHY: `apply_jdot_jrank_floor` above OVERRIDES the configured
    // blackbox ibp_dot/ibp_rank upward to the max dots/rank carried by
    // the targets.  For AMFSystem::build_boundary's recursive sub-family
    // reduce the targets are REGION-EXPANDED boundary integrands (dots
    // and rank up to the expansion depth, 3 on the production nonplanar
    // family), so the
    // Masters-mode `select_mandatory_recursively` seeds at
    // (r: n_top+3, s: 3, d: 3) over the sub-family's whole sector tree.
    // On a production 22-position family (12-prop top sector
    // 4095) that is a 4,571,222-integral mandatory closure ->
    // std::bad_alloc inside Kira's Select-equations at <1 G RSS
    // (measured).  Neither
    // AMFLOW_PREHEAT_DOT (read only in ibp::diffeq) nor
    // AMFLOW_REDUCE_DCAP (Reduce-mode reduce_sectors line only) reaches
    // this Masters-mode seed.
    //
    // FIX (env-gated, default OFF = byte-identical staging):
    // when the call is marked `boundary_stage` (set ONLY by
    // build_boundary), AMFLOW_BOUNDARY_PREHEAT_DOT /
    // AMFLOW_BOUNDARY_PREHEAT_RANK CAP the Masters-mode preheat seed
    // (min with the auto-floor; never raised).  The target_reduce stage
    // KEEPS the floored (rank, dot) so every dotted/ranked target stays
    // inside Kira's seeded system (reduce.cpp:885-class correctness).
    // Master-detection safety: the Route-B precedent (ibp::diffeq on the 3-loop crossed box,
    // mandatory 3584 -> 22 at the pinned dot, basis unchanged) plus this
    // function's reduce-side master UNION (`n_added` block below), which
    // appends any Reduce-mode master the capped preheat enumeration
    // missed.  If the capped preheat returns ZERO masters we re-run it
    // UNCAPPED before concluding "boundary reduces to zero" — a capped
    // seed must never silently drop a boundary family.
    // Cache note: the Masters-mode `d:` is result-determining and IS
    // hashed (ibp_cache.cpp jobs-canonicalization comment), so capped
    // and uncapped preheats key separately — no poisoning, but a
    // relaunch that drops the env goes cold (same pitfall class as the
    // AMFLOW_REDUCE_DCAP row).
    ReduceOptions opts_preheat = opts_eff;
    bool preheat_capped = false;
    if (opts.boundary_stage) {
        auto env_cap = [](const char* k) -> long {
            const char* v = std::getenv(k);
            if (!v || !*v) return -1;
            try { return std::stol(v); } catch (...) { return -1; }
        };
        const long dcap = env_cap("AMFLOW_BOUNDARY_PREHEAT_DOT");
        const long scap = env_cap("AMFLOW_BOUNDARY_PREHEAT_RANK");
        if (dcap >= 0 && dcap < opts_preheat.ibp_dot) {
            std::cerr << "[amflow][boundary-staging] fam=" << fc.family
                      << " masters_preheat dot seed CAPPED "
                      << opts_preheat.ibp_dot << " -> " << dcap
                      << " (AMFLOW_BOUNDARY_PREHEAT_DOT; target_reduce"
                         " keeps the floored seed)\n";
            opts_preheat.ibp_dot = dcap;
            preheat_capped = true;
        }
        if (scap >= 0 && scap < opts_preheat.ibp_rank) {
            std::cerr << "[amflow][boundary-staging] fam=" << fc.family
                      << " masters_preheat rank seed CAPPED "
                      << opts_preheat.ibp_rank << " -> " << scap
                      << " (AMFLOW_BOUNDARY_PREHEAT_RANK; target_reduce"
                         " keeps the floored seed)\n";
            opts_preheat.ibp_rank = scap;
            preheat_capped = true;
        }
    }

    KiraConfig cfg = make_kira_config(fc, top_pattern, opts_eff);
    KiraConfig cfg_preheat =
        preheat_capped ? make_kira_config(fc, top_pattern, opts_preheat)
                       : cfg;
    auto parse_ctx = make_parse_ctx(cfg, opts, out.red_ctx.ctx);

    // Mirror upstream `BlackBoxReduce` (Kira/interface.m): a two-Kira-call
    // sequence at the same `(rank, dot)`.  First call (`IBPSystem`) is
    // `select_mandatory_recursively` for a sector-wide master
    // enumeration; second call (`AnalyticReduction`) is
    // `select_mandatory_list` for the specific targets.  Upstream
    // shares one `$ReductionDirectory` between the two calls; our
    // Kira 2.x version refuses to do that because the Masters-mode
    // `run_initiate: masters` doesn't register the `-s` numeric
    // substitutions, so the subsequent Reduce-mode call sees a stale
    // auxiliary state and Aborts with `Kira::update_auxiliary_file:
    // Last Kira run set 0 variables to numeric values, this time you
    // request N.`  We therefore run the two calls in *separate*
    // subdirectories (`masters_preheat/` and `target_reduce/`) — both
    // at the same `(rank, dot)`, so Kira returns consistent master
    // lists — and bridge them with the SubsetQ check below mirroring
    // upstream's
    // `If[!SubsetQ[masters, str], Abort["inconsistent masters from Kira"]]`.
    const std::string masters_dir = opts.work_dir + "/masters_preheat";
    const std::string reduce_dir  = opts.work_dir + "/target_reduce";

    // === Step 1: Masters-mode preheat (sector-wide enumeration) ===
    // Uses `cfg_preheat` (== `cfg` unless the boundary-staging cap is
    // active) so ONLY this stage's seed shrinks; the target_reduce
    // below stays on the floored `cfg`.
    fs::remove_all(masters_dir);
    fs::create_directories(masters_dir);
    kira_write_config(cfg_preheat, masters_dir);
    kira_write_preferred(preferred, fc, masters_dir);
    kira_write_jobs(cfg_preheat, masters_dir, KiraReductionMode::Masters);
    kira_run_cached(cfg_preheat, masters_dir, opts,
                    /*run_split_export=*/false,
                    KiraReductionMode::Masters);

    out.masters = kira_read_masters(fc, masters_dir);
    if (out.masters.empty() && preheat_capped) {
        // A capped preheat must NOT convert "basis exists" into
        // "boundary reduces to zero" (the caller skips the family on
        // empty masters).  Re-run once at the UNCAPPED floor to rule
        // out seed truncation before accepting the empty enumeration.
        std::cerr << "[amflow][boundary-staging] fam=" << fc.family
                  << " capped preheat returned 0 masters — re-running"
                     " UNCAPPED to rule out seed truncation\n";
        fs::remove_all(masters_dir);
        fs::create_directories(masters_dir);
        kira_write_config(cfg, masters_dir);
        kira_write_preferred(preferred, fc, masters_dir);
        kira_write_jobs(cfg, masters_dir, KiraReductionMode::Masters);
        kira_run_cached(cfg, masters_dir, opts,
                        /*run_split_export=*/false,
                        KiraReductionMode::Masters);
        out.masters = kira_read_masters(fc, masters_dir);
    }
    if (out.masters.empty()) {
        // Empty masters_mma means Kira identified the sector as scaleless
        // (or otherwise had no integrals to enumerate at the chosen
        // (rank, dot)).  Mirror upstream `AnalyticReduction`
        // (`Kira/interface.m:461`): `If[Length[masters]===0, Return[{{},{}}]]`.
        // Caller (AMFSystem::build_boundary) treats empty `red.masters` as
        // "boundary integrand reduces to zero, skip family".
        out.rules.reserve(targets.size());
        for (std::size_t i = 0; i < targets.size(); ++i) out.rules.emplace_back();
        return out;
    }
    sort_integrals_like_amflow(out.masters);

    // Snapshot the canonical (preheat) master set for the SubsetQ check.
    auto sector_master_keys = std::set<std::string>();
    for (const auto& m : out.masters) {
        std::ostringstream key;
        key << m.family();
        for (long v : m.indices()) key << "|" << v;
        sector_master_keys.insert(key.str());
    }

    // === Step 2: Reduce-mode (target reduction in a sibling subdir) ===
    //
    // **Preferred file**: mirror MMA `AnalyticReduction`, which reuses
    // the preferred file IBPSystem wrote (the *input* `preferred`).
    // Do NOT write the sorted preheat output as preferred here; it
    // changes Kira's elimination order and breaks the resulting
    // reduction rule structure.
    // AMFLOW_KIRA_FFSAVE=1: preserve an input-identical existing dir so
    // a re-fired Kira resumes its FireFly ff_save state (see the
    // ffsave_prepare_kira_dir block above).  Knob OFF (default) takes
    // the stock remove_all+recreate+write path, byte-identical.
    const bool ffsave_resumed = ffsave_prepare_kira_dir(
        reduce_dir, ffsave_fingerprint(cfg),
        [&](const std::string& d) {
            kira_write_config(cfg, d);
            kira_write_preferred(preferred, fc, d);  // mirror MMA: reuse input preferred
            kira_write_targets(targets, fc, d);
            kira_write_jobs(cfg, d, KiraReductionMode::Reduce);
        });
    // Layer-A IBP cache: hash inputs → on HIT skip Kira; on MISS run +
    // (under FireFly split) export kira2math now so the cached entry
    // carries `kira_target.m`, then store.  See `kira_run_cached`.
    kira_run_cached(cfg, reduce_dir, opts, /*run_split_export=*/true,
                    KiraReductionMode::Reduce, ffsave_resumed);

    // === Master-basis UNION fix (distinct-mass AnalyticReduction bug) ===
    //
    // Historically this was a SubsetQ guard that THREW if the
    // Reduce-mode masters file was not ⊆ the Masters-mode sector
    // enumeration, mirroring upstream
    // `If[!SubsetQ[masters, str], Abort["inconsistent masters from Kira"]]`.
    //
    // For EQUAL-mass (and most other) families the two Kira passes
    // (`select_mandatory_recursively` preheat vs `select_mandatory_list`
    // reduce) enumerate the SAME master basis, so the guard never
    // fired.  But for DISTINCT masses (e.g. the {4,1} 4-loop banana,
    // masses^2 = (1,1,1,1,9)) the recursive preheat heuristic
    // genuinely MISSES masters that the actual target reduction needs:
    // the preheat finds 25 sector-31 masters while the reduce pass
    // needs 26 (an extra dot-3 D5 master like j[...,1,1,1,1,4,...] that
    // is genuinely irreducible because the distinct D5 mass breaks the
    // symmetry the preheat heuristic relies on).  This is a real basis
    // disagreement, NOT Kira nondeterminism or a fixable config flag
    // (see checkpoints/amflow_diagnosis.json).
    //
    // The principled fix is to take out.masters = UNION(preheat, reduce)
    // and build the reduction over the consistent superset basis.  Any
    // Reduce-mode master absent from the preheat enumeration is appended
    // to out.masters (preserving the preheat order for the original
    // masters so equal-mass behaviour is byte-identical).  When the
    // reduce masters are already ⊆ preheat (equal-mass et al.), the
    // union adds nothing and this takes EXACTLY the old code path.
    auto reduce_masters = kira_read_masters(fc, reduce_dir);
    {
        std::size_t n_added = 0;
        for (const auto& m : reduce_masters) {
            std::ostringstream key;
            key << m.family();
            for (long v : m.indices()) key << "|" << v;
            if (sector_master_keys.insert(key.str()).second) {
                // Not in the preheat enumeration -> union it in.
                out.masters.push_back(m);
                ++n_added;
            }
        }
        if (n_added > 0) {
            // Re-sort the unioned master set with the same (no-op) hook
            // the preheat path used, so any future non-trivial ordering
            // stays consistent between the two code paths.
            sort_integrals_like_amflow(out.masters);
        }
    }

    auto raw_rules = kira_read_target_table(fc, reduce_dir);

    std::map<std::string, std::size_t> master_index;
    for (std::size_t i = 0; i < out.masters.size(); ++i) {
        std::ostringstream key;
        key << out.masters[i].family();
        for (long v : out.masters[i].indices()) {
            key << "|" << v;
        }
        master_index[key.str()] = i;
    }

    out.rules.reserve(targets.size());
    for (std::size_t i = 0; i < targets.size(); ++i) {
        out.rules.emplace_back();
    }

    auto target_key = [](const qft::JIntegral& j) {
        std::ostringstream k;
        k << j.family();
        for (long v : j.indices()) k << "|" << v;
        return k.str();
    };
    std::map<std::string, std::size_t> target_index;
    for (std::size_t i = 0; i < targets.size(); ++i) {
        target_index[target_key(targets[i])] = i;
    }

    std::size_t n_rhs_union = 0;
    for (const auto& kr : raw_rules) {
        auto it = target_index.find(target_key(kr.lhs));
        if (it == target_index.end()) {
            continue;
        }
        std::size_t i = it->second;
        for (const auto& [coef_str, rhs_j] : kr.rhs) {
            // Mirror Kira/interface.m:486 (AnalyticReduction): every
            // J integral on the RHS of a Kira reduction must be in
            // the master list.  Upstream Aborts with "wrong reduction
            // results from Kira."
            //
            // EXTENSION of the master-basis UNION fix above (PM eikonal
            // boundary family pm4_K3p0_b0_r0_f0): Kira's two
            // passes can disagree one step FURTHER — a target rule's RHS
            // can reference an integral that is in NEITHER masters file
            // (the reduce pass found extra relations killing its sector's
            // preheat masters, but the dangling reference was not in
            // `select_mandatory_list`, so its own rule was never
            // exported).  Such an integral is a well-defined integral of
            // the family; unioning it into the master basis is redundant
            // but correct (downstream boundary setup computes its value
            // directly, exactly as for any other master).  Cap the number
            // of appends so genuinely corrupt reductions still throw.
            if (master_index.find(target_key(rhs_j)) == master_index.end()) {
                constexpr std::size_t kMaxRhsUnion = 16;
                if (n_rhs_union >= kMaxRhsUnion) {
                    std::ostringstream msg;
                    msg << "ibp::black_box_reduce: Kira returned an RHS J "
                           "integral that is not in the master list (target="
                        << target_key(kr.lhs) << ", offending rhs="
                        << target_key(rhs_j)
                        << ") and the RHS-union cap (" << kMaxRhsUnion
                        << ") is exhausted.  This indicates an inconsistent "
                           "Kira reduction.";
                    throw std::runtime_error(msg.str());
                }
                std::cerr << "[amflow][ibp] WARN black_box_reduce: RHS "
                             "integral not in master list; unioning it into "
                             "the master basis (target=" << target_key(kr.lhs)
                          << ", rhs=" << target_key(rhs_j) << ")\n";
                master_index[target_key(rhs_j)] = out.masters.size();
                out.masters.push_back(rhs_j);
                ++n_rhs_union;
            }
            DerivTerm dt;
            dt.coef = parse_kira_coef(cfg, opts, out.red_ctx.ctx,
                                      parse_ctx, coef_str);
            dt.integ = rhs_j;
            out.rules[i].push_back(std::move(dt));
        }
    }

    for (std::size_t i = 0; i < targets.size(); ++i) {
        if (!out.rules[i].empty()) continue;
        auto mit = master_index.find(target_key(targets[i]));
        if (mit != master_index.end()) {
            DerivTerm dt;
            dt.coef = Mfrac::one(out.red_ctx.ctx);
            dt.integ = targets[i];
            out.rules[i].push_back(std::move(dt));
        }
    }

    return out;
}

DiffeqResult
diffeq(const qft::FamilyConfig& fc,
       const std::vector<qft::JIntegral>& jpreferred,
       const std::vector<std::string>& vars,
       const std::vector<int>& top_pattern,
       const ReduceOptions& opts) {
    DiffeqResult out;
    out.red_ctx = make_reduction_context(fc);
    out.vars = vars;

    if (jpreferred.empty()) return out;

    ReduceOptions opts_eff =
        apply_jdot_jrank_floor(opts, {&jpreferred}, true);

    // === Route B: minimal-basis preheat seed override ===
    //
    // The auto-floor above sets the Masters-mode preheat seed to
    // `max_dot(jpreferred)+1` / `max_rank(jpreferred)`.  For the the 3-loop crossed-box
    // crossed-box aux-mass family this floors to (r=11, d=5), and Kira's
    // `select_mandatory_recursively` then enumerates a 3584-integral
    // mandatory closure of the whole (rank≤11, dot≤5) box — which it
    // reduces SYMBOLICALLY over (eta,d) via Fermat (FireFly off), with
    // per-coefficient eta-polynomials up to ~eta^49.  That is the 680GB
    // blow-up: the runaway is the *mandatory-list size driven by the
    // recursive box*, not the genuine master basis (which is the finite
    // 67-master parent set the η-injection leaves invariant).
    //
    // The genuine eta-injected sub-sector master basis is recovered at
    // the minimal seed that still enumerates every jpreferred master.
    // The +1 dot padding (`dot_plus_one=true`) is pure headroom: the
    // jpreferred masters themselves are the basis, so seeding the preheat
    // at exactly `max_dot(jpreferred)` (no +1) keeps the enumeration a
    // superset of jpreferred while collapsing the mandatory closure to
    // the minimal basis (crossed-box top: 22 mandatory at the genuine dot vs
    // 3584 at the padded box).  The `preferred_ ⊆ sortedmasters` guard in
    // AMFSystem::build_diffeq still verifies completeness.
    //
    // `AMFLOW_PREHEAT_DOT` / `AMFLOW_PREHEAT_RANK` pin the preheat seed
    // explicitly (minimal-basis seeding).  DEFAULT is the auto-floor
    // (`max_dot+1`, byte-identical to the prior path) — the +1 padding is
    // retained by default because the Masters-mode `select_mandatory_
    // recursively` enumeration is the one Kira call whose master DETECTION
    // we must not risk truncating (the downstream `preferred_ ⊆
    // sortedmasters` guard would abort).  For the 3-loop crossed-box the
    // minimal seed that still recovers the full 67-master basis is
    // `AMFLOW_PREHEAT_DOT=4` (= max_dot, no padding): the mandatory
    // closure collapses from 3584 (d=5) to a far smaller set while the
    // basis is unchanged.  Pin it per-family from the run config rather
    // than lowering the default for every family.
    {
        auto envl = [](const char* k, long dflt) -> long {
            const char* v = std::getenv(k);
            if (!v || !*v) return dflt;
            try { return std::stol(v); } catch (...) { return dflt; }
        };
        opts_eff.ibp_dot  = envl("AMFLOW_PREHEAT_DOT",  opts_eff.ibp_dot);
        opts_eff.ibp_rank = envl("AMFLOW_PREHEAT_RANK", opts_eff.ibp_rank);
    }

    KiraConfig cfg = make_kira_config(fc, top_pattern, opts_eff);

    // === AMFLOW_DIFFEQ_NUMD activation (diffeq leg ONLY) ===
    //
    // The numeric-d substitution is scoped to ibp::diffeq — the leg
    // whose Reduce-mode FireFly reconstruction runs over symbolic
    // {d, eta} (a measured symbolic-reconstruction wall).  ibp::reduce (boundary /
    // physical reductions) is left byte-identical.  Consistency with
    // the run's eps is enforced (a) here when the caller pinned a
    // numeric `eps`, and (b) at pipeline startup against eps_samples in
    // black_box_amflow_single (the silent-wrong-d killer).
    const std::string diffeq_numd = diffeq_numd_env();
    if (!diffeq_numd.empty()) {
        if (auto it = opts.numeric_values.find("eps");
            it != opts.numeric_values.end()) {
            fmpq_t eps_q, d_q, numd_q;
            fmpq_init(eps_q); fmpq_init(d_q); fmpq_init(numd_q);
            const bool eps_ok =
                fmpq_set_str(eps_q, it->second.c_str(), 10) == 0;
            const bool numd_ok =
                fmpq_set_str(numd_q, diffeq_numd.c_str(), 10) == 0;
            bool match = false;
            if (eps_ok && numd_ok) {
                fmpq_canonicalise(eps_q);
                fmpq_canonicalise(numd_q);
                fmpq_set_si(d_q, 4, 1);
                fmpq_t two_eps; fmpq_init(two_eps);
                fmpq_mul_si(two_eps, eps_q, 2);
                fmpq_sub(d_q, d_q, two_eps);
                fmpq_clear(two_eps);
                match = fmpq_equal(d_q, numd_q) != 0;
            }
            fmpq_clear(eps_q); fmpq_clear(d_q); fmpq_clear(numd_q);
            if (!match) {
                throw std::runtime_error(
                    "AMFLOW_DIFFEQ_NUMD='" + diffeq_numd +
                    "' does not equal 4 - 2*eps for this run's eps='" +
                    it->second + "' (or one of them failed to parse).  "
                    "Refusing to bake a wrong d into the reduction "
                    "(silent-wrong-d killer).");
            }
        }
        cfg.diffeq_numd = diffeq_numd;
    }

    auto parse_ctx = make_parse_ctx(cfg, opts, out.red_ctx.ctx);

    // Mirror upstream `BlackBoxDiffeq` (Kira/interface.m): a two-Kira-call
    // sequence at the same `(rank, dot)`.  Same separate-subdir
    // pattern as in `ibp::reduce` above (see the matching block for
    // the rationale).  Preheat is the sector-wide master enumeration;
    // the second call reduces the libp_deriv-derived integral list
    // against those masters.  Both Kira runs use `opts_eff`'s
    // `(rank, dot)`; re-flooring over `all_ints` for the inner call
    // would break master-count agreement with upstream.
    const std::string masters_dir = opts.work_dir + "/masters_preheat";
    const std::string reduce_dir  = opts.work_dir + "/target_reduce";

    // === Step 1: Masters-mode preheat (sector-wide enumeration) ===
    fs::remove_all(masters_dir);
    fs::create_directories(masters_dir);
    kira_write_config(cfg, masters_dir);
    kira_write_preferred(jpreferred, fc, masters_dir);
    kira_write_jobs(cfg, masters_dir, KiraReductionMode::Masters);
    kira_run_cached(cfg, masters_dir, opts, /*run_split_export=*/false,
                    KiraReductionMode::Masters);

    auto masters = kira_read_masters(fc, masters_dir);
    if (masters.empty()) {
        throw std::runtime_error("ibp::diffeq: no masters from Kira");
    }
    sort_integrals_like_amflow(masters);
    out.sortedmasters = masters;

    auto int_key = [](const qft::JIntegral& j) {
        std::ostringstream k; k << j.family();
        for (long v : j.indices()) k << "|" << v;
        return k.str();
    };

    // Snapshot the preheat master set for the SubsetQ check below.
    std::set<std::string> preheat_master_keys;
    for (const auto& m : masters) preheat_master_keys.insert(int_key(m));

    // Compute libp_deriv per master per requested variable.  Mirrors
    // upstream `der = ComputeDerivative[masters, #]&/@vars`.
    //
    // Two memory/perf optimisations applied here:
    //   (a) reproject each per-term Mfrac to the narrow red_ctx
    //       ({eta, d}) immediately after libp_deriv returns and BEFORE
    //       `simplify_terms` runs.  Numeric substitution alone is not
    //       enough — without reprojection the Mfrac stays on `fc.ctx`
    //       (10+ variables for multi-mass 3L families) and FLINT runs
    //       full-ctx multivariate-GCD on every subsequent operation.
    //   (b) hoist `libp_denoms_deriv` out of the per-master loop: `dd`
    //       depends only on (fc, var), so compute it once per `var`
    //       instead of `masters.size()` times.  Per-master `libp_deriv`
    //       is then inlined using the shared `dd`, mirroring
    //       `src/ibp/libp_deriv.cpp::libp_deriv`'s body, and
    //       `neg_ak_cache[ak]` amortises `Mfrac::from_si` allocations.
    long N_props = (long)fc.propagators_after_conservation.size();

    // Per-var derivative state (`dd`, n_d_local, neg_ak_cache), kept
    // alive past the initial preheat-master loop so that the
    // master-basis UNION fix below (see the matching block after the
    // Reduce-mode SubsetQ check) can compute libp_deriv rows for any
    // masters the Reduce pass discovers but the preheat enumeration
    // missed.  The per-master derivative body is factored into
    // `compute_der_row` so both the preheat masters and the unioned
    // ones go through identical code.
    struct VarDerivState {
        LibpDenomsDerivResult dd;
        long n_d_local = 0;
        std::map<long, Mfrac> neg_ak_cache;
    };
    std::vector<VarDerivState> var_state(vars.size());

    auto compute_der_row =
        [&](VarDerivState& vs, const qft::JIntegral& m) -> std::vector<DerivTerm> {
        if ((long)m.n_indices() != N_props) {
            throw std::invalid_argument(
                "ibp::diffeq: master indices length != #propagators");
        }
        std::vector<DerivTerm> raw;
        for (long k = 0; k < N_props; ++k) {
            long ak = m.indices()[(std::size_t)k];
            if (ak == 0) continue;

            auto cache_it = vs.neg_ak_cache.find(ak);
            if (cache_it == vs.neg_ak_cache.end()) {
                cache_it = vs.neg_ak_cache.emplace(
                    ak, Mfrac::from_si(out.red_ctx.ctx, -ak)).first;
            }
            const Mfrac& neg_ak = cache_it->second;

            for (long jp = 0; jp < vs.n_d_local; ++jp) {
                const Mfrac& cf = vs.dd.coef[(std::size_t)k][(std::size_t)jp];
                if (cf.is_zero()) continue;
                std::vector<long> new_idx = m.indices();
                new_idx[(std::size_t)k] += 1;
                if (jp < N_props) new_idx[(std::size_t)jp] -= 1;
                qft::JIntegral shifted(m.family(), new_idx);
                DerivTerm dt;
                dt.integ = shifted;
                dt.coef = neg_ak.clone();
                dt.coef *= cf;
                raw.push_back(std::move(dt));
            }

            const Mfrac& ck = vs.dd.constant[(std::size_t)k];
            if (!ck.is_zero()) {
                std::vector<long> new_idx = m.indices();
                new_idx[(std::size_t)k] += 1;
                qft::JIntegral shifted(m.family(), new_idx);
                DerivTerm dt;
                dt.integ = shifted;
                dt.coef = neg_ak.clone();
                dt.coef *= ck;
                raw.push_back(std::move(dt));
            }
        }
        return simplify_terms(std::move(raw));
    };

    std::vector<std::vector<std::vector<DerivTerm>>> der(vars.size());
    for (std::size_t v = 0; v < vars.size(); ++v) {
        var_state[v].dd = libp_denoms_deriv(fc, vars[v],
                                            opts.numeric_values,
                                            out.red_ctx.ctx);
        var_state[v].n_d_local = (long)var_state[v].dd.completede.size();

        der[v].reserve(masters.size());
        for (const auto& m : masters) {
            der[v].push_back(compute_der_row(var_state[v], m));
        }
    }

    // Collect the unique integrals appearing in any derivative term,
    // plus the masters themselves (so the Reduce-mode call's target
    // list covers the basis).  Mirrors upstream
    // `integrals = Cases[der, _?(Head[#]===Symbol["j"]&), Infinity] //
    // DeleteDuplicates`.
    std::vector<qft::JIntegral> all_ints;
    std::set<std::string> all_ints_seen;
    {
        for (const auto& vec : der) {
            for (const auto& terms : vec) {
                for (const auto& t : terms) {
                    std::string k = int_key(t.integ);
                    if (all_ints_seen.insert(k).second) all_ints.push_back(t.integ);
                }
            }
        }
        // Do not add `masters` to `all_ints`.  MMA's
        // `DifferentialEquation` (`Kira/interface.m:500`) only sends
        // `Cases[der, j[...], Infinity] // DeleteDuplicates` to
        // `AnalyticReduction`.  Adding preheat masters here inflates
        // the Reduce target list and Kira's memory footprint without
        // gain — identity rules for masters-of-themselves are filled
        // in by the loop below (after raw_rules is parsed).
    }

    // === Master-basis UNION pre-seed (distinct-mass AnalyticReduction
    //     bug) ===
    //
    // For distinct masses the Reduce pass enumerates masters the
    // recursive preheat misses (see the UNION block after the Reduce
    // run for the full rationale).  Those extra masters get appended to
    // the basis, and the DE matrix needs a CORRECT derivative row for
    // each of them — i.e. every integral appearing in their libp_deriv
    // must be reduced by the Reduce pass.  Since the Reduce target list
    // `all_ints` is built only from the preheat masters' derivatives,
    // the extra masters' own derivative integrals would otherwise be
    // unreduced (no rule -> silently-dropped terms -> wrong row).
    //
    // The caller (`AMFSystem::build_diffeq`) passes the parent's full
    // master basis as `jpreferred` — for the {4,1} banana this already
    // includes the extra dot-3 D5 master.  We therefore pre-seed
    // `all_ints` with the derivative integrals of every `jpreferred`
    // master that is NOT already a preheat master, guaranteeing the
    // Reduce pass reduces everything the unioned DE rows reference.
    // For equal-mass (and any family where jpreferred ⊆ preheat) this
    // adds nothing and `all_ints` is byte-identical to the old path.
    for (const auto& jp : jpreferred) {
        if (preheat_master_keys.count(int_key(jp))) continue;
        for (std::size_t v = 0; v < vars.size(); ++v) {
            for (const auto& t : compute_der_row(var_state[v], jp)) {
                std::string k = int_key(t.integ);
                if (all_ints_seen.insert(k).second) all_ints.push_back(t.integ);
            }
        }
    }

    // === Step 2: Reduce-mode (target reduction in a sibling subdir) ===
    // Same `(rank, dot)` as preheat — NO re-floor over `all_ints`.
    // The re-floor used to be applied implicitly by the nested
    // `reduce()` call's own `apply_jdot_jrank_floor(opts, {targets,
    // preferred}, false)`; bypassing the nested call and writing the
    // Reduce yaml directly here preserves `opts_eff`'s `(rank, dot)`.
    //
    // **Preferred file**: mirrors MMA `AnalyticReduction`
    // (Kira/interface.m:458-489), which **does not call `Preferred[]`**
    // -- it reuses the file IBPSystem wrote (the *input* `preferred`).
    // Writing the sorted-preheat output (228 masters) here instead of
    // the original input (172 masters) changes Kira's elimination
    // order and produces a different reduction rule (the dotted
    // self-referencing form with 49+ RHS masters seen on pentabox).
    // That mismatched rule structure cascades through the DE matrix
    // construction and breaks pentabox top-sector master values by
    // ~10^12.  Pentagon/wzbox/hexagon/mercedes/sunset aren't sensitive
    // because their preheat masters happen to coincide with the input
    // preferred ordering closely enough.
    // Re-floor the REDUCE pass over `all_ints` (the actual integrals being
    // reduced), NOT just over `jpreferred`.  The derivative of a dot-d
    // master is a dot-(d+1) integral; the preheat/jpreferred floor
    // (`max dots in jpreferred + 1`) is generally LOWER than the dots
    // present in the derivative target list.  When the Reduce pass is
    // seeded below the dot of a target, Kira cannot express that target in
    // the master basis and instead PROMOTES it to a spurious sub-sector
    // master (e.g. fam[1,1,0,2,1,1], then fam[1,1,0,3,1,1], ...).  Its
    // derivative is then one dot higher again, so the closure loop below
    // would chase an ever-rising dot tower (observed: dot 2,3,...,19,
    // never converging) — the crossed-box sectors 59/55/31 "UNION fix
    // incomplete" failure.  Seeding the Reduce pass at the actual target
    // (rank,dot) lets Kira reduce these dotted sub-sector integrals to the
    // genuine masters, so no spurious master is created and the closure
    // converges.  The preheat (Masters) enumeration is left at the
    // jpreferred floor so the canonical master basis is unchanged for
    // equal-mass families; only the Reduce *seeding* is widened, and the
    // UNION fix already reconciles any master-set difference.
    //
    // === Dot-tower PRE-SEED (the crossed-box sec59/55/31 root fix) ===
    //
    // For the crossed aux-mass families the closure loop below would,
    // with only the `all_ints`-derived floor, raise the Reduce-pass dot
    // by exactly +1 per pass: each pass promotes a dot-(d) sub-sector
    // integral to a spurious master, whose derivative is dot-(d+1), which
    // re-floors the next pass +1 higher, and so on (observed tower
    // r = 9,10,...,21, never closing within kMaxClosurePasses).  The list
    // of "mandatory masters" GROWS monotonically (122->127->132->...) — the
    // signature of under-seeding at depth, not of a genuinely infinite
    // basis.
    //
    // Root fix: seed the Reduce pass at a dot floor RAISED by an explicit
    // headroom up front, so Kira reduces the whole dotted sub-sector tower
    // to its genuine (finite) master set in a SINGLE pass — no spurious
    // promotions, so the closure loop discovers no new extras and the basis
    // closes immediately.  `dot_headroom` is added to BOTH the rank and dot
    // floors (the auxiliary-mass DE can also raise numerator/ISP rank on the
    // crossed sub-sectors).  For equal-mass / planar families the closure
    // loop below never runs its body (`reduce_masters ⊆ preheat` after the
    // first pass), so a non-zero headroom only WIDENS the single Reduce
    // call's (rank,dot) — it cannot change the master basis (the canonical
    // masters are still detected by the preheat at the original floor; the
    // wider Reduce just reduces more dotted targets to that same basis).  We
    // therefore keep `dot_headroom = 0` for the first pass (byte-identical
    // to the prior path when no extra is discovered) and only raise it inside
    // the closure loop, where extras have already been seen.
    auto run_reduce_pass = [&](long dot_headroom) {
        ReduceOptions reduce_opts =
            apply_jdot_jrank_floor(opts_eff, {&all_ints}, false);
        reduce_opts.ibp_dot  += dot_headroom;
        reduce_opts.ibp_rank += dot_headroom;
        KiraConfig reduce_cfg = make_kira_config(fc, top_pattern, reduce_opts);
        // AMFLOW_DIFFEQ_NUMD: propagate the (validated) substitution
        // value onto every Reduce pass of this diffeq leg.  Empty when
        // the env is unset — byte-identical stock behaviour.
        reduce_cfg.diffeq_numd = cfg.diffeq_numd;
        // AMFLOW_KIRA_FFSAVE=1: preserve an input-identical existing dir
        // so a re-fired Kira resumes its FireFly ff_save state — THE
        // monolith wall above (this eta-DE `target_reduce` is the dir the
        // relaunch used to wipe at Done 8857/8946).  A closure pass with a
        // GROWN target list / dot headroom writes different inputs and
        // therefore takes the FRESH path automatically.  Knob OFF
        // (default) = stock remove_all+recreate+write, byte-identical.
        const bool ffsave_resumed = ffsave_prepare_kira_dir(
            reduce_dir, ffsave_fingerprint(reduce_cfg),
            [&](const std::string& d) {
                kira_write_config(reduce_cfg, d);
                kira_write_preferred(jpreferred, fc, d);  // mirror MMA: reuse input preferred
                kira_write_targets(all_ints, fc, d);
                kira_write_jobs(reduce_cfg, d, KiraReductionMode::Reduce);
            });
        // NOTE on the FireFly split (reduce_cfg.split_kira2math): the
        // FireFly reduce commits coefficients to kira.db but emits NO
        // <target>.m in-process.  Without caching we defer the standalone
        // kira2math export to AFTER the closure loop (so a multi-pass
        // closure does not pay the export cost per pass).  WITH caching
        // active, `kira_run_cached` runs the export NOW so the cached
        // entry carries `kira_target.m` — on a later HIT there is no
        // kira.db to export from.  The (rare) extra per-pass export cost
        // is paid once per cache-miss closure pass and never again.
        const bool cache_on = !ibp_cache_resolve_root(opts.ibp_cache_dir).empty();
        kira_run_cached(reduce_cfg, reduce_dir, opts,
                        /*run_split_export=*/cache_on,
                        KiraReductionMode::Reduce, ffsave_resumed);
        return kira_read_masters(fc, reduce_dir);
    };

    auto reduce_masters = run_reduce_pass(/*dot_headroom=*/0);

    // === Master-basis UNION fix (distinct-mass AnalyticReduction bug) ===
    //
    // Historically this was a SubsetQ guard that THREW when a
    // Reduce-mode master was absent from the Masters-mode preheat
    // enumeration, mirroring upstream
    // `If[!SubsetQ[masters, str], Abort["inconsistent masters from Kira"]]`.
    //
    // For DISTINCT masses (e.g. the {4,1} 4-loop banana, masses^2 =
    // (1,1,1,1,9)) the recursive preheat (`select_mandatory_recursively`)
    // genuinely MISSES masters that the target reduction
    // (`select_mandatory_list`) needs — the broken symmetry of the
    // distinct D5 mass makes an extra dot-3 D5 master irreducible.  This
    // is a real basis disagreement, not Kira nondeterminism (see
    // checkpoints/amflow_diagnosis.json).  We therefore take
    // masters = UNION(preheat, reduce) and build the DE over the
    // consistent superset basis.
    //
    // Because the DE matrix is built over `masters` (as both rows and
    // RHS columns) and `der[]` is computed per master, any unioned
    // master must ALSO get its libp_deriv row computed here — `der[]`
    // for the preheat masters was filled above, before this Reduce
    // pass discovered the extras.  We append the extra masters to
    // `masters`/`out.sortedmasters` (preserving the preheat order for
    // the original masters, so equal-mass behaviour is byte-identical),
    // extend `preheat_master_keys`, and call `compute_der_row` for each
    // added master per var.  When the reduce masters are already ⊆
    // preheat (equal-mass et al.), nothing is appended and this takes
    // EXACTLY the old code path (union == preheat => no behavioural
    // change).
    //
    // === Reduce-discovered-extra CLOSURE (the 3-loop crossed-box gap) ===
    //
    // The original pre-seed (above) only added the derivative integrals
    // of the *jpreferred* non-preheat masters to `all_ints`.  But the
    // Reduce pass can itself discover an extra master that is in NEITHER
    // the preheat enumeration NOR `jpreferred` — e.g. for the the 3-loop crossed-box
    // crossed double box, sectors 59/55/31 turn up a Reduce-only dot-2
    // sub-sector master `j[fam,5,1,0,2,1,1,...]` whose own derivative
    // `j[fam,5,1,0,3,1,1,...]` (dot-3) was never on the Reduce target
    // list.  `append_master` then THREW ("UNION fix incomplete") because
    // that derivative integral was unreduced — the DE row would be
    // wrong.  The previous code ran the Reduce pass exactly ONCE, so it
    // could not cover integrals that only become required once the extra
    // master is known.
    //
    // The complete fix is a fixed-point CLOSURE with a DOT-TOWER PRE-SEED.
    // After each Reduce pass, compute the libp_deriv of every newly-
    // discovered extra master (without appending it yet), collect any of
    // its derivative integrals that are not yet in `all_ints`, add them to
    // the Reduce target list, and RE-RUN the Reduce pass — but seed that
    // re-run at a dot floor RAISED BY A GROWING HEADROOM (see the
    // `run_reduce_pass` comment), so the entire dotted sub-sector tower is
    // reduced to its genuine finite master set in O(1) passes instead of
    // climbing one dot per pass.  Each iteration can only add masters in a
    // finite (rank,dot)-bounded set, and `all_ints` grows monotonically, so
    // the loop terminates.  When no extra master is discovered (equal-mass
    // et al.), the loop runs the body zero extra times and the basis is
    // byte-identical to the single-pass path.
    //
    // HONEST CONVERGENCE GUARD: we track the running master count.  A finite
    // closure shows the count STABILIZING (it stops growing once the headroom
    // covers the tower).  If the count keeps growing pass after pass even as
    // the headroom is raised well past the observed tower depth, the basis is
    // genuinely non-finite (pathological aux-mass DE) and we throw a distinct,
    // honest error rather than silently looping or capping at a wrong value.
    {
        // Headroom schedule is overridable via env so the seeding cost vs
        // pass-count tradeoff can be tuned per family without a rebuild
        // (crossed sec59's eta-power tower makes a large single high-seed
        // Kira pass memory-pathological; a gentler step trades more passes
        // for a tractable per-pass IBP system).
        auto envl = [](const char* k, long dflt) -> long {
            const char* v = std::getenv(k);
            if (!v || !*v) return dflt;
            try { return std::stol(v); } catch (...) { return dflt; }
        };
        // Route B: minimal-basis seeding favours a GENTLE headroom climb.
        // The genuine eta-injected master basis is finite and equals the
        // jpreferred set; a Reduce-discovered "extra" is a dotted
        // sub-sector integral that the current seed could not express in
        // that basis, so it needs at most a few more dots of headroom —
        // NOT a jump to a (rank+6, dot+6) box whose `select_mandatory_*`
        // closure is the 3584-integral, eta^49-coefficient, 680GB system.
        // A step of 1 reduces each Reduce pass to the SMALLEST seed that
        // covers one more dot of the tower, so the per-pass IBP system
        // stays minimal; the basis still closes in O(tower-depth) passes
        // (the count STABILISES — see the honest convergence guard below)
        // and the peak per-pass memory is bounded by the genuine basis,
        // not the recursive box.  Overridable via env for tuning.
        const int  kMaxClosurePasses = (int)envl("AMFLOW_CLOSURE_MAX_PASSES", 60);
        const long kHeadroomStep     = envl("AMFLOW_CLOSURE_HEADROOM_STEP", 1);
        const long kMaxHeadroom      = envl("AMFLOW_CLOSURE_HEADROOM_MAX", 36);
        int  pass = 0;
        long headroom = 0;
        std::size_t prev_master_count = reduce_masters.size();
        int  stagnant_growth_passes = 0;     // consecutive passes where count grew at max headroom
        for (;;) {
            // Newly-discovered extras = reduce_masters not yet in the
            // preheat/union set.  (Do NOT mutate preheat_master_keys here;
            // the authoritative append happens after closure.)
            std::vector<qft::JIntegral> new_extras;
            for (const auto& m : reduce_masters) {
                if (!preheat_master_keys.count(int_key(m)))
                    new_extras.push_back(m);
            }
            if (new_extras.empty()) break;

            // Collect derivative integrals of the extras that the current
            // Reduce target list does not yet cover.
            std::vector<qft::JIntegral> missing;
            for (const auto& m : new_extras) {
                for (std::size_t v = 0; v < vars.size(); ++v) {
                    for (const auto& t : compute_der_row(var_state[v], m)) {
                        const std::string k = int_key(t.integ);
                        if (all_ints_seen.insert(k).second) {
                            all_ints.push_back(t.integ);
                            missing.push_back(t.integ);
                        }
                    }
                }
            }

            if (missing.empty()) {
                // Every extra's derivatives are already reduced — closure
                // reached; the extras get appended below.
                break;
            }

            if (++pass > kMaxClosurePasses) {
                std::ostringstream msg;
                msg << "ibp::diffeq: UNION closure did not converge after "
                    << kMaxClosurePasses << " Reduce passes ("
                    << missing.size() << " integrals still uncovered, e.g. "
                    << int_key(missing.front())
                    << "; final master count " << reduce_masters.size()
                    << ", headroom " << headroom
                    << ").  Reduce-discovered master basis is not closing.";
                throw std::runtime_error(msg.str());
            }

            // Pre-seed the dot/rank tower: jump the headroom ahead by
            // kHeadroomStep each pass so one Reduce call reduces the whole
            // sub-sector dot tower instead of climbing +1 dot per pass.
            if (headroom < kMaxHeadroom) {
                headroom += kHeadroomStep;
                if (headroom > kMaxHeadroom) headroom = kMaxHeadroom;
            } else {
                // Already at the seeding ceiling.  If the master count is
                // STILL growing pass-on-pass here, the basis is genuinely
                // non-finite for this family — bail out honestly.
                if (reduce_masters.size() > prev_master_count) {
                    if (++stagnant_growth_passes >= 3) {
                        std::ostringstream msg;
                        msg << "ibp::diffeq: UNION closure master basis "
                               "does NOT stabilize even at seeding headroom "
                            << headroom << " (master count still growing: "
                            << prev_master_count << " -> "
                            << reduce_masters.size()
                            << " over " << stagnant_growth_passes
                            << " passes at the ceiling).  The auxiliary-mass "
                               "DE for this family has an apparently non-finite "
                               "(pathological) dotted master tower; no finite "
                               "fixed point exists in this reducer.";
                        throw std::runtime_error(msg.str());
                    }
                } else {
                    stagnant_growth_passes = 0;
                }
            }
            prev_master_count = reduce_masters.size();

            // Re-run Reduce over the enlarged target list, seeded with the
            // raised headroom so the dotted sub-sector tower is pre-covered;
            // this can both reduce the previously-missing integrals AND
            // surface further extra masters, which the next iteration handles.
            reduce_masters = run_reduce_pass(headroom);
        }
    }

    // Append a master to the union basis: extend `masters`/
    // `out.sortedmasters`, compute its libp_deriv DE row per var, and
    // verify every derivative integral was a Reduce target (else the
    // DE-row build below would silently drop the unreduced term and
    // produce a wrong row — fail loudly instead).  After the closure loop
    // above this guard can no longer fire for Reduce-discovered extras;
    // it remains as a correctness invariant.  Shared by both the
    // Reduce-discovered union and the `jpreferred` union below so the
    // two go through identical (audited) code.
    auto append_master = [&](const qft::JIntegral& m) {
        masters.push_back(m);
        out.sortedmasters.push_back(m);
        for (std::size_t v = 0; v < vars.size(); ++v) {
            std::vector<DerivTerm> row = compute_der_row(var_state[v], m);
            for (const auto& t : row) {
                const std::string k = int_key(t.integ);
                if (!all_ints_seen.count(k)) {
                    throw std::runtime_error(
                        "ibp::diffeq: UNION fix incomplete — appended "
                        "master '" + m.to_string() + "' has a "
                        "derivative integral (" + k + ") that was not "
                        "reduced (absent from the Reduce target list). "
                        "The extra master is not covered by jpreferred; "
                        "the DE row would be incorrect.");
                }
            }
            der[v].push_back(std::move(row));
        }
    };

    {
        std::vector<qft::JIntegral> added;
        for (const auto& m : reduce_masters) {
            if (preheat_master_keys.insert(int_key(m)).second) {
                added.push_back(m);
            }
        }
        if (!added.empty()) {
            // Re-sort with the same (no-op) hook the preheat path used.
            sort_integrals_like_amflow(added);
            for (const auto& m : added) append_master(m);
        }
    }

    // === jpreferred master-basis UNION (distinct-mass amfsystem-assembly
    //     SubsetQ) ===
    //
    // The caller `AMFSystem::build_diffeq` enforces `preferred_ ⊆
    // sortedmasters` (mirrors MMA `AMFSystemCombineSolution`'s SubsetQ);
    // `solve_one_eps` then extracts each preferred value from the solved
    // DE solution at `pref_to_sorted_[i]`.  So every `jpreferred` master
    // is by definition a row of the DE and MUST appear in
    // `out.sortedmasters`.
    //
    // For distinct masses the parent's BlackBoxReduce master list
    // (`jpreferred`) can contain a master — e.g. the {4,1} banana
    // numerator/ISP master `[1,1,1,1,1,0,0,0,0,-1,0,-1,0,0]` — that is
    // present in NEITHER the eta-injected preheat enumeration NOR the
    // Reduce pass's discovered masters.  The reduce_masters union above
    // therefore never adds it, and the amfsystem SubsetQ guard fires.
    // (It is a genuinely different master from the dot-3 D5 master the
    // reduce_masters union targets: this one carries negative ISP
    // indices.)
    //
    // A `jpreferred` master is a master by construction, so we union it
    // into the basis with its own libp_deriv DE row.  The pre-seed block
    // above already added the derivative integrals of every non-preheat
    // `jpreferred` master to the Reduce target list (`all_ints`), so each
    // appended row's terms are reduced and the `append_master` guard
    // holds.  For equal-mass (and any family where jpreferred ⊆ preheat)
    // nothing is appended here — `preheat_master_keys` already contains
    // every jpreferred — so this is a strict no-op and the basis is
    // byte-identical to the pre-patch path.
    {
        std::vector<qft::JIntegral> added;
        for (const auto& jp : jpreferred) {
            if (preheat_master_keys.insert(int_key(jp)).second) {
                added.push_back(jp);
            }
        }
        if (!added.empty()) {
            sort_integrals_like_amflow(added);
            for (const auto& m : added) append_master(m);
        }
    }

    // Kira 3.1 FireFly -> kira2math export split (the the 3-loop crossed-box export fault
    // fix; see KiraConfig::split_kira2math).  The final Reduce pass
    // committed its reconstructed coefficients to reduce_dir/kira.db but
    // — under the split — emitted no `results/<fam>/<target>.m`.  Export
    // it now in a FRESH kira process (sharing the committed kira.db) so
    // `kira_read_target_table` below has a populated .m.  A fresh process
    // dodges the in-process kira2math-after-firefly SIG6/bad_alloc.  Run
    // once, on the final reduce_dir state only.
    //
    // Under the IBP cache, `run_reduce_pass` already exported (or restored
    // from cache) `kira_target.m`; skip the export when it is already
    // present and non-empty so a cache HIT does not try to kira2math from
    // a non-existent kira.db.
    if (cfg.split_kira2math) {
        std::error_code ec;
        fs::path tgt_m = fs::path(reduce_dir) / "results" / fc.family
                       / "kira_target.m";
        if (!fs::exists(tgt_m, ec) || fs::file_size(tgt_m, ec) == 0) {
            kira_run(cfg, reduce_dir, opts.log_file, "jobs_kira2math.yaml");
        }
    }

    auto raw_rules = kira_read_target_table(fc, reduce_dir);

    std::map<std::string, std::size_t> master_index;
    for (std::size_t i = 0; i < masters.size(); ++i) {
        master_index[int_key(masters[i])] = i;
    }

    // Parse `raw_rules` (target_table.m) into a lhs-keyed map of
    // DerivTerm lists.  Identity rule for any `all_ints` element that
    // is itself a master and didn't appear as an explicit rule LHS.
    std::map<std::string, std::vector<DerivTerm>> rule_by_lhs;
    for (const auto& kr : raw_rules) {
        std::vector<DerivTerm> dts;
        dts.reserve(kr.rhs.size());
        for (const auto& [coef_str, rhs_j] : kr.rhs) {
            // Mirror Kira/interface.m:486 (AnalyticReduction RHS-in-master
            // guard); also the C++ D6 fix.  Now that the SubsetQ guard
            // above passes, RHS-not-in-master can only fire on a Kira
            // internal inconsistency.
            if (master_index.find(int_key(rhs_j)) == master_index.end()) {
                std::ostringstream msg;
                msg << "ibp::diffeq: Kira returned an RHS J integral "
                       "that is not in the master list (target="
                    << int_key(kr.lhs) << ", offending rhs="
                    << int_key(rhs_j)
                    << ").  This indicates an inconsistent Kira reduction.";
                throw std::runtime_error(msg.str());
            }
            DerivTerm dt;
            dt.coef = parse_kira_coef(cfg, opts, out.red_ctx.ctx,
                                      parse_ctx, coef_str);
            dt.integ = rhs_j;
            dts.push_back(std::move(dt));
        }
        rule_by_lhs[int_key(kr.lhs)] = std::move(dts);
    }
    for (const auto& integ : all_ints) {
        auto k = int_key(integ);
        if (rule_by_lhs.count(k)) continue;
        auto mit = master_index.find(k);
        if (mit != master_index.end()) {
            DerivTerm dt;
            dt.coef = Mfrac::one(out.red_ctx.ctx);
            dt.integ = integ;
            // `vector<DerivTerm>{std::move(dt)}` would go through
            // initializer_list which copies — DerivTerm holds a
            // move-only Mfrac, so the copy is deleted.  Build the
            // single-element vector via push_back instead.
            std::vector<DerivTerm> identity_rule;
            identity_rule.push_back(std::move(dt));
            rule_by_lhs[k] = std::move(identity_rule);
        }
    }

    // Batched lcm-sum per matrix entry, mirroring MMA's
    // `Together[der/.j->red]` at `Kira/interface.m:505`.  Pairwise
    // `out.diffeq[v][row][col] += product` would canonicalise each
    // addition (every step ~3× operand-size temporaries inside FLINT);
    // collecting all (num_i, den_i) pairs first, computing
    // `L = lcm(den_i)` once, scaling, summing, then canonicalising
    // once cuts that to a single canonicalisation per entry.
    auto batched_sum = [](std::vector<Mfrac>& terms,
                          const std::shared_ptr<MpolyContext>& ctx) -> Mfrac {
        if (terms.empty()) return Mfrac::zero(ctx);
        if (terms.size() == 1) return std::move(terms.front());
        // L = lcm of all denominators.
        // lcm(a, b) = a / gcd(a, b) * b  (uses exact_divide for the
        // division step, since gcd divides a by definition).
        Mpoly L = terms[0].denominator();
        for (std::size_t i = 1; i < terms.size(); ++i) {
            Mpoly g = Mpoly::gcd(L, terms[i].denominator());
            Mpoly L_div_g(ctx);
            if (!Mpoly::exact_divide(L_div_g, L, g)) {
                throw std::runtime_error(
                    "ibp::diffeq batched_sum: gcd did not divide lcm "
                    "(should be impossible)");
            }
            L = L_div_g * terms[i].denominator();
        }
        // N_total = sum_i (num_i * (L / den_i)).
        Mpoly N_total = Mpoly::zero(ctx);
        for (auto& t : terms) {
            Mpoly scale(ctx);
            if (!Mpoly::exact_divide(scale, L, t.denominator())) {
                throw std::runtime_error(
                    "ibp::diffeq batched_sum: term denominator does not "
                    "divide lcm (should be impossible)");
            }
            N_total += t.numerator() * scale;
        }
        return Mfrac(std::move(N_total), std::move(L));  // canonicalises once
    };

    out.diffeq.reserve(vars.size());
    for (std::size_t v = 0; v < vars.size(); ++v) {
        out.diffeq.emplace_back();
        out.diffeq[v].reserve(masters.size());
        for (std::size_t row = 0; row < masters.size(); ++row) {
            // Collect all products per column for this row, deferring
            // canonicalisation to one batched_sum call per column.
            std::map<std::size_t, std::vector<Mfrac>> per_col_products;
            for (const auto& dt : der[v][row]) {
                auto rit = rule_by_lhs.find(int_key(dt.integ));
                if (rit == rule_by_lhs.end()) continue;
                for (const auto& mt : rit->second) {
                    auto mit = master_index.find(int_key(mt.integ));
                    if (mit == master_index.end()) continue;
                    Mfrac product = dt.coef.clone();
                    product *= mt.coef;
                    per_col_products[mit->second].push_back(std::move(product));
                }
            }
            out.diffeq[v].emplace_back();
            out.diffeq[v][row].reserve(masters.size());
            for (std::size_t col = 0; col < masters.size(); ++col) {
                auto it = per_col_products.find(col);
                if (it == per_col_products.end()) {
                    out.diffeq[v][row].push_back(Mfrac::zero(out.red_ctx.ctx));
                } else {
                    out.diffeq[v][row].push_back(
                        batched_sum(it->second, out.red_ctx.ctx));
                }
            }
        }
    }

    return out;
}

}  // namespace amflow::ibp
