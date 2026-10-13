// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Anthropic, PBC. Added by the BootLoops fork of AMFlow.cpp; see NOTICE and PATCHES.md.
// ibp::ibp_cache — content-addressed cache for Kira IBP-reduction results.
//
// PROBLEM: `solve_integrals` (via the AMFSystem recursion) re-runs the
// FULL Kira IBP reduction for every (kinematic point × ε-sample × sub-
// system), even though the symbolic reduction is identical — only the
// substituted numeric kinematics differ, and those are baked into the
// Kira input YAML before the reducer is launched.  For a 5-point production
// family this is ~4 h/point; a 79-point oracle farm is multi-week.
//
// FIX (two layers, both keyed on a content hash of the Kira *input*
// files in a working directory):
//
//   * Layer A — work_dir results cache.  Hash {integralfamilies.yaml,
//     kinematics.yaml, jobs.yaml, jobs_kira2math.yaml, preferred,
//     target} → 16-hex key.  Under `<cache_root>/ibp/<key>/` store the
//     Kira `results/<fam>/` outputs (`masters`, `kira_target.m`).
//     Before forking Kira, look up the key; on HIT, copy the cached
//     results into the work_dir and SKIP the kira_run.  On MISS, run
//     Kira and STORE.  This makes RESTARTS free and de-duplicates the
//     per-ε reductions (same kinematic point, multiple ε values
//     currently re-reduce because solve_integrals wipes the work_dir
//     per call).
//
//   * Layer B — symbolic-kinematics reduction.  When `symbolic_ibp` is
//     on (see ReduceOptions / KiraConfig), `kira_yaml.cpp` writes the
//     kinematics.yaml WITHOUT numeric substitution, so the input hash
//     is identical for EVERY kinematic point and Layer A's cache
//     becomes a single shared reduction across the whole oracle farm.
//     The Kira output then carries rational functions of {d, η, s12,
//     s23, …}; `reduce.cpp` parses them on a wide MpolyContext and
//     `substitute_and_narrow`s the per-point numeric_values back to
//     {η, d}.  Symbolic-kinematics Kira is slower per-reduction (more
//     variables → bigger Fermat/FireFly) but runs ONCE.  Gated behind
//     `AMFLOW_SYMBOLIC_IBP=1` / `"symbolic_ibp": true` because some
//     families blow up symbolically.
//
// The cache root is resolved (in priority order) from:
//   1. `ReduceOptions::ibp_cache_dir` (surfaced as
//      `amf_options.blackbox.ibp_cache_dir` in the JSON);
//   2. env `AMFLOW_IBP_CACHE`;
//   3. empty → caching disabled (byte-identical to the pre-cache path).

#ifndef AMFLOW_IBP_IBP_CACHE_HPP
#define AMFLOW_IBP_IBP_CACHE_HPP

#include <string>

namespace amflow::ibp {

// Resolve the cache root: explicit `opt_root` if non-empty, else
// env `AMFLOW_IBP_CACHE`, else "".  Empty return ⇒ caching disabled.
std::string ibp_cache_resolve_root(const std::string& opt_root);

// FNV-1a 64-bit hash over the canonical Kira input files in `dir`
// (see header comment for the file list).  Missing files contribute
// nothing; the file *name* is folded into the hash so `target` vs
// no-`target` (Masters vs Reduce) keys diverge even when the other
// inputs coincide.  `jobs.yaml` is canonicalized before hashing (see
// `canonicalize_jobs_yaml`) so result-neutral env knobs (DCAP,
// n_thread) don't perturb the key.  Returns a 16-hex-char lowercase
// string.
//
// `extra_salt` (default "" = byte-identical legacy key) is folded into
// the hash AFTER the input files.  Used by the AMFLOW_DIFFEQ_NUMD
// numeric-d substitution (KiraConfig::diffeq_numd): the substitution
// rewrites tmp/SYSTEM files that are NOT hashed inputs, so without the
// salt a numeric-d reduction would collide with the symbolic-d entry
// for the same inputs (the eps-independent-key collision
// class) and poison the cache in both directions.
std::string ibp_cache_key(const std::string& dir,
                          const std::string& extra_salt = "");

// Canonicalize jobs.yaml content for cache-key hashing: strip the
// `, d: <int>` search-bound cap from the `reduce:` seed-spec line
// (the `AMFLOW_REDUCE_DCAP` pitfall) and any `n_thread`/`threads:`
// parallelism lines.  The `d:` inside `select_mandatory_recursively`
// (Masters mode, RESULT-determining) is preserved.  Exposed for unit
// testing; `ibp_cache_key` calls this internally.
std::string canonicalize_jobs_yaml(const std::string& raw);

// Look up `<cache_root>/ibp/<key>/results/<family>/` and, if it
// contains a non-empty `masters` and/or `kira_target.m`, copy them
// into `<dir>/results/<family>/`.  Returns true on HIT (something was
// restored).  Logs a `[amflow][ibp-cache] HIT …` line on stderr.
bool ibp_cache_lookup(const std::string& cache_root,
                      const std::string& key,
                      const std::string& family,
                      const std::string& dir);

// Copy `<dir>/results/<family>/{masters,kira_target.m}` into the cache
// at `<cache_root>/ibp/<key>/results/<family>/`.  Uses a write-then-
// rename so a concurrent lookup never sees a half-written entry.  Logs
// `[amflow][ibp-cache] STORED …`.  No-op (with a warning) if `dir` has
// no results to store.
void ibp_cache_store(const std::string& cache_root,
                     const std::string& key,
                     const std::string& family,
                     const std::string& dir);

// === Stampede lock ===
//
// When N farm jobs launched together all reach the same (s,t)-
// independent IBP stage, all N MISS the cache simultaneously and all N
// run kira — N−1 wasted (measured on a 14-cold-job farm: jobs × 10
// η-only boundary stages × 5-30 s ≈ 1 h aggregate waste).  This RAII
// guard takes an advisory `flock()` on `<cache_root>/ibp/<key>.lock`
// around the MISS→run-kira→STORE path so the first job computes and
// the rest block, then re-check the cache → HIT, skip kira.
//
// Construction acquires the lock (LOCK_EX).  A non-blocking probe
// (LOCK_NB) is tried first; if it fails the lock logs
//   `[amflow][ibp-cache] WAIT key=... (another job computing)`
// and polls until acquired or `AMFLOW_IBP_LOCK_TIMEOUT_S` (default
// 7200) elapses, after which it proceeds anyway (a crashed holder must
// not deadlock the farm).  `waited()` reports whether the constructor
// had to wait — caller uses it to decide whether to re-check the cache.
// Destruction releases.  Gated behind `AMFLOW_IBP_CACHE_LOCK` (default
// ON; set to 0 to disable).  When disabled or `cache_root` empty, the
// guard is a no-op.
class IbpCacheLock {
  public:
    IbpCacheLock(const std::string& cache_root, const std::string& key);
    ~IbpCacheLock();
    IbpCacheLock(const IbpCacheLock&) = delete;
    IbpCacheLock& operator=(const IbpCacheLock&) = delete;

    // True iff the constructor had to wait on another holder (i.e. the
    // initial LOCK_NB probe failed).  Caller should re-lookup the cache
    // when this is true — the entry was very likely just populated.
    bool waited() const noexcept { return waited_; }

    // True iff the lock is actually held (false when disabled, when
    // cache_root is empty, or after a timeout fallthrough).
    bool held() const noexcept { return fd_ >= 0 && held_; }

    // Is the stampede lock enabled?  Default ON; AMFLOW_IBP_CACHE_LOCK=0
    // disables.
    static bool enabled();

  private:
    int  fd_     = -1;
    bool held_   = false;
    bool waited_ = false;
};

}  // namespace amflow::ibp

#endif  // AMFLOW_IBP_IBP_CACHE_HPP
