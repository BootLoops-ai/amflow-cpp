// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Anthropic, PBC. Added by the BootLoops fork of AMFlow.cpp; see NOTICE and PATCHES.md.
// pipeline::system_checkpoint — persist completed per-system DE solves.
//
// PROBLEM: a multi-hour 3-loop η-flow (15-propagator self-energies, bananas, …) walks a tree of
// `AMFSystem` nodes; each non-ending node solves its η-deformed differential
// equation at every ε-sample (`AMFSystem::solve` → `solve_one_eps` →
// `ode::amflow`), and each ending node looks up / Kira-reduces its boundary.
// When the process dies mid-flow (an aborted reduction, a stall, a machine
// restart) the
// WHOLE flow restarts from scratch — every already-completed system is
// re-solved.  The Kira *reduction* is already protected by the content-
// addressed `ibp_cache` (survives a crash, re-keys correctly); the expensive
// piece that was NOT persisted is the per-system numeric solve — the master
// values on the ε-grid.
//
// FIX: a content-addressed checkpoint, one file per completed AMFSystem.
// Before a node solves, it computes a stable key from its identity
// (family + preferred masters + etac + ε-grid + working precision +
// direction + ending flag) and looks for `<root>/system_<key>.ckpt`.  On a
// HIT the node's `solutions_` (master values per ε, plus the SingleMass /
// Cutkosky `global_*` remap) are restored losslessly and the solve is
// SKIPPED — including the recursion into its sub-systems, since a node's
// solution already incorporates them.  On a MISS the node solves normally
// and, on success, writes the checkpoint with a write-then-rename so a
// concurrent / interrupted run never reads a half-written file.
//
// LOSSLESS: acb balls are serialized via FLINT `arb_dump_str` /
// `arb_load_str` (exact midpoint+radius, base-16), so a resumed flow is
// bit-identical to an uninterrupted one — the checkpoint is a cache, never
// an approximation.
//
// The checkpoint root is resolved (in priority order) from:
//   1. env `AMFLOW_SYSTEM_CHECKPOINT` (explicit opt-in directory);
//   2. `<ibp_cache_dir>/syssolve`     — co-located with the reduction cache
//      (which already persists across runs and re-keys on config), so
//      enabling `AMFLOW_IBP_CACHE` transparently also resumes solves;
//   3. empty → checkpointing disabled (byte-identical to the pre-checkpoint
//      path).
//
// Disabled by default unless one of those roots resolves non-empty.

#ifndef AMFLOW_PIPELINE_SYSTEM_CHECKPOINT_HPP
#define AMFLOW_PIPELINE_SYSTEM_CHECKPOINT_HPP

#include <string>
#include <vector>

#include "amflow/numeric/acb_wrap.hpp"

namespace amflow::pipeline {

struct AMFSystemSolution;  // fwd (defined in amfsystem.hpp)

// Resolve the checkpoint root.  `ibp_cache_dir` is the resolved ibp cache
// root (may be empty).  Returns "" when checkpointing is disabled.
std::string system_checkpoint_resolve_root(const std::string& ibp_cache_dir);

// Build the stable per-system content key (16 hex chars).  Folds in the
// family name, the preferred-master keys (order-preserving — the solution
// vector is indexed by this order), the η-coefficient vector, the ε-grid
// (lossless), the working precision in bits, the path direction, the
// ending flag, and — read from the process state at call time — every
// numeric option (chop_pre, rationalize_pre, d0, the expansion orders, the
// running options) and the solve-affecting environment overrides listed in
// kSolveEnvKnobs (src/pipeline/system_checkpoint.cpp: AMFLOW_DIFFEQ_NUMD,
// AMFLOW_SPARSE_CHOP_DIGITS, AMFLOW_NO_SPARSE_CHOP, AMFLOW_RESONANCE_CHOP,
// AMFLOW_X00_PIVOT_GUARD, AMFLOW_FORCE_ENDING_DEPTH, AMFLOW_FIXEDPOINT_PROBE,
// ...; "" when unset).  Debug dumps, parallelism and eps-grid shaping
// variables are not folded (they do not change a system's solution).  Two
// invocations of the same config produce the same key; a change to any
// folded input changes the key.
std::string system_checkpoint_key(const std::string& family,
                                   const std::vector<std::string>& preferred_keys,
                                   const std::vector<int>& etac,
                                   const std::vector<numeric::AcbValue>& epslist,
                                   long working_prec_bits,
                                   int direction,
                                   bool is_ending);

// Restore `out` from `<root>/system_<key>.ckpt` if it exists and parses.
// Returns true on a HIT (out fully populated).  Logs a HIT line on stderr.
// A malformed / truncated file (e.g. from a kill mid-write of an older
// format) is treated as a MISS, not an error.
bool system_checkpoint_load(const std::string& root,
                            const std::string& key,
                            std::vector<AMFSystemSolution>& out);

// Persist `sol` to `<root>/system_<key>.ckpt` via write-then-rename.
// No-op when `root` is empty.  Logs a STORED line on stderr.
void system_checkpoint_store(const std::string& root,
                             const std::string& key,
                             const std::vector<AMFSystemSolution>& sol);

}  // namespace amflow::pipeline

#endif  // AMFLOW_PIPELINE_SYSTEM_CHECKPOINT_HPP
