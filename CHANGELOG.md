# Changelog

> 中文版本: [CHANGELOG_zh.md](CHANGELOG_zh.md)

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added
- Oracle suite at 228 triplets across L = 1, 2, 3, 4
  (1L=52, 2L=53, 3L=64, 4L=59).  Sub-master reuse of parent MMA
  caches enables cheap 4L diversity.
- Kira yaml emission aligned with MMA `interface.m:126`
  (`[Propagator, 0]` form, closed-form rank-1 propagator
  reconstruction, raw caller-order preferred list).
- Opt-in `AMFLOW_STAGE_TIMING` env var: emits `[kira_time]` per
  Kira call and `[total_time]` at amflow_cli exit, lets external
  tools split wallclock into Kira-subprocess vs amflow-internal
  time.
- CMake install rules now include `amflow_cli`: after
  `sudo cmake --install build` the binary lands at
  `${CMAKE_INSTALL_PREFIX}/bin/amflow_cli` (default `/usr/local/bin/`)
  alongside the library and headers, so the CLI is reachable from
  any working directory without needing to invoke `./build/src/cli/`.

### Fixed
- Pentabox 2L 5-leg Jordan eigenvector normalization:
  `fmpq_mat_nullspace_exact` rescales basis vectors so the last
  non-zero entry is 1, matching Mathematica's convention.
- `build_boundary` projection on multi-mass topologies:
  `to_complete_explicit` rank-filter no longer drops mass-bearing
  propagators; numeric substitution restored.
- Multi-mass 3L corner-master divergence: `analyze_block`
  extend / Gather / Complement now match MMA; default
  `sparse_chop_digits = max(chop_pre, working_pre - 40)` keeps acb
  noise budget below legitimate values.

### Changed
- Complex-valued numeric kinematics is out-of-scope (previously
  "deferred"); JSON dispatcher rejection is unchanged.


## [1.1.0] — 2026-05-12

Post-v1.0 audit-driven correctness pass + Phase 3 oracle diversity
expansion.  Together they close all silent-wrong-result paths
surfaced by the post-release line-level MMA parity audit
([`docs/AUDIT_MMA_PARITY.md`](docs/AUDIT_MMA_PARITY.md)) and open
the project's five "we-didn't-think-of-this" oracle diversity axes.

After v1.1.0 the audit table reads
**86 🟢 / 0 🟡 / 7 🔴 (6 fixed + 1 D5 deferred indefinitely) /
17 ⚪**, and the oracle suite covers all 5 diversity axes
(loop number L ≤ 4, ≥ 3 kinematic invariants, multi-cut Cutkosky,
mixed-mass, and ε-extremes) with rel ≤ ~10⁻³⁰ across the new
oracles.

### Added — Phase 3 oracle diversity expansion
- 3.F **L=4 banana oracle** (`tools/bench/banana_4loop_eps001_*`):
  4-loop equal-mass banana sunrise at `psq=-3, msq=1, eps=1/1000`.
  Surfaced and locked audit divergence **D7** (see Fixed below).
  Matches MMA at rel 7.2 × 10⁻³¹ / 1.6 × 10⁻³⁰ (C++ 399 s vs MMA
  385 s).
- 3.G **multi-invariant electroweak box oracle**
  (`tools/bench/ewbox_1loop_eps001_*`): 1-loop electroweak box with
  alternating W/Z masses and 4 distinct kinematic invariants
  {`s, t, mWsq, mZsq`}.  Matches at rel ≤ 3.8 × 10⁻³⁰ (C++ 39 s
  vs MMA 62 s).  WW bubble above threshold gives physical imag
  part, locking cross-mass branch-cut handling.
- 3.H **multi-cut Cutkosky oracle** (`tools/bench/cutbanana_4L_eps001_*`):
  4-loop massless cutbanana with all 5 internal lines on-shell
  (5-particle Cutkosky cut).  Matches at rel ≤ 2.2 × 10⁻³⁰
  (C++ 141 s vs MMA 149 s).  5-level subsystem recursion verified
  end-to-end.
- 3.I **mixed-mass oracle** (`tools/bench/bn3mix_eps001_*`):
  3-loop banana with 1 W-massive + 3 massless internal propagators
  — first 3-loop mixed-mass case.  Matches at rel ≤ 2.9 × 10⁻³⁰
  (C++ 56 s vs MMA 96 s).  Exercises scaleless-sub-sector
  detection under mixed mass.
- 3.J **ε-extremes oracles** (`tools/bench/cutbubble_1L_eps{2,10000}_*`):
  cutbubble at eps = 1/2 (D = 3) matches at rel 1.9 × 10⁻⁶⁴ on
  the exact value 1/8; at eps = 10⁻⁴ matches at rel 1.0 × 10⁻³².
  An initial run at eps = 1 (D = 2) revealed an upstream MMA
  AMFlow limitation (DESolver returns partially-symbolic
  `(1/2π) Im[DESolver\`Private\`variables[1, 1]]`); eps = 1/2 was
  used as the next-most-extreme rational that yields a clean
  numeric — documented in the bench source comment.

### Fixed (alignment with upstream)
- `RunningOptions::run_length` default raised from 200 to **1000** to
  match upstream `RunLength = 1000` (`AMFlow.m:259`,
  `DESolver.m:97`).  Prevents premature `RunUnit` aborts on crowded
  pole landscapes.  (Audit D1.)
- `GlobalOptions::rationalize_pre` default raised from 20 to **100**
  to match upstream `RationalizePre = 100`.  Removes the silent
  precision-narrowing on rationalization steps along the contour.
  (Audit D2.)
- `ibp::black_box_reduce` now throws `std::runtime_error` if Kira's
  reduction returns an RHS J-integral that is not in the master
  list, mirroring upstream `Kira/interface.m:486` `Abort`.  Prior
  behaviour silently dropped the offending row.  (Audit D6.)
- **D3 — Tradition-with-cut boundary projection.**  v1.0 silently
  dropped the parent system's `Cut` when constructing the boundary
  sub-family.  The Phase 1A patch converted that to a loud abort.
  **Phase 1B implements the full projection** mirroring upstream
  `ReduceBoundary` (`AMFlow.m:790-803`): each parent cut prop is
  passed through the bare `region.transform.map`; each fam.prop is
  matched against the transformed cut props modulo
  `reduced_replacement`; the new `sub_cut` is built and asserted to
  preserve `Count[cut, 1]`; the result is passed to the
  `qft::FamilyConfig::build` of the sub-family.  Backed by a new
  oracle benchmark (`tradcut_phase_2L_eps001_*`) — 2-loop
  Tradition-with-cut from upstream `examples/automatic_phasespace`,
  matches MMA at relative error 2.68 × 10⁻³⁰.  (Audit D3.)
- `branch_to_loop` (`src/qft/region.cpp`) gains a defensive assert
  that `det(A) = ±1` (constant), where `A` is the loop-redefinition
  matrix.  This catches any future relaxation of the `branch_momenta`
  unit-leading-loop-coefficient precondition before the missing
  `|Det|^(4-2eps)` Jacobian factor (upstream `AMFlow.m:731-732`)
  silently produces wrong boundary integrands.  (Audit D4.)
- **D7 — dual-Kira-call master-count divergence.**  `ibp::reduce`
  and `ibp::diffeq` now mirror upstream's `BlackBoxReduce` /
  `BlackBoxDiffeq` two-call pattern: a Masters-mode preheat call
  (sector-wide enumeration via `select_mandatory_recursively`)
  followed by a Reduce-mode call (`select_mandatory_list` for the
  specific targets), both at the same `(rank, dot)`.  Each call
  runs in its own subdirectory (`<work_dir>/masters_preheat/` and
  `<work_dir>/target_reduce/`) because the Kira 2.x release
  refuses to share a `$ReductionDirectory` between the two calls
  (Masters-mode `run_initiate: masters` doesn't register the `-s`
  numeric substitutions, and the subsequent Reduce-mode call
  aborts with `Kira::update_auxiliary_file: Last Kira run set 0
  variables to numeric values, this time you request N`).
  `ibp::diffeq` also no longer calls nested `reduce()` (which
  would re-floor `(rank, dot)` over the derivative integrals — the
  L=4 banana failure mode); it inlines the Reduce-mode Kira
  invocation at `opts_eff`'s `(rank, dot)`, mirroring upstream's
  `AnalyticReduction` which inherits `IBPRank`/`IBPDot` globals
  from `IBPSystem`.  Both functions add a SubsetQ guard on the
  Reduce-mode master file against the Masters-mode sector
  enumeration, mirroring upstream's
  `If[!SubsetQ[masters, str], Abort["inconsistent masters from
  Kira"]]`.  Locked by the L=4 banana oracle.  (Audit D7.)

### Added
- [`docs/AUDIT_MMA_PARITY.md`](docs/AUDIT_MMA_PARITY.md): line-level
  upstream-parity audit report.  86 🟢 verified / 0 🟡 unverified /
  7 🔴 (6 fully fixed, 1 D5 deferred indefinitely) / 17 ⚪ not
  ported.  (Started at 65 🟢 / 21 🟡 / 6 🔴 in v1.0.0; the post-v1.0
  pass closed every 🟡, fixed 6 of the 7 🔴, and surfaced D7 as
  the seventh 🔴 — now fixed too.)
- [`docs/ROADMAP.md`](docs/ROADMAP.md) §"Performance benchmarking":
  wall-clock baseline vs MMA on the 12 oracle benchmarks.  Median
  speedup ~2× on Kira-light benches; converges to ~1× on Kira-heavy
  benches.
- [`docs/ROADMAP.md`](docs/ROADMAP.md): post-v1.0 development plan
  (Phase 1 implementation completeness, Phase 2 oracle expansion,
  Phase 3 ongoing diversification).
- [`tools/bench/run_perf_audit.sh`](tools/bench/run_perf_audit.sh):
  shell driver to reproduce the perf run.
- New oracle benchmarks under `tools/bench/`:
  - `tradcut_phase_2L_eps001_*` — first Tradition-with-cut parity
    test (D3 acceptance gate).
  - `banana_4loop_eps001_*` — Phase 3.F first L=4 oracle (D7 gate).
  - `ewbox_1loop_eps001_*` — Phase 3.G multi-invariant oracle.
  - `cutbanana_4L_eps001_*` — Phase 3.H 5-particle Cutkosky oracle.
  - `bn3mix_eps001_*` — Phase 3.I 3-loop mixed-mass oracle.
  - `cutbubble_1L_eps{2,10000}_*` — Phase 3.J ε-extremes oracles.
- New ending scheme: `EndingScheme::Trivial`, auto-appended to the
  user's `ending_schemes` list as a final fallback (mirror of
  upstream `AMFlow.m:1034`).
- New `solve_integrals` single-eps fast path: when
  `numeric_values["eps"]` is supplied, skip the Laurent fit and
  return the integral evaluated at that eps (mirror of
  `AMFlow.m:1364-1374`).
- New per-system path direction: `AMFSystem::setup()` computes the
  η-touching loops' prescription consensus and overrides the global
  `run_direction` for that system's ODE solve (mirror of
  `AMFlow.m:981-991`).
- Cutkosky setup now validates that all phase-volume component
  masses are non-negative after `Numeric` substitution; raises
  otherwise (mirror of `AMFlow.m:1050`).
- `apply_blackbox_options` (`src/api/run_json.cpp`) now actively
  rejects the complex-numeric object form
  `{"re":..,"im":..}` in `amf_options.blackbox.numeric_values`
  with a clear "not implemented" error pointing to audit divergence
  D5.  Prior behaviour would have failed deeper in the parser with an
  unclear message; the explicit rejection makes the feature gap
  visible at the input boundary.

### Known limitations (deferred indefinitely)
- **D5 — Kira `ComplexMode` / imaginary-numeric pipeline is not
  implemented.**  Users who need to evaluate at numeric kinematics
  with non-zero imaginary parts cannot do so via this port; the
  upstream filters such values through `IBPRule` / `CompensateRule`
  but the C++ port treats `numeric_values` as a flat real-valued
  map.  All 12 oracle benchmarks use purely-real numerics, so this
  surface is untested.  **Reclassified out of scope** (see Unreleased
  §Changed) — see [`docs/AUDIT_MMA_PARITY.md`](docs/AUDIT_MMA_PARITY.md) §D5.
## [1.0.0] — 2026-05-08

First public release.  C++17 reimplementation of the
auxiliary-mass-flow algorithm originally introduced in the Mathematica
package AMFlow by Liu and Ma
(<https://gitlab.com/multiloop-pku/amflow>; Comput. Phys. Commun. 283
(2023) 108565).  All algorithmic credit belongs to the upstream
authors; this project contributes only the C++ implementation and a
numerical-parity test/benchmark harness.

### Added
- Domain-oriented architecture under the `amflow::` namespace:
  `numeric`, `algebra`, `ode`, `qft`, `ibp`, `pipeline`, `api`, `cli`.
- Public headers under `include/amflow/<domain>/` with strict ABI/PIMPL
  boundaries between domains.
- ODE solver — reimplementation of the upstream `DESolver.m` covering
  the workflows exercised by the oracle benchmarks; matches the
  reference at relative precision ~1e-30 across all 12 oracle cases.
- AMFlow + Kira pipeline implementing the core algorithms of upstream
  `AMFlow.m`: `amflow`, `black_box_amflow`, and `solve_integrals` modes.
- `amflow_cli` — JSON-driven driver around `amflow::api::run_json`,
  exposing all three top-level modes.
- 500+ GoogleTest cases, including 12 sampled-parity oracle benchmarks
  under `tools/bench/` (`eps = 1/1000`).
- Mathematica reference drivers under `tools/math_ref/` and
  `tools/bench/` (require a local clone of upstream AMFlow; see
  `reference/README.md`), plus committed JSON reference outputs for
  regression tests.
- CMake install / export rules; downstream projects can use
  `find_package(AMFlowCpp)` and link `AMFlowCpp::amflow`.

### Notes
- The upstream Mathematica AMFlow source is **not** vendored.  See
  `reference/README.md` for how to obtain it locally to regenerate
  benchmark reference data.

### Out of scope
- `SolveIntegralsGaugeLink`, HQET / SCET / Wilson-line workflows — see
  [`docs/AUDIT_MMA_PARITY.md`](docs/AUDIT_MMA_PARITY.md).

[1.1.0]: https://github.com/chang18/amflow-cpp/releases/tag/v1.1
[1.0.0]: https://github.com/chang18/amflow-cpp/releases/tag/v1.0
