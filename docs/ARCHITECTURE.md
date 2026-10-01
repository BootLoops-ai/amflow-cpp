# Architecture

> 中文版本: [ARCHITECTURE_zh.md](ARCHITECTURE_zh.md)

This document describes the domain-oriented architecture of the AMFlow
C++ port.

---

## Two ways to use the library

**Use mode 1 — direct ODE solve.** The user already has a rational ODE
matrix `A(η)` and an asymptotic boundary specification at η = ∞. The
goal is to evaluate `I(η = 0)`. This is the path implemented entirely
within the `ode` domain, mirroring `amflow[de, bc]` in `DESolver.m`:

```
    dI(η)/dη  =  A(η) · I(η)                 (rational matrix)
    I_i(η)   ~  Σ_k v_{i,k} · η^{μ_{i,k}}    near η = ∞
                                               (BoundarySpec)
```

**Use mode 2 — full AMFlow workflow.** The user supplies a *family*
(loop momenta, legs, replacement rules, propagators) and a list of
target master integrals; the library produces per-ε Laurent expansions
of every master. The `pipeline` domain wraps `ode` + `qft` + `ibp` to
make this happen automatically: Kira IBP, differential equation
construction, region decomposition, asymptotic boundary integrand
construction, recursive sub-system building.  See
[`AUDIT_MMA_PARITY.md`](AUDIT_MMA_PARITY.md) for the validated parity
surface.

---

## Domain DAG

```
                          ┌──────────────┐
                          │   numeric    │   options / log / acb wrap /
                          │              │   rational helpers
                          └──────┬───────┘
                                 │
                          ┌──────▼───────┐
                          │   algebra    │   mpoly, mpoly_matrix,
                          │              │   context migration
                          └─┬─────────┬──┘
                            │         │
                  ┌─────────▼─┐   ┌───▼────────┐
                  │    ode    │   │    qft     │   family_config / family_uf
                  │           │   │            │   topology / amfmode /
                  │           │   │            │   region / boundary / vacuum /
                  │           │   │            │   jintegral / dense_q
                  └─────┬─────┘   └─┬────────┬─┘
                        │           │        │
                        │   ┌───────▼────┐   │
                        │   │    ibp     │   │   libp_deriv, kira (yaml/run/
                        │   │            │   │   parse), reduce
                        │   └─────┬──────┘   │
                        │         │          │
                        │  ┌──────▼──────────▼┐
                        │  │     pipeline    │  factorize, AMFSystem,
                        │  │                 │  black_box_amflow{,_single},
                        │  │                 │  GenerateNumericalConfig,
                        │  │                 │  FitEps, SolveIntegrals
                        │  └────────┬────────┘
                        │           │
                        │     ┌─────▼─────┐
                        │     │    api    │  JSON dispatcher (run_json)
                        │     └─────┬─────┘
                        │           │
                        │     ┌─────▼─────┐
                        │     │    cli    │  amflow_cli executable
                        │     └───────────┘
                        │
                       (mode-1 entry: `ode::amflow(de, bcs)`
                        is reachable directly from api/cli too)
```

- Each domain is a sub-namespace of `amflow::` and a directory under
  `include/amflow/<domain>/` + `src/<domain>/`.
- `numeric → algebra` is the foundation; both the `ode` chain and the
  `qft → ibp → pipeline` chain build on top.
- A unit-test target `amflow_tests` covers all domains; the cli
  integration test runs `amflow_cli` as a subprocess.

---

## Mode-1 data flow (`ode` domain only)

```
   ┌─────────────────┐   inf_to_regular   ┌─────────────────┐
   │  asymptotic BC  │ ─────────────────▶ │  numerical BC   │
   │   at η = ∞      │   (ode::inf)       │   at η = x_0    │
   └─────────────────┘                    └────────┬────────┘
                                                   │ regular_run
                                                   │  (ode::regular)
                                                   ▼
                                          ┌─────────────────┐
                                          │  numerical BC   │
                                          │   at η = x_n    │
                                          └────────┬────────┘
                                                   │ solve_asy_exp
                                                   │  (ode::zero, ode::normalize,
                                                   │   ode::jordan, ode::qqbar)
                                                   ▼
                                          ┌─────────────────┐
                                          │  asymptotic     │
                                          │  expansion at   │
                                          │  η = 0          │
                                          └────────┬────────┘
                                                   │ pick_zero_solution
                                                   ▼
                                              I(η = 0)
```

The contour `[x_0, x_1, ..., x_n]` is a chain of regular complex points
joining "near infinity" to "near zero" along a path that avoids all
poles, generated by `ode::path`'s `run_eta`.

The public entry is `amflow::ode::amflow(de, bcs, prec)`.

## Mode-2 data flow (`pipeline` domain wrapping `qft` + `ibp` + `ode`)

```
   ┌──────────────────────────────────────────────────────────┐
   │  user inputs                                             │
   │    * qft::FamilyConfig                                   │
   │    * preferred masters: List<qft::JIntegral>             │
   │    * numeric values for invariants                       │
   │    * pipeline::AMFSystemOptions  (modes, ending schemes) │
   └─────────────────────────┬────────────────────────────────┘
                             │  pipeline::amf_systems_setup
                             ▼
        ┌──────────────────────────────────────────────────┐
        │  pipeline::AMFSystem (root)                      │
        │   if ending → look up qft::vacuum / explicit BC  │
        │   else:                                          │
        │     a. qft::analyze_top_sector                   │
        │     b. qft::amf_position(modes) → η_c            │
        │     c. inject η into propagators                 │
        │     d. ibp::diffeq                               │
        │           ├── ibp::reduce → masters list         │
        │           ├── ibp::libp_deriv → ∂Master/∂η       │
        │           └── ibp::kira_run analytic reduction   │
        │           ⇒ dM/dη = A(η, ε) · M                  │
        │     e. qft::find_all_region                      │
        │     f. for each non-zero region:                 │
        │           qft::region_power      → μ(ε)          │
        │           qft::boundary_integrands → integrand   │
        │           qft::boundary_integrals → per-family   │
        │           recurse: build child AMFSystem on each │
        │                     boundary family              │
        └─────────────────────────┬────────────────────────┘
                                  │  pipeline::amf_systems_solution(epslist)
                                  ▼
        ┌──────────────────────────────────────────────────┐
        │  for each ε in epslist:                          │
        │   bottom-up traversal of the AMF tree:           │
        │    for ending: lookup qft::vacuum / explicit     │
        │    for non-ending:                               │
        │      build BC by combining child solutions       │
        │      via the recorded RegionBoundary tables      │
        │      → ode::amflow(de_at_eps, BC)                │
        │      → master values for this node               │
        └─────────────────────────┬────────────────────────┘
                                  ▼
                  per-ε master values for the root preferred[]
                  (Laurent expansion follows via pipeline::FitEps,
                   user-facing entry: pipeline::solve_integrals)
```

### How `pipeline` reuses `ode`

Each non-ending AMF node ultimately calls `ode::amflow(de_at_eps, BC)`
(the mode-1 entry). `pipeline` only **builds** the per-ε ODE matrix and
boundary spec, and **decides** when to stop recursing.

---

## What each domain owns

| Domain | Headers (selected) | Responsibilities |
|--------|---------------------|------------------|
| `numeric` | `options.hpp`, `acb_wrap.hpp`, `log.hpp`, `rational.hpp` | Working precision, RAII wrappers around acb/arb, rational-function & rational-matrix types, `get_poles` |
| `algebra` | `mpoly.hpp`, `mpoly_matrix.hpp`, `context_migration.hpp` | Multivariate ℤ-polynomials & ℚ rationals over named contexts; matrix det/adjugate/inverse; cross-context migration |
| `ode` | `blocks.hpp`, `sparse.hpp`, `asy.hpp`, `regular.hpp`, `inf.hpp`, `zero.hpp`, `path.hpp`, `amflow.hpp`, `qqbar.hpp`, `jordan.hpp`, `normalize.hpp`, `acb_rational.hpp` | Block analysis, sparse Gauss, asymptotic-expansion algebra, regular & singular point recurrences, `NormalizeMat`/Jordan/`Calcx00`, contour `RunEta`, top-level `amflow()` driver |
| `qft` | `family_config.hpp`, `family_uf.hpp`, `jintegral.hpp`, `topology.hpp`, `amfmode.hpp`, `region.hpp`, `boundary.hpp`, `vacuum.hpp`, `dense_q.hpp` | Family description, Symanzik U/F, sector algebra, η-injection strategies, region decomposition, boundary-integrand construction, vacuum closed-form table |
| `ibp` | `libp_deriv.hpp`, `kira.hpp`, `reduce.hpp` | IBP derivatives, Kira yaml writers / process driver / output parser, `reduce` and `diffeq` (Kira-only entry) |
| `pipeline` | `amfsystem.hpp`, `factorize.hpp`, `solve_integrals.hpp` | `AMFSystem` recursion tree, ending schemes (Tradition/Cutkosky/SingleMass), `factorize_family`, `BlackBoxAMFlow{,Single}`, `FitEps`, `SolveIntegrals` |
| `api` | `run_json.hpp` | JSON dispatcher (mode = `amflow` / `solve_integrals` / `black_box_amflow`); takes `nlohmann::json` directly |
| `cli` | (none) | Thin `main.cpp` that hands argv → `api::run_json` |

---

## Where to start reading code

| Goal | Entry points |
|---|---|
| Mode-1 pipeline (ode domain) | `src/ode/amflow.cpp::amflow()`, `src/ode/inf.cpp::calc_inf()`, `src/ode/zero.cpp::calcx00()` (read side-by-side with upstream `diffeq_solver/DESolver.m:856-917` from <https://gitlab.com/multiloop-pku/amflow>) |
| Mode-2 pipeline (pipeline domain) | `include/amflow/pipeline/amfsystem.hpp` header comment, then `src/pipeline/amfsystem.cpp::AMFSystem::setup()`/`solve()`, `src/pipeline/factorize.cpp::factorize_family()`, `src/ibp/reduce.cpp::diffeq()` |
| Data types | `numeric/options.hpp`, `numeric/rational.hpp`, `ode/asy.hpp`, `algebra/mpoly.hpp`, `qft/family_config.hpp`, `qft/jintegral.hpp` |
| Invoking from a host app | `src/cli/main.cpp` (one-screen JSON wrapper); `tests/test_*.cpp` (each domain has its own test file); `tests/test_cli_integration.cpp` (subprocess JSON example) |
