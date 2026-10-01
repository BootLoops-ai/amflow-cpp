> **BootLoops fork.** This repository is the BootLoops project's patched fork of
> [AMFlow.cpp](https://github.com/chang18/amflow-cpp) (MIT; Zenodo
> [10.5281/zenodo.20087172](https://doi.org/10.5281/zenodo.20087172)), the C++17
> reimplementation of the auxiliary-mass-flow method written and maintained by
> the AMFlow.cpp contributors (maintainer: GitHub
> [\@chang18](https://github.com/chang18)): a snapshot of upstream `main` at
> commit `006a275` (2026-05-22; v1.1.0 plus the [Unreleased] CHANGELOG entries)
> with the BootLoops project's modifications applied; upstream's history and
> releases live at the link above. We are grateful to the AMFlow.cpp maintainer:
> the 228-case Mathematica-oracle benchmark and the GoogleTest suite that ship
> with upstream are what made every change in this fork checkable. The
> auxiliary-mass-flow method is due to Xiao Liu, Yan-Qing Ma and their
> collaborators Chen-Yu Wang, Zhi-Feng Liu, Wei Tao, Peng Zhang and Rui-Jun
> Huang, and the original Mathematica package
> [AMFlow](https://gitlab.com/multiloop-pku/amflow) (MIT) is by Xiao Liu and
> Yan-Qing Ma (AMFlow 2.0 with Rui-Jun Huang); please cite their papers, listed
> with the other scientific credits in
> [PATCHES.md](PATCHES.md#scientific-credit-and-thanks), together with upstream's
> `CITATION.cff`, whenever results obtained with this fork are published. The IBP
> reductions this engine launches are done by Kira (created by Philipp
> Maierhöfer, Johann Usovitsch and Peter Uwer; developed today by Fabian Lange,
> Johann Usovitsch and Zihao Wu, with Jonas Klappert on Kira 2.0) with FireFly
> (Jonas Klappert, Sven Yannick Klein, Fabian Lange) and Fermat (Robert H.
> Lewis); see the same section for references. What changed relative to that
> base and why: [PATCHES.md](PATCHES.md); build and run notes for the additions:
> [BUILD.md](BUILD.md); our validation harness for the engine:
> [bootloops-wrappers/](bootloops-wrappers/). Upstream license and notices
> (`LICENSE`, `NOTICE`, `CITATION.cff`) are unchanged; our modifications are
> released under the same MIT license, Copyright (c) 2026 Anthropic, PBC. Created by
> Matthew D. Schwartz; code for the modifications written by Claude (Anthropic) under
> his supervision. This is not an officially supported Anthropic product; it is
> maintained by Matthew D. Schwartz (https://www.bootloops.ai).
>
> The "Project provenance" section, `AGENTS.md`, `CONTRIBUTING.md`, `NOTICE`,
> `CITATION.cff`, the CI badge and the maintainer contact below are upstream
> AMFlow.cpp's and describe upstream's own development (upstream states its C++
> source was developed primarily by AI coding agents under its maintainer's
> direction); they do not describe the BootLoops modifications. What this fork
> vouches for is its own delta ([PATCHES.md](PATCHES.md)) and the acceptance
> evidence it re-ran: the GoogleTest suite, the 228-case Mathematica-oracle
> parity harness and the engine-independent closed-form/Bessel anchors in
> [bootloops-wrappers/](bootloops-wrappers/). Questions or bugs concerning the
> BootLoops modifications: open an issue on this repository (or on the main
> BootLoops repository, `bootloops`, published beside this one by the same
> organization); the CI badge, issue links and maintainer address below
> are upstream's and apply to unmodified AMFlow.cpp only.

## Maintenance, reporting and security

**Maintenance.** This repository is maintained by Matthew D. Schwartz, not by
Anthropic. It is not an officially supported Anthropic product, and Anthropic
does not provide support, updates or fixes for it.

**Reporting issues.** Please report bugs and security problems through this
repository's GitHub issues.

**Security considerations.** Treat input files from others as code. These are
research tools meant to be run locally on inputs you trust. Many of them
evaluate the contents of their input files (JSON, YAML, `.m`, `.ms`, `.jl`,
pickle and similar), so a file received from someone else can run arbitrary
commands on your machine. Only run files you wrote yourself or got from a
source you trust, or run them in a sandbox or container. The integrity checks
and certificates in this repository guard against accidents. They are not a
security boundary.

# AMFlow.cpp

[![CI](https://github.com/chang18/amflow-cpp/actions/workflows/ci.yml/badge.svg)](https://github.com/chang18/amflow-cpp/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Version](https://img.shields.io/badge/version-1.1.0-blue.svg)](CHANGELOG.md)
[![C++17](https://img.shields.io/badge/C%2B%2B-17-blue.svg)](https://en.cppreference.com/w/cpp/17)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20087172.svg)](https://doi.org/10.5281/zenodo.20087172)

> 中文版本: [README_zh.md](README_zh.md)

A C++17 reimplementation of the auxiliary-mass-flow algorithm for
multi-loop Feynman integrals.

## Scientific attribution

This project re-implements the algorithms presented in the Mathematica
package **AMFlow** by Xiao Liu and Yan-Qing Ma:

- Upstream Mathematica package: <https://gitlab.com/multiloop-pku/amflow>
- Algorithm paper: X. Liu and Y.-Q. Ma, *AMFlow: A Mathematica package
  for Feynman integrals computation via auxiliary mass flow*, Comput.
  Phys. Commun. **283** (2023) 108565,
  [doi:10.1016/j.cpc.2022.108565](https://doi.org/10.1016/j.cpc.2022.108565).

This repository contributes **only** the C++17 implementation and a
numerical-parity test/benchmark harness against the Mathematica
reference; **all algorithmic credit belongs to the AMFlow authors**.
We claim no algorithmic novelty.

If you use this software in published research, please cite **both** the
upstream paper above (algorithm) and this repository (implementation).
The repository is archived on Zenodo with DOI
[10.5281/zenodo.20087172](https://doi.org/10.5281/zenodo.20087172)
(concept DOI; always resolves to the latest release).  Machine-readable
metadata is in [`CITATION.cff`](CITATION.cff); see [`NOTICE`](NOTICE)
for full attribution.

> **Note** — the upstream Mathematica source is **not** vendored.
> Documentation in this repo refers to symbols and line numbers in
> upstream files (`AMFlow.m`, `diffeq_solver/DESolver.m`,
> `ibp_interface/Kira/interface.m`); to follow those references and
> to regenerate Mathematica reference data, clone the upstream into
> `reference/amflow-master/` (gitignored).  See
> [`reference/README.md`](reference/README.md).

## Project provenance

This C++17 implementation was developed primarily by **AI coding
agents** under the direction of the maintainer, with every behavioural
change gated by numerical-parity verification against the upstream
Mathematica reference.  The 228 oracle benchmarks under `tools/bench/`
(1L=52, 2L=53, 3L=64, 4L=59) and the 549-case GoogleTest suite are
the contract that the AI-led implementation has to honour for any
change to land.

Maintainer / contact: **3250800970@qq.com** (please open a GitHub
issue first when possible; use email for inquiries that don't fit a
public issue).

## Layout

Domain-oriented architecture:

```
amflow::numeric -> amflow::algebra -> {amflow::ode,
                                       amflow::qft -> {amflow::ibp,
                                                       amflow::pipeline -> amflow::api -> amflow::cli}}
```

Public headers under `include/amflow/<domain>/`; implementation under
`src/<domain>/`.

## Working rules

- If a workflow succeeds in upstream `AMFlow.m` and the corresponding
  branch is already ported here, a C++ mismatch is a parity bug.
- Do not invent extra correctness criteria beyond Mathematica parity.
- Sampled single-point parity at `eps = 1/1000` is the primary
  numerical contract; see `tools/bench/`.
- All Mathematica-computed benchmark values are kept committed as
  regression standards (raw caches are regenerable; see
  `tools/bench/README.md`).

## Status

| Area | State |
|---|---|
| Core ODE solver (port of upstream `DESolver.m`) | Implemented; 228/228 oracle cases match MMA reference at `rel ~1e-30` (rel ~1e-10 on the pentabox 2L corner where the integral's intrinsic cancellation horizon dominates) |
| AMFlow + Kira pipeline (upstream `AMFlow.m` core algorithms) | Implemented for the covered workflows |
| Top-level entries (`amflow`, `black_box_amflow`, `solve_integrals`) | Implemented; exposed via `amflow_cli` JSON modes |
| Line-level MMA parity audit | Closed; every public symbol has at least one oracle exercising it. See [`docs/AUDIT_MMA_PARITY.md`](docs/AUDIT_MMA_PARITY.md) |
| Wall-clock vs MMA | Benchmarking deferred until performed on a dedicated machine; see [`docs/ROADMAP.md`](docs/ROADMAP.md) §"Performance benchmarking" |
| `SolveIntegralsGaugeLink`, HQET / SCET / Wilson lines | Out of scope (intentional) |
| Complex-valued numeric kinematics | Out of scope (not supported); JSON dispatcher rejects `{"re":..,"im":..}` form. See [`docs/FAQ.md`](docs/FAQ.md) "What's *not* implemented?" |

See [`docs/AUDIT_MMA_PARITY.md`](docs/AUDIT_MMA_PARITY.md) for the
line-level upstream-parity audit (full divergence inventory) and
[`tools/bench/`](tools/bench/) for the oracle bench triplets.

## Quick start

On Ubuntu, the packaged dependencies needed for the default build are:

```bash
sudo apt-get update
sudo apt-get install -y build-essential cmake pkg-config libflint-dev \
  libgtest-dev nlohmann-json3-dev
```

```bash
cmake -S . -B build -DAMFLOW_BUILD_DRIVER=ON
cmake --build build -j32
ctest --test-dir build --output-on-failure -j 4

# Install to /usr/local/{bin,lib,include} (default CMake prefix).
# Puts `amflow_cli` on PATH and exposes the library + headers for
# downstream find_package().  Needs sudo for system-wide; pass
# --prefix "$HOME/.local" instead for a user-local install.
sudo cmake --install build

# Raw ODE example (no IBP).  Without the install step, run the binary
# from build/ directly: ./build/src/cli/amflow_cli examples/power_law.json
amflow_cli examples/power_law.json

# Full AMFlow examples (require Kira + Fermat at runtime)
amflow_cli examples/box1_black_box_amflow.json
amflow_cli examples/bubble_solve_integrals.json
amflow_cli examples/box1_solve_integrals.json
```

To use as a CMake dependency once installed:

```cmake
find_package(AMFlowCpp 1.1 REQUIRED)
target_link_libraries(my_target PRIVATE AMFlowCpp::amflow)
```

## Mathematica reference generation

The repository includes Wolfram drivers for regenerating committed
reference data under `tests/data/math_ref/` and sampled benchmark data
under `tools/bench/`.  These drivers expect a local clone of the
upstream Mathematica AMFlow at `reference/amflow-master/`
(gitignored — see [`reference/README.md`](reference/README.md) for the
two acceptable layouts).

```bash
# Either clone upstream directly into reference/amflow-master/ ...
git clone https://gitlab.com/multiloop-pku/amflow.git reference/amflow-master

# ... or symlink to an existing clone:
# ln -s /path/to/your/clone reference/amflow-master

cd tools/math_ref
AMF_REF_CFG=$PWD/cfg_bubble_1L.wl math -script run_amflow_kira.wl
```

For current parity benchmarks at `eps = 1/1000` see
[`tools/bench/`](tools/bench/) and
[`docs/AUDIT_MMA_PARITY.md`](docs/AUDIT_MMA_PARITY.md).

## Documentation map

**For users** — read in this order:

1. [`docs/USER_GUIDE.md`](docs/USER_GUIDE.md) — step-by-step walk-through of one no-IBP example and one full Feynman-integral example.
2. [`docs/JSON_SCHEMA.md`](docs/JSON_SCHEMA.md) — every field of every input/output, for all three CLI modes.
3. [`docs/FAQ.md`](docs/FAQ.md) — common build / runtime / parity questions.

**For contributors** — additionally:

| Document | What it covers |
|---|---|
| [`CONTRIBUTING.md`](CONTRIBUTING.md) | TL;DR contribution rules + maintainer contact |
| [`docs/CONTRIBUTING.md`](docs/CONTRIBUTING.md) | Detailed dev workflow: tests, benchmarks, references, debugging, env-var trace flags |
| [`AGENTS.md`](AGENTS.md) | AI-coding-agent onboarding (the AI-led development contract) |
| [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md) | Domain DAG and end-to-end data flow |
| [`docs/INVARIANTS.md`](docs/INVARIANTS.md) | Project-wide rules that must stay true (precision, RAII, ODE-boundary traps) |
| [`docs/REFERENCE_MAP.md`](docs/REFERENCE_MAP.md) | Upstream Mathematica symbol → C++ symbol mapping |
| [`docs/AUDIT_MMA_PARITY.md`](docs/AUDIT_MMA_PARITY.md) | Live parity status against upstream MMA AMFlow (per-row divergence catalog + closure log) |
| [`docs/ROADMAP.md`](docs/ROADMAP.md) | Active development plan (forward-looking — oracle expansion) |
| [`reference/README.md`](reference/README.md) | How to clone the upstream MMA AMFlow locally for reference data regeneration |
| [`tools/bench/README.md`](tools/bench/README.md) | Sampled-parity benchmark conventions and regeneration |

## Dependencies

| Component | Why |
|---|---|
| C++17 compiler | Core build |
| CMake >= 3.16 | Build system |
| FLINT / Arb | Rational, polynomial, and arbitrary-precision arithmetic |
| GoogleTest | Test suite |
| nlohmann/json | Driver JSON I/O |
| Kira + Fermat | IBP reduction backend (runtime; tests gracefully skip when absent) |
| Mathematica / Wolfram Engine + upstream AMFlow clone | Optional; only needed to regenerate reference data |

## Contributing

See [`CONTRIBUTING.md`](CONTRIBUTING.md) for the GitHub-facing entry
point and [`docs/CONTRIBUTING.md`](docs/CONTRIBUTING.md) for the full
project workflow.  Behavioural changes must be backed by upstream
AMFlow parity evidence and focused C++ regression tests.

For questions, bug reports, or collaboration inquiries, open a GitHub
issue at <https://github.com/chang18/amflow-cpp/issues>, or contact
the maintainer at <3250800970@qq.com>.

## License

MIT — see [`LICENSE`](LICENSE) and [`NOTICE`](NOTICE).  The upstream
Mathematica AMFlow is also MIT-licensed; this project does not bundle
any upstream source.
