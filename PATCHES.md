> **Layout note.** This file is the BootLoops project's record of its changes
> to AMFlow.cpp. Path names below map onto this repository as follows:
> `amflow-cpp/` = the repository root; `wrappers/` = `bootloops-wrappers/`;
> `upgrades/kira` = the BootLoops project's separate Kira fork (the sibling
> repository `kira`, carrying the SQLite transaction fix also
> referenced from `BUILD.md`), not part of this repository. This fork is part of
> the BootLoops release; the main repository, `bootloops` (MIT), is published
> beside it by the same organization, does not include this fork, and points here
> from its `upgrades/ENGINES.md`. The `LICENSE` and `NOTICE` at
> this repository's root are upstream AMFlow.cpp's, unchanged. Questions or bugs
> concerning the BootLoops modifications: open an issue on this repository (or
> the main BootLoops repository); upstream's issue tracker and maintainer
> address apply to unmodified AMFlow.cpp only.

# amflow-cpp fork + AMFlow-port wrapper suite

## What this is

`amflow-cpp/` is our fork of **AMFlow.cpp** (v1.1.0 + `main`@`006a275`;
<https://github.com/chang18/amflow-cpp>, MIT; Zenodo
doi:10.5281/zenodo.20087172), a C++17 reimplementation of the
auxiliary-mass-flow method written and maintained by the AMFlow.cpp
contributors (maintainer: GitHub \@chang18). AMFlow.cpp re-implements, and is
validated case by case against, the Mathematica package **AMFlow** by Xiao Liu
and Yan-Qing Ma (<https://gitlab.com/multiloop-pku/amflow>, MIT). All
algorithmic credit for the method belongs to the AMFlow authors, as upstream's
`NOTICE` says; all credit for the C++ engine, its architecture and its parity
harness belongs to the AMFlow.cpp contributors. The upstream `LICENSE` and
`NOTICE` files ship intact inside `amflow-cpp/`; our modifications are released
under the same MIT license.

`wrappers/` is our own validation harness for the engine (original code, no
upstream lineage): closed-form (mpmath) and Bessel-moment anchors, cross-checked
externally against pySecDec, plus the runner scripts we used to audit the
fork's numerics.

The engine serves as the IBP/numeric oracle for the BootLoops program —
1-loop through 4-loop families, including 15-propagator 3-loop topologies at
100+ requested digits. That is well outside the operating point upstream
targets (parity with the original AMFlow package at its default precision), and
each change below addresses a need that arose only at that scale. None of them
reflects a defect in upstream's design; several are candidates to offer back
(see the end of this file).

## Scientific credit and thanks

Nothing in this fork is new physics. The methods it implements or extends are:

- **Auxiliary mass flow.** X. Liu, Y.-Q. Ma and C.-Y. Wang, "A systematic and
  efficient method to compute multi-loop master integrals", Phys. Lett. B 779
  (2018) 353 [arXiv:1711.09572]; X. Liu and Y.-Q. Ma, "Multiloop corrections
  for collider processes using auxiliary mass flow", Phys. Rev. D 105 (2022)
  L051503 [arXiv:2107.01864]; Z.-F. Liu and Y.-Q. Ma, "Determining Feynman
  integrals with only input from linear algebra", Phys. Rev. Lett. 129 (2022)
  222001 [arXiv:2201.11637]; X. Liu, Y.-Q. Ma, W. Tao and P. Zhang,
  "Calculation of Feynman loop integration and phase-space integration via
  auxiliary mass flow", Chin. Phys. C 45 (2021) 013115 [arXiv:2009.07987]. The
  package: X. Liu and Y.-Q. Ma, "AMFlow: a Mathematica package for Feynman
  integrals computation via auxiliary mass flow", Comput. Phys. Commun. 283
  (2023) 108565 [arXiv:2201.11669]; and now R.-J. Huang, X. Liu and Y.-Q. Ma,
  "AMFlow 2.0: significant algorithmic and software improvements for Feynman
  integral evaluation" [arXiv:2607.08477], which ships the AMFlow authors' own
  high-performance C++ differential-equation solver. Users who want the
  reference implementation should use AMFlow itself.
- **Linear (eikonal/worldline) propagators.** The auxiliary-mass-flow
  treatment of integrals with linear propagators was developed by Z.-F. Liu
  and Y.-Q. Ma, "Automatic computation of Feynman integrals containing linear
  propagators via auxiliary mass flow", Phys. Rev. D 105 (2022) 074003
  [arXiv:2201.11636] (an auxiliary quadratic term for each linear propagator,
  removed again through the differential equations), and is available in the
  Mathematica AMFlow; that paper is the originating and reference treatment of
  this integral class and should be cited whenever the additions below are
  used. AMFlow.cpp lists the class as out of scope; the product-of-tadpoles
  vacuum, the half-η guard and the experimental eikonal region engine added
  here (§3) handle a subset of eikonal/worldline families inside this engine
  and claim nothing beyond that paper. The η→∞ region analysis is the method
  of expansion by regions, M. Beneke and V. A. Smirnov, Nucl. Phys. B 522
  (1998) 321 [hep-ph/9711391], as used for the η→∞ boundary conditions in the
  auxiliary-mass-flow papers above. The closed-form eikonal vacua are the
  standard HQET integrals: one-loop, see A. G. Grozin, *Heavy Quark Effective
  Theory*, Springer Tracts Mod. Phys. 201 (Springer, 2004); two-loop, D. J.
  Broadhurst and A. G. Grozin, "Two-loop renormalization of the effective
  field theory of a static quark", Phys. Lett. B 267 (1991) 105
  [hep-ph/9908362]; two-velocity (cusp) integrals, G. P. Korchemsky and A. V.
  Radyushkin, "Renormalization of the Wilson loops beyond the leading order",
  Nucl. Phys. B 283 (1987) 342.
- **Maximal-cut ending (Mercedes, experimental).** Baikov representation: P.
  A. Baikov, Phys. Lett. B 385 (1996) 404 [hep-ph/9603267] and Nucl. Instrum.
  Meth. A 389 (1997) 347 [hep-ph/9611449]; cuts in Baikov variables: H.
  Frellesvig and C. G. Papadopoulos, JHEP 04 (2017) 083 [arXiv:1701.07356];
  maximal cuts as boundary data for banana-type differential equations: A.
  Primo and L. Tancredi, Nucl. Phys. B 921 (2017) 316 [arXiv:1704.05465]. The
  Gram-determinant monomials were generated with SymPy (A. Meurer et al.,
  PeerJ Comput. Sci. 3 (2017) e103).
- **Vacuum boundary values.** The single-mass vacuum table extended in §3 is
  X. Liu and Y.-Q. Ma's (Mathematica AMFlow; see Phys. Rev. D 99 (2019) 071501
  [arXiv:1801.10523] and the AMF papers above); the all-loop banana closed
  form is the iterated one-loop bubble/tadpole result found in any textbook
  (e.g. V. A. Smirnov, *Analytic Tools for Feynman Integrals*, Springer Tracts
  Mod. Phys. 250 (Springer, 2012)).

**Engines we drive, and thanks.** Every reduction in this engine is performed
by **Kira** — created by P. Maierhöfer, J. Usovitsch and P. Uwer, "Kira — a
Feynman integral reduction program", Comput. Phys. Commun. 230 (2018) 99
[arXiv:1705.05610]; J. Klappert, F. Lange, P. Maierhöfer and J. Usovitsch,
"Integral reduction with Kira 2.0 and finite field methods", Comput. Phys.
Commun. 266 (2021) 108024 [arXiv:2008.06494]; F. Lange, J. Usovitsch and Z. Wu,
"Kira 3: integral reduction with efficient seeding and optimized equation
selection", Comput. Phys. Commun. 322 (2026) 109999 [arXiv:2505.20197] — using
**FireFly** for finite-field reconstruction — J. Klappert and F. Lange, Comput.
Phys. Commun. 247 (2020) 106951 [arXiv:1904.00009]; J. Klappert, S. Y. Klein
and F. Lange, Comput. Phys. Commun. 264 (2021) 107968 [arXiv:2004.01463] — and
**Fermat**, the computer algebra system by Robert H. Lewis (Fordham University,
<http://home.bway.net/lewis/>), for polynomial algebra. The robustness items in
§2 are about running these programs far outside their usual envelope; that
they run there at all is to their authors' credit, and we thank the Kira,
FireFly and Fermat authors for making them available. Our Kira-side patch
lives in the separate Kira fork (the sibling repository `kira`).
Arbitrary-precision arithmetic throughout is FLINT/Arb (W. Hart, D. Harvey, S.
Pancratz, F. Johansson, A. Ahlbäck and the FLINT team, <https://flintlib.org>;
Arb, now part of FLINT 3: F. Johansson, "Arb: efficient arbitrary-precision
midpoint-radius interval arithmetic", IEEE Trans. Comput. 66 (2017) 1281
[arXiv:1611.02831]); JSON I/O is nlohmann/json (Niels Lohmann, MIT); the tests
use GoogleTest (Google, BSD-3-Clause); jemalloc is by Jason Evans and
contributors. Please cite Kira, FireFly and the AMFlow papers above in any
publication that uses numbers produced with this fork.

## What changed vs upstream, and why

The changes over upstream v1.1.0, grouped by theme:

### 1. Precision correctness

The headline fix is the **exponent-separation solver fix**: the eta→0
matching solve classifies solution behaviors by exponent, and has to tell
`eta^(n + a*eps)` apart from `eta^n`. Upstream uses fixed chop/rationalization
tolerances, which are correct at the default precision it targets and match
the original AMFlow package there. At the much finer eps grids we request,
`a*eps` can fall below a fixed tolerance, the two behaviors are then classified
together, and the Laurent coefficients come out wrong without an error being
raised. The fix scales `chop_pre` (and the rationalization thresholds)
with the eps-grid scale so the separation margin tracks the request. In the
same set: an underflow-safe zero-imaginary test (supports
`rationalize_pre > 308`, where `10^-p` underflows a double), never chopping
exact eps samples, and integer/exponent resonance-tolerance tests in the
`calcx00` boundary-matching path (a resonant boundary at eps=0 previously
produced a spurious "interval too wide" failure — now either solved correctly or
reported honestly as a structural divergence).

### 2. Robustness of the embedded Kira/FireFly reduce

- **SIG6/SIG11 → FireFly auto-staging with a dot-cap.** On large eta-deformed
  top-sector reduces the Fermat back-substitution path aborts (SIGABRT — an
  internal assertion at ~35-45 GB peak in our runs, or std::bad_alloc under a
  memory cap; SIGSEGV on some families); the fork catches the abort and re-stages the reduce through
  FireFly finite-field reconstruction *with a capped dot-degree mandatory
  list*. Measured on the failing families: FireFly alone still stalls; the
  dot-cap plus FireFly together clears them.
- **CPU-progress stall watchdog.** A monolithic-eta FireFly reduce can
  stall indefinitely at zero CPU (nothing dies, so a blocking wait hangs
  indefinitely). The launcher now watches aggregate CPU time of the child
  process tree; only a continuous zero-progress window (default 180 s,
  `AMFLOW_KIRA_STALL_SECS`) declares a stall, kills the tree, and re-stages
  through the auto-staging path. Runtime is never the trigger.
- **Real child-error surfacing.** A dead Kira child used to report a bare
  "KILLED BY SIGNAL n"; the launcher now surfaces the log tail (e.g.
  `std::bad_alloc` under a memory cap vs an internal Kira assert), which is
  the difference between "raise the cap" and "report the fault".
- **`AMFLOW_KIRA_FFSAVE=1` — resumable FireFly reduces.** Upstream relaunches
  start each reduce from a clean directory — the safe default — so an
  interrupted multi-hour FireFly reconstruction restarts from the first probe.
  With the knob on, matching inputs preserve the
  dir and FireFly resumes from its saved states; mismatches go fresh; an
  unresumable save (killed pre-shift) is detected and cleared so a bad save
  cannot wedge the job. Default off, byte-identical to stock.
- **jemalloc auto-preload** into the Kira child (see `amflow-cpp/BUILD.md`),
  plus a Kira-side SQLite transaction fix carried in our Kira fork
  (the sibling repository `kira`).

### 3. Boundary and vacuum extensions

- **Massless-vacuum boundary auto-complete**: auto-completes the ISP basis
  and factorizes the eta→inf vacua of fully-massless multi-loop graphs, with
  the single-mass "banana" closed form generalized from the original five
  hard-coded entries to every loop order (value-checked against the original
  table in the test suite).
- **(L,L) product-of-tadpoles vacuum.** For families with propagators linear
  in the loop momenta (worldline/eikonal propagators; the class whose
  auxiliary-mass-flow treatment is due to Z.-F. Liu and Y.-Q. Ma,
  arXiv:2201.11636; see the credit section above): the eta→inf hard
  region collapses linear propagators to constants, leaving an L-loop L-line
  product-of-tadpoles master, `(-Gamma(-1+eps))^L`, that Kira does not
  factorize away; it is now recognized directly so the recursion terminates.
- **`AMFMode::Bubble`** (2-propagator sub-loop eta-exclusion) as an
  additional eta-placement mode.
- **Vacuum hard-zero fix**: caps the resonance chop precision so a
  legitimate exact-zero vacuum boundary is not misclassified
  (regression-tested).
- **Eikonal region engine** (`AMFLOW_EIKONAL_REGIONS`, experimental): region
  finding and boundary construction for eikonal/worldline kinematics;
  implemented for the one- and two-velocity classes up to three loops listed
  in `amfsystem.cpp`; other configurations stop with an explicit error
  (closed forms: the standard HQET one- and two-loop integrals,
  Broadhurst–Grozin 1991, and cusp integrals, Korchemsky–Radyushkin 1987;
  references in the credit section above).
- **Mercedes maximal-cut ending** (`AMFLOW_MERCEDES_CUT=<family>`, experimental,
  research use): substitutes the Baikov maximal cut of the equal-mass K4
  vacuum for the top-sector ending masters of one named 15-index family; the
  overall normalization relative to the AMFlow measure is supplied by the
  user via `AMFLOW_MERCEDES_CONV_RE/IM` and is not validated in this tree
  (Baikov 1996; Frellesvig–Papadopoulos 2017; Primo–Tancredi 2017).

### 4. Scale and operations

- **DE checkpoint/resume** for the per-system eta-flow solves (crash
  recovery on multi-hour runs). The checkpoint key covers every numeric
  option (`rationalize_pre`, `chop_pre`, the expansion orders, ...) and the
  values of the solve-affecting environment overrides it lists
  (`AMFLOW_DIFFEQ_NUMD`, `AMFLOW_SPARSE_CHOP_DIGITS`,
  `AMFLOW_FORCE_ENDING_DEPTH`, ...; see `kSolveEnvKnobs` in
  `src/pipeline/system_checkpoint.cpp`), so a re-run with any of them
  changed re-solves instead of replaying the stored solution.
- **eps-grid parallelism** and **mode:diffeq** with numeric-d substitution
  (`AMFLOW_DIFFEQ_NUMD`) inside the diffeq target reduce. With the
  substitution active the per-eps solve uses the exact `(4 - P/Q)/2` as its
  eps, matching the DE tables, which carry `d = P/Q` exactly; the sample is
  not rationalized in that mode. Every other rationalization (region
  exponents, path points, poles, the boundary-matching matrix) recovers the
  simplest rational within `10^-rationalize_pre` rather than truncating to
  decimal digits, so an exponent `a + b*eps` is exact and leaves no residual
  for the sparse solver to pivot on.
- **Node-tree parallelism** (`amf_options.node_parallel` /
  `AMFLOW_NODE_PARALLEL`): a non-ending node's sibling sub-system solves
  are mutually independent and are dispatched to a thread pool, for wide
  shallow AMF trees / short eps grids where the eps-grid lever alone
  leaves cores idle (the `tree_stats` diagnostic reports the headroom).
  Default off; bit-identical to the serial path (regression-tested,
  including concurrent ending-system Kira solves).
- **IBP cache** with load-adaptive scheduling (measured on a 14-job
  cold-start farm where redundant identical Kira runs cost ~2x wall).
- **Boundary-table JSON dump** (`AMFLOW_BOUNDARY_DUMP`) for offline
  inspection of the eta=inf boundary tables.
- Config hygiene: unrecognized JSON keys are never silently ignored: a stderr
  warning names where the key would be read (fatal under
  `AMFLOW_STRICT_CONFIG=1`); per-family Kira
  seeding knobs surfaced in the JSON schema (see `amflow-cpp/docs/JSON_SCHEMA.md` and `amflow-cpp/BUILD.md`).

### 5. Repository housekeeping

- `.github/workflows/docker.yml` runs only on `workflow_dispatch` and never
  logs in to a registry or pushes an image (`push: false`, no `packages: write`):
  the image would bundle third-party engines (Kira, FireFly, GiNaC/CLN, FLINT)
  under their own licenses, so nothing is published automatically from this
  fork. Upstream's repository still carries its own version of the workflow.
- `.gitignore` ignores `bootloops-wrappers/results/` (the harness's default
  output directory).
- `docs/FAQ.md` / `docs/FAQ_zh.md`: the parity line reads 228 of 228 (the
  three 4L `sunset_bubble` cases were closed upstream in `4a292ab`/`d7b093d`;
  upstream's FAQ still said 225).
- Editorial pass over our own comments and messages (plain wording for run
  interruptions, relaunches and observed engine failures; no internal run
  identifiers).
- Upstream's Mathematica benchmark reference dump
  `tools/bench/moller_scalar_qed/mma_refs/` (~32,600 files, ~0.6 GB) is not
  included in this snapshot; it is reference output only, not needed to build
  or test, and remains available in the upstream repository.

## Verification

The tree builds clean and passes the full test suite — **579/579 tests**
(upstream's plus ours, including the Mathematica-reference oracle
comparisons).

## The wrapper suite (`wrappers/`)

Our audit harness, usable as an acceptance gate for any build of the engine:

- `run_tier1.sh` — re-runs upstream's 228 oracle bench cases under
  `amflow-cpp/tools/bench/` and compares against the committed Mathematica
  references; on `006a275` in our environment (2026-06-11) 224/228 pass —
  three `pentabox_2L` cases abort in the eta→0 residue step and `wzbox_2L`
  under-resolves at the committed truncation orders; upstream reports 228/228
  (README status table; `docs/AUDIT_MMA_PARITY.md`); the four differences are
  specific to our environment and under investigation.
- `gen_tier23_configs.py` + `run_tier23.sh` — independent closed-form
  anchors: massless bubble (Gamma closed form), massless box (2F1/polylog
  closed form), massive bubble incl. above-threshold imaginary part, and the
  equal-mass sunrise (Bessel-moment anchor), at 30/100/300 goal digits.
- `ref_closed_forms.py`, `ref_sunrise_bessel.py` — the anchor generators
  (mpmath, engine-independent).
- `compare_orders.py` — order-by-order comparison of two result JSONs.
- `run_tier5.sh` — cost-vs-digits and cost-vs-eps-order scaling fits.

Point them at a build via `AMFLOW_REPO` (defaults to `../amflow-cpp`).

## Candidates to offer upstream

- The precision fixes, headed by the exponent-separation solver fix.
- SIG6/SIG11 → FireFly auto-staging + the stall watchdog.
- `AMFLOW_KIRA_FFSAVE` resumable reduces.

## License

MIT throughout. Upstream amflow-cpp is Copyright (c) 2026 AMFlow.cpp
contributors and is distributed under the MIT License; its `LICENSE` and
`NOTICE` files ship unmodified inside `amflow-cpp/`, and its per-file SPDX
headers are kept. The auxiliary-mass-flow method it implements is due to Xiao
Liu and Yan-Qing Ma (AMFlow, MIT); the original Mathematica package is not
bundled. The fork modifications described above and everything under
`wrappers/` are Copyright (c) 2026 Anthropic, PBC and are released under
the same MIT License; they were created by Matthew D. Schwartz, with the code
written by Claude (Anthropic) under his supervision. This is
not an officially supported Anthropic product; it is maintained by Matthew D. Schwartz
(https://www.bootloops.ai). See also `LICENSE` and `NOTICE` at
this repository's root (upstream's, unchanged); the main BootLoops repository
carries its own `THIRD_PARTY.md` listing this fork as an external engine.
