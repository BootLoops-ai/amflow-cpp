> **Layout note.** This file is the BootLoops project's record of its changes
> to AMFlow.cpp. Path names below map onto this repository as follows:
> `amflow-cpp/` = the repository root; `wrappers/` = `bootloops-wrappers/`;
> `upgrades/kira` = the BootLoops project's separate Kira fork (the sibling
> repository `kira`, carrying the SQLite transaction fix also
> referenced from `BUILD.md`), not part of this repository. This fork is part of
> the BootLoops release; the main repository, `bootloops` (MIT), is published
> beside it by the same organization, does not include this fork, and points here
> from its `upgrades/ENGINES.md`. The `LICENSE` at this repository's root is
> upstream AMFlow.cpp's, unchanged. `NOTICE` is upstream's (its first 35 lines,
> unchanged) followed by sections that this fork appended and marked as its
> own. This repository does not carry upstream's commit history: it starts from
> a snapshot of upstream `main` at commit `006a275`, and upstream's history and
> tags are in upstream's repository. Questions or bugs
> concerning the BootLoops modifications: open an issue on this repository (or
> the main BootLoops repository); upstream's issue tracker and maintainer
> address apply to unmodified AMFlow.cpp only.

# amflow-cpp fork + AMFlow-port wrapper suite

## What this is

`amflow-cpp/` is our fork of **AMFlow.cpp** (v1.1.0 + `main`@`006a275`;
<https://github.com/chang18/amflow-cpp>, MIT; Zenodo
doi:10.5281/zenodo.20087172), the C++17 port of the Mathematica package
**AMFlow** by Xiao Liu and Yan-Qing Ma
(<https://gitlab.com/multiloop-pku/amflow>, MIT), written and maintained by
the AMFlow.cpp contributors (maintainer: GitHub \@chang18). Upstream calls
itself a reimplementation; its `docs/REFERENCE_MAP.md` maps each Mathematica
symbol to its C++ counterpart, with line numbers that track AMFlow commit
`efda1db` (version 1.2), "the snapshot taken when this port was performed",
and its source comments cite `AMFlow.m`, `DESolver.m` and `Kira/interface.m`
by line. AMFlow.cpp is validated case by case against AMFlow. All
algorithmic credit for the method belongs to the AMFlow authors, as upstream's
`NOTICE` says; all credit for the C++ engine, its architecture and its parity
harness belongs to the AMFlow.cpp contributors. The derivative routines of
AMFlow's Kira interface, which AMFlow.cpp ports in `src/ibp/libp_deriv.cpp`,
are T. Peraro's LiteIBP (see "Scientific credit and thanks" and the License
section). The upstream `LICENSE` ships unchanged inside `amflow-cpp/`, and
upstream's `NOTICE` text ships unchanged as the first 35 lines of `NOTICE`,
with our sections appended below it; our modifications are released under the
same MIT license, except as stated under "LiteIBP" in `NOTICE`.

The BootLoops additions were written with the Mathematica source of AMFlow 2.0
(R.-J. Huang, X. Liu and Y.-Q. Ma; tag 2.0, commit `2600551`) open. Where an
addition follows that source, the entry below says so.

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
- **Derivatives of integrals with respect to invariants (LiteIBP).** The
  derivative routines of AMFlow's Kira interface are from `LiteIBP.m` by
  Tiziano Peraro, part of FiniteFlow MathTools
  (<https://github.com/peraro/finiteflow-mathtools>, `packages/LiteIBP.m`;
  T. Peraro, "FiniteFlow: multivariate functional reconstruction using finite
  fields and dataflow graphs" [arXiv:1905.08019]). AMFlow says so in its
  source, `ibp_interface/Kira/interface.m`: "compute derivatives (code from
  LiteIBP.m provided by T.Peraro)" (line 280 at commit `efda1db`, line 400 at
  tag 2.0). Upstream AMFlow.cpp ported that block as
  `src/ibp/libp_deriv.cpp` and `include/amflow/ibp/libp_deriv.hpp`. The
  momentum-derivative chain rule that this fork added to those files (§6)
  takes the place of LiteIBP's `LIBPDerivivative` (`LiteIBP.m` line 565 at
  commit `f77e2b5`; `Kira/interface.m` line 383 at `efda1db`, line 503 at tag
  2.0). LiteIBP builds on LiteRed by R. N. Lee.
- **Boundary matching at η→0.** The matching-equation count and the
  rationalization precision of the matching system restored in §1 are those
  of `diffeq_solver/DESolver.m` in AMFlow 2.0 (Huang, Liu and Ma), read from
  its source.
- **Reduction to Fuchsian form.** The kernel-shear rescue in §3 uses the step
  that its source comment names the Lee/Barkatou right-kernel shear: the
  reduction to Fuchsian form of R. N. Lee [arXiv:1411.0911] and the Moser
  reduction of Barkatou and Pfluegel.

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

The main changes over upstream v1.1.0, grouped by theme. The themed entries
describe what was added and why; they are not a line-by-line record. The
complete list of files this fork added or changed is in §7, and the exact
difference is what `git diff` reports against upstream commit `006a275`.

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

Two of the changes in the `calcx00` path follow the Mathematica source of
AMFlow 2.0, `diffeq_solver/DESolver.m` (Huang, Liu and Ma, MIT), which was read
for them (three blocks of `src/ode/zero.cpp`, 47 non-blank lines of code and
comment):

- **Matching-equation count.** The number of match equations kept is
  `compensate = n_var - (number of raw pairing equations)`, the rule of
  `DESolver.m` lines 1011-1014 at tag 2.0 (`compensate = Length[allvar] -
  Length[pairtoequ]`), in place of the rank-filtered count of the C++ port.
- **Rationalization precision of the matching system.** The boundary-matching
  matrix is rationalized at `rationalize_pre` digits, as `AMFRationalize` does
  in `DESolver.m` (line 135 at tag 2.0), in place of `working_pre`.

The rank-selected fallback for resonant blocks that follows them in
`zero.cpp` has no counterpart in the Mathematica package. A comment in
`src/ode/inf.cpp` quotes `DESolver.m` as well (the `getboundary` null check).

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
- **Kira 3 seed truncation, `kira3_truncate_sp`**
  (`amf_options.blackbox.kira3_truncate_sp` / `AMFLOW_KIRA3_TRUNCATE_SP`):
  when on, the Reduce and Masters jobs written for Kira carry a `truncate_sp`
  block with `l` set to the number of loops. This is an option of Kira 3
  (Lange, Usovitsch and Wu, arXiv:2505.20197); the fork only writes it into
  the job file. Default off, set per family (see `amflow-cpp/BUILD.md`). The
  companion knob `kira3_sectorwise` is also default off.
- **Symbolic-kinematics reduction, `symbolic_ibp`**
  (`amf_options.blackbox.symbolic_ibp` / `AMFLOW_SYMBOLIC_IBP`): the
  kinematics file is written without numeric substitution, Kira reduces over
  the symbolic invariants once, and the per-point numbers are substituted
  into the parsed result. With the IBP cache (§4) every kinematic point then
  shares one reduction. Default off.

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
- **Half-integer coset lattice** (`AMFLOW_ODE_HALF_COSETS=1`, with
  `AMFLOW_BC_PARTIAL_ZEROFILL=1` for partly supplied boundary data;
  `src/ode/inf.cpp`): the boundary content of eikonal/worldline families puts
  the members of one coupled η→∞ block on two cosets of the integer power
  lattice, with exponents differing by half-integers. The block equation is
  rewritten in `t = sqrt(x)`, which puts both cosets on one integer lattice,
  and the existing boundary-order routine is applied at doubled depth.
  Default off. The Mathematica AMFlow treats linear propagators by a
  different route (an auxiliary quadratic term, arXiv:2201.11636).
- **Kernel-shear rescue for the reduction to Fuchsian form**
  (`AMFLOW_TOFUCHSIAN_KERNEL_RESCUE=1`; `src/ode/normalize.cpp`,
  `try_kernel_shear_rescue`): when the projector search of
  `to_fuchsian_local` stops with "irreducible Poincare rank", a projector
  built from the right kernel of the leading matrix (then from the left
  kernel) is tried, and accepted only if it lowers the Poincare rank, or
  keeps it and lowers the rank of the leading matrix. Method: R. N. Lee
  [arXiv:1411.0911]; Barkatou and Pfluegel (Moser reduction). Default off;
  behavior is unchanged unless the knob is set and the projector search
  fails.
- **Fixed-point and cut-boundary switches.** Five names that a downstream
  driver, PMflow (part of the BootLoops toolkit), uses to close the boundary
  of families where the eta→inf recursion does not end on its own. PMflow
  sets `AMFLOW_ETAC_OVERRIDE` and `AMFLOW_DUMP_EPS_GRID`, reads the file
  written under `AMFLOW_FIXEDPOINT_PROBE`, and recognizes the two guard
  errors in the run's output. The three variables are unset by default;
  the two guards are on by default.
  - `AMFLOW_FIXEDPOINT_PROBE=<file.json>`: at the top system, skip the
    boundary recursion, write the matrix `A(eps)` of `M(0) = A(eps)·v` to
    the file, and exit without solving the targets.
  - `AMFLOW_ETAC_OVERRIDE=<list>`: replace the automatic eta placement of
    the top system with a comma-separated per-propagator list.
  - `AMFLOW_DUMP_EPS_GRID`: print the eps grid, master values, reduction
    coefficients and target values of the black-box solve to stderr
    (output only).
  - `GRAVITYFLOW_FIXEDPOINT`: the name of the error raised when a boundary
    family reproduces its parent's propagators with no smaller top sector;
    `AMFLOW_ALLOW_FIXEDPOINT_RECURSION=1` suppresses it.
  - `GRAVITYFLOW_CUTREGION`: the name of the error raised when a family with
    a cut quadratic propagator loses every eta→inf region, which would give
    false zeros; `AMFLOW_ALLOW_CUT_FALSE_ZERO=1` suppresses it.

  The two `GRAVITYFLOW_` names are error tags, not environment variables.
  `docs/CONTRIBUTING.md` lists the values each switch takes and whether it
  enters the checkpoint key.

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

### 6. Derivatives with respect to kinematic invariants

- **Momentum-derivative chain rule** (`src/ibp/libp_deriv.cpp`:
  `build_kinematic_chain`, `decompose_chain`, `accumulate_chain_for_prop` and
  the helpers `legleg_key` and `invert_mfrac_matrix_local`, 253 non-blank
  lines; declared in `include/amflow/ibp/libp_deriv.hpp`; tested in
  `tests/test_ibp_libp_deriv.cpp`). Upstream's `libp_deriv` differentiates
  with respect to symbols that appear directly in the propagators (η, masses).
  Its header said that the momentum-derivative path of the Mathematica code,
  `LIBPDerivivative`, was intentionally not implemented. This fork implements
  that path, for invariants defined through replacement rules (for example
  `s` defined by `(p1+p2)^2 -> s`). The routine whose place it takes, and the
  block of AMFlow's Kira interface it plugs into, are T. Peraro's LiteIBP as
  carried in AMFlow (see the credit section above). The fork's routine
  follows `LIBPDerivivative` in what it computes and where it plugs in. It
  obtains the chain-rule coefficients differently: in closed form from the
  inverse of the external Gram matrix, `C = 1/2 G^-1 J`, where LiteIBP sets up
  a linear system for the derivative operators and solves it. For a variable
  that appears in no replacement rule the added term is zero and the earlier
  code path is taken unchanged. The License section states the license status
  of LiteIBP.

### 7. Files added and files changed by this fork

Compared with upstream `main` at commit `006a275` by content hash. Upstream's
Mathematica benchmark reference dump (§5) is left out of the comparison.

**Added by this fork (21 files).** These files are the BootLoops
project's. The C++ files among them carry the line "Copyright (c) 2026
Anthropic, PBC. Added by the BootLoops fork of AMFlow.cpp" under the SPDX
line; the scripts under `bootloops-wrappers/` carry the Anthropic copyright
header.

```
BUILD.md
PATCHES.md
bootloops-wrappers/README.md
bootloops-wrappers/compare_orders.py
bootloops-wrappers/gen_tier23_configs.py
bootloops-wrappers/ref_closed_forms.py
bootloops-wrappers/ref_sunrise_bessel.py
bootloops-wrappers/run_tier1.sh
bootloops-wrappers/run_tier23.sh
bootloops-wrappers/run_tier5.sh
examples/dbox_massless_solve_integrals.json
include/amflow/ibp/ibp_cache.hpp
include/amflow/ibp/load_adaptive.hpp
include/amflow/pipeline/system_checkpoint.hpp
src/ibp/ibp_cache.cpp
src/ibp/numd_subst.cpp
src/pipeline/system_checkpoint.cpp
tests/test_ibp_cache.cpp
tests/test_ibp_kira_stall.cpp
tests/test_pipeline_numd_frame.cpp
tests/test_pipeline_system_checkpoint.cpp
```

**Upstream files changed by this fork (73 files).** These are
upstream AMFlow.cpp's files with BootLoops modifications. Upstream's SPDX
line is kept and no other header was added, so the file itself does not show
which lines are ours; `git diff 006a275` against upstream's repository does.
`src/ibp/libp_deriv.cpp` and `include/amflow/ibp/libp_deriv.hpp` are in this
list and carry a provenance note (LiteIBP, §6).

```
.github/workflows/docker.yml
.gitignore
CHANGELOG.md
CHANGELOG_zh.md
CMakeLists.txt
NOTICE
README.md
README_zh.md
cmake/AMFlowCppConfig.cmake.in
docker/Dockerfile
docs/CONTRIBUTING.md
docs/CONTRIBUTING_zh.md
docs/DOCKER.md
docs/DOCKER_zh.md
docs/FAQ.md
docs/FAQ_zh.md
docs/JSON_SCHEMA.md
docs/JSON_SCHEMA_zh.md
include/amflow/ibp/kira.hpp
include/amflow/ibp/libp_deriv.hpp
include/amflow/ibp/reduce.hpp
include/amflow/numeric/acb_wrap.hpp
include/amflow/numeric/options.hpp
include/amflow/ode/amflow.hpp
include/amflow/ode/asy.hpp
include/amflow/ode/inf.hpp
include/amflow/ode/regular.hpp
include/amflow/ode/sparse.hpp
include/amflow/ode/zero.hpp
include/amflow/pipeline/amfsystem.hpp
include/amflow/qft/amfmode.hpp
include/amflow/qft/vacuum.hpp
src/algebra/mpoly.cpp
src/algebra/numeric_subst.cpp
src/api/run_json.cpp
src/ibp/CMakeLists.txt
src/ibp/kira_run.cpp
src/ibp/kira_yaml.cpp
src/ibp/libp_deriv.cpp
src/ibp/reduce.cpp
src/numeric/acb_wrap.cpp
src/numeric/matrix.cpp
src/numeric/options.cpp
src/ode/amflow.cpp
src/ode/asy.cpp
src/ode/inf.cpp
src/ode/normalize.cpp
src/ode/path.cpp
src/ode/qqbar.cpp
src/ode/regular.cpp
src/ode/sparse.cpp
src/ode/zero.cpp
src/pipeline/CMakeLists.txt
src/pipeline/amfsystem.cpp
src/pipeline/solve_integrals.cpp
src/qft/amfmode.cpp
src/qft/boundary.cpp
src/qft/complete.cpp
src/qft/findregion.cpp
src/qft/region.cpp
src/qft/topology.cpp
src/qft/vacuum.cpp
tests/CMakeLists.txt
tests/test_amflow_amfsystem.cpp
tests/test_amflow_amfsystem_mma_ref.cpp
tests/test_amflow_solve_integrals.cpp
tests/test_api_run_json.cpp
tests/test_ibp_kira.cpp
tests/test_ibp_libp_deriv.cpp
tests/test_numeric_acb_wrap.cpp
tests/test_ode_sparse.cpp
tests/test_ode_zero.cpp
tests/test_qft_vacuum.cpp
```

**Unchanged.** The remaining 884 files are byte-identical to upstream at
`006a275`.

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

Upstream amflow-cpp is Copyright (c) 2026 AMFlow.cpp
contributors and is distributed under the MIT License; its `LICENSE` ships
unmodified inside `amflow-cpp/`, its `NOTICE` text ships unmodified as the
first 35 lines of `NOTICE`, and its per-file SPDX
headers are kept. The auxiliary-mass-flow method it implements is due to Xiao
Liu and Yan-Qing Ma; AMFlow.cpp is a port of their Mathematica package AMFlow
(MIT; AMFlow 2.0 with Rui-Jun Huang), whose license text is reproduced in
`NOTICE`. The Mathematica package itself is not
bundled. The fork modifications described above and everything under
`wrappers/` are Copyright (c) 2026 Anthropic, PBC and are released under
the same MIT License, except as stated in the next paragraph; they were
created by Matthew D. Schwartz, with the code
written by Claude (Anthropic) under his supervision. This is
not an officially supported Anthropic product; it is maintained by Matthew D. Schwartz
(https://www.bootloops.ai). See also `LICENSE` (upstream's, unchanged) and
`NOTICE` (upstream's text followed by the sections this fork appended) at
this repository's root; the main BootLoops repository
carries its own `THIRD_PARTY.md` listing this fork as an external engine.

**LiteIBP.** `src/ibp/libp_deriv.cpp` and `include/amflow/ibp/libp_deriv.hpp`
derive from `LiteIBP.m` by Tiziano Peraro
(<https://github.com/peraro/finiteflow-mathtools>, `packages/LiteIBP.m`, read
at commit `f77e2b5`), through the block of AMFlow's
`ibp_interface/Kira/interface.m` that AMFlow marks "code from LiteIBP.m
provided by T.Peraro": upstream AMFlow.cpp ported that block, and this fork
added the momentum-derivative chain rule (§6) in the place of LiteIBP's
`LIBPDerivivative`. The finiteflow-mathtools repository carries no license
file, at that commit or at any of its 19 commits, and `LiteIBP.m` states no
license terms. The Anthropic copyright and the MIT grant on the fork
modifications do not cover the portion of these two files that derives from
LiteIBP.
