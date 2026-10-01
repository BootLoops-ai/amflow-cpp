# FAQ

> 中文版本: [FAQ_zh.md](FAQ_zh.md)

Common questions for users and packagers of AMFlow.cpp.  If your
question isn't here, please [open a GitHub issue](https://github.com/chang18/amflow-cpp/issues)
or contact the maintainer at <3250800970@qq.com>.

---

## Build and install

### Q. CMake fails with `Could not find module 'flint' in version >= 3.0`.

Your system FLINT is too old.  Symbols this project needs
(`fmpz_mpoly_q_set_str_pretty`, `fmpz_mpoly_q_evaluate_acb`, …) appeared
in FLINT 3.x.  Notably **Ubuntu 22.04 ships FLINT 2.8.4**, which won't
work.  Options:

- **Build FLINT 3.x from source** (~5–10 min):
  ```bash
  curl -fsSL https://github.com/flintlib/flint/releases/download/v3.4.0/flint-3.4.0.tar.gz \
       | tar -xz
  cd flint-3.4.0
  ./bootstrap.sh && ./configure --prefix=/opt/flint --disable-static
  make -j$(nproc) && sudo make install
  echo "/opt/flint/lib" | sudo tee /etc/ld.so.conf.d/flint.conf
  sudo ldconfig
  export PKG_CONFIG_PATH=/opt/flint/lib/pkgconfig:$PKG_CONFIG_PATH
  ```
  Then re-run `cmake -S . -B build` — pkg-config will see the new FLINT.
- **Use Ubuntu 24.04 / Fedora / Arch / Conda-forge**, which package
  FLINT 3.x.  See `.github/workflows/ci.yml` for a reproducible CI
  recipe.

### Q. Build succeeds but linker complains about missing `acb_*` / `arb_*` / `flint_*` symbols at run time.

Same root cause as above — multiple FLINT installations on the same
system, with the build picking one and the runtime loader picking
another.  Make sure `LD_LIBRARY_PATH` (or `/etc/ld.so.conf.d/*`)
exposes the same FLINT that pkg-config selected at configure time.

### Q. Can I build without Kira?

Yes.  Kira and Fermat are **runtime** dependencies of `solve_integrals`
and `black_box_amflow` only.  `cmake --build build -j` works fine
without them.  The 548-case test suite runs without Kira too — Kira-
dependent tests auto-skip when the `kira` binary isn't on `$PATH`.

The raw `amflow` mode (the ODE engine) doesn't need Kira at all — see
[`examples/power_law.json`](../examples/power_law.json) and
[`examples/constant.json`](../examples/constant.json) for runnable
examples.

### Q. How do I install the library and the `amflow_cli` binary system-wide?

After a successful build:

```bash
sudo cmake --install build                       # default prefix /usr/local
# or, for a user-local install (no sudo):
cmake --install build --prefix "$HOME/.local"
```

This installs four things under the chosen prefix:

- `bin/amflow_cli` — the JSON-driven CLI on `$PATH`.
- `lib/libamflow.a` — the static library.
- `include/amflow/...` — public headers per domain.
- `lib/cmake/AMFlowCpp/` — CMake config for downstream `find_package`.

For a custom prefix:

```bash
cmake --install build --prefix /some/prefix
```

### Q. How do I make CMake's `find_package(AMFlowCpp)` work for a downstream project?

In your downstream `CMakeLists.txt`:

```cmake
list(APPEND CMAKE_PREFIX_PATH /some/prefix)   # only if non-standard
find_package(AMFlowCpp 1.1 REQUIRED)
target_link_libraries(my_target PRIVATE AMFlowCpp::amflow)
```

The exported target name is `AMFlowCpp::amflow` (note the lowercase
`amflow`); the namespace is `AMFlowCpp::`.

---

## Running

### Q. `amflow_cli` exits with an error before producing JSON output.

The CLI wraps every error in a single line written to **stderr**.
Common ones:

| Message starts with | Meaning |
|---|---|
| `error: cannot open <path>` | Wrong input filename |
| `error: failed to parse JSON: ...` | Input file is not valid JSON |
| `unknown mode '...'` | `mode` field is not `"amflow"`, `"solve_integrals"`, or `"black_box_amflow"` |
| `missing required array '...'` | A required key is absent |
| `<key> must be an array / a string / ...` | Wrong JSON type |

Re-run with `2>&1` redirected to inspect the error.

### Q. The Kira step hangs forever.

Symptoms: `amflow_cli` is sitting at high CPU but never returns.

Most common causes:

- **Fermat path wrong.**  Set
  `amf_options.blackbox.fermat_executable` (or `FERMATPATH` env var)
  to the *directory* containing `fer64`, not the binary itself.
- **Kira can't reach Fermat.**  Run Kira directly on the dumped
  workdir to see the real error:
  ```bash
  cd /your/work_dir
  kira jobs.yaml
  ```
- **Workdir on a slow filesystem.**  Kira does heavy I/O.  Use a
  local SSD path, not a network mount.

To dump the Kira workdir for inspection (it's normally cleaned up):

```bash
AMFLOW_KIRA_DUMP=1 ./build/src/cli/amflow_cli your_input.json
```

### Q. Result is `NaN` or `inf`.

See [`docs/CONTRIBUTING.md`](CONTRIBUTING.md) §"Common debugging
strategies".  Usual culprits: working precision was clamped to zero, a
boundary specification was malformed, a Kira reduction returned
incomplete data.

### Q. Result is finite but disagrees with Mathematica AMFlow by a constant factor.

Likely culprits:

- **Wrong sign convention** in `propagators`.  This project uses the
  Minkowski-`+i0` convention and writes propagators as `momentum² −
  mass²`.  If your reference uses Euclidean signature the sign flips.
- **Numeric values not assigned** to all kinematic invariants.  Every
  invariant in `replacement` plus every bare mass in `propagators`
  must appear in `amf_options.blackbox.numeric_values`.
- **`D0` mismatch.**  The default space-time dimension is 4.  If you're
  in a non-integer dimension scheme, set
  `options.d0 = "<rational>"` to match.

### Q. How do I increase / decrease working precision?

`options.working_pre = N` sets working precision to N decimal digits
(default 100).  `options.chop_pre = M` sets the chop tolerance to
`10^-M` digits (default 20).  Both are in **decimal** at the API; the
implementation converts to bits internally.

For more on the precision policy see [`docs/INVARIANTS.md`](INVARIANTS.md) §1.

### Q. Can I compute in dimensions other than 4?

Yes.  The host space-time dimension is `D = D₀ - 2ε`, where `D₀` is
configurable via `options.d0` (rational `"p/q"` or integer; default
`"4"`):

```jsonc
"options": {
  "d0": "7/3"      // → D = 7/3 - 2ε
}
```

This is supported in `solve_integrals` and `black_box_amflow` (mirrors
upstream AMFlow.m's `D0` option).  The engine internally evaluates the
family at `ε + (4 - D₀)/2` and fits the Laurent expansion against the
user-facing grid, so `D₀` cancels cleanly at the API.

Validated at the oracle bench
[`tools/bench/box1_d0_7_3_solve_integrals_cpp.json`](../tools/bench/box1_d0_7_3_solve_integrals_cpp.json)
(D₀ = 7/3, matches Mathematica AMFlow to ~30 sig digits).
See [`docs/USER_GUIDE.md`](USER_GUIDE.md) §4.4 for a worked example
and [`AUDIT_MMA_PARITY.md`](AUDIT_MMA_PARITY.md) for the parity numbers.

---

## Numerical agreement with Mathematica AMFlow

### Q. The C++ result and the Mathematica result differ in their last few digits.  Bug?

Almost certainly **not** a bug.  AMFlow is a numerical algorithm; the
exact least-significant digits depend on internal floating-point
choices.  This project's parity contract is:

- **Sampled single-point comparison** at `eps = 1/1000` (or sometimes
  `eps = 1/100`) must match upstream to relative error
  `~10^-30` after subtracting working-precision noise.

If you see mismatches *bigger* than ~10⁻²⁵ at `working_pre = 100`,
that's worth investigating.  See [`AUDIT_MMA_PARITY.md`](AUDIT_MMA_PARITY.md)
and [`tools/bench/`](../tools/bench/) for the list of validated
families and observed relative errors.

### Q. How do I regenerate the committed Mathematica reference data?

You'll need:

- A working **Mathematica / Wolfram Engine** install.
- A **local clone of upstream AMFlow** at `reference/amflow-master/`
  (this project does not vendor it).  See
  [`reference/README.md`](../reference/README.md) for the two
  acceptable layouts (clone-into or symlink).

Then for benchmark caches:

```bash
cd tools/bench
math -script <name>_mma.wl     # writes <name>_mma_reference.json
```

Or for the larger reference set:

```bash
cd tools/math_ref
AMF_REF_CFG=$PWD/cfg_bubble_1L.wl math -script run_amflow_kira.wl
```

Full details: [`tools/bench/README.md`](../tools/bench/README.md).

---

## Project scope and contributing

### Q. What's *not* implemented?

**Out of scope (intentional, won't be added):**
- `SolveIntegralsGaugeLink`
- HQET / SCET / Wilson-line workflows
- IBP backends other than Kira (no FIRE / LiteRed / FiniteFlow / Blade)
- **Complex-valued numeric kinematics.**  Upstream MMA filters
  imaginary-part numerics through `IBPRule` / `CompensateRule` (apply
  real parts to Kira, substitute imaginary parts after).  This C++
  port treats `amf_options.blackbox.numeric_values` as a flat
  real-valued map; the C++ algebra layer is over Q (FLINT
  `fmpz_mpoly_q_t`) and cannot carry complex coefficients.  The JSON
  dispatcher rejects the object form `{"re":..,"im":..}` with a
  clear error.  Workaround: provide only purely-real numeric values.
  See [`docs/AUDIT_MMA_PARITY.md`](AUDIT_MMA_PARITY.md).

**See also:**
- [`docs/AUDIT_MMA_PARITY.md`](AUDIT_MMA_PARITY.md) — line-level
  upstream-parity audit (full divergence inventory).
- [`docs/REFERENCE_MAP.md`](REFERENCE_MAP.md) — Mathematica → C++
  symbol mapping (with explicit "not ported" entries).

### Q. Why "AI-developed" — should I trust the code?

This project's C++17 source was developed primarily by AI coding
agents under the maintainer's direction.  The trustworthiness contract
is **numerical parity**:

- 228 of 228 oracle benchmarks under `tools/bench/` match upstream
  Mathematica AMFlow at `rel ~ 10⁻³⁰` (rel ~ 10⁻¹⁰ on the pentabox 2L
  5-leg corner, where the integral's intrinsic cancellation horizon
  bounds the achievable precision at fixed `WorkingPre`); the three 4L
  2-leg `sunset_bubble` cases were closed upstream in 4a292ab/d7b093d.
- 548 GoogleTest cases gate the public surface.
- Every commit is gated by both, in CI on every push.

The tests (not the code) are the contract.  If the tests are right,
the code is right; if the code drifts, the tests catch it.  See
[`AGENTS.md`](../AGENTS.md) for how AI agents are instructed to
develop / verify changes.

### Q. I want to add a new family / port a missing upstream feature / fix a parity bug.

Read these in order:

1. [`AGENTS.md`](../AGENTS.md) — working rules and parity contract.
2. [`docs/CONTRIBUTING.md`](CONTRIBUTING.md) — practical workflow,
   debugging strategies, env vars.
3. [`docs/INVARIANTS.md`](INVARIANTS.md) — non-obvious traps in the
   precision / RAII / boundary logic.
4. [`docs/REFERENCE_MAP.md`](REFERENCE_MAP.md) — find your upstream
   function's C++ counterpart.

Then open a PR.  CI must stay green; behavioural changes must be
backed by upstream parity evidence.

---

## Citing

### Q. How should I cite this software?

Cite **both** the upstream paper (algorithm) and this repository
(implementation):

```bibtex
@software{amflow_cpp,
  title  = {{AMFlow.cpp}: A C++17 reimplementation of the AMFlow algorithm},
  author = {AMFlow.cpp contributors},
  year   = {2026},
  url    = {https://github.com/chang18/amflow-cpp},
  doi    = {10.5281/zenodo.20087172}
}

@article{liu_ma_amflow_2023,
  title   = {{AMFlow}: A {Mathematica} package for {Feynman} integrals computation via auxiliary mass flow},
  author  = {Liu, Xiao and Ma, Yan-Qing},
  journal = {Computer Physics Communications},
  volume  = {283},
  pages   = {108565},
  year    = {2023},
  doi     = {10.1016/j.cpc.2022.108565}
}
```

The Zenodo DOI `10.5281/zenodo.20087172` is the **concept DOI** —
always resolves to the latest release.  Per-version DOIs are listed
in [`CITATION.cff`](../CITATION.cff).
