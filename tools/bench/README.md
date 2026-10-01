# Sampled Benchmark Assets

> 中文版本: [README_zh.md](README_zh.md)

This directory stores sampled `eps = 1/1000` parity inputs, Mathematica
reference values, and C++ benchmark configs.

Mathematica AMFlow is the reference implementation. Keep every
Mathematica value and raw cache produced for these benchmarks.

## Compare Existing Outputs

`compare_sampled.py` compares a C++ sampled output JSON against the
committed Mathematica reference JSON. It reads the requested/top/master
reference entries and reports real/imaginary absolute and relative
errors per integral.

```bash
tools/bench/compare_sampled.py \
  tools/bench/tt_2loop_box_eps001_mma_reference.json \
  /tmp/amflow_bench_tt_2loop_box_cpp_out.json

tools/bench/compare_sampled.py \
  tools/bench/doublebox_sv_eps001_mma_reference.json \
  /tmp/amflow_bench_doublebox_sv_eps001_cpp_out.json
```

Defaults:

- relative tolerance: `1e-25`
- absolute tolerance: `1e-40`

Use `--json-out <path>` to write a machine-readable comparison report.

## Run Long C++ Benchmarks

These commands run the sampled C++ benchmark path. They are long-running
benchmark commands, not default `ctest` entries.

The commands below assume `amflow_cli` is on `$PATH` (after
`sudo cmake --install build`).  If you haven't installed, substitute
`./build/src/cli/amflow_cli` for each `amflow_cli`.

```bash
amflow_cli \
  tools/bench/tt_2loop_box_eps001_black_box_amflow_cpp.json \
  /tmp/amflow_bench_tt_2loop_box_cpp_out.json

amflow_cli \
  tools/bench/doublebox_sv_eps001_black_box_amflow_cpp.json \
  /tmp/amflow_bench_doublebox_sv_eps001_cpp_out.json

amflow_cli \
  tools/bench/tt_higher_rank_eps001_black_box_amflow_cpp.json \
  /tmp/amflow_bench_tt_higher_rank_eps001_cpp_out.json

amflow_cli \
  tools/bench/xbox_2loop_eps001_black_box_amflow_cpp.json \
  /tmp/amflow_bench_xbox_2loop_eps001_cpp_out.json

amflow_cli \
  tools/bench/tt_cutkosky_probe_eps001_black_box_amflow_cpp.json \
  /tmp/amflow_bench_tt_cutkosky_probe_eps001_cpp_out.json

AMFLOW_DEBUG_SCHEME=1 amflow_cli \
  tools/bench/cutbubble_1L_eps001_black_box_amflow_cpp.json \
  /tmp/amflow_bench_cutbubble_1L_eps001_cpp_out.json

AMFLOW_DEBUG_SCHEME=1 amflow_cli \
  tools/bench/cutsunrise_2L_eps001_black_box_amflow_cpp.json \
  /tmp/amflow_bench_cutsunrise_2L_eps001_cpp_out.json

AMFLOW_DEBUG_SCHEME=1 amflow_cli \
  tools/bench/cutbanana_3L_eps001_black_box_amflow_cpp.json \
  /tmp/amflow_bench_cutbanana_3L_eps001_cpp_out.json
```

After a run, compare the output with `compare_sampled.py` before
updating benchmark status fields.

## Sampled Benchmark Inventory

The live inventory lives in the files themselves: each bench is a
triplet `<name>_{cpp.json, mma.wl, mma_reference.json}` (sub-master
oracles reuse parent caches and have only `cpp.json + mma_reference.json`).
Use `ls tools/bench/*_mma_reference.json` for the full list (228
entries as of 2026-05-17).

Wallclock numbers are intentionally not maintained here; see
[`docs/ROADMAP.md`](../../docs/ROADMAP.md) §"Performance benchmarking".

## Notes

- `tt_cutkosky_probe` only differs from `tt_2loop_box` in the
  `EndingScheme` priority list. With this kinematics (s = 30, t = -10/3,
  msq = 1) the boundary integrals are handled by `Tradition`; the
  presence of `Cutkosky` in the list is exercised at the protocol level
  but does not perturb numerics.
- `cutbubble_1L` is the *opposite* counterpart to `tt_cutkosky_probe`:
  the family already has `cut = {1, 1}` at the top level, so the top
  sector satisfies `phase_volume_q` (|var| = L+1, all cuts == 1) and
  `AMFSystemSetupMaster` chooses `Cutkosky` for the root system. Run
  with `AMFLOW_DEBUG_SCHEME=1` to confirm the firing line
  `[scheme] Cutkosky fired: family=cutbubble phase_loop_num=1 ...` is
  printed. The eps^0 leading coefficient is exactly `1/(8*Pi)`, the
  massless 2-particle phase-space volume.
- `cutsunrise_2L` extends the Cutkosky probe to `L = 2`: 5-prop basis
  (3 actual + 2 ISPs) with `cut = {1, 1, 1, 0, 0}`, target
  `j[cutsunrise, 1, 1, 1, 0, 0]`. `phase_volume_q` triggers because
  `|var| = 3 = L+1` with all live cuts == 1, so `AMFSystemSetupMaster`
  picks Cutkosky on the top family. The firing line reads
  `[scheme] Cutkosky fired: family=cutsunrise phase_loop_num=2 top_position={0,1,2}`.
  This bench is the strongest available check on the
  `(Pi^(2-eps)*(2*Pi)^(2*eps-4))^L * (-1)^(L+1)` Cutkosky prefactor at
  `L > 1` (the sign factor flips relative to the 1-loop case) and on
  `phase_loop_num` propagation through `AMFSystemSetupMaster`.
- `cutbanana_3L` pushes the Cutkosky probe to `L = 3`: 9-prop basis (4
  actual + 5 ISPs) with `cut = {1, 1, 1, 1, 0, 0, 0, 0, 0}`, target
  `j[cutbanana, 1, 1, 1, 1, 0, 0, 0, 0, 0]`. After cut clearing the
  eta-flow runs as a 3-loop banana with eta on `l1^2` (1 massive + 3
  massless legs); subsystem recursion goes `system 1 -> 2 -> 3 ->
  trivial`. Firing line:
  `[scheme] Cutkosky fired: family=cutbanana phase_loop_num=3 top_position={0,1,2,3}`.
  Together with `cutbubble_1L` (sign +1) and `cutsunrise_2L` (sign -1),
  this completes the L=1/2/3 sign-pattern check on `(-1)^(L+1)` (back
  to +1 at L=3) and exercises the deepest recursion in the bench
  inventory. **Requires `BlackBoxDot = 5`** (matching the equal-mass
  `banana_3loop`): with `BlackBoxDot = 1` the boundary Reduce step
  fails with `master count differs between Masters and Reduce calls`.
  Both the MMA `*_mma.wl` and C++ `*_cpp.json` configs set this floor
  explicitly — AMFlow.m's `Max[$BlackBoxDot, JDot/@...]` pre-compute
  cannot raise it on its own when the input integrals all have
  `JDot = 0` (every positive index is exactly 1).
- `doublebox2m` is the regression guard for D12 (see
  `docs/AUDIT_MMA_PARITY.md`).  Interleaved cross-loop mass placement
  gives the top-sector differential-equation matrix a complex-conjugate
  pole pair (`η² + η + 12 = 0`, `|η| ≈ 3.46`), tightening the
  Frobenius series convergence radius for the 4 top-sector masters.
  The committed C++ config sets `working_pre=200, x_order=400,
  extra_x_order=480` — roughly doubled vs the all-massless / block-mass
  defaults — to converge.  At the older defaults `(160, 200, 240)` C++
  produced a sign-flipped Re and `O(0.1)` spurious Im; the elevated
  parameters are the recommended baseline for any topology with
  complex-pole pairs near the NegIm contour.
- `banana_3loop` requires `BlackBoxDot = 5` for both the Masters and
  Reduce passes — at lower dot the Reduce step cannot cover the master
  set returned by Masters, and AMFlow's protocol-level sanity check
  refuses to proceed. The committed configs (`*_cpp.json` and
  `*_mma.wl`) set `dot = 5`. Pak symmetry detection in Kira
  (`search_symmetry_relations`) must also be enabled for this family;
  with symmetry off some sub-banana sectors do not collapse and the
  Reduce master count drifts again.

## Regenerating Mathematica caches (`mma_refs/`)

Raw Mathematica caches under `mma_refs/` are large (~80 MB) and **not
committed**. They can be regenerated from the `*_mma.wl` scripts in
this directory. For each bench triplet `<name>_{cpp.json, mma.wl,
mma_reference.json}`:

```bash
cd tools/bench
math -script <name>_mma.wl
```

Each `*_mma.wl` writes its raw cache directory under `mma_refs/<name>_mma_cache/`
and the summarised reference values into `<name>_mma_reference.json` (which
**is** committed).

If `tools/math_ref/cache/` is also missing (likewise gitignored), use
`tools/math_ref/run_amflow_kira.wl` driven by a `cfg_*.wl` file to
regenerate.
