# Contributing & Debugging Guide

> 中文版本: [CONTRIBUTING_zh.md](CONTRIBUTING_zh.md)

This document is for practical project work: build, test, benchmark,
reference generation, and parity debugging.

---

## Workflow checklist

When you change code:

1. Build your active tree:
   ```bash
   cmake --build build -j32
   ```
2. Run the full local suite:
   ```bash
   ctest --test-dir build --output-on-failure -j 4
   ```
3. If a test fails, read the assertion message first. Most numeric tests
   already print the compared values.
4. Update the live docs when project status changes:
   `README.md`, `docs/AUDIT_MMA_PARITY.md`, `CHANGELOG.md`.

When you add or fix functionality:

1. Put the code in the correct domain module under `src/<domain>/`
   (numeric / algebra / ode / qft / ibp / pipeline / api / cli).
2. Add or extend the matching test file under `tests/`.
3. If the change alters public surface or project status, update docs in
   the same turn.
4. If the change closes or exposes a parity issue, preserve the
   Mathematica side of the comparison.

## Current validation policy

- Mathematica AMFlow is the spec.
- The immediate parity screen is sampled single-point values at
  `eps = 1/1000`.
- Do not make broader fitted/Laurent coverage the first move when a
  sampled parity comparison is open.
- Preserve every Mathematica output used to establish parity.

## Running tests

Single test:

```bash
./build/tests/amflow_tests \
  --gtest_filter='AmfsystemTest.*'
```

Run a single domain group:

```bash
./build/tests/amflow_tests \
  --gtest_filter='*Numeric*:*Qft*'
```

List tests:

```bash
ctest --test-dir build -N
```

## Generating Mathematica reference data

The Wolfram driver lives at `tools/math_ref/run_amflow_kira.wl` and
reads its config from `AMF_REF_CFG`.

```bash
cd tools/math_ref
AMF_REF_CFG=$PWD/cfg_bubble_1L.wl math -script run_amflow_kira.wl
```

The committed direct reference JSON files live under
`tests/data/math_ref/`. Sampled parity benchmark references live under
`tools/bench/`.

## Sampled parity workflow

Use this workflow for current family bring-up:

1. Choose integral targets and fix `eps = 1/1000`.
2. Generate or preserve Mathematica sampled values.
3. Save summarized values in `tools/bench/*_mma_reference.json`.
4. Preserve raw Mathematica caches in `tools/bench/mma_refs/` whenever
   they contain useful intermediate values.
5. Run the C++ sampled path on exactly the same family, replacements,
   numeric values, and targets.
6. If C++ mismatches, compare the first wrong direct subintegral or
   child-system value before touching fit/Laurent code.
7. If C++ matches, add or extend a focused regression when runtime is
   reasonable.

## Common debugging strategies

### Result is NaN

Usually one of:

- `prec` was omitted or became zero before an Arb call;
- a normalization path received invalid precision or a singular
  transform;
- the path generator produced a degenerate step;
- a Kira reduction returned incomplete data and downstream code accepted
  it blindly.

First checks:

```bash
./build/tests/amflow_tests \
  --gtest_filter='OdeAmflowTest.*'
./build/tests/amflow_tests \
  --gtest_filter='AmfsystemTest.BubblePipeline_RunsToCompletion'
```

### Result is finite but wrong by a simple factor or sign

Usually one of:

- 1-based vs 0-based translation;
- reverse coefficient indexing around the sparse system;
- missing prefactor propagation in a recursive boundary chain;
- wrong numeric substitution into family replacements before setup.

Check `docs/INVARIANTS.md` first.

### All masters come back as zero or the boundary list is empty

Look for:

- missing pattern-zero boundary entries;
- wrong `deinf` construction;
- wrong `determine_boundary_order` input type;
- ending/global value stitching problems.

### Recursion does not shrink or setup loops forever

Check:

- `ending_q` parity with `AMFlow.m`;
- whether child families truly drop loops or propagators;
- whether a child setup is being built from the fully substituted family
  data instead of a partly symbolic one.

### `MpolyContext mismatch`

You mixed symbolic objects created in different contexts. Re-parse via
string round-trip into the target context. Do not compare contexts by
name; the runtime rule is shared ownership identity.

### Kira errors or hangs

Check:

- absolute `kira_path` and `fermat_path`;
- `FERMATPATH` directory contents;
- dumped Kira workdir via `AMFLOW_KIRA_DUMP=1`;
- generated `integralfamilies.yaml`, `jobs.yaml`, `kinematics.yaml`,
  `preferred`, and `target`.

## Debug flags

Runtime environment variables — the canonical list lives in
[`include/amflow/numeric/log.hpp`](../include/amflow/numeric/log.hpp);
add new trace categories there when introducing them. Currently:

| Env var | Effect |
|---|---|
| `AMFLOW_KIRA_DUMP=1` | Keep the Kira workdir after a run |
| `AMFLOW_JORDAN_DUMP=...` | Dump Jordan-form intermediate matrices to file |
| `AMFLOW_NO_SPARSE_CHOP=1` | Disable the sparse-row chop heuristic |
| `AMFLOW_SPARSE_CHOP_DIGITS=N` | Override sparse-chop digit count (default = `numeric::chop_pre()`) |
| `AMFLOW_DEBUG_BC=1` | Emit boundary-condition and `solve_one_eps` diagnostics |
| `AMFLOW_DEBUG_DBO=1` | Emit `determine_block_boundary_order` diagnostics |
| `AMFLOW_DEBUG_PRECISE=1` | Verbose precision tracing in normalize / Jordan paths |
| `AMFLOW_DEBUG_REGULAR=1` | Emit `regular_run` step-by-step trace |
| `AMFLOW_DEBUG_STAGES=1` | Emit per-stage system traces (Inf→Reg, RegRun, Zero) |
| `AMFLOW_DEBUG_CT_PERM=1` | Emit `ct_perm` Cutkosky permutation trace |
| `AMFLOW_DEBUG_SCHEME=1` | Emit `[scheme] ...` EndingScheme dispatch trace from `amf_system_setup_master` (Cutkosky / SingleMass firing) |
| `AMFLOW_TRACE_LAYER6=1` | Fine-grained `src/ode/inf.cpp` trace (boundary-Taylor / Frobenius diagnostics) |

## Adding tests

Use the existing pattern:

```cpp
class Layer5Test : public ::testing::Test {
protected:
    void SetUp() override { set_default_options(); }
    void TearDown() override { set_default_options(); }
};
```

Always reset options in tests. Many failures that look numerical are
just leaked global settings.

For numeric checks, prefer helpers already used in the suite rather than
string comparison. Relax tolerances only when the workflow genuinely
accumulates more numerical error.

## Things to avoid

- Modifying the upstream Mathematica AMFlow to match C++ behaviour
  (the upstream is the spec; we conform to it, not the other way
  around).
- Replacing parity work with weaker tests or looser tolerances.
- Mixing symbolic contexts without re-parsing.
- Excluding a ported Mathematica path from parity checks instead of
  debugging the mismatch.
