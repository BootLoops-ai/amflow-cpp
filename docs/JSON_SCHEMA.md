# `amflow_cli` JSON Schema Reference

> 中文版本: [JSON_SCHEMA_zh.md](JSON_SCHEMA_zh.md)

`amflow_cli` consumes a single JSON document and emits a single JSON
document.  This page documents every field of every input/output, for
all three top-level modes.

The C++ entry point is `amflow::api::run_json(const json&)`
([`include/amflow/api/run_json.hpp`](../include/amflow/api/run_json.hpp));
this document and the implementation in
[`src/api/run_json.cpp`](../src/api/run_json.cpp) must agree.

> **Convention** — keys marked **required** must be present.  Keys
> marked **optional** have a documented default.  Where a field
> accepts multiple aliases, the first listed is preferred.

---

## Top-level shape

```jsonc
{
  "mode":        "amflow" | "solve_integrals" | "black_box_amflow",
  "options":     { ... },             // optional global options
  // ... mode-specific keys (see below)
}
```

If `mode` is omitted, `"amflow"` is assumed.  Output always contains
the keys `mode` (echoed), `options` (effective options after applying
`options`), and `result` (mode-specific shape).

---

## Mode `amflow` — raw ODE solve

Mirrors `AMFlow[de, bc]` in upstream `DESolver.m`.  Solves
`dI/dη = A(η)·I` with an asymptotic boundary specification at η = ∞
and returns the value at η = 0.  No IBP, no Kira, no family.

### Input

| Key          | Required | Type   | Meaning |
|---           |---       |---     |---      |
| `mode`       | optional | string | Must be `"amflow"` (or omitted) |
| `matrix`     | **required** | 2-D array | The rational ODE matrix `A(η)`, square `N × N`; see [Rational entries](#rational-entries) |
| `boundaries` | **required** | array of arrays | `boundaries[i]` is the BC list for row `i`; each entry is `{ "mu": <complex>, "value": <complex> }` |
| `options`    | optional | object | See [Global options](#global-options) |

`boundaries.length` must equal `matrix.length`.

### Output

```jsonc
{
  "mode":   "amflow",
  "result": [ <complex>, <complex>, ... ],   // length N
  "options": { /* effective options */ }
}
```

`result[i]` is the value `I_i(η = 0)`.

### Minimal example

```jsonc
// 1×1 system, dI/dη = (1/η)·I, boundary I ~ η^1 at infinity → I(0) = 1
{
  "options": { "x_order": 30, "extra_x_order": 10, "silent_mode": true },
  "matrix":     [[ { "num": ["1"], "den": ["0", "1"] } ]],
  "boundaries": [[ { "mu": { "re": "1", "im": "0" },
                    "value": { "re": "1", "im": "0" } } ]],
  "mode": "amflow"
}
```

(Run this with `amflow_cli examples/power_law.json` — or `./build/src/cli/amflow_cli examples/power_law.json` if you haven't installed the CLI.)

---

## Mode `solve_integrals` — Laurent expansion in ε

Mirrors `SolveIntegrals[jints, goal, order]` in upstream `AMFlow.m`.
The full pipeline: family → ε-grid → Kira IBP → AMFSystem recursion →
sampled values → Laurent fit.

### Input

| Key            | Required | Type   | Meaning |
|---             |---       |---     |---      |
| `mode`         | optional | string | Must be `"solve_integrals"` |
| `family`       | **required** | object | The integral family — see [Family object](#family-object) |
| `integrals` (alias: `jints`, `targets`) | **required** | array | Target integrals — see [Integral entry](#integral-entry) |
| `goal_digits` (alias: `goal`) | **required** | integer | Decimal digits of accuracy in the Laurent coefficients |
| `target_eps_power` | **preferred** | integer | The physical highest ε power K you want (e.g. `4` for "through ε⁴").  Internally sets `eps_order = K + 2L`.  Use this instead of `eps_order`. |
| `eps_order` (alias: `order`)  | required if `target_eps_power` absent | integer | (Legacy AMFlow convention.)  Number of ε orders fitted above the algorithmic leading −2L (L = loop count). Highest fitted ε power is `ε^(eps_order − 2L)`. Example: 1-loop with `eps_order=4` fits up to ε² (orders ε⁻² … ε²). To request "result up to ε^K", pass `eps_order = K + 2L`. Leading orders that vanish are trimmed from the output. **Pitfall**: a 2-loop with `eps_order=4` only reaches ε⁰; an `eps_order` below 2L is auto-corrected to 2L with a warning (set `AMFLOW_EPS_ORDER_NO_BUMP=1` to keep the requested under-order), and if every fitted order vanishes the manufactured `0` is a hard error (set `AMFLOW_ALLOW_TRUNCATED_ZERO=1` to demote). |
| `work_dir`     | optional | string  | Filesystem directory used for Kira intermediate files (default: a tmp dir) |
| `amf_options`  | optional | object | See [AMF options](#amf-options); also carries `numeric_values` |
| `options`      | optional | object | See [Global options](#global-options) |

`amf_options.blackbox.numeric_values` typically supplies the numeric
values of every kinematic invariant (e.g. `"s": "100"`, `"t": "-1"`).

### Output

```jsonc
{
  "mode":   "solve_integrals",
  "result": [
    {
      "integral":      { "family": "...", "indices": [...] },
      "leading_order": -2,                      // smallest ε power
      "coefficients":  [
        { "order": -2, "value": { "re": "...", "im": "..." } },
        { "order": -1, "value": { ... } },
        { "order":  0, "value": { ... } },
        ...
      ]
    },
    ...
  ],
  "options": { /* effective options */ }
}
```

Each entry is the Laurent expansion `Σ_k coefficients[k].value · ε^(leading_order + k)`.

---

## Mode `black_box_amflow` — sampled, no Laurent fit

Mirrors `BlackBoxAMFlow[jints, epslist, dir]`.  Returns each target
integral evaluated at every supplied ε sample, **without** fitting a
Laurent expansion.  Use this for parity comparison at single ε values
or when you need control over the ε grid yourself.

### Input

| Key | Required | Type | Meaning |
|---|---|---|---|
| `mode`         | optional | string | Must be `"black_box_amflow"` |
| `family`       | **required** | object | See [Family object](#family-object) |
| `integrals` (alias: `jints`, `targets`) | **required** | array | See [Integral entry](#integral-entry) |
| `eps_samples` (alias: `epslist`) | **required** | array of complex | The ε values to evaluate at; rationals like `"1/100"` are accepted |
| `work_dir`     | optional | string | Same as in `solve_integrals` |
| `amf_options`  | optional | object | See [AMF options](#amf-options) |
| `options`      | optional | object | See [Global options](#global-options) |

### Output

```jsonc
{
  "mode":   "black_box_amflow",
  "result": [
    {
      "integral": { "family": "...", "indices": [...] },
      "samples":  [
        { "eps":   { "re": "1/100", "im": "0" },
          "value": { "re": "...",   "im": "..." } },
        ...
      ]
    },
    ...
  ],
  "options": { /* effective options */ }
}
```

---

## Family object

Used by `solve_integrals` and `black_box_amflow`.

| Key            | Required | Type | Meaning |
|---             |---       |---   |---      |
| `name` (alias: `family`) | **required** | string | Family name (used as `j[<name>, …]`) |
| `loops`        | **required** | string[] | Loop momenta (e.g. `["l1", "l2"]`) |
| `legs`         | **required** | string[] | External momenta (e.g. `["p1", "p2", "p3", "p4"]`) |
| `propagators`  | **required** | string[] | Propagator denominators as algebraic strings (e.g. `"l^2 - msq"`) |
| `conservation` | optional | object | Map of momentum → expression for momentum conservation rewrites (e.g. `{"p4": "-p1 - p2 - p3"}`) |
| `replacement`  | optional | object | Map of scalar product → invariant (e.g. `{"(p1+p2)^2": "s"}`); used to rewrite kinematics in invariant variables |
| `cut`          | optional | int[] | One entry per propagator, `1` = on shell / cut, `0` = uncut |
| `prescription` | optional | int[] | One entry per propagator; `+1` / `-1` are the two ±i0 sign choices |

**Cut and prescription must each have length equal to `propagators`** if
present.

### Family example

```jsonc
"family": {
  "name": "box1",
  "loops":   ["l"],
  "legs":    ["p1", "p2", "p3", "p4"],
  "conservation": { "p4": "-p1 - p2 - p3" },
  "replacement":  {
    "p1^2": "0", "p2^2": "0", "p3^2": "0", "p4^2": "0",
    "(p1+p2)^2": "s",
    "(p1+p3)^2": "t"
  },
  "propagators": [
    "l^2",
    "(l + p1)^2",
    "(l + p1 + p2)^2",
    "(l + p1 + p2 + p4)^2"
  ]
}
```

---

## Integral entry

Each entry of `integrals` / `jints` / `targets` is one of:

```jsonc
// (a) bare indices array — uses the outer family
{ "indices": [1, 0, 1, 0] }

// (b) override family
{ "family": "box1", "indices": [1, 0, 1, 0] }
```

`indices` length must equal `propagators` length of the resolved
family.  Negative indices encode irreducible scalar products
(e.g. `[-2, 1, 1, 1]` = ISP rank 2 in slot 0).

---

## Global options

Top-level `"options"` object.  All fields are **optional**; defaults
are the upstream AMFlow defaults.

| Key              | Type           | Default | Mirrors upstream |
|---               |---             |---      |---               |
| `working_pre`    | int (digits)   | 100     | `WorkingPre`     |
| `chop_pre`       | int (digits)   | 20      | `ChopPre`        |
| `rationalize_pre`| int (digits)   | 100     | `RationalizePre` |
| `silent_mode`    | bool           | false   | `SilentMode`     |
| `d0`             | rational ("p/q") or int | `"4"`  | `D0` — host space-time dimension |
| `x_order`        | int            | 100     | `XOrder`         |
| `extra_x_order`  | int            | 20      | `ExtraXOrder`    |
| `learn_x_order`  | int            | -1      | `LearnXOrder`    |
| `test_x_order`   | int            | 5       | `TestXOrder`     |
| `run_radius`     | int            | 2       | `RunRadius`      |
| `run_length`     | int            | 1000    | `RunLength`      |
| `run_candidate`  | int            | 10      | `RunCandidate`   |
| `run_direction`  | `"Re"`/`"Im"`/`"NegRe"`/`"NegIm"` | `"NegIm"` | `RunDirection` |

Precision values are **decimal digits** at the API; conversion to
binary bits happens internally
(see [`docs/INVARIANTS.md`](INVARIANTS.md) §1).

---

## AMF options

`"amf_options"` is the per-recursion-tree configuration for the AMFSystem
pipeline used by `solve_integrals` / `black_box_amflow`.  All fields are
optional.

| Key | Type | Default | Meaning |
|---|---|---|---|
| `amf_modes`     | string[] | `["Prescription", "Mass", "Propagator"]` | η-injection candidate modes (mirrors `AMFMode`) |
| `ending_schemes`| string[] (`"Tradition"`/`"Cutkosky"`/`"SingleMass"`) | `["Tradition", "Cutkosky", "SingleMass"]` | Ending schemes tried in order |
| `cache_root`    | string   | empty | Root for AMF subsystem caches; per-node directories are created underneath |
| `direction`     | string   | upstream default | Override path direction for this run |
| `max_recursion_depth` | int | unlimited | Hard cap on AMF subsystem recursion |
| `blackbox` (alias: `bb`) | object | (see below) | Kira/Fermat/IBP knobs |

### `amf_options.blackbox`

| Key | Type | Default | Meaning |
|---|---|---|---|
| `numeric_values`   | object map<string,string\|number> | `{}` | Numeric values for kinematic invariants — *required for non-trivial families* (e.g. `{"s": "100", "t": "-1"}`).  **Real scalars only**; the complex object form `{"re":..,"im":..}` is rejected with a clear error (out of scope). |
| `ibp_rank`         | int | upstream default | `BlackBoxRank` floor |
| `ibp_dot`          | int | upstream default | `BlackBoxDot` floor |
| `n_thread`         | int | upstream default | `NThread` for Kira |
| `integral_order`   | int | 5 | Kira `integral_ordering` |
| `kira_executable`  | string | `/usr/local/bin/kira` | Path to `kira` |
| `fermat_executable`| string | `/usr/share/Ferl7/fer64` | Path to `fer64` |
| `work_dir`         | string | per-mode default | Override work dir at this level |
| `log_file`         | string | `<kira_working_dir>/kira_stdout.log` | Where Kira stdout/stderr is written (per-run) |
| `run_firefly`      | bool | `false` | Route the Reduce-mode IBP solve through FireFly finite-field reconstruction (the memory-light fix for aux-mass+η coefficient swell) |
| `firefly_bunch`    | int | `0` | FireFly bunch size |
| `kira3_truncate_sp`| bool | `false` | Kira 3 heuristic seed truncation `truncate_sp {l=#loops}`.  **Per-family**: set `true` for the 3-loop crossed-box / two-mass 5-point-style nonplanar tops (collapses 473→39 masters and avoids the abort we observed in the Fermat path on that family); leave `false` for box1-style families (it destabilizes the closure loop, r:5→r:48).  Also overridable via `AMFLOW_KIRA3_TRUNCATE_SP=0/1`. |
| `kira3_sectorwise` | bool | `false` | Kira 3 `run_triangular: sectorwise` + `iterative_reduction: sectorwise`.  Did not terminate on the aux-mass DE reduce in our tests; opt-in only.  Also via `AMFLOW_KIRA3_SECTORWISE=0/1`. |

The launcher also LD_PRELOADs `libjemalloc.so.2` into the Kira child by default
(belt-and-suspenders; Kira is now also *linked* with jemalloc).  Set
`AMFLOW_KIRA_NO_JEMALLOC=1` to disable, or `AMFLOW_JEMALLOC_PATH` / `MALLOC_CONF`
to override.

---

## Rational entries

Used in `mode: amflow` (matrix entries).  Two equivalent forms:

```jsonc
// (a) polynomial only, denominator implicitly 1
{ "coeffs": ["c0", "c1", "c2", ...] }

// (b) explicit numerator / denominator polynomials in η
{ "num": ["a0", "a1", ...], "den": ["b0", "b1", ...] }
```

Each coefficient is parsed as an `fmpq` (rational); strings like
`"3/7"`, `"-1"`, integers, and decimals are all accepted.  Index `i`
in the array is the coefficient of `η^i`.

Example:

```jsonc
// 1 / η  →  num = 1, den = 0 + 1·η = η
{ "num": ["1"], "den": ["0", "1"] }

// (3 + 5η - η²) / 1
{ "coeffs": ["3", "5", "-1"] }
```

---

## Complex values

Used everywhere a single complex number appears (boundary `mu` /
`value`, `eps_samples` entries, output).  Two equivalent forms:

```jsonc
// (a) explicit real / imag parts (each parsed as fmpq or arb decimal)
{ "re": "1/2", "im": "0" }

// (b) bare real scalar — imaginary part = 0
"1/2"
"-1.5"
3
```

Output complex values always use form (a).

---

## Config hygiene — unrecognized keys warn loudly

Every config object the dispatcher parses (top level per mode, `options`,
`amf_options`, `amf_options.blackbox`/`bb`, `family`, per-integral entries)
checks its keys against the set the parsing code actually reads.  A key that
is not read is **never silently ignored**:

- **Default**: one stderr warning per (object, key) per process; parsing and
  the run continue unchanged (legacy configs keep running, stdout output is
  byte-identical).  If the stray key is a known knob of *another* object, the
  warning names where it WOULD be read, e.g.

  ```
  [amflow][config] WARNING: 'eps_parallel' found in options — not read; use amf_options.eps_parallel or AMFLOW_EPS_PARALLEL
  ```

  For keys unknown anywhere, the warning lists the accepted keys of the
  object instead.

- **`AMFLOW_STRICT_CONFIG=1`** (any non-empty value other than `0`):
  upgrades every such warning to a fatal `std::runtime_error` at parse time,
  before any compute or work-dir creation.

- Keys beginning with `_` are reserved for comments and never flagged.

Motivation: a knob in the wrong JSON object (the `eps_parallel`-in-`options`
case) previously parsed cleanly, did nothing, and wasted a long serial run.
When adding a knob to the code, update the accepted-key list at the parse
site **and** `config_key_locations()` in `src/api/run_json.cpp` in the same
edit.

## Error handling

The dispatcher throws `std::runtime_error` on any schema or parse
failure; `amflow_cli` catches and writes the exception text to stderr,
exiting non-zero.  Typical errors:

| Message starts with | Cause |
|---|---|
| `unknown mode '...'`           | `mode` is not one of the three known values |
| `missing required array '...'` | Required key missing |
| `... must be an array`         | Wrong JSON type for an array-shaped key |
| `... has N indices, expected M`| Integral indices length ≠ propagators length |
| `boundaries length must match matrix size; got N vs M` | Boundary list cardinality mismatch in `mode: amflow` |
| `failed to parse fmpq from '...'` | Coefficient string isn't a valid rational |

---

## See also

- [`docs/USER_GUIDE.md`](USER_GUIDE.md) — step-by-step walkthrough.
- [`docs/REFERENCE_MAP.md`](REFERENCE_MAP.md) — Mathematica → C++ symbol map (which upstream functions back which JSON modes).
- [`examples/`](../examples) — five runnable JSON inputs.
- [`src/api/run_json.cpp`](../src/api/run_json.cpp) — implementation; this document and that file must stay in sync.
