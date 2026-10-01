# User Guide

> 中文版本: [USER_GUIDE_zh.md](USER_GUIDE_zh.md)

This guide walks through using **AMFlow.cpp** end-to-end:

1. What it computes and when to use it
2. Picking a mode
3. A no-IBP example you can run in 30 seconds (no Kira required)
4. A full Feynman-integral example (requires Kira + Fermat)
5. Reading the output
6. Where to go next

For the algorithm itself, refer to the upstream paper:
[Liu & Ma, *Comput. Phys. Commun.* **283** (2023) 108565](https://doi.org/10.1016/j.cpc.2022.108565).

---

## 1. What does AMFlow.cpp compute?

Given a Feynman integral `I(ε)` (or a list of master integrals), it
returns the value as an ε-Laurent series

$$
I(\varepsilon) = \sum_{k \ge k_0} c_k \, \varepsilon^k
$$

to a user-requested precision, by the **auxiliary mass flow** algorithm
of Liu & Ma.  Internally:

1. The integrals are reduced to a closed set of master integrals
   (this project ships a Kira backend).
2. A differential equation `dM/dη = A(η, ε) · M` is constructed in an
   auxiliary mass `η`.
3. The equation is solved from η = ∞ (where boundary conditions are
   trivial: vacuum bubbles, single-mass tadpoles) down to η = 0
   (the physical limit).
4. Sub-systems are recursively built and solved at each region of the
   asymptotic expansion.
5. Per-ε numerical values are fitted into a Laurent expansion.

You don't need to know any of step 2-5 to use the tool, but if you
want a deeper picture see [`docs/ARCHITECTURE.md`](ARCHITECTURE.md).

---

## 2. Three usage modes

`amflow_cli` consumes a single JSON document.  The top-level `"mode"`
key picks one of three workflows:

| Mode | When to use | What you supply | What you get back |
|---|---|---|---|
| `amflow`           | You already have a rational ODE matrix `A(η)` and an asymptotic boundary spec at η = ∞.  Useful for testing the ODE engine in isolation, or for inputs computed by other tools. | A square rational matrix and a list of `(μ, value)` boundary pairs per row | The vector `I(η = 0)` |
| `solve_integrals`  | You have a Feynman family + target integrals + ε-precision goal, and you want a Laurent expansion. | Family config (loops, legs, propagators, kinematics), target integrals, `goal_digits`, `eps_order` | Each target as `Σ c_k ε^k` |
| `black_box_amflow` | Same as above but you supply the ε grid yourself (no Laurent fit). | Family config + targets + `eps_samples` | Each target evaluated at every ε sample |

Most users want **`solve_integrals`** — it does the full job.  Use
`black_box_amflow` for parity comparisons against Mathematica AMFlow at
single ε values.  Use `amflow` (the raw mode) only if you have the ODE
matrix already in hand.

For complete field-level reference see
[`docs/JSON_SCHEMA.md`](JSON_SCHEMA.md).

> **Out of scope** — **complex-valued numeric kinematics** are not
> supported.  Supply only real values in
> `amf_options.blackbox.numeric_values`; the dispatcher rejects the
> object form `{"re":..,"im":..}` with a clear error.  Intentional
> non-goals are listed in
> [`docs/AUDIT_MMA_PARITY.md`](AUDIT_MMA_PARITY.md); user-facing
> non-goals in [`docs/FAQ.md`](FAQ.md) "What's *not* implemented?".

---

## 3. 30-second example — no IBP, no Kira

The file [`examples/power_law.json`](../examples/power_law.json)
encodes a 1×1 ODE that admits a closed-form solution:

```math
\frac{dI(\eta)}{d\eta} = \frac{1}{\eta} I(\eta), \qquad
I(\eta) \sim \eta^1 \;\; (\eta \to \infty).
```

The exact answer is `I(η) = η`, so `I(0) = 0`.  The JSON:

```jsonc
{
  "options": {
    "x_order":       30,
    "extra_x_order": 10,
    "silent_mode":   true
  },
  "matrix":     [[ { "num": ["1"], "den": ["0", "1"] } ]],
  "boundaries": [[ { "mu":    { "re": "1", "im": "0" },
                     "value": { "re": "1", "im": "0" } } ]],
  "mode": "amflow"
}
```

What every key means:

- `mode: "amflow"` — raw ODE solve, no family, no Kira.
- `matrix` — a `1 × 1` array.  The single entry is the rational
  function `1/η` written as numerator polynomial `[1]` (= 1) divided
  by denominator polynomial `[0, 1]` (= η).  The array index is the
  coefficient of η^index, so `[0, 1]` means `0 + 1·η`.
- `boundaries` — one boundary list per row.  The outer array has one
  entry (the only row); the inner list has one term:
  `value · η^μ = 1 · η^1`, asserting the leading behaviour at η = ∞.
- `options.x_order = 30` — keep 30 terms in series expansions.
- `options.silent_mode = true` — suppress progress logging.

Build the project (assumes you've already done `cmake -S . -B build`):

```bash
cmake --build build -j32
```

Run:

```bash
./build/src/cli/amflow_cli examples/power_law.json
```

> Tip — `sudo cmake --install build` puts `amflow_cli` on `$PATH`
> (default prefix `/usr/local/`), so all subsequent commands in this
> guide can drop the `./build/src/cli/` prefix and just say
> `amflow_cli ...`.  See [`FAQ.md`](FAQ.md) "How do I install the
> library and the `amflow_cli` binary system-wide?" for the full
> install workflow.

Output (abridged):

```jsonc
{
  "mode": "amflow",
  "result": [
    { "re": "0", "im": "0" }
  ],
  "options": { "x_order": 30, ... }
}
```

The `result` array has one entry — the value `I(η = 0) = 0` — exactly
matching the closed form.

[`examples/constant.json`](../examples/constant.json) is even
simpler: dI/dη = 0, BC `I ~ 42·η^0` → output `42`.

You now have a working end-to-end run with the C++ ODE engine.

---

## 4. Full example — `solve_integrals` with IBP

`solve_integrals` is the headline workflow.  It needs:

1. A Feynman family description.
2. A list of target master integrals.
3. Precision goals (`goal_digits`, `eps_order`).
4. Numeric values for every kinematic invariant.

### 4.1 Prerequisite — install Kira and Fermat

The IBP backend is **Kira**.  Without it, `solve_integrals` and
`black_box_amflow` cannot run.

- Kira: <https://kira.hepforge.org/> — install to `/usr/local/bin/kira`
  (or override the path via `amf_options.blackbox.kira_executable`).
- Fermat: <https://home.bway.net/lewis/> — install to
  `/usr/share/Ferl7/fer64` (or override via
  `amf_options.blackbox.fermat_executable`).

Both tools are *runtime* dependencies, not build dependencies.  The
build itself does not link against either; the project shells out to
them when an IBP step is needed.

### 4.2 The example: 1-loop massive bubble

The bubble integral

```math
\mathcal{B}(s) = \int \frac{\mathrm{d}^d \ell}{(2\pi)^d} \;
                 \frac{1}{(\ell^2 - m^2)\, ((\ell - p)^2 - m^2)},
\quad s = p^2
```

is the simplest example of a non-trivial Feynman integral with a known
analytic answer.

[`examples/bubble_solve_integrals.json`](../examples/bubble_solve_integrals.json):

```jsonc
{
  "mode": "solve_integrals",
  "options": {
    "chop_pre":        20,
    "silent_mode":     true
  },
  "family": {
    "name":  "bubblefam",
    "loops": ["l"],
    "legs":  ["p"],
    "replacement": { "p^2": "s" },
    "propagators": [
      "l^2 - msq",
      "(l - p)^2 - msq"
    ]
  },
  "integrals":   [ { "indices": [1, 1] } ],
  "goal_digits": 25,
  "eps_order":   2,
  "work_dir":    "/tmp/desolver_example_bubble",
  "amf_options": {
    "blackbox": {
      "numeric_values": { "s": "9", "msq": "1" }
    }
  }
}
```

Key-by-key:

- **`family`** — describes the topology.
  - `loops` / `legs` — the loop momentum `l` and external `p`.
  - `replacement` — kinematic invariants.  Here `p² = s`.
  - `propagators` — two scalar denominators in terms of `l`, `p`,
    and the symbolic mass parameter `msq`.
- **`integrals`** — the targets.  Each entry has `indices` of length
  equal to `propagators`.  `[1, 1]` is the master integral with
  unit power on each denominator.
- **`goal_digits: 25`** — fit each ε-coefficient to ~25 decimal digits.
- **`eps_order: 2`** — for this 1-loop bubble (L=1) it fits up to `ε^(2 − 2L) = ε^0`, so the result will contain ε⁻¹ and ε⁰ coefficients (the formally-leading ε⁻² vanishes for this finite bubble and is trimmed). General rule: pass `eps_order = K + 2L` if you want the result up to ε^K.
- **`work_dir`** — Kira intermediate files go here.  Reusable across
  runs (cache reuse is fingerprint-guarded).
- **`amf_options.blackbox.numeric_values`** — assign numeric values
  to *every* symbol that appeared in `replacement` and to every
  symbolic mass that appeared in `propagators`.  Here `s = 9` and
  `msq = 1`.

### 4.3 Running it

```bash
./build/src/cli/amflow_cli examples/bubble_solve_integrals.json
```

Expected runtime: a few seconds on a modern desktop.  The Kira workdir
at `/tmp/desolver_example_bubble/` will accumulate yaml files,
intermediate results, and a `kira_target.m` per IBP pass.

Output (abridged):

```jsonc
{
  "mode": "solve_integrals",
  "result": [
    {
      "integral": { "family": "bubblefam", "indices": [1, 1] },
      "leading_order": 0,
      "coefficients": [
        { "order": 0, "value": { "re": "...", "im": "..." } },
        { "order": 1, "value": { "re": "...", "im": "..." } },
        { "order": 2, "value": { "re": "...", "im": "..." } }
      ]
    }
  ],
  "options": { ... }
}
```

The Laurent expansion is

`B(s, ε) ≈ c_0 + c_1 · ε + c_2 · ε²`

where `c_k = result[0].coefficients[k].value`.

### 4.4 Computing in non-4 spacetime dimensions

By default the host space-time dimension is `D = 4 - 2ε` (so `ε → 0`
recovers the physical four-dimensional value).  To work in a different
dimension scheme — for instance `D = 7/3 - 2ε` — set
`options.d0 = "<rational>"`:

```jsonc
{
  "mode": "solve_integrals",
  "options": {
    "d0": "7/3",          // any rational "p/q" or integer; default "4"
    "rationalize_pre": 100
  },
  "family":     { ... },
  "integrals":  [ ... ],
  "goal_digits": 30,
  "eps_order":   4,
  "amf_options": { "blackbox": { "numeric_values": { "s": "100", "t": "-1" } } }
}
```

What happens internally: the user-facing ε grid stays user-facing,
but the engine evaluates the family at `ε + (4 − D₀)/2` (mirroring
upstream `AMFlow.m:1342`/`1351`) and fits the Laurent expansion against
the original grid (mirroring `AMFlow.m:1356`).  This keeps `D₀`
invisible at the API while letting it cancel cleanly inside the
algorithm.

A committed end-to-end oracle benchmark
[`tools/bench/box1_d0_7_3_solve_integrals_cpp.json`](../tools/bench/box1_d0_7_3_solve_integrals_cpp.json)
solves a 1-loop box at `D₀ = 7/3, s = 100, t = -1` and matches
Mathematica AMFlow to ~30 significant digits across orders 0..2.

### 4.5 Adapting to your own integral

To solve a different integral, edit:

1. **`family.propagators`** — your topology.
2. **`family.replacement`** — your kinematics (everything that's not a
   loop momentum and isn't `i0`).
3. **`integrals`** — your targets, one entry per master.
4. **`amf_options.blackbox.numeric_values`** — concrete numbers for
   every invariant in `replacement` and every bare mass that appeared
   in `propagators` but not in `replacement`.
5. **`goal_digits`** / **`eps_order`** — set as you wish; cost grows
   roughly linearly in `goal_digits` and exponentially with the number
   of internal masses + invariants.

Larger families ([`examples/box1_solve_integrals.json`](../examples/box1_solve_integrals.json)
is a 1-loop box with two invariants `s, t`) follow the same template.

---

## 5. Reading the output

Every run produces a top-level object:

```jsonc
{
  "mode":    "<echoed mode name>",
  "result":  <mode-specific shape>,
  "options": <effective options after applying input.options>
}
```

Each complex number is encoded as `{ "re": "<decimal string>",
"im": "<decimal string>" }`.  Decimal strings preserve exact precision
out to the working precision (`options.working_pre` digits).  Parsing
with Python:

```python
import json, mpmath
mpmath.mp.dps = 100        # match working_pre

out = json.load(open("/tmp/bubble_out.json"))
for entry in out["result"]:
    integral = entry["integral"]
    coefficients = entry["coefficients"]
    print(f"\n{integral['family']}{integral['indices']}:")
    for c in coefficients:
        order = c["order"]
        val = mpmath.mpc(c["value"]["re"], c["value"]["im"])
        print(f"  ε^{order}: {val}")
```

For details on every output field, see
[`docs/JSON_SCHEMA.md`](JSON_SCHEMA.md).

---

## 6. Going further

| Want to... | Read |
|---|---|
| Look up a specific input/output field | [`docs/JSON_SCHEMA.md`](JSON_SCHEMA.md) |
| Understand the algorithm and library design | [`docs/ARCHITECTURE.md`](ARCHITECTURE.md) |
| Check parity status and validated families | [`docs/AUDIT_MMA_PARITY.md`](AUDIT_MMA_PARITY.md) |
| Find which Mathematica function backs a C++ entry | [`docs/REFERENCE_MAP.md`](REFERENCE_MAP.md) |
| Solve a build / runtime problem | [`docs/FAQ.md`](FAQ.md) |
| Regenerate Mathematica reference data | [`reference/README.md`](../reference/README.md) and [`tools/bench/README.md`](../tools/bench/README.md) |
| Contribute | [`CONTRIBUTING.md`](../CONTRIBUTING.md) and [`docs/CONTRIBUTING.md`](CONTRIBUTING.md) |
| Cite this software | [`CITATION.cff`](../CITATION.cff) |

For questions not covered here or in the FAQ, open an issue at
<https://github.com/chang18/amflow-cpp/issues> or contact the
maintainer at <3250800970@qq.com>.
