# MMA Parity Audit

> 中文版本: [AUDIT_MMA_PARITY_zh.md](AUDIT_MMA_PARITY_zh.md)

Status reference for the C++17 port against upstream Mathematica
AMFlow (snapshot at commit `efda1db` of
<https://gitlab.com/multiloop-pku/amflow>).

## Status

- **228 of 228 oracle benches under [`tools/bench/`](../tools/bench/)
  match upstream at rel ~10⁻³⁰**, except where the integral's
  intrinsic cancellation horizon limits precision.
- All identified semantic divergences are resolved in code or
  intentionally out of scope (see §3 below).
- No outstanding unverified branches — every public symbol port has
  at least one oracle exercising it, or a unit test pinning its
  behaviour.

---

## 1. Methodology

For each upstream file, an audit pass classified every public symbol
into one of:

- **🟢 verified** — semantics match upstream, oracle covers the path.
- **🟡 unverified** — semantics look correct but no oracle exercises
  the path.  Latent risk; resolve by adding an oracle or unit test.
- **🔴 divergent** — actual semantic difference.  Resolve in code or
  document as out-of-scope.
- **⚪ not ported** — explicitly out of scope (gauge link, HQET, SCET,
  Wilson, non-Kira IBP backends, WSL plumbing).

Three audit passes ran in parallel — one for each upstream file
domain (AMFlow.m core, Kira/interface.m, DESolver.m + dependencies)
— with line-by-line cross-referencing against the port.

To run a future audit pass (after upstream advances or new symbols
are introduced):

1. Pin the upstream commit hash.
2. For each public symbol, locate the corresponding C++ port.
3. Walk the call graph and assert oracle coverage.
4. Either land a new oracle or unit test, or document the gap.

---

## 2. Not ported (intentional, by design)

These upstream features are **deliberately absent** from the C++
port; no plan to add them.

- `SolveIntegralsGaugeLinkSingle` / `SolveIntegralsGaugeLink`
  (gauge-link / HQET / SCET / Wilson lines).
- `ExpandGaugeX`, `GenerateSquare` (linear-propagator reshaping for
  gauge links).
- `IBPReducer = "FiniteFlow+LiteRed"` and other non-Kira backends
  (`FIRE`, `LiteRed`, `Blade`) — port is Kira-only by design.
- `WSL` Windows compatibility plumbing and `BinaryFormat` Windows
  newline handling.
- `$PermutationOption` (Kira yaml field, never set in upstream
  benches).
- `FireFly` / `Mixed` / `NoFactorScan` Kira reduction modes (only
  `Masters` and `Kira` modes are used).
- `install.m` discovery scripts (handled by CMake).
- `ComplexMode` / `CompensateRule` real-only filter for complex
  numeric kinematics.  All numeric_values must be real scalars; the
  object form `{"re":..,"im":..}` is rejected explicitly in
  `api::run_json::apply_blackbox_options`.

---

## 3. Out-of-scope user-visible behaviour

Complex-valued numeric kinematics (see ⚪ above) is the one
intentional gap with user-facing semantics: the JSON parser rejects
the complex-form object loudly rather than truncating to the real
part or returning a wrong answer.

---

## 4. Upstream drift watch

Diffs against the pinned reference snapshot (`efda1db`):

- **New `UseCache` and `SkipReduction` options** at `AMFlow.m:259` —
  AMFSystem persistence caching/skip toggles.  Not ported; not
  blocking on any current functionality.
- **Source-comment line-number drift** of ~30 lines from upstream's
  growth since the port (~3 %).  Pin future MMA citations to a
  specific upstream commit hash to keep references stable.

All upstream-side changes since the port are **additions**; no
upstream behaviour we already ported has changed.

---

## 5. Forward-looking work

See [`docs/ROADMAP.md`](ROADMAP.md) for diversity-driven oracle
expansion and other planned work.
