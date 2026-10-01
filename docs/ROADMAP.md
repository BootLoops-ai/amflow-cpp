# Roadmap

> 中文版本: [ROADMAP_zh.md](ROADMAP_zh.md)

Forward-looking development plan.  For per-release narrative see
[`CHANGELOG.md`](../CHANGELOG.md); for MMA-parity status see
[`docs/AUDIT_MMA_PARITY.md`](AUDIT_MMA_PARITY.md).

---

## Active practice — diversity-driven oracle expansion

The MMA-parity audit is closed: every public symbol port has at
least one oracle exercising it.  Going forward the quality loop is
**oracle diversification** — build new benchmarks along axes the
existing suite under-covers, and fix anything that diverges from MMA
inline.

The 228 committed oracle triplets currently cover:

| Diversity axis | Current coverage |
|---|---|
| Loop number          | L = 1, 2, 3, 4 |
| Kinematic invariants | 1 (`s`), 2 (`s, t`), 4 (`s, t, mWsq, mZsq`) |
| External legs        | 2, 3, 4, 5, 6 (pentagon 1L, pentagon 1L massive, pentabox 2L, hexagon 1L) |
| Mass configurations  | all-massless, all-equal, mixed at L = 1 / 2 / 3 / 4 |
| Cutkosky cuts        | 2-, 3-, 4-, 5-particle cuts |
| ε extremes           | `1/2`, `1/100`, `1/1000`, `1/10000` |
| Sector size          | ≤ 8 propagators (pentabox 2L corner) |

Target additions: higher-loop multi-invariant cases, larger sectors
(≥ 10 propagators), and any topology a downstream user reports.

Each new oracle is roughly one day of work (construct family → MMA
reference → C++ run → commit triplet).  When a new oracle exposes a
divergence, fix it inline; the bug becomes a parity-validated case
going forward.

---

## Out of scope — won't be implemented

Deliberate non-goals.  The JSON entry-point rejects them with a
clear error rather than silently mishandling.

- **`SolveIntegralsGaugeLink`**, HQET / SCET / Wilson-line workflows.
- **IBP backends other than Kira** — no FIRE / LiteRed / FiniteFlow /
  Blade.
- **Complex-valued numeric kinematics**.  The C++ algebra layer is
  over Q (FLINT `fmpz_mpoly_q_t`) and does not carry a complex
  coefficient ring.  The JSON dispatcher rejects the object form
  `{"re":..,"im":..}` in `amf_options.blackbox.numeric_values`;
  provide only purely-real numeric values.  See
  [`docs/FAQ.md`](FAQ.md) and
  [`docs/AUDIT_MMA_PARITY.md`](AUDIT_MMA_PARITY.md).

---

## Update policy

- Each new oracle triplet (cpp.json + mma.wl + mma_reference.json)
  is added to `run_perf_audit.sh` rotation on landing.
- When a new oracle exposes a real divergence, fix in the same
  session if practical and document the relevant assumption inline
  in code.
- Do not silently drop items that turn out to be infeasible —
  document the reason where the user will see it (FAQ or AUDIT).

---

## Performance benchmarking

Wall-clock benchmarking against upstream MMA AMFlow is **intentionally
deferred** until after the project has completed broader correctness
validation and the measurements can be run on a dedicated, load-isolated
machine.  When that pass is performed, the numbers will be added here.
