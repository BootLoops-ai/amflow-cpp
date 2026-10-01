# AGENTS.md — onboarding for AI coding agents

> 中文版本: [AGENTS_zh.md](AGENTS_zh.md)

This file is the canonical onboarding document for any AI coding
agent picking up work on this repo. It captures the project-level
rules and conventions a new collaborator must internalise before
touching code or docs. Deeper documentation lives in [`docs/`](docs/)
and [`docs/AUDIT_MMA_PARITY.md`](docs/AUDIT_MMA_PARITY.md); this file points at it and lists the
load-bearing rules.

---

## Project at a glance

- This repo is **AMFlow.cpp** — a domain-oriented C++17 reimplementation
  of the Mathematica package AMFlow
  (<https://gitlab.com/multiloop-pku/amflow>), the auxiliary-mass-flow
  framework for multi-loop Feynman integrals.  The upstream Mathematica
  source is **not** vendored — see [`reference/README.md`](reference/README.md).
- Layout: domain sub-namespaces `amflow::<numeric|algebra|ode|qft|ibp|pipeline|api|cli>::`,
  sources under `src/<domain>/` and `include/amflow/<domain>/`.
- Mathematica AMFlow is the only behavioural specification. Any C++
  mismatch on a ported branch is a parity bug.
- **This codebase was developed primarily by AI coding agents** (this
  document is the onboarding contract).  The maintainer
  (<3250800970@qq.com>) directs the work but does not line-review;
  the parity contract — the 228 oracle benchmarks under `tools/bench/`
  and the 548-case GoogleTest suite — is the gate every change must
  pass.

## Where to start reading

| Goal | Start here |
|---|---|
| Project intro, build instructions | [`README.md`](README.md) |
| Current parity status, validated families, benchmark inventory | [`docs/AUDIT_MMA_PARITY.md`](docs/AUDIT_MMA_PARITY.md) |
| Domain DAG and data-flow overview | [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md) |
| Rules that must stay true across changes | [`docs/INVARIANTS.md`](docs/INVARIANTS.md) |
| MMA symbol → C++ symbol mapping | [`docs/REFERENCE_MAP.md`](docs/REFERENCE_MAP.md) |
| Adding tests / benches / references; debugging | [`docs/CONTRIBUTING.md`](docs/CONTRIBUTING.md) |

---

## Working rules

All current work is functional / parity work (or, very rarely,
behaviour-preserving code-health work).

### Functional / parity work — primary mode

Triggered by: closing a parity gap, porting a missing branch from
`AMFlow.m`, fixing a wrong numerical result, adding a previously
unported feature.

**Rules:**

- **Source of truth is the upstream Mathematica `*.m` files**
  (<https://gitlab.com/multiloop-pku/amflow>). Read the relevant
  Mathematica passage before designing a fix; cite the upstream
  file + line range in the commit / PR description.
- **Strict mirror — no shortcuts.** Do not substitute closed-form
  special cases, disable features, narrow scope, or paper over a
  divergence. Any non-mirror fix is itself a second divergence under
  a different shape.
- **"Existing tests pass" ≠ "aligned with MMA."** A green ctest only
  proves the previously-exercised paths still work. Validate closure
  against the AMFlow.m line range, not against downstream test
  status alone.
- Breaking long-standing port-level invariants (loop vs. leg
  classification, FamilyConfig shape, etc.) is *already authorised*
  when the alternative is a non-mirror shortcut.

### Code-health / behaviour-preserving work — secondary mode

Triggered by: structural improvements that **must not** alter
numerical output (renames, file splits, dead-code removal,
documentation fixes inside source comments).

**Rules:**

- **Source of truth is the test suite + committed sampled benches.**
  AMFlow.m source-line correspondence is *not* required during the
  edit, but the behaviour anchored by the tests must be preserved
  byte-for-byte.
- **Contract every commit must respect:**
  - `ctest --test-dir build --output-on-failure -j 4` stays green;
  - every committed sampled-bench triplet listed in [`docs/AUDIT_MMA_PARITY.md`](docs/AUDIT_MMA_PARITY.md)
    stays numerically green when re-run;
  - public APIs in `include/amflow/<domain>/` keep their signatures
    unless the change is a deliberate, documented break.
- **One concern per commit.** Don't bundle "rename + split + rewrite"
  — diffs must stay bisectable.
- If a structural change would alter a test expectation or a bench
  output, **stop**. That is functional work, not code-health work.

### Telling the modes apart

| Question | Functional | Code-health |
|---|---|---|
| Could this change alter numerical output? | Yes | No |
| What pins behaviour? | AMFlow.m line correspondence | Test/bench suite |
| Dividing test | Cite an AMFlow.m line range | All tests + benches still pass |

If you cannot tell which mode you are in, ask the user before
proceeding.

---

## Operational conventions

### Build & test

```bash
cmake -S . -B build -DAMFLOW_BUILD_DRIVER=ON
cmake --build build -j32
ctest --test-dir build --output-on-failure -j 4
```

- **Use `-j32` for parallel builds**, not `$(nproc)`. The host shares
  cores with other workloads; saturating CPU starves them.
- **NEVER rebuild only `amflow_cli` and trust the `ctest` you ran
  afterwards.** Running `cmake --build build --target amflow_cli`
  rebuilds the CLI but leaves `amflow_tests` stale; `ctest` will then
  exercise the **old** test binary against the **new** library and
  silently report green when the actual behaviour has changed.
  Either rebuild without `--target` (default builds everything
  including `amflow_tests`), or pass `--target amflow_tests` before
  running `ctest`.  When in doubt, run the bare
  `cmake --build build -j32`.
- Bench commands: see [`tools/bench/README.md`](tools/bench/README.md).
  Each committed sampled bench has a triplet:
  `*_mma.wl` (regenerator) + `*_cpp.json` (driver input) +
  `*_mma_reference.json` (reference values to compare against).

### Dialogue language

- All user-facing dialogue replies must be in **Chinese**
  (本项目所有对话回复一律使用中文).
- Code, comments, commit messages, and file contents (including
  this file and all `docs/*.md`) stay in **English**.

### Doc-sync philosophy

- When changing project state, sync the **whole** `docs/` tree.
- **Trim closed-issue history aggressively.** Resolved bugs,
  abandoned approaches, and "back when X was broken" narratives
  belong in git history, not in living docs.
- Deep diagnostic recipes (specific debug-env-var combos, gotcha
  lists for one-off failure modes) belong in agent-private memory
  or commit messages, not in project docs.

### Bench WL operational gotcha

When writing a top-level `tools/bench/*_solve_integrals_mma.wl`
that mirrors AMFlow.m's `SolveIntegrals` body inline (so a
committed cache directory backs `BlackBoxAMFlow`), three Private-
context symbols MUST be qualified — bare references silently leak
unbound `Global`` symbols and produce hangs / wrong leading order /
collapsed coefficients:

| Bare ref (`Global`) | Must use | Where AMFlow binds it |
|---|---|---|
| `$D0` | `` `AMFlow`Private`$D0` `` | `AMFlow.m:262` (in `Begin["`Private`"]`) |
| `Loop` | `AMFlowInfo["Loop"]` | `AMFlow.m:192` |
| `$Eps` | `` `AMFlow`Private`$Eps` `` (resolves to `Symbol["Global`eps"]`) | `AMFlow.m:181` |

Canonical fix pattern lives in
[`tools/bench/box1_d0_7_3_solve_integrals_mma.wl`](tools/bench/box1_d0_7_3_solve_integrals_mma.wl).

### Acting with care

- Default to **transparently asking** before any hard-to-reverse
  action: force-push, `git reset --hard`, force-amend of published
  commits, dropping tables, killing processes, etc.
- Do NOT skip git hooks (`--no-verify`, `--no-gpg-sign`) unless the
  user explicitly authorises.
- Do NOT commit unless the user explicitly asks. Showing a diff and
  awaiting approval is the safe default.
- Investigate root causes; don't reach for destructive operations
  as shortcuts to make obstacles go away.

---

## What you should NOT do

- **Don't substitute shortcuts** for genuine MMA-aligned fixes
  during functional work.
- **Don't change a test expectation or a bench output** during a
  code-health commit. If you think you must, you are doing
  functional work.
- **Don't run `cmake --build` without `-j32`** unless you have a
  reason that overrides the host-load convention.
- **Don't bundle multiple concerns** into one commit.
- **Don't bypass git hooks** without explicit user authorisation.
- **Don't commit autonomously.** Show the diff, await user
  approval.
- **Don't add speculative abstractions** for hypothetical future
  needs. Three similar lines beats a premature abstraction.
- **Don't add docstrings/comments to code you didn't change.**
  Only comment where the logic isn't self-evident.
- **Don't guess at recent project state** — `git log` and
  `docs/AUDIT_MMA_PARITY.md` are authoritative; check them before
  asserting.
