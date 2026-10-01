# Future Optimization Proposals

> 中文版本: [future_optimization_proposals_zh.md](future_optimization_proposals_zh.md)

Forward-looking proposals only.  Each is **unverified and unimplemented**.
Drop entries once accepted-and-implemented or rejected; do not retain
their post-mortems here.

---

## Reuse `sectormappings/` across `system_X_diffeq` runs of the same family

Kira's `sectormappings/` tables (sector symmetries, IBP & LI relations,
trivial-sector detection) depend on the propagator family's structural
form (loop momenta + external momenta + masses) but **not** on which
propagator carries the auxiliary mass η.  So the tables built during
`amf/system_0_diffeq/masters_preheat/` could be reused by
`amf/system_1_diffeq/`, `system_2_diffeq/`, ... instead of having
each Kira invocation rebuild them from scratch.

**Where the change would land**:
- `src/ibp/kira_yaml.cpp::kira_write_config` — initializes the
  sectormappings dir per Kira call.
- `src/pipeline/amfsystem.cpp::AMFSystem::build_diffeq` — per-system
  loop where successive `system_X_diffeq` runs are launched.

For `system_0`: materialize sectormappings as today.
For `system_1+`: point Kira at the existing sectormappings dir
(symlink, copy, or Kira config), bypassing rebuild.

**Before implementing, confirm**:
1. **Correctness**: the η-independence claim.  Propagator lists across
   `system_X_diffeq` differ only by which propagator has `-eta` added.
   Plausible that symmetry/sector relations are unaffected, but
   double-check with a Kira-internals reader.
2. **Whether MMA already does this implicitly**: MMA reuses a single
   Kira working dir, so it may get the speedup for free.  The C++
   sibling-dir pattern (introduced as a Kira-2.x workaround at
   `reduce.cpp` — `masters_preheat/` and `target_reduce/` separated)
   may be the only thing preventing it on C++ side.  Revisit whether
   the Kira-2.x bug still applies on the current Kira version; if
   not, dir-reuse is the smaller fix.

Per `docs/ROADMAP.md` §"Performance benchmarking": perf-only changes
wait until after correctness validation and a dedicated benchmark host.
