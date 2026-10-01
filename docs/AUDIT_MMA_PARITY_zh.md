# MMA Parity Audit

> 英文版本: [AUDIT_MMA_PARITY.md](AUDIT_MMA_PARITY.md)

C++17 port 相对上游 Mathematica AMFlow 的状态参考 (快照在 <https://gitlab.com/multiloop-pku/amflow> 的 commit `efda1db`).

## 状态

- **[`tools/bench/`](../tools/bench/) 下 228 / 228 oracle bench 与上游一致到 rel ~10⁻³⁰**, 除非积分本身的取消尺度限制了精度.
- 所有识别出来的语义 divergence 都已在代码里解决, 或有意 out of scope (见下文 §3).
- 没有未验证的分支 — 每个公开符号 port 都至少有一个 oracle 覆盖, 或有一个 unit test 钉住其行为.

---

## 1. 方法论

对每个上游文件, audit pass 把每个公开符号归类为:

- **🟢 verified** — 语义与上游匹配, oracle 覆盖路径.
- **🟡 unverified** — 语义看起来正确但没有 oracle 走过该路径. 潜在风险; 通过加 oracle 或 unit test 解决.
- **🔴 divergent** — 实际语义差异. 在代码里解决或文档化为 out-of-scope.
- **⚪ not ported** — 显式 out of scope (gauge link, HQET, SCET, Wilson, 非 Kira IBP 后端, WSL 适配).

三个 audit pass 并发跑 — 每个上游文件 domain 一个 (AMFlow.m core, Kira/interface.m, DESolver.m + dependencies) — 逐行对照 port.

未来 audit pass (上游推进或新符号出现后) 的流程:

1. 钉上游 commit hash.
2. 对每个公开符号, 找到对应的 C++ port.
3. 走 call graph, 断言 oracle 覆盖.
4. 要么 land 一个新 oracle / unit test, 要么文档化 gap.

---

## 2. Not ported (有意, 按设计)

下列上游功能**有意**在 C++ port 里缺失; 没有加入计划.

- `SolveIntegralsGaugeLinkSingle` / `SolveIntegralsGaugeLink` (gauge-link / HQET / SCET / Wilson lines).
- `ExpandGaugeX`, `GenerateSquare` (gauge link 的 linear-propagator reshape).
- `IBPReducer = "FiniteFlow+LiteRed"` 以及其它非 Kira 后端 (`FIRE`, `LiteRed`, `Blade`) — port 设计上是 Kira-only.
- `WSL` 的 Windows 适配以及 `BinaryFormat` 的 Windows 换行处理.
- `$PermutationOption` (Kira yaml 字段, 上游 bench 从未设置).
- `FireFly` / `Mixed` / `NoFactorScan` Kira reduction mode (只用 `Masters` 和 `Kira` mode).
- `install.m` 发现脚本 (由 CMake 处理).
- `ComplexMode` / `CompensateRule` 复数 kinematic 的实部过滤. 所有 numeric_values 必须是实数标量; 对象形式 `{"re":..,"im":..}` 在 `api::run_json::apply_blackbox_options` 显式拒绝.

---

## 3. 面向用户可见的 out-of-scope 行为

复数值数值 kinematics (见 ⚪) 是面向用户语义的唯一有意 gap: JSON parser 高声拒绝复数对象形式, 而不是截到实部或返回错答案.

---

## 4. 上游漂移监控

相对钉住的参考快照 (`efda1db`) 的 diff:

- **新 `UseCache` 与 `SkipReduction` 选项** 在 `AMFlow.m:259` — AMFSystem 持久化 caching / 跳约化的开关. 未 port; 不阻塞当前任何功能.
- **源注释行号漂移** ~30 行, 因 port 之后上游增长 (~3 %). 未来 MMA 引用应钉到特定 commit hash 保持引用稳定.

所有 port 之后上游侧的变化都是**新增**; 我们已 port 的上游行为没变.

---

## 5. 前向工作

diversity 驱动的 oracle 扩展和其它计划工作见 [`docs/ROADMAP.md`](ROADMAP.md).
