# 路线图

> 英文版本: [ROADMAP.md](ROADMAP.md)

前向开发计划. 每个 release 的叙事见 [`CHANGELOG.md`](../CHANGELOG.md); MMA-parity 状态见 [`docs/AUDIT_MMA_PARITY.md`](AUDIT_MMA_PARITY.md).

---

## 活跃实践 — diversity 驱动的 oracle 扩展

MMA-parity audit 已 close: 每个公开符号 port 都至少有一个 oracle 覆盖. 后续质量循环是 **oracle diversification** — 沿现有 suite 覆盖不足的轴构造新 benchmark, 跟 MMA 不一致的就 inline 修.

228 个 committed oracle 三元组目前覆盖:

| Diversity 轴 | 当前覆盖 |
|---|---|
| Loop 数         | L = 1, 2, 3, 4 |
| Kinematic invariant | 1 (`s`), 2 (`s, t`), 4 (`s, t, mWsq, mZsq`) |
| 外腿数         | 2, 3, 4, 5, 6 (pentagon 1L, pentagon 1L massive, pentabox 2L, hexagon 1L) |
| 质量配置        | 全 massless, 全等, 在 L = 1 / 2 / 3 / 4 上混合 |
| Cutkosky 切    | 2-, 3-, 4-, 5-particle cut |
| ε 极值         | `1/2`, `1/100`, `1/1000`, `1/10000` |
| Sector size    | ≤ 8 propagator (pentabox 2L 角落) |

未来增量目标: 更高 loop 的多 invariant case, 更大 sector (≥ 10 propagator), 以及下游用户报告的任何 topology.

每个新 oracle 大约一天工作量 (构造 family → MMA 参考 → C++ 运行 → commit 三元组). 新 oracle 暴露 divergence 时, 当场 inline 修复; 这个 bug 之后就成 parity-validated 的 case.

---

## 不在 scope — 不会实现

有意的非目标. JSON 入口会用清晰错误拒绝, 而不是静默错处理.

- **`SolveIntegralsGaugeLink`**, HQET / SCET / Wilson-line 工作流.
- **Kira 以外的 IBP 后端** — 不上 FIRE / LiteRed / FiniteFlow / Blade.
- **复数值数值 kinematics**. C++ algebra 层是 Q 上的 (FLINT `fmpz_mpoly_q_t`), 没有复系数环. JSON dispatcher 拒绝 `amf_options.blackbox.numeric_values` 里的 `{"re":..,"im":..}` 对象形式; 只提供纯实数. 见 [`docs/FAQ.md`](FAQ.md) 与 [`docs/AUDIT_MMA_PARITY.md`](AUDIT_MMA_PARITY.md).

---

## 更新策略

- 每个新 oracle 三元组 (cpp.json + mma.wl + mma_reference.json) landing 时加进 `run_perf_audit.sh` 轮转.
- 新 oracle 暴露真的 divergence 时, 能在同一 session 内修就修, 并把相关假设 inline 写到代码里.
- 不要静默丢掉证实不可行的项目 — 在用户能看见的地方 (FAQ 或 AUDIT) 文档化原因.

---

## 性能 benchmark

相对上游 MMA AMFlow 的 wallclock benchmark **有意延迟**, 等项目完成更广泛的正确性验证, 并且测量能在专用的, 负载隔离的机器上运行之后再做. 到时跑了这一 pass, 数据会加到本文件这一段.
