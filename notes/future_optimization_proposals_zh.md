# 未来优化提议

> 英文版本: [future_optimization_proposals.md](future_optimization_proposals.md)

只放前向提议. 每一条**未验证, 未实现**. 一旦被采纳并实现或被拒绝, 就删掉; 不在这里保留事后分析.

---

## 在同 family 的 `system_X_diffeq` 之间复用 `sectormappings/`

Kira 的 `sectormappings/` 表 (sector symmetry, IBP & LI 关系, trivial-sector 检测) 依赖 propagator family 的结构形式 (loop 动量 + 外动量 + 质量), **但**不依赖哪个 propagator 携带 auxiliary mass η. 所以 `amf/system_0_diffeq/masters_preheat/` 期间建出来的表, 可以被 `amf/system_1_diffeq/`, `system_2_diffeq/`, ... 复用, 而不是让每次 Kira 调用都从头重建.

**改动会落在哪**:
- `src/ibp/kira_yaml.cpp::kira_write_config` — 每次 Kira 调用初始化 sectormappings 目录.
- `src/pipeline/amfsystem.cpp::AMFSystem::build_diffeq` — 每个 system 的循环, 顺序启动 `system_X_diffeq`.

`system_0`: 像今天一样物化 sectormappings.
`system_1+`: 让 Kira 指向已有的 sectormappings 目录 (symlink, copy, 或 Kira 配置), 绕过重建.

**实现前要确认**:
1. **正确性**: η 独立性这个论断. 各 `system_X_diffeq` 间 propagator 列表只差在哪条 propagator 加了 `-eta`. symmetry / sector 关系不受影响是合理的, 但请用读过 Kira internals 的人 double-check.
2. **MMA 是否隐式已经这样做了**: MMA 复用单一 Kira 工作目录, 可能免费拿到这个加速. C++ 用同级目录模式 (作为 Kira-2.x 临时方案在 `reduce.cpp` — `masters_preheat/` 与 `target_reduce/` 拆开) 可能是 C++ 侧唯一阻止这件事的因素. 复查 Kira-2.x 的 bug 在当前 Kira 版本上是否还成立; 若已不成立, 目录复用是更小的修复.

按 `docs/ROADMAP.md` §"性能 benchmark": 纯性能改动等到正确性验证和专用 benchmark 主机之后.
