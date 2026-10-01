# 更新日志

> 英文版本: [CHANGELOG.md](CHANGELOG.md)

本项目所有重要的变更都记录在这里.

格式基于 [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), 遵循 [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added
- Oracle suite 已扩到 228 三元组, 覆盖 L = 1, 2, 3, 4 (1L=52, 2L=53, 3L=64, 4L=59). sub-master 复用父 MMA cache 让 4L diversity 成本低.
- Kira yaml 输出对齐 MMA `interface.m:126` (`[Propagator, 0]` 形式, 闭式 rank-1 propagator 重建, 原 caller-order preferred list).
- Opt-in `AMFLOW_STAGE_TIMING` 环境变量: 每次 Kira 调用发 `[kira_time]`, amflow_cli 退出时发 `[total_time]`, 让外部工具能把 wallclock 拆成 Kira-subprocess 与 amflow 内部时间.
- CMake install rules 现在包含 `amflow_cli`: 跑过 `sudo cmake --install build` 后, 二进制落在 `${CMAKE_INSTALL_PREFIX}/bin/amflow_cli` (默认 `/usr/local/bin/`), 与 library 和 headers 一起, CLI 在任何工作目录下都能直接调用, 不再需要 `./build/src/cli/`.

### Fixed
- Pentabox 2L 5-leg Jordan eigenvector normalization: `fmpq_mat_nullspace_exact` rescale basis vector, 使最后非零项 = 1, 匹配 Mathematica 约定.
- multi-mass 拓扑上 `build_boundary` projection: `to_complete_explicit` rank-filter 不再丢带 mass 的 propagator; 数值代入已恢复.
- multi-mass 3L corner-master divergence: `analyze_block` 的 extend / Gather / Complement 现在与 MMA 对齐; 默认 `sparse_chop_digits = max(chop_pre, working_pre - 40)` 把 acb 噪声预算压在合法值之下.

### Changed
- 复数值数值 kinematics 改为 out-of-scope (之前是 "deferred"); JSON dispatcher 拒绝行为不变.


## [1.1.0] — 2026-05-12

post-v1.0 的 audit 驱动正确性 pass + Phase 3 oracle diversity 扩展. 这俩一起关闭了 post-release 行级 MMA parity audit ([`docs/AUDIT_MMA_PARITY.md`](docs/AUDIT_MMA_PARITY.md)) 暴露的所有静默错值路径, 开了项目的五个 "我们当时没想到" oracle diversity 轴.

v1.1.0 之后 audit 表为 **86 🟢 / 0 🟡 / 7 🔴 (6 fixed + 1 D5 无限期 deferred) / 17 ⚪**, oracle 套件覆盖了全部 5 个 diversity 轴 (loop number L ≤ 4, ≥ 3 kinematic invariants, multi-cut Cutkosky, mixed-mass, ε-extremes), 新 oracle 上 rel ≤ ~10⁻³⁰.

### Added — Phase 3 oracle diversity 扩展
- 3.F **L=4 banana oracle** (`tools/bench/banana_4loop_eps001_*`): 4-loop 等质量 banana sunrise, `psq=-3, msq=1, eps=1/1000`. 暴露并锁住了 audit divergence **D7** (见下 Fixed). 与 MMA 一致到 rel 7.2 × 10⁻³¹ / 1.6 × 10⁻³⁰ (C++ 399 s vs MMA 385 s).
- 3.G **多 invariant electroweak box oracle** (`tools/bench/ewbox_1loop_eps001_*`): 1-loop electroweak box, 交替 W/Z 质量, 4 个独立 kinematic invariant {`s, t, mWsq, mZsq`}. 一致到 rel ≤ 3.8 × 10⁻³⁰ (C++ 39 s vs MMA 62 s). 阈值之上的 WW bubble 给出物理虚部, 锁住跨质量 branch-cut 处理.
- 3.H **多 cut Cutkosky oracle** (`tools/bench/cutbanana_4L_eps001_*`): 4-loop massless cutbanana, 全部 5 条内线 on-shell (5-particle Cutkosky cut). 一致到 rel ≤ 2.2 × 10⁻³⁰ (C++ 141 s vs MMA 149 s). 端到端验证了 5 级 subsystem 递归.
- 3.I **混合质量 oracle** (`tools/bench/bn3mix_eps001_*`): 3-loop banana 1 个 W-massive + 3 个 massless 内 propagator — 第一个 3-loop mixed-mass case. 一致到 rel ≤ 2.9 × 10⁻³⁰ (C++ 56 s vs MMA 96 s). 锻炼了混合质量下的 scaleless 子 sector 检测.
- 3.J **ε 极值 oracle** (`tools/bench/cutbubble_1L_eps{2,10000}_*`): cutbubble eps = 1/2 (D = 3) 时与精确值 1/8 一致到 rel 1.9 × 10⁻⁶⁴; eps = 10⁻⁴ 时一致到 rel 1.0 × 10⁻³². 最初 eps = 1 (D = 2) 跑出来揭示了一个上游 MMA AMFlow 的限制 (DESolver 返回部分符号化的 `(1/2π) Im[DESolver\`Private\`variables[1, 1]]`); 改用 eps = 1/2 作为次极端的 rational, 拿到干净数值 — bench 源注释里有记录.

### Fixed (对齐上游)
- `RunningOptions::run_length` 默认从 200 调到 **1000**, 匹配上游 `RunLength = 1000` (`AMFlow.m:259`, `DESolver.m:97`). 防止 pole 密集时 `RunUnit` 过早 abort. (Audit D1.)
- `GlobalOptions::rationalize_pre` 默认从 20 调到 **100**, 匹配上游 `RationalizePre = 100`. 移除了 contour 上 rationalization 步的静默精度收窄. (Audit D2.)
- `ibp::black_box_reduce` 在 Kira 约化返回的 RHS J-integral 不在 master 列表里时, 现在抛 `std::runtime_error`, 镜像上游 `Kira/interface.m:486` 的 `Abort`. 之前会静默丢这一行. (Audit D6.)
- **D3 — Tradition-with-cut boundary projection.** v1.0 在构造 boundary 子 family 时静默丢了父系统的 `Cut`. Phase 1A patch 把它改成 loud abort. **Phase 1B 实现了完整 projection**, 镜像上游 `ReduceBoundary` (`AMFlow.m:790-803`): 每个父 cut prop 过 bare `region.transform.map`; 每个 fam.prop 用 `reduced_replacement` 模匹配 transformed cut prop; 构造新的 `sub_cut`, 断言 `Count[cut, 1]` 不变; 结果送给子 family 的 `qft::FamilyConfig::build`. 一个新 oracle benchmark (`tradcut_phase_2L_eps001_*`) 顶住 — 来自上游 `examples/automatic_phasespace` 的 2-loop Tradition-with-cut, 与 MMA 一致到相对误差 2.68 × 10⁻³⁰. (Audit D3.)
- `branch_to_loop` (`src/qft/region.cpp`) 加了一个防御性 assert: `det(A) = ±1` (常数), `A` 是 loop-redefinition 矩阵. 这样能在未来 `branch_momenta` 的 unit-leading-loop-coefficient 前提条件被放松前, 先抓到缺失的 `|Det|^(4-2eps)` Jacobian 因子 (上游 `AMFlow.m:731-732`), 避免静默生成错误 boundary integrand. (Audit D4.)
- **D7 — dual-Kira-call master-count divergence.** `ibp::reduce` 和 `ibp::diffeq` 现在镜像上游 `BlackBoxReduce` / `BlackBoxDiffeq` 的两次调用模式: 一次 Masters-mode 预热调用 (通过 `select_mandatory_recursively` 做 sector 级枚举), 再一次 Reduce-mode 调用 (针对具体目标用 `select_mandatory_list`), 两者用相同的 `(rank, dot)`. 每次调用各自在子目录 (`<work_dir>/masters_preheat/` 和 `<work_dir>/target_reduce/`) 跑, 因为 Kira 2.x release 拒绝两次调用共享 `$ReductionDirectory` (Masters-mode `run_initiate: masters` 不注册 `-s` 数值替换, 后续 Reduce-mode 调用会 abort: `Kira::update_auxiliary_file: Last Kira run set 0 variables to numeric values, this time you request N`). `ibp::diffeq` 也不再嵌套调 `reduce()` (那会在导数积分上重新 floor `(rank, dot)` — L=4 banana 的失败模式); 它在 `opts_eff` 的 `(rank, dot)` 上 inline 一次 Reduce-mode Kira 调用, 镜像上游 `AnalyticReduction` 从 `IBPSystem` 继承 `IBPRank`/`IBPDot` 全局. 两个函数都在 Reduce-mode master 文件上加了 SubsetQ guard 对照 Masters-mode sector 枚举, 镜像上游的 `If[!SubsetQ[masters, str], Abort["inconsistent masters from Kira"]]`. 由 L=4 banana oracle 锁定. (Audit D7.)

### Added
- [`docs/AUDIT_MMA_PARITY.md`](docs/AUDIT_MMA_PARITY.md): 行级上游 parity audit 报告. 86 🟢 verified / 0 🟡 unverified / 7 🔴 (6 fully fixed, 1 D5 无限期 deferred) / 17 ⚪ not ported. (起步在 v1.0.0 时是 65 🟢 / 21 🟡 / 6 🔴; post-v1.0 pass 关闭了所有 🟡, 修了 7 个 🔴 里的 6 个, 还暴露出 D7 作为第七个 🔴 — 现已也修.)
- [`docs/ROADMAP.md`](docs/ROADMAP.md) §"性能 benchmark": 12 oracle benchmark 上相对 MMA 的 wallclock baseline. Kira 轻的 bench 中位数加速 ~2×; Kira 重的 bench 收敛到 ~1×.
- [`docs/ROADMAP.md`](docs/ROADMAP.md): post-v1.0 开发计划 (Phase 1 实现完整性, Phase 2 oracle 扩展, Phase 3 持续多样化).
- [`tools/bench/run_perf_audit.sh`](tools/bench/run_perf_audit.sh): 重现 perf run 的 shell driver.
- `tools/bench/` 下的新 oracle benchmark:
  - `tradcut_phase_2L_eps001_*` — 第一个 Tradition-with-cut parity 测试 (D3 验收闸).
  - `banana_4loop_eps001_*` — Phase 3.F 第一个 L=4 oracle (D7 闸).
  - `ewbox_1loop_eps001_*` — Phase 3.G 多 invariant oracle.
  - `cutbanana_4L_eps001_*` — Phase 3.H 5-particle Cutkosky oracle.
  - `bn3mix_eps001_*` — Phase 3.I 3-loop mixed-mass oracle.
  - `cutbubble_1L_eps{2,10000}_*` — Phase 3.J ε-extremes oracle.
- 新 ending scheme: `EndingScheme::Trivial`, 自动追加到用户 `ending_schemes` 列表末尾作为最终 fallback (镜像上游 `AMFlow.m:1034`).
- `solve_integrals` 新增单 eps 快速路径: 若 `numeric_values["eps"]` 已提供, 跳过 Laurent 拟合, 直接返回在该 eps 上求值的积分 (镜像 `AMFlow.m:1364-1374`).
- 每个 system 的路径方向: `AMFSystem::setup()` 计算 η-接触的 loop 的 prescription 共识, 覆盖全局 `run_direction` 用于该 system 的 ODE 求解 (镜像 `AMFlow.m:981-991`).
- Cutkosky setup 现在验证 phase-volume 分量质量在 `Numeric` 代入后都非负; 否则抛错 (镜像 `AMFlow.m:1050`).
- `apply_blackbox_options` (`src/api/run_json.cpp`) 现在主动拒绝 `amf_options.blackbox.numeric_values` 里的复数对象形式 `{"re":..,"im":..}`, 用清晰的 "not implemented" 错误信息指向 audit divergence D5. 之前的行为会在 parser 深处用不清楚的信息失败; 显式拒绝让功能 gap 在输入边界处可见.

### 已知限制 (无限期 deferred)
- **D5 — Kira `ComplexMode` / 虚部数值流水线没实现.** 在需要在复数 kinematics (非零虚部) 上求值的用户无法通过本 port 完成; 上游通过 `IBPRule` / `CompensateRule` 过滤这种值, 但 C++ port 把 `numeric_values` 当成纯实数的平直 map. 12 个 oracle benchmark 都用纯实数, 这块 surface 未测. **已重新分类为 out of scope** (见 Unreleased §Changed) — 详见 [`docs/AUDIT_MMA_PARITY.md`](docs/AUDIT_MMA_PARITY.md) §D5.

## [1.0.0] — 2026-05-08

首次公开发布. auxiliary-mass-flow 算法 (上游为 Liu & Ma 的 Mathematica 包 AMFlow, <https://gitlab.com/multiloop-pku/amflow>; Comput. Phys. Commun. 283 (2023) 108565) 的 C++17 重实现. 所有算法上的归属属于上游作者; 本项目仅贡献 C++ 实现与数值 parity 测试 / benchmark 框架.

### Added
- 按 domain 组织的架构, 都在 `amflow::` namespace 下: `numeric`, `algebra`, `ode`, `qft`, `ibp`, `pipeline`, `api`, `cli`.
- 公开头文件在 `include/amflow/<domain>/`, domain 之间是严格的 ABI/PIMPL 边界.
- ODE solver — 上游 `DESolver.m` 的重实现, 覆盖 oracle benchmark 中走过的工作流; 12 个 oracle case 全部相对参考一致到 ~1e-30.
- AMFlow + Kira pipeline — 实现上游 `AMFlow.m` 的核心算法: `amflow`, `black_box_amflow`, `solve_integrals` 三种 mode.
- `amflow_cli` — `amflow::api::run_json` 的 JSON driver, 暴露三种顶层 mode.
- 500+ 个 GoogleTest case, 含 `tools/bench/` 下 12 个采样 parity oracle benchmark (`eps = 1/1000`).
- `tools/math_ref/` 和 `tools/bench/` 下的 Mathematica reference driver (需要上游 AMFlow 的本地 clone, 见 `reference/README.md`), 以及 committed 的 JSON reference 输出用作回归测试.
- CMake install / export 规则; 下游项目可以用 `find_package(AMFlowCpp)` 并链接 `AMFlowCpp::amflow`.

### Notes
- 上游 Mathematica AMFlow 源码**未**随仓库分发. 见 `reference/README.md` 怎么本地拿到以便重新生成 benchmark reference 数据.

### Out of scope
- `SolveIntegralsGaugeLink`, HQET / SCET / Wilson-line 工作流 — 见 [`docs/AUDIT_MMA_PARITY.md`](docs/AUDIT_MMA_PARITY.md).

[1.1.0]: https://github.com/chang18/amflow-cpp/releases/tag/v1.1
[1.0.0]: https://github.com/chang18/amflow-cpp/releases/tag/v1.0
