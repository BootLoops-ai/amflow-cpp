# 架构

> 英文版本: [ARCHITECTURE.md](ARCHITECTURE.md)

本文档描述 AMFlow C++ port 的按 domain 组织的架构.

---

## 库的两种用法

**用法 1 — 直接 ODE 求解.** 用户已经有 rational ODE 矩阵 `A(η)` 和 η = ∞ 处的渐近边界规约, 目标是 `I(η = 0)` 的值. 这条路径完全在 `ode` domain 内实现, 镜像 `DESolver.m` 的 `amflow[de, bc]`:

```
    dI(η)/dη  =  A(η) · I(η)                 (rational matrix)
    I_i(η)   ~  Σ_k v_{i,k} · η^{μ_{i,k}}    near η = ∞
                                               (BoundarySpec)
```

**用法 2 — 完整 AMFlow 工作流.** 用户提供一个 *family* (圈动量, 腿, replacement 规则, propagator) 和一组目标 master integrals; 库自动产出每个 master 的 ε-Laurent 展开. `pipeline` domain 把 `ode` + `qft` + `ibp` 串起来, 完成: Kira IBP, 微分方程构造, region decomposition, 渐近 boundary integrand 构造, 递归子系统构建. 已验证的 parity surface 见 [`AUDIT_MMA_PARITY.md`](AUDIT_MMA_PARITY.md).

---

## Domain DAG

```
                          ┌──────────────┐
                          │   numeric    │   options / log / acb wrap /
                          │              │   rational helpers
                          └──────┬───────┘
                                 │
                          ┌──────▼───────┐
                          │   algebra    │   mpoly, mpoly_matrix,
                          │              │   context migration
                          └─┬─────────┬──┘
                            │         │
                  ┌─────────▼─┐   ┌───▼────────┐
                  │    ode    │   │    qft     │   family_config / family_uf
                  │           │   │            │   topology / amfmode /
                  │           │   │            │   region / boundary / vacuum /
                  │           │   │            │   jintegral / dense_q
                  └─────┬─────┘   └─┬────────┬─┘
                        │           │        │
                        │   ┌───────▼────┐   │
                        │   │    ibp     │   │   libp_deriv, kira (yaml/run/
                        │   │            │   │   parse), reduce
                        │   └─────┬──────┘   │
                        │         │          │
                        │  ┌──────▼──────────▼┐
                        │  │     pipeline    │  factorize, AMFSystem,
                        │  │                 │  black_box_amflow{,_single},
                        │  │                 │  GenerateNumericalConfig,
                        │  │                 │  FitEps, SolveIntegrals
                        │  └────────┬────────┘
                        │           │
                        │     ┌─────▼─────┐
                        │     │    api    │  JSON dispatcher (run_json)
                        │     └─────┬─────┘
                        │           │
                        │     ┌─────▼─────┐
                        │     │    cli    │  amflow_cli executable
                        │     └───────────┘
                        │
                       (用法 1 入口: `ode::amflow(de, bcs)`
                        也可直接从 api/cli 调到)
```

- 每个 domain 是 `amflow::` 下的子 namespace, 同时是 `include/amflow/<domain>/` + `src/<domain>/` 的目录.
- `numeric → algebra` 是基础; `ode` 链和 `qft → ibp → pipeline` 链都建在上面.
- 单测 target `amflow_tests` 覆盖所有 domain; CLI 集成测试以子进程方式跑 `amflow_cli`.

---

## 用法 1 的数据流 (仅 `ode` domain)

```
   ┌─────────────────┐   inf_to_regular   ┌─────────────────┐
   │  asymptotic BC  │ ─────────────────▶ │  numerical BC   │
   │   at η = ∞      │   (ode::inf)       │   at η = x_0    │
   └─────────────────┘                    └────────┬────────┘
                                                   │ regular_run
                                                   │  (ode::regular)
                                                   ▼
                                          ┌─────────────────┐
                                          │  numerical BC   │
                                          │   at η = x_n    │
                                          └────────┬────────┘
                                                   │ solve_asy_exp
                                                   │  (ode::zero, ode::normalize,
                                                   │   ode::jordan, ode::qqbar)
                                                   ▼
                                          ┌─────────────────┐
                                          │  asymptotic     │
                                          │  expansion at   │
                                          │  η = 0          │
                                          └────────┬────────┘
                                                   │ pick_zero_solution
                                                   ▼
                                              I(η = 0)
```

contour `[x_0, x_1, ..., x_n]` 是连接 "近 ∞" 到 "近 0" 的一串 regular 复数点链, 避开所有极点, 由 `ode::path` 的 `run_eta` 生成.

公开入口是 `amflow::ode::amflow(de, bcs, prec)`.

## 用法 2 的数据流 (`pipeline` 包装 `qft` + `ibp` + `ode`)

```
   ┌──────────────────────────────────────────────────────────┐
   │  用户输入                                                 │
   │    * qft::FamilyConfig                                   │
   │    * preferred masters: List<qft::JIntegral>             │
   │    * invariant 数值                                       │
   │    * pipeline::AMFSystemOptions  (mode, ending scheme)   │
   └─────────────────────────┬────────────────────────────────┘
                             │  pipeline::amf_systems_setup
                             ▼
        ┌──────────────────────────────────────────────────┐
        │  pipeline::AMFSystem (root)                      │
        │   ending 节点 → 查 qft::vacuum / 显式 BC          │
        │   非 ending:                                     │
        │     a. qft::analyze_top_sector                   │
        │     b. qft::amf_position(modes) → η_c            │
        │     c. 把 η 注入 propagators                      │
        │     d. ibp::diffeq                               │
        │           ├── ibp::reduce → masters 列表          │
        │           ├── ibp::libp_deriv → ∂Master/∂η       │
        │           └── ibp::kira_run analytic reduction   │
        │           ⇒ dM/dη = A(η, ε) · M                  │
        │     e. qft::find_all_region                      │
        │     f. 对每个非零 region:                          │
        │           qft::region_power      → μ(ε)          │
        │           qft::boundary_integrands → integrand   │
        │           qft::boundary_integrals → per-family   │
        │           递归: 在每个 boundary family 上构建子    │
        │                 AMFSystem                         │
        └─────────────────────────┬────────────────────────┘
                                  │  pipeline::amf_systems_solution(epslist)
                                  ▼
        ┌──────────────────────────────────────────────────┐
        │  对每个 ε in epslist:                              │
        │   自底向上遍历 AMF 树:                              │
        │    ending: 查 qft::vacuum / 显式                   │
        │    非 ending:                                     │
        │      通过记录的 RegionBoundary 表把子解组合成 BC    │
        │      → ode::amflow(de_at_eps, BC)                │
        │      → 本节点的 master values                      │
        └─────────────────────────┬────────────────────────┘
                                  ▼
                  root preferred[] 的 per-ε master 值
                  (Laurent 展开通过 pipeline::FitEps 拿到;
                   用户面向的入口是 pipeline::solve_integrals)
```

### `pipeline` 怎么复用 `ode`

每个非 ending AMF 节点最终都调 `ode::amflow(de_at_eps, BC)` (即用法 1 入口). `pipeline` 只**构造**每个 ε 的 ODE 矩阵和边界 spec, 并**决定**什么时候停止递归.

---

## 每个 domain 负责什么

| Domain | 头文件 (节选) | 职责 |
|--------|---------------------|------------------|
| `numeric` | `options.hpp`, `acb_wrap.hpp`, `log.hpp`, `rational.hpp` | 工作精度, acb/arb 的 RAII 封装, rational 函数 / rational 矩阵类型, `get_poles` |
| `algebra` | `mpoly.hpp`, `mpoly_matrix.hpp`, `context_migration.hpp` | 命名 context 上的多变量 ℤ-多项式与 ℚ rational; 矩阵 det/adjugate/inverse; 跨 context 迁移 |
| `ode` | `blocks.hpp`, `sparse.hpp`, `asy.hpp`, `regular.hpp`, `inf.hpp`, `zero.hpp`, `path.hpp`, `amflow.hpp`, `qqbar.hpp`, `jordan.hpp`, `normalize.hpp`, `acb_rational.hpp` | block 分析, 稀疏 Gauss, 渐近展开代数, regular / 奇点递推, `NormalizeMat`/Jordan/`Calcx00`, contour `RunEta`, 顶层 `amflow()` 驱动 |
| `qft` | `family_config.hpp`, `family_uf.hpp`, `jintegral.hpp`, `topology.hpp`, `amfmode.hpp`, `region.hpp`, `boundary.hpp`, `vacuum.hpp`, `dense_q.hpp` | family 描述, Symanzik U/F, sector 代数, η 注入策略, region decomposition, boundary integrand 构造, vacuum 闭式表 |
| `ibp` | `libp_deriv.hpp`, `kira.hpp`, `reduce.hpp` | IBP 导数, Kira yaml writer / 进程驱动 / 输出 parser, `reduce` 与 `diffeq` (Kira-only 入口) |
| `pipeline` | `amfsystem.hpp`, `factorize.hpp`, `solve_integrals.hpp` | `AMFSystem` 递归树, ending scheme (Tradition/Cutkosky/SingleMass), `factorize_family`, `BlackBoxAMFlow{,Single}`, `FitEps`, `SolveIntegrals` |
| `api` | `run_json.hpp` | JSON dispatcher (mode = `amflow` / `solve_integrals` / `black_box_amflow`); 直接接受 `nlohmann::json` |
| `cli` | (无) | 轻薄的 `main.cpp`, 把 argv 转给 `api::run_json` |

---

## 从哪开始读代码

| 目标 | 入口 |
|---|---|
| 用法 1 pipeline (ode domain) | `src/ode/amflow.cpp::amflow()`, `src/ode/inf.cpp::calc_inf()`, `src/ode/zero.cpp::calcx00()` (与上游 `diffeq_solver/DESolver.m:856-917` 来自 <https://gitlab.com/multiloop-pku/amflow> 对照阅读) |
| 用法 2 pipeline (pipeline domain) | `include/amflow/pipeline/amfsystem.hpp` 头文件注释, 然后 `src/pipeline/amfsystem.cpp::AMFSystem::setup()`/`solve()`, `src/pipeline/factorize.cpp::factorize_family()`, `src/ibp/reduce.cpp::diffeq()` |
| 数据类型 | `numeric/options.hpp`, `numeric/rational.hpp`, `ode/asy.hpp`, `algebra/mpoly.hpp`, `qft/family_config.hpp`, `qft/jintegral.hpp` |
| 从宿主应用调用 | `src/cli/main.cpp` (一屏的 JSON 包装); `tests/test_*.cpp` (每个 domain 自己的测试文件); `tests/test_cli_integration.cpp` (子进程 JSON 示例) |
