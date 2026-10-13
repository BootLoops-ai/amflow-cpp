> **BootLoops fork.** This repository is the BootLoops project's patched fork of
> [AMFlow.cpp](https://github.com/chang18/amflow-cpp) (MIT; Zenodo
> [10.5281/zenodo.20087172](https://doi.org/10.5281/zenodo.20087172)), the C++17
> port, written and maintained by the AMFlow.cpp contributors (maintainer: GitHub
> [\@chang18](https://github.com/chang18)), of the Mathematica package AMFlow by
> Xiao Liu and Yan-Qing Ma: a snapshot of upstream `main` at
> commit `006a275` (2026-05-22; v1.1.0 plus the [Unreleased] CHANGELOG entries)
> with the BootLoops project's modifications applied. This repository does not
> carry upstream's commit history; upstream's history and
> releases live at the link above. We are grateful to the AMFlow.cpp maintainer:
> the 228-case Mathematica-oracle benchmark and the GoogleTest suite that ship
> with upstream are what made every change in this fork checkable. The
> auxiliary-mass-flow method is due to Xiao Liu, Yan-Qing Ma and their
> collaborators Chen-Yu Wang, Zhi-Feng Liu, Wei Tao, Peng Zhang and Rui-Jun
> Huang, and the original Mathematica package
> [AMFlow](https://gitlab.com/multiloop-pku/amflow) (MIT) is by Xiao Liu and
> Yan-Qing Ma (AMFlow 2.0 with Rui-Jun Huang). The BootLoops additions were
> written with the Mathematica source of AMFlow 2.0 open. The derivative
> routines of AMFlow's Kira interface, which AMFlow.cpp ports in
> `src/ibp/libp_deriv.cpp`, are from LiteIBP by Tiziano Peraro
> ([finiteflow-mathtools](https://github.com/peraro/finiteflow-mathtools));
> that repository carries no license file (see [NOTICE](NOTICE) and the License
> section of [PATCHES.md](PATCHES.md#license)).
> Please cite the AMFlow authors' papers, listed
> with the other scientific credits in
> [PATCHES.md](PATCHES.md#scientific-credit-and-thanks), together with upstream's
> `CITATION.cff`, whenever results obtained with this fork are published. The IBP
> reductions this engine launches are done by Kira (created by Philipp
> Maierhöfer, Johann Usovitsch and Peter Uwer; developed today by Fabian Lange,
> Johann Usovitsch and Zihao Wu, with Jonas Klappert on Kira 2.0) with FireFly
> (Jonas Klappert, Sven Yannick Klein, Fabian Lange) and Fermat (Robert H.
> Lewis); see the same section for references. What changed relative to that
> base and why: [PATCHES.md](PATCHES.md); build and run notes for the additions:
> [BUILD.md](BUILD.md); our validation harness for the engine:
> [bootloops-wrappers/](bootloops-wrappers/). Upstream's `LICENSE` and
> `CITATION.cff` are unchanged. `NOTICE` is upstream's text, unchanged, followed
> by sections that this fork appended and marked as its own. Our modifications are
> released under the same MIT license, except as stated under "LiteIBP" in `NOTICE`,
> Copyright (c) 2026 Anthropic, PBC. Created by
> Matthew D. Schwartz; code for the modifications written by Claude (Anthropic) under
> his supervision. This is not an officially supported Anthropic product; it is
> maintained by Matthew D. Schwartz (https://www.bootloops.ai).
> The provenance section, `AGENTS_zh.md`, `CONTRIBUTING_zh.md`, the first
> 35 lines of `NOTICE`,
> `CITATION.cff`, the CI badge and the maintainer contact below are upstream
> AMFlow.cpp's and describe upstream's own development; they do not describe the
> BootLoops modifications. Questions or bugs concerning the BootLoops
> modifications: open an issue on this repository; the issue links and
> maintainer address below apply to unmodified AMFlow.cpp only. The upstream
> Chinese README follows unchanged.

# AMFlow.cpp

[![CI](https://github.com/chang18/amflow-cpp/actions/workflows/ci.yml/badge.svg)](https://github.com/chang18/amflow-cpp/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Version](https://img.shields.io/badge/version-1.1.0-blue.svg)](CHANGELOG.md)
[![C++17](https://img.shields.io/badge/C%2B%2B-17-blue.svg)](https://en.cppreference.com/w/cpp/17)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20087172.svg)](https://doi.org/10.5281/zenodo.20087172)

> 英文版本: [README.md](README.md)

多圈 Feynman 积分的 auxiliary mass flow (辅助质量流) 算法的 C++17 重实现.

## 科学引用归属

本项目是 Xiao Liu 与 Yan-Qing Ma 提出的 Mathematica 软件包 **AMFlow** 中算法的重实现:

- 上游 Mathematica 软件包: <https://gitlab.com/multiloop-pku/amflow>
- 算法论文: X. Liu and Y.-Q. Ma, *AMFlow: A Mathematica package
  for Feynman integrals computation via auxiliary mass flow*, Comput.
  Phys. Commun. **283** (2023) 108565,
  [doi:10.1016/j.cpc.2022.108565](https://doi.org/10.1016/j.cpc.2022.108565).

本仓库**仅**贡献 C++17 实现以及对 Mathematica 参考实现做数值一致性验证的测试 / benchmark 框架; **算法的原创归属完全属于 AMFlow 作者**, 本项目不主张任何算法上的新意.

若在已发表研究中使用本软件, 请**同时**引用上游论文 (算法) 和本仓库 (实现). 本仓库已存档于 Zenodo, DOI 为 [10.5281/zenodo.20087172](https://doi.org/10.5281/zenodo.20087172) (concept DOI, 始终指向最新发布版). 机读元数据见 [`CITATION.cff`](CITATION.cff); 完整归属信息见 [`NOTICE`](NOTICE).

> **注意** — 上游 Mathematica 源代码**未**随本仓库分发. 项目文档中引用上游文件 (`AMFlow.m`, `diffeq_solver/DESolver.m`, `ibp_interface/Kira/interface.m`) 时直接给出符号名与行号; 要跟随这些引用并重新生成 Mathematica 参考数据, 请将上游 clone 到 `reference/amflow-master/` (该路径已 gitignored). 详见 [`reference/README.md`](reference/README.md).

## 项目背景

本 C++17 实现主要由 **AI 编程 agent** 在维护者指导下开发, 每一项行为变更都以"与上游 Mathematica 参考实现的数值一致性"作为闸门. `tools/bench/` 下的 228 个 oracle benchmark (1L=52, 2L=53, 3L=64, 4L=59) 与 549-case GoogleTest 测试套是 AI 主导实现必须遵守的契约, 任何 landing 的改动都要满足.

维护者 / 联系方式: **3250800970@qq.com** (能用 GitHub issue 时优先开 issue; 邮件用于不适合公开 issue 的咨询).

## 代码结构

按 domain 组织的架构:

```
amflow::numeric -> amflow::algebra -> {amflow::ode,
                                       amflow::qft -> {amflow::ibp,
                                                       amflow::pipeline -> amflow::api -> amflow::cli}}
```

公开头文件位于 `include/amflow/<domain>/`, 实现位于 `src/<domain>/`.

## 工作规则

- 若一个工作流在上游 `AMFlow.m` 上能跑通, 而对应分支已经在本项目 port 完成, 那么 C++ 端不一致就是一致性 bug.
- 不要发明 Mathematica 一致性之外的额外正确性判据.
- 单点采样一致性 (`eps = 1/1000`) 是主要的数值契约, 见 `tools/bench/`.
- 所有 Mathematica 算出的 benchmark 数值都作为回归基准 committed 在仓库里 (原始 cache 可重新生成, 见 `tools/bench/README.md`).

## 状态

| 模块 | 状态 |
|---|---|
| 核心 ODE solver (上游 `DESolver.m` 的 port) | 已实现; 228/228 oracle case 与 MMA 参考一致, 相对精度 ~1e-30 (pentabox 2L 角落 case 因积分本身的取消尺度限制, 一致到 ~1e-10) |
| AMFlow + Kira pipeline (上游 `AMFlow.m` 的核心算法) | 在已覆盖的工作流上完成实现 |
| 顶层 entry (`amflow`, `black_box_amflow`, `solve_integrals`) | 已实现; 通过 `amflow_cli` 的 JSON mode 暴露 |
| 行级 MMA 一致性审计 | 已 close; 每个公开符号都至少被一个 oracle 覆盖. 详见 [`docs/AUDIT_MMA_PARITY.md`](docs/AUDIT_MMA_PARITY.md) |
| 相对 MMA 的 wallclock 性能 | 延迟到专用机器运行后再做正式 benchmark, 详见 [`docs/ROADMAP.md`](docs/ROADMAP.md) §"性能 benchmark" |
| `SolveIntegralsGaugeLink`, HQET / SCET / Wilson lines | 不在 scope 内 (明确决定) |
| 复数值数值 kinematics | 不在 scope 内 (不支持); JSON dispatcher 会拒绝 `{"re":..,"im":..}` 形式. 见 [`docs/FAQ.md`](docs/FAQ.md) "What's *not* implemented?" |

行级上游一致性审计 (完整 divergence inventory) 见 [`docs/AUDIT_MMA_PARITY.md`](docs/AUDIT_MMA_PARITY.md); oracle bench 三元组见 [`tools/bench/`](tools/bench/).

## 快速开始

Ubuntu 上, 默认构建所需的系统包:

```bash
sudo apt-get update
sudo apt-get install -y build-essential cmake pkg-config libflint-dev \
  libgtest-dev nlohmann-json3-dev
```

```bash
cmake -S . -B build -DAMFLOW_BUILD_DRIVER=ON
cmake --build build -j32
ctest --test-dir build --output-on-failure -j 4

# 安装到 /usr/local/{bin,lib,include} (CMake 默认 prefix). 这一步
# 把 amflow_cli 放进 PATH, 同时把 library + headers 暴露给下游
# find_package(). 系统级安装要 sudo; 若想装到用户目录, 改用
# --prefix "$HOME/.local".
sudo cmake --install build

# 纯 ODE 示例 (不需要 IBP). 不装的话直接用 build/ 下的二进制:
# ./build/src/cli/amflow_cli examples/power_law.json
amflow_cli examples/power_law.json

# 完整 AMFlow 示例 (运行时需要 Kira + Fermat)
amflow_cli examples/box1_black_box_amflow.json
amflow_cli examples/bubble_solve_integrals.json
amflow_cli examples/box1_solve_integrals.json
```

作为 CMake 依赖使用 (安装后):

```cmake
find_package(AMFlowCpp 1.1 REQUIRED)
target_link_libraries(my_target PRIVATE AMFlowCpp::amflow)
```

## Mathematica 参考数据生成

仓库提供 Wolfram driver, 用于重新生成 committed 在 `tests/data/math_ref/` 下的参考数据以及 `tools/bench/` 下的采样 benchmark 数据. 这些 driver 需要在 `reference/amflow-master/` 路径下有一份本地的上游 Mathematica AMFlow clone (该路径已 gitignored, 两种合法布局见 [`reference/README.md`](reference/README.md)).

```bash
# 方式一: 直接 clone 上游到 reference/amflow-master/
git clone https://gitlab.com/multiloop-pku/amflow.git reference/amflow-master

# 方式二: 用 symlink 指向已有的 clone
# ln -s /path/to/your/clone reference/amflow-master

cd tools/math_ref
AMF_REF_CFG=$PWD/cfg_bubble_1L.wl math -script run_amflow_kira.wl
```

`eps = 1/1000` 处的当前一致性 benchmark 见 [`tools/bench/`](tools/bench/) 与 [`docs/AUDIT_MMA_PARITY.md`](docs/AUDIT_MMA_PARITY.md).

## 文档导览

**面向使用者** — 按以下顺序阅读:

1. [`docs/USER_GUIDE.md`](docs/USER_GUIDE.md) — 一份"无 IBP 例子"加一份"完整 Feynman 积分例子"的逐步教程.
2. [`docs/JSON_SCHEMA.md`](docs/JSON_SCHEMA.md) — 三种 CLI mode 全部输入 / 输出字段的完整定义.
3. [`docs/FAQ.md`](docs/FAQ.md) — 常见 build / 运行时 / 一致性问题.

**面向贡献者** — 额外补充:

| 文档 | 涵盖内容 |
|---|---|
| [`CONTRIBUTING.md`](CONTRIBUTING.md) | 贡献规则简介 + 维护者联系方式 |
| [`docs/CONTRIBUTING.md`](docs/CONTRIBUTING.md) | 详细开发流程: 测试, benchmark, 参考数据, 调试, env-var trace 开关 |
| [`AGENTS.md`](AGENTS.md) | AI 编程 agent 的 onboarding (AI 主导开发的契约) |
| [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md) | Domain DAG 与端到端数据流 |
| [`docs/INVARIANTS.md`](docs/INVARIANTS.md) | 项目级必须维持的不变式 (精度, RAII, ODE 边界陷阱等) |
| [`docs/REFERENCE_MAP.md`](docs/REFERENCE_MAP.md) | 上游 Mathematica 符号 → C++ 符号映射 |
| [`docs/AUDIT_MMA_PARITY.md`](docs/AUDIT_MMA_PARITY.md) | 与上游 MMA AMFlow 的活跃一致性状态 (逐行 divergence 目录 + 闭环日志) |
| [`docs/ROADMAP.md`](docs/ROADMAP.md) | 活跃开发计划 (前向 — oracle 扩展) |
| [`reference/README.md`](reference/README.md) | 如何 clone 上游 MMA AMFlow 到本地以重生成参考数据 |
| [`tools/bench/README.md`](tools/bench/README.md) | 采样一致性 benchmark 的约定与再生方式 |

## 依赖

| 组件 | 用途 |
|---|---|
| C++17 编译器 | 核心构建 |
| CMake >= 3.16 | 构建系统 |
| FLINT / Arb | rational, 多项式, 任意精度算术 |
| GoogleTest | 测试套 |
| nlohmann/json | Driver 的 JSON I/O |
| Kira + Fermat | IBP 约化后端 (运行时依赖; 测试在缺失时会优雅 skip) |
| Mathematica / Wolfram Engine + 上游 AMFlow clone | 可选; 仅用于重新生成参考数据 |

## 参与贡献

GitHub 入口见 [`CONTRIBUTING.md`](CONTRIBUTING.md), 完整开发流程见 [`docs/CONTRIBUTING.md`](docs/CONTRIBUTING.md). 行为性变更必须有"与上游 AMFlow 一致"的证据以及聚焦的 C++ 回归测试支撑.

如有问题, bug 反馈或合作咨询, 请到 <https://github.com/chang18/amflow-cpp/issues> 开 issue, 或联系维护者 <3250800970@qq.com>.

## 许可

MIT — 见 [`LICENSE`](LICENSE) 与 [`NOTICE`](NOTICE). 上游 Mathematica AMFlow 同样是 MIT 许可; 本项目不打包任何上游源码.
