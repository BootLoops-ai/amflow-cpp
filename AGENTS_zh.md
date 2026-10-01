# AGENTS.md — AI 编程 agent 的 onboarding

> 英文版本: [AGENTS.md](AGENTS.md)

本文件是任何 AI 编程 agent 接手本仓库工作时的规范 onboarding 文档. 它捕获项目级规则和新协作者动代码或文档之前必须内化的约定. 更深的文档在 [`docs/`](docs/) 与 [`docs/AUDIT_MMA_PARITY.md`](docs/AUDIT_MMA_PARITY.md); 本文件指向那里, 并列出承重规则.

---

## 项目一览

- 本仓库是 **AMFlow.cpp** — Mathematica 包 AMFlow (<https://gitlab.com/multiloop-pku/amflow>, 多圈 Feynman 积分的 auxiliary-mass-flow 框架) 的 domain-oriented C++17 重实现. 上游 Mathematica 源**未**打包 — 见 [`reference/README.md`](reference/README.md).
- 布局: 按 domain 划分的子 namespace `amflow::<numeric|algebra|ode|qft|ibp|pipeline|api|cli>::`, 源代码在 `src/<domain>/` 与 `include/amflow/<domain>/`.
- Mathematica AMFlow 是唯一的行为规范. 已 port 分支上 C++ 不一致就是 parity bug.
- **本 codebase 主要由 AI 编程 agent 开发** (本文档是 onboarding 契约). 维护者 (<3250800970@qq.com>) 指导工作但不做行级 review; parity 契约 — `tools/bench/` 下的 228 oracle benchmark 与 548 case GoogleTest 测试套 — 是每个改动必须过的闸.

## 从哪开始读

| 目标 | 入口 |
|---|---|
| 项目介绍, 构建指引 | [`README.md`](README.md) |
| 当前 parity 状态, 已验证 family, benchmark 清单 | [`docs/AUDIT_MMA_PARITY.md`](docs/AUDIT_MMA_PARITY.md) |
| Domain DAG 与数据流总览 | [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md) |
| 跨改动必须保持的规则 | [`docs/INVARIANTS.md`](docs/INVARIANTS.md) |
| MMA 符号 → C++ 符号映射 | [`docs/REFERENCE_MAP.md`](docs/REFERENCE_MAP.md) |
| 加测试 / bench / 参考数据; 调试 | [`docs/CONTRIBUTING.md`](docs/CONTRIBUTING.md) |

---

## 工作规则

当前所有工作都是 functional / parity 工作 (或极少数情况下, 保留行为的代码健康工作).

### Functional / parity 工作 — 主要模式

触发条件: 关闭 parity gap, port `AMFlow.m` 里缺的分支, 修错的数值结果, 加之前没 port 的功能.

**规则:**

- **真理来源是上游 Mathematica `*.m` 文件** (<https://gitlab.com/multiloop-pku/amflow>). 设计修复前读相关 Mathematica 段落; 在 commit / PR 描述里引用上游 文件 + 行号范围.
- **严格镜像 — 不走捷径.** 不要替换成 closed-form 特例, 不要关功能, 不要缩窄 scope, 不要在 divergence 上糊纸. 任何非镜像修复本身就是另一形态的第二个 divergence.
- **"现有测试通过" ≠ "与 MMA 对齐."** 绿 ctest 只证明之前走过的路径还能跑. 验证 closure 要对照 AMFlow.m 行号范围, 不要单看下游测试状态.
- 当替代方案是非镜像捷径时, 打破长存的 port 级不变式 (loop vs leg 分类, FamilyConfig 形状等) *已经预先授权*.

### 代码健康 / 保留行为工作 — 次要模式

触发条件: **不得**改变数值输出的结构性改进 (重命名, 文件拆分, 死代码移除, 源注释里的文档修复).

**规则:**

- **真理来源是测试套 + committed 采样 bench.** 编辑过程中*不*要求 AMFlow.m 源行对应, 但测试锚定的行为必须 byte-for-byte 保留.
- **每个 commit 都要遵守的契约:**
  - `ctest --test-dir build --output-on-failure -j 4` 保持绿;
  - [`docs/AUDIT_MMA_PARITY.md`](docs/AUDIT_MMA_PARITY.md) 列出的每个 committed 采样 bench 三元组重跑时仍数值上绿;
  - `include/amflow/<domain>/` 中的公开 API 保持签名, 除非改动是有意, 有文档的破坏.
- **一个 commit 一个关注点.** 不要打包 "重命名 + 拆分 + 重写" — diff 必须保持可 bisect.
- 如果结构性改动会改变测试期望或 bench 输出, **停**. 那是 functional 工作, 不是代码健康工作.

### 区分两种模式

| 问题 | Functional | 代码健康 |
|---|---|---|
| 这个改动会改变数值输出吗? | 会 | 不会 |
| 什么钉住行为? | AMFlow.m 行号对应 | 测试 / bench 套件 |
| 分界测试 | 引用一个 AMFlow.m 行号范围 | 所有测试 + bench 仍通过 |

不知道你处于哪个模式时, 先问用户再推进.

---

## 操作约定

### 构建 & 测试

```bash
cmake -S . -B build -DAMFLOW_BUILD_DRIVER=ON
cmake --build build -j32
ctest --test-dir build --output-on-failure -j 4
```

- **用 `-j32` 并行构建**, 不是 `$(nproc)`. 宿主与别的 workload 共享核; 占满 CPU 会饿死它们.
- **永远不要只重构建 `amflow_cli` 然后信任之后跑的 `ctest`.** 跑 `cmake --build build --target amflow_cli` 只重建 CLI 但留下旧的 `amflow_tests`; 然后 `ctest` 用**旧**测试二进制对**新**库, 行为已经改了却静默报绿. 要么不用 `--target` 跑 (默认构建所有, 含 `amflow_tests`), 要么在跑 `ctest` 之前用 `--target amflow_tests`. 拿不准时, 跑裸 `cmake --build build -j32`.
- bench 命令: 见 [`tools/bench/README.md`](tools/bench/README.md). 每个 committed 采样 bench 有一个三元组: `*_mma.wl` (重新生成) + `*_cpp.json` (driver 输入) + `*_mma_reference.json` (对比的参考值).

### 对话语言

- 所有面向用户的对话回复必须用**中文** (本项目所有对话回复一律使用中文).
- 代码, 注释, commit 信息, 文件内容 (含本文件和所有 `docs/*.md` 英文版) 保持**英文**.

### 文档同步原则

- 改项目状态时, 同步**整个** `docs/` 树.
- **删历史问题要狠.** 已解决的 bug, 放弃的方案, "X 坏的时候" 这种叙事都属于 git 历史, 不属于 living docs.
- 深度诊断 recipe (特定 debug-env-var 组合, 一次性失败模式的坑位列表) 属于 agent 私人记忆或 commit 信息, 不属于项目文档.

### Bench WL 操作坑

写一份顶层 `tools/bench/*_solve_integrals_mma.wl` inline 镜像 AMFlow.m 的 `SolveIntegrals` 主体 (这样 committed cache 目录支撑 `BlackBoxAMFlow`) 时, 三个 Private-context 符号必须 qualify — 裸引用静默泄漏成未绑定 `Global`` 符号, 产生卡死 / 错的 leading order / 折叠系数:

| 裸引用 (`Global`) | 必须用 | AMFlow 在哪绑它 |
|---|---|---|
| `$D0` | `` `AMFlow`Private`$D0` `` | `AMFlow.m:262` (在 `Begin["`Private`"]` 内) |
| `Loop` | `AMFlowInfo["Loop"]` | `AMFlow.m:192` |
| `$Eps` | `` `AMFlow`Private`$Eps` `` (resolve 到 `Symbol["Global`eps"]`) | `AMFlow.m:181` |

规范修复模式见 [`tools/bench/box1_d0_7_3_solve_integrals_mma.wl`](tools/bench/box1_d0_7_3_solve_integrals_mma.wl).

### 谨慎行事

- 默认对任何**难撤销**操作 (force-push, `git reset --hard`, 强制 amend 已发布 commit, 删表, 杀进程等) **公开询问**再做.
- 不要跳过 git hook (`--no-verify`, `--no-gpg-sign`), 除非用户明确授权.
- 用户没明确要求时不要 commit. 给 diff 等批准是安全默认.
- 调研根本原因; 不要把破坏性操作当作 "让障碍消失" 的捷径.

---

## 你**不应该**做的事

- **不要在 functional 工作里用捷径替代**真正与 MMA 对齐的修复.
- **不要在代码健康 commit 里改测试期望或 bench 输出.** 觉得必须改时, 你做的就是 functional 工作.
- **不要不带 `-j32` 跑 `cmake --build`**, 除非有理由覆盖 host-load 约定.
- **不要把多个关注点打包**进一个 commit.
- **没有明确用户授权不要绕过 git hook.**
- **不要自主 commit.** 给 diff, 等用户批准.
- **不要为假设的未来需求加投机抽象.** 三行类似代码胜过一个过早的抽象.
- **不要给没改的代码加 docstring/注释.** 只在逻辑不显然处加注释.
- **不要猜近期项目状态** — `git log` 与 `docs/AUDIT_MMA_PARITY.md` 是权威; 断言前先查.
