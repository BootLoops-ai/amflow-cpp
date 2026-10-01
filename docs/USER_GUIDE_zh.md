# 用户指南

> 英文版本: [USER_GUIDE.md](USER_GUIDE.md)

本指南端到端地讲解如何使用 **AMFlow.cpp**:

1. 它算什么, 什么时候用
2. 怎么选 mode
3. 30 秒可跑的无 IBP 例子 (不需要 Kira)
4. 完整的 Feynman 积分例子 (需要 Kira + Fermat)
5. 怎么读输出
6. 接下来读什么

算法本身请参考上游论文:
[Liu & Ma, *Comput. Phys. Commun.* **283** (2023) 108565](https://doi.org/10.1016/j.cpc.2022.108565).

---

## 1. AMFlow.cpp 算的是什么?

给定一个 Feynman 积分 `I(ε)` (或一组 master integrals), 它通过 Liu & Ma 提出的 **auxiliary mass flow** 算法, 返回 ε 的 Laurent 级数

$$
I(\varepsilon) = \sum_{k \ge k_0} c_k \, \varepsilon^k
$$

到用户指定的精度. 内部流程:

1. 把积分约化到一组封闭的 master integrals (本项目用 Kira 作为后端).
2. 在 auxiliary mass `η` 上构造微分方程 `dM/dη = A(η, ε) · M`.
3. 从 η = ∞ (边界条件平凡: vacuum bubbles, single-mass tadpoles) 求解到 η = 0 (物理极限).
4. 在渐近展开的每个 region 递归构造并求解子系统.
5. 在每个 ε 上拿到数值后, 拟合成 Laurent 展开.

使用工具不需要懂 step 2-5, 但想了解全貌见 [`docs/ARCHITECTURE.md`](ARCHITECTURE.md).

---

## 2. 三种使用 mode

`amflow_cli` 接受单一 JSON 文档. 顶层的 `"mode"` 键选择三种工作流之一:

| Mode | 什么时候用 | 你提供什么 | 拿回什么 |
|---|---|---|---|
| `amflow`           | 你已经有了 ODE 系数矩阵 `A(η)` (rational) 和 η = ∞ 处的渐近边界条件, 想隔离测试 ODE 引擎, 或输入来自其它工具. | 一个方形 rational 矩阵 + 每行一组 `(μ, value)` 边界对 | 向量 `I(η = 0)` |
| `solve_integrals`  | 你有 Feynman family + 目标积分 + ε 精度目标, 要 Laurent 展开. | family 配置 (loops, legs, propagators, kinematics) + 目标积分 + `goal_digits` + `eps_order` | 每个目标的 `Σ c_k ε^k` |
| `black_box_amflow` | 同上, 但你自己指定 ε 网格 (不做 Laurent 拟合). | family 配置 + 目标 + `eps_samples` | 每个目标在每个 ε 上的数值 |

大多数用户应该用 **`solve_integrals`** — 它做完整的工作. `black_box_amflow` 用于在单 ε 上与 Mathematica AMFlow 做 parity 对比. `amflow` (raw mode) 只在你已经手上有 ODE 矩阵时用.

字段级完整参考见 [`docs/JSON_SCHEMA.md`](JSON_SCHEMA.md).

> **不在 scope 内** — **复数值数值 kinematics** 不支持. `amf_options.blackbox.numeric_values` 只允许实数; dispatcher 会用清晰错误信息拒绝 `{"re":..,"im":..}` 对象形式. 完整的有意非目标见 [`docs/AUDIT_MMA_PARITY.md`](AUDIT_MMA_PARITY.md); 面向用户的非目标见 [`docs/FAQ.md`](FAQ.md) "What's *not* implemented?".

---

## 3. 30 秒例子 — 无 IBP, 无 Kira

[`examples/power_law.json`](../examples/power_law.json) 编码了一个有解析解的 1×1 ODE:

```math
\frac{dI(\eta)}{d\eta} = \frac{1}{\eta} I(\eta), \qquad
I(\eta) \sim \eta^1 \;\; (\eta \to \infty).
```

精确解是 `I(η) = η`, 所以 `I(0) = 0`. JSON 文件:

```jsonc
{
  "options": {
    "x_order":       30,
    "extra_x_order": 10,
    "silent_mode":   true
  },
  "matrix":     [[ { "num": ["1"], "den": ["0", "1"] } ]],
  "boundaries": [[ { "mu":    { "re": "1", "im": "0" },
                     "value": { "re": "1", "im": "0" } } ]],
  "mode": "amflow"
}
```

每个键的含义:

- `mode: "amflow"` — raw ODE 求解, 不涉及 family / Kira.
- `matrix` — 一个 `1 × 1` 数组. 唯一一项是 rational 函数 `1/η`, 分子多项式写成 `[1]` (= 1), 分母多项式写成 `[0, 1]` (= η). 数组 index 表示 η^index 的系数, 所以 `[0, 1]` 代表 `0 + 1·η`.
- `boundaries` — 每行一个边界列表. 外层数组一项 (唯一一行); 内层列表一项: `value · η^μ = 1 · η^1`, 表示 η = ∞ 处的 leading 行为.
- `options.x_order = 30` — 级数展开保留 30 项.
- `options.silent_mode = true` — 抑制 progress 输出.

构建项目 (假定你已 `cmake -S . -B build`):

```bash
cmake --build build -j32
```

运行:

```bash
./build/src/cli/amflow_cli examples/power_law.json
```

> 提示 — 跑过 `sudo cmake --install build` 后 `amflow_cli` 就进了 `$PATH` (默认 prefix `/usr/local/`), 本指南后续命令都可以省略 `./build/src/cli/` 前缀, 直接写 `amflow_cli ...`. 完整安装流程见 [`FAQ_zh.md`](FAQ_zh.md) "怎么把 library 和 `amflow_cli` 二进制装到系统?".

输出 (节选):

```jsonc
{
  "mode": "amflow",
  "result": [
    { "re": "0", "im": "0" }
  ],
  "options": { "x_order": 30, ... }
}
```

`result` 数组一项 — 数值 `I(η = 0) = 0` — 与解析解吻合.

[`examples/constant.json`](../examples/constant.json) 更简单: `dI/dη = 0`, 边界 `I ~ 42·η^0` → 输出 `42`.

至此你已经跑通一个端到端的 C++ ODE 引擎调用.

---

## 4. 完整例子 — `solve_integrals` 配合 IBP

`solve_integrals` 是主打工作流, 需要:

1. 一份 Feynman family 描述.
2. 一组目标 master integrals.
3. 精度目标 (`goal_digits`, `eps_order`).
4. 每个 kinematic invariant 的数值.

### 4.1 前置依赖 — 安装 Kira 和 Fermat

IBP 后端用 **Kira**. 没有它的话, `solve_integrals` 与 `black_box_amflow` 都不能跑.

- Kira: <https://kira.hepforge.org/> — 安装到 `/usr/local/bin/kira` (或通过 `amf_options.blackbox.kira_executable` 重写路径).
- Fermat: <https://home.bway.net/lewis/> — 安装到 `/usr/share/Ferl7/fer64` (或通过 `amf_options.blackbox.fermat_executable` 重写).

两者都是**运行时**依赖, 不是构建依赖. 构建本身不链接它们; 项目在需要 IBP 时通过子进程 shell 调用.

### 4.2 示例: 1-loop massive bubble

bubble 积分

```math
\mathcal{B}(s) = \int \frac{\mathrm{d}^d \ell}{(2\pi)^d} \;
                 \frac{1}{(\ell^2 - m^2)\, ((\ell - p)^2 - m^2)},
\quad s = p^2
```

是带解析解的最简单非平凡 Feynman 积分例子.

[`examples/bubble_solve_integrals.json`](../examples/bubble_solve_integrals.json):

```jsonc
{
  "mode": "solve_integrals",
  "options": {
    "chop_pre":        20,
    "silent_mode":     true
  },
  "family": {
    "name":  "bubblefam",
    "loops": ["l"],
    "legs":  ["p"],
    "replacement": { "p^2": "s" },
    "propagators": [
      "l^2 - msq",
      "(l - p)^2 - msq"
    ]
  },
  "integrals":   [ { "indices": [1, 1] } ],
  "goal_digits": 25,
  "eps_order":   2,
  "work_dir":    "/tmp/desolver_example_bubble",
  "amf_options": {
    "blackbox": {
      "numeric_values": { "s": "9", "msq": "1" }
    }
  }
}
```

逐键解释:

- **`family`** — 拓扑描述.
  - `loops` / `legs` — 圈动量 `l` 和外动量 `p`.
  - `replacement` — kinematic invariants. 这里 `p² = s`.
  - `propagators` — 两个标量分母, 用 `l`, `p` 和符号化的 mass 参数 `msq` 写出.
- **`integrals`** — 目标积分. 每一项有长度等于 `propagators` 的 `indices`. `[1, 1]` 是每个分母都 1 次的 master integral.
- **`goal_digits: 25`** — 每个 ε 系数拟合到 ~25 位精度.
- **`eps_order: 2`** — 对这个 1-loop bubble (L=1) 来说拟到 `ε^(2 − 2L) = ε⁰`, 结果含 ε⁻¹ 和 ε⁰ 两个系数 (形式上的 ε⁻² leading 对此有限 bubble 为零, 被 trim 掉). 通用规则: 想要结果到 ε^K 就传 `eps_order = K + 2L`.
- **`work_dir`** — Kira 中间文件放在这里. 多次运行可复用 (cache 命中靠指纹守住).
- **`amf_options.blackbox.numeric_values`** — 给*所有*在 `replacement` 里出现过的符号和*所有*在 `propagators` 里出现过的符号化 mass 都赋数值. 这里 `s = 9`, `msq = 1`.

### 4.3 运行

```bash
./build/src/cli/amflow_cli examples/bubble_solve_integrals.json
```

预期用时: 现代桌面机几秒. Kira 的 `/tmp/desolver_example_bubble/` workdir 会堆积 yaml 文件, 中间结果, 以及每轮 IBP 一个 `kira_target.m`.

输出 (节选):

```jsonc
{
  "mode": "solve_integrals",
  "result": [
    {
      "integral": { "family": "bubblefam", "indices": [1, 1] },
      "leading_order": 0,
      "coefficients": [
        { "order": 0, "value": { "re": "...", "im": "..." } },
        { "order": 1, "value": { "re": "...", "im": "..." } },
        { "order": 2, "value": { "re": "...", "im": "..." } }
      ]
    }
  ],
  "options": { ... }
}
```

Laurent 展开:

`B(s, ε) ≈ c_0 + c_1 · ε + c_2 · ε²`

其中 `c_k = result[0].coefficients[k].value`.

### 4.4 在非 4 维时空里计算

默认宿主时空维度 `D = 4 - 2ε` (`ε → 0` 时恢复物理四维). 要在别的维度方案里工作 — 比如 `D = 7/3 - 2ε` — 设 `options.d0 = "<rational>"`:

```jsonc
{
  "mode": "solve_integrals",
  "options": {
    "d0": "7/3",          // 任何 rational "p/q" 或整数; 默认 "4"
    "rationalize_pre": 100
  },
  "family":     { ... },
  "integrals":  [ ... ],
  "goal_digits": 30,
  "eps_order":   4,
  "amf_options": { "blackbox": { "numeric_values": { "s": "100", "t": "-1" } } }
}
```

内部发生的事: 用户面向的 ε 网格不变, 但引擎在 `ε + (4 − D₀)/2` 上对 family 求值 (镜像上游 `AMFlow.m:1342`/`1351`), 然后把 Laurent 展开拟合回原始网格 (镜像 `AMFlow.m:1356`). 这样 `D₀` 在 API 层不可见, 而在算法内部干净抵消.

仓库 committed 一个端到端的 oracle benchmark [`tools/bench/box1_d0_7_3_solve_integrals_cpp.json`](../tools/bench/box1_d0_7_3_solve_integrals_cpp.json), 在 `D₀ = 7/3, s = 100, t = -1` 上解一个 1-loop box, 与 Mathematica AMFlow 在 order 0..2 上吻合到 ~30 位有效数字.

### 4.5 改成你自己的积分

要解别的积分, 改这几处:

1. **`family.propagators`** — 你的拓扑.
2. **`family.replacement`** — 你的 kinematics (所有不是 loop 动量也不是 `i0` 的量).
3. **`integrals`** — 你的目标, 一个 master 一项.
4. **`amf_options.blackbox.numeric_values`** — 给 `replacement` 里的每个 invariant 和在 `propagators` 里出现但不在 `replacement` 里的每个裸 mass 都赋具体数值.
5. **`goal_digits`** / **`eps_order`** — 按需设. 代价大约随 `goal_digits` 线性增长, 随内部 mass + invariant 数量指数增长.

更大的 family ([`examples/box1_solve_integrals.json`](../examples/box1_solve_integrals.json) 是一个 1-loop box, 带两个 invariants `s, t`) 套同样的模板.

---

## 5. 怎么读输出

每次运行都产生一个顶层对象:

```jsonc
{
  "mode":    "<echo 的 mode 名>",
  "result":  <按 mode 不同的结构>,
  "options": <input.options 应用后的有效 options>
}
```

每个复数编码为 `{ "re": "<十进制字符串>", "im": "<十进制字符串>" }`. 十进制字符串保留到 working precision (`options.working_pre` 位) 的精确精度. 用 Python 解析:

```python
import json, mpmath
mpmath.mp.dps = 100        # 与 working_pre 对齐

out = json.load(open("/tmp/bubble_out.json"))
for entry in out["result"]:
    integral = entry["integral"]
    coefficients = entry["coefficients"]
    print(f"\n{integral['family']}{integral['indices']}:")
    for c in coefficients:
        order = c["order"]
        val = mpmath.mpc(c["value"]["re"], c["value"]["im"])
        print(f"  ε^{order}: {val}")
```

每个输出字段的细节见 [`docs/JSON_SCHEMA.md`](JSON_SCHEMA.md).

---

## 6. 接下来读什么

| 你想... | 读 |
|---|---|
| 查具体的输入 / 输出字段 | [`docs/JSON_SCHEMA.md`](JSON_SCHEMA.md) |
| 了解算法和库的设计 | [`docs/ARCHITECTURE.md`](ARCHITECTURE.md) |
| 看 parity 状态和已验证的 family | [`docs/AUDIT_MMA_PARITY.md`](AUDIT_MMA_PARITY.md) |
| 查某个 Mathematica 函数对应的 C++ 实现 | [`docs/REFERENCE_MAP.md`](REFERENCE_MAP.md) |
| 解决 build / 运行时问题 | [`docs/FAQ.md`](FAQ.md) |
| 重新生成 Mathematica 参考数据 | [`reference/README.md`](../reference/README.md) 与 [`tools/bench/README.md`](../tools/bench/README.md) |
| 参与贡献 | [`CONTRIBUTING.md`](../CONTRIBUTING.md) 与 [`docs/CONTRIBUTING.md`](CONTRIBUTING.md) |
| 引用本软件 | [`CITATION.cff`](../CITATION.cff) |

本指南和 FAQ 都没覆盖的问题, 请到 <https://github.com/chang18/amflow-cpp/issues> 开 issue, 或联系维护者 <3250800970@qq.com>.
