# `amflow_cli` JSON Schema 参考

> 英文版本: [JSON_SCHEMA.md](JSON_SCHEMA.md)

`amflow_cli` 接受单一 JSON 文档作为输入, 输出单一 JSON 文档. 本页给出三种顶层 mode 中所有输入 / 输出字段的完整定义.

C++ 入口是 `amflow::api::run_json(const json&)` ([`include/amflow/api/run_json.hpp`](../include/amflow/api/run_json.hpp)); 本文档与 [`src/api/run_json.cpp`](../src/api/run_json.cpp) 的实现必须保持一致.

> **约定** — 标 **required** 的键必须存在. 标 **optional** 的键有 documented 的默认值. 字段接受多个别名时, 列出的第一个是首选.

---

## 顶层结构

```jsonc
{
  "mode":        "amflow" | "solve_integrals" | "black_box_amflow",
  "options":     { ... },             // optional global options
  // ... 按 mode 不同的其它键 (见下文)
}
```

`mode` 缺省时按 `"amflow"` 处理. 输出总是含 `mode` (回显), `options` (应用输入 `options` 后的有效值), `result` (按 mode 不同的结构).

---

## Mode `amflow` — raw ODE 求解

镜像上游 `DESolver.m` 的 `AMFlow[de, bc]`. 用 η = ∞ 处的渐近边界规约求解 `dI/dη = A(η)·I`, 返回 η = 0 处的值. 不涉及 IBP, Kira, family.

### 输入

| 键          | 必填 | 类型   | 含义 |
|---          |---       |---     |---      |
| `mode`       | optional | string | 必须是 `"amflow"` (或缺省) |
| `matrix`     | **required** | 二维数组 | rational ODE 矩阵 `A(η)`, 方阵 `N × N`; 见 [Rational 项](#rational-项) |
| `boundaries` | **required** | 数组的数组 | `boundaries[i]` 是第 `i` 行的 BC 列表; 每项 `{ "mu": <complex>, "value": <complex> }` |
| `options`    | optional | object | 见 [全局 options](#全局-options) |

`boundaries.length` 必须等于 `matrix.length`.

### 输出

```jsonc
{
  "mode":   "amflow",
  "result": [ <complex>, <complex>, ... ],   // 长度 N
  "options": { /* 有效 options */ }
}
```

`result[i]` 是数值 `I_i(η = 0)`.

### 最小示例

```jsonc
// 1×1 系统, dI/dη = (1/η)·I, 边界 I ~ η^1 at infinity → I(0) = 1
{
  "options": { "x_order": 30, "extra_x_order": 10, "silent_mode": true },
  "matrix":     [[ { "num": ["1"], "den": ["0", "1"] } ]],
  "boundaries": [[ { "mu": { "re": "1", "im": "0" },
                    "value": { "re": "1", "im": "0" } } ]],
  "mode": "amflow"
}
```

(用 `amflow_cli examples/power_law.json` 运行 — 没装 CLI 时改用 `./build/src/cli/amflow_cli examples/power_law.json`.)

---

## Mode `solve_integrals` — ε 的 Laurent 展开

镜像上游 `AMFlow.m` 的 `SolveIntegrals[jints, goal, order]`. 完整 pipeline: family → ε 网格 → Kira IBP → AMFSystem 递归 → 采样数值 → Laurent 拟合.

### 输入

| 键            | 必填 | 类型   | 含义 |
|---             |---       |---     |---      |
| `mode`         | optional | string | 必须是 `"solve_integrals"` |
| `family`       | **required** | object | 积分 family — 见 [Family 对象](#family-对象) |
| `integrals` (别名: `jints`, `targets`) | **required** | 数组 | 目标积分 — 见 [Integral 项](#integral-项) |
| `goal_digits` (别名: `goal`) | **required** | int | Laurent 系数的精度 (十进制位数) |
| `eps_order` (别名: `order`)  | **required** | int | 拟合在算法 leading 阶 −2L 之上的 ε 阶数 (L = 圈数). 最高拟合 ε 幂为 `ε^(eps_order − 2L)`. 例: 1-loop 用 `eps_order=4` 拟到 ε² (阶 ε⁻² … ε²). 想要 "结果到 ε^K", 传 `eps_order = K + 2L`. 输出会 trim 掉 leading 处为零的阶. |
| `work_dir`     | optional | string  | Kira 中间文件所在目录 (默认: 临时目录) |
| `amf_options`  | optional | object | 见 [AMF options](#amf-options); 也承载 `numeric_values` |
| `options`      | optional | object | 见 [全局 options](#全局-options) |

`amf_options.blackbox.numeric_values` 通常用来给每个 kinematic invariant 赋值 (例如 `"s": "100"`, `"t": "-1"`).

### 输出

```jsonc
{
  "mode":   "solve_integrals",
  "result": [
    {
      "integral":      { "family": "...", "indices": [...] },
      "leading_order": -2,                      // 最小 ε 次幂
      "coefficients":  [
        { "order": -2, "value": { "re": "...", "im": "..." } },
        { "order": -1, "value": { ... } },
        { "order":  0, "value": { ... } },
        ...
      ]
    },
    ...
  ],
  "options": { /* 有效 options */ }
}
```

每一项即 Laurent 展开 `Σ_k coefficients[k].value · ε^(leading_order + k)`.

---

## Mode `black_box_amflow` — 单 ε 采样, 不做 Laurent 拟合

镜像上游 `BlackBoxAMFlow[jints, epslist, dir]`. 返回每个目标积分在每个用户给定 ε 上的数值, **不**做 Laurent 拟合. 适用于在单 ε 上做与 Mathematica AMFlow 的 parity 对比, 或你想自己控制 ε 网格的场景.

### 输入

| 键 | 必填 | 类型 | 含义 |
|---|---|---|---|
| `mode`         | optional | string | 必须是 `"black_box_amflow"` |
| `family`       | **required** | object | 见 [Family 对象](#family-对象) |
| `integrals` (别名: `jints`, `targets`) | **required** | 数组 | 见 [Integral 项](#integral-项) |
| `eps_samples` (别名: `epslist`) | **required** | 复数数组 | 要求值的 ε 值列表; 接受 `"1/100"` 这样的 rational |
| `work_dir`     | optional | string | 同 `solve_integrals` |
| `amf_options`  | optional | object | 见 [AMF options](#amf-options) |
| `options`      | optional | object | 见 [全局 options](#全局-options) |

### 输出

```jsonc
{
  "mode":   "black_box_amflow",
  "result": [
    {
      "integral": { "family": "...", "indices": [...] },
      "samples":  [
        { "eps":   { "re": "1/100", "im": "0" },
          "value": { "re": "...",   "im": "..." } },
        ...
      ]
    },
    ...
  ],
  "options": { /* 有效 options */ }
}
```

---

## Family 对象

`solve_integrals` 和 `black_box_amflow` 共用.

| 键            | 必填 | 类型 | 含义 |
|---             |---       |---   |---      |
| `name` (别名: `family`) | **required** | string | family 名 (作为 `j[<name>, …]` 使用) |
| `loops`        | **required** | string[] | 圈动量 (例如 `["l1", "l2"]`) |
| `legs`         | **required** | string[] | 外动量 (例如 `["p1", "p2", "p3", "p4"]`) |
| `propagators`  | **required** | string[] | propagator 分母, 用代数字符串写 (例如 `"l^2 - msq"`) |
| `conservation` | optional | object | 动量守恒重写规则, 动量 → 表达式 (例如 `{"p4": "-p1 - p2 - p3"}`) |
| `replacement`  | optional | object | scalar product → invariant 映射 (例如 `{"(p1+p2)^2": "s"}`); 用来把 kinematics 改写为 invariant 变量 |
| `cut`          | optional | int[] | 每个 propagator 一项, `1` = on shell / 切, `0` = 未切 |
| `prescription` | optional | int[] | 每个 propagator 一项; `+1` / `-1` 是两个 ±i0 符号选择 |

如果指定 `cut` 和 `prescription`, **它们的长度必须等于 `propagators`**.

### Family 示例

```jsonc
"family": {
  "name": "box1",
  "loops":   ["l"],
  "legs":    ["p1", "p2", "p3", "p4"],
  "conservation": { "p4": "-p1 - p2 - p3" },
  "replacement":  {
    "p1^2": "0", "p2^2": "0", "p3^2": "0", "p4^2": "0",
    "(p1+p2)^2": "s",
    "(p1+p3)^2": "t"
  },
  "propagators": [
    "l^2",
    "(l + p1)^2",
    "(l + p1 + p2)^2",
    "(l + p1 + p2 + p4)^2"
  ]
}
```

---

## Integral 项

`integrals` / `jints` / `targets` 的每一项写成下面之一:

```jsonc
// (a) 仅 indices 数组 — 使用外层 family
{ "indices": [1, 0, 1, 0] }

// (b) 覆盖 family
{ "family": "box1", "indices": [1, 0, 1, 0] }
```

`indices` 长度必须等于解析后的 family 的 `propagators` 长度. 负数 indices 编码 irreducible scalar products (例如 `[-2, 1, 1, 1]` 表示在 slot 0 上 ISP rank 为 2).

---

## 全局 options

顶层 `"options"` 对象. 所有字段都是 **optional**; 默认值与上游 AMFlow 默认值一致.

| 键              | 类型           | 默认 | 镜像上游 |
|---               |---             |---      |---               |
| `working_pre`    | int (digits)   | 100     | `WorkingPre`     |
| `chop_pre`       | int (digits)   | 20      | `ChopPre`        |
| `rationalize_pre`| int (digits)   | 100     | `RationalizePre` |
| `silent_mode`    | bool           | false   | `SilentMode`     |
| `d0`             | rational ("p/q") 或 int | `"4"`  | `D0` — 宿主时空维度 |
| `x_order`        | int            | 100     | `XOrder`         |
| `extra_x_order`  | int            | 20      | `ExtraXOrder`    |
| `learn_x_order`  | int            | -1      | `LearnXOrder`    |
| `test_x_order`   | int            | 5       | `TestXOrder`     |
| `run_radius`     | int            | 2       | `RunRadius`      |
| `run_length`     | int            | 1000    | `RunLength`      |
| `run_candidate`  | int            | 10      | `RunCandidate`   |
| `run_direction`  | `"Re"`/`"Im"`/`"NegRe"`/`"NegIm"` | `"NegIm"` | `RunDirection` |

精度值在 API 层是**十进制位**, 实现内部转换成 binary bit (见 [`docs/INVARIANTS.md`](INVARIANTS.md) §1).

---

## AMF options

`"amf_options"` 是 `solve_integrals` / `black_box_amflow` 用的 AMFSystem pipeline 的每棵递归树的配置. 所有字段都是 optional.

| 键 | 类型 | 默认 | 含义 |
|---|---|---|---|
| `amf_modes`     | string[] | `["Prescription", "Mass", "Propagator"]` | η 注入候选 mode (镜像 `AMFMode`) |
| `ending_schemes`| string[] (`"Tradition"`/`"Cutkosky"`/`"SingleMass"`) | `["Tradition", "Cutkosky", "SingleMass"]` | 按顺序尝试的 ending scheme |
| `cache_root`    | string   | empty | AMF 子系统 cache 的根目录; 每个节点的子目录在下面创建 |
| `direction`     | string   | 上游默认 | 本次运行的路径方向覆盖 |
| `max_recursion_depth` | int | 无限 | AMF 子系统递归深度硬限制 |
| `blackbox` (别名: `bb`) | object | (见下) | Kira / Fermat / IBP 旋钮 |

### `amf_options.blackbox`

| 键 | 类型 | 默认 | 含义 |
|---|---|---|---|
| `numeric_values`   | object map<string,string\|number> | `{}` | kinematic invariant 的数值 — *非平凡 family 必填* (例如 `{"s": "100", "t": "-1"}`). **仅接受实数标量**; 复数对象形式 `{"re":..,"im":..}` 会被显式错误拒绝 (out of scope). |
| `ibp_rank`         | int | 上游默认 | `BlackBoxRank` 下限 |
| `ibp_dot`          | int | 上游默认 | `BlackBoxDot` 下限 |
| `n_thread`         | int | 上游默认 | Kira 的 `NThread` |
| `integral_order`   | int | 5 | Kira `integral_ordering` |
| `kira_executable`  | string | `/usr/local/bin/kira` | `kira` 可执行文件路径 |
| `fermat_executable`| string | `/usr/share/Ferl7/fer64` | `fer64` 可执行文件路径 |
| `work_dir`         | string | 按 mode 不同的默认值 | 在本层覆盖 work dir |
| `log_file`         | string | `/tmp/amflow_cli_kira.log` | Kira stdout / stderr 追加写入的位置 |

---

## Rational 项

用于 `mode: amflow` (矩阵元素). 两种等价形式:

```jsonc
// (a) 仅多项式, 分母隐式为 1
{ "coeffs": ["c0", "c1", "c2", ...] }

// (b) 显式 numerator / denominator, 都是 η 的多项式
{ "num": ["a0", "a1", ...], "den": ["b0", "b1", ...] }
```

每个系数被解析为 `fmpq` (rational); `"3/7"`, `"-1"`, 整数, 小数都可接受. 数组 index `i` 表示 `η^i` 的系数.

示例:

```jsonc
// 1 / η  →  num = 1, den = 0 + 1·η = η
{ "num": ["1"], "den": ["0", "1"] }

// (3 + 5η - η²) / 1
{ "coeffs": ["3", "5", "-1"] }
```

---

## 复数值

凡是出现单一复数的地方 (boundary `mu` / `value`, `eps_samples` 项, 输出) 都用以下两种等价形式之一:

```jsonc
// (a) 显式 real / imag (各自解析为 fmpq 或 arb decimal)
{ "re": "1/2", "im": "0" }

// (b) 纯实数 — 虚部 = 0
"1/2"
"-1.5"
3
```

输出端的复数总是用 (a) 形式.

---

## 错误处理

Dispatcher 在任何 schema 或解析错误上抛 `std::runtime_error`; `amflow_cli` 捕获并把异常文本写到 stderr, 非零退出. 典型错误:

| 起始信息 | 原因 |
|---|---|
| `unknown mode '...'`           | `mode` 不是三个已知值之一 |
| `missing required array '...'` | 必填键缺失 |
| `... must be an array`         | 期望数组的键类型错了 |
| `... has N indices, expected M`| Integral indices 长度与 propagator 长度不符 |
| `boundaries length must match matrix size; got N vs M` | `mode: amflow` 下 boundary 列表基数不符 |
| `failed to parse fmpq from '...'` | 系数字符串不是合法 rational |

---

## 另见

- [`docs/USER_GUIDE.md`](USER_GUIDE.md) — 逐步教程.
- [`docs/REFERENCE_MAP.md`](REFERENCE_MAP.md) — Mathematica → C++ 符号映射 (哪些上游函数支撑哪些 JSON mode).
- [`examples/`](../examples) — 五个可运行的 JSON 输入.
- [`src/api/run_json.cpp`](../src/api/run_json.cpp) — 实现; 本文档与该文件必须保持同步.
