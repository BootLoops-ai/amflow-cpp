# 贡献与调试指南

> 英文版本: [CONTRIBUTING.md](CONTRIBUTING.md)

本文档面向实际工程工作: 构建, 测试, benchmark, 参考数据生成, parity 调试.

---

## 工作流 checklist

改代码时:

1. 构建当前工作树:
   ```bash
   cmake --build build -j32
   ```
2. 跑完整本地测试:
   ```bash
   ctest --test-dir build --output-on-failure -j 4
   ```
3. 测试失败时, 先读 assertion 信息. 多数数值测试已经把比较值打出来了.
4. 项目状态变化时同步更新 live docs: `README.md`, `docs/AUDIT_MMA_PARITY.md`, `CHANGELOG.md`.

加 / 改功能时:

1. 把代码放到 `src/<domain>/` 下对应的 domain 模块 (numeric / algebra / ode / qft / ibp / pipeline / api / cli).
2. 在 `tests/` 下添加或扩展对应测试文件.
3. 如果改动改了公开 surface 或项目状态, 同一次提交里同步更新 docs.
4. 如果改动关闭或暴露了 parity 问题, 保留比较中的 Mathematica 一侧.

## 当前验证策略

- Mathematica AMFlow 是规范.
- 当前的 parity 筛子是 `eps = 1/1000` 处的单点采样值对比.
- 当一个采样 parity 对比还开着时, 不要先去做更广的 fit / Laurent 覆盖.
- 用于确立 parity 的每一份 Mathematica 输出都要保留.

## 运行测试

单测:

```bash
./build/tests/amflow_tests \
  --gtest_filter='AmfsystemTest.*'
```

单 domain 组:

```bash
./build/tests/amflow_tests \
  --gtest_filter='*Numeric*:*Qft*'
```

列测试:

```bash
ctest --test-dir build -N
```

## 生成 Mathematica 参考数据

Wolfram driver 在 `tools/math_ref/run_amflow_kira.wl`, 从 `AMF_REF_CFG` 读配置.

```bash
cd tools/math_ref
AMF_REF_CFG=$PWD/cfg_bubble_1L.wl math -script run_amflow_kira.wl
```

直接参考 JSON 文件 committed 在 `tests/data/math_ref/`. 采样 parity benchmark 参考在 `tools/bench/`.

## 采样 parity 工作流

新 family 拉起来时用这个流程:

1. 选积分目标, 固定 `eps = 1/1000`.
2. 生成或保留 Mathematica 采样值.
3. 把汇总的值存到 `tools/bench/*_mma_reference.json`.
4. 凡是含中间值的 Mathematica 原始 cache, 都保留到 `tools/bench/mma_refs/`.
5. 在完全相同的 family / replacement / 数值 / 目标上跑 C++ 采样路径.
6. 如果 C++ 不一致, 先比较第一处错的直接子积分 / 子系统数值, 不要先改 fit / Laurent 代码.
7. 如果 C++ 一致, 在合理运行时间内添加或扩展聚焦的回归测试.

## 常见调试策略

### 结果是 NaN

通常是其中之一:

- `prec` 在某个 Arb 调用前被省略或变成 0;
- 某个 normalization 路径收到了无效精度或奇异变换;
- 路径生成器产生了 degenerate step;
- Kira 约化返回不完整数据, 下游代码盲信了.

第一波检查:

```bash
./build/tests/amflow_tests \
  --gtest_filter='OdeAmflowTest.*'
./build/tests/amflow_tests \
  --gtest_filter='AmfsystemTest.BubblePipeline_RunsToCompletion'
```

### 结果是有限值, 但差一个简单因子或符号

通常是其中之一:

- 1-based vs 0-based 的转换;
- 稀疏系统周围的反向系数索引;
- 递归 boundary 链里缺了 prefactor 传播;
- setup 之前对 family replacement 的 numeric 代入做错了.

先看 `docs/INVARIANTS.md`.

### 所有 master 都回零, 或 boundary 列表为空

找:

- 缺失的 pattern-zero boundary 项;
- `deinf` 构造错;
- `determine_boundary_order` 入参类型错;
- ending / global value 拼接问题.

### 递归不收缩或 setup 死循环

检查:

- `ending_q` 与 `AMFlow.m` 的 parity;
- 子 family 是否真的减少了 loop 或 propagator;
- 子 setup 是否从完全代入的 family 数据构造, 而不是半符号的.

### `MpolyContext mismatch`

你混用了不同 context 创建的符号对象. 用字符串往返重新解析到目标 context. 不要按名字比较 context; 运行时规则是共享所有权同一性.

### Kira 报错或卡住

检查:

- `kira_path` 与 `fermat_path` 是不是绝对路径;
- `FERMATPATH` 目录里的内容;
- 用 `AMFLOW_KIRA_DUMP=1` dump 的 Kira workdir;
- 生成出来的 `integralfamilies.yaml`, `jobs.yaml`, `kinematics.yaml`, `preferred`, `target`.

## 调试开关

运行时环境变量 — 权威清单在 [`include/amflow/numeric/log.hpp`](../include/amflow/numeric/log.hpp); 引入新 trace category 时去那里加. 当前:

| Env var | 效果 |
|---|---|
| `AMFLOW_KIRA_DUMP=1` | 运行后保留 Kira workdir |
| `AMFLOW_JORDAN_DUMP=...` | dump Jordan 形式中间矩阵到文件 |
| `AMFLOW_NO_SPARSE_CHOP=1` | 关掉稀疏行 chop 启发 |
| `AMFLOW_SPARSE_CHOP_DIGITS=N` | 覆盖稀疏 chop 位数 (默认 = `numeric::chop_pre()`) |
| `AMFLOW_DEBUG_BC=1` | 发 boundary-condition 与 `solve_one_eps` 诊断 |
| `AMFLOW_DEBUG_DBO=1` | 发 `determine_block_boundary_order` 诊断 |
| `AMFLOW_DEBUG_PRECISE=1` | normalize / Jordan 路径里的精细精度追踪 |
| `AMFLOW_DEBUG_REGULAR=1` | 发 `regular_run` 逐步追踪 |
| `AMFLOW_DEBUG_STAGES=1` | 发逐 stage 系统追踪 (Inf→Reg, RegRun, Zero) |
| `AMFLOW_DEBUG_CT_PERM=1` | 发 `ct_perm` Cutkosky permutation 追踪 |
| `AMFLOW_DEBUG_SCHEME=1` | 发 `[scheme] ...` EndingScheme dispatch 追踪 (`amf_system_setup_master` 里 Cutkosky / SingleMass 触发) |
| `AMFLOW_TRACE_LAYER6=1` | `src/ode/inf.cpp` 的细粒度追踪 (boundary-Taylor / Frobenius 诊断) |

## 加测试

按现有模式:

```cpp
class Layer5Test : public ::testing::Test {
protected:
    void SetUp() override { set_default_options(); }
    void TearDown() override { set_default_options(); }
};
```

测试里一定要 reset options. 很多看起来像数值问题的失败, 其实是泄露的全局设置.

数值检查优先用测试套已有的辅助函数, 别拿字符串比较. 只在工作流确实积累了更多数值误差时才放宽容差.

## 不要做的事

- 改上游 Mathematica AMFlow 去配合 C++ 行为 (上游是规范, 我们顺它, 不是反过来).
- 用更弱的测试或更松的容差替换 parity 工作.
- 没重新解析就跨 context 混用符号.
- 不调查不一致就把一个 ported 的 Mathematica 路径从 parity 检查里排除掉.
