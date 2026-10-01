# FAQ

> 英文版本: [FAQ.md](FAQ.md)

AMFlow.cpp 的使用者与打包者常见问题. 如果你的问题不在这里, 请到 [GitHub issue](https://github.com/chang18/amflow-cpp/issues) 反馈, 或联系维护者 <3250800970@qq.com>.

---

## 构建与安装

### Q. CMake 报错 `Could not find module 'flint' in version >= 3.0`.

你的系统 FLINT 太老了. 本项目用到的符号 (`fmpz_mpoly_q_set_str_pretty`, `fmpz_mpoly_q_evaluate_acb`, …) 是 FLINT 3.x 才有的. 注意 **Ubuntu 22.04 装的是 FLINT 2.8.4**, 不能用. 选项:

- **源码构建 FLINT 3.x** (~5–10 min):
  ```bash
  curl -fsSL https://github.com/flintlib/flint/releases/download/v3.4.0/flint-3.4.0.tar.gz \
       | tar -xz
  cd flint-3.4.0
  ./bootstrap.sh && ./configure --prefix=/opt/flint --disable-static
  make -j$(nproc) && sudo make install
  echo "/opt/flint/lib" | sudo tee /etc/ld.so.conf.d/flint.conf
  sudo ldconfig
  export PKG_CONFIG_PATH=/opt/flint/lib/pkgconfig:$PKG_CONFIG_PATH
  ```
  然后重新跑 `cmake -S . -B build` — pkg-config 会找到新的 FLINT.
- **改用 Ubuntu 24.04 / Fedora / Arch / Conda-forge**, 这些发行版都打包了 FLINT 3.x. CI 的可复现配方见 `.github/workflows/ci.yml`.

### Q. 构建成功但运行时 linker 报缺失 `acb_*` / `arb_*` / `flint_*` 符号.

跟上面同源 — 系统上有多份 FLINT, 构建用了一份, 运行时 loader 选了另一份. 确保 `LD_LIBRARY_PATH` (或 `/etc/ld.so.conf.d/*`) 暴露的是 configure 时 pkg-config 选中的同一份 FLINT.

### Q. 没装 Kira 能构建吗?

可以. Kira 和 Fermat 只是 `solve_integrals` 和 `black_box_amflow` 的**运行时**依赖. `cmake --build build -j` 在没它们的环境里也能跑通. 548 个 GoogleTest case 在没 Kira 时也能跑 — 依赖 Kira 的 test 在 `kira` 二进制不在 `$PATH` 时会自动 skip.

Raw `amflow` mode (ODE 引擎) 完全不需要 Kira — 可跑示例见 [`examples/power_law.json`](../examples/power_law.json) 和 [`examples/constant.json`](../examples/constant.json).

### Q. 怎么把 library 和 `amflow_cli` 二进制装到系统?

构建成功后:

```bash
sudo cmake --install build                       # 默认 prefix /usr/local
# 用户级安装 (不需要 sudo):
cmake --install build --prefix "$HOME/.local"
```

会在选定的 prefix 下装四样东西:

- `bin/amflow_cli` — JSON 驱动的 CLI, 进 `$PATH`.
- `lib/libamflow.a` — 静态库.
- `include/amflow/...` — 按 domain 划分的公开头文件.
- `lib/cmake/AMFlowCpp/` — 下游 `find_package` 用的 CMake config.

自定义 prefix:

```bash
cmake --install build --prefix /some/prefix
```

### Q. 怎么让下游项目的 `find_package(AMFlowCpp)` 找到本项目?

下游的 `CMakeLists.txt`:

```cmake
list(APPEND CMAKE_PREFIX_PATH /some/prefix)   # 仅非标准路径需要
find_package(AMFlowCpp 1.1 REQUIRED)
target_link_libraries(my_target PRIVATE AMFlowCpp::amflow)
```

导出的 target 名是 `AMFlowCpp::amflow` (注意 `amflow` 是小写), namespace 是 `AMFlowCpp::`.

---

## 运行

### Q. `amflow_cli` 在产生 JSON 输出之前就退出报错了.

CLI 把每个错误包装成单行写到 **stderr**. 常见错误:

| 起始信息 | 含义 |
|---|---|
| `error: cannot open <path>` | 输入文件名错 |
| `error: failed to parse JSON: ...` | 输入文件不是合法 JSON |
| `unknown mode '...'` | `mode` 字段不是 `"amflow"`, `"solve_integrals"`, `"black_box_amflow"` 之一 |
| `missing required array '...'` | 缺少必填键 |
| `<key> must be an array / a string / ...` | JSON 类型错了 |

重跑时把 stderr 一起重定向 (`2>&1`) 看具体错误.

### Q. Kira 步卡死永远不返回.

症状: `amflow_cli` 持续高 CPU 但从不返回.

最常见原因:

- **Fermat 路径错.** `amf_options.blackbox.fermat_executable` (或 `FERMATPATH` 环境变量) 应该指向**目录** (含 `fer64`), 不是二进制文件本身.
- **Kira 找不到 Fermat.** 直接在 dump 出来的 workdir 上跑 Kira 看真正的错:
  ```bash
  cd /your/work_dir
  kira jobs.yaml
  ```
- **workdir 在慢盘上.** Kira IO 很重. 用本地 SSD 路径, 不要用网络挂载.

要 dump Kira workdir 以便检查 (默认会清理掉):

```bash
AMFLOW_KIRA_DUMP=1 ./build/src/cli/amflow_cli your_input.json
```

### Q. 结果是 `NaN` 或 `inf`.

见 [`docs/CONTRIBUTING.md`](CONTRIBUTING.md) §"Common debugging strategies". 常见原因: 工作精度被夹到 0, 边界规约写坏, 或 Kira 约化返回不完整数据.

### Q. 结果是有限值, 但跟 Mathematica AMFlow 差一个常数因子.

可能原因:

- **`propagators` 符号约定错.** 本项目用 Minkowski-`+i0` 约定, propagator 写成 `momentum² − mass²`. 如果你的参考是 Euclidean signature, 符号翻转.
- **没给所有 kinematic invariant 赋数值.** 在 `replacement` 里出现的每个 invariant + 在 `propagators` 里出现的每个裸 mass 都必须出现在 `amf_options.blackbox.numeric_values`.
- **`D0` 不匹配.** 默认时空维度是 4, 如果你在非整数维方案下, 设 `options.d0 = "<rational>"` 匹配上.

### Q. 怎么调高 / 调低工作精度?

`options.working_pre = N` 设 working precision 为 N 位十进制 (默认 100). `options.chop_pre = M` 设 chop tolerance 为 `10^-M` 位 (默认 20). 两者都是 API 层的**十进制**; 实现内部会转成 bit.

精度策略详见 [`docs/INVARIANTS.md`](INVARIANTS.md) §1.

### Q. 能在 4 以外的维度算吗?

可以. 宿主时空维度是 `D = D₀ - 2ε`, `D₀` 通过 `options.d0` 可配 (rational `"p/q"` 或整数; 默认 `"4"`):

```jsonc
"options": {
  "d0": "7/3"      // → D = 7/3 - 2ε
}
```

`solve_integrals` 和 `black_box_amflow` 都支持 (镜像上游 AMFlow.m 的 `D0` 选项). 引擎内部在 `ε + (4 - D₀)/2` 上对 family 求值, 然后拟合 Laurent 展开回用户面向的网格, 所以 `D₀` 在 API 层干净抵消.

在 oracle bench [`tools/bench/box1_d0_7_3_solve_integrals_cpp.json`](../tools/bench/box1_d0_7_3_solve_integrals_cpp.json) 上验证过 (D₀ = 7/3, 与 Mathematica AMFlow 吻合到 ~30 位有效数字). 工作示例见 [`docs/USER_GUIDE.md`](USER_GUIDE.md) §4.4, parity 数据见 [`AUDIT_MMA_PARITY.md`](AUDIT_MMA_PARITY.md).

---

## 与 Mathematica AMFlow 的数值一致性

### Q. C++ 结果与 Mathematica 结果在末尾几位不一样. 是 bug 吗?

几乎肯定**不是** bug. AMFlow 是个数值算法, 最末几位精确数字依赖于内部浮点选择. 本项目的一致性契约是:

- **单点采样对比** 在 `eps = 1/1000` (有时 `eps = 1/100`) 上, 在扣除 working-precision 噪声后, 相对误差与上游一致到 `~10^-30`.

如果你看到 `working_pre = 100` 下大于 ~10⁻²⁵ 的差异, 这值得调查. 已验证的 family 列表和观测到的相对误差见 [`AUDIT_MMA_PARITY.md`](AUDIT_MMA_PARITY.md) 和 [`tools/bench/`](../tools/bench/).

### Q. 怎么重新生成 committed 的 Mathematica 参考数据?

你需要:

- 一份能用的 **Mathematica / Wolfram Engine** 安装.
- 一份**上游 AMFlow 的本地 clone**, 放在 `reference/amflow-master/` (本项目不打包上游). 两种合法布局 (clone-into 或 symlink) 见 [`reference/README.md`](../reference/README.md).

然后, benchmark cache:

```bash
cd tools/bench
math -script <name>_mma.wl     # 写入 <name>_mma_reference.json
```

或更大的 reference set:

```bash
cd tools/math_ref
AMF_REF_CFG=$PWD/cfg_bubble_1L.wl math -script run_amflow_kira.wl
```

完整细节见 [`tools/bench/README.md`](../tools/bench/README.md).

---

## 项目 scope 与贡献

### Q. 哪些功能**没**实现?

**不在 scope 内 (有意决定, 不会加):**
- `SolveIntegralsGaugeLink`
- HQET / SCET / Wilson-line 工作流
- Kira 以外的 IBP 后端 (不上 FIRE / LiteRed / FiniteFlow / Blade)
- **复数值数值 kinematics.** 上游 MMA 通过 `IBPRule` / `CompensateRule` 过滤虚部 (实部传 Kira, 虚部在之后代入). 本 C++ port 把 `amf_options.blackbox.numeric_values` 当作平直的实数 map; C++ algebra 层是 Q 上的 (FLINT `fmpz_mpoly_q_t`), 不能携带复数系数. JSON dispatcher 用清晰错误信息拒绝 `{"re":..,"im":..}` 对象形式. 绕过方式: 只提供纯实数. 见 [`docs/AUDIT_MMA_PARITY.md`](AUDIT_MMA_PARITY.md).

**另见:**
- [`docs/AUDIT_MMA_PARITY.md`](AUDIT_MMA_PARITY.md) — 行级上游一致性审计 (完整 divergence inventory).
- [`docs/REFERENCE_MAP.md`](REFERENCE_MAP.md) — Mathematica → C++ 符号映射 (含显式的 "not ported" 条目).

### Q. 为什么是 "AI 开发的" — 我能信任代码吗?

本项目的 C++17 源码主要由 AI 编程 agent 在维护者指导下开发. 可信度契约是**数值一致性**:

- `tools/bench/` 下 228/228 oracle benchmark 与上游 Mathematica AMFlow 在 `rel ~ 10⁻³⁰` 上吻合 (pentabox 2L 5-leg 角落 case 一致到 `rel ~ 10⁻¹⁰`, 因为该积分的内在取消尺度在固定 `WorkingPre` 下限制了可达精度); 此前延迟的 3 个 4L 2-leg `sunset_bubble` case 已在上游 4a292ab/d7b093d 中关闭.
- 548 个 GoogleTest case 把住公开 surface 的闸.
- 每个 commit 都通过这两道闸 (CI 在 push 时跑).

测试 (不是代码) 是契约. 测试对的话, 代码就是对的; 代码漂移的话, 测试能抓到. AI agent 的开发 / 验证指令见 [`AGENTS.md`](../AGENTS.md).

### Q. 我想加新 family / port 缺失的上游功能 / 修一个 parity bug.

按以下顺序读:

1. [`AGENTS.md`](../AGENTS.md) — 工作规则与 parity 契约.
2. [`docs/CONTRIBUTING.md`](CONTRIBUTING.md) — 实际工作流, 调试策略, env 变量.
3. [`docs/INVARIANTS.md`](INVARIANTS.md) — 精度 / RAII / 边界逻辑里不显眼的陷阱.
4. [`docs/REFERENCE_MAP.md`](REFERENCE_MAP.md) — 查上游函数的 C++ 对应.

然后开 PR. CI 必须保持绿色; 行为性变更必须有上游 parity 证据支撑.

---

## 引用

### Q. 怎么引用本软件?

**同时**引用上游论文 (算法) 与本仓库 (实现):

```bibtex
@software{amflow_cpp,
  title  = {{AMFlow.cpp}: A C++17 reimplementation of the AMFlow algorithm},
  author = {AMFlow.cpp contributors},
  year   = {2026},
  url    = {https://github.com/chang18/amflow-cpp},
  doi    = {10.5281/zenodo.20087172}
}

@article{liu_ma_amflow_2023,
  title   = {{AMFlow}: A {Mathematica} package for {Feynman} integrals computation via auxiliary mass flow},
  author  = {Liu, Xiao and Ma, Yan-Qing},
  journal = {Computer Physics Communications},
  volume  = {283},
  pages   = {108565},
  year    = {2023},
  doi     = {10.1016/j.cpc.2022.108565}
}
```

Zenodo DOI `10.5281/zenodo.20087172` 是 **concept DOI** — 始终指向最新发布. 每个版本的 DOI 列表见 [`CITATION.cff`](../CITATION.cff).
