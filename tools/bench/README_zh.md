# 采样 Benchmark 资产

> 英文版本: [README.md](README.md)

本目录存放 `eps = 1/1000` 的采样 parity 输入, Mathematica 参考值, C++ benchmark 配置.

Mathematica AMFlow 是参考实现. 保留这些 benchmark 产生的每个 Mathematica 值与原始 cache.

## 对比已有输出

`compare_sampled.py` 把 C++ 采样输出 JSON 与 committed Mathematica 参考 JSON 对比. 它读取 requested / top / master 参考项, 报告每个 integral 的实部 / 虚部绝对与相对误差.

```bash
tools/bench/compare_sampled.py \
  tools/bench/tt_2loop_box_eps001_mma_reference.json \
  /tmp/amflow_bench_tt_2loop_box_cpp_out.json

tools/bench/compare_sampled.py \
  tools/bench/doublebox_sv_eps001_mma_reference.json \
  /tmp/amflow_bench_doublebox_sv_eps001_cpp_out.json
```

默认:

- 相对容差: `1e-25`
- 绝对容差: `1e-40`

用 `--json-out <path>` 把对比报告写成 machine-readable 文件.

## 跑长时 C++ Benchmark

下面命令跑采样 C++ benchmark 路径. 它们是长时 benchmark 命令, 不是默认 `ctest` 项.

下面命令假设 `amflow_cli` 已在 `$PATH` 上 (跑过 `sudo cmake --install build`). 没装的话, 把每处 `amflow_cli` 替换为 `./build/src/cli/amflow_cli` 即可.

```bash
amflow_cli \
  tools/bench/tt_2loop_box_eps001_black_box_amflow_cpp.json \
  /tmp/amflow_bench_tt_2loop_box_cpp_out.json

amflow_cli \
  tools/bench/doublebox_sv_eps001_black_box_amflow_cpp.json \
  /tmp/amflow_bench_doublebox_sv_eps001_cpp_out.json

amflow_cli \
  tools/bench/tt_higher_rank_eps001_black_box_amflow_cpp.json \
  /tmp/amflow_bench_tt_higher_rank_eps001_cpp_out.json

amflow_cli \
  tools/bench/xbox_2loop_eps001_black_box_amflow_cpp.json \
  /tmp/amflow_bench_xbox_2loop_eps001_cpp_out.json

amflow_cli \
  tools/bench/tt_cutkosky_probe_eps001_black_box_amflow_cpp.json \
  /tmp/amflow_bench_tt_cutkosky_probe_eps001_cpp_out.json

AMFLOW_DEBUG_SCHEME=1 amflow_cli \
  tools/bench/cutbubble_1L_eps001_black_box_amflow_cpp.json \
  /tmp/amflow_bench_cutbubble_1L_eps001_cpp_out.json

AMFLOW_DEBUG_SCHEME=1 amflow_cli \
  tools/bench/cutsunrise_2L_eps001_black_box_amflow_cpp.json \
  /tmp/amflow_bench_cutsunrise_2L_eps001_cpp_out.json

AMFLOW_DEBUG_SCHEME=1 amflow_cli \
  tools/bench/cutbanana_3L_eps001_black_box_amflow_cpp.json \
  /tmp/amflow_bench_cutbanana_3L_eps001_cpp_out.json
```

跑完后, 在更新 benchmark 状态字段前用 `compare_sampled.py` 比较输出.

## 采样 Benchmark 清单

活清单在文件本身: 每个 bench 是 `<name>_{cpp.json, mma.wl, mma_reference.json}` 三元组 (sub-master oracle 复用父 cache, 只有 `cpp.json + mma_reference.json`). 用 `ls tools/bench/*_mma_reference.json` 取完整列表 (截至 2026-05-17 共 228 项).

Wallclock 数字有意不在这里维护; 见 [`docs/ROADMAP.md`](../../docs/ROADMAP.md) §"性能 benchmark".

## 备注

- `tt_cutkosky_probe` 与 `tt_2loop_box` 只在 `EndingScheme` 优先级列表上不同. 在这个 kinematic (s = 30, t = -10/3, msq = 1) 下 boundary integral 由 `Tradition` 处理; 列表里 `Cutkosky` 的存在在协议层被锻炼但不扰动数值.
- `cutbubble_1L` 是 `tt_cutkosky_probe` 的*反向*对照: family 在顶层已有 `cut = {1, 1}`, 所以顶 sector 满足 `phase_volume_q` (|var| = L+1, all cut == 1), `AMFSystemSetupMaster` 在 root 系统上选 `Cutkosky`. 用 `AMFLOW_DEBUG_SCHEME=1` 跑, 应该看到 firing line `[scheme] Cutkosky fired: family=cutbubble phase_loop_num=1 ...` 被打印. eps^0 leading 系数恰好是 `1/(8*Pi)`, 即无质量 2-particle phase-space volume.
- `cutsunrise_2L` 把 Cutkosky probe 扩到 `L = 2`: 5-prop 基 (3 真 + 2 ISP), `cut = {1, 1, 1, 0, 0}`, target `j[cutsunrise, 1, 1, 1, 0, 0]`. `phase_volume_q` 触发因为 `|var| = 3 = L+1` 且所有活 cut == 1, 所以 `AMFSystemSetupMaster` 在顶 family 选 Cutkosky. firing line: `[scheme] Cutkosky fired: family=cutsunrise phase_loop_num=2 top_position={0,1,2}`. 这个 bench 是当前最强的 `(Pi^(2-eps)*(2*Pi)^(2*eps-4))^L * (-1)^(L+1)` Cutkosky prefactor 在 `L > 1` 下的检查 (符号因子相对 1-loop case 翻转), 也是 `phase_loop_num` 经 `AMFSystemSetupMaster` 传播的检查.
- `cutbanana_3L` 把 Cutkosky probe 推到 `L = 3`: 9-prop 基 (4 真 + 5 ISP), `cut = {1, 1, 1, 1, 0, 0, 0, 0, 0}`, target `j[cutbanana, 1, 1, 1, 1, 0, 0, 0, 0, 0]`. cut 清空后 eta-flow 作为 3-loop banana 跑 (eta 在 `l1^2` 上, 1 个 massive + 3 个无质量 leg); 子系统递归 `system 1 -> 2 -> 3 -> trivial`. firing line: `[scheme] Cutkosky fired: family=cutbanana phase_loop_num=3 top_position={0,1,2,3}`. 与 `cutbubble_1L` (符号 +1) 和 `cutsunrise_2L` (符号 -1) 一起, 完成了 L=1/2/3 上 `(-1)^(L+1)` 的符号规律检查 (L=3 时回到 +1), 也是 bench 清单里最深的递归. **需要 `BlackBoxDot = 5`** (匹配等质量 `banana_3loop`): `BlackBoxDot = 1` 时 boundary Reduce 步会失败, 报 `master count differs between Masters and Reduce calls`. MMA `*_mma.wl` 与 C++ `*_cpp.json` 配置都显式设了这个 floor — 当输入 integral 全部 `JDot = 0` 时 (每个正 index 都是 1), AMFlow.m 的 `Max[$BlackBoxDot, JDot/@...]` 预计算无法自己抬高它.
- `doublebox2m` 是 D12 的回归守 (见 `docs/AUDIT_MMA_PARITY.md`). 交错的跨 loop mass 放置给顶 sector 微分方程矩阵带来一对共轭复极点 (`η² + η + 12 = 0`, `|η| ≈ 3.46`), 收紧了顶 sector 4 个 master 的 Frobenius 级数收敛半径. committed 的 C++ 配置设 `working_pre=200, x_order=400, extra_x_order=480` — 相对全无质量 / block-mass 默认大约翻倍 — 才收敛. 在旧默认 `(160, 200, 240)` 下 C++ 产生反向的 Re 和 `O(0.1)` 量级的虚假 Im; 抬高的参数是任何 NegIm contour 附近有共轭复极点的拓扑的推荐 baseline.
- `banana_3loop` 在 Masters 和 Reduce 两轮都需要 `BlackBoxDot = 5` — 低 dot 时 Reduce 步无法覆盖 Masters 返回的 master 集合, AMFlow 的协议级 sanity 检查拒绝继续. committed 配置 (`*_cpp.json` 与 `*_mma.wl`) 设 `dot = 5`. Kira 的 Pak symmetry 检测 (`search_symmetry_relations`) 也必须为该 family 开启; 关掉 symmetry 时一些 sub-banana sector 不会塌缩, Reduce master count 又会漂.

## 重新生成 Mathematica cache (`mma_refs/`)

`mma_refs/` 下的原始 Mathematica cache 较大 (~80 MB), **未** committed. 它们能从本目录的 `*_mma.wl` 脚本重生. 对每个 bench 三元组 `<name>_{cpp.json, mma.wl, mma_reference.json}`:

```bash
cd tools/bench
math -script <name>_mma.wl
```

每个 `*_mma.wl` 把原始 cache 目录写到 `mma_refs/<name>_mma_cache/`, 把汇总参考值写到 `<name>_mma_reference.json` (这个**确实** committed).

若 `tools/math_ref/cache/` 也缺 (同样 gitignored), 用 `tools/math_ref/run_amflow_kira.wl` 加 `cfg_*.wl` 配置文件来重生.
