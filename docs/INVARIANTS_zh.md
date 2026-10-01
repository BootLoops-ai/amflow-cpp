# 项目级不变式

> 英文版本: [INVARIANTS.md](INVARIANTS.md)

这些约定**不是**每个 layer 各自的 — 它们在整个 codebase 都要遵守. 违反其中任何一条都会静默产生错误结果.

> _注_: 下面 "Layer N" 标签指 v1.0 port 时的遗留模块组织. 当前源码树按 **domain** 组织 (`numeric` → `algebra` → `ode`/`qft`/`ibp`/`pipeline`/`api`/`cli`) 在 `src/<domain>/` 下; 当前布局见 [`ARCHITECTURE.md`](ARCHITECTURE.md).

---

## 1. API 层是十进制精度, 内部是二进制位

面向用户的 option 用**十进制位**:

```cpp
GlobalOptions g;
g.working_pre = 100;   // 表示 100 位十进制
```

FLINT/Arb 内部用**二进制 bit**. 转换通过 `Precision` value type:

```cpp
amflow::numeric::Precision::from_decimal_digits(int);
```

公式 `bits = ceil(digits · log2(10))`, 夹住 ≥ 53.

**规则:** 凡是接受 `prec` 参数的 Arb 调用, **每一次**都必须使用从 `Precision::from_decimal_digits` (或本项目其它 `Precision`-typed helper) 派生的值. 不要硬编码.

---

## 2. FLINT RAII 包装是 move-only; rational/poly 值可复制

`AcbValue`, `AcbVector`, `AcbMatrix`, `AcbPoly`, `FmpqPoly`, `RationalFunction`, `AsyTerm` 等包装是 *move-only*.

```cpp
AcbValue a; a.set_si(7);
AcbValue b = std::move(a);   // OK
AcbValue c = a;              // 编译错 (deleted copy ctor)
AcbValue d = a.clone();      // 显式深拷贝
```

原因:

- 意外复制一个百万 bit 的 `acb_t` 既看不见又昂贵.
- 我们希望 `std::vector<AcbValue>` 只通过 move 增长.

**Caveat:** `RationalComplex` *可以*复制, 因为只持有两个 `fmpq_t` (便宜). `BlockEquation`, `BlockEquationNum`, `AsyExpansion` 等因为持有 move-only 字段, 也继承了 move-only 性.

需要 `vector<AcbValue> v = something()` 时, 用 `std::move(...)` 或调辅助函数 `clone_vec(...)` (按需在各 layer 定义).

`std::vector<AcbValue>::assign(N, AcbValue())` **不**编译, 因为 `assign(n, value)` 需要复制. 改用:

```cpp
v.clear();
v.resize(N);   // 值初始化; 调 N 次 AcbValue() 默认 ctor
```

这里类型故意不同: `amflow::numeric::RationalPolynomial`, `RationalFunction`, `RationalFunctionMatrix` 是可复制的值对象. 这个边界是给精确代数所有权的, 让普通 C++ 值语义在更高级 domain model 里更简单. 任意精度 ball 包装当复制会隐藏大的精度成本时, 仍然保持 move-only.

---

## 3. Layer 3 / 6 / 7 的反向系数索引

`construct_matrix(dx, ax, totalorder)` 布置一个 `(T+2)·N` 列的稀疏矩阵. 列 index `k` 对应**未知积分级数的第 *(T+1-k)* 个系数**, **不是**第 *k* 个.

也就是:

```
                    column index k:  0           1           2         …    T+1
   represents the coefficient of:    f_{T+1}     f_T         f_{T-1}   …    f_0
```

这跟 Mathematica `DetermineBlockBoundaryOrder` 里 `fid[n]` 公式的约定一致:

```mathematica
fid[n_] := { ... [Mod[n-1, |new|]+1], ExtraXOrder + 1 - Quotient[n - 1, |new|] }
```

Layer 5 (`calcx1x2`) **不**用这个约定, 因为它用直接 in-place 递推, 从不碰稀疏列. 只有 Layer 3 风格的 routine (`construct_matrix`, `sparse_gaussian`, Layer 6 和 Layer 7 的 caller) 用.

要改一侧约定, 必须同时改另一侧.

---

## 4. Mathematica 1-based, C++ 0-based — `nh[All, n+1]` 变成 `nh[i][n]`

任何 .m 公式里写:

```mathematica
nh[[All, n+1]]   (* 1-based, 所以是第 (n+1) 个元素 = 0-based index n *)
```

C++ 译过来是 `nh[i][n]`, **不是** `nh[i][n+1]`.

Layer-5 覆盖很快能抓到这类翻译错误; 不要特例处理 index 算术.

---

## 5. Chop tolerance 是十进制, 应用在*中点*上

`AMFChop[x, 10^-ChopPre]` 把幅度低于 `10^-ChopPre` 的数清零. 我们的等价物:

```cpp
bool acb_is_chop_zero(acb_srcptr x, int chop_digits);
bool AcbValue::is_chop_zero(int chop_digits = chop_pre()) const;
```

两者都看复数球的**中点**. 它们*不*查 radius. 这跟 .m 语义对齐 — .m 是在 `N[..., WorkingPre]` (没误差球) 之后做 Chop.

如果你需要一个 interval-aware 的 "肯定为零" 测试, 直接用 `acb_contains_zero(...)` 构造; *不要*改 `is_chop_zero` 的语义.

---

## 6. Layer-7 Jordan decomposition 在 Q 上精确

`NormalizeMat` 用精确 rational residue matrix:

```text
RationalMatrix -> fmpq_mat residue -> jordan_decomposition_exact
```

如果特征多项式在 Q 上有 degree > 1 的不可约因子, exact-Q Jordan 路径必须明确失败. 当前代数 fallback 可能仍用精确 FLINT `qqbar` 数据做 leading 指数和 integer-floor 判定, 但公开的 `NormalizeMat()` 不能在下游表示只支持 rational 时, 静默假装完整代数 Jordan 基底存在. `CalcZero` 允许用单独的内部 `normalize_mat_for_calc_zero()` 路径, 因为该路径也把对应的 block-local 代数旋转一路带到 `calcx00()`.

---

## 7. `acb_mat_t` 是数组 typedef, 不是 class

```c
typedef acb_mat_struct acb_mat_t[1];
```

含义:

- `new acb_mat_t` **不**编译. 用 `new acb_mat_struct`.
- `std::vector<acb_mat_t>` 不编译 (数组不可赋值). 用 `std::vector<acb_mat_struct*>` (heap 指针) 或 `std::vector<acb_mat_struct>` 加手动初始化.
- 函数接受 `acb_mat_t mat` (数组退化成 `acb_mat_struct*`). 传变量名或 struct 指针; *不要*用 `&mat` (那是 `acb_mat_struct(*)[1]`, 错的类型).

同样适用于 `acb_t`, `arb_t`, `fmpq_t`, `fmpz_t`, `acb_poly_t`, `fmpq_poly_t`.

---

## 8. `fmpq_cmp_si(x, c)` 取一个整数, 不是分数

每个人都会被这个坑一次.

```c
int fmpq_cmp_si(const fmpq_t x, slong c);   // 比较 x 与 c (整数)
```

**没有** `fmpq_cmp_si(x, p, q)` 这种 "x 是不是等于 p/q?" 的版本. 用 helper

```cpp
bool fmpq_eq_pq(const fmpq* x, long p, long q);
```

(在测试里局部声明), 或者构造临时 `fmpq_t` 调 `fmpq_equal`.

---

## 9. `fmpq_mat_t` 遵循同样的数组-typedef 规则

```c
typedef fmpq_mat_struct fmpq_mat_t[1];
```

`jordan_decomposition_exact` 这类输出参数 API 因此接受预初始化的 `fmpq_mat_t` 对象. 不要直接返回或 move `fmpq_mat_t`; 在局部用 RAII 包装, 或从 caller 传入输出矩阵.

---

## 10. 精确 Jordan block 计数用 kernel-tower 差

`src/ode/jordan.cpp` 里, 对一个特征值:

```text
V_s = ker((A - lambda I)^s)
r[s] = dim(V_s) - dim(V_{s-1})
blocks_of_size[s] = r[s] - r[s+1]
```

在最大深度处, `r[depth+1]` 有意是 `0`, 因为没有比已发现的最大链更大的 block. 跟 `dim(V_{depth+1}) = dim(V_depth)` 混淆会破坏精确 block 计数.

Layer 7 对角 ToFuchsian 遵循同样的 exact-Q 规则: `ReduceL0`, `FindProjector`, `Balance`, `InvBalance` 在 rational 矩阵上操作. 不要把 projector 循环走 Arb 特征向量或数值 rank 测试; 不可约 leading Poincare rank 是精确代数错误.

---

## 11. `AcbRationalFunction` **不**做约分

跟 Layer 1 的 `RationalFunction` (维持 `gcd(num, den) = 1` 和 `den` monic) 不同, Layer 7 的 `AcbRationalFunction` (acb 系数 rational function) **不**约分. Acb 没有精确 GCD, 约分 ill-defined.

操作盲目堆叠 `num` 和 `den`. 通分清理只在 η = 0 求值时 (`value_at_zero()`) 发生一次, 它先剥离 leading η prefix, 再除常数.

需要检查结果时, 通过 `acb_rational_matrix_value_at_zero(...)` 落到 `acb_mat_t`.

---

## 12. `MpolyContext` 经 `shared_ptr` 共享; 跨 context 操作抛错

每个 `Mpoly` / `Mfrac` / `MpolyMatrix` / `MfracMatrix` 都带一个 `shared_ptr<MpolyContext>` 标识它所在的变量列表.

```cpp
auto ctxA = MpolyContext::make({"eta", "eps"});
auto ctxB = MpolyContext::make({"eta", "eps", "D[1]", "D[2]"});

Mpoly p = Mpoly::parse(ctxA, "eta + eps");
Mpoly q = Mpoly::parse(ctxB, "D[1]");

Mpoly r = p + q;   // 抛 std::runtime_error: ctx mismatch
```

要桥接两个 context, **通过字符串往返**:

```cpp
Mpoly p_in_B = Mpoly::parse(ctxB, p.to_string());   // ok
```

这是项目的纪律:

- Layer 12 携带 family 的 "base" context.
- Layer 13–14 一般用 `D[1] ... D[k]`, eta, region scale 扩展它 — 创建新 `MpolyContext`, 重新解析从父来的东西.
- Layer 16 的 child AMFSystem (尤其 SingleMass 子节点) 拥有自己的 context; 父在传 propagator 进来之前必须重新解析.

混用*碰巧*变量名相同的 context 没有警告 — 运行时检查用 `shared_ptr` identity. 这是有意的: 同名的两个 context 可能变量顺序不同, 不能静默重解释指数向量.

---

## 13. §3 的反向列约定在 Layer 14 的 boundary 表里要保留

Layer 14e (`boundary_integrals`) 最终把 boundary term 喂给 Layer 6 (`calc_inf` → `calc_taylor`). Term 存为 `(sub_index, coefficient)` tuple, `sub_index` 指 *sub-master* (子 family 上的 master integral).

Layer 16 的 `solve_one_eps` 在解时消费这个列表; 父 family 的 `BoundaryEntry { mu, value }` 必须遵守与 §3 相同的反向索引约定:

- `mu` 是 η 指数 (复数 rational), 直接传给 Layer 6 的 `read_bcs`.
- `value` 是在匹配 sub_index 上求值的子解的线性组合, 加到对应的 `BoundaryEntry`.

如果你改了 `construct_matrix` 的列方向约定 (§3), 必须同时翻转 Layer 16 从 `RegionBoundary::terms` 填 `BoundarySpec` 的顺序. 两者通过 `read_bcs` → `calc_taylor` → `sparse_gaussian` 耦合.

---

## 14. Kira 子进程契约

Layer 15c (`kira_run`) shell 出 Kira 二进制. 契约:

1. **只接受绝对路径.** `KiraConfig::kira_path` 和 `KiraConfig::fermat_path` 必须绝对. 我们在 `kira_run` 入口验证; 相对路径会静默失败, 因为 `fork()` 继承 caller 的 cwd, 然后 `execve` 重置它.
2. **`FERMATPATH` 环境变量.** 我们总是在子进程 env 里设 `FERMATPATH = cfg.fermat_path`. Kira 启动时用这个变量探测 Fermat; 缺失的话 Kira 在产生任何输出前 abort.
3. **Workdir 是 `cfg.work_dir`.** 我们在子进程 `execve` 前 `chdir(cfg.work_dir)`; 所有 yaml 文件已经预先写在那里.
4. **`BlackBoxReduce` 两轮.** 第一轮 `KiraJobMode::Masters` (只写 master 列表); 第二轮 `KiraJobMode::Reduce` (写 target 表). 两轮共享 workdir.
5. **Workdir 复用按 request 形状守护.** `ibp::black_box_reduce` / `black_box_diffeq` 路径在 Kira 跑成功后写 `.amflow_kira_fingerprint`. 指纹覆盖 family 形状, kinematics, job mode, 数值代入, IBP option, dimension base, preferred master, target. 下一个对同 workdir 的请求若匹配, 且预期的 Kira 输出文件齐全, cache 路径会直接解析这些文件不再启动 Kira. 请求不同或输出不全时, 只清掉 Kira 生成的产物 (`config`, `tmp`, `results`, `sectormappings`, `firefly_saves`, `preferred`, `target`, `jobs.yaml`) 加旧指纹, 再写新输入. 这样避免 Kira 辅助文件在不兼容的 numeric-variable 形状之间复用, 同时保留 caller 拥有的父目录.
6. **输出解析.** `kira_read_masters(cfg)` 解析 `<workdir>/results/<family>/masters`. `kira_read_target_table(cfg)` 经 `kira_parse_expression` 解析 `<workdir>/results/<family>/kira_target.m`, 支持 MMA 前缀语法和经 `fmpz_set_str` 解析大整数.

预飞失败 (缺二进制, 缺 yaml, 错 family 名) 抛 `std::runtime_error` 带 stderr. Post-Kira 失败 (Kira 跑了但返回非零) 抛 `std::runtime_error` 带 exit code 和 stderr 末 4 KiB.

调试 Kira 时, 设 env `AMFLOW_KIRA_DUMP=1` (或 `BlackBoxOptions::dump_kira_io = true`) — 跑完不会清 workdir.

---

## 15. SingleMass 子节点里 `global_values` 优先于 `master_values`

Layer 16 的 `AMFSystem::solve_one_eps` 从子节点读 boundary 贡献时, 优先顺序:

```cpp
const auto& src = child.solutions[eps_index];
const auto& boundary_values =
    src.global_values.empty() ? src.master_values
                              : src.global_values;
```

为什么: SingleMass 子节点的 `master_values` 是 *sub-master integral* (factorization 后每个连通分量一个) 的值. loop 重定义产生的 Γ-函数 prefactor 在 `solve_one_eps` 中应用, 单独存入 `global_values`. 父系统在拼接子系统时期待 *带 prefactor 的*值.

无条件读 `master_values` 会丢 Γ prefactor, 全局误差量级 `Γ(1+ε)^L`. **不要**简化这个条件分支, 除非你也在 call site 传播 prefactor.

---

## 16. `BuildTaylor` 必须在 `NHEquations` *之前*应用对角修正

Mathematica 的 `BuildTaylor[mat, ini]` 返回**完全修正后**的矩阵: `mat'[i,j] = η^{-ini[i]} · mat[i,j] · η^{ini[j]}` 若 `i ≠ j`, `mat'[i,i] = mat[i,i] - ini[i]/η` 若 `i = j` — 两者都通过 `Together` 保持 entry 为约简的 rational function. 这个修正后的矩阵再传给 `CalcTaylor`, `CalcTaylor` 内部对其调 `NHEquations`.

正确的 C++ 实现在 `src/ode/inf.cpp`:

- `build_taylor_symbolic(mat, int_offsets, ini_rat)` 在 `RationalFunction` (即 `fmpq`) 形式下做完整 BuildTaylor (非对角的 `η^(int_offsets[j] - int_offsets[i])` 移位*和*对角的 `mat[i,i] - ini_rat[i] / η` 减项), 通过现有 operator overload 应用 `Together`.
- `prepare_taylor_system(m_pure, …)` 然后在已修正的矩阵上跑 `nh_equations(…, Taylor)`. `PoincareRank`, `dx`, `ax` 都来自单一一致的计算.
- 任何代码路径都不能在 `NHEquations` *之后*从 `axexp[i][i][k+rank]` 数值地减去 `ini[i]·dx_inner[k]`.

**未来的具体不变式:** 任何未来对带非零 `ini` offset 矩阵跑 `nh_equations(..., Taylor)` 的代码路径, 必须**符号化地**预合成 `mat[i,i] - ini[i]/η` 减项. 不能有数值的事后修正. 如果未来用例需要复数值 `ini`, 扩展符号表示, 不要重新引入这个分裂.

---

## 17. AMFSystem 树的 `solve()` 自底向上; 永不 mutate 已访问的子节点

每个非 ending `AMFSystem` 节点持有:

- `region_boundaries` — `RegionBoundary { region, mu, terms }` 列表, `terms` 按 `sub_index` 引用子节点.
- `differential_eq_matrix` — 每个 eps 上的 ODE 矩阵 (`setup()` 通过 `BlackBoxDiffeq` 一次性建).
- `children` — `vector<AMFSystem>` (每个 region 每个 family 槽一个).

`solve(epslist)` 做:

```text
for each eps in epslist:
   1. for each child in children:
        递归: child.solve_one_eps(eps)
   2. 从 region_boundaries + 子解 build BoundarySpec
   3. 在 eps 上求值 differential_eq_matrix
   4. 调 amflow(de_at_eps, BC) → master_values
   5. (若 SingleMass 父) 应用 Γ prefactor → global_values
```

不变式: **一旦子节点的 `solve_one_eps(eps)` 返回, 该子节点在该 eps 上的解就固定了**. 不要从父节点 mutate 子节点的 solutions 数组. 多父原则上可以共享一个子节点 (DAG, 不是严格树), 所以 mutate 会破坏其它父节点. 当前结构是严格树, 但我们保持纪律, 让结构以后能变成 DAG 而不破坏任何东西.

加 caching, 并行 eps 求解, 或任何跨树 memoization 时这很关键: 自然位置是在子节点上, 对每个 (子节点, eps) 对正好写一次.

---

## 18. `EndingScheme` 在 `AMFSystemOptions::ending_schemes` 里的顺序很重要

`AMFSystemOptions::ending_schemes` 是 `vector<EndingScheme>`, 列出每个节点*尝试*的策略, 按顺序. 第一个返回 `ending_q == true` 的胜出.

默认顺序: `{Tradition, SingleMass}`. 这意味着能用 vacuum lookup 终结的节点就那样终结; 不行就试 SingleMass; 都不行就递归 (= Kira + 导数 + 子 AMFSystem).

如果你 reorder — 例如 `{SingleMass, Tradition}` — 同 family 会产生不同 (但数学等价) 的 AMF 树, 运行时代价 / 数值精度可能戏剧性变化. 调数值误差时, **先**强制单一 ending scheme (只 `{Tradition}` 或只 `{SingleMass}`) 定位问题.

---

## 19. `factor_poly(Taylor, rank)` 字面使用 `η^(rank+1)` — 不夹

朴素 port 里 `nh_equations` 用的 per-block `factor` 会是

```
factor = η^(rank + 1)    for Taylor mode
```

其中 `rank = poincare_rank(mat)` 可能是 `-1` (矩阵在 η=0 处 regular).

**正确做法:** `rank = -1` 时 `factor = η⁰ = 1`. Mathematica 字面这么做, 我们的 `factor_poly` 现在也一样.

**反模式:** 用 `if (power < 1) power = 1;` 夹住, 即使矩阵 regular 也强制 `factor = η`, 给 `dx` 灌一个虚假的 `η` 因子, 重新触发 INVARIANTS §16 旨在防止的列移位问题.

保留 `factor = η^max(0, rank+1)` 语义 (或等价的 — `rank+1 ≤ 0` 时多项式必须是 `1`), 别再引入夹.

---

## 20. `determine_boundary_order` 接受 rational `power` — 永不 round

`include/amflow/inf.hpp` 暴露两个 overload:

```cpp
std::vector<long>
determine_boundary_order(const RationalMatrix& mat,
                         const std::vector<long>& power);

std::vector<long>
determine_boundary_order(const RationalMatrix& mat,
                         const std::vector<RationalFunction>& power_q);
```

**`build_boundary` 用 rational overload.** Mathematica 的 `DetermineBoundaryOrder` 接受任意 `ini` (含非整数 rational, 如 `1 - ε` 在 `ε = 1/100` 处求值得到 `-99/100`). round 到最近的 `long` 给出不同答案:

```
MMA  DetermineBoundaryOrder[deinf, {-99/100}]  = {0}
MMA  DetermineBoundaryOrder[deinf, {-1}]       = {-1}
```

如果你的代码路径算出 rational `npat`, **不要** round. 把 `RationalFunction::from_fmpq(q)` 传给 rational overload.

C++ 入口是 `ode::determine_boundary_order(...)`, 接受 `std::vector<numeric::RationalFunction>` 指数. 把 AMFlow boundary order 接到 execution 层时保持 rational exponent 路径.

---

## 21. `deinf = -1/η² · de(1/η)` — 两半都要

`RationalFunction::substitute_inverse_eta()` **只**实现 `r(1/η)`. 它*不*乘 `-1/η²`.

Layer 6 的 helper `make_deinf(de)` 把两步包起来:

```cpp
auto subbed = de.substituted_inverse_eta();
RationalFunction prefactor = RationalFunction::monomial(-2, neg_one);
out(i, j) = subbed(i, j) * prefactor;
```

Layer 16 (`AMFSystem::build_boundary`) 里同样的 `deinf` 从每个 per-ε 求值后的 `Mfrac` inline 构造; 该代码也必须显式应用 `-1/η²` 因子. 漏掉它会让 `deinf` 在 η 上是多项式, 而 MMA 在 η=0 处有简单极点 — `determine_boundary_order` 看到完全不同的矩阵, 返回胡乱的 border, 非零 sector region 的 `rb.integrand_terms` 最后会空.

---

## 22. `solve_one_eps` 里的 Pattern-zero BC 補

Mathematica 的 `AMFSystemSolution` 以 (AMFlow.m line 1149) 结束:

```mathematica
bc = MapThread[Join[#1,#2]&,
               {bc, Transpose[Thread /@ Thread[(pattern/.epsrule) -> 0]]}];
```

每个 pattern group 的 per-master μ 值必须出现在 `bc_sorted[i]` 中, 值 `0` — 即使没有 boundary integrand 贡献匹配该 μ. 没有 pattern-zero 项的话, 一个 region 全 `border = -1` 的子系统会得到空的 `bc_sorted`, 然后 `amflow(de, {}, prec)` 返回恒为零.

我们的实现把 pattern 值存在 `AMFSystem::bc_pattern_` (在 `build_boundary` 里由 `boundary_pattern(all_powers)` 结果填), 在调 `amflow(de, bc_sorted, prec)` **之前**, 在 `solve_one_eps` 里追加 `(μ_i → 0)` 项. 不要跳过这一步.

---

## 23. `BoundaryFamily::terms` 布局是 `[input][order][laporta_term]`

`boundary_integrals` 的输出 — 文档在 `include/amflow/boundary.hpp` — 有三层嵌套维度:

```cpp
std::vector<std::vector<std::vector<LaportaTerm>>> terms;
//  └─ input_i ─┴─ order_j ─┴─ laporta_term_k ─┘
```

其中:

- `input_i` 遍历输入 integral (= 从 `build_boundary` 调用时的 `diff_masters` 条目);
- `order_j` 在 `0 .. border[input_i]` 上遍历 (`border[i] = -1` 时内层 vector 长度 0);
- 每个最内 entry 是 `LaportaTerm` = `{ indices, coef }`.

Mathematica 的 `BoundaryIntegrals` 等效地做 `DynamicPartition[flat_apart_pieces, partition]` 产生这个形状; 内层 `GatherBy` + `Plus` 聚合把同 `(input_i, order_j)` 共享 family 签名的所有 PFD 片段合成单个 `Mfrac` 再过 `LaportaIntegrals`.

**反模式:** 按平直 `pfd_per_input[flat_idx]` 布局 (`flat_idx` 把 `(input, order)` 当作单个序列, 内层是 per-PFD-piece) 索引 `terms`, 会让 `terms.size()` 不等于 `diff_masters.size()`, 导致 `build_boundary` 循环要么跳过条目要么写到错的槽. 平直布局不匹配 MMA 输出形状, 也没有下游消费者.

---

## 24. Layer 14 契约只支持 standard-QFT 语义

C++ port 只覆盖 standard QFT (二次 propagator) — 线性 propagator / gauge-link / Wilson-line 工作流不在 scope (见 `docs/ROADMAP.md`). 两个 Layer-14 入口强制这个契约, 在非 standard-QFT 形状的输入上**必须抛错**:

- `branch_momenta(family)` 要求每项 `coe1 == 1`; `1/2 * coe2` Wilson-line 折半分支没实现.
- `sp_list_to_dlist_symbol(family, sp_list)` 拒绝双线性系数不是常整数的 sp-list 项.

两者在 linear-propagator / gauge-link 输入上抛错. 不要放宽这些检查去 "支持" 混合 quadratic/linear family — pipeline 其它部分 (Layer 13 candidate selection, Layer 14 region decomposition, Layer 16 ending detection) 没有 gauge-link 分支.
