# Mathematica ↔ C++ 符号对照

> 英文版本: [REFERENCE_MAP.md](REFERENCE_MAP.md)

本文档把三份上游 Mathematica 源文件的符号映射到对应的 C++ 实现. 上游在 <https://gitlab.com/multiloop-pku/amflow>; 本项目不打包 (见 [`reference/README.md`](../reference/README.md)).

- 上游 `diffeq_solver/DESolver.m` (ODE 引擎, Layer 0–10). 见下文 Part I.
- 上游 `AMFlow.m` (Feynman 积分流水线, Layer 11–16). 见下文 Part II.
- 上游 `ibp_interface/Kira/interface.m` 和小的 Kira-glue 辅助. 见下文 Part III.

下面引用的行号对应上游 commit [`efda1db`](https://gitlab.com/multiloop-pku/amflow/-/commit/efda1db) ("add_version_number_1.2"), 即本 port 进行时的快照. 上游可能已往后推 (任何漂移见 [`AUDIT_MMA_PARITY.md`](AUDIT_MMA_PARITY.md) §5), 所以单个行号在更新的 revision 上可能对不上.
**有疑问请按符号名搜索.** 源码里的 `// Mirrors AMFlow.m:NNN` 注释也遵循同样约定.

需要翻译一个 Mathematica 表达式或追查为何结果不同时, 把本图当作导航辅助.

约定:

- 下面 "Layer N" 标签是 v1.0 port 时的**遗留模块名**. 当前源码树按 **domain** 组织在 `src/<domain>/` 下 (`numeric`, `algebra`, `ode`, `qft`, `ibp`, `pipeline`, `api`, `cli`); 当前组织见 [`docs/ARCHITECTURE.md`](ARCHITECTURE.md). Layer 6/7 ↔ `src/ode/`, Layer 11-16 ↔ `src/qft/` + `src/ibp/` + `src/pipeline/`, Layer 17 ↔ `src/cli/`.
- "—" 表示该函数没有 C++ 对应, 要么是 Mathematica 习惯用法 (例如 `MatrixDensity` 用于诊断), 要么是已重构掉的内部辅助.

---

# Part I — `DESolver.m` (Layer 0–10)

## `DESolver.m` 顶部 — Options 块 (lines 62-132)

| .m 符号 | C++ |
|---|---|
| `WorkingPre`, `ChopPre`, `SilentMode`, `RationalizePre` | `GlobalOptions` 字段 |
| `XOrder`, `ExtraXOrder`, `LearnXOrder`, `TestXOrder` | `ExpansionOptions` 字段 |
| `RunRadius`, `RunLength`, `RunCandidate`, `RunDirection` | `RunningOptions` 字段 |
| `SetGlobalOptions[opt]` | `set_global_options(GlobalOptions)` |
| `SetExpansionOptions[opt]` | `set_expansion_options(ExpansionOptions)` |
| `SetRunningOptions[opt]` | `set_running_options(RunningOptions)` |
| `SetWorkingPrecision[p]` | `working_prec_bits()`, `decimal_digits_to_bits(int)` |
| `AMFN[a]` | 隐式 (我们全程用任意精度) |
| `AMFChop[a]` | `acb_is_chop_zero(...)`, `AcbValue::is_chop_zero(...)` |
| `AMFPrint[...]` | `log_line(...)` |
| `AMFRationalize[a]` | 局部 helper (例如 `arb_to_rational`, `acb_real_to_fmpq`) |
| `SetDefaultOptions[]` | `set_default_options()` |

Block 作用域 `Block[{XOrder=..., SilentMode=True}, ...]` → `GlobalScope` RAII 类型 (commit() 装入, dtor 还原).

---

## Helper (lines 142-186)

| .m 符号 | C++ |
|---|---|
| `MatrixDensity[matrix]` | `RationalMatrix::density()` |
| `MatrixDegree[matrix]` | `RationalMatrix::max_numerator_degree()`, `max_denominator_degree()` |
| `ListPlotComplex` | — (诊断用) |
| `FitEps` | — (后处理工具) |
| `EvaluateExpansion[exp, x0]` | `evaluate_expansion(exp, x0)` |
| `EvaluateAsymptoticExpansion[mu -> exp, x0]` | `evaluate_asy_term(mu, exp, x0)` |
| `PickElement[list, order]` | `pick_element(list, order)` |
| `PickList[list, order]` | `pick_list(list, order)` |
| `PickMat[mat0, order]` | `pick_mat(mat, order)` |
| `TimesPoly[f, poly]` | `times_poly(f, poly)` |
| `InversePoly[f, poly]` | `inverse_poly(f, poly)` |
| `RationalExpansion[num, de, expansion]` | `rational_expansion(num, den, expansion)` |
| `MapRationalExpansion[nummap, demap, expansions]` | `map_rational_expansion(...)` |
| `ToStringInput[exp]` | — (仅 Mathematica) |

---

## 方程分析 (lines 189-304)

| .m 符号 | C++ |
|---|---|
| `AnalyzeBlockOld[mat]` | `analyze_block_old(mat)` |
| `AnalyzeBlock[mat]`, `AnalyzeBlock0[mat]` | `analyze_block(mat)` |
| `SubBlockID[mat]` | `sub_block_ids(mat, blocks)` |
| `UnionBehavior[behs]` | 局部 `union_behavior(...)` in `ode::asymptotic_behavior` |
| `UnionBehavior2[behs]` | 局部 `union_behavior2(...)` in `ode::asymptotic_behavior` |
| `AsymptoticBehavior[mat]` | `ode::asymptotic_behavior(mat, chop_digits, prec_bits)`; CalcZero 旋转 fallback: `ode::asymptotic_behavior_from_rotations(mat, rotations, chop_digits, prec_bits)` |
| `PoincareRank[mat]` | `poincare_rank(mat)` |
| `NHEquations[mat, mode]` | `nh_equations(mat, EquationMode)` (还有一个接受 BlockPartition 的 overload) |
| `ToNum[poly]` | `to_num(poly)` |
| `NHEquationsNum[nheqs]` | `nh_equations_num(eqs)` |

---

## 对角 block 的 normalization (lines 307-446)

| .m 符号 | C++ |
|---|---|
| `PBar[P]` | (内联; `I - P`) |
| `Balance[P]` | 局部 `balance(P)` in `src/ode/normalize.cpp`, 用于 ToFuchsian projector path |
| `InvBalance[P]` | 局部 `inv_balance(P)`, 用于 ToFuchsian projector path |
| `ReduceL0[L0, r, lblock]` | 局部 `reduce_l0_exact(...)`, 在 `fmpq_mat_t` 上 |
| `DynamicPartition[l, p]` | — (按需内联) |
| `FindProjector[cp, cp1]` | 局部 `find_projector_exact(...)`, 在 `fmpq_mat_t` 上 |
| `ToFuchsian[mat]` | 局部 `to_fuchsian_local(...)` in `src/ode/normalize.cpp`; 不可约的 leading rank 抛错 |
| `JordanDecomposition[(mat eta) /. eta -> 0]` | `jordan_decomposition_exact(...)` in `src/ode/jordan.cpp`, 用于 `NormalizeMat` rational 残数 |
| `NormalEigen[feps]` | 局部 `normal_eigen(...)` in `src/ode/inf.cpp`; Layer 7 整数 floor 在精确 Jordan 之后用精确 `fmpq_floor_si(...)` |
| `NormalizeEigenQ[mat]` | 局部 `normalize_eigen_q(...)` in `src/ode/normalize.cpp` |
| `ShearingTransformation[mat]` | 局部 `shearing_transformation(...)` in `src/ode/normalize.cpp` |
| `NormalizeEigen[mat]` | 局部 `normalize_eigen(...)` in `src/ode/normalize.cpp` |
| `NormalizeDiagonal[mat]` | 局部 `normalize_diagonal(...)` in `src/ode/normalize.cpp` |

---

## 非对角 block 的 normalization (lines 449-490)

| .m 符号 | C++ |
|---|---|
| `SolveOffDiagonal[a0, b0, c0, p]` | `internal::solve_off_diagonal_fmpq(...)` in `src/ode/normalize.cpp` |
| `ToFuchsianGlobal[mat]` | 局部 `to_fuchsian_global(...)` in `src/ode/normalize.cpp` |

---

## 归一化 Fuchsian 形式 (lines 493-515)

| .m 符号 | C++ |
|---|---|
| `NormalizeMat[mat]` | `normalize_mat(mat)` |

---

## 积分 contour (lines 519-586)

| .m 符号 | C++ |
|---|---|
| `AllFactors[exp_List]` | 在 `get_poles(...)` 内 |
| `Zeros[poly]` | 在 `get_poles(...)` 内 |
| `GetPoles[mat]` | `ode::singular_points(de, path_options)` |
| `FirstStep[polelist]` | `ode::first_step(poles, PathOptions)` |
| `LastStep[polelist]` | `ode::last_step(poles, PathOptions)` |
| `RunUnit[polelist]` | `ode::run_unit(poles, PathOptions)` |
| `RunSegment[polelist, ini, fin]` | `ode::run_segment(poles, ini, fin, PathOptions)` |
| `RunEtaDirection[polelist, 1]` | `run_eta_direction_positive(poles)` |
| `RunEtaDirection[polelist, direction]` | `run_eta_direction(poles, direction)` |
| `RunEtaDirection[polelist, mode_String]` | `ode::run_eta_direction(poles, PathDirectionMode, PathOptions)` |
| `RunEta[polelist]` | `ode::run_eta(poles, mode, PathOptions)` / `ode::PathPlan` |

---

## 稀疏求解 (lines 589-680)

| .m 符号 | C++ |
|---|---|
| `Sparsify[matrix]` | `sparsify_matrix(...)`, `sparsify_row(...)` |
| `SparseTranslate[sp, n]` | `sparse_translate(...)` |
| `ConstructMatrix[dx, ax, totalorder]` | `construct_matrix(...)` |
| `SplitSystem[sys, m]` | `split_system(...)` |
| `PlusSparse[sp1, sp2]` | `plus_sparse(...)` |
| `ScalarSparse[s, sp]` | `scalar_sparse(...)` |
| `ForwardSparseGaussian[sparse, vset]` | `forward_sparse_gaussian(...)` |
| `SparseGaussian[{dim, sparselist}, nh]` | `sparse_gaussian(...)` |
| `NSparse[sp]` | `chop_sparse(...)` |

---

## 边界条件 (lines 684-731)

| .m 符号 | C++ |
|---|---|
| `BuildTaylor[mat, ini]` | `ode::build_taylor_matrix(mat, initial_exponents)` |
| `DetermineBlockBoundaryOrder[mat, power]` | `determine_block_boundary_order(...)` |
| `DetermineBoundaryOrder[mat, power]` | `ode::determine_boundary_order(...)` |
| `ReverseBCS[bcs]` | `ode::reverse_boundary_exponents(...)` |
| `UnionBCS[bcs]` | `ode::union_boundary_exponents(...)` |
| `ReadBCS[bcs, region]` | `ode::read_boundary_region(...)` |

---

## 展开 (lines 734-795)

| .m 符号 | C++ |
|---|---|
| `CalcTaylor[mat, bc]` | `ode::calc_taylor(...)`, 用在 `BuildTaylor` 后的矩阵上 |
| `CalcInf[de, bcs]` | `ode::solve_at_infinity(...)` |

---

## 正则点 (lines 798-849)

| .m 符号 | C++ |
|---|---|
| `ExpandNHEquationsNum[nheqn, x0]` | `expand_nh_equations_num(nheqn, x0)` |
| `Calcx1x2[nheqn, bc, x0]` | `calcx1x2(...)` (返回 `TaylorCoefficients`) |
| `CalcRun[de, bc, run]` | `calc_run(...)` (接受 `RationalMatrix` 或 `vector<BlockEquationNum>` 的 overload) |

---

## 零点 (lines 852-917)

| .m 符号 | C++ |
|---|---|
| `Calcx00[nheq, nheqn, bc, x0, behavior]` | C++ 边界: `ZeroRecurrenceTable`, `ZeroCoefficientStore`, `build_zero_lower_sources(...)`, `build_zero_recurrence_block(...)`, `solve_zero_boundary_block(...)`, `ZeroBlockRotation`, `solve_zero_blocks(...)` |
| `FindLogPower[exp]` | `find_log_power(exp)` |
| `LearnFromRuleS[rules]` | `learn_from_rule_s(asy)` |
| `LearnFromRuleSAll[ruleslist, blocks]` | `learn_from_rule_s_all(asy_list, blocks)` |
| `ExtendExpansion[exp, n]` | `extend_expansion(exp, n)` |
| `PlusExpansion[exps]` | `plus_expansion(exps)` |
| `RescaleExpansion[exp, order]` | `rescale_expansion(exp, order)` |
| `PlusRuleS[rules0]` | `plus_rule_set(region)` |
| `UnionRuleS[rules]` | `union_rule_set(asy)` |
| `PSTimesRuleS[ps, rules]` | `ps_times_rule_set(ps, asy)` |
| `PSMapRuleS[psmap, ruleslist]` | `ps_map_rule_set(psmap, asy_list)` |
| `ToPS[poly]` | `to_power_series(rf)`, `to_power_series_matrix(mat)` |
| `CalcZero[de, bc, x0]` | C++ 路径: `ode::solve_at_zero(...)`; CalcZero 代数-Jordan fallback: `normalize_mat_for_calc_zero(...)`, `ZeroBlockRotation`, `asymptotic_behavior_from_rotations(...)`, `solve_zero_blocks(...)` 的 rotation overload |
| `TimesAsyExp[rational, asyexp]` | `times_asy_exp(rational, asy)` |
| `PlusAsyExp[asyexplist]` | `plus_asy_exp(list)` |

---

## 系统句柄 (lines 989-1080)

| .m 符号 | C++ |
|---|---|
| `LoadSystem[sysid, de, bc, point]` | C++ 用 value entry point 代替 mutable handle |
| `ClearSystem[sysid]` | `System` 析构 (RAII) |
| `DE[sysid]` | `System::de()` |
| `BC[sysid]` | `System::bc_singular()` / `System::bc_regular()` |
| `P[sysid]` | `System::point()` (或 `point_at_infinity()`) |
| `AsyExp[sysid]` | `System::asy_exp()` |
| `InfToRegular[sysid, x0]` | C++ 内部步骤, 在 `ode::solve_ode(...)` 里 |
| `RegularRun[sysid, run]` | `ode::solve_ode_from_regular(...)` |
| `RegularInterpolation[sysid, samples]` | `ode::regular_interpolation(de, boundary, point, samples, options)` |
| `SolveAsyExp[sysid]` | C++ 内部 `ode::solve_at_zero(...)` 步骤 |

---

## DESolver.m 里的 AMFlow (lines 1091-1138)

| .m 符号 | C++ |
|---|---|
| `PickZeroRuleS[rules]` | `ode::pick_zero_solution(asy, chop, prec)` |
| `$InternalSystem` | `amflow(...)` 内部局部的 `System` 实例 |
| `AMFlow[de, bc]` | `ode::solve_ode(de, boundaries, options)`, 带一个显式 poles overload |

---

# Part II — `AMFlow.m` (Layer 11–16)

.m 源里这部分引入 Feynman 积分概念 — family, `JIntegral`, top sector, region, eta 注入, SingleMass / Cutkosky ending scheme, 以及 `BlackBoxAMFlow` 顶层 driver. C++ 端重组成 Layer 11–16 (按 layer 的细节见 [`docs/ARCHITECTURE.md`](ARCHITECTURE.md)).

## 多变量代数原语

这些没有直接的 .m 对应 — Mathematica 通过内置符号引擎处理. C++ 项目要显式实现, 因为我们用 FLINT/Arb.

| 概念 (.m 习惯用法) | C++ |
|---|---|
| `{eta, eps, x[1], ..., x[k]}` 上的符号多项式 | `Mpoly` (Layer 11) |
| 同变量上的符号 rational | `Mfrac` (Layer 11) |
| `Together[...]`, 自动 gcd 消去 | 在 `Mfrac` 里隐式 (每个操作自动 reduce) |
| 多项式 / rational 矩阵的 `MatrixForm[...]` | `MpolyMatrix`, `MfracMatrix` (Layer 11) |
| `Det[m]`, `Inverse[m]`, `Cofactor[m]` | `MpolyMatrix::det()` / `MfracMatrix::inverse()` / `MpolyMatrix::adjugate()` |
| `m[[i,j]] /. var -> value` | `Mfrac::evaluate_to_acb({{name, value}, ...}, prec)` |
| `Variables[expr]` | `MpolyContext::variable_names()` |

## Family 与积分类型 (AMFlow.m, 文件顶端全局)

| .m 符号 | C++ |
|---|---|
| `AMFlowInfo[family]` (全局 association) | `FamilyConfig` 不可变对象 (Layer 12b) |
| `AMFlowInfo[family]["Loops"]` | `FamilyConfig::loops()` |
| `AMFlowInfo[family]["Externals"]` | `FamilyConfig::legs()` |
| `AMFlowInfo[family]["Conservation"]` | `FamilyConfig::conservation()` |
| `AMFlowInfo[family]["Replacement"]` | `FamilyConfig::replacement()` |
| `AMFlowInfo[family]["Propagators"]` | `FamilyConfig::propagators()` |
| `AMFlowInfo[family]["Cut"]` | `FamilyConfig::cut()` |
| `AMFlowInfo[family]["Prescription"]` | `FamilyConfig::prescription()` |
| `j[family, n_1, ..., n_k]` | `JIntegral { family, indices }` (Layer 12a) |
| `Sector[j[...]]` | `sector_of(JIntegral)` |
| `SortIntegrals[list]` | `sort_integrals(list)` |
| `GetTopSector[list]` | `get_top_sector(list)` |
| `SplitTarget[list]` | `split_target(list)` |
| `SymanzikU[family]`, `SymanzikF[family]` | `evaluate_uf(family, ctx) -> UFResult` (Layer 12c) |
| `ABCofU[a, family]` | `evaluate_abc(family, alpha)` (Layer 12c) |
| `AnalyzeTopology[family]` | `analyze_topology(family)` (Layer 12d) |
| `ZeroSectorQ[family, sector]` | `zero_sector_q(family, sector)` (Layer 12d) |

## Eta 注入 (AMFlow.m: AMFCandidate / AMFPosition / AMFEtaC)

| .m 符号 | C++ |
|---|---|
| `AMFCandidate[family, mode]` | `amf_candidate(family, mode)` (Layer 13) |
| `AMFPosition[family, mode]` | `amf_position(family, mode)` |
| `AMFEtaC[family, position]` | `amf_eta_c(family, position)` |
| `AllPossiblePosition[family, mode]` | `all_possible_position(family, mode)` |
| `AnalyzeTopSector[family, jint]` | `analyze_top_sector(family, jint)` |
| Mode 字符串 `"Mass"`, `"Propagator"`, `"Branch"`, `"Loop"`, `"All"` | `enum class AMFMode { Prescription, Mass, Propagator, Branch, Loop, All }` |
| `EndingQ[component]` | `ending_q(component)` (还有 `vacuum_q`, `single_mass_q`, `phase_volume_q`) |
| `TopSectorComponent[id, props, loops, masses]` | `TopSectorComponentInfo { component_id, propagator_indices, loop_subspace, mass_count, ... }` |

## Region decomposition 与 boundary integrand (AMFlow.m: 中段)

| .m 符号 | C++ |
|---|---|
| `BranchMomenta[family]` | `branch_momenta(family)` (Layer 14a) |
| `BranchToLoop[family]` | `branch_to_loop(family)` |
| `RegionRule[family, mode]` | `region_rule(family, mode)` |
| `FindAllRegion[family, mode]` | `find_all_region(family, mode)` |
| `ZeroRegionQ[family, region]` | `zero_region_q(family, region)` |
| `RegionPower[family, jint, region, eps]` | `region_power(family, jint, region, eps)` |
| `LoopTransform[matrix]` | `LoopTransform { matrix }` (struct in `region.hpp`) |
| `SquaredDenominators[family]` | `squared_denominators(family)` (Layer 14a) |
| `ToSquareAll[family, expr]` | `to_square_all(family, expr)` |
| `ToCompleteExplicit[family, expr]` | `to_complete_explicit(family, expr)` |
| `SPListToDListSymbol[family, sp_list]` | `sp_list_to_dlist_symbol(family, sp_list)` |
| 符号 `D[i]` 分母基 | `DListContext` (Layer 14a) |
| `ApartOneVar[r, var]` | `apart_one_var(rational, var)` (Layer 14e) |
| `ApartRationals[list, vars]` | `apart_rationals(list, vars)` |
| `BoundaryPattern[family, region]` | `qft::boundary_power_patterns(...)` |
| `BoundaryIntegrands[family, region]` | `boundary_integrands(family, region)` |
| `LaportaIntegrals[integrand]` | `amflow::laporta_boundary_terms(...)` / `analyze_boundary_pfd_term(...)` / `make_boundary_integrand_terms_from_pfd(...)` |
| `BoundaryIntegrals[family, region, jint_list]` | `boundary_integrals(family, region, jint_list) -> vector<BoundaryTerm>` |
| `Vacuum[L, n, eps]` (5-entry 表) | `vacuum(L, n, eps)` (Layer 14f) |

## IBP 导数辅助 (AMFlow.m + DESolver.m bridge)

| .m 符号 | C++ |
|---|---|
| `LIBPDeriv[jint, family, var]` | `libp_deriv(jint, family, var)` (Layer 15a) |
| `ComputeDerivative[jint, family, var]` | `compute_derivative(jint, family, var)` |
| `DenomsDeriv[...]` (内部辅助) | `LibpDenomsDerivResult` 容器 |

## SingleMass / Cutkosky / FactorizeFamily (AMFlow.m: ending scheme)

| .m 符号 | C++ |
|---|---|
| `FactorizeFamily[family, preferred]` | `factorize_family(family, preferred)` (Layer 16b) |
| `FactorizedComponent[...]` | `FactorizedComponent { sub_family, sub_preferred, removed_propagator, loop_redefinition_matrix, gamma_prefactor }` |
| `RREFTransform[matrix]` (loop 重定义) | 在 `factorize_family(...)` 内 |
| `SingleMassEnding[family, ...]` | `AMFSystem::single_mass_setup(...)` 内的分支 (Layer 16) |
| `CutkoskyEnding[family, ...]` | `amf_system_setup_master(...)` 里的 `EndingScheme::Cutkosky` 分支; 清空 cut/prescription, 应用 phase-space prefactor 与 `Im[...]`; 含 `CutkoskyPrefactor` + `evaluate_cutkosky_prefactor` + `apply_cutkosky_prefactor` + `cutkosky_prefactor_for_request` + `apply_cutkosky_top_level_result` |
| Γ-函数 prefactor 相乘 | 在 `AMFSystem::solve(...)` 里组合 `master_values` → `global_values` 时应用; 含 `SingleMassPrefactor` + `evaluate_single_mass_prefactor` + `apply_single_mass_child_prefactors` / `apply_single_mass_prefactors` + `run_single_mass_execution_steps` |

## AMFSystem driver (AMFlow.m: AMFSystem / AMFSystemsSetup / AMFSystemsSolution)

| .m 符号 | C++ |
|---|---|
| `AMFSystem[family, preferred, opts]` | `AMFSystem` 类 (Layer 16) |
| `AMFSystemSetup[...]` | `AMFSystem::setup()`; 经 eta 注入 DE setup `make_eta_injected_kinematics` + `make_differential_equation_request` |
| `AMFSystemsSetup[...]` | `amf_systems_setup(family, preferred, opts)` (自由函数) |
| `AMFSystemSolution[..., epslist]` | `AMFSystem::solve(epslist)` (request/result 类型见 `include/amflow/pipeline/amfsystem.hpp`) |
| `AMFSystemsSolution[..., epslist]` | `amf_systems_solution(root, epslist)` (装配 / cache / API 表面见 `include/amflow/pipeline/amfsystem.hpp` 与 `include/amflow/api/`) |
| `RegionBoundary[region, mu, terms]` | `RegionBoundary { region, mu, terms }` |
| `BoundaryTerm[sub_index, coefficient]` | `BoundaryMasterTerm` + `BoundaryExpansionTable` + `run_reduction_step` + `run_boundary_expansion_step` + `make_boundary_evaluation_point` + `BoundaryContributionTable` + `run_boundary_contribution_step` + `BoundaryVectorProjection` handoff |
| `EndingMode -> "Tradition"` | `EndingScheme::Tradition` |
| `EndingMode -> "SingleMass"` | `EndingScheme::SingleMass` |
| `EndingMode -> "Cutkosky"` | `EndingScheme::Cutkosky` |
| `MasterValues[...]` | `AMFSystemSolution::master_values`; 普通 Tradition root 上 `AmflowResult::values()` |
| `GlobalValues[...]` (post-prefactor) | `AMFSystemSolution::global_values`; Cutkosky / SingleMass 装配后的 `AmflowResult::values()` |

## 顶层 driver (AMFlow.m: BlackBoxAMFlow / SolveIntegrals / FitEps)

| .m 符号 | C++ 状态 |
|---|---|
| `BlackBoxAMFlow[family, jints, epslist, dir]` | `black_box_amflow(fc, jints, eps_samples, opts, work_dir)`; 经采样 facade `api::try_black_box_amflow(request)` 加 CLI adapter `cli::parse_black_box_amflow_json_payload(...)`, `cli::try_run_black_box_amflow_json(...)`, `cli::run_black_box_amflow_json(...)`, `cli::run_refactor_cli_json(...)`, `cli::black_box_amflow_output_to_json(...)` |
| `BlackBoxAMFlowSingle[family, jints, epslist, dir]` | `black_box_amflow_single(fc, jints, eps_samples, opts, work_dir)` |
| `GenerateNumericalConfig[goal, order]` | `generate_numerical_config(loop_count, goal_digits, eps_order)`; `generate_ordinary_numerical_config(loop_count, goal_digits, eps_order)` |
| `FitEps[values, epslist, order]` | `fit_eps(eps_samples, values, leading_order)`; `fit_laurent_coefficients(eps_samples, values, leading_order, precision_bits, chop_digits)` |
| `SolveIntegrals[jints, goal, order]` | `solve_integrals(fc, jints, goal_digits, eps_order, opts, work_dir)`; 经普通 boundary `ordinary_solve_integrals(request)` over `solve_ordinary_integrals_laurent(request)`; API facade `api::solve_integrals(request)`; CLI adapter `cli::parse_solve_integrals_json_payload(...)`, `cli::parse_kinematic_family_json(...)`, `cli::plan_solve_integrals_subtree(...)`, `cli::plan_solve_integrals_subtree_shallow(...)`, `cli::plan_solve_integrals_json_child_specs(...)`, `cli::try_run_planned_solve_integrals(...)`, `cli::try_run_solve_integrals_json(...)`, `cli::run_solve_integrals_json(...)`, `cli::run_refactor_cli_json(...)`, `cli::solve_integrals_output_to_json(...)` |
| `"D0" -> 4` (`Options[SetAMFOptions]`, AMFlow.m:262) | `GlobalOptions::d0` (rational 字符串, 默认 `"4"`) |
| `epslist = epslist0 + (4-$D0)*1/2` (AMFlow.m:1342 / 1351) | `d0_eps_shift_fmpq(fmpq_t)` + `solve_integrals` 内部的 shift; `shift_epsilon_samples_for_dimension(...)` |
| `GenerateSquare[prop, symbol]` | **未 port** (gauge-link) |
| `SolveIntegralsGaugeLink[jints, goal, order]` | **未 port** — HQET / SCET / Wilson lines 对本 port out of scope |

这些在 C++ 端暴露在 Layer 17. `amflow_cli` 覆盖 raw `amflow`, 采样 `black_box_amflow`, `solve_integrals` 三种 JSON mode. `src/cli/` 下的 CLI 树提供 solve-integrals JSON payload/result, kinematic-family 拼接, 计划执行, 顶层 JSON runner, JSON 子树 spec adapter; `amflow_cli` 是该 runner 之上的独立二进制 adapter. JSON runner 把 `amf_options.ending_schemes`, `amf_options.max_recursion_depth`, 以及支持的 `amf_options.amf_modes` / `eta_modes` 在 plan 之前映射到 `amflow::AmflowPlanOptions`.

---

# Part III — `ibp_interface/Kira/interface.m` 与辅助

这是包裹 Kira 子进程调用的小 Mathematica 模块. C++ 端对应的在 Layer 15 (`kira.hpp`, `blackbox.hpp`).

## YAML 与 target 写入

| .m 符号 | C++ |
|---|---|
| `KiraConfig` (association) | `KiraConfig` struct (Layer 15) |
| `WriteIntegralfamilies[family, cfg]` | `kira_write_integralfamilies(family, cfg)` |
| `WriteKinematics[family, cfg]` | `kira_write_kinematics(family, cfg)` |
| `WriteJobs[cfg, "Masters"]` | `kira_write_jobs(cfg, KiraJobMode::Masters)` |
| `WriteJobs[cfg, "Reduce"]` | `kira_write_jobs(cfg, KiraJobMode::Reduce)` |
| `WritePreferred[masters, cfg]` | `kira_write_preferred(masters, cfg)` |
| `WriteTarget[targets, cfg]` | `kira_write_target(targets, cfg)` |

## 子进程与 parser

| .m 符号 | C++ |
|---|---|
| `RunKira[cfg]` (一般是 `Run["kira jobs.yaml"]`) | `kira_run(cfg)` (fork / execve / wait) |
| `ReadMasters[cfg]` | `kira_read_masters(cfg)` |
| `ReadTargetTable[cfg]` (parse `kira_target.m`) | `kira_read_target_table(cfg)` |
| Mathematica `ToExpression["Plus[Times[...], j[...]]"]` | `kira_parse_expression(str, dlist_ctx)` |
| `kira_target.m` 里的大整数 | 经 `fmpz_set_str` 解析 |

## 高层约化包装

| .m 符号 | C++ |
|---|---|
| `BlackBoxReduce[family, targets, cfg]` | `ibp::black_box_reduce` (Layer 15d) 加 `amflow::run_reduction_step`, 用于 boundary-task 执行 |
| `BlackBoxDiffeq[family, masters, var, cfg]` | `BlackBoxDiffeq(family, masters, var, cfg) -> BlackBoxDiffeqResult` |
| `ReductionContext[family, masters, targets]` | `ReductionContext` struct |
| `BlackBoxOptions` (cfg flag 的 association) | `BlackBoxOptions` struct (传给 KiraConfig) |
| `$KiraDumpIO` flag | env `AMFLOW_KIRA_DUMP=1` (与 `BlackBoxOptions::dump_kira_io`) |
| `$KiraPath`, `$FermatPath` | `KiraConfig::kira_path`, `KiraConfig::fermat_path` |
