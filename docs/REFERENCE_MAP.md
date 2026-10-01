# Mathematica ↔ C++ Cross-Reference

> 中文版本: [REFERENCE_MAP_zh.md](REFERENCE_MAP_zh.md)

This document maps the symbols of the three upstream Mathematica
source files to their C++ counterparts.  Upstream lives at
<https://gitlab.com/multiloop-pku/amflow>; this project does not
vendor it (see [`reference/README.md`](../reference/README.md)).

- Upstream `diffeq_solver/DESolver.m` (the ODE engine, Layers 0–10).
  Covered in Part I below.
- Upstream `AMFlow.m` (the Feynman-integral pipeline, Layers 11–16).
  Covered in Part II below.
- Upstream `ibp_interface/Kira/interface.m` and the small Kira-glue
  helpers.  Covered in Part III below.

Line numbers cited below track upstream commit
[`efda1db`](https://gitlab.com/multiloop-pku/amflow/-/commit/efda1db)
("add_version_number_1.2"), the snapshot taken when this port was
performed.  Upstream may have moved on since (see
[`AUDIT_MMA_PARITY.md`](AUDIT_MMA_PARITY.md) §5 for any drift), so
individual line numbers may be off in newer revisions.
**Search by symbol name when in doubt.**  In-source `// Mirrors
AMFlow.m:NNN` comments follow the same convention.

Use this map as a navigation aid when you need to translate a
Mathematica expression or trace why a result differs.

Conventions:

- "Layer N" labels below are **legacy module names** from the v1.0
  port.  The current source tree is organized by **domain** under
  `src/<domain>/` (`numeric`, `algebra`, `ode`, `qft`, `ibp`,
  `pipeline`, `api`, `cli`); see
  [`docs/ARCHITECTURE.md`](ARCHITECTURE.md) for the current
  organization.  Layer 6/7 ↔ `src/ode/`, Layer 11-16 ↔
  `src/qft/` + `src/ibp/` + `src/pipeline/`, Layer 17 ↔ `src/cli/`.
- "—" means the function has no C++ counterpart, either because it is
  a Mathematica idiom (e.g. `MatrixDensity` for diagnostics) or because
  it is an internal helper that has been refactored.

---

# Part I — `DESolver.m` (Layers 0–10)

## Top of `DESolver.m` — Options block (lines 62-132)

| .m symbol | C++ |
|---|---|
| `WorkingPre`, `ChopPre`, `SilentMode`, `RationalizePre` | `GlobalOptions` fields |
| `XOrder`, `ExtraXOrder`, `LearnXOrder`, `TestXOrder` | `ExpansionOptions` fields |
| `RunRadius`, `RunLength`, `RunCandidate`, `RunDirection` | `RunningOptions` fields |
| `SetGlobalOptions[opt]` | `set_global_options(GlobalOptions)` |
| `SetExpansionOptions[opt]` | `set_expansion_options(ExpansionOptions)` |
| `SetRunningOptions[opt]` | `set_running_options(RunningOptions)` |
| `SetWorkingPrecision[p]` | `working_prec_bits()`, `decimal_digits_to_bits(int)` |
| `AMFN[a]` | implicit (we work in arbitrary precision throughout) |
| `AMFChop[a]` | `acb_is_chop_zero(...)`, `AcbValue::is_chop_zero(...)` |
| `AMFPrint[...]` | `log_line(...)` |
| `AMFRationalize[a]` | local helpers (e.g. `arb_to_rational`, `acb_real_to_fmpq`) |
| `SetDefaultOptions[]` | `set_default_options()` |

Block scope `Block[{XOrder=..., SilentMode=True}, ...]` →
`GlobalScope` RAII type (commit() to install, dtor restores).

---

## Helpers (lines 142-186)

| .m symbol | C++ |
|---|---|
| `MatrixDensity[matrix]` | `RationalMatrix::density()` |
| `MatrixDegree[matrix]` | `RationalMatrix::max_numerator_degree()`, `max_denominator_degree()` |
| `ListPlotComplex` | — (diagnostic) |
| `FitEps` | — (post-processing utility) |
| `EvaluateExpansion[exp, x0]` | `evaluate_expansion(exp, x0)` |
| `EvaluateAsymptoticExpansion[mu -> exp, x0]` | `evaluate_asy_term(mu, exp, x0)` |
| `PickElement[list, order]` | `pick_element(list, order)` |
| `PickList[list, order]` | `pick_list(list, order)` |
| `PickMat[mat0, order]` | `pick_mat(mat, order)` |
| `TimesPoly[f, poly]` | `times_poly(f, poly)` |
| `InversePoly[f, poly]` | `inverse_poly(f, poly)` |
| `RationalExpansion[num, de, expansion]` | `rational_expansion(num, den, expansion)` |
| `MapRationalExpansion[nummap, demap, expansions]` | `map_rational_expansion(...)` |
| `ToStringInput[exp]` | — (Mathematica-only) |

---

## Equation analysis (lines 189-304)

| .m symbol | C++ |
|---|---|
| `AnalyzeBlockOld[mat]` | `analyze_block_old(mat)` |
| `AnalyzeBlock[mat]`, `AnalyzeBlock0[mat]` | `analyze_block(mat)` |
| `SubBlockID[mat]` | `sub_block_ids(mat, blocks)` |
| `UnionBehavior[behs]` | local `union_behavior(...)` in `ode::asymptotic_behavior` |
| `UnionBehavior2[behs]` | local `union_behavior2(...)` in `ode::asymptotic_behavior` |
| `AsymptoticBehavior[mat]` | `ode::asymptotic_behavior(mat, chop_digits, prec_bits)`; CalcZero rotation fallback: `ode::asymptotic_behavior_from_rotations(mat, rotations, chop_digits, prec_bits)` |
| `PoincareRank[mat]` | `poincare_rank(mat)` |
| `NHEquations[mat, mode]` | `nh_equations(mat, EquationMode)` (also overload taking BlockPartition) |
| `ToNum[poly]` | `to_num(poly)` |
| `NHEquationsNum[nheqs]` | `nh_equations_num(eqs)` |

---

## Diagonal blocks normalisation (lines 307-446)

| .m symbol | C++ |
|---|---|
| `PBar[P]` | (inlined; `I - P`) |
| `Balance[P]` | local `balance(P)` in `src/ode/normalize.cpp` for the ToFuchsian projector path |
| `InvBalance[P]` | local `inv_balance(P)` for the ToFuchsian projector path |
| `ReduceL0[L0, r, lblock]` | local `reduce_l0_exact(...)` over `fmpq_mat_t` |
| `DynamicPartition[l, p]` | — (inlined where needed) |
| `FindProjector[cp, cp1]` | local `find_projector_exact(...)` over `fmpq_mat_t` |
| `ToFuchsian[mat]` | local `to_fuchsian_local(...)` in `src/ode/normalize.cpp`; irreducible leading rank throws |
| `JordanDecomposition[(mat eta) /. eta -> 0]` | `jordan_decomposition_exact(...)` in `src/ode/jordan.cpp` for `NormalizeMat` rational residues |
| `NormalEigen[feps]` | local `normal_eigen(...)` in `src/ode/inf.cpp`; Layer 7 integer floors use exact `fmpq_floor_si(...)` after exact Jordan |
| `NormalizeEigenQ[mat]` | local `normalize_eigen_q(...)` in `src/ode/normalize.cpp` |
| `ShearingTransformation[mat]` | local `shearing_transformation(...)` in `src/ode/normalize.cpp` |
| `NormalizeEigen[mat]` | local `normalize_eigen(...)` in `src/ode/normalize.cpp` |
| `NormalizeDiagonal[mat]` | local `normalize_diagonal(...)` in `src/ode/normalize.cpp` |

---

## Off-diagonal blocks normalisation (lines 449-490)

| .m symbol | C++ |
|---|---|
| `SolveOffDiagonal[a0, b0, c0, p]` | `internal::solve_off_diagonal_fmpq(...)` in `src/ode/normalize.cpp` |
| `ToFuchsianGlobal[mat]` | local `to_fuchsian_global(...)` in `src/ode/normalize.cpp` |

---

## Normalised Fuchsian form (lines 493-515)

| .m symbol | C++ |
|---|---|
| `NormalizeMat[mat]` | `normalize_mat(mat)` |

---

## Integration contour (lines 519-586)

| .m symbol | C++ |
|---|---|
| `AllFactors[exp_List]` | inside `get_poles(...)` |
| `Zeros[poly]` | inside `get_poles(...)` |
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

## Sparse solve (lines 589-680)

| .m symbol | C++ |
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

## Boundary conditions (lines 684-731)

| .m symbol | C++ |
|---|---|
| `BuildTaylor[mat, ini]` | `ode::build_taylor_matrix(mat, initial_exponents)` |
| `DetermineBlockBoundaryOrder[mat, power]` | `determine_block_boundary_order(...)` |
| `DetermineBoundaryOrder[mat, power]` | `ode::determine_boundary_order(...)` |
| `ReverseBCS[bcs]` | `ode::reverse_boundary_exponents(...)` |
| `UnionBCS[bcs]` | `ode::union_boundary_exponents(...)` |
| `ReadBCS[bcs, region]` | `ode::read_boundary_region(...)` |

---

## Expand (lines 734-795)

| .m symbol | C++ |
|---|---|
| `CalcTaylor[mat, bc]` | `ode::calc_taylor(...)` for post-`BuildTaylor` matrices |
| `CalcInf[de, bcs]` | `ode::solve_at_infinity(...)` |

---

## Regular point (lines 798-849)

| .m symbol | C++ |
|---|---|
| `ExpandNHEquationsNum[nheqn, x0]` | `expand_nh_equations_num(nheqn, x0)` |
| `Calcx1x2[nheqn, bc, x0]` | `calcx1x2(...)` (returns `TaylorCoefficients`) |
| `CalcRun[de, bc, run]` | `calc_run(...)` (overloads taking `RationalMatrix` or `vector<BlockEquationNum>`) |

---

## Zero (lines 852-917)

| .m symbol | C++ |
|---|---|
| `Calcx00[nheq, nheqn, bc, x0, behavior]` | C++ boundary: `ZeroRecurrenceTable`, `ZeroCoefficientStore`, `build_zero_lower_sources(...)`, `build_zero_recurrence_block(...)`, `solve_zero_boundary_block(...)`, `ZeroBlockRotation`, `solve_zero_blocks(...)` |
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
| `CalcZero[de, bc, x0]` | C++ path: `ode::solve_at_zero(...)`; CalcZero algebraic-Jordan fallback: `normalize_mat_for_calc_zero(...)`, `ZeroBlockRotation`, `asymptotic_behavior_from_rotations(...)`, rotation overload of `solve_zero_blocks(...)` |
| `TimesAsyExp[rational, asyexp]` | `times_asy_exp(rational, asy)` |
| `PlusAsyExp[asyexplist]` | `plus_asy_exp(list)` |

---

## System handle (lines 989-1080)

| .m symbol | C++ |
|---|---|
| `LoadSystem[sysid, de, bc, point]` | C++ uses value entry points instead of a mutable handle |
| `ClearSystem[sysid]` | `System` destructor (RAII) |
| `DE[sysid]` | `System::de()` |
| `BC[sysid]` | `System::bc_singular()` / `System::bc_regular()` |
| `P[sysid]` | `System::point()` (or `point_at_infinity()`) |
| `AsyExp[sysid]` | `System::asy_exp()` |
| `InfToRegular[sysid, x0]` | C++ internal step in `ode::solve_ode(...)` |
| `RegularRun[sysid, run]` | `ode::solve_ode_from_regular(...)` |
| `RegularInterpolation[sysid, samples]` | `ode::regular_interpolation(de, boundary, point, samples, options)` |
| `SolveAsyExp[sysid]` | C++ internal `ode::solve_at_zero(...)` step |

---

## AMFlow inside DESolver.m (lines 1091-1138)

| .m symbol | C++ |
|---|---|
| `PickZeroRuleS[rules]` | `ode::pick_zero_solution(asy, chop, prec)` |
| `$InternalSystem` | local `System` instance inside `amflow(...)` |
| `AMFlow[de, bc]` | `ode::solve_ode(de, boundaries, options)` with an explicit-poles overload |

---

# Part II — `AMFlow.m` (Layers 11–16)

This part of the .m source introduces the Feynman-integral concepts —
families, `JIntegral`s, top sectors, regions, eta injection,
SingleMass / Cutkosky ending schemes, and the `BlackBoxAMFlow`
top-level driver.  The C++ side reorganises this into Layers 11–16
(see [`docs/ARCHITECTURE.md`](ARCHITECTURE.md) for per-layer detail).

## Multivariate algebra primitives

These have no direct .m counterpart — Mathematica handles them via the
built-in symbolic engine.  The C++ project must implement them
explicitly because we use FLINT/Arb.

| Concept (.m idiom) | C++ |
|---|---|
| Symbolic polynomial in `{eta, eps, x[1], ..., x[k]}` | `Mpoly` (Layer 11) |
| Symbolic rational in same variables | `Mfrac` (Layer 11) |
| `Together[...]`, automatic gcd cancellation | implicit in `Mfrac` (every op auto-reduces) |
| `MatrixForm[...]` of polynomial / rational matrices | `MpolyMatrix`, `MfracMatrix` (Layer 11) |
| `Det[m]`, `Inverse[m]`, `Cofactor[m]` | `MpolyMatrix::det()` / `MfracMatrix::inverse()` / `MpolyMatrix::adjugate()` |
| `m[[i,j]] /. var -> value` | `Mfrac::evaluate_to_acb({{name, value}, ...}, prec)` |
| `Variables[expr]` | `MpolyContext::variable_names()` |

## Family and integral types (AMFlow.m, top-of-file globals)

| .m symbol | C++ |
|---|---|
| `AMFlowInfo[family]` (global association) | `FamilyConfig` immutable object (Layer 12b) |
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

## Eta-injection (AMFlow.m: AMFCandidate / AMFPosition / AMFEtaC)

| .m symbol | C++ |
|---|---|
| `AMFCandidate[family, mode]` | `amf_candidate(family, mode)` (Layer 13) |
| `AMFPosition[family, mode]` | `amf_position(family, mode)` |
| `AMFEtaC[family, position]` | `amf_eta_c(family, position)` |
| `AllPossiblePosition[family, mode]` | `all_possible_position(family, mode)` |
| `AnalyzeTopSector[family, jint]` | `analyze_top_sector(family, jint)` |
| Mode strings `"Mass"`, `"Propagator"`, `"Branch"`, `"Loop"`, `"All"` | `enum class AMFMode { Prescription, Mass, Propagator, Branch, Loop, All }` |
| `EndingQ[component]` | `ending_q(component)` (also `vacuum_q`, `single_mass_q`, `phase_volume_q`) |
| `TopSectorComponent[id, props, loops, masses]` | `TopSectorComponentInfo { component_id, propagator_indices, loop_subspace, mass_count, ... }` |

## Region decomposition and boundary integrands (AMFlow.m: middle section)

| .m symbol | C++ |
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
| Symbolic `D[i]` denominator basis | `DListContext` (Layer 14a) |
| `ApartOneVar[r, var]` | `apart_one_var(rational, var)` (Layer 14e) |
| `ApartRationals[list, vars]` | `apart_rationals(list, vars)` |
| `BoundaryPattern[family, region]` | `qft::boundary_power_patterns(...)` |
| `BoundaryIntegrands[family, region]` | `boundary_integrands(family, region)` |
| `LaportaIntegrals[integrand]` | `amflow::laporta_boundary_terms(...)` / `analyze_boundary_pfd_term(...)` / `make_boundary_integrand_terms_from_pfd(...)` |
| `BoundaryIntegrals[family, region, jint_list]` | `boundary_integrals(family, region, jint_list) -> vector<BoundaryTerm>` |
| `Vacuum[L, n, eps]` (5-entry table) | `vacuum(L, n, eps)` (Layer 14f) |

## IBP-derivative helpers (AMFlow.m + DESolver.m bridge)

| .m symbol | C++ |
|---|---|
| `LIBPDeriv[jint, family, var]` | `libp_deriv(jint, family, var)` (Layer 15a) |
| `ComputeDerivative[jint, family, var]` | `compute_derivative(jint, family, var)` |
| `DenomsDeriv[...]` (internal helper) | `LibpDenomsDerivResult` container |

## SingleMass / Cutkosky / FactorizeFamily (AMFlow.m: ending schemes)

| .m symbol | C++ |
|---|---|
| `FactorizeFamily[family, preferred]` | `factorize_family(family, preferred)` (Layer 16b) |
| `FactorizedComponent[...]` | `FactorizedComponent { sub_family, sub_preferred, removed_propagator, loop_redefinition_matrix, gamma_prefactor }` |
| `RREFTransform[matrix]` (loop redefinition) | inside `factorize_family(...)` |
| `SingleMassEnding[family, ...]` | branch in `AMFSystem::single_mass_setup(...)` (Layer 16) |
| `CutkoskyEnding[family, ...]` | `EndingScheme::Cutkosky` branch in `amf_system_setup_master(...)`; clears cut/prescription, applies phase-space prefactor and `Im[...]`; including `CutkoskyPrefactor` + `evaluate_cutkosky_prefactor` + `apply_cutkosky_prefactor` + `cutkosky_prefactor_for_request` + `apply_cutkosky_top_level_result` |
| Γ-function prefactor multiplications | applied inside `AMFSystem::solve(...)` when combining `master_values` → `global_values`; including `SingleMassPrefactor` + `evaluate_single_mass_prefactor` + `apply_single_mass_child_prefactors` / `apply_single_mass_prefactors` + `run_single_mass_execution_steps` |

## AMFSystem driver (AMFlow.m: AMFSystem / AMFSystemsSetup / AMFSystemsSolution)

| .m symbol | C++ |
|---|---|
| `AMFSystem[family, preferred, opts]` | `AMFSystem` class (Layer 16) |
| `AMFSystemSetup[...]` | `AMFSystem::setup()`; via eta-injected DE setup `make_eta_injected_kinematics` + `make_differential_equation_request` |
| `AMFSystemsSetup[...]` | `amf_systems_setup(family, preferred, opts)` (free function) |
| `AMFSystemSolution[..., epslist]` | `AMFSystem::solve(epslist)` (see `include/amflow/pipeline/amfsystem.hpp` for the request/result types) |
| `AMFSystemsSolution[..., epslist]` | `amf_systems_solution(root, epslist)` (see `include/amflow/pipeline/amfsystem.hpp` and `include/amflow/api/` for assembly / cache / API surface) |
| `RegionBoundary[region, mu, terms]` | `RegionBoundary { region, mu, terms }` |
| `BoundaryTerm[sub_index, coefficient]` | `BoundaryMasterTerm` + `BoundaryExpansionTable` + `run_reduction_step` + `run_boundary_expansion_step` + `make_boundary_evaluation_point` + `BoundaryContributionTable` + `run_boundary_contribution_step` + `BoundaryVectorProjection` handoff |
| `EndingMode -> "Tradition"` | `EndingScheme::Tradition` |
| `EndingMode -> "SingleMass"` | `EndingScheme::SingleMass` |
| `EndingMode -> "Cutkosky"` | `EndingScheme::Cutkosky` |
| `MasterValues[...]` | `AMFSystemSolution::master_values`; `AmflowResult::values()` on ordinary Tradition roots |
| `GlobalValues[...]` (post-prefactor) | `AMFSystemSolution::global_values`; `AmflowResult::values()` after Cutkosky / SingleMass assembly |

## Top-level driver (AMFlow.m: BlackBoxAMFlow / SolveIntegrals / FitEps)

| .m symbol | C++ status |
|---|---|
| `BlackBoxAMFlow[family, jints, epslist, dir]` | `black_box_amflow(fc, jints, eps_samples, opts, work_dir)`; via sampled facade `api::try_black_box_amflow(request)` plus CLI adapters `cli::parse_black_box_amflow_json_payload(...)`, `cli::try_run_black_box_amflow_json(...)`, `cli::run_black_box_amflow_json(...)`, `cli::run_refactor_cli_json(...)`, and `cli::black_box_amflow_output_to_json(...)` |
| `BlackBoxAMFlowSingle[family, jints, epslist, dir]` | `black_box_amflow_single(fc, jints, eps_samples, opts, work_dir)` |
| `GenerateNumericalConfig[goal, order]` | `generate_numerical_config(loop_count, goal_digits, eps_order)`; `generate_ordinary_numerical_config(loop_count, goal_digits, eps_order)` |
| `FitEps[values, epslist, order]` | `fit_eps(eps_samples, values, leading_order)`; `fit_laurent_coefficients(eps_samples, values, leading_order, precision_bits, chop_digits)` |
| `SolveIntegrals[jints, goal, order]` | `solve_integrals(fc, jints, goal_digits, eps_order, opts, work_dir)`; via ordinary boundary `ordinary_solve_integrals(request)` over `solve_ordinary_integrals_laurent(request)`; API facade `api::solve_integrals(request)`; CLI adapters `cli::parse_solve_integrals_json_payload(...)`, `cli::parse_kinematic_family_json(...)`, `cli::plan_solve_integrals_subtree(...)`, `cli::plan_solve_integrals_subtree_shallow(...)`, `cli::plan_solve_integrals_json_child_specs(...)`, `cli::try_run_planned_solve_integrals(...)`, `cli::try_run_solve_integrals_json(...)`, `cli::run_solve_integrals_json(...)`, `cli::run_refactor_cli_json(...)`, and `cli::solve_integrals_output_to_json(...)` |
| `"D0" -> 4` (`Options[SetAMFOptions]`, AMFlow.m:262) | `GlobalOptions::d0` (rational string, default `"4"`) |
| `epslist = epslist0 + (4-$D0)*1/2` (AMFlow.m:1342 / 1351) | `d0_eps_shift_fmpq(fmpq_t)` + internal shift inside `solve_integrals`; `shift_epsilon_samples_for_dimension(...)` |
| `GenerateSquare[prop, symbol]` | **not ported** (gauge-link) |
| `SolveIntegralsGaugeLink[jints, goal, order]` | **not ported** — HQET / SCET / Wilson lines are out of scope for this port |

These are exposed on the C++ side in Layer 17.  `amflow_cli` covers
raw `amflow`, sampled `black_box_amflow`, and `solve_integrals` JSON
modes.  The CLI tree under `src/cli/` provides solve-integrals JSON
payload/result, kinematic-family sourcing, planned execution, a
top-level JSON runner, and JSON child-subtree spec adapters, with
`amflow_cli` being the standalone binary adapter over that runner.
The JSON runner maps `amf_options.ending_schemes`,
`amf_options.max_recursion_depth`, and the supported
`amf_options.amf_modes` / `eta_modes` into
`amflow::AmflowPlanOptions` before planning.

---

# Part III — `ibp_interface/Kira/interface.m` and helpers

This is the small Mathematica module that wraps Kira sub-process
invocation.  In C++ the equivalent lives in Layer 15 (`kira.hpp`,
`blackbox.hpp`).

## YAML and target writers

| .m symbol | C++ |
|---|---|
| `KiraConfig` (association) | `KiraConfig` struct (Layer 15) |
| `WriteIntegralfamilies[family, cfg]` | `kira_write_integralfamilies(family, cfg)` |
| `WriteKinematics[family, cfg]` | `kira_write_kinematics(family, cfg)` |
| `WriteJobs[cfg, "Masters"]` | `kira_write_jobs(cfg, KiraJobMode::Masters)` |
| `WriteJobs[cfg, "Reduce"]` | `kira_write_jobs(cfg, KiraJobMode::Reduce)` |
| `WritePreferred[masters, cfg]` | `kira_write_preferred(masters, cfg)` |
| `WriteTarget[targets, cfg]` | `kira_write_target(targets, cfg)` |

## Sub-process and parser

| .m symbol | C++ |
|---|---|
| `RunKira[cfg]` (typically `Run["kira jobs.yaml"]`) | `kira_run(cfg)` (fork / execve / wait) |
| `ReadMasters[cfg]` | `kira_read_masters(cfg)` |
| `ReadTargetTable[cfg]` (parses `kira_target.m`) | `kira_read_target_table(cfg)` |
| Mathematica `ToExpression["Plus[Times[...], j[...]]"]` | `kira_parse_expression(str, dlist_ctx)` |
| Big integers in `kira_target.m` | parsed via `fmpz_set_str` |

## High-level reduction wrappers

| .m symbol | C++ |
|---|---|
| `BlackBoxReduce[family, targets, cfg]` | `ibp::black_box_reduce` (Layer 15d) plus `amflow::run_reduction_step` for boundary-task execution |
| `BlackBoxDiffeq[family, masters, var, cfg]` | `BlackBoxDiffeq(family, masters, var, cfg) -> BlackBoxDiffeqResult` |
| `ReductionContext[family, masters, targets]` | `ReductionContext` struct |
| `BlackBoxOptions` (association of cfg flags) | `BlackBoxOptions` struct (passes through to KiraConfig) |
| `$KiraDumpIO` flag | env `AMFLOW_KIRA_DUMP=1` (and `BlackBoxOptions::dump_kira_io`) |
| `$KiraPath`, `$FermatPath` | `KiraConfig::kira_path`, `KiraConfig::fermat_path` |
