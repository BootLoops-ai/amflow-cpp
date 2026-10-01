# Project-Wide Invariants

> 中文版本: [INVARIANTS_zh.md](INVARIANTS_zh.md)

These conventions are *not* per-layer — they are honoured everywhere in
the codebase.  Violating any of them will silently produce wrong results.

> _Note_: "Layer N" labels below refer to the legacy module organization
> used during the v1.0 port.  The current source tree is organized by
> **domain** (`numeric` → `algebra` → `ode`/`qft`/`ibp`/`pipeline`/`api`/`cli`)
> under `src/<domain>/`; see [`ARCHITECTURE.md`](ARCHITECTURE.md) for the
> current layout.

---

## 1. Precision is decimal at the API, binary internally

User-facing options use **decimal digits**:

```cpp
GlobalOptions g;
g.working_pre = 100;   // means 100 decimal digits
```

FLINT/Arb uses **binary bits** internally.  Conversion goes through
the `Precision` value type:

```cpp
amflow::numeric::Precision::from_decimal_digits(int);
```

with the formula `bits = ceil(digits · log2(10))`, clamped to ≥ 53.

**Rule:** *every* call into Arb that takes a `prec` argument must use a
value derived from `Precision::from_decimal_digits` (or one of the
project's other `Precision`-typed helpers).  Never hard-code.

---

## 2. FLINT RAII wrappers are move-only; rational/poly values are copyable

The wrappers `AcbValue`, `AcbVector`, `AcbMatrix`, `AcbPoly`, `FmpqPoly`,
`RationalFunction`, `AsyTerm` are *move-only*.

```cpp
AcbValue a; a.set_si(7);
AcbValue b = std::move(a);   // OK
AcbValue c = a;              // compile error (deleted copy ctor)
AcbValue d = a.clone();      // explicit deep copy
```

Reasons:

- An accidental copy of a million-bit `acb_t` would be invisible and
  expensive.
- We want `std::vector<AcbValue>` to grow only through moves.

**Caveat:** `RationalComplex` *is* copyable because it only holds two
`fmpq_t` (cheap).  `BlockEquation`, `BlockEquationNum`, `AsyExpansion`,
etc. inherit move-only-ness because they hold move-only fields.

When you need `vector<AcbValue> v = something()`, use `std::move(...)`
or call helper `clone_vec(...)` (defined per-layer where needed).

`std::vector<AcbValue>::assign(N, AcbValue())` does **not** compile
because `assign(n, value)` requires copy.  Use this idiom instead:

```cpp
v.clear();
v.resize(N);   // value-initialised; calls AcbValue() default ctor N times
```

The exact types intentionally differ here:
`amflow::numeric::RationalPolynomial`, `RationalFunction`, and
`RationalFunctionMatrix` are copyable value objects.  That boundary is
for exact algebraic ownership where ordinary C++ value semantics make
the higher-level domain model simpler.  Arbitrary-precision ball
wrappers may still stay move-only when copying would hide large
precision-dependent costs.

---

## 3. Reverse coefficient indexing in Layer 3 / 6 / 7

`construct_matrix(dx, ax, totalorder)` lays out an `(T+2)·N` column
sparse matrix.  Column index `k` corresponds to **the *(T+1-k)*-th
coefficient** of the unknown integral series, *not* the *k*-th.

That is:

```
                    column index k:  0           1           2         …    T+1
   represents the coefficient of:    f_{T+1}     f_T         f_{T-1}   …    f_0
```

This matches the Mathematica convention used in `DetermineBlockBoundaryOrder`'s
`fid[n]` formula:

```mathematica
fid[n_] := { ... [Mod[n-1, |new|]+1], ExtraXOrder + 1 - Quotient[n - 1, |new|] }
```

Layer 5 (`calcx1x2`) does **not** use this convention because it uses a
direct in-place recurrence and never touches sparse columns.  Only Layer
3-style routines (`construct_matrix`, `sparse_gaussian`, callers in Layer
6 and Layer 7) do.

If you change either side of the convention, change both.

---

## 4. Mathematica is 1-based, C++ is 0-based — `nh[All, n+1]` becomes `nh[i][n]`

Whenever a .m formula has

```mathematica
nh[[All, n+1]]   (* 1-based, so this is the (n+1)-th element = 0-based index n *)
```

the C++ translation is `nh[i][n]`, **not** `nh[i][n+1]`.

Layer-5 coverage catches this class of translation mistake quickly; do
not special-case the index arithmetic.

---

## 5. Chop tolerance is decimal, applied to *midpoints*

`AMFChop[x, 10^-ChopPre]` zeroes any number with magnitude below
`10^-ChopPre`.  Our equivalents:

```cpp
bool acb_is_chop_zero(acb_srcptr x, int chop_digits);
bool AcbValue::is_chop_zero(int chop_digits = chop_pre()) const;
```

Both look at the **midpoint** of the complex ball.  They do *not*
consult the radius.  This matches the .m semantics, which apply Chop
after `N[..., WorkingPre]` (which has no error ball).

If you need an interval-aware "definitely zero" test, build it from
`acb_contains_zero(...)` directly; do *not* mutate `is_chop_zero`.

---

## 6. Layer-7 Jordan decomposition is exact over Q

`NormalizeMat` uses exact rational residue matrices:

```text
RationalMatrix -> fmpq_mat residue -> jordan_decomposition_exact
```

If the characteristic polynomial has an irreducible factor of degree
greater than 1 over Q, the exact-Q Jordan path must fail clearly.  The
current algebraic fallback may still use exact FLINT `qqbar` data for
leading exponents and integer-floor decisions, but public
`NormalizeMat()` must not silently pretend that a full algebraic Jordan
basis exists when its downstream representation is rational-only.
`CalcZero` is allowed to use the separate internal
`normalize_mat_for_calc_zero()` path, because that path also carries the
corresponding block-local algebraic rotations all the way through
`calcx00()`.

---

## 7. `acb_mat_t` is an array typedef, not a class

```c
typedef acb_mat_struct acb_mat_t[1];
```

Implications:

- `new acb_mat_t` does **not** compile.  Use `new acb_mat_struct`.
- `std::vector<acb_mat_t>` does not compile (arrays not assignable).
  Use `std::vector<acb_mat_struct*>` (heap pointers) or
  `std::vector<acb_mat_struct>` plus manual initialisation.
- Functions take `acb_mat_t mat` (array decays to `acb_mat_struct*`).
  Pass either the variable name or the struct pointer; do *not* take
  `&mat` (that gives `acb_mat_struct(*)[1]`, the wrong type).

The same applies to `acb_t`, `arb_t`, `fmpq_t`, `fmpz_t`, `acb_poly_t`,
`fmpq_poly_t`.

---

## 8. `fmpq_cmp_si(x, c)` takes one integer, not a fraction

This catches everyone once.

```c
int fmpq_cmp_si(const fmpq_t x, slong c);   // compares x against c (integer)
```

There is **no** `fmpq_cmp_si(x, p, q)` for "is x equal to p/q?".  Use
the helper

```cpp
bool fmpq_eq_pq(const fmpq* x, long p, long q);
```

(declared locally in tests) or build a temporary `fmpq_t` and call
`fmpq_equal`.

---

## 9. `fmpq_mat_t` follows the same array-typedef rule

```c
typedef fmpq_mat_struct fmpq_mat_t[1];
```

Output-parameter APIs such as `jordan_decomposition_exact` therefore
take pre-initialized `fmpq_mat_t` objects.  Do not try to return or move
`fmpq_mat_t` directly; wrap it in RAII locally or pass output matrices
from the caller.

---

## 10. Exact Jordan block counts use kernel-tower differences

In `src/ode/jordan.cpp`, for one eigenvalue:

```text
V_s = ker((A - lambda I)^s)
r[s] = dim(V_s) - dim(V_{s-1})
blocks_of_size[s] = r[s] - r[s+1]
```

At the maximum depth, `r[depth+1]` is intentionally `0`, because there
are no blocks larger than the largest discovered chain.  Confusing this
with `dim(V_{depth+1}) = dim(V_depth)` breaks the exact block count.

Layer 7 diagonal ToFuchsian follows the same exact-Q rule: `ReduceL0`,
`FindProjector`, `Balance`, and `InvBalance` operate on rational
matrices.  Do not route the projector loop through Arb eigenvectors or
numeric rank tests; irreducible leading Poincare rank is an exact
algebraic error.

---

## 11. `AcbRationalFunction` does NOT reduce

Unlike Layer 1's `RationalFunction` (which keeps `gcd(num, den) = 1` and
`den` monic), Layer 7's `AcbRationalFunction` (acb-coefficient rational
function) **does not** reduce.  Acb has no exact GCD so reduction is
ill-defined.

Operations build up `num` and `den` blindly.  Common-denominator
clean-up happens once, when evaluating at η = 0 via
`value_at_zero()`, which strips a leading η prefix and then divides
constants.

If you need to inspect a result, materialise to `acb_mat_t` via
`acb_rational_matrix_value_at_zero(...)`.

---

## 12. `MpolyContext` is shared by `shared_ptr`; cross-context ops throw

Every `Mpoly` / `Mfrac` / `MpolyMatrix` / `MfracMatrix` carries a
`shared_ptr<MpolyContext>` identifying the variable list it lives in.

```cpp
auto ctxA = MpolyContext::make({"eta", "eps"});
auto ctxB = MpolyContext::make({"eta", "eps", "D[1]", "D[2]"});

Mpoly p = Mpoly::parse(ctxA, "eta + eps");
Mpoly q = Mpoly::parse(ctxB, "D[1]");

Mpoly r = p + q;   // throws std::runtime_error: ctx mismatch
```

To bridge two contexts, **round-trip through string**:

```cpp
Mpoly p_in_B = Mpoly::parse(ctxB, p.to_string());   // ok
```

This is the project's discipline:

- Layer 12 carries the family's "base" context.
- Layer 13–14 typically extend it with `D[1] ... D[k]`, eta, region
  scales — they create a new `MpolyContext` and re-parse whatever
  came from the parent.
- Layer 16 child AMFSystems (especially SingleMass children) own
  their own context; the parent must re-parse propagators before
  passing them in.

There is no warning if you mix contexts that *happen* to have
identical variable names — the runtime check uses `shared_ptr`
identity.  This is intentional: two contexts with the same names
might still have different variable orderings, and we must not
silently re-interpret an exponent vector.

---

## 13. The reverse-column convention of §3 is preserved by Layer 14's boundary tables

Layer 14e (`boundary_integrals`) eventually feeds boundary terms to
Layer 6 (`calc_inf` → `calc_taylor`).  The terms are stored as
`(sub_index, coefficient)` tuples where `sub_index` indexes a
*sub-master* (a master integral on a daughter family).

When this list is consumed at solve time by Layer 16's `solve_one_eps`,
the resulting `BoundaryEntry { mu, value }` for the parent family must
respect the same reverse indexing convention as §3:

- `mu` is the η-exponent (a complex rational), passed directly to
  Layer 6's `read_bcs`.
- `value` is the linear combination of child solutions evaluated at
  the matching sub_index, then summed into the appropriate
  `BoundaryEntry`.

If you ever change the column-direction convention in `construct_matrix`
(§3), you must also flip the order in which Layer 16 fills
`BoundarySpec` from `RegionBoundary::terms`.  The two are coupled
through `read_bcs` → `calc_taylor` → `sparse_gaussian`.

---

## 14. Kira sub-process contract

Layer 15c (`kira_run`) shells out to the Kira binary.  The contract:

1. **Absolute paths only.**  `KiraConfig::kira_path` and
   `KiraConfig::fermat_path` must be absolute.  We validate this at
   `kira_run` entry; relative paths fail silently because `fork()`
   inherits the caller's cwd then `execve` resets it.
2. **`FERMATPATH` env var.**  We always set `FERMATPATH = cfg.fermat_path`
   in the child env.  Kira's startup probes for Fermat using this
   variable; if it's missing, Kira aborts before producing any output.
3. **Workdir is `cfg.work_dir`.**  We `chdir(cfg.work_dir)` in the
   child before `execve`; all yaml files are pre-written there.
4. **Two-pass for `BlackBoxReduce`.**  First pass uses
   `KiraJobMode::Masters` (writes only the master list); second pass
   uses `KiraJobMode::Reduce` (writes the target table).  The two passes
   share a workdir.
5. **Workdir reuse is request-shape guarded.**  The
   `ibp::black_box_reduce` / `black_box_diffeq` path writes
   `.amflow_kira_fingerprint` after a successful Kira run.
   The fingerprint covers family shape, kinematics, job mode, numeric
   substitutions, IBP options, dimension base, preferred masters, and
   targets.  If the next request for the same workdir matches and the
   expected Kira output files are present, the cached path parses those
   files without launching Kira again.  If the request differs or outputs
   are incomplete, only Kira-generated artifacts (`config`, `tmp`,
   `results`, `sectormappings`, `firefly_saves`, `preferred`, `target`,
   and `jobs.yaml`) plus the old fingerprint are cleared before new
   inputs are written.  This avoids Kira auxiliary-file reuse across
   incompatible numeric-variable shapes while preserving caller-owned
   parent directories.
6. **Output parsing.**  `kira_read_masters(cfg)` parses
   `<workdir>/results/<family>/masters`.  `kira_read_target_table(cfg)`
   parses `<workdir>/results/<family>/kira_target.m` via
   `kira_parse_expression`, which supports MMA prefix syntax and
   big-integer parsing through `fmpz_set_str`.

Pre-flight failures (missing binaries, missing yaml files, bad family
name) raise `std::runtime_error` with stderr captured.  Post-Kira
failures (Kira ran but returned non-zero) raise `std::runtime_error`
with the exit code and the last 4 KiB of stderr.

If you need to debug Kira, set env `AMFLOW_KIRA_DUMP=1` (or
`BlackBoxOptions::dump_kira_io = true`) — the workdir will not be
cleaned up after the run.

---

## 15. `global_values` takes precedence over `master_values` in
   SingleMass children

In Layer 16's `AMFSystem::solve_one_eps`, when reading boundary
contributions from a child node, the order of preference is:

```cpp
const auto& src = child.solutions[eps_index];
const auto& boundary_values =
    src.global_values.empty() ? src.master_values
                              : src.global_values;
```

Why: SingleMass children's `master_values` are the values of the
*sub-master integrals* (one per connected component after factorisation).
The Γ-function prefactor that arises from the loop redefinition is
applied in `solve_one_eps` and stored separately in `global_values`.
The parent expects the *prefactor-included* values when stitching
sub-systems together.

Reading `master_values` unconditionally drops the Γ prefactor and
produces an across-the-board error of magnitude `Γ(1+ε)^L`.
**Do not** simplify the conditional unless you also propagate the
prefactor at the call site.

---

## 16. `BuildTaylor` must apply the diagonal correction *before* `NHEquations`

Mathematica's `BuildTaylor[mat, ini]` returns the **fully-corrected**
matrix `mat'[i,j] = η^{-ini[i]} · mat[i,j] · η^{ini[j]}` if `i ≠ j`,
and `mat'[i,i] = mat[i,i] - ini[i]/η` if `i = j` — both via `Together`
to keep entries as reduced rational functions.  This corrected matrix
is then passed to `CalcTaylor`, which internally invokes `NHEquations`
on it.

The correct C++ implementation lives in `src/ode/inf.cpp`:

- `build_taylor_symbolic(mat, int_offsets, ini_rat)` does the full
  BuildTaylor (off-diagonal `η^(int_offsets[j] - int_offsets[i])`
  shifts *and* diagonal `mat[i,i] - ini_rat[i] / η` subtraction) in
  `RationalFunction` (i.e. `fmpq`) form, applying `Together` via the
  existing operator overloads.
- `prepare_taylor_system(m_pure, …)` then runs `nh_equations(…, Taylor)`
  on the already-corrected matrix.  `PoincareRank`, `dx`, and `ax` all
  come out of a single consistent computation.
- No code path may subtract `ini[i]·dx_inner[k]` from
  `axexp[i][i][k+rank]` numerically *after* `NHEquations`.

**Concrete invariant going forward:** any future code path that
computes `nh_equations(..., Taylor)` on a matrix with a non-zero
`ini` offset MUST pre-compose the `mat[i,i] - ini[i]/η` subtraction
symbolically.  There must be no numeric post-correction step.  If a
new use case needs complex-valued `ini`, extend the symbolic
representation rather than re-introducing the split.

---

## 17. AMFSystem tree's `solve()` is bottom-up; never mutate visited
    children

Each non-ending `AMFSystem` node holds:

- `region_boundaries` — a list of `RegionBoundary { region, mu,
  terms }`, where `terms` references children by `sub_index`.
- `differential_eq_matrix` — the per-eps ODE matrix (built once by
  `setup()` via `BlackBoxDiffeq`).
- `children` — `vector<AMFSystem>` (one per region per family slot).

`solve(epslist)` does:

```text
for each eps in epslist:
   1. for each child in children:
        recurse: child.solve_one_eps(eps)
   2. build BoundarySpec from region_boundaries + child solutions
   3. evaluate differential_eq_matrix at this eps
   4. call amflow(de_at_eps, BC) → master_values
   5. (if SingleMass parent) apply Γ prefactor → global_values
```

The invariant: **once a child's `solve_one_eps(eps)` returns, the
child's solution for that eps is fixed**.  Do not mutate a child's
solutions array from the parent.  Multiple parents may share a child
in principle (DAG, not strict tree), so mutation would corrupt other
parents.  In practice today the structure is a strict tree, but we
preserve the discipline so the structure can later become a DAG
without breaking anything.

This matters when adding caching, parallel-eps solving, or any kind
of memoisation across the tree: the natural place is on the child
node, written exactly once per (child, eps) pair.

---

## 18. `EndingScheme` order matters in `AMFSystemOptions::ending_schemes`

`AMFSystemOptions::ending_schemes` is a `vector<EndingScheme>` listing
the strategies to *try* for each node, in order.  The first one that
returns `ending_q == true` wins.

Default order: `{Tradition, SingleMass}`.  This means a node that *can*
be terminated by a vacuum lookup will be terminated that way; if not,
SingleMass is tried; if neither works, the node recurses (meaning Kira
+ derivative + sub-AMFSystems).

If you reorder — e.g. `{SingleMass, Tradition}` — the same family will
produce different (but mathematically equivalent) AMF trees, and the
runtime cost / numerical accuracy may shift dramatically.  When
debugging numerical errors, **first** try forcing a single ending
scheme (`{Tradition}` only or `{SingleMass}` only) to localise the
problem.

---

## 19. `factor_poly(Taylor, rank)` uses `η^(rank+1)` literally — no clamp

In a naive port the per-block `factor` used by `nh_equations` would be

```
factor = η^(rank + 1)    for Taylor mode
```

where `rank = poincare_rank(mat)` may be `-1` (matrix regular at η=0).

**Correct:** for `rank = -1`, `factor = η⁰ = 1`.  Mathematica does this
literally and our `factor_poly` now does the same.

**Anti-pattern:** clamping with `if (power < 1) power = 1;` forces
`factor = η` even when the matrix is regular, inflating `dx` by a
spurious `η` factor and re-triggering the column-shift problem that
INVARIANTS §16 is meant to prevent.

Preserve `factor = η^max(0, rank+1)` semantics (or equivalent — the
polynomial must be `1` when `rank+1 ≤ 0`) and never re-introduce the
clamp.

---

## 20. `determine_boundary_order` accepts rational `power` — never round

`include/amflow/inf.hpp` exposes two overloads:

```cpp
std::vector<long>
determine_boundary_order(const RationalMatrix& mat,
                         const std::vector<long>& power);

std::vector<long>
determine_boundary_order(const RationalMatrix& mat,
                         const std::vector<RationalFunction>& power_q);
```

**Use the rational overload from `build_boundary`.**  Mathematica's
`DetermineBoundaryOrder` accepts any `ini` (including non-integer
rationals like `-99/100` that arise when evaluating `1 - ε` at
`ε = 1/100`).  Rounding to the nearest `long` gives a different
answer:

```
MMA  DetermineBoundaryOrder[deinf, {-99/100}]  = {0}
MMA  DetermineBoundaryOrder[deinf, {-1}]       = {-1}
```

If your code path computes a rational `npat`, **do not** round.  Pass
`RationalFunction::from_fmpq(q)` into the rational overload.

The C++ entry point is `ode::determine_boundary_order(...)`, which
accepts `std::vector<numeric::RationalFunction>` exponents.  Keep the
rational exponent path when wiring AMFlow boundary orders into the
execution layer.

---

## 21. `deinf = -1/η² · de(1/η)` — both halves are required

`RationalFunction::substitute_inverse_eta()` implements `r(1/η)`
**only**.  It does *not* multiply by `-1/η²`.

In Layer 6 the helper `make_deinf(de)` wraps both steps:

```cpp
auto subbed = de.substituted_inverse_eta();
RationalFunction prefactor = RationalFunction::monomial(-2, neg_one);
out(i, j) = subbed(i, j) * prefactor;
```

In Layer 16 (`AMFSystem::build_boundary`) the same `deinf` is built
inline from each per-ε-evaluated `Mfrac`; that code MUST also apply
the `-1/η²` factor explicitly.  Forgetting it makes `deinf` a
polynomial in η where MMA has a simple pole at η=0 —
`determine_boundary_order` then sees an entirely different matrix and
returns nonsense borders, and `rb.integrand_terms` ends up empty for
non-zero-sector regions.

---

## 22. Pattern-zero BC補 in `solve_one_eps`

Mathematica's `AMFSystemSolution` ends with (AMFlow.m line 1149):

```mathematica
bc = MapThread[Join[#1,#2]&,
               {bc, Transpose[Thread /@ Thread[(pattern/.epsrule) -> 0]]}];
```

Every pattern group's per-master μ value must appear in `bc_sorted[i]`
with value `0` — even if no boundary-integrand contribution matches
that μ.  Without the pattern-zero entries, a sub-system whose regions
all have `border = -1` ends up with an empty `bc_sorted`, and
`amflow(de, {}, prec)` returns identically zero.

Our implementation stores pattern values in `AMFSystem::bc_pattern_`
(filled in `build_boundary` from the result of `boundary_pattern(all_powers)`)
and appends `(μ_i → 0)` entries in `solve_one_eps` **before** calling
`amflow(de, bc_sorted, prec)`.  Do not skip this step.

---

## 23. `BoundaryFamily::terms` layout is `[input][order][laporta_term]`

Output of `boundary_integrals` — documented at
`include/amflow/boundary.hpp` — has three nested dimensions:

```cpp
std::vector<std::vector<std::vector<LaportaTerm>>> terms;
//  └─ input_i ─┴─ order_j ─┴─ laporta_term_k ─┘
```

where:

- `input_i` ranges over input integrals (= `diff_masters` entries
  when called from `build_boundary`);
- `order_j` ranges over `0 .. border[input_i]` (so when `border[i] = -1`
  the inner vector has length 0);
- each innermost entry is a `LaportaTerm` = `{ indices, coef }`.

Mathematica's `BoundaryIntegrals` effectively does
`DynamicPartition[flat_apart_pieces, partition]` to produce exactly
this shape; the inner `GatherBy` + `Plus` aggregation merges all PFD
pieces of the same `(input_i, order_j)` that share a family
signature into a single `Mfrac` before passing through
`LaportaIntegrals`.

**Anti-pattern:** indexing `terms` by the flat
`pfd_per_input[flat_idx]` layout (where `flat_idx` iterates over
`(input, order)` as a single sequence and per-PFD-piece as the inner
index) means `terms.size()` is not equal to `diff_masters.size()`,
which causes the `build_boundary` loop to either skip entries or write
to the wrong slot.  The flat layout does not match MMA's output shape
and has no downstream consumer.

---

## 24. Layer 14 is contracted to standard-QFT-only semantics

The C++ port covers standard QFT (quadratic propagators) only — linear
propagators / gauge-link / Wilson-line workflows are out of scope (see
`docs/ROADMAP.md`).  Two Layer-14 entry points enforce this contract
and **must throw** on inputs outside the standard-QFT shape:

- `branch_momenta(family)` requires `coe1 == 1` for every entry; the
  `1/2 * coe2` Wilson-line halving branch is not implemented.
- `sp_list_to_dlist_symbol(family, sp_list)` rejects sp-list entries
  whose bilinear coefficient is not a constant integer.

Both throw on linear-propagator / gauge-link inputs.  Do not relax
these checks to "support" mixed quadratic/linear families — the rest
of the pipeline (Layer 13 candidate selection, Layer 14 region
decomposition, Layer 16 ending detection) has no gauge-link branch.
