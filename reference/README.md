# reference/ — pointer to upstream Mathematica AMFlow

This C++ project is a reimplementation of the algorithms in the
Mathematica package **AMFlow** by Xiao Liu and Yan-Qing Ma:

- Upstream repository: <https://gitlab.com/multiloop-pku/amflow>
- Algorithm paper: Liu & Ma, *AMFlow: A Mathematica package for Feynman
  integrals computation via auxiliary mass flow*, Comput. Phys. Commun.
  **283** (2023) 108565,
  [doi:10.1016/j.cpc.2022.108565](https://doi.org/10.1016/j.cpc.2022.108565).

The upstream Mathematica source is **not** vendored in this repository.
Documentation and code comments in this repo refer to symbols and line
numbers in the upstream files (`AMFlow.m`,
`diffeq_solver/DESolver.m`, `ibp_interface/Kira/interface.m`); to follow
those references and to regenerate Mathematica reference data with
`tools/bench/*.wl` or `tools/math_ref/*.wl`, you need a local copy of
upstream AMFlow at this path:

```
reference/amflow-master/
```

## How to populate `reference/amflow-master/`

Pick one of the following — both are fine for development.

### Option 1 — clone upstream into this repo (gitignored)

```bash
git clone https://gitlab.com/multiloop-pku/amflow.git reference/amflow-master
```

`reference/amflow-master/` is in `.gitignore`, so the cloned tree will
not be tracked by this repository.

### Option 2 — symlink to an existing clone elsewhere

```bash
git clone https://gitlab.com/multiloop-pku/amflow.git ~/src/amflow
ln -s ~/src/amflow reference/amflow-master
```

The symlink itself is gitignored.

## Notes

- The line numbers cited in this project's documentation track the
  upstream `master` branch at the time of the v1.0.0 release.  If
  upstream has moved on, search by symbol name rather than by line
  number.
- All algorithmic credit belongs to the AMFlow authors above; this
  repository contributes only the C++17 reimplementation and a
  numerical-parity test/benchmark harness.
