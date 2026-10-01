# Contributing

> 中文版本: [CONTRIBUTING_zh.md](CONTRIBUTING_zh.md)

AMFlow.cpp is a parity-driven C++17 reimplementation of the AMFlow
algorithm.  Before making behavioural changes, read the
project-specific workflow in
[`docs/CONTRIBUTING.md`](docs/CONTRIBUTING.md).

The short version:

- Upstream Mathematica AMFlow
  (<https://gitlab.com/multiloop-pku/amflow>) is the behavioural
  specification.
- Keep Mathematica reference outputs and raw caches when they establish
  or protect parity.
- Add focused tests for changed domains (`tests/test_<domain>_*.cpp`).
- Update `README.md` and `docs/AUDIT_MMA_PARITY.md` when project status changes.
- Do not relax numerical tolerances or remove parity checks to make a
  C++ mismatch pass.

For local verification, use:

```bash
cmake -S . -B build -DAMFLOW_BUILD_DRIVER=ON
cmake --build build -j
ctest --test-dir build --output-on-failure -j 4
```

Tests that require Kira/Fermat are skipped when those tools are not
available. Mathematica / Wolfram Engine is only required when
regenerating reference data.

## About this implementation

The C++17 source in this repository was developed primarily by AI
coding agents under the direction of the maintainer.  Every
behavioural change is gated by numerical-parity verification against
the upstream Mathematica reference; the 228 oracle benchmarks under
`tools/bench/` and the 548-case GoogleTest suite are the contract
that any change must honour.

## Contact

For bug reports and feature requests, please open a GitHub issue at
<https://github.com/chang18/amflow-cpp/issues>.

For collaboration inquiries, citation questions, or anything that
doesn't fit a public issue, contact the maintainer at
<3250800970@qq.com>.
