## Summary

<!-- One short paragraph: what changes, and why. -->

## Parity / correctness

- [ ] No behavioural change (refactor / docs / build only)
- [ ] Behavioural change — covered by an existing test
- [ ] Behavioural change — adds a new test
- [ ] Updates Mathematica reference data — explain why below

<!-- If this fixes a parity bug against upstream AMFlow
(https://gitlab.com/multiloop-pku/amflow), link the offending family /
workflow and the rel-error before/after. -->

## Build & test

- [ ] `cmake --build build` clean (no new warnings)
- [ ] `ctest --test-dir build` green (note any deliberately skipped
      tests, e.g. Kira-dependent)
- [ ] If touching benchmarks: ran the relevant `tools/bench/*` driver
      and updated reference JSONs intentionally

## Notes for reviewers

<!-- Tricky parts, design decisions you'd like a second opinion on. -->
