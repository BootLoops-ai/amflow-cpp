> **Layout note for this repository.** In the BootLoops fork of AMFlow.cpp this
> harness lives at `bootloops-wrappers/` inside the engine repository, so the
> engine checkout is the parent directory, not `../amflow-cpp`. Build the engine
> at the repository root (`cmake -S . -B build -DCMAKE_BUILD_TYPE=Release &&
> make -C build -j$(nproc)`), then run the scripts from this directory with
> `export AMFLOW_REPO="$(git rev-parse --show-toplevel)"` (also the default).
> The body below is copied unchanged from the main BootLoops repository except
> for the paths and this License section.

# AMFlow-port wrapper suite

Validation harness for the `amflow-cpp` engine shipped alongside (original
code, MIT; see the License section below). Build the engine first:

```sh
cd ..
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
make -C build -j$(nproc)
```

Then, from this directory:

| script | what it does |
|---|---|
| `run_tier1.sh [jobs]` | re-runs upstream's 228 oracle bench cases in `../tools/bench/` and compares vs upstream's committed AMFlow-package reference values; on `006a275` in our environment (2026-06-11) 224/228 pass — three `pentabox_2L` cases abort in the eta→0 residue step and `wzbox_2L` under-resolves at the committed truncation orders; upstream reports 228/228 (README status table; `docs/AUDIT_MMA_PARITY.md`); the four differences are specific to our environment and under investigation |
| `gen_tier23_configs.py` | generates the tier-2/3 solve_integrals configs (bubble/box/massive-bubble/sunrise at 30/100/300 goal digits) |
| `run_tier23.sh` | runs the tier-2/3 configs |
| `ref_closed_forms.py --out refs.json` | engine-independent anchors: Gamma / 2F1-polylog closed forms (mpmath) |
| `ref_sunrise_bessel.py` | equal-mass sunrise anchor via Bessel moments (elliptic case) |
| `compare_orders.py a.json b.json` | order-by-order Laurent-coefficient comparison |
| `run_tier5.sh` | cost-vs-digits and cost-vs-eps-order scaling fits |

Environment: `AMFLOW_REPO` points at the engine checkout (default: the parent
directory of this one); the tier scripts write results under
`$AMFLOW_PORT_ROOT/results/` (default: this directory). Tier-1/2/3 need a
`kira` binary on PATH for the reduction legs; the closed-form anchors need
only python3 + mpmath.

## Credits

The tier-1 leg re-runs the 228-case oracle benchmark built by the AMFlow.cpp
contributors against X. Liu and Y.-Q. Ma's Mathematica AMFlow; upstream
reports 228/228 on its reference setup (README status table;
`docs/AUDIT_MMA_PARITY.md`), and the four cases that differ in our environment
(table above) are recorded as environment observations, not as claims against
upstream. The closed-form anchors follow G. Duplančić and B. Nižić, Eur. Phys.
J. C 20 (2001) 357 [hep-ph/0006249] (massless box) and, for the equal-mass
sunrise, S. Groote, J. G. Körner and A. A. Pivovarov, Annals Phys. 322 (2007)
2374 [hep-ph/0506286] and D. H. Bailey, J. M. Borwein, D. Broadhurst and M. L.
Glasser, J. Phys. A 41 (2008) 205203 [arXiv:0801.0891]. External cross-checks
used pySecDec — S. Borowka, G. Heinrich, S. Jahn, S. P. Jones, M. Kerner, J.
Schlenk and T. Zirke, Comput. Phys. Commun. 222 (2018) 313 [arXiv:1703.09692]
(and later versions by the pySecDec collaboration). Numerics use mpmath
(Fredrik Johansson and the mpmath contributors); the runners use GNU parallel
by Ole Tange, whose author asks that it be cited (O. Tange, *GNU Parallel
2018*, doi:10.5281/zenodo.1146014; or O. Tange, "GNU Parallel — the
command-line power tool", ;login: The USENIX Magazine, February 2011, 42-47).
We thank all of them.

## License

The files in bootloops-wrappers/ are Copyright (c) 2026 Anthropic, PBC and released under the MIT License, the same terms as
the LICENSE at this repository's root (which is upstream AMFlow.cpp's; the root NOTICE is upstream's text followed by sections this fork appended). Created by Matthew D. Schwartz;
code written by Claude (Anthropic) under his supervision. This is not an officially supported Anthropic product; it is
maintained by Matthew D. Schwartz (https://www.bootloops.ai). Third-party programs these scripts invoke — Kira (Maierhöfer,
Usovitsch, Uwer; Klappert, Lange, Wu), GNU parallel (O. Tange), mpmath (F. Johansson et al.) — keep their own licenses.
