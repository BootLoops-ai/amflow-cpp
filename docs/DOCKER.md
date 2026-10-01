# Docker

> 中文版本: [DOCKER_zh.md](DOCKER_zh.md)

`ghcr.io/chang18/amflow-cpp` bundles **FLINT 3.4 + Kira 3.1 + FireFly +
amflow_cli** in a single Linux image so you don't need to compile the
dependency stack yourself.  This page documents how to use it.

> **What this image does *not* contain:**
> - **Fermat** (closed-source, free for academic use, but
>   non-redistributable per its license).  You must bind-mount your
>   local Fermat install at runtime — see ["Fermat: bind-mount
>   required"](#fermat-bind-mount-required) below.
> - **Mathematica** (closed-source, commercial).  AMFlow.cpp does not
>   need MMA at runtime; MMA is only used by maintainers to generate
>   reference oracle values.

---

## Quick start

```bash
# 1. Pull the image (or `docker build`, see "Building from source")
docker pull ghcr.io/chang18/amflow-cpp:latest

# 2. Run a sample bench: 1-loop cutbubble at eps = 1/100
docker run --rm \
    -v "$HOME/Ferl7:/usr/share/Ferl7:ro" \
    -v "$(pwd):/work" \
    ghcr.io/chang18/amflow-cpp:latest \
    /opt/amflow-cpp/tools/bench/cutbubble_1L_eps001_black_box_amflow_cpp.json \
    /work/out.json

# 3. Inspect the result
cat out.json | python3 -m json.tool | head -40
```

The first run pulls a multi-hundred-MB image; subsequent runs reuse the
cached layers.

---

## What's bundled

| Component | Version | Path inside image | License |
|---|---|---|---|
| `amflow_cli` | 1.1.0 | `/usr/local/bin/amflow_cli` | MIT |
| FLINT | 3.4.0 | `/opt/flint/{bin,include,lib}` | LGPL-3 |
| Kira | 3.1 | `/opt/kira/{bin,lib,share}` + `/usr/local/bin/kira` (symlink) | GPL-3 |
| FireFly | `kira-2` branch | `/opt/kira/lib/libfirefly.so.*` | (FireFly project license) |
| Oracle benches | from v1.1.0 | `/opt/amflow-cpp/tools/bench/` | (project license) |

The image is built on top of `ubuntu:22.04`; the runtime stage ships
without `-dev` packages but keeps every shared library the binaries
need (`libginac11`, `libcln6`, `libyaml-cpp0.7`, `libjemalloc2`,
`libgmp10`, `libmpfr6`, `zlib1g`, `libstdc++6`).

> Kira is built with `-Dkyotocabinet=false`.  Kyoto Cabinet is an
> optional pyRed key-value store backend that AMFlow.cpp doesn't use;
> Kira 3.1's `keyvaluedb.h` also has an include-order quirk
> (`#include <kcpolydb.h>` placed before `#include "pyred/config.h"`
> in the kyoto-cabinet branch) that breaks the compile when the macro
> guarding it isn't yet visible, so we disable it.

---

## Fermat: bind-mount required

Robert Lewis's Fermat is free for academic use but the license does
not permit redistribution.  AMFlow.cpp's IBP backend (Kira) needs it
for multivariate rational algebra, so the image **must** see your
local Fermat install at runtime.

The default Fermat lookup path inside the image is
`/usr/share/Ferl7/fer64` (this matches the hard-coded default in
[`src/api/run_json.cpp`](../src/api/run_json.cpp)).  Two mount styles
work:

**Style A — match the default path** (recommended, no env var needed):

```bash
docker run --rm \
    -v /your/local/Ferl7:/usr/share/Ferl7:ro \
    -v "$(pwd):/work" \
    ghcr.io/chang18/amflow-cpp:latest \
    input.json /work/out.json
```

**Style B — custom mount + `FERMATPATH` env var:**

```bash
docker run --rm \
    -v /custom/path:/opt/fermat:ro \
    -e FERMATPATH=/opt/fermat/fer64 \
    -v "$(pwd):/work" \
    ghcr.io/chang18/amflow-cpp:latest \
    input.json /work/out.json
```

If Fermat is missing, the entrypoint exits with code 2 and a clear
error message before Kira is even invoked.

You can grab Fermat from <https://home.bway.net/lewis/>; copy the
`Ferl7/` directory anywhere on your host and point the bind-mount at
it.

---

## Mount semantics

- `-v /your/Ferl7:/usr/share/Ferl7:ro` — Fermat install (read-only OK).
- `-v "$(pwd):/work"` — your working directory.  The container's
  `WORKDIR` is `/work`, so any **relative** path in a bench JSON is
  resolved here.  Outputs you write to `/work/...` appear on the host
  immediately.
- Input file path: pass **absolute paths** (either
  `/work/your-input.json` if you mounted it, or
  `/opt/amflow-cpp/tools/bench/...` for a bundled oracle).

`work_dir` inside the bench JSON is interpreted **inside the
container**.  The shipped oracle benches use `/tmp/...` which is
ephemeral.  If you want persistent intermediate Kira state, point
`work_dir` at a path under `/work`.

---

## Running the bundled oracle benches

The image ships every oracle bench from the v1.1.0 release.  List
them with:

```bash
docker run --rm ghcr.io/chang18/amflow-cpp:latest \
    ls /opt/amflow-cpp/tools/bench
```

A 4-loop banana example (~6 minutes, comparable to MMA):

```bash
docker run --rm \
    -v "$HOME/Ferl7:/usr/share/Ferl7:ro" \
    -v "$(pwd):/work" \
    ghcr.io/chang18/amflow-cpp:latest \
    /opt/amflow-cpp/tools/bench/banana_4loop_eps001_black_box_amflow_cpp.json \
    /work/banana_4loop_out.json
```

To verify against the committed MMA reference, run the comparison
script on the host (or inside the container with Python 3):

```bash
python3 tools/bench/compare_sampled.py \
    tools/bench/banana_4loop_eps001_mma_reference.json \
    banana_4loop_out.json
```

---

## Building from source

Build the image yourself when you want to pin a specific FLINT/Kira
version or modify `Dockerfile`:

```bash
git clone https://github.com/chang18/amflow-cpp
cd amflow-cpp
docker build -f docker/Dockerfile -t amflow-cpp:dev .
```

The build runs three big stages back-to-back:

1. **FLINT 3.4 from source** (~5–10 min).  Ubuntu 22.04 ships FLINT
   2.8.4 which is too old; we fetch the tarball from GitHub and
   compile with `--disable-static` to keep the image small.
2. **Kira 3.1 + FireFly (subproject)** (~10–15 min).  `meson setup`
   auto-pulls FireFly from `gitlab.com/firefly-library/firefly`
   (revision `kira-2`) via `subprojects/firefly.wrap`.  We set
   `-Dflint=true` so FireFly shares the FLINT 3.4 install instead of
   pulling the 2.8.4 fallback.
3. **amflow_cli** (~2 min).  Standard `cmake + ninja`.

Total cold-build wall-clock is ~20–30 min depending on hardware.

---

## Troubleshooting

### `ERROR: Fermat executable not found at "..."`

The entrypoint couldn't find Fermat.  Double-check your `-v` mount
(see ["Fermat: bind-mount required"](#fermat-bind-mount-required)).
Inside the container, the file `${FERMATPATH:-/usr/share/Ferl7/fer64}`
must exist and be `chmod +x`.

### `kira::ExceptionInternal: ...` or Kira crashes

Most often a **Fermat version mismatch** — the image was built against
some FireFly+Kira combo that expects a specific Fermat protocol.  If
your local Fermat is from a different vintage, try `Ferl7` (the
project's tested baseline) before reporting.

### Output JSON is empty / `permission denied` on output

Docker writes container output as `root` by default, which can hit
write-permission issues on the host bind mount.  Two fixes:

- **Run as your user:** add `--user "$(id -u):$(id -g)"` to the
  `docker run` flags.  The image is small enough that root vs user
  doesn't matter; this is purely about ownership of the result file.
- **Mount with explicit write permission and let root own the
  output**, then `sudo chown` afterwards.

### `unauthorized: authentication required` when pulling

The image is hosted on GitHub Container Registry (`ghcr.io`), which
is public.  If you hit auth errors, your docker client may have a
stale token — `docker logout ghcr.io` typically resolves it.

---

## Limitations

- **Linux/amd64 only** at the moment.  arm64 will require either
  cross-compilation (FLINT and Kira both build cleanly on ARM, but
  CI hasn't been wired up) or native-arm runners.
- **No Mathematica integration in-container.**  If you want to
  regenerate `*_mma_reference.json` against upstream MMA, run those
  steps on the host with your own Mathematica install — the C++ side
  doesn't need MMA and the Docker image deliberately doesn't either.
- **Image size is ~1 GB.**  Most of that is Kira + FireFly compiled
  artefacts; FLINT is ~120 MB.  Slimming will need either a deeper
  multi-stage prune of unused FireFly object files or a Distroless
  base.

---

## See also

- [`docs/USER_GUIDE.md`](USER_GUIDE.md) — JSON input format, oracle
  catalogue, full CLI options.
- [`docs/FAQ.md`](FAQ.md) — non-Docker build/install questions.
- [`tools/bench/README.md`](../tools/bench/README.md) — oracle bench
  layout and re-generation flow.
