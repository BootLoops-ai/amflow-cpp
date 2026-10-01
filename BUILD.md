# BUILD.md — BootLoops tooling build notes for amflow-cpp

Operational build/run notes for running this AMFlow.cpp port as the
IBP/numeric oracle for the BootLoops elliptic-bootstrap program. These are *deployment* notes; they do not change the algorithm
(see `README.md` / `docs/ARCHITECTURE.md` for the science).

## (a) jemalloc requirement for Kira

With the default glibc allocator, the heap can fragment during Kira's
"Select equations recursively" phase on our largest systems, enough that a
single contiguous allocation throws `std::bad_alloc` with memory still free
in aggregate. jemalloc — a documented Kira build option (`-Djemalloc=true`;
the Kira README reports it "often by more than 20%" faster in their
experience when FireFly is used) — fragments far less and (with the decay
knobs below) returns dirty pages to the OS, lowering peak RSS.

Build Kira with jemalloc linked in:

```sh
apt-get install libjemalloc-dev
# in the Kira meson build:
meson setup build -Djemalloc=true
```

Confirm with `kira --version`; the build configuration line then reports
`with jemalloc`.

In addition, `src/ibp/kira_run.cpp` (the launcher that forks the Kira
child) **auto-`LD_PRELOAD`s the system libjemalloc** into the Kira child
as a belt-and-suspenders so that older binaries / Fermat / FireFly
children inherit it too, and sets a tuned `MALLOC_CONF`
(`background_thread:true,dirty_decay_ms:5000,muzzy_decay_ms:5000,narenas:32`).

Controls (all read from the environment by the launcher):

- `AMFLOW_KIRA_NO_JEMALLOC=1` — disable the auto-`LD_PRELOAD` entirely.
- `AMFLOW_JEMALLOC_PATH=/path/to/libjemalloc.so.2` — override the default
  `/lib/x86_64-linux-gnu/libjemalloc.so.2`.
- A pre-existing `LD_PRELOAD` or `MALLOC_CONF` in the environment is
  respected and not overwritten.

The preload only takes effect if the resolved `libjemalloc.so.2` is
readable; otherwise the launcher silently proceeds without it.

## (b) Kira-side SQLite patch

The Kira fork carrying the SQLite transaction fix is a separate repository
(the BootLoops Kira fork, GPL-3.0-or-later; Patch 1 in its PATCHES.md; linked
from upgrades/ENGINES.md in the main BootLoops repository). Nothing from it is
included here; this port only launches the built `kira` binary as a
subprocess.

That patch guards `DataBase::begin_transaction()` /
`DataBase::commit_transaction()` (in `src/kira/dataBase.cpp`) with
`sqlite3_get_autocommit()`: a nested `BEGIN` on an already-open
connection is skipped, and a `COMMIT` with no active transaction becomes
a no-op (rather than printing "cannot commit - no transaction is active"
on essentially every commit during long reductions). The same
`PATCHES.md` documents the jemalloc build described in (a). Use that
fork's binary as the `kira` invoked by this port.

## (c) Per-family run knobs

### `kira3_truncate_sp`

Kira 3 seeding lever, surfaced as an explicit per-family JSON knob under
`amf_options.blackbox.kira3_truncate_sp`. It **defaults
OFF** because `truncate_sp` regresses on `box1`. It is a per-family
decision:

| family            | `kira3_truncate_sp` |
|-------------------|---------------------|
| 3-loop crossed-box / two-mass 5-point top | `true`  (required) |
| box1              | `false` (regresses) |

(The sibling lever `kira3_sectorwise` also defaults OFF — it hangs on the
aux-mass DE reduce.)

### `eps_order = K + 2L` convention

The AMFlow convention is `eps_order = K + 2L`, i.e. the highest fitted
power is `eps^K` with `K = eps_order - 2L` for an `L`-loop integral.
Consequences:

- The **finite part** `eps^0` requires `eps_order >= 2L`.
- "Through `eps^K`" requires `eps_order = K + 2L`.

So e.g. `eps_order = 4` on a 2-loop family fits only up to `eps^0`.
`solve_integrals` prints the fitted power range; if `eps_order < 2L` (the
finite part would not be fitted) it warns and auto-corrects `eps_order` to
`2L` (set `AMFLOW_EPS_ORDER_NO_BUMP=1` to keep the requested under-order);
and if every fitted order of an integral vanishes, so that truncation would
return a manufactured `0 @ order 0`, it stops with an error naming the
`eps_order` needed (set `AMFLOW_ALLOW_TRUNCATED_ZERO=1` to demote that to a
warning when only the pole part is wanted).

## (d) AMFLOW_KIRA_FFSAVE — resumable embedded-Kira FireFly reduces

The embedded Kira reduce runs FireFly with Kira's own save mechanism
always on (`reductionFF.cpp` hardcodes `bool save = true`; `ff_save/`
accumulates per-function saved states during the run).  Upstream
(re)launches clear the reduce work dirs (a clean-start default), so that
state does not survive a relaunch and an interrupted multi-hour FireFly
reconstruction restarts from probe zero.

`AMFLOW_KIRA_FFSAVE=1` makes the embedded lifecycle RESUMABLE
(default unset/0 = the upstream clean start, the identical code path to
upstream; the ON path has no automated test in this tree):

- `black_box_amflow_single` keeps the `reduce/` + `amf/` trees across a
  relaunch; `ibp::reduce`/`ibp::diffeq` byte-compare the Kira input
  files they are about to write (plus a `.amflow_ffsave.job` sidecar:
  numeric_values incl. the argv-only `-sd` eps, kira executable)
  against the existing dir.  Match ⇒ PRESERVE: Kira re-runs in place
  and FireFly `resume_from_saved_state()` fast-forwards past done
  functions (`Loading saved states` / `Resuming in prime field`).
  Mismatch ⇒ FRESH; if the old dir held live `ff_save` states it is
  moved to `<dir>.ffsave_stale` (ONE slot) for salvage/diagnosis
  instead of being deleted.
- Resume-engagement guard: a save interrupted in FireFly's pre-shift
  window is unresumable (FireFly hard-fails `Shift file not found!`);
  a non-signal Kira failure inside a PRESERVED dir clears the state
  and retries ONCE fresh, so a bad save can never wedge the job.
- `kira_stdout.log` is opened O_APPEND (not O_TRUNC) so an interrupted
  attempt's streamed FireFly `Done: N / M` counters survive the relaunch.

What survives where after an external kill (knob ON): `ff_save/`
(states/probes/validation.gz/tags/shift/anchor_points — resumable iff
`shift` exists), `tmp/<fam>/SYSTEM_*.gz` (initiate output, reused on
rerun), `results/kira.db` (checkpoint), `sectormappings/`,
`kira_stdout.log` + `firefly.log` (append-forensic).  NOT covered:
`build_boundary` sub-family reduces (pid+counter-unique temp dir names
never collide across launches, so there is nothing to resume), and the
Fermat→FireFly SIG6 auto-staging rewrite of jobs.yaml in a preserved
dir (next relaunch sees mismatched inputs and goes fresh).
