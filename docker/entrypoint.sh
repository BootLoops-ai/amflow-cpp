#!/usr/bin/env bash
#
# AMFlow.cpp Docker entrypoint.
#
# 1.  Quick-help / version bypass — no Fermat needed for printing usage.
# 2.  Fermat detection — non-redistributable, must be bind-mounted by
#     the user.  We look at $FERMATPATH first (if user set it), then at
#     the hard-coded default that src/api/run_json.cpp expects.
# 3.  exec amflow_cli "$@".

set -euo pipefail

FERMAT_DEFAULT=/usr/share/Ferl7/fer64

print_help() {
    cat <<'EOF'
AMFlow.cpp Docker image — C++17 multi-loop Feynman integral computation.

This image bundles FLINT 3.4 + Kira 3.1 + FireFly + amflow_cli.
Fermat (closed-source, free academic use, non-redistributable) is the
only piece you must provide at runtime via bind-mount.

Usage:
  docker run \
      -v /path/to/Ferl7:/usr/share/Ferl7:ro \
      -v "$(pwd)":/work \
      ghcr.io/chang18/amflow-cpp:latest \
      /opt/amflow-cpp/tools/bench/cutbubble_1L_eps001_black_box_amflow_cpp.json \
      /work/out.json

  - The first -v mounts Fermat (whose absolute path inside the image
    matches the hard-coded default that amflow_cli expects).
  - The second -v exposes your working directory as /work; outputs go
    there, and any bench JSON you pass by relative path is resolved
    from /work.
  - The 20 bundled oracle benches live at /opt/amflow-cpp/tools/bench/
    inside the image — reference them by absolute path as shown above.
    Run with 'ls /opt/amflow-cpp/tools/bench' to list them.

For full documentation:
  https://github.com/chang18/amflow-cpp/blob/main/docs/DOCKER.md

Fermat home page:
  https://home.bway.net/lewis/
EOF
}

# --- 1. Help bypass (no Fermat check) ---
case "${1:-}" in
    ""|-h|--help|help)
        print_help
        exit 0
        ;;
esac

# --- 2. Fermat detection ---
FERMAT="${FERMATPATH:-$FERMAT_DEFAULT}"
if [ ! -x "$FERMAT" ]; then
    cat >&2 <<EOF
ERROR: Fermat executable not found at "$FERMAT".

AMFlow.cpp's IBP backend (Kira) needs Fermat for multivariate algebra.
Fermat is closed-source (free for academic use, not redistributable),
so this image does NOT bundle it.  Mount your local install at runtime:

  docker run \\
      -v /path/to/Ferl7:/usr/share/Ferl7:ro \\
      -v "\$(pwd)":/work \\
      <this-image> input.json output.json

Or override the lookup path via FERMATPATH if your mount point differs:

  docker run \\
      -v /custom/path/Ferl7:/opt/fermat:ro \\
      -e FERMATPATH=/opt/fermat/fer64 \\
      -v "\$(pwd)":/work \\
      <this-image> input.json output.json

Fermat: https://home.bway.net/lewis/
Run with 'help' (or no args) for full usage.
EOF
    exit 2
fi

# --- 3. Pass through to amflow_cli ---
exec /usr/local/bin/amflow_cli "$@"
