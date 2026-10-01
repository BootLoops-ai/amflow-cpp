#!/usr/bin/env bash
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code
# written by Claude (Anthropic) under his supervision.
# Tier 2/3: run the generated solve_integrals configs, with timing.
set -u
ROOT="${AMFLOW_PORT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
REPO=${AMFLOW_REPO:-$ROOT/..}
CLI=$REPO/build/src/cli/amflow_cli
CFG=$ROOT/results/tier23/configs
OUT=$ROOT/results/tier23
JOBS=${1:-4}
ONLY=${2:-}
mkdir -p "$OUT/out" "$OUT/time"

run_one() {
    cfg=$1
    ROOT="${AMFLOW_PORT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
    REPO=${AMFLOW_REPO:-$ROOT/..}
CLI=$REPO/build/src/cli/amflow_cli
    OUT=$ROOT/results/tier23
    name=$(basename "$cfg" .json)
    [ -s "$OUT/out/${name}.json" ] && { echo "$name SKIP(done)"; return; }
    /usr/bin/time -v -o "$OUT/time/${name}.time" \
        "$CLI" "$cfg" "$OUT/out/${name}.json" > "$OUT/out/${name}.log" 2>&1
    rc=$?
    wall=$(grep "Elapsed (wall" "$OUT/time/${name}.time" 2>/dev/null | awk '{print $NF}')
    rss=$(grep "Maximum resident" "$OUT/time/${name}.time" 2>/dev/null | awk '{print $NF}')
    echo "$name rc=$rc wall=$wall rss_kb=$rss"
}
export -f run_one

ls "$CFG"/*${ONLY}*.json | parallel -j "$JOBS" run_one {} | tee -a "$OUT/runlog.txt"
