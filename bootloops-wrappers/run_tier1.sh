#!/usr/bin/env bash
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code
# written by Claude (Anthropic) under his supervision.
# Tier 1: re-run upstream's 228 oracle bench cases of amflow-cpp and compare vs
# upstream's committed AMFlow-package reference values. On 006a275 in our environment
# (2026-06-11) 224/228 passed: three pentabox_2L cases abort in the eta->0
# residue step and wzbox_2L under-resolves at the committed truncation orders.
#
# Usage: run_tier1.sh [jobs]   (default 20 concurrent cases; configs use n_thread=4 for Kira)
set -u
ROOT="${AMFLOW_PORT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
REPO=${AMFLOW_REPO:-$ROOT/..}
CLI=$REPO/build/src/cli/amflow_cli
BENCH=$REPO/tools/bench
OUT=${TIER1_OUT:-$ROOT/results/tier1}
JOBS=${1:-20}
mkdir -p "$OUT/out" "$OUT/cmp" "$OUT/work" "$OUT/time"

run_case() {
    name=$1
    ROOT="${AMFLOW_PORT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
    REPO=${AMFLOW_REPO:-$ROOT/..}
    CLI=$REPO/build/src/cli/amflow_cli
    BENCH=$REPO/tools/bench
    OUT=${TIER1_OUT:-$ROOT/results/tier1}
    cfg=$BENCH/${name}_black_box_amflow_cpp.json
    [ -f "$cfg" ] || cfg=$BENCH/${name}_cpp.json
    ref=$BENCH/${name}_mma_reference.json
    # rewrite work_dir to a private one (configs hardcode /tmp paths)
    python3 - "$cfg" "$OUT/work/$name" > "$OUT/work/${name}.json" <<'EOF'
import json, sys
cfg = json.load(open(sys.argv[1]))
cfg["work_dir"] = sys.argv[2]
json.dump(cfg, sys.stdout)
EOF
    /usr/bin/time -v -o "$OUT/time/${name}.time" \
        "$CLI" "$OUT/work/${name}.json" "$OUT/out/${name}.json" \
        > "$OUT/out/${name}.log" 2>&1
    rc=$?
    if [ $rc -ne 0 ]; then
        echo "$name RUN_FAIL rc=$rc"
        return
    fi
    if [ -f "$ref" ]; then
        python3 "$BENCH/compare_sampled.py" "$ref" "$OUT/out/${name}.json" \
            --json-out "$OUT/cmp/${name}.cmp.json" > "$OUT/cmp/${name}.cmp.txt" 2>&1
        crc=$?
        wall=$(grep "Elapsed (wall" "$OUT/time/${name}.time" | awk '{print $NF}')
        rss=$(grep "Maximum resident" "$OUT/time/${name}.time" | awk '{print $NF}')
        if [ $crc -eq 0 ]; then echo "$name PASS wall=$wall rss_kb=$rss";
        else echo "$name CMP_FAIL wall=$wall rss_kb=$rss"; fi
    else
        echo "$name NO_REF"
    fi
}
export -f run_case

# resume-aware: skip cases already recorded in summary.txt
touch "$OUT/summary.txt"
ls "$BENCH"/*_mma_reference.json | sed 's/.*\///;s/_mma_reference\.json//' | sort | \
    comm -23 - <(awk '{print $1}' "$OUT/summary.txt" | sort) | \
    parallel -j "$JOBS" run_case {} | tee -a "$OUT/summary.txt"
echo "DONE: $(grep -c ' PASS ' "$OUT/summary.txt") pass / $(grep -c 'FAIL' "$OUT/summary.txt") fail"
