#!/usr/bin/env bash
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code
# written by Claude (Anthropic) under his supervision.
# Tier 5: cost-vs-digits and cost-vs-eps-depth scaling, run SEQUENTIALLY on a
# quiet machine (patched build). Configs from gen_tier23_configs-era t5*.json.
set -u
ROOT="${AMFLOW_PORT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
REPO=${AMFLOW_REPO:-$ROOT/..}
CLI=$REPO/build/src/cli/amflow_cli
CFG=$ROOT/results/tier23/configs
OUT=$ROOT/results/tier5
mkdir -p "$OUT/out" "$OUT/time" "$OUT/work"

for cfg in "$CFG"/t5dig_box0_g30_o8.json "$CFG"/t5dig_box0_g100_o8.json "$CFG"/t5dig_box0_g300_o8.json \
           "$CFG"/t5dig_sunrise_g30_o8.json "$CFG"/t5dig_sunrise_g100_o8.json "$CFG"/t5dig_sunrise_g300_o8.json \
           "$CFG"/t5depth_sunrise_g60_o4.json "$CFG"/t5depth_sunrise_g60_o6.json "$CFG"/t5depth_sunrise_g60_o8.json "$CFG"/t5depth_sunrise_g60_o10.json; do
    name=$(basename "$cfg" .json)
    [ -s "$OUT/out/${name}.json" ] && { echo "$name SKIP"; continue; }
    /usr/bin/time -v -o "$OUT/time/${name}.time" \
        "$CLI" "$cfg" "$OUT/out/${name}.json" > "$OUT/out/${name}.log" 2>&1
    rc=$?
    wall=$(grep "Elapsed (wall" "$OUT/time/${name}.time" | awk '{print $NF}')
    rss=$(grep "Maximum resident" "$OUT/time/${name}.time" | awk '{print $NF}')
    echo "$name rc=$rc wall=$wall rss_kb=$rss"
done | tee -a "$OUT/runlog.txt"

# scaling fits
python3 - <<'EOF'
import re, math, os
def secs(w):
    p=[float(x) for x in w.split(':')]
    return p[0]*60+p[1] if len(p)==2 else p[0]*3600+p[1]*60+p[2]
rows={}
for line in open(os.environ.get('AMFLOW_PORT_ROOT','.')+'/results/tier5/runlog.txt'):
    m=re.match(r'(t5\w+) rc=0 wall=([\d:.]+) rss_kb=(\d+)', line)
    if m: rows[m.group(1)]=(secs(m.group(2)), int(m.group(3)))
def fit(tag, xs, names):
    pts=[(x, rows[n][0]) for x,n in zip(xs,names) if n in rows]
    if len(pts)<2: print(f"{tag}: insufficient data"); return
    import statistics
    lx=[math.log(p[0]) for p in pts]; ly=[math.log(p[1]) for p in pts]
    n=len(pts); b=(n*sum(a*c for a,c in zip(lx,ly))-sum(lx)*sum(ly))/(n*sum(a*a for a in lx)-sum(lx)**2)
    print(f"{tag}: exponent ~ {b:.2f}   points: {[(p[0], round(p[1],1)) for p in pts]}")
fit("box0 time vs digits", [30,100,300], ['t5dig_box0_g30_o8','t5dig_box0_g100_o8','t5dig_box0_g300_o8'])
fit("sunrise time vs digits", [30,100,300], ['t5dig_sunrise_g30_o8','t5dig_sunrise_g100_o8','t5dig_sunrise_g300_o8'])
fit("sunrise time vs eps_order", [4,6,8,10], [f't5depth_sunrise_g60_o{o}' for o in (4,6,8,10)])
for k,(t,r) in sorted(rows.items()): print(f"  {k}: {t:.1f}s rss={r/1024:.0f}MB")
EOF
