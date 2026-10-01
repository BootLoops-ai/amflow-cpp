#!/usr/bin/env bash
# Run all committed C++ oracle benchmarks once each, recording wall-clock time.
# Used for the v1.0 performance audit.  Outputs CSV: name,seconds,exit_code.

set -u
cd "$(dirname "$0")/../.."

OUT_CSV="${1:-/tmp/amflow_perf_results.csv}"
echo "name,seconds,exit_code" > "$OUT_CSV"

CLI="./build/src/cli/amflow_cli"
[ -x "$CLI" ] || { echo "missing $CLI; run cmake --build build -j first" >&2; exit 2; }

# Each bench gets a fresh tmp work_dir so Kira-cache state doesn't bias
# adjacent runs.  We respect the work_dir baked into the json input
# unless that is missing (which it shouldn't be for committed benches).

BENCHES=(
  tools/bench/doublebox_4mass_2L_eps001_black_box_amflow_cpp.json
  tools/bench/xbox_2L_3mass_eps001_black_box_amflow_cpp.json
  tools/bench/vtx2_2L_4mass_eps001_black_box_amflow_cpp.json
  tools/bench/sunset_2L_4mass_eps001_black_box_amflow_cpp.json
  tools/bench/box_1L_4mass_alt_eps001_black_box_amflow_cpp.json
  tools/bench/triangle_1L_2mass_alt_eps001_black_box_amflow_cpp.json
  tools/bench/banana_3loop_eps001_black_box_amflow_cpp.json
  tools/bench/banana_4L_2mass_eps001_black_box_amflow_cpp.json
  tools/bench/banana_4L_mixed_black_box_amflow_cpp.json
  tools/bench/banana_4loop_eps001_black_box_amflow_cpp.json
  tools/bench/bn3_4mass_3L_eps001_black_box_amflow_cpp.json
  tools/bench/hexagon_1L_2mass_eps001_black_box_amflow_cpp.json
  tools/bench/hexagon_1L_3mass_eps001_black_box_amflow_cpp.json
  tools/bench/hexagon_1L_4mass_eps001_black_box_amflow_cpp.json
  tools/bench/hexagon_1L_5mass_eps001_black_box_amflow_cpp.json
  tools/bench/hexagon_1L_6mass_eps001_black_box_amflow_cpp.json
  tools/bench/pentagon_1L_2mass_eps001_black_box_amflow_cpp.json
  tools/bench/pentagon_1L_4mass_eps001_black_box_amflow_cpp.json
  tools/bench/pentagon_1L_5mass_eps001_black_box_amflow_cpp.json
  tools/bench/bn3mix_eps001_black_box_amflow_cpp.json
  tools/bench/box_1L_2mass_eps001_black_box_amflow_cpp.json
  tools/bench/box_1L_3mass_eps001_black_box_amflow_cpp.json
  tools/bench/box_1L_4mass_eps001_black_box_amflow_cpp.json
  tools/bench/box1_d0_7_3_solve_integrals_cpp.json
  tools/bench/box1_multitarget_black_box_amflow_cpp.json
  tools/bench/box1_single_solve_integrals_cpp.json
  tools/bench/bubble_1L_largemass_eps001_black_box_amflow_cpp.json
  tools/bench/bubble_1L_diffmass_black_box_amflow_cpp.json
  tools/bench/cutbanana_3L_eps001_black_box_amflow_cpp.json
  tools/bench/cutbanana_4L_eps001_black_box_amflow_cpp.json
  tools/bench/cutbubble_1L_eps001_black_box_amflow_cpp.json
  tools/bench/cutbubble_1L_eps2_black_box_amflow_cpp.json
  tools/bench/cutbubble_1L_eps10000_black_box_amflow_cpp.json
  tools/bench/cutsunrise_2L_eps001_black_box_amflow_cpp.json
  tools/bench/dotted_3L_banana_black_box_amflow_cpp.json
  tools/bench/doublebox_sv_eps001_black_box_amflow_cpp.json
  tools/bench/ewbox_1loop_eps001_black_box_amflow_cpp.json
  tools/bench/sunset_2L_3mass_eps001_black_box_amflow_cpp.json
  tools/bench/sunset_2L_2mass_eps001_black_box_amflow_cpp.json
  tools/bench/sunrise_2L_2mass_eps001_black_box_amflow_cpp.json
  tools/bench/sunrise_2L_threshold_black_box_amflow_cpp.json
  tools/bench/sunset_2L_onshell_black_box_amflow_cpp.json
  tools/bench/tadbubble_2L_black_box_amflow_cpp.json
  tools/bench/triangle_1L_2mass_black_box_amflow_cpp.json
  tools/bench/triangle_1L_3mass_alt_eps001_black_box_amflow_cpp.json
  tools/bench/triangle_1L_3mass_black_box_amflow_cpp.json
  tools/bench/tt_2loop_box_eps001_black_box_amflow_cpp.json
  tools/bench/tt_cutkosky_probe_eps001_black_box_amflow_cpp.json
  tools/bench/tradcut_phase_2L_eps001_black_box_amflow_cpp.json
  tools/bench/tt_higher_rank_eps001_black_box_amflow_cpp.json
  tools/bench/vacuum_2L_2mass_eps001_black_box_amflow_cpp.json
  tools/bench/vacuum_2L_sunrise_black_box_amflow_cpp.json
  tools/bench/vtx2_2L_1mass_eps001_black_box_amflow_cpp.json
  tools/bench/vtx2_2L_2mass_eps001_black_box_amflow_cpp.json
  tools/bench/vtx2_2loop_vertex_eps001_black_box_amflow_cpp.json
  tools/bench/vtx2_2loop_vertex_masters_only_eps001_cpp.json
  tools/bench/xbox_2L_2mass_eps001_black_box_amflow_cpp.json
  tools/bench/xbox_2L_eps10_black_box_amflow_cpp.json
  tools/bench/xbox_2loop_eps001_black_box_amflow_cpp.json
  tools/bench/hexagon_1L_2mass_eps01_black_box_amflow_cpp.json
  tools/bench/hexagon_1L_3mass_eps01_black_box_amflow_cpp.json
  tools/bench/doublebox_3mass_alt_2L_eps001_black_box_amflow_cpp.json
  tools/bench/pentagon_1L_2mass_eps01_black_box_amflow_cpp.json
  tools/bench/triangle_1L_2mass_eps01_black_box_amflow_cpp.json
  tools/bench/bubble_1L_diffmass_eps01_black_box_amflow_cpp.json
  tools/bench/hexagon_1L_5mass_eps01_black_box_amflow_cpp.json
  tools/bench/box_1L_4mass_eps10_black_box_amflow_cpp.json
  tools/bench/pentagon_1L_5mass_eps01_black_box_amflow_cpp.json
  tools/bench/triangle_1L_1mass_alt_eps001_black_box_amflow_cpp.json
  tools/bench/box_1L_1mass_alt_eps001_black_box_amflow_cpp.json
  tools/bench/pentagon_1L_1mass_eps001_black_box_amflow_cpp.json
  tools/bench/hexagon_1L_1mass_eps001_black_box_amflow_cpp.json
  tools/bench/sunset_2L_2mass_eps01_black_box_amflow_cpp.json
  tools/bench/doublebox_2L_2mass_eps01_black_box_amflow_cpp.json
  tools/bench/vtx2_2L_2mass_eps01_black_box_amflow_cpp.json
  tools/bench/pentabox_2L_2mass_eps001_black_box_amflow_cpp.json
  tools/bench/pentabox_2L_1mass_eps001_black_box_amflow_cpp.json
  tools/bench/sunset_2L_3mass_eps01_black_box_amflow_cpp.json
  tools/bench/sunset_2L_4mass_eps01_black_box_amflow_cpp.json
  tools/bench/vtx2_2L_3mass_eps01_black_box_amflow_cpp.json
  tools/bench/vtx2_2L_4mass_eps01_black_box_amflow_cpp.json
  tools/bench/xbox_2L_3mass_eps01_black_box_amflow_cpp.json
  tools/bench/mercedes_3L_1mass_eps01_black_box_amflow_cpp.json
  tools/bench/mercedes_3L_2mass_eps01_black_box_amflow_cpp.json
  tools/bench/bn3_1mass_3L_eps001_black_box_amflow_cpp.json
  tools/bench/bn3_2mass_3L_eps001_black_box_amflow_cpp.json
  tools/bench/bn3mix_3mass_3L_eps001_black_box_amflow_cpp.json
  tools/bench/banana_3L_2mass_eps001_black_box_amflow_cpp.json
  tools/bench/banana_3L_eqmass_eps01_black_box_amflow_cpp.json
  tools/bench/banana_4L_1mass_eps001_black_box_amflow_cpp.json
  tools/bench/banana_4L_2mass_alt_eps001_black_box_amflow_cpp.json
  tools/bench/banana_4L_2mass_eps01_black_box_amflow_cpp.json
  tools/bench/banana_4L_1mass_eps01_black_box_amflow_cpp.json
  tools/bench/bn3mix_eqmass_3L_eps001_black_box_amflow_cpp.json
  tools/bench/xbox_2L_4mass_eps001_black_box_amflow_cpp.json
  tools/bench/banana_3L_eqmass_eps0001_black_box_amflow_cpp.json
  tools/bench/banana_3L_1mass_eps01_black_box_amflow_cpp.json
  tools/bench/bn3_2mass_alt_3L_eps001_black_box_amflow_cpp.json
  tools/bench/mercedes_3L_3mass_eps01_black_box_amflow_cpp.json
  tools/bench/doublebox_2L_3mass_eps01_black_box_amflow_cpp.json
  tools/bench/banana_4L_2mass_alt2_eps001_black_box_amflow_cpp.json
  tools/bench/mercedes_3L_2mass_alt_eps001_black_box_amflow_cpp.json
  tools/bench/mercedes_3L_1mass_eps0001_black_box_amflow_cpp.json
  tools/bench/mercedes_3L_massless_eps001_black_box_amflow_cpp.json
  tools/bench/banana_3L_eqmass_psq4_eps001_black_box_amflow_cpp.json
  tools/bench/banana_3L_2mass_eps01_black_box_amflow_cpp.json
  tools/bench/bn3mix_4mass_3L_eps001_black_box_amflow_cpp.json
  tools/bench/banana_3L_2mass_alt_eps001_black_box_amflow_cpp.json
  tools/bench/bn3_2mass_eps01_3L_black_box_amflow_cpp.json
  tools/bench/banana_3L_eqmass_psqQuarter_eps001_black_box_amflow_cpp.json
  tools/bench/bn3mix_2mass_eps01_3L_black_box_amflow_cpp.json
  tools/bench/bn3_1mass_eps01_3L_black_box_amflow_cpp.json
  tools/bench/mercedes_3L_massless_eps01_black_box_amflow_cpp.json
  tools/bench/banana_3L_4mass_eps001_black_box_amflow_cpp.json
  tools/bench/banana_3L_2mass_eps01_v2_black_box_amflow_cpp.json
  tools/bench/banana_3L_eqmass_psq10_eps001_black_box_amflow_cpp.json
  tools/bench/bn3_2mass_psq4_3L_eps001_black_box_amflow_cpp.json
  tools/bench/banana_3L_eqmass_psq2_eps001_black_box_amflow_cpp.json
  tools/bench/mercedes_3L_2mass_alt_eps01_black_box_amflow_cpp.json
  tools/bench/banana_4L_2mass_eps0001_black_box_amflow_cpp.json
  tools/bench/mercedes_3L_eqmass_3main_eps001_black_box_amflow_cpp.json
  tools/bench/banana_3L_2mass_psqQuarter_eps001_black_box_amflow_cpp.json
  tools/bench/doublebox_2L_4mass_eps01_black_box_amflow_cpp.json
  tools/bench/mercedes_3L_2mass_eps0001_black_box_amflow_cpp.json
  tools/bench/pentabox_2L_3mass_eps001_black_box_amflow_cpp.json
  tools/bench/banana_4L_2mass_psq4_eps001_black_box_amflow_cpp.json
  tools/bench/pentabox_2L_1mass_eps01_black_box_amflow_cpp.json
  tools/bench/mercedes_3L_3mass_eps0001_black_box_amflow_cpp.json
  tools/bench/mercedes_3L_3mass_psq4_eps001_black_box_amflow_cpp.json
  tools/bench/bn3_4mass_eps01_3L_black_box_amflow_cpp.json
  tools/bench/mercedes_3L_2mass_psq2_eps001_black_box_amflow_cpp.json
  tools/bench/banana_4L_2mass_psq1Quarter_eps001_black_box_amflow_cpp.json
  tools/bench/banana_4L_2mass_psq1Half_eps001_black_box_amflow_cpp.json
  tools/bench/mercedes_3L_eqmass_psq2_eps001_black_box_amflow_cpp.json
  tools/bench/banana_3L_eqmass_psq6_eps001_black_box_amflow_cpp.json
  tools/bench/mercedes_3L_1mass_psq8_eps001_black_box_amflow_cpp.json
  tools/bench/banana_3L_eqmass_psq8_eps001_black_box_amflow_cpp.json
  tools/bench/mercedes_3L_1mass_psq2_eps001_black_box_amflow_cpp.json
  tools/bench/mercedes_3L_1mass_psq4_eps001_black_box_amflow_cpp.json
  tools/bench/mercedes_3L_eqmass_psq4_eps001_black_box_amflow_cpp.json
  tools/bench/banana_4L_2mass_psq2_eps001_black_box_amflow_cpp.json
  tools/bench/banana_4L_1mass_psq4_eps001_black_box_amflow_cpp.json
  tools/bench/banana_4L_1mass_psq1Half_eps001_black_box_amflow_cpp.json
  tools/bench/mercedes_3L_3mass_psq8_eps001_black_box_amflow_cpp.json
  tools/bench/banana_3L_2mass_psq2_eps001_black_box_amflow_cpp.json
  tools/bench/mercedes_3L_3mass_psq2_eps001_black_box_amflow_cpp.json
  tools/bench/banana_3L_1mass_psq4_eps001_black_box_amflow_cpp.json
  tools/bench/banana_3L_2mass_psq8_eps001_black_box_amflow_cpp.json
  tools/bench/banana_4L_2mass_psq6_eps001_black_box_amflow_cpp.json
  tools/bench/mercedes_3L_2mass_psq4_eps001_black_box_amflow_cpp.json
  tools/bench/banana_3L_eqmass_psq5_eps001_black_box_amflow_cpp.json
  tools/bench/banana_4L_1mass_psq2_eps001_black_box_amflow_cpp.json
  tools/bench/banana_3L_1mass_psq2_eps001_black_box_amflow_cpp.json
  tools/bench/banana_3L_eqmass_psq3_eps001_black_box_amflow_cpp.json
  tools/bench/bn3mix_2mass_alt_3L_eps001_black_box_amflow_cpp.json
  tools/bench/banana_4L_1mass_psq8_eps001_black_box_amflow_cpp.json
  tools/bench/banana_4L_1mass_psq16_eps001_black_box_amflow_cpp.json
  tools/bench/banana_4L_2mass_alt2_psq2_eps001_black_box_amflow_cpp.json
  tools/bench/banana_4L_2mass_alt_psq4_eps001_black_box_amflow_cpp.json
  tools/bench/sunset_bubble_4L_2leg_eqmass_eps001_black_box_amflow_cpp.json
  tools/bench/sunset_bubble_4L_1mass_l1_eps001_black_box_amflow_cpp.json
  tools/bench/sunset_bubble_4L_1mass_l4_eps001_black_box_amflow_cpp.json
  tools/bench/sunset_bubble_4L_1mass_top_eps001_black_box_amflow_cpp.json
  tools/bench/sunset_bubble_4L_2mass_bub_core_eps001_black_box_amflow_cpp.json
  tools/bench/sunset_bubble_4L_2leg_alt_eqmass_eps001_black_box_amflow_cpp.json
  tools/bench/sunset_bubble_4L_alt1_eqmass_eps001_black_box_amflow_cpp.json
  tools/bench/cutbanana_4L_1mass_eps001_black_box_amflow_cpp.json
  tools/bench/cutbanana_4L_2mass_eps001_black_box_amflow_cpp.json
  tools/bench/banana_4loop_sub1_eps001_black_box_amflow_cpp.json
  tools/bench/banana_4loop_sub2_eps001_black_box_amflow_cpp.json
  tools/bench/banana_4loop_sub3_eps001_black_box_amflow_cpp.json
  tools/bench/banana_4loop_sub4_eps001_black_box_amflow_cpp.json
  tools/bench/banana_4loop_sub5_eps001_black_box_amflow_cpp.json
  tools/bench/banana_4loop_sub6_eps001_black_box_amflow_cpp.json
  tools/bench/banana_4L_1mass_sub1_eps001_black_box_amflow_cpp.json
  tools/bench/banana_4L_1mass_sub2_eps001_black_box_amflow_cpp.json
  tools/bench/banana_4L_1mass_sub3_eps001_black_box_amflow_cpp.json
  tools/bench/banana_4L_1mass_sub4_eps001_black_box_amflow_cpp.json
  tools/bench/banana_4L_1mass_sub5_eps001_black_box_amflow_cpp.json
  tools/bench/banana_4L_1mass_sub6_eps001_black_box_amflow_cpp.json
  tools/bench/banana_4L_2mass_sub1_eps001_black_box_amflow_cpp.json
  tools/bench/banana_4L_2mass_sub2_eps001_black_box_amflow_cpp.json
  tools/bench/banana_4L_2mass_sub3_eps001_black_box_amflow_cpp.json
  tools/bench/banana_4L_2mass_sub4_eps001_black_box_amflow_cpp.json
  tools/bench/banana_4L_2mass_sub5_eps001_black_box_amflow_cpp.json
  tools/bench/banana_4L_2mass_sub6_eps001_black_box_amflow_cpp.json
  tools/bench/banana_4L_mixed_sub1_eps001_black_box_amflow_cpp.json
  tools/bench/banana_4L_mixed_sub2_eps001_black_box_amflow_cpp.json
  tools/bench/banana_4L_mixed_sub3_eps001_black_box_amflow_cpp.json
  tools/bench/banana_4L_mixed_sub4_eps001_black_box_amflow_cpp.json
  tools/bench/banana_4L_mixed_sub5_eps001_black_box_amflow_cpp.json
  tools/bench/banana_4L_mixed_sub6_eps001_black_box_amflow_cpp.json
  tools/bench/cutbanana_4L_sub1_eps001_black_box_amflow_cpp.json
  tools/bench/cutbanana_4L_sub2_eps001_black_box_amflow_cpp.json
  tools/bench/cutbanana_4L_sub3_eps001_black_box_amflow_cpp.json
  tools/bench/sunset_bubble_4L_alt2_eqmass_eps001_black_box_amflow_cpp.json
  tools/bench/sunset_2L_1mass_eps001_black_box_amflow_cpp.json
  tools/bench/vtx2_2L_1mass_eps01_black_box_amflow_cpp.json
)

for input in "${BENCHES[@]}"; do
    name="$(basename "$input" .json)"
    [ -f "$input" ] || { echo "  skip (missing): $name" >&2; continue; }
    out="/tmp/amflow_perf_${name}_out.json"
    echo "[$(date +%H:%M:%S)] running: $name" >&2
    start=$(date +%s.%N)
    "$CLI" "$input" "$out" >/dev/null 2>"/tmp/amflow_perf_${name}_err.log"
    rc=$?
    end=$(date +%s.%N)
    secs=$(python3 -c "print(f'{$end - $start:.2f}')")
    echo "$name,$secs,$rc" >> "$OUT_CSV"
    echo "  -> $secs s (exit $rc)" >&2
done

echo "" >&2
echo "results written to $OUT_CSV" >&2
echo "" >&2
column -t -s',' "$OUT_CSV"
