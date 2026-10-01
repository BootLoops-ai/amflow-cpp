#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code
# written by Claude (Anthropic) under his supervision.
"""Compare amflow-cpp solve_integrals output against analytic anchors,
per epsilon order, reporting TRUE error vs CLAIMED ball radius.

Usage:
  compare_orders.py <engine_out.json> <anchor.json> <anchor_case_key> [--loops L]

Anchor formats supported:
  - ref_closed_forms.py output: anchor["cases"][key]["orders"][order] = {re, im}
  - ref_sunrise_bessel.py output: anchor["orders"][order] = {re, im}

Output: one line per eps order with
  true_err  (|engine_mid - anchor|),
  radius    (engine ball radius, max of re/im),
  honest    (true_err <= radius ?),
  digits    (correct decimal digits = -log10(true rel err)).
"""
import json, sys, re, argparse
from mpmath import mp, mpf, mpc, mpmathify

mp.dps = 400


def parse_ball(s):
    """'[1.23 +/- 4.5e-31]' or plain '1.23' -> (mid, rad)."""
    s = str(s).strip()
    m = re.match(r"\[(.*?)\s*\+/-\s*(.*?)\]", s)
    if m:
        return mpmathify(m.group(1)), mpmathify(m.group(2))
    return mpmathify(s), mpf(0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("engine_out")
    ap.add_argument("anchor")
    ap.add_argument("anchor_key", nargs="?", default=None)
    ap.add_argument("--target", type=int, default=0, help="result index")
    ap.add_argument("--re-only", action="store_true",
                    help="compare real parts only (anchor Im untrusted; "
                         "integral known real). Engine Im is still reported "
                         "as |Im| in the true_err column context.")
    args = ap.parse_args()

    eng = json.load(open(args.engine_out))
    anc = json.load(open(args.anchor))
    if "cases" in anc:
        orders_ref = anc["cases"][args.anchor_key]["orders"]
    else:
        orders_ref = anc["orders"]

    res = eng["result"][args.target]
    lead = res["leading_order"]
    rows = []
    for c in res["coefficients"]:
        k = str(c["order"])
        if k not in orders_ref:
            rows.append((c["order"], None))
            continue
        mid_re, rad_re = parse_ball(c["value"]["re"])
        mid_im, rad_im = parse_ball(c["value"]["im"])
        ref = mpc(mpmathify(orders_ref[k]["re"]), mpmathify(orders_ref[k]["im"]))
        if args.re_only:
            err = abs(mid_re - mp.re(ref)) + abs(mid_im)  # engine Im must be ~0
        else:
            err = abs(mpc(mid_re, mid_im) - ref)
        rad = max(rad_re, rad_im)
        scale = abs(ref) if abs(ref) > 0 else mpf(1)
        digits = float(-mp.log10(err / scale)) if err > 0 else 999.0
        honest = bool(err <= rad * mp.sqrt(2) + mpf(10) ** (-mp.dps + 20))
        rows.append((c["order"], dict(
            true_err=mp.nstr(err, 3), radius=mp.nstr(rad, 3),
            honest=honest, correct_digits=round(digits, 1))))

    print(f"{'order':>6} {'true_err':>12} {'radius':>12} {'honest':>7} {'digits':>8}")
    out = {"rows": []}
    for o, r in rows:
        if r is None:
            print(f"{o:>6} (no anchor)")
            continue
        print(f"{o:>6} {r['true_err']:>12} {r['radius']:>12} {str(r['honest']):>7} {r['correct_digits']:>8}")
        out["rows"].append({"order": o, **r})
    jpath = args.engine_out.replace(".json", ".orders_cmp.json")
    json.dump(out, open(jpath, "w"), indent=1)


if __name__ == "__main__":
    main()
