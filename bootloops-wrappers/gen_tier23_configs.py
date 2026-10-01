#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code
# written by Claude (Anthropic) under his supervision.
"""Generate Tier-2/3 solve_integrals configs for amflow-cpp.

Cases (all rational Euclidean kinematics unless noted):
  bubble0  massless bubble        s=-1            1L  anchor: Gamma closed form
  box0     massless on-shell box  s=-3, t=-7      1L  anchor: 2F1 closed form
  bubblem  equal-mass bubble      s=-1 / s=5, m=1 1L  anchor: 1D log-moments
                                  (s=5 is ABOVE threshold: Im-part sign check)
  sunrise  equal-mass sunrise     p^2=-1, m^2=1   2L  anchor: Bessel moments (elliptic)

Variants: goal_digits in {30,100,300}; 'auto' = only goal_digits set (tests
automatic precision propagation), 'manual' = working_pre/x_order raised by hand.
"""
import json, os, sys, itertools

OUT = os.environ.get("TIER23_CONFIGS", "results/tier23/configs")
WORK = os.environ.get("TIER23_WORK", "results/tier23/work")
os.makedirs(OUT, exist_ok=True)

FAMILIES = {
    "bubble0": {
        "family": {
            "name": "bubble0", "loops": ["l"], "legs": ["p"],
            "replacement": {"p^2": "s"},
            "propagators": ["l^2", "(l - p)^2"],
        },
        "numeric": {"s": "-1"}, "targets": [[1, 1]], "loops": 1,
    },
    "box0": {
        "family": {
            "name": "box0", "loops": ["l"], "legs": ["p1", "p2", "p3", "p4"],
            "conservation": {"p4": "-p1 - p2 - p3"},
            "replacement": {"p1^2": "0", "p2^2": "0", "p3^2": "0", "p4^2": "0",
                            "(p1 + p2)^2": "s", "(p1 + p3)^2": "t"},
            "propagators": ["l^2", "(l + p1)^2", "(l + p1 + p2)^2",
                            "(l + p1 + p2 + p4)^2"],
        },
        "numeric": {"s": "-3", "t": "-7"}, "targets": [[1, 1, 1, 1]], "loops": 1,
    },
    "bubblem_eu": {
        "family": {
            "name": "bubblem", "loops": ["l"], "legs": ["p"],
            "replacement": {"p^2": "s"},
            "propagators": ["l^2 - msq", "(l - p)^2 - msq"],
        },
        "numeric": {"s": "-1", "msq": "1"}, "targets": [[1, 1]], "loops": 1,
    },
    "bubblem_phys": {
        "family": {
            "name": "bubblem", "loops": ["l"], "legs": ["p"],
            "replacement": {"p^2": "s"},
            "propagators": ["l^2 - msq", "(l - p)^2 - msq"],
        },
        "numeric": {"s": "5", "msq": "1"}, "targets": [[1, 1]], "loops": 1,
    },
    "sunrise": {
        "family": {
            "name": "sunrise", "loops": ["l1", "l2"], "legs": ["p"],
            "conservation": {},
            "replacement": {"p^2": "s"},
            "propagators": ["l1^2 - msq", "l2^2 - msq", "(l1 + l2 - p)^2 - msq",
                            "(l1 + p)^2", "(l2 + p)^2"],
        },
        "numeric": {"s": "-1", "msq": "1"}, "targets": [[1, 1, 1, 0, 0]], "loops": 2,
    },
}

# eps_order: highest fitted power = eps_order - 2L. We want through eps^4
# minimum (eps^6 on 1L cases for depth scaling).
PLAN = [
    # (case, goal, eps_order, variant, extra_options)
    *[("bubble0", g, 8, "auto", {}) for g in (30, 100, 300)],
    *[("box0", g, 8, "auto", {}) for g in (30, 100, 300)],
    ("bubblem_eu", 100, 8, "auto", {}),
    ("bubblem_phys", 100, 8, "auto", {}),
    *[("sunrise", g, 8, "auto", {}) for g in (30, 100, 300)],
    # manual variants: raise dials in case 'auto' saturates
    ("bubble0", 300, 8, "manual",
     {"working_pre": 450, "x_order": 1200, "extra_x_order": 200, "rationalize_pre": 450}),
    ("box0", 300, 8, "manual",
     {"working_pre": 450, "x_order": 1200, "extra_x_order": 200, "rationalize_pre": 450}),
    ("sunrise", 300, 8, "manual",
     {"working_pre": 450, "x_order": 1200, "extra_x_order": 200, "rationalize_pre": 450}),
]


def make(case, goal, eps_order, variant, extra):
    f = FAMILIES[case]
    name = f"{case}_g{goal}_o{eps_order}_{variant}"
    cfg = {
        "mode": "solve_integrals",
        "options": {"silent_mode": True, **extra},
        "family": f["family"],
        "integrals": [{"indices": t} for t in f["targets"]],
        "goal_digits": goal,
        "eps_order": eps_order,
        "work_dir": f"{WORK}/{name}",
        "amf_options": {"blackbox": {"numeric_values": f["numeric"],
                                     "n_thread": 8}},
    }
    path = os.path.join(OUT, name + ".json")
    with open(path, "w") as fh:
        json.dump(cfg, fh, indent=1)
    return path


if __name__ == "__main__":
    for spec in PLAN:
        print(make(*spec))
