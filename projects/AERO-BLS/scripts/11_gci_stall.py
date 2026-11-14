#!/usr/bin/env python3
"""Three-level grid-convergence (GCI) analysis at stall, NACA 0012, alpha=15, Re=1e6.

Celik et al. (2008) procedure: observed order p from the three solutions,
Richardson extrapolation to the continuum value, grid convergence indices.
"""
import csv
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

CELLS = {"coarse": 16100, "medium": 33740, "fine": 51360}
CASES = {
    "coarse": "naca0012_aoa15_coarse",
    "medium": "naca0012_aoa15",
    "fine": "naca0012_aoa15_fine",
}
FS = 1.25


def load(case):
    with open(os.path.join(ROOT, "data", "results", "summary.csv")) as f:
        for r in csv.DictReader(f):
            if r["case"] == case:
                return float(r["Cl"]), float(r["Cd"])
    raise SystemExit(f"case {case} not in summary.csv")


def gci_trio(phi_fine, phi_med, phi_coarse, n_fine, n_med, n_coarse):
    r21 = (n_fine / n_med) ** 0.5
    r32 = (n_med / n_coarse) ** 0.5
    e21 = phi_med - phi_fine
    e32 = phi_coarse - phi_med
    ratio = e32 / e21 if e21 != 0 else float("nan")
    if e21 == 0 or e32 == 0 or ratio <= 0:
        return {"r21": r21, "r32": r32, "eps21": e21, "eps32": e32,
                "ratio": ratio, "p": None, "oscillatory": ratio < 0}
    p = math.log(abs(ratio)) / math.log(r21)
    phi_ext = (r21 ** p * phi_fine - phi_med) / (r21 ** p - 1)
    gci21 = FS * abs(e21 / phi_fine) / (r21 ** p - 1)
    gci32 = FS * abs(e32 / phi_med) / (r32 ** p - 1)
    return {"r21": r21, "r32": r32, "eps21": e21, "eps32": e32,
            "ratio": ratio, "p": p, "phi_ext": phi_ext,
            "gci21": gci21, "gci32": gci32, "oscillatory": False}


def main():
    out = []
    for qty in ("Cl", "Cd"):
        vals = {}
        for lvl, case in CASES.items():
            cl, cd = load(case)
            vals[lvl] = cl if qty == "Cl" else cd
        res = gci_trio(vals["fine"], vals["medium"], vals["coarse"],
                       CELLS["fine"], CELLS["medium"], CELLS["coarse"])
        out.append((qty, vals, res))
    for qty, vals, res in out:
        print(f"--- {qty} at alpha=15 (NACA 0012, Re=1e6) ---")
        print(f"  coarse  ({CELLS['coarse']:>6} cells): {vals['coarse']:.4f}")
        print(f"  medium  ({CELLS['medium']:>6} cells): {vals['medium']:.4f}")
        print(f"  fine    ({CELLS['fine']:>6} cells): {vals['fine']:.4f}")
        print(f"  r32={res['r32']:.3f}  r21={res['r21']:.3f}")
        print(f"  eps32={res['eps32']:+.4f}  eps21={res['eps21']:+.4f}  "
              f"ratio={res['ratio']:+.3f}")
        if res["p"] is None:
            print(f"  NON-MONOTONE (ratio<0): observed order undefined; "
                  f"family oscillatory or non-asymptotic at stall")
        else:
            print(f"  observed order p={res['p']:.2f}")
            print(f"  Richardson phi_ext={res['phi_ext']:.4f}")
            print(f"  GCI_med={res['gci32']*100:.1f}%  GCI_fine={res['gci21']*100:.1f}%")


if __name__ == "__main__":
    main()
