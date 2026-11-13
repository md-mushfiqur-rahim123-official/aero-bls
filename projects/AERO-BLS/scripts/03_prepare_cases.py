#!/usr/bin/env python3
"""Prepare one OpenFOAM case from base_case for a given airfoil/AoA/Re.

Substitutes freestream velocity, turbulence inlet quantities, solver
directions, and kinematic viscosity into the templated case files.
"""
import argparse
import json
import math
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

NU_AIR = 1.5e-5      # kinematic viscosity of air [m^2/s]
CHORD = 1.0          # chord length [m]
TU = 0.001           # freestream turbulence intensity
LT = 0.01            # turbulence length scale [m] (0.01 c)
C_MU = 0.09
SPAN = 0.1           # 2D slab thickness [m]


def freestream_turbulence(u_inf):
    k_inf = 1.5 * (TU * u_inf) ** 2
    omega_inf = math.sqrt(k_inf) / (C_MU ** 0.25 * LT)
    nut_inf = k_inf / omega_inf
    return k_inf, omega_inf, nut_inf


def substitute(path, mapping):
    with open(path) as f:
        text = f.read()
    for key, val in mapping.items():
        text = text.replace(key, val)
    with open(path, "w") as f:
        f.write(text)


def prepare(airfoil, alpha, re, case_dir, mesh_path, end_time=4000):
    u_inf = re * NU_AIR / CHORD
    alpha_rad = math.radians(alpha)
    ux, uy = u_inf * math.cos(alpha_rad), u_inf * math.sin(alpha_rad)
    lift = (-math.sin(alpha_rad), math.cos(alpha_rad))
    drag = (math.cos(alpha_rad), math.sin(alpha_rad))
    k_inf, omega_inf, nut_inf = freestream_turbulence(u_inf)

    if os.path.exists(case_dir):
        shutil.rmtree(case_dir)
    shutil.copytree(os.path.join(ROOT, "openfoam_cases", "base_case"), case_dir)

    mapping = {
        "__UX__": f"{ux:.8f}",
        "__UY__": f"{uy:.8f}",
        "__KINF__": f"{k_inf:.8e}",
        "__OMEGAINF__": f"{omega_inf:.8e}",
        "__NUTINF__": f"{nut_inf:.8e}",
        "__NU__": f"{NU_AIR:.8e}",
        "__MAGUINF__": f"{u_inf:.6f}",
        "__LIFTX__": f"{lift[0]:.8f}",
        "__LIFTY__": f"{lift[1]:.8f}",
        "__DRAGX__": f"{drag[0]:.8f}",
        "__DRAGY__": f"{drag[1]:.8f}",
        "__AREF__": f"{CHORD * SPAN:.6f}",
        "__ENDTIME__": str(end_time),
    }
    for fname in ("U", "p", "k", "omega", "nut"):
        substitute(os.path.join(case_dir, "0", fname), mapping)
    substitute(os.path.join(case_dir, "constant", "transportProperties"), mapping)
    substitute(os.path.join(case_dir, "system", "controlDict"), mapping)

    os.makedirs(os.path.join(case_dir, "constant", "polyMesh"), exist_ok=True)
    with open(os.path.join(case_dir, "mesh_source.txt"), "w") as f:
        f.write(os.path.abspath(mesh_path) + "\n")
    info = {"airfoil": airfoil, "alpha_deg": alpha, "reynolds": re,
            "u_inf": u_inf, "nu": NU_AIR, "k_inf": k_inf,
            "omega_inf": omega_inf, "nut_inf": nut_inf,
            "chord": CHORD, "span": SPAN, "aref": CHORD * SPAN,
            "mesh": os.path.abspath(mesh_path)}
    with open(os.path.join(case_dir, "case_info.json"), "w") as f:
        json.dump(info, f, indent=2)
    print(f"[case] {case_dir}: NACA {airfoil}, AoA={alpha} deg, "
          f"Re={re:.2e}, U={u_inf:.2f} m/s")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--airfoil", required=True)
    ap.add_argument("--alpha", type=float, required=True)
    ap.add_argument("--re", type=float, default=1e6)
    ap.add_argument("--case", required=True)
    ap.add_argument("--mesh", required=True)
    ap.add_argument("--end-time", type=int, default=4000)
    args = ap.parse_args()
    prepare(args.airfoil, args.alpha, args.re, args.case, args.mesh, args.end_time)


if __name__ == "__main__":
    main()
