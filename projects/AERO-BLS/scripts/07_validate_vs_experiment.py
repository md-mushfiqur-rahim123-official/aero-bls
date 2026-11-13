#!/usr/bin/env python3
"""Validate CFD results against NASA/NACA experimental data.

Primary dataset: NACA 0012 from the NASA Langley Turbulence Modeling
Resource (TMR, https://tmbwg.github.io/turbmodels/naca0012_val.html):
  - Ladson (NASA TM 4074, 1988): Cl/Cd, Re=6e6, transition tripped
  - Abbott & von Doenhoff (1959): Cl/Cd, Re=6e6, free transition (digitized)
  - Gregory & O'Reilly (R&M 3726, 1970): Cl (Re=2.88e6) and Cp (Re=2.88e6)
  - CFL3D SST reference results (code-to-code verification), Re=6e6

Secondary: key metrics for NACA 2412 / 4412 transcribed from
Abbott & von Doenhoff, "Theory of Wing Sections" (1959), Re=6e6,
clearly labelled as approximate transcriptions.

Writes:
  data/validation/validation_table.csv   per-alpha CFD vs experiment (0012)
  data/validation/cp_comparison_<alpha>.csv
  data/validation/metrics.json           MAE / stall / Cd_min summary
"""
import csv
import glob
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
EXP = os.path.join(ROOT, "data", "experimental")
VALID = os.path.join(ROOT, "data", "validation")

# Key metrics transcribed from Abbott & von Doenhoff (1959), Re=6e6.
ABBOTT_METRICS = {
    "naca0012": {"alpha0_deg": 0.0, "cl_max": 1.50, "alpha_clmax_deg": 16.0},
    "naca2412": {"alpha0_deg": -2.1, "cl_max": 1.60, "alpha_clmax_deg": 16.0},
    "naca4412": {"alpha0_deg": -4.1, "cl_max": 1.67, "alpha_clmax_deg": 16.0},
}


def parse_tecplot(path):
    """Parse TMR .dat files: # comments, variables line, zone lines."""
    zones = []
    current = None
    columns = []
    for line in open(path):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("variables"):
            columns = [c.strip().strip('"') for c in line.split("=", 1)[1].split(",")]
            continue
        if line.lower().startswith("zone"):
            current = {"name": line, "data": []}
            zones.append(current)
            continue
        parts = line.replace(",", " ").split()
        try:
            vals = [float(p) for p in parts]
        except ValueError:
            continue
        if current is None:
            current = {"name": "default", "data": []}
            zones.append(current)
        current["data"].append(vals)
    for z in zones:
        z["data"] = np.array(z["data"])
    return zones, columns


def load_summary():
    path = os.path.join(ROOT, "data", "results", "summary.csv")
    with open(path) as f:
        return list(csv.DictReader(f))


def interp_exp(alpha, exp_alpha, exp_cl):
    order = np.argsort(exp_alpha)
    return np.interp(alpha, exp_alpha[order], exp_cl[order])


def main():
    os.makedirs(VALID, exist_ok=True)
    summary = load_summary()

    # ---------------- NACA 0012 experimental datasets ----------------
    ladson_z, _ = parse_tecplot(os.path.join(EXP, "CLCD_Ladson_expdata.dat"))
    ladson = ladson_z[0]["data"]  # zone "80 grit": tripped, Re=6e6
    abbott_cl_z, _ = parse_tecplot(os.path.join(EXP, "0012.abbottdata.cl.dat"))
    abbott_cd_z, _ = parse_tecplot(os.path.join(EXP, "0012.abbottdata.cd.dat"))
    gregory_cl_z, _ = parse_tecplot(os.path.join(EXP, "CL_Gregory_expdata.dat"))
    cfl3d_z, _ = parse_tecplot(os.path.join(EXP, "n0012clcd_cfl3d_sst.dat"))
    cfl3d = cfl3d_z[0]["data"]

    cfd0012 = sorted(
        [r for r in summary
         if r["airfoil"] == "naca0012" and float(r["reynolds"]) == 1e6
         and r["case"] == f"naca0012_aoa{int(float(r['alpha_deg']))}"],
        key=lambda r: float(r["alpha_deg"]))
    if not cfd0012:
        raise SystemExit("no NACA 0012 Re=1e6 cases in summary.csv")

    rows = []
    for r in cfd0012:
        a = float(r["alpha_deg"])
        cl, cd = float(r["Cl"]), float(r["Cd"])
        row = {
            "airfoil": "naca0012", "alpha_deg": a,
            "Cl_CFD": cl, "Cd_CFD": cd,
            "Cl_Ladson": interp_exp(a, ladson[:, 0], ladson[:, 1]),
            "Cd_Ladson": interp_exp(a, ladson[:, 0], ladson[:, 2]),
            "Cl_Abbott": interp_exp(a, abbott_cl_z[0]["data"][:, 0], abbott_cl_z[0]["data"][:, 1]),
            "Cd_Abbott": interp_exp(a, abbott_cd_z[0]["data"][:, 0], abbott_cd_z[0]["data"][:, 1]),
            "Cl_Gregory": interp_exp(a, gregory_cl_z[0]["data"][:, 0], gregory_cl_z[0]["data"][:, 1]),
        }
        row["Cl_err_Ladson_pct"] = 100 * (row["Cl_CFD"] - row["Cl_Ladson"]) / abs(row["Cl_Ladson"]) if row["Cl_Ladson"] else np.nan
        row["Cd_err_Ladson_pct"] = 100 * (row["Cd_CFD"] - row["Cd_Ladson"]) / row["Cd_Ladson"]
        m = np.isclose(cfl3d[:, 0], a)
        if m.any():
            row["Cl_CFL3D_SST"] = float(cfl3d[m][0, 1])
            row["Cd_CFL3D_SST"] = float(cfl3d[m][0, 2])
        rows.append(row)

    cols = ["airfoil", "alpha_deg", "Cl_CFD", "Cl_Ladson", "Cl_err_Ladson_pct",
            "Cl_Abbott", "Cl_Gregory", "Cd_CFD", "Cd_Ladson", "Cd_err_Ladson_pct",
            "Cd_Abbott", "Cl_CFL3D_SST", "Cd_CFL3D_SST"]
    with open(os.path.join(VALID, "validation_table.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for row in rows:
            w.writerow({k: (f"{row[k]:.6g}" if isinstance(row.get(k), float) else row.get(k, ""))
                        for k in cols})

    # ---------------- summary metrics ----------------
    pre_stall = [r for r in rows if r["alpha_deg"] <= 14.22]
    attached = [r for r in rows if r["alpha_deg"] <= 10.0]
    cl_mae = np.mean([abs(r["Cl_CFD"] - r["Cl_Ladson"]) for r in pre_stall])
    sig = [r for r in pre_stall if abs(r["Cl_Ladson"]) >= 0.3]
    cl_mae_pct = np.mean([abs(r["Cl_err_Ladson_pct"]) for r in sig]) if sig else np.nan
    cd_mae_pct = np.mean([abs(r["Cd_err_Ladson_pct"]) for r in pre_stall])
    cl_mae_abbott = np.mean([abs(r["Cl_CFD"] - r["Cl_Abbott"]) for r in pre_stall])
    cl_mae_gregory = np.mean([abs(r["Cl_CFD"] - r["Cl_Gregory"]) for r in pre_stall])
    cd_mae_ladson_abs = np.mean([abs(r["Cd_CFD"] - r["Cd_Ladson"]) for r in pre_stall])

    cl_curve = [(float(r["alpha_deg"]), float(r["Cl"])) for r in cfd0012]
    alpha_stall_cfd = max(cl_curve, key=lambda p: p[1])[0]
    cl_max_cfd = max(c for _, c in cl_curve)
    ladson_pre = ladson[ladson[:, 0] <= 17.2]
    i = np.argmax(ladson_pre[:, 1])
    alpha_stall_exp, cl_max_exp = ladson_pre[i, 0], ladson_pre[i, 1]

    cfl3d_pts = []
    for row in rows:
        if "Cl_CFL3D_SST" not in row:
            continue
        cl_abs = row["Cl_CFD"] - row["Cl_CFL3D_SST"]
        cd_pct = 100 * (row["Cd_CFD"] - row["Cd_CFL3D_SST"]) / row["Cd_CFL3D_SST"]
        cl_pct = (100 * cl_abs / row["Cl_CFL3D_SST"]
                  if abs(row["Cl_CFL3D_SST"]) >= 0.3 else np.nan)
        cfl3d_pts.append({
            "alpha": row["alpha_deg"],
            "cl_diff_abs": cl_abs,
            "cl_diff_pct": cl_pct,
            "cd_diff_pct": cd_pct,
        })

    # ---------------- cambered airfoils: Abbott metrics ----------------
    camber_rows = []
    for airfoil, met in ABBOTT_METRICS.items():
        cases = sorted([r for r in summary if r["airfoil"] == airfoil
                        and float(r["reynolds"]) == 1e6 and "serr" not in r["case"]],
                       key=lambda r: float(r["alpha_deg"]))
        if not cases:
            continue
        cl_by_a = {float(r["alpha_deg"]): float(r["Cl"]) for r in cases}
        cl0_cfd = cl_by_a.get(0.0, np.nan)
        a_stall_cfd, clmax_cfd = max(cl_by_a.items(), key=lambda kv: kv[1])
        camber_rows.append({
            "airfoil": airfoil,
            "Cl0_CFD": cl0_cfd, "Cl0_Abbott_transcribed": cl0_cfd - 0.1096 * met["alpha0_deg"],
            "alpha0_Abbott_transcribed": met["alpha0_deg"],
            "Clmax_CFD": clmax_cfd, "Clmax_Abbott_transcribed": met["cl_max"],
            "alpha_stall_CFD": a_stall_cfd, "alpha_stall_Abbott_transcribed": met["alpha_clmax_deg"],
        })
    with open(os.path.join(VALID, "camber_airfoil_metrics.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(camber_rows[0].keys()))
        w.writeheader()
        w.writerows(camber_rows)

    # ---------------- Cp comparison (Gregory, upper surface) ----------------
    cp_mae = {}
    cp_files = {0: None, 10: None, 15: None}
    for a in cp_files:
        hits = glob.glob(os.path.join(ROOT, "data", "results", "cp", f"naca0012_aoa{a}.csv"))
        if hits:
            cp_files[a] = hits[0]
    greg_cp_z, _ = parse_tecplot(os.path.join(EXP, "CP_Gregory_expdata.dat"))
    for zone, a in zip(greg_cp_z, (0, 10, 15)):
        if cp_files[a] is None:
            continue
        cfd = list(csv.DictReader(open(cp_files[a])))
        up = [(float(r["x_over_c"]), float(r["Cp"])) for r in cfd if r["surface"] == "upper"]
        up.sort()
        xs = np.array([u[0] for u in up])
        cs = np.array([u[1] for u in up])
        ex = zone["data"][:, 0]
        ec = zone["data"][:, 1]
        order = np.argsort(xs)
        exp_at_cfd = np.interp(xs[order], ex, ec)
        out = os.path.join(VALID, f"cp_comparison_aoa{a}.csv")
        with open(out, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["x_over_c", "Cp_CFD_upper", "Cp_Gregory_interp_upper", "Cp_diff"])
            for xi, ci, ei in zip(xs[order], cs[order], exp_at_cfd):
                w.writerow([f"{xi:.6f}", f"{ci:.6f}", f"{ei:.6f}", f"{ci - ei:.6f}"])
        mae = np.mean(np.abs(cs[order] - exp_at_cfd))
        cp_mae[a] = float(mae)
        print(f"[cp] aoa{a}: upper-surface Cp MAE vs Gregory = {mae:.4f}")

    metrics = {
        "cp_mae_gregory": cp_mae,
        "naca0012_ladson": {
            "cl_mae_pre_stall": cl_mae,
            "cl_mae_pre_stall_pct": cl_mae_pct,
            "cl_mae_attached": float(np.mean([abs(r["Cl_CFD"] - r["Cl_Ladson"]) for r in attached])),
            "cl_mae_attached_pct": float(np.mean([abs(r["Cl_err_Ladson_pct"]) for r in attached
                                                  if abs(r["Cl_Ladson"]) >= 0.3])),
            "cd_mae_pre_stall_pct": cd_mae_pct,
            "cl_mae_abbott": cl_mae_abbott,
            "cl_mae_gregory": cl_mae_gregory,
            "cd_mae_ladson_abs": cd_mae_ladson_abs,
            "alpha_stall_CFD_deg": alpha_stall_cfd,
            "alpha_stall_Ladson_deg": float(alpha_stall_exp),
            "cl_max_CFD": cl_max_cfd,
            "cl_max_Ladson": float(cl_max_exp),
        },
        "cfl3d_sst_verification": cfl3d_pts,
        "camber_airfoils": camber_rows,
        "notes": [
            "CFD at Re=1e6, fully turbulent k-omega SST; Ladson data Re=6e6, tripped.",
            "Fully turbulent CFD overpredicts drag vs tripped experiment; expected bias.",
            "Gregory Cp data Re=2.88e6; TMR recommends it for surface-pressure validation.",
        ],
    }
    with open(os.path.join(VALID, "metrics.json"), "w") as f:
        json.dump(metrics, f, indent=2)

    print(f"[done] validation table + metrics -> {VALID}")
    print(f"  Cl MAE (pre-stall, vs Ladson): {cl_mae:.4f} ({cl_mae_pct:.1f}%)")
    print(f"  Cd MAE (pre-stall, vs Ladson): {cd_mae_pct:.1f}%")
    print(f"  stall: CFD alpha={alpha_stall_cfd:.1f} deg (Clmax={cl_max_cfd:.3f}), "
          f"Ladson alpha={alpha_stall_exp:.1f} deg (Clmax={cl_max_exp:.3f})")


if __name__ == "__main__":
    main()
