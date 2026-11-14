#!/usr/bin/env python3
"""Post-process the URANS verification cases (pimpleFoam, k-omega SST).

For each *-urans case:
  - reads postProcessing/forceCoeffs/*/coefficient.dat (stage1 + stage2 history)
  - averages the production window (t >= T_PRODUCTION) -> mean +/- std
  - compares against the steady simpleFoam sibling (summary.csv) and the
    Reynolds-matched Sheldahl & Klimas (SAND80-2114) experiment

Writes:
  data/results/urans_summary.csv
  latex/figures/fig15_urans_history.png
"""
import csv
import glob
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RESULTS = os.path.join(ROOT, "data", "results")
EXP = os.path.join(ROOT, "data", "experimental")
FIGS = os.path.join(ROOT, "latex", "figures")

T_PRODUCTION = 0.15
SANDIA = os.path.join(EXP, "sheldahl1981_naca0012_re360k.dat")

plt.rcParams.update({
    "font.size": 11, "axes.grid": True, "grid.alpha": 0.3,
    "figure.constrained_layout.use": True, "savefig.dpi": 300,
})


def read_coefficients(case):
    hits = glob.glob(os.path.join(case, "postProcessing", "forceCoeffs",
                                  "*", "coefficient*.dat"))
    if not hits:
        return None
    rows = []
    for path in hits:
        for line in open(path):
            if line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 8:
                rows.append((float(parts[0]), float(parts[1]),
                             float(parts[4]), float(parts[7])))
    if not rows:
        return None
    rows.sort(key=lambda r: r[0])
    dedup = {}
    for r in rows:
        dedup[round(r[0], 6)] = r
    rows = [dedup[k] for k in sorted(dedup)]
    arr = np.array(rows)
    return arr[:, 0], arr[:, 1], arr[:, 2], arr[:, 3]


def load_steady():
    path = os.path.join(RESULTS, "summary.csv")
    if not os.path.exists(path):
        return {}
    with open(path) as f:
        return {r["case"]: r for r in csv.DictReader(f)}


def load_sandia():
    if not os.path.exists(SANDIA):
        return {}
    out = {}
    for line in open(SANDIA):
        if line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) >= 3:
            out[round(float(parts[0]), 2)] = (float(parts[1]), float(parts[2]))
    return out


def main():
    steady = load_steady()
    sandia = load_sandia()
    rows = []
    cases = sorted(glob.glob(os.path.join(ROOT, "openfoam_cases", "*_urans")))
    cases = [c for c in cases if os.path.isdir(c)]
    for case in cases:
        name = os.path.basename(case)
        info_path = os.path.join(case, "case_info.json")
        if not os.path.exists(info_path):
            continue
        data = read_coefficients(case)
        if data is None:
            print(f"[skip] {name}: no forceCoeffs output")
            continue
        t, cd, cl, cm = data
        t_end = float(t[-1])
        if t_end < T_PRODUCTION:
            print(f"[skip] {name}: history ends at t={t_end:.3f} "
                  f"(< production window start {T_PRODUCTION})")
            continue
        m = t >= T_PRODUCTION
        cl_m, cl_s = float(cl[m].mean()), float(cl[m].std())
        cd_m, cd_s = float(cd[m].mean()), float(cd[m].std())
        cm_m = float(cm[m].mean())

        info = json.load(open(info_path))
        alpha = float(info["alpha_deg"])
        sib = steady.get(name[:-len("_urans")], {})
        st_cl = float(sib["Cl"]) if sib.get("Cl") else np.nan
        st_cd = float(sib["Cd"]) if sib.get("Cd") else np.nan

        exp_cl = exp_cd = np.nan
        for da in (alpha, round(alpha), int(alpha)):
            if da in sandia:
                exp_cl, exp_cd = sandia[da]
                break

        row = {
            "case": name, "airfoil": info["airfoil"], "alpha_deg": alpha,
            "reynolds": info["reynolds"],
            "Cl_mean": cl_m, "Cl_std": cl_s,
            "Cd_mean": cd_m, "Cd_std": cd_s, "Cm_mean": cm_m,
            "t_start": float(t[m][0]), "t_end": t_end, "n_samples": int(m.sum()),
            "steady_Cl": st_cl, "steady_Cd": st_cd,
            "exp_Cl": exp_cl, "exp_Cd": exp_cd,
            "Cl_minus_steady": cl_m - st_cl if not np.isnan(st_cl) else np.nan,
            "Cd_minus_steady": cd_m - st_cd if not np.isnan(st_cd) else np.nan,
            "Cl_err_vs_exp": (cl_m - exp_cl) / abs(exp_cl) if not np.isnan(exp_cl) else np.nan,
            "Cd_err_vs_exp": (cd_m - exp_cd) / abs(exp_cd) if not np.isnan(exp_cd) else np.nan,
        }
        rows.append(row)
        print(f"[ok] {name}: Cl={cl_m:.4f}+/-{cl_s:.4f}  Cd={cd_m:.5f}+/-{cd_s:.5f} "
              f"(t={t[m][0]:.2f}-{t_end:.2f}, n={int(m.sum())})  "
              f"steady Cl={st_cl:.4f} Cd={st_cd:.5f}  "
              f"exp Cl={exp_cl:.3f} Cd={exp_cd:.4f}")

    if not rows:
        print("[done] no completed URANS cases yet")
        return

    keys = list(rows[0].keys())
    with open(os.path.join(RESULTS, "urans_summary.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    print(f"[done] {len(rows)} URANS cases -> data/results/urans_summary.csv")

    os.makedirs(FIGS, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.6))
    colors = {"12": "tab:blue", "15": "tab:red"}
    for row in rows:
        case = os.path.join(ROOT, "openfoam_cases", row["case"])
        t, cd, cl, _ = read_coefficients(case)
        u_inf = float(row["reynolds"]) * 1.5e-5
        t = t * u_inf
        c = colors.get(str(int(row["alpha_deg"])), None)
        a = int(row["alpha_deg"])
        axes[0].plot(t, cl, lw=0.8, color=c, label=rf"$\alpha={a}^\circ$")
        axes[0].axhline(row["Cl_mean"], color=c, ls="--", lw=1.0)
        axes[1].plot(t, cd, lw=0.8, color=c)
        axes[1].axhline(row["Cd_mean"], color=c, ls="--", lw=1.0)
        if not np.isnan(row["exp_Cl"]):
            axes[0].axhline(row["exp_Cl"], color=c, ls=":", lw=1.2,
                            alpha=0.8)
    axes[0].set_xlabel(r"$t^*$ ($c/U_\infty$; one unit = one flow pass)")
    axes[0].set_ylabel(r"$C_L$")
    axes[0].set_title("URANS lift history (dashed: mean, dotted: Sheldahl exp.)")
    axes[1].set_xlabel(r"$t^*$ ($c/U_\infty$; one unit = one flow pass)")
    axes[1].set_ylabel(r"$C_D$")
    axes[1].set_title("URANS drag history (dashed: mean)")
    axes[0].legend(frameon=False, fontsize=9)
    fig.savefig(os.path.join(FIGS, "fig15_urans_history.png"))
    print(f"[fig] {FIGS}/fig15_urans_history.png")


if __name__ == "__main__":
    main()
