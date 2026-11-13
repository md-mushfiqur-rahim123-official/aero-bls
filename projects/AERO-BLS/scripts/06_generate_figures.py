#!/usr/bin/env python3
"""Generate the 9 publication figures from data/results and data/validation.

Outputs to latex/figures/*.png (300 dpi) and echoes the list for the report.
Field contours (fig6, fig7) require foamToVTK output in data/fields/<case>/;
they are skipped gracefully when absent.
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
VALID = os.path.join(ROOT, "data", "validation")
EXP = os.path.join(ROOT, "data", "experimental")
FIGS = os.path.join(ROOT, "latex", "figures")

plt.rcParams.update({
    "font.size": 11, "axes.grid": True, "grid.alpha": 0.3,
    "figure.constrained_layout.use": True, "savefig.dpi": 300,
})
COLORS = {"naca0012": "tab:blue", "naca2412": "tab:orange", "naca4412": "tab:green"}
AIRFOILS = ["naca0012", "naca2412", "naca4412"]


def load_summary():
    with open(os.path.join(RESULTS, "summary.csv")) as f:
        return list(csv.DictReader(f))


EXCLUDE_SUBSTRINGS = ("serr", "coarse", "fine", "farfield", "long",
                      "_re1e5", "_re3e5", "tub")


def core_cases(summary, airfoil, reynolds=1e6):
    rows = [r for r in summary
            if r["airfoil"] == airfoil and float(r["reynolds"]) == reynolds
            and all(s not in r["case"] for s in EXCLUDE_SUBSTRINGS)]
    return sorted(rows, key=lambda r: float(r["alpha_deg"]))


def variant_cases(summary, airfoil, reynolds=1e6):
    """Cases for one device variant (baseline/serrated/tubercle) at one Re."""
    rows = [r for r in summary
            if r["airfoil"] == airfoil and float(r["reynolds"]) == reynolds
            and all(s not in r["case"] for s in ("coarse", "fine", "farfield", "long"))]
    return sorted(rows, key=lambda r: float(r["alpha_deg"]))


def series(rows, key):
    a = [float(r["alpha_deg"]) for r in rows]
    v = [float(r[key]) for r in rows]
    return np.array(a), np.array(v)


def parse_tecplot(path):
    zones = []
    current = None
    for line in open(path):
        line = line.strip()
        if not line or line.startswith("#") or line.startswith("variables"):
            continue
        if line.lower().startswith("zone"):
            current = {"data": []}
            zones.append(current)
            continue
        try:
            vals = [float(p) for p in line.replace(",", " ").split()]
        except ValueError:
            continue
        if current is None:
            current = {"data": []}
            zones.append(current)
        current["data"].append(vals)
    for z in zones:
        z["data"] = np.array(z["data"])
    return zones


def fig0_mesh():
    """Zoomed mesh view from the medium Gmsh file."""
    msh = os.path.join(ROOT, "mesh", "naca0012_medium.msh")
    if not os.path.exists(msh):
        print("[skip] fig0_mesh: no mesh file")
        return
    with open(msh) as f:
        lines = f.read().split("\n")
    nodes = {}
    quads = []
    mode = None
    for i, line in enumerate(lines):
        s = line.strip()
        if s == "$Nodes":
            mode = "nodes"; n_expected = int(lines[i + 1]); idx = i + 2; count = 0
            continue
        if s == "$Elements":
            mode = "elements"; n_expected = int(lines[i + 1]); idx = i + 2; count = 0
            continue
        if s.startswith("$"):
            mode = None
            continue
        if mode == "nodes" and count < n_expected:
            parts = s.split()
            if len(parts) == 4:
                nodes[int(parts[0])] = (float(parts[1]), float(parts[2]))
                count += 1
        elif mode == "elements" and count < n_expected:
            parts = s.split()
            if len(parts) == 6 and parts[1] == "3":
                quads.append([int(p) for p in parts[2:]])
            elif len(parts) == 9 and parts[1] == "5":
                quads.append([int(p) for p in parts[5:]])
            count += 1
    fig, ax = plt.subplots(figsize=(9, 3.2))
    for q in quads:
        pts = [nodes[n] for n in q if n in nodes]
        if len(pts) < 3:
            continue
        xs = [p[0] for p in pts] + [pts[0][0]]
        ys = [p[1] for p in pts] + [pts[0][1]]
        if min(xs) > 1.3 or max(xs) < -0.15 or min(ys) > 0.35 or max(ys) < -0.35:
            continue
        ax.plot(xs, ys, color="0.55", lw=0.25)
    ax.set_xlim(-0.15, 1.3)
    ax.set_ylim(-0.35, 0.35)
    ax.set_aspect("equal")
    ax.set_xlabel("$x/c$")
    ax.set_ylabel("$y/c$")
    ax.set_title("Polar O-grid near the NACA 0012 (medium density, 33{,}740 cells)")
    fig.savefig(os.path.join(FIGS, "fig0_mesh.png"))
    plt.close(fig)


def fig1_polar(summary):
    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    for af in AIRFOILS:
        rows = core_cases(summary, af)
        if not rows:
            continue
        cd, cl = series(rows, "Cd")[1], series(rows, "Cl")[1]
        ax.plot(cd, cl, "o-", color=COLORS[af], label=af.upper())
    ax.set_xlabel(r"$C_D$")
    ax.set_ylabel(r"$C_L$")
    ax.set_title(r"Lift-drag polar, $Re=10^6$ (2D RANS, k-$\omega$ SST)")
    ax.legend()
    fig.savefig(os.path.join(FIGS, "fig1_cl_cd_polar.png"))
    plt.close(fig)


def fig2_cl_alpha(summary):
    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    for af in AIRFOILS:
        rows = core_cases(summary, af)
        if not rows:
            continue
        a, cl = series(rows, "Cl")
        ax.plot(a, cl, "o-", color=COLORS[af], label=af.upper())
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xlabel(r"$\alpha$ [deg]")
    ax.set_ylabel(r"$C_L$")
    ax.set_title(r"Lift curve, $Re=10^6$")
    ax.legend()
    fig.savefig(os.path.join(FIGS, "fig2_cl_vs_alpha.png"))
    plt.close(fig)


def fig3_cd_alpha(summary):
    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    for af in AIRFOILS:
        rows = core_cases(summary, af)
        if not rows:
            continue
        a, cd = series(rows, "Cd")
        ax.plot(a, cd, "o-", color=COLORS[af], label=af.upper())
    ax.set_yscale("log")
    ax.set_xlabel(r"$\alpha$ [deg]")
    ax.set_ylabel(r"$C_D$")
    ax.set_title(r"Drag rise, $Re=10^6$")
    ax.legend()
    fig.savefig(os.path.join(FIGS, "fig3_cd_vs_alpha.png"))
    plt.close(fig)


def fig4_efficiency(summary):
    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    for af in AIRFOILS:
        rows = core_cases(summary, af)
        if not rows:
            continue
        a, l_over_d = series(rows, "Cl_Cd")
        ax.plot(a, l_over_d, "o-", color=COLORS[af], label=af.upper())
    ax.set_xlabel(r"$\alpha$ [deg]")
    ax.set_ylabel(r"$C_L/C_D$")
    ax.set_title(r"Aerodynamic efficiency, $Re=10^6$")
    ax.legend()
    fig.savefig(os.path.join(FIGS, "fig4_cl_cd_ratio.png"))
    plt.close(fig)


def fig5_cp(summary):
    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    alphas = [0, 5, 10, 15, 18]
    cmap = plt.cm.viridis(np.linspace(0, 0.9, len(alphas)))
    for c, a in zip(cmap, alphas):
        path = os.path.join(RESULTS, "cp", f"naca0012_aoa{a}.csv")
        if not os.path.exists(path):
            continue
        rows = list(csv.DictReader(open(path)))
        for surf, ls in (("upper", "-"), ("lower", "--")):
            pts = [(float(r["x_over_c"]), float(r["Cp"])) for r in rows if r["surface"] == surf]
            pts.sort()
            if pts:
                ax.plot([p[0] for p in pts], [p[1] for p in pts], ls,
                        color=c, lw=1.4, label=f"$\\alpha$={a} deg" if surf == "upper" else None)
    ax.invert_yaxis()
    ax.set_xlabel("$x/c$")
    ax.set_ylabel(r"$C_p$")
    ax.set_title(r"NACA 0012 surface pressure, $Re=10^6$")
    ax.legend(ncol=2, fontsize=9)
    fig.savefig(os.path.join(FIGS, "fig5_cp_distribution.png"))
    plt.close(fig)


def read_vtk_fields(case):
    """Read latest foamToVTK output for a case; return dict with x, y and cell fields."""
    try:
        import meshio
    except ImportError:
        return None
    hits = sorted(glob.glob(os.path.join(ROOT, "data", "fields", case, "**", "*.vtu"),
                            recursive=True)
                  + glob.glob(os.path.join(ROOT, "data", "fields", case, "**", "*.vtk"),
                              recursive=True))
    if not hits:
        return None
    mesh = meshio.read(hits[-1])
    pts = mesh.points
    ctr = pts[mesh.cells_dict.get("hexahedron", np.zeros((0, 8)))].mean(axis=1) \
        if "hexahedron" in mesh.cells_dict else None
    if ctr is None:
        return None
    out = {"x": ctr[:, 0], "y": ctr[:, 1]}
    for name, data in mesh.cell_data.items():
        for arr in data:
            if arr.shape[0] == ctr.shape[0]:
                out.setdefault(name, arr)
    if "UMag" not in out and "U" in out:
        U = np.asarray(out["U"])
        out["UMag"] = np.linalg.norm(U[:, :min(3, U.shape[1])], axis=1) if U.ndim > 1 else U
    return out


def contour_figure(field_name, fname, title, clabel, case, log_scale=False):
    fields = read_vtk_fields(case)
    if fields is None or field_name not in fields:
        print(f"[skip] {fname}: no field data for {case} (run foamToVTK)")
        return False
    from scipy.interpolate import griddata
    x, y = fields["x"], fields["y"]
    v = fields[field_name]
    if log_scale:
        v = np.log10(np.maximum(v, 1e-12))
    xi = np.linspace(-0.2, 1.6, 400)
    yi = np.linspace(-0.6, 0.6, 300)
    X, Y = np.meshgrid(xi, yi)
    V = griddata((x, y), v, (X, Y), method="linear")
    fig, ax = plt.subplots(figsize=(7.6, 3.4))
    cf = ax.contourf(X, Y, V, levels=40, cmap="viridis")
    ax.plot([], [])  # placeholder keeps layout stable
    ax.set_xlim(-0.2, 1.6)
    ax.set_ylim(-0.6, 0.6)
    ax.set_aspect("equal")
    ax.set_xlabel("$x/c$")
    ax.set_ylabel("$y/c$")
    ax.set_title(title)
    cbar = fig.colorbar(cf, ax=ax)
    cbar.set_label(clabel)
    fig.savefig(os.path.join(FIGS, fname))
    plt.close(fig)
    return True


def fig6_velocity(case):
    ok = contour_figure("UMag", "fig6_velocity_contour.png",
                        f"Velocity magnitude, {case.replace('_', ' ')}",
                        r"$|U|$ [m/s]", case)
    if ok:
        print(f"[ok] fig6 from {case}")


def fig7_tke(case):
    ok = contour_figure("k", "fig7_turbulent_kinetic_energy.png",
                        f"Turbulent kinetic energy, {case.replace('_', ' ')}",
                        r"$\log_{10} k$ [m$^2$/s$^2$]", case, log_scale=True)
    if ok:
        print(f"[ok] fig7 from {case}")


def fig8_cf(summary):
    fig, ax = plt.subplots(figsize=(6.8, 4.4))
    alphas = [0, 5, 10, 15, 18]
    cmap = plt.cm.viridis(np.linspace(0, 0.9, len(alphas)))
    for c, a in zip(cmap, alphas):
        path = os.path.join(RESULTS, "cf", f"naca0012_aoa{a}.csv")
        if not os.path.exists(path):
            continue
        rows = list(csv.DictReader(open(path)))
        for surf, ls in (("upper", "-"), ("lower", "--")):
            pts = [(float(r["x_over_c"]), float(r["Cf"])) for r in rows if r["surface"] == surf]
            pts.sort()
            if pts:
                ax.plot([p[0] for p in pts], [p[1] for p in pts], ls, color=c, lw=1.4,
                        label=f"$\\alpha$={a} deg" if surf == "upper" else None)
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xlabel("$x/c$")
    ax.set_ylabel(r"$C_f$")
    ax.set_title(r"NACA 0012 skin friction (solid: upper, dashed: lower)")
    ax.legend(ncol=2, fontsize=9)
    fig.savefig(os.path.join(FIGS, "fig8_wall_shear_stress.png"))
    plt.close(fig)


def fig10_cp_validation():
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.6), sharey=True)
    greg = parse_tecplot(os.path.join(EXP, "CP_Gregory_expdata.dat"))
    for ax, (zone, a) in zip(axes, zip(greg, (0, 10, 15))):
        path = os.path.join(RESULTS, "cp", f"naca0012_aoa{a}.csv")
        if os.path.exists(path):
            rows = list(csv.DictReader(open(path)))
            pts = sorted((float(r["x_over_c"]), float(r["Cp"]))
                         for r in rows if r["surface"] == "upper")
            ax.plot([p[0] for p in pts], [p[1] for p in pts], "-", color="tab:red",
                    lw=1.6, label="CFD (Re=1e6)")
        ex = zone["data"]
        ax.plot(ex[:, 0], ex[:, 1], "o", color="k", ms=3.5, mfc="none",
                label="Gregory exp. (Re=2.88e6)")
        ax.invert_yaxis()
        ax.set_title(f"$\\alpha = {a}^\\circ$", fontsize=11)
        ax.set_xlabel("$x/c$")
    axes[0].set_ylabel(r"$C_p$ (upper surface)")
    axes[0].legend(fontsize=8)
    fig.suptitle("NACA 0012 surface pressure: CFD vs Gregory & O'Reilly", fontsize=12)
    fig.savefig(os.path.join(FIGS, "fig10_cp_validation.png"))
    plt.close(fig)


def fig9_validation():
    fig, ax = plt.subplots(figsize=(5.8, 4.6))
    ladson = parse_tecplot(os.path.join(EXP, "CLCD_Ladson_expdata.dat"))[0]["data"]
    greg = parse_tecplot(os.path.join(EXP, "CL_Gregory_expdata.dat"))[0]["data"]
    abbott = parse_tecplot(os.path.join(EXP, "0012.abbottdata.cl.dat"))[0]["data"]
    ax.plot(ladson[:, 0], ladson[:, 1], "k-", lw=1.5,
            label="Ladson exp., Re=6e6, tripped")
    ax.plot(greg[:, 0], greg[:, 1], "k--", lw=1.2, label="Gregory exp., Re=2.88e6")
    ax.plot(abbott[:, 0], abbott[:, 1], "k:", lw=1.2, label="Abbott exp., Re=6e6")
    cfl3d = parse_tecplot(os.path.join(EXP, "n0012clcd_cfl3d_sst.dat"))[0]["data"]
    ax.plot(cfl3d[:, 0], cfl3d[:, 1], "s", color="grey", ms=7,
            label="CFL3D SST (NASA TMR)")
    path = os.path.join(RESULTS, "summary.csv")
    rows = core_cases(load_summary(), "naca0012")
    a, cl = series(rows, "Cl")
    a5, cl5 = series(rows, "Cl")
    err = np.array([float(r["Cl_std"]) for r in rows])
    ax.errorbar(a, cl, yerr=err, fmt="o-", color="tab:red", capsize=3,
                label="This work: OpenFOAM SST, Re=1e6")
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xlabel(r"$\alpha$ [deg]")
    ax.set_ylabel(r"$C_L$")
    ax.set_title("NACA 0012: CFD vs experiment")
    ax.legend(fontsize=8)
    fig.savefig(os.path.join(FIGS, "fig9_validation.png"))
    plt.close(fig)


def fig11_re_sweep(summary):
    """Lift and drag curves at three Reynolds numbers, one column per airfoil."""
    fig, axes = plt.subplots(2, 3, figsize=(10.8, 6.2), sharex=True)
    re_list = [(1e5, r"10^5"), (3e5, r"3\times10^5"), (1e6, r"10^6")]
    for col, af in enumerate(AIRFOILS):
        for row, (key, ylab) in enumerate((("Cl", r"$C_L$"), ("Cd", r"$C_D$"))):
            ax = axes[row, col]
            for re, lab in re_list:
                rows = [r for r in summary
                        if r["airfoil"] == af and float(r["reynolds"]) == re
                        and all(s not in r["case"] for s in
                                ("serr", "coarse", "fine", "farfield", "long", "tub"))]
                if not rows:
                    continue
                a, v = series(sorted(rows, key=lambda r: float(r["alpha_deg"])), key)
                ax.plot(a, v, "o-", ms=3.5, label=rf"$Re={lab}$")
            ax.set_title(af.upper(), fontsize=10)
            if row == 1:
                ax.set_xlabel(r"$\alpha$ [deg]")
            if col == 0:
                ax.set_ylabel(ylab)
            if key == "Cd":
                ax.set_yscale("log")
            if row == 0 and col == 2:
                ax.legend(fontsize=8)
    fig.suptitle(r"Reynolds-number sensitivity, 2D RANS k-$\omega$ SST")
    fig.savefig(os.path.join(FIGS, "fig11_re_sweep.png"))
    plt.close(fig)


def fig12_tubercle(summary):
    """Passive-device comparison: baseline vs serrated TE vs tubercle LE surrogate."""
    devices = [("naca0012", "Baseline"),
               ("naca0012s", "Serrated TE"),
               ("naca0012t", "Tubercle LE (2D surrogate)")]
    fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.2))
    for af, lab in devices:
        rows = variant_cases(summary, af)
        if not rows:
            continue
        a, cl = series(rows, "Cl")
        _, cd = series(rows, "Cd")
        axes[0].plot(a, cl, "o-", ms=3.5, label=lab)
        axes[1].plot(a, cd, "o-", ms=3.5, label=lab)
    axes[0].axhline(0, color="k", lw=0.5)
    axes[0].set_xlabel(r"$\alpha$ [deg]")
    axes[0].set_ylabel(r"$C_L$")
    axes[0].set_title(r"Lift curve, $Re=10^6$")
    axes[0].legend(fontsize=8)
    axes[1].set_yscale("log")
    axes[1].set_xlabel(r"$\alpha$ [deg]")
    axes[1].set_ylabel(r"$C_D$")
    axes[1].set_title(r"Drag polar rise, $Re=10^6$")
    axes[1].legend(fontsize=8)
    fig.savefig(os.path.join(FIGS, "fig12_tubercle.png"))
    plt.close(fig)


def main():
    os.makedirs(FIGS, exist_ok=True)
    summary = load_summary()
    fig0_mesh()
    fig1_polar(summary)
    fig2_cl_alpha(summary)
    fig3_cd_alpha(summary)
    fig4_efficiency(summary)
    fig5_cp(summary)
    fig6_velocity("naca0012_aoa15")
    fig7_tke("naca0012_aoa15")
    fig8_cf(summary)
    fig9_validation()
    fig10_cp_validation()
    fig11_re_sweep(summary)
    fig12_tubercle(summary)
    made = sorted(os.path.basename(p) for p in glob.glob(os.path.join(FIGS, "*.png")))
    print("[done] figures:", ", ".join(made))


if __name__ == "__main__":
    main()
