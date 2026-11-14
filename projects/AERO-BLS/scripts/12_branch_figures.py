import csv, glob, re
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

rows = list(csv.DictReader(open("data/results/summary.csv")))
idx = {r["case"]: r for r in rows}

# ---- Figure A: separation location vs grid refinement ----
series = {
    r"$\alpha=5^\circ$":  [("coarse", 16100, "naca0012_aoa5_coarse"), ("medium", 33740, "naca0012_aoa5"), ("fine", 51360, "naca0012_aoa5_fine")],
    r"$\alpha=15^\circ$": [("coarse", 16100, "naca0012_aoa15_coarse"), ("medium", 33740, "naca0012_aoa15"), ("fine", 51360, "naca0012_aoa15_fine")],
    r"$\alpha=18^\circ$": [("medium", 33740, "naca0012_aoa18"), ("fine", 51360, "naca0012_aoa18_fine")],
}
fig, ax = plt.subplots(figsize=(6.4, 3.6))
colors = {"$\\alpha=5^\\circ$": "#2c7fb8", "$\\alpha=15^\\circ$": "#d95f02", "$\\alpha=18^\\circ$": "#7570b3"}
for label, pts in series.items():
    xs, ys, txts = [], [], []
    for name, n, case in pts:
        r = idx[case]
        sep = r["sep_upper"].strip()
        xs.append(n)
        if sep == "":
            ys.append(1.06); txts.append("attached")
        else:
            ys.append(float(sep)); txts.append(f"{float(sep):.2f}")
    ax.plot(xs, ys, "o-", color=colors[label], label=label, clip_on=False)
    for j, (x, y, t) in enumerate(zip(xs, ys, txts)):
        dy = 9 if j % 2 == 0 else -14
        ax.annotate(t, (x, y), textcoords="offset points", xytext=(12, dy), ha="left", fontsize=9, color=colors[label],
                    bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85), zorder=5)
ax.set_xscale("log")
ax.set_xticks([16100, 33740, 51360]); ax.set_xticklabels(["16,100\n(coarse)", "33,740\n(medium)", "51,360\n(fine)"], fontsize=8)
ax.set_xlabel("Grid cells")
ax.set_ylabel("Upper-surface separation $x_{sep}/c$")
ax.set_ylim(-0.05, 1.2)
ax.axhline(1.0, color="grey", lw=0.6, ls=":")
ax.text(16500, 1.09, "trailing edge (attached)", fontsize=7, color="grey")
ax.legend(fontsize=8, loc="center left")
ax.set_title("Grid refinement selects different separation branches at stall", fontsize=9)
fig.tight_layout()
fig.savefig("latex/figures/figA_branch_map.png", dpi=200)
print("figA done")

# ---- Figure B: sliding-window means (non-stationarity) ----
def load_urans(case):
    ts, cls, cds = [], [], []
    for f in sorted(glob.glob(f"openfoam_cases/{case}/postProcessing/forceCoeffs/*/coefficient*.dat")):
        for line in open(f):
            if line.startswith("#") or not line.strip():
                continue
            p = line.split()
            try:
                t = float(p[0]); cd = float(p[1]); cl = float(p[4])
            except (ValueError, IndexError):
                continue
            if t >= 0.05:
                ts.append(t); cls.append(cl); cds.append(cd)
    order = sorted(range(len(ts)), key=lambda i: ts[i])
    return [ts[i] for i in order], [cls[i] for i in order], [cds[i] for i in order]

U_INF = 5.4
fig, axes = plt.subplots(1, 2, figsize=(7.6, 3.0))
steady = {"12": (1.039, 0.040), "15": (1.058, 0.075)}
exp = {"12": (0.483, 0.022), "15": (0.331, 0.190)}
for ang, case in [("12", "naca0012tr_aoa12_re360k_urans"), ("15", "naca0012tr_aoa15_re360k_urans")]:
    ts, cls, cds = load_urans(case)
    tstar = [t * U_INF for t in ts]
    W = 0.10 * U_INF  # 0.1 s window in convective units
    centers, mcl, mcd = [], [], []
    t0, t1 = tstar[0], tstar[-1]
    c = t0 + W / 2
    while c + W / 2 <= t1:
        w = [i for i in range(len(tstar)) if c - W/2 <= tstar[i] <= c + W/2]
        if len(w) >= 5:
            centers.append(c)
            mcl.append(sum(cls[i] for i in w) / len(w))
            mcd.append(sum(cds[i] for i in w) / len(w))
        c += W / 4
    color = "#d95f02" if ang == "15" else "#2c7fb8"
    axes[0].plot(centers, mcl, "o-", ms=3, color=color, label=fr"$\alpha={ang}^\circ$")
    axes[1].plot(centers, mcd, "o-", ms=3, color=color)
    axes[0].axhline(steady[ang][0], color=color, ls="--", lw=0.8, alpha=0.6)
    axes[0].axhline(exp[ang][0], color=color, ls=":", lw=0.8, alpha=0.6)
    axes[1].axhline(steady[ang][1], color=color, ls="--", lw=0.8, alpha=0.6)
    axes[1].axhline(exp[ang][1], color=color, ls=":", lw=0.8, alpha=0.6)
axes[0].set_ylabel("Windowed mean $C_L$"); axes[1].set_ylabel("Windowed mean $C_D$")
for a in axes:
    a.set_xlabel(r"Window centre $t^*$ ($c/U_\infty$)")
axes[0].set_title("lift: drift away from steady (dashed)", fontsize=9)
axes[1].set_title("drag: steady (dashed) and experiment (dotted) bracketed", fontsize=9)
axes[0].legend(fontsize=8)
fig.tight_layout()
fig.savefig("latex/figures/figB_window_drift.png", dpi=200)
print("figB done")
