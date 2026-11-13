#!/usr/bin/env python3
"""Post-process all OpenFOAM cases: forces, Cp, Cf, separation, convergence.

Writes:
  data/results/summary.csv          one row per case (Cl, Cd, Cm, Cdp, Cdf, convergence)
  data/results/cp/<case>.csv        Cp(x/c) on upper/lower surface
  data/results/cf/<case>.csv        Cf(x/c) + separation location per surface
  data/results/convergence/<case>.csv  residual + force history
"""
import glob
import json
import os
import re

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RESULTS = os.path.join(ROOT, "data", "results")


def read_coefficients(case):
    path = os.path.join(case, "postProcessing", "forceCoeffs", "0", "coefficient.dat")
    if not os.path.exists(path):
        hits = glob.glob(os.path.join(case, "postProcessing", "forceCoeffs", "*", "coefficient.dat"))
        if not hits:
            return None
        path = hits[-1]
    t, cd, cl, cm = [], [], [], []
    for line in open(path):
        if line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) >= 8:
            t.append(float(parts[0]))
            cd.append(float(parts[1]))
            cl.append(float(parts[4]))
            cm.append(float(parts[7]))
    arr = np.array([t, cd, cl, cm])
    return arr if arr.size else None


def pressure_viscous_split(log_path):
    """Parse the last forceCoeffs console block for pressure/viscous split."""
    if not os.path.exists(log_path):
        return np.nan, np.nan
    cdp, cdf = np.nan, np.nan
    for line in open(log_path):
        m = re.match(r"\s+Cd:\s+(\S+)\s+(\S+)\s+(\S+)", line)
        if m:
            cdp, cdf = float(m.group(2)), float(m.group(3))
    return cdp, cdf


def read_raw_surface(path):
    rows = []
    for line in open(path):
        if line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) >= 4:
            rows.append([float(v) for v in parts])
    return np.array(rows)


def surface_quantities(case, u_inf, alpha_deg=0.0):
    """Return cp_table, cf_table, separation dict.

    Surface files are POINT_DATA polylines ordered TE -> upper -> LE ->
    lower -> TE with consecutive duplicates; Cf uses the local wall tangent
    oriented downstream (positive dot with the freestream direction).
    """
    hits = sorted(glob.glob(os.path.join(case, "postProcessing", "surfaces", "*")),
                  key=lambda p: float(os.path.basename(p)) if os.path.basename(p).isdigit() else 0)
    if not hits:
        return None, None, None
    latest = hits[-1]
    pfile = os.path.join(latest, "p_airfoil.raw")
    wfile = os.path.join(latest, "wallShearStress_airfoil.raw")
    if not os.path.exists(pfile):
        return None, None, None
    d = read_raw_surface(pfile)
    keep = np.ones(len(d), dtype=bool)
    keep[1:] = np.hypot(np.diff(d[:, 0]), np.diff(d[:, 1])) > 1e-12
    d = d[keep]
    x, y, p = d[:, 0], d[:, 1], d[:, 3]
    q = 0.5 * u_inf ** 2
    cp = p / q
    i_le = int(np.argmin(x))
    if len(d) < 10 or i_le in (0, len(d) - 1):
        return None, None, None
    upper_first = y[1] >= y[-2]
    segs = [("upper", 0, i_le + 1), ("lower", i_le, len(d))]
    if not upper_first:
        segs = [("lower", 0, i_le + 1), ("upper", i_le, len(d))]
    cp_rows = []
    for name, i0, i1 in segs:
        for xi, yi, ci in zip(x[i0:i1], y[i0:i1], cp[i0:i1]):
            cp_rows.append((name, xi, yi, ci))
    cf_rows, sep = [], {"upper": None, "lower": None}
    if not os.path.exists(wfile):
        return cp_rows, cf_rows, sep
    w = read_raw_surface(wfile)[keep]
    taux, tauy = w[:, 3], w[:, 4]
    u_dir = np.array([np.cos(np.radians(alpha_deg)), np.sin(np.radians(alpha_deg))])
    for name, i0, i1 in segs:
        xs, ys = x[i0:i1], y[i0:i1]
        tx = np.gradient(xs)
        ty = np.gradient(ys)
        ln = np.hypot(tx, ty)
        good = ln > 1e-12
        tx, ty, ln = tx[good], ty[good], ln[good]
        sx, sy, wx2, wy2 = xs[good], ys[good], taux[i0:i1][good], tauy[i0:i1][good]
        flip = (tx * u_dir[0] + ty * u_dir[1]) < 0
        tx[flip], ty[flip] = -tx[flip], -ty[flip]
        cfv = -(wx2 * tx + wy2 * ty) / (ln * q)
        for xi, ci in zip(sx, cfv):
            cf_rows.append((name, xi, ci))
        o = np.argsort(sx)
        sx_s, cfv_s = sx[o], cfv[o]
        body = sx_s < 0.995
        sx_s, cfv_s = sx_s[body], cfv_s[body]
        rev = cfv_s < -1e-6
        if rev.any():
            run = np.convolve(rev.astype(int), [1, 1], mode="valid") == 2
            if run.any():
                sep[name] = float(sx_s[np.where(run)[0][0]])
    return cp_rows, cf_rows, sep


def final_residuals(case):
    hits = glob.glob(os.path.join(case, "postProcessing", "residuals", "*", "residuals.dat"))
    if not hits:
        return {}
    lines = [l for l in open(hits[-1]) if not l.startswith("#")]
    if not lines:
        return {}
    parts = lines[-1].split()
    names = [l.strip() for l in open(hits[-1]) if l.startswith("#") and "Initial residual" in l]
    # residual file header: # Time  p  U  k  omega (initial residuals)
    header = None
    for l in open(hits[-1]):
        if l.startswith("#") and "Time" in l:
            header = l.strip("# \n").split()
    out = {}
    if header and len(parts) == len(header):
        for name, val in zip(header[1:], parts[1:]):
            out[name] = float(val)
    return out


def main():
    os.makedirs(os.path.join(RESULTS, "cp"), exist_ok=True)
    os.makedirs(os.path.join(RESULTS, "cf"), exist_ok=True)
    os.makedirs(os.path.join(RESULTS, "convergence"), exist_ok=True)
    summary = []
    for case in sorted(glob.glob(os.path.join(ROOT, "openfoam_cases", "*"))):
        if not os.path.isdir(case) or os.path.basename(case) == "base_case":
            continue
        name = os.path.basename(case)
        info_path = os.path.join(case, "case_info.json")
        if not os.path.exists(info_path):
            continue
        info = json.load(open(info_path))
        coeffs = read_coefficients(case)
        if coeffs is None:
            print(f"[skip] {name}: no forceCoeffs output")
            continue
        t, cd, cl, cm = coeffs
        n_tail = max(1, int(0.2 * len(t)))
        cl_m, cl_s = cl[-n_tail:].mean(), cl[-n_tail:].std()
        cd_m, cd_s = cd[-n_tail:].mean(), cd[-n_tail:].std()
        cm_m = cm[-n_tail:].mean()
        cdp, cdf = pressure_viscous_split(os.path.join(ROOT, "logs", name + ".log"))
        res = final_residuals(case)
        cp_rows, cf_rows, sep = surface_quantities(case, info["u_inf"], info["alpha_deg"])
        if cp_rows:
            with open(os.path.join(RESULTS, "cp", name + ".csv"), "w") as f:
                f.write("surface,x_over_c,y_over_c,Cp\n")
                for r in cp_rows:
                    f.write("%s,%.6f,%.6f,%.6f\n" % r)
        if cf_rows:
            with open(os.path.join(RESULTS, "cf", name + ".csv"), "w") as f:
                f.write("surface,x_over_c,Cf\n")
                for r in cf_rows:
                    f.write("%s,%.6f,%.6f\n" % r)
        with open(os.path.join(RESULTS, "convergence", name + ".csv"), "w") as f:
            f.write("time,Cd,Cl,Cm\n")
            for ti, ci, li, mi in zip(t, cd, cl, cm):
                f.write("%.1f,%.6f,%.6f,%.6f\n" % (ti, ci, li, mi))
        summary.append({
            "case": name, "airfoil": info["airfoil"], "alpha_deg": info["alpha_deg"],
            "reynolds": info["reynolds"], "Cl": cl_m, "Cl_std": cl_s,
            "Cd": cd_m, "Cd_std": cd_s, "Cm": cm_m,
            "Cd_pressure": cdp, "Cd_viscous": cdf,
            "Cl_Cd": cl_m / cd_m if cd_m else np.nan,
            "sep_upper": sep["upper"] if sep else None,
            "sep_lower": sep["lower"] if sep else None,
            "iters": len(t), "residuals": json.dumps(res),
        })
        print(f"[ok] {name}: Cl={cl_m:.4f} Cd={cd_m:.5f} "
              f"sep_up={sep['upper'] if sep else None}")
    keys = sorted({k for s in summary for k in s})
    with open(os.path.join(RESULTS, "summary.csv"), "w") as f:
        f.write(",".join(keys) + "\n")
        for s in summary:
            f.write(",".join("" if s.get(k) is None else str(s.get(k)) for k in keys) + "\n")
    print(f"[done] {len(summary)} cases -> data/results/summary.csv")


if __name__ == "__main__":
    main()
