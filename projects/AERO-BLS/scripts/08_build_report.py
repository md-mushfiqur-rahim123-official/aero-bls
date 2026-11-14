#!/usr/bin/env python3
"""Build the LaTeX report: generate numbers.tex macros from results,
render the appendix convergence figure, and compile main.tex with tectonic.
"""
import csv
import glob
import json
import os
import re
import subprocess
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
LATEX = os.path.join(ROOT, "latex")
RESULTS = os.path.join(ROOT, "data", "results")
VALID = os.path.join(ROOT, "data", "validation")
TECTONIC = os.path.join(os.environ.get("MAMBA_ROOT_PREFIX", "/workspace/tools/mamba/root"),
                        "envs", "main", "bin", "tectonic")
if not os.path.exists(TECTONIC):
    TECTONIC = "/workspace/tools/env/bin/tectonic"


def load_summary():
    with open(os.path.join(RESULTS, "summary.csv")) as f:
        return list(csv.DictReader(f))


def get(rows, airfoil, alpha, suffix="", reynolds=1e6, key="Cl"):
    for r in rows:
        if (r["airfoil"] == airfoil and abs(float(r["alpha_deg"]) - alpha) < 0.01
                and float(r["reynolds"]) == reynolds and r["case"].endswith(suffix)
                and ("serr" not in r["case"] if suffix == "" else True)):
            return float(r[key])
    return None


def fmt(v, nd=4):
    return f"{v:.{nd}f}" if v is not None else "n/a"


def sep_text(case_rows, airfoil, alpha):
    for r in case_rows:
        if r["airfoil"] == airfoil and abs(float(r["alpha_deg"]) - alpha) < 0.01:
            v = r.get("sep_upper", "")
            if v and v not in ("None", ""):
                return f"{float(v):.2f}"
            return r"$>0.99$ (attached)"
    return "n/a"


def main():
    summary = load_summary()
    metrics = {}
    if os.path.exists(os.path.join(VALID, "metrics.json")):
        metrics = json.load(open(os.path.join(VALID, "metrics.json")))

    m = {}
    lad = metrics.get("naca0012_ladson", {})
    m["ClMAEpct"] = f"{lad.get('cl_mae_attached_pct', lad.get('cl_mae_pre_stall_pct', 0)):.1f}" if lad else "n/a"
    m["ClMaeLadson"] = f"{lad.get('cl_mae_attached', lad.get('cl_mae_pre_stall', 0)):.4f}" if lad else "n/a"
    m["ClMAEPreStallPct"] = f"{lad.get('cl_mae_pre_stall_pct', 0):.1f}" if lad else "n/a"
    m["ClMaeLadsonPreStall"] = f"{lad.get('cl_mae_pre_stall', 0):.4f}" if lad else "n/a"
    m["ClMaeAbbott"] = f"{lad.get('cl_mae_abbott', 0):.4f}" if lad else "n/a"
    m["ClMaeGregory"] = f"{lad.get('cl_mae_gregory', 0):.4f}" if lad else "n/a"
    m["CdMaeLadsonAbs"] = f"{lad.get('cd_mae_ladson_abs', 0):.4f}" if lad else "n/a"
    cfl = metrics.get("cfl3d_sst_verification", [])
    if cfl:
        cl_pct = [abs(p["cl_diff_pct"]) for p in cfl
                  if p.get("cl_diff_pct") is not None and not np.isnan(p["cl_diff_pct"])
                  and p["alpha"] <= 10]
        cl_abs = [abs(p["cl_diff_abs"]) for p in cfl]
        m["CflThreeDClPct"] = f"{max(cl_pct):.1f}" if cl_pct else f"{max(cl_abs):.3f} abs"
        deep = [p for p in cfl if p["alpha"] == 15]
        if deep:
            d = deep[0]
            m["CflThreeDClFifteen"] = (f"{d['cl_diff_pct']:.0f}" if d.get("cl_diff_pct") is not None
                                       and not np.isnan(d["cl_diff_pct"]) else f"{d['cl_diff_abs']:.2f} abs")
        else:
            m["CflThreeDClFifteen"] = "n/a"
    else:
        m["CflThreeDClPct"] = "n/a"
        m["CflThreeDClFifteen"] = "n/a"

    cp_mae = metrics.get("cp_mae_gregory", {})
    m["CpMaeZero"] = f"{cp_mae.get('0', float('nan')):.3f}" if "0" in cp_mae else "n/a"
    m["CpMaeTen"] = f"{cp_mae.get('10', float('nan')):.3f}" if "10" in cp_mae else "n/a"
    m["CpMaeFifteen"] = f"{cp_mae.get('15', float('nan')):.3f}" if "15" in cp_mae else "n/a"

    def is_core(r):
        name = r["case"]
        return (float(r["reynolds"]) == 1e6 and name == f"{r['airfoil']}_aoa{int(float(r['alpha_deg']))}")
    core = [r for r in summary if is_core(r)]
    core.sort(key=lambda r: float(r["alpha_deg"]))
    cl0 = get(core, "naca0012", 0)
    m["ClZeroZeroOneTwoZero"] = fmt(cl0) if cl0 is not None else "n/a"
    m["SepAten"] = sep_text(core, "naca0012", 10)
    m["SepAfifteen"] = sep_text(core, "naca0012", 15)

    def ldmax(af):
        vals = [float(r["Cl_Cd"]) for r in core if r["airfoil"] == af and r["Cl_Cd"]]
        return fmt(max(vals), 1) if vals else "n/a"
    m["LdMaxZeroZeroOneTwo"] = ldmax("naca0012")
    m["LdMaxFourFourOneTwo"] = ldmax("naca4412")

    def grid_case(suffix):
        for r in summary:
            if r["case"] == f"naca0012_aoa5{suffix}":
                return r
        return None
    coarse, medium, fine = grid_case("_coarse"), grid_case(""), grid_case("_fine")
    sizes = {"_coarse": 16100, "": 33740, "_fine": 51360}
    m["NumCoarse"], m["NumMedium"], m["NumFine"] = "16{,}100", "33{,}740", "51{,}360"
    m["ClCoarse"] = fmt(float(coarse["Cl"])) if coarse else "n/a"
    m["ClMedium"] = fmt(float(medium["Cl"])) if medium else "n/a"
    m["ClFine"] = fmt(float(fine["Cl"])) if fine else "n/a"
    if coarse and medium and fine:
        dcl = 100 * abs(float(fine["Cl"]) - float(medium["Cl"])) / abs(float(medium["Cl"]))
        dcd = 100 * abs(float(fine["Cd"]) - float(medium["Cd"])) / float(medium["Cd"])
        m["ClGridDelta"] = f"{dcl:.2f}\\%"
        m["CdGridDelta"] = f"{dcd:.2f}\\%"
    else:
        m["ClGridDelta"] = m["CdGridDelta"] = "n/a"

    stall_coarse = next((r for r in summary if r["case"] == "naca0012_aoa15_coarse"), None)
    stall_medium = next((r for r in summary if r["case"] == "naca0012_aoa15"), None)
    stall_fine = next((r for r in summary if r["case"] == "naca0012_aoa15_fine"), None)
    if stall_coarse:
        m["ClStallCoarse"] = fmt(float(stall_coarse["Cl"]))
        m["CdStallCoarse"] = fmt(float(stall_coarse["Cd"]))
        v = stall_coarse.get("sep_upper", "")
        m["SepStallCoarse"] = f"{float(v):.2f}" if v and v not in ("None", "") else "n/a"
    else:
        m["ClStallCoarse"] = m["CdStallCoarse"] = m["SepStallCoarse"] = "n/a"
    if stall_coarse and stall_medium and stall_fine:
        eps21 = float(stall_medium["Cl"]) - float(stall_fine["Cl"])
        eps32 = float(stall_coarse["Cl"]) - float(stall_medium["Cl"])
        m["StallClRatio"] = f"{eps32 / eps21:+.2f}" if eps21 != 0 else "n/a"
        eps21d = float(stall_medium["Cd"]) - float(stall_fine["Cd"])
        eps32d = float(stall_coarse["Cd"]) - float(stall_medium["Cd"])
        m["StallCdRatio"] = f"{eps32d / eps21d:+.2f}" if eps21d != 0 else "n/a"
    else:
        m["StallClRatio"] = m["StallCdRatio"] = "n/a"

    sep_rows = []
    for r in core:
        if r["airfoil"] != "naca0012":
            continue
        sep_rows.append(f"{float(r['alpha_deg']):.0f} & {fmt(float(r['Cl']))} & "
                        f"{fmt(float(r['Cd']))} & {fmt(float(r['Cl_Cd']), 1)} & "
                        f"{sep_text(core, 'naca0012', float(r['alpha_deg']))} \\\\")
    m["SepTableRows"] = "\n    ".join(sep_rows)

    re_rows = []
    re_cl, re_cd = {}, {}
    for r in summary:
        if r["airfoil"] == "naca0012" and abs(float(r["alpha_deg"]) - 10) < 0.01 \
                and "serr" not in r["case"] and "coarse" not in r["case"] and "fine" not in r["case"]:
            if int(float(r["reynolds"])) in re_cl:
                continue
            re_rows.append((int(float(r["reynolds"])), r))
            re_cl[int(float(r["reynolds"]))] = float(r["Cl"])
            re_cd[int(float(r["reynolds"]))] = float(r["Cd"])
    re_rows.sort(key=lambda p: p[0])
    def re_sep(r):
        v = r.get("sep_upper", "")
        return f"{float(v):.2f}" if v and v not in ("None", "") else "$>0.99$ (attached)"
    re_label = {100000.0: "$10^5$", 300000.0: "$3\\times10^5$", 360000.0: "$3.6\\times10^5$", 1000000.0: "$10^6$"}
    m["ReTableRows"] = "\n    ".join(
        f"{re_label.get(float(re), f'{re:.0e}')} & {fmt(float(r['Cl']))} & {fmt(float(r['Cd']))} & "
        f"{fmt(float(r['Cl_Cd']), 1)} & {re_sep(r)} \\\\"
        for re, r in re_rows)
    if 1e5 in re_cl and 1e6 in re_cl:
        m["ReClDropPct"] = f"{100 * (re_cl[1e6] - re_cl[int(1e5)]) / re_cl[1e6]:.0f}"
        m["ReCdRisePct"] = f"{100 * (re_cd[int(1e5)] - re_cd[1e6]) / re_cd[1e6]:.0f}"
    else:
        m["ReClDropPct"] = m["ReCdRisePct"] = "n/a"
    if 3e5 in re_cl and 1e6 in re_cl:
        m["ReClDropThreePct"] = f"{100 * (re_cl[1e6] - re_cl[int(3e5)]) / re_cl[1e6]:.0f}"
        m["ReCdRiseThreePct"] = f"{100 * (re_cd[int(3e5)] - re_cd[1e6]) / re_cd[1e6]:.0f}"
    else:
        m["ReClDropThreePct"] = m["ReCdRiseThreePct"] = "n/a"
    sep_re1e5 = next((r for r in summary if r["case"] == "naca0012_aoa10_re1e5"), None)
    v = sep_re1e5.get("sep_upper", "") if sep_re1e5 else ""
    m["SepReLow"] = f"{float(v):.2f}" if v and v not in ("None", "") else "n/a"
    m["LdReSix"] = f"{re_cd.get(1e6, np.nan) and re_cl.get(1e6, 0) / re_cd.get(1e6, 1):.1f}" if 1e6 in re_cd else "n/a"
    m["LdReFive"] = f"{re_cl.get(int(1e5), 0) / re_cd.get(int(1e5), 1):.1f}" if 1e5 in re_cd else "n/a"

    gridhigh_rows = []
    for label, mesh_label, case in (("Medium", "4000 iters", "naca0012_aoa10"),
                                    ("Fine", "4000 iters", "naca0012_aoa10_fine"),
                                    ("Medium", "12000 iters", "naca0012_aoa10_long"),
                                    ("Medium, $R=50c$", "4000 iters", "naca0012_aoa10_farfield50"),
                                    ("Medium", "$\\alpha=12^\\circ$", "naca0012_aoa12"),
                                    ("Fine", "$\\alpha=12^\\circ$", "naca0012_aoa12_fine"),
                                    ("Medium", "$\\alpha=15^\\circ$", "naca0012_aoa15"),
                                    ("Fine", "$\\alpha=15^\\circ$", "naca0012_aoa15_fine"),
                                    ("Coarse", "$\\alpha=15^\\circ$", "naca0012_aoa15_coarse")):
        r = next((x for x in summary if x["case"] == case), None)
        if r:
            v = r.get("sep_upper", "")
            sep = f"{float(v):.2f}" if v and v not in ("None", "") else "attached"
            gridhigh_rows.append(f"{label} & {mesh_label} & {fmt(float(r['Cl']))} & {fmt(float(r['Cd']))} & {sep} \\\\")
    m["GridHighRows"] = "\n    ".join(gridhigh_rows) if gridhigh_rows else "--- & --- & --- & --- \\\\"

    serr_rows = []
    serr_drops = {}
    for a in (5, 10, 15):
        base = next((r for r in summary if r["case"] == f"naca0012_aoa{a}"), None)
        serr = next((r for r in summary if r["case"] == f"naca0012s_aoa{a}"), None)
        if base and serr:
            bl, bc = float(base["Cl"]), float(base["Cd"])
            sl, sc = float(serr["Cl"]), float(serr["Cd"])
            serr_rows.append(f"{a} & {fmt(bl)} & {fmt(sl)} & {fmt(bc)} & {fmt(sc)} \\\\")
            serr_drops[a] = (100 * (bc - sc) / bc, 100 * (bl - sl) / bl)
    m["SerrTableRows"] = "\n    ".join(serr_rows)
    m["SerrCdDropFive"] = f"{serr_drops[5][0]:.1f}" if 5 in serr_drops else "n/a"
    m["SerrClDropFive"] = f"{serr_drops[5][1]:.1f}" if 5 in serr_drops else "n/a"
    m["SerrCdDropFifteen"] = f"{serr_drops[15][0]:.1f}" if 15 in serr_drops else "n/a"
    m["SerrClRiseFifteen"] = f"{-serr_drops[15][1]:.1f}" if 15 in serr_drops else "n/a"
    m["SerrCdIncrFive"] = f"{-serr_drops[5][0]:.1f}" if 5 in serr_drops else "n/a"

    re_sweep_rows = []
    sweep_re_cl = {}
    for af in ("naca0012", "naca2412", "naca4412"):
        for re_val in (1e5, 3e5, 1e6):
            cl_by_a = {}
            for r in summary:
                if (r["airfoil"] == af and float(r["reynolds"]) == re_val
                        and all(s not in r["case"] for s in
                                ("serr", "coarse", "fine", "farfield", "long", "tub"))):
                    cl_by_a[int(float(r["alpha_deg"]))] = float(r["Cl"])
            if not cl_by_a:
                continue
            sweep_re_cl[(af, re_val)] = cl_by_a
            cells = " & ".join(fmt(cl_by_a.get(a)) for a in (0, 5, 10, 15, 18))
            re_sweep_rows.append(f"{af.upper()} & {re_label.get(re_val, f'{re_val:.0e}')} & {cells} \\\\")
    m["ReSweepTableRows"] = "\n    ".join(re_sweep_rows) or "--- \\\\"

    def sweep_delta(af, alpha, re_low=1e5, re_high=1e6):
        hi = sweep_re_cl.get((af, re_high), {}).get(alpha)
        lo = sweep_re_cl.get((af, re_low), {}).get(alpha)
        if hi is None or lo is None or abs(hi) < 1e-3:
            return None
        return 100 * (hi - lo) / abs(hi)
    d12 = sweep_delta("naca0012", 10)
    d24 = sweep_delta("naca2412", 10)
    d44 = sweep_delta("naca4412", 10)
    m["ReSweepClDropTenZeroZeroOneTwo"] = f"{d12:.1f}" if d12 is not None else "n/a"
    m["ReSweepClDropTenTwoFourOneTwo"] = f"{d24:.1f}" if d24 is not None else "n/a"
    m["ReSweepClDropTenFourFourOneTwo"] = f"{d44:.1f}" if d44 is not None else "n/a"

    tub_rows = []
    tub_stats = {}
    for a in (0, 5, 10, 15, 18):
        base = next((r for r in summary if r["case"] == f"naca0012_aoa{a}"), None)
        tub = next((r for r in summary if r["case"] == f"naca0012t_aoa{a}"), None)
        if base and tub:
            bl, bc = float(base["Cl"]), float(base["Cd"])
            tl, tc = float(tub["Cl"]), float(tub["Cd"])
            if abs(bl) < 0.02:
                dcl_txt = "--"
                dcl = None
            else:
                dcl = 100 * (tl - bl) / abs(bl)
                dcl_txt = f"{dcl:.1f}\\%"
            dcd = 100 * (tc - bc) / bc
            tub_rows.append(f"{a} & {fmt(bl)} & {fmt(tl)} & {dcl_txt} & "
                            f"{fmt(bc)} & {fmt(tc)} & {dcd:.1f}\\% \\\\")
            tub_stats[a] = (dcl, dcd)
    m["TubTableRows"] = "\n    ".join(tub_rows) or "--- \\\\"
    m["TubCdPctFive"] = f"{tub_stats[5][1]:.1f}" if 5 in tub_stats else "n/a"
    m["TubClPctFive"] = f"{tub_stats[5][0]:.1f}" if 5 in tub_stats and tub_stats[5][0] is not None else "n/a"
    m["TubCdPctTen"] = f"{tub_stats[10][1]:.1f}" if 10 in tub_stats else "n/a"
    m["TubCdPctFifteen"] = f"{tub_stats[15][1]:.1f}" if 15 in tub_stats else "n/a"
    m["TubClPctTen"] = f"{tub_stats[10][0]:.1f}" if 10 in tub_stats and tub_stats[10][0] is not None else "n/a"
    m["TubClPctFifteen"] = f"{tub_stats[15][0]:.1f}" if 15 in tub_stats and tub_stats[15][0] is not None else "n/a"
    tub15 = next((r for r in summary if r["case"] == "naca0012t_aoa15"), None)
    if tub15:
        v = tub15.get("sep_upper", "")
        m["TubSepFifteen"] = f"{float(v):.2f}" if v and v not in ("None", "") else "$>0.99$ (attached)"
    else:
        m["TubSepFifteen"] = "n/a"

    def tr_get(af, alpha, key):
        for r in summary:
            if (r["airfoil"] == af + "tr" and abs(float(r["alpha_deg"]) - alpha) < 0.01
                    and float(r["reynolds"]) == 1e5):
                v = r.get(key, "")
                return float(v) if v and v not in ("None", "") else None
        return None

    def ft_get(af, alpha, key):
        for r in summary:
            if (r["airfoil"] == af and abs(float(r["alpha_deg"]) - alpha) < 0.01
                    and float(r["reynolds"]) == 1e5
                    and all(s not in r["case"] for s in
                            ("serr", "coarse", "fine", "farfield", "long", "tub"))):
                v = r.get(key, "")
                return float(v) if v and v not in ("None", "") else None
        return None

    tr_rows = []
    tr_stats = {}
    for af, tag in (("naca0012", "ZeroZeroOneTwo"), ("naca2412", "TwoFourOneTwo"),
                    ("naca4412", "FourFourOneTwo")):
        for alpha in (0, 5, 10, 15, 18):
            cl_ft, cl_tr = ft_get(af, alpha, "Cl"), tr_get(af, alpha, "Cl")
            cd_ft, cd_tr = ft_get(af, alpha, "Cd"), tr_get(af, alpha, "Cd")
            sep_ft, sep_tr = ft_get(af, alpha, "sep_upper"), tr_get(af, alpha, "sep_upper")
            if cl_ft is None or cl_tr is None:
                continue
            dcl = 100 * (cl_tr - cl_ft) / abs(cl_ft) if abs(cl_ft) > 0.02 else None
            dcd = 100 * (cd_tr - cd_ft) / cd_ft if cd_ft else None
            def cell(v, nd=4):
                return f"{v:.{nd}f}" if v is not None else "--"
            def sep_cell(v):
                if v is None:
                    return "attached"
                if abs(v) < 0.005:
                    return "0.00 (LE)"
                return f"{v:.2f}"
            tr_rows.append(
                f"{af.upper()} & {alpha} & {cell(cl_ft)} & {cell(cl_tr)} & "
                f"{cell(dcl, 1)} & {cell(cd_ft)} & {cell(cd_tr)} & {cell(dcd, 1)} & "
                f"{sep_cell(sep_ft)} & {sep_cell(sep_tr)} \\\\")
            tr_stats[(af, alpha)] = (dcl, dcd)
    m["TrTableRows"] = "\n    ".join(tr_rows) or "--- \\\\"
    for af, tag in (("naca0012", "ZeroZeroOneTwo"), ("naca2412", "TwoFourOneTwo"),
                    ("naca4412", "FourFourOneTwo")):
        for alpha, aname in ((0, "Zero"), (5, "Five"), (10, "Ten"),
                             (15, "Fifteen"), (18, "Eighteen")):
            dcl, dcd = tr_stats.get((af, alpha), (None, None))
            m[f"TrClPct{aname}{tag}"] = f"{dcl:.1f}" if dcl is not None else "n/a"
            m[f"TrCdPct{aname}{tag}"] = f"{dcd:.1f}" if dcd is not None else "n/a"

    def sandia_mae(model, key):
        mpath = os.path.join(ROOT, "data", "validation", "metrics.json")
        if os.path.exists(mpath):
            met = json.load(open(mpath))
            v = met.get("low_re_sandia", {}).get(model, {}).get(key)
            return float(v) if v is not None else None
        return None

    ft_cl = sandia_mae("FT", "cl_mae_prestall")
    tr_cl = sandia_mae("TR", "cl_mae_prestall")
    ft_cd = sandia_mae("FT", "cd_mae_prestall")
    tr_cd = sandia_mae("TR", "cd_mae_prestall")
    m["SandiaClMaeFT"] = f"{ft_cl:.4f}" if ft_cl is not None else "n/a"
    m["SandiaClMaeTR"] = f"{tr_cl:.4f}" if tr_cl is not None else "n/a"
    m["SandiaCdMaeFT"] = f"{ft_cd:.4f}" if ft_cd is not None else "n/a"
    m["SandiaCdMaeTR"] = f"{tr_cd:.4f}" if tr_cd is not None else "n/a"
    if ft_cd and tr_cd:
        m["SandiaCdErrDropPct"] = f"{100 * (ft_cd - tr_cd) / ft_cd:.0f}"
    else:
        m["SandiaCdErrDropPct"] = "n/a"

    def get_at(af, alpha, re, key, tr=False):
        for r in summary:
            if (r["airfoil"] == af + ("tr" if tr else "")
                    and abs(float(r["alpha_deg"]) - alpha) < 0.01
                    and float(r["reynolds"]) == re
                    and all(s not in r["case"] for s in
                            ("serr", "coarse", "fine", "farfield", "long", "tub"))):
                v = r.get(key, "")
                return float(v) if v and v not in ("None", "") else None
        return None

    m["SandiaClTenExp"] = "0.981"
    m["SandiaClTwelveExp"] = "0.483"
    cl12_ft = get_at("naca0012", 12, 3.6e5, "Cl")
    m["SandiaClTwelveCFD"] = f"{cl12_ft:.2f}" if cl12_ft is not None else "n/a"
    m["SandiaCdFifteenExp"] = "0.190"
    cd15_ft = get_at("naca0012", 15, 3.6e5, "Cd")
    m["SandiaCdFifteenCFD"] = f"{cd15_ft:.3f}" if cd15_ft is not None else "n/a"

    cl_ft0 = get_at("naca4412", 0, 3e5, "Cl")
    cl_tr0 = get_at("naca4412", 0, 3e5, "Cl", tr=True)
    cd_ft0 = get_at("naca4412", 0, 3e5, "Cd")
    cd_tr0 = get_at("naca4412", 0, 3e5, "Cd", tr=True)
    if cl_ft0 and cl_tr0:
        m["TrClPctZeroFourFourOneTwoReThr"] = f"{100 * (cl_tr0 - cl_ft0) / abs(cl_ft0):.1f}"
    if cd_ft0 and cd_tr0:
        m["TrCdDropZeroFourFourOneTwoReThr"] = f"{100 * (cd_tr0 - cd_ft0) / cd_ft0:.1f}"

    allforce_rows = []
    for af in ("naca0012", "naca2412", "naca4412"):
        for r in core:
            if r["airfoil"] != af:
                continue
            allforce_rows.append(
                f"{af.upper()} & {float(r['alpha_deg']):.0f} & {fmt(float(r['Cl']))} & "
                f"{fmt(float(r['Cd']))} & {fmt(float(r['Cl_Cd']), 1)} \\\\")
    m["AllForcesRows"] = "\n    ".join(allforce_rows) or "--- \\\\"

    def mesh_quality(logname):
        path = os.path.join(ROOT, "logs", logname + ".log")
        if not os.path.exists(path):
            return None
        txt = open(path, errors="ignore").read()
        no = re.search(r"non-orthogonality Max:\s*([\d.]+)\s+average:\s*([\d.]+)", txt)
        sk = re.search(r"Max skewness\s*=\s*([\d.]+)", txt)
        if not no or not sk:
            return None
        return no.group(1), no.group(2), sk.group(1)
    for dens, case in (("Coarse", "naca0012_aoa5_coarse"),
                       ("Medium", "naca0012_aoa5"),
                       ("Fine", "naca0012_aoa5_fine")):
        q = mesh_quality(case)
        m[f"NonOrth{dens}"] = q[0] if q else "n/a"
        m[f"NonOrthAvg{dens}"] = q[1] if q else "n/a"
        m[f"Skew{dens}"] = q[2] if q else "n/a"
    try:
        cl_med, cl_fine = float(m["ClMedium"]), float(m["ClFine"])
        m["GciPct"] = f"{100 * abs(cl_fine - cl_med) / abs(cl_fine):.2f}"
    except (KeyError, TypeError, ValueError):
        m["GciPct"] = "n/a"

    with open(os.path.join(LATEX, "numbers.tex"), "w") as f:
        for k, v in m.items():
            f.write(f"\\newcommand{{\\{k}}}{{{v}}}\n")

    conv_fig = os.path.join(LATEX, "figures", "appendix_convergence.png")
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    ax2 = ax.twinx()
    for case, color in (("naca0012_aoa5", "tab:blue"), ("naca0012_aoa15", "tab:red")):
        path = os.path.join(RESULTS, "convergence", case + ".csv")
        if os.path.exists(path):
            rows = list(csv.DictReader(open(path)))
            t = [float(r["time"]) for r in rows]
            cl = [float(r["Cl"]) for r in rows]
            cd = [float(r["Cd"]) for r in rows]
            bad = [i for i in range(len(t))
                   if abs(cl[i]) >= 3 * abs(cl[-1]) or abs(cd[i]) >= 3 * abs(cd[-1])]
            start = (bad[-1] + 1) if bad else 0
            ax.plot(t[start:], cl[start:], color=color, label=f"{case} $C_L$")
            ax2.plot(t[start:], cd[start:], "--", color=color, alpha=0.6,
                     label=f"{case} $C_D$")
    ax.set_xlabel("Iteration")
    ax.set_ylabel("$C_L$")
    ax2.set_ylabel("$C_D$")
    ax.set_title("Force convergence history")
    handles = ax.get_legend_handles_labels()[0] + ax2.get_legend_handles_labels()[0]
    labels = ax.get_legend_handles_labels()[1] + ax2.get_legend_handles_labels()[1]
    ax.legend(handles, labels, fontsize=8, loc="center right")
    fig.savefig(conv_fig)
    plt.close(fig)

    if "--no-compile" in sys.argv:
        print("[done] numbers.tex + appendix figure written")
        return
    env = dict(os.environ)
    result = subprocess.run([TECTONIC, "main.tex"], cwd=LATEX, capture_output=True, text=True)
    print(result.stdout[-3000:])
    if result.returncode != 0:
        print(result.stderr[-3000:])
        sys.exit(result.returncode)
    print(f"[done] {LATEX}/main.pdf")


if __name__ == "__main__":
    import numpy as np
    main()
