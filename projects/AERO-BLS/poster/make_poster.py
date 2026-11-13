#!/usr/bin/env python3
"""Build the A0 landscape conference poster (poster/poster.html).

Static sections use the publication figures; the interactive velocity
section is a Plotly contour built from the foamToVTK field data when
available (falls back to a static figure otherwise).
"""
import csv
import glob
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
FIGS = os.path.join(ROOT, "latex", "figures")
RESULTS = os.path.join(ROOT, "data", "results")


def load_summary():
    path = os.path.join(RESULTS, "summary.csv")
    if not os.path.exists(path):
        return []
    with open(path) as f:
        return list(csv.DictReader(f))


def plotly_velocity():
    """Return (html_fragment, used_plotly)."""
    case = "naca0012_aoa15"
    hits = (sorted(glob.glob(os.path.join(ROOT, "data", "fields", case, "**", "*.vtu"),
                             recursive=True))
            + sorted(glob.glob(os.path.join(ROOT, "data", "fields", case, "**", "*.vtk"),
                               recursive=True)))
    if not hits:
        return ("<p>Field data not available; see figure fig6 in the paper.</p>", False)
    try:
        import meshio
        import numpy as np
    except ImportError:
        return ("<p>meshio not installed; see figure fig6 in the paper.</p>", False)
    mesh = meshio.read(hits[-1])
    if "hexahedron" not in mesh.cells_dict:
        return ("<p>Unexpected field format.</p>", False)
    pts = mesh.points[mesh.cells_dict["hexahedron"]].mean(axis=1)
    umag = None
    for name, arrs in mesh.cell_data.items():
        for arr in arrs:
            if name == "UMag" and arr.shape[0] == pts.shape[0]:
                umag = arr
    if umag is None and "U" in mesh.cell_data:
        U = mesh.cell_data["U"][0]
        umag = np.linalg.norm(U[:, :3], axis=1)
    if umag is None:
        return ("<p>No velocity data.</p>", False)
    x, y, v = pts[:, 0].tolist(), pts[:, 1].tolist(), umag.tolist()
    data = {
        "type": "contour", "x": x, "y": y, "z": v,
        "colorscale": "Viridis", "contours": {"coloring": "heatmap"},
        "line": {"width": 0},
    }
    layout = {
        "title": "Velocity magnitude, NACA 0012 at 15 deg (drag to rotate, scroll to zoom)",
        "xaxis": {"title": "x/c", "range": [-0.2, 1.6]},
        "yaxis": {"title": "y/c", "scaleanchor": "x", "range": [-0.6, 0.6]},
        "width": 1150, "height": 520, "margin": {"l": 60, "r": 20, "t": 50, "b": 50},
    }
    return (f"<div id='velplot'></div>"
            f"<script>Plotly.newPlot('velplot', [{json.dumps(data)}], {json.dumps(layout)});</script>",
            True)


def val_table():
    path = os.path.join(ROOT, "data", "validation", "validation_table.csv")
    if not os.path.exists(path):
        return "<p>Validation table not yet generated.</p>"
    rows = list(csv.DictReader(open(path)))
    head = ("<tr><th>α [°]</th><th>C<sub>L</sub> CFD</th><th>C<sub>L</sub> Ladson</th>"
            "<th>err [%]</th><th>C<sub>D</sub> CFD</th><th>C<sub>D</sub> Ladson</th></tr>")
    body = ""
    for r in rows:
        body += (f"<tr><td>{float(r['alpha_deg']):.0f}</td>"
                 f"<td>{float(r['Cl_CFD']):.3f}</td>"
                 f"<td>{float(r['Cl_Ladson']):.3f}</td>"
                 f"<td>{float(r['Cl_err_Ladson_pct']):.1f}</td>"
                 f"<td>{float(r['Cd_CFD']):.4f}</td>"
                 f"<td>{float(r['Cd_Ladson']):.4f}</td></tr>")
    return f"<table>{head}{body}</table>"


def main():
    plotly_html, used_plotly = plotly_velocity()
    metrics = json.load(open(os.path.join(ROOT, "data", "validation", "metrics.json")))
    cl_mae_pct = metrics["naca0012_ladson"]["cl_mae_pre_stall_pct"]
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>AERO-BLS: Turbulent Boundary Layer Separation on NACA Airfoils</title>
<script src="https://cdn.plot.ly/plotly-2.32.0.min.js"></script>
<style>
  @page {{ size: A0 landscape; margin: 0; }}
  * {{ box-sizing: border-box; }}
  body {{ margin: 0; font-family: 'Helvetica Neue', Arial, sans-serif;
         background: #0e2a47; color: #10202e; }}
  .poster {{ width: 1189mm; height: 841mm; padding: 24mm; display: grid;
            grid-template-columns: 1.05fr 1fr 1fr; grid-template-rows: auto 1fr 1fr;
            gap: 10mm; background: linear-gradient(160deg,#0e2a47 0%,#123a5e 100%); }}
  .banner {{ grid-column: 1 / 4; background: #ffffff; border-radius: 6mm;
            padding: 10mm 14mm; display: flex; align-items: center; gap: 14mm; }}
  .banner h1 {{ font-size: 52pt; margin: 0; color: #0e2a47; }}
  .banner p {{ font-size: 20pt; margin: 4mm 0 0; color: #345; }}
  .card {{ background: #ffffff; border-radius: 6mm; padding: 8mm 10mm;
          overflow: hidden; }}
  .card h2 {{ font-size: 28pt; margin: 0 0 5mm; color: #0e2a47;
             border-bottom: 1.2mm solid #1f6fb2; padding-bottom: 2mm; }}
  .card p, .card li {{ font-size: 16pt; line-height: 1.35; }}
  .card img {{ max-width: 100%; }}
  ul {{ margin: 2mm 0; padding-left: 7mm; }}
  table {{ border-collapse: collapse; width: 100%; font-size: 13pt; }}
  th, td {{ border: 0.3mm solid #b9c6d2; padding: 1.6mm 2.4mm; text-align: center; }}
  th {{ background: #e8f0f7; }}
  .eq {{ font-size: 17pt; background: #f2f7fb; padding: 3mm 5mm; border-radius: 3mm; }}
  .foot {{ grid-column: 1 / 4; color: #cfe2f3; font-size: 13pt; text-align: center; }}
  .foot a {{ color: #cfe2f3; }}
</style>
</head>
<body>
<div class="poster">
  <div class="banner">
    <div>
      <h1>Turbulent Boundary Layer Separation &amp; Passive Drag Reduction on NACA Airfoils</h1>
      <p>A validated 2D RANS study with the k–ω SST closure — OpenFOAM v2412 · Re = 10<sup>5</sup>–10<sup>6</sup> · α = 0–18°</p>
    </div>
  </div>

  <div class="card">
    <h2>1 · Motivation &amp; Physics</h2>
    <p>Separation onset limits maximum lift and drives the drag rise near stall on
       small-UAS and wind-turbine sections.</p>
    <div class="eq">∂U<sub>i</sub>/∂x<sub>i</sub> = 0 &nbsp;&nbsp;·&nbsp;&nbsp;
       U<sub>j</sub>∂U<sub>i</sub>/∂x<sub>j</sub> = −(1/ρ)∂p/∂x<sub>i</sub> + ∂/∂x<sub>j</sub>[ν(∂U<sub>i</sub>+∂U<sub>j</sub>) − u′<sub>i</sub>u′<sub>j</sub>]</div>
    <ul>
      <li><b>Closure:</b> k–ω SST (Menter) — best-in-class for adverse pressure gradients; γ–Re<sub>θt</sub> repeats for transition physics</li>
      <li><b>Profiles:</b> NACA 0012 · 2412 · 4412 (analytic, 638 pts)</li>
      <li><b>Mesh:</b> polar O-grid, y⁺ ≈ 0.4, farfield 20c, 34k cells</li>
      <li><b>Solver:</b> steady SIMPLEC, fully turbulent, Tu = 0.1%</li>
    </ul>
  </div>

  <div class="card">
    <h2>2 · Method Pipeline</h2>
    <ul>
      <li>Analytic NACA 4-digit geometry → serrated / tubercle variants</li>
      <li>Structured O-grid generator → Gmsh v2.2 → gmshToFoam</li>
      <li>91 scripted OpenFOAM cases (21 core + full Re sweep + grid/serration/tubercle + 25 γ–Re<sub>θt</sub> transition repeats)</li>
      <li>Forces, Cp, Cf, separation from wall shear reversal</li>
      <li>Validation: NASA TMR NACA 0012 (Ladson, Gregory, Abbott) + CFL3D SST</li>
    </ul>
    <img src="../latex/figures/fig5_cp_distribution.png" alt="Cp distributions">
  </div>

  <div class="card">
    <h2>3 · Validation vs NASA Data</h2>
    {val_table()}
    <p style="font-size:13pt">CFD: Re = 10<sup>6</sup>, fully turbulent. Ladson: Re = 6×10<sup>6</sup>, tripped.
       Pre-stall C<sub>L</sub> MAE ≈ {cl_mae_pct:.1f}%.</p>
    <img src="../latex/figures/fig9_validation.png" alt="Validation">
  </div>

  <div class="card">
    <h2>4 · Separation &amp; Stall</h2>
    <img src="../latex/figures/fig8_wall_shear_stress.png" alt="Cf distributions">
    <p>Separation detected from sustained C<sub>f</sub> &lt; 0:
       attached to α = 10°, x/c ≈ 0.20 at α = 15°, with a five-fold drag rise.</p>
  </div>

  <div class="card">
    <h2>5 · Interactive Velocity Field (α = 15°)</h2>
    {plotly_html}
  </div>

  <div class="card">
    <h2>6 · Conclusions</h2>
    <ul>
      <li>Pre-stall lift matches tripped experiment within ≈ few % (MAE); drag biased high by the fully turbulent assumption — quantified and consistent.</li>
      <li>Camber orders maximum lift: 4412 &gt; 2412 &gt; 0012; stall softens with camber.</li>
      <li>Lowering Re from 10<sup>6</sup> to 10<sup>5</sup> cuts lift and raises drag at fixed incidence; the stall-angle shift is muted by the fully turbulent closure.</li>
      <li>γ–Re<sub>θt</sub> transition repeats (25 cases, identical grids): attached-range drag −11–22% at Re = 10<sup>5</sup>, laminar separation bubbles at x/c = 0.06–0.75 that the fully turbulent model cannot produce; near stall the closures diverge (ΔC<sub>L</sub> −32% to +30%) as bubble burst reshapes separation.</li>
      <li>2D serrated TE <i>increases</i> attached-flow drag (+8.6% at α = 5°) but cuts drag −13.6% at α = 15° via base-pressure recovery; the 2D tubercle surrogate is a null result — both devices need spanwise (3D) physics.</li>
    </ul>
    <img src="../latex/figures/fig2_cl_vs_alpha.png" alt="Lift curves">
    <img src="../latex/figures/fig13_transition.png" alt="Transition vs fully turbulent">
  </div>

  <div class="foot">AERO-BLS — fully scripted reproducible workflow (scripts 00–08) · OpenFOAM v2412 ·
    Data: NASA Turbulence Modeling Resource, NASA TM-4074, R&amp;M 3726, Abbott &amp; von Doenhoff (1959)</div>
</div>
</body>
</html>
"""
    out = os.path.join(HERE, "poster.html")
    with open(out, "w") as f:
        f.write(html)
    print(f"[done] {out} (plotly={used_plotly})")


if __name__ == "__main__":
    main()
