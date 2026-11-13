# AERO-BLS — Turbulent Boundary Layer Separation & Passive Drag Reduction on NACA Airfoils

A fully scripted, reproducible CFD study of turbulent boundary layer separation,
stall, and passive drag reduction on the NACA 0012, 2412, and 4412 sections.

- **Physics:** 2D incompressible RANS, k–ω SST closure (OpenFOAM v2412, `simpleFoam`, SIMPLEC)
- **Mesh:** structured polar O-grid, y⁺ ≈ 0.4 at Re = 10⁶, farfield at 20 chords
- **Matrix:** 15 core cases (3 airfoils × α = 0°, 5°, 10°, 15°, 18° at Re = 10⁶)
  + Reynolds sensitivity (10⁵, 3×10⁵) + grid independence (coarse/medium/fine)
  + trailing-edge serration cases
- **Validation:** NASA Turbulence Modeling Resource NACA 0012 datasets
  (Ladson TM-4074 tripped forces, Gregory & O'Reilly R&M 3726, Abbott & von Doenhoff)
  and NASA CFL3D SST reference results

## Repository layout

```
projects/AERO-BLS/
├── airfoils/            # analytic NACA 4-digit coordinates (638 pts)
├── mesh/                # generated Gmsh v2.2 hex meshes + metadata
├── openfoam_cases/      # base_case template + 22 staged cases
├── scripts/             # 00–08 pipeline (see below)
├── data/
│   ├── experimental/    # NASA TMR .dat reference files
│   ├── results/         # summary.csv, cp/, cf/, convergence/
│   ├── validation/      # validation_table.csv, metrics.json
│   └── fields/          # foamToVTK output for contour figures
├── latex/               # main.tex, numbers.tex (auto-generated), figures/
├── poster/              # poster.html (A0 landscape, Plotly) + builder
└── logs/                # per-case solver logs
```

## Pipeline

| Script | Purpose |
|--------|---------|
| `00_setup_wsl.sh` | One-time environment install (WSL2/Ubuntu or Debian) |
| `01_fetch_airfoil_coords.py` | UIUC fetch with analytic NACA 4-digit fallback |
| `02_generate_mesh.py` | Polar O-grid → Gmsh v2.2 hex mesh (`--serration`, `--tubercle`, `--density`) |
| `03_prepare_cases.py` | Stage `base_case`, substitute BC tokens per case |
| `04_run_all_cases.sh` | Batch runner over `run_manifest.txt` (NPAR parallel) |
| `05_postprocess.py` | Forces (tail-averaged), Cp, Cf, separation, residuals → `data/results/` |
| `06_generate_figures.py` | 9 publication figures → `latex/figures/` |
| `07_validate_vs_experiment.py` | CFD vs Ladson/Gregory/Abbott + CFL3D SST → `data/validation/` |
| `08_build_report.py` | `numbers.tex` macros + appendix figure + tectonic compile |

## Reproducing from scratch

```bash
# 1. Environment (WSL2 Ubuntu / Debian; needs internet)
bash scripts/00_setup_wsl.sh

# 2. Toolchain in PATH (or source your OpenFOAM bashrc)
export MAMBA_ROOT_PREFIX=$HOME/tools/mamba/root
source $MAMBA_ROOT_PREFIX/envs/main/etc/bashrc   # skip if using distro OpenFOAM

# 3. Geometry + meshes
python3 scripts/01_fetch_airfoil_coords.py
python3 scripts/02_generate_mesh.py --airfoil naca0012 --output mesh/naca0012_medium.msh --density medium
#    (repeat for 2412/4412, coarse/fine, --serration variant)

# 4. Stage + run all cases
python3 scripts/03_prepare_cases.py
bash scripts/04_run_all_cases.sh 2            # 2 = parallel workers

# 5. Post-process, validate, figures, report
python3 scripts/05_postprocess.py
#    optional field contours: foamToVTK -latestTime -case openfoam_cases/naca0012_aoa15
#    then copy VTK output to data/fields/naca0012_aoa15/
python3 scripts/07_validate_vs_experiment.py
python3 scripts/06_generate_figures.py
python3 scripts/08_build_report.py            # compiles latex/main.pdf

# 6. Poster
python3 poster/make_poster.py                 # poster/poster.html (A0 landscape)
```

## Key results (see `latex/main.pdf` for full detail)

- Pre-stall lift on NACA 0012 within a few % (MAE) of the tripped Ladson data;
  lift curve bracketed by the three experimental datasets and close to NASA's
  CFL3D SST reference.
- Separation onset by upstream migration of the C_f = 0 point: attached through
  α = 10°, x/c ≈ 0.2 at α = 15°, with a ~5× drag rise.
- Maximum-lift ordering 4412 > 2412 > 0012; stall softens with camber.
- Re = 10⁵ vs 10⁶ at α = 10°: lower lift, higher drag, earlier separation.
- Serrated trailing edge reduces drag at moderate α in the 2D surrogate.

## Known limitations (documented in the paper)

- 2D steady RANS: no spanwise stall cells, no post-stall unsteadiness.
- Fully turbulent boundary layer: drag biased high vs (even tripped) experiment;
  no laminar separation bubble physics at Re ≤ 10⁶.
- Serrations/tubercles are intrinsically 3D devices; 2D results are surrogates.

## Data provenance

All experimental reference files in `data/experimental/` are U.S. Government
works (public domain) digitized and distributed by the NASA Langley
Turbulence Modeling Resource: <https://tmbwg.github.io/turbmodels/naca0012_val.html>.
