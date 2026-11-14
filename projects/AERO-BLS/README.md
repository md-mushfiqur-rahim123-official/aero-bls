# AERO-BLS — Transition-Sensitive RANS Modelling of Boundary-Layer Separation on NACA Airfoils at Low Reynolds Number

A fully scripted, reproducible CFD study of boundary-layer separation, transition,
and stall on the NACA 0012, 2412, and 4412 sections, with the passive-device
question retained as an explicit **2D negative control** (Appendix D of the paper).

- **Physics:** 2D incompressible RANS, k–ω SST and γ–Re_θt (Langtry–Menter)
  closures (OpenFOAM v2412); steady SIMPLEC (`simpleFoam`) + a URANS
  verification subset (`pimpleFoam`) for the near-stall unsteady check
- **Mesh:** structured polar O-grid, y⁺ ≈ 0.4 at Re = 10⁶, farfield at 20 chords
  (50-chord verification mesh for the URANS subset)
- **Matrix:** 21 core cases (3 airfoils × α = 0°, 5°, 10°, 12°, 15°, 18°, 20° at Re = 10⁶)
  + full Reynolds sweep (Re = 10⁵, 3×10⁵ for **all three** airfoils, α = 0–18°)
  + grid independence (coarse/medium/fine at α = 5°, full family at the stall
    angle α = 15°, fine repeat at α = 18°) + farfield/convergence verification
  + trailing-edge serration and leading-edge tubercle (2D surrogate) cases
  + transition-sensitive repeats with the γ–Re_θt (Langtry–Menter) closure:
    NACA 0012/2412/4412 × α = 0–18° × Re = 10⁵/3×10⁵/10⁶ (full factorial, 45 cases),
    warm-started from the converged fully turbulent siblings
  + low-Reynolds free-transition validation at Re = 3.6×10⁵ vs Sandia SAND80-2114
    (Sheldahl & Klimas): NACA 0012 fully turbulent + γ–Re_θt, α = 0–15° (10 cases)
    — the **primary (Reynolds-matched) validation** of the study
  + URANS verification subset (pimpleFoam, k–ω SST, 50-chord farfield,
    33,740-cell mesh): NACA 0012 γ–Re_θt-warm-started, α = 12°/15° at
    Re = 3.6×10⁵, two-stage integration (dt = 10⁻⁵ spin-up, Co ≤ 20 production)
  — **125 cases total** (123 steady + 2 URANS)
- **Validation:** primary = Reynolds-matched Sandia SAND80-2114 (Re = 3.6×10⁵);
  NASA Turbulence Modeling Resource datasets (Ladson TM-4074 tripped forces,
  Gregory & O'Reilly R&M 3726, Abbott & von Doenhoff) and NASA CFL3D SST
  reference results are used **qualitatively only** (Re mismatch 3–6×)
- **Reviewer response:** point-by-point response letter in `docs/reviewer_response.md`

## Repository layout

```
projects/AERO-BLS/
├── airfoils/            # analytic NACA 4-digit coordinates (638 pts)
├── mesh/                # generated Gmsh v2.2 hex meshes + metadata
├── openfoam_cases/      # base_case templates + 123 staged cases (incl. URANS subset)
├── scripts/             # 00–10 pipeline (see below)
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
| `06_generate_figures.py` | 16 publication figures → `latex/figures/` |
| `07_validate_vs_experiment.py` | CFD vs Sandia (primary, Re-matched) + Ladson/Gregory/Abbott (qualitative) + CFL3D SST → `data/validation/` |
| `08_build_report.py` | `numbers.tex` macros + appendix figure + tectonic compile |
| `09_run_urans.sh` | URANS verification subset runner (stage 1 spin-up + stage 2 Co-limited production) |
| `10_postprocess_urans.py` | URANS production-window averaging → `data/results/urans_summary.csv` + fig15 |

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
#    (repeat for 2412/4412, coarse/fine, --serration/--tubercle variants)
#    NOTE: cambered airfoils need the O-grid centre on the camber line to keep
#    mesh quality acceptable:  add --centre-y 0.019 (2412) / --centre-y 0.0375 (4412)

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

- **Primary validation (Re-matched, Re = 3.6×10⁵, free transition):** the
  γ–Re_θt closure cuts attached-range drag error vs Sandia SAND80-2114 by 33%
  relative to fully turbulent k–ω SST (Cd MAE 0.0102 → 0.0068); lift MAE ≈ 0.04
  for both closures over the attached range.
- Steady RANS lift on NACA 0012 within 5–6% of the (higher-Re) tripped Ladson
  data over the attached range — used qualitatively only, given the 3–6×
  Reynolds mismatch documented in the paper.
- Separation onset by upstream migration of the C_f = 0 point: attached through
  α = 10°–12°, x/c ≈ 0.2 at α = 15° on the production grid, with a ~5× drag rise.
  Fine-grid repeats (α = 12°, 15°) show the medium grid over-separates at high
  lift: the fine grid stays attached at α = 15° (C_L 1.24 vs 0.99), so modelled
  onset is grid-sensitive — reported as a bracketed range, not a point value.
- URANS verification subset (α = 12°/15°, Re = 3.6×10⁵): time-averaged loads
  and their fluctuation envelopes quantify how much of the steady-RANS/steady-
  experiment discrepancy at near-stall incidence is unsteady in origin
  (`data/results/urans_summary.csv`, fig15).
- Maximum-lift ordering 4412 (1.476 @ 15°) > 2412 (1.351 @ 18°) > 0012
  (plateau ≈ 1.05 beyond 12°); stall softens with camber.
- Full Reynolds sweep (10⁵, 3×10⁵, all three airfoils): attached-range lift
  drops 6–9% at α = 10° (ordered by camber), but near stall the effect is
  dramatic — at Re = 10⁵ all three sections are fully separated by α = 18°
  (separation at the leading edge, C_L max falls 30–45%).
- Passive devices, 2D negative control (Appendix D): serrated TE **increases**
  drag +8.6% at α = 5° (attached) and reduces it −13.6% at α = 15° (separated,
  base-pressure recovery); tubercle LE surrogate is penalised at every
  incidence (drag up to +231%, lift up to −51%). Both device classes depend on
  spanwise (3D) physics, so these 2D results are reported only as negative
  controls demonstrating that 2D RANS cannot assess them.

## Known limitations (documented in the paper)

- 2D RANS: no spanwise stall cells; URANS subset is 2D as well (verification,
  not a 3D stall prediction).
- Fully turbulent baseline: drag biased high vs (even tripped) experiment;
  the γ–Re_θt subset addresses laminar-run/LSB physics at Re ≤ 10⁶.
- Post-stall bubble bursting remains outside RANS validity; URANS results near
  stall are qualitative (Co ≤ 20).

## Data provenance

All experimental reference files in `data/experimental/` are U.S. Government
works (public domain) digitized and distributed by the NASA Langley
Turbulence Modeling Resource: <https://tmbwg.github.io/turbmodels/naca0012_val.html>.
