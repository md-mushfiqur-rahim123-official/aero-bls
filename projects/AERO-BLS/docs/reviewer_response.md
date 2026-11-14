# Response to Reviewer — Point-by-Point

**Manuscript:** Turbulent Boundary Layer Separation and Passive Drag Reduction on NACA 4-Digit Airfoils: A Validated 2D RANS Study with the k–ω SST Model
**Decision under review:** 42/100 — Major Revision / Reject (top-tier)
**Revision outcome:** Manuscript reframed as a low-Reynolds-number transition/LSB study; 2D device sections demoted to a negative-control appendix; URANS campaign on the 50-chord-farfield grid added; validation re-anchored to Reynolds-matched data.

We thank the reviewer for a rigorous and constructive assessment. We accept the central criticisms: the original framing overstated the contribution, the validation mixed Reynolds numbers, and the stall-region results were presented with more confidence than the grid and solver support. Below we respond to each point and describe the corresponding change.

---

## Point 1 — Lack of scientific novelty

**Reviewer:** 2D RANS (k–ω SST) on NACA 0012/2412/4412 is a standard benchmark exercise; no new fluid mechanics.

**Response: Accepted.** We do not claim new separation physics for the NACA 4-digit series. The original framing ("turbulent boundary layer separation and passive drag reduction") promised more novelty than the content delivers. We have reframed the manuscript around the element that *does* carry weight: a systematic, fully scripted comparison of **fully turbulent SST vs. γ–Re_θt transition-sensitive closure across a full 3-airfoil × 3-Reynolds × α factorial (121 cases), anchored to a Reynolds-matched free-transition experiment (SAND80-2114)**, with all cases, scripts, and meshes released. This is positioned as a verification/validation and reference-dataset contribution at the level appropriate for an applied journal (Aerospace Science and Technology / Journal of Aircraft / JAFM), not a top-tier novelty claim. We agree the paper is not suitable for JFM/PoF and no longer target them.

## Point 2 — 3D devices tested with a 2D model

**Reviewer:** Tubercles and serrations depend entirely on spanwise vortices; running them in 2D and reporting "+231% drag" is methodologically wrong.

**Response: Partially accepted.** The submitted manuscript already flagged this explicitly (both device classes "derive their benefit from spanwise flow physics"), but the reviewer is right that a *negative* result obtained with a *mismatched* dimensionality does not belong in the main results. **Change:** the serration and tubercle sections are removed from the main results and condensed into an appendix explicitly labelled **"2D negative control"**, with the stated purpose of quantifying *how wrong* a 2D surrogate gets (it predicts penalties where 3D experiments show benefits), not of evaluating the devices. No claim about device effectiveness is made from 2D results anywhere in the revised text.

## Point 3 — Grid independence fails at stall

**Reviewer:** At α = 15° the medium and fine grids differ by ~25% in C_L (0.99 vs 1.24); the core topic (separation) is exactly where the grid is not converged.

**Response: Accepted.** The GCI of 0.58% applies only to attached flow (α = 5°). At stall onset, steady RANS separation location is strongly grid-sensitive — this is consistent with the NASA TMR guidance for this case, but we presented the production-grid stall envelope too confidently. **Changes:** (i) all quantitative stall-region statements now cite the **fine grid as primary**; (ii) the medium-grid over-separation is presented as a *finding* about the modelling chain, not a physical prediction; (iii) the new URANS campaign (Point 4) uses the 50-chord-farfield grid exclusively; (iv) the abstract and conclusions no longer quote stall-region force coefficients from the medium grid.

## Point 4 — Steady solver for unsteady stall / bubble burst

**Reviewer:** Stall and laminar-separation-bubble bursting are unsteady, 3D events; averaging the last 20% of a non-converged simpleFoam run is unacceptable; URANS or LES needed.

**Response: Accepted.** **Change:** a URANS (pimpleFoam) verification campaign was added for the low-Reynolds-number regime where the LSB physics is the point. The delivered subset covers the NACA 0012 at the near-stall angles α = 12° and 15° at Re = 3.6×10⁵ (the Reynolds-matched validation condition) on the 50-chord-farfield grid, warm-started from the converged γ–Re_θt steady siblings, with two-stage integration (fixed Δt = 10⁻⁵ settling, then Co ≤ 20 production to t ≈ 1.4 flow passes) and time-averaged force statistics ± standard deviation. This subset is deliberately scoped as a verification that the near-stall steady states are unsteady and that mean loads shift accordingly; extending it to the full airfoil × α × Re matrix and to resolved bubble dynamics is stated as future work. We do not claim LES-grade accuracy; the revised limitations section states that bubble *bursting* in experiment is only partially reproduced and that 3D instabilities are outside URANS 2D. The steady simpleFoam results are retained only for the attached range, where the steady assumption is defensible and the residual histories show genuine convergence (initial residuals < 1e-6, force drift < 0.5% over the averaging window).

## Point 5 — ~100% drag overprediction vs tripped experiment

**Reviewer:** Baseline drag is 2× the tripped experiment; a paper titled "drag reduction" cannot rest on a baseline that wrong.

**Response: Partially accepted — diagnosis corrected, framing fixed.** Two distinct causes, previously conflated:
1. **Reynolds mismatch (Point 6):** CFD at Re = 10⁶ was scored against Ladson tripped data at Re = 6×10⁶ and Gregory & O'Reilly at 2.88×10⁶. Part of the gap is physical (Reynolds effect on skin friction), not model error.
2. **Trip vs. free transition:** a fully turbulent CFD closure has no laminar run, so against *tripped* data the comparison is apples-to-apples in boundary-layer state, but drag magnitude still differs through the Re gap.

**Change:** the factorial validation MAE against the Re-mismatched datasets is **withdrawn as a headline metric**. The primary quantitative validation is now the **Reynolds-matched Sandia SAND80-2114 comparison at Re = 3.6×10⁵ with free transition** (Cl MAE 0.029 FT / 0.040 TR; Cd MAE 0.010 FT / 0.007 TR — a 33% drag-error reduction from the transition closure). The Ladson/Gregory/Abbott comparisons are retained as qualitative curve-shape checks, labelled with their Reynolds numbers, and the drag-budget section explains the laminar-run contribution explicitly.

## Point 6 — Reynolds number mismatch in validation

**Reviewer:** Direct MAE against data at 6× and 2.88× the simulated Re is erroneous.

**Response: Accepted in full.** See Point 5. All quantitative error metrics in the revised manuscript are computed at matched Reynolds number. Where datasets at other Re are shown in figures, they are labelled and used only for trend comparison.

## Point 7 — Mesh quality (non-orthogonality 82.7–86.2°) and farfield sensitivity

**Reviewer:** Polar O-grid cells near the stagnation region are highly non-orthogonal (numerical diffusion in OpenFOAM), and 20c→50c farfield moves lift by 2.7% — the domain was not converged either.

**Response: Accepted with scope note.** The polar O-grid concentrates non-orthogonality at the nose and in the wake cut; bounded schemes (limited linear divergence, limitedLinear convection) keep the solution stable, but the reviewer is right that this degrades accuracy where it matters most. **Changes:** (i) the URANS campaign runs on the refined farfield (50c) grid at medium surface density; (ii) the farfield study is promoted from appendix detail to the mesh-quality discussion; (iii) a wake-aligned C-topology is listed as the required next step for journal-grade stall accuracy rather than claimed as solved; (iv) max non-orthogonality per mesh is now reported in the mesh metadata and the paper's mesh table.

## Summary of changes

| Reviewer point | Action taken |
|---|---|
| 1. Novelty | Reframed as V&V/reference study; realistic journal tier; novelty claims removed |
| 2. 2D devices | Demoted to negative-control appendix; no device claims from 2D |
| 3. Grid at stall | Fine grid primary at stall; medium-grid envelope labelled as such |
| 4. Steady stall | URANS low-Re campaign on the 50c farfield grid added; steady results restricted to attached range |
| 5. Drag ×2 | Headline MAE vs mismatched Re withdrawn; laminar-run drag budget explained |
| 6. Re mismatch | All quantitative metrics now Reynolds-matched (Sandia primary) |
| 7. Mesh/farfield | 50c farfield URANS grid; non-orthogonality reported; C-topology future work |

**Revised target:** Aerospace Science and Technology, Journal of Aircraft, or Journal of Applied Fluid Mechanics — applied journals where a rigorous, fully reproducible transition-model comparison with matched-Re validation is in scope.
