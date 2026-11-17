# Rebuilt coupling code — design choices and first findings (2026-11-17)

This is a clean re-implementation written from the Methods and the released geometry. It is NOT the
author's original code (which is on the author's home PC). Numbers from it will differ from the
submitted figures and the manuscript must be regenerated from it.

## Layout
- `hepcoupling/exchange.py` typed record I(i,k) -> O(i,k), declared units, closure residual
- `hepcoupling/network.py`  geometry loading, flow law, control-volume volume and residence time
- `hepcoupling/adapters.py` SBMLAdapter (libRoadRunner), CellMLAdapter (re-implementation)
- `hepcoupling/solver.py`   plug-flow marching, computed oxygen, lobule outlets, well-mixed baseline
- `smoke_test.py`           closure, flow balance, extraction, eps_h, well-mixed comparison

## Explicit assumptions (every one is an argument; all need author confirmation)
| Item | Default | Status |
|---|---|---|
| 1 geometry unit | 0.5 mm | ASSUMED so portal-central paths are ~0.3-0.6 mm |
| sinusoid radius | 5 um | typical, to cite |
| mean velocity | 300 um/s | typical, to cite |
| flow law | Poiseuille, common radius: q_i ∝ 1/L_i within a lobule | `uniform` also available |
| oxygen | COMPUTED, first-order k_O2 = -ln(0.75)/1.267 s (nominal 25%; actual mean path extraction ~20% because the full-scaffold mean residence time is 0.99 s) | demonstrator |
| inlet O2 | portal 6.5, artery 8.5 mol/m3, artery share 0.25 | to cite |
| inlet APAP | 0.07 mol/m3 = 70 uM | from the manuscript's 0.07 |
| CYP2E1 scaling | OFF by default; `cyp_o2_half` switches on f(O2) as an external override of Vmax_2E1_APAP | author to supply original form |
| hepatocyte:blood volume ratio phi | 1 (SBML compartment = control-volume blood) | rescales reaction time |
| transported species | all six SBML species advected with the flow | simplification, state in paper |
| flux J | integral of stoichiometry x reaction rates (Simpson), independent of x_out | closure is not by construction |
| CellML | re-implementation, O2 clamped, transported O2 mapped by ratio to a reference | unit placeholders unresolved |

## First findings (8 lobules, 32 paths, 16 control volumes each)
1. Mean path residence time 1.3 s; single-pass APAP extraction 0.19% (max 0.30%).
   Referee 1's estimate is confirmed: with the published SBML parameters at 70 uM the total
   elimination is ~1.1e-4 mM/s (CYP2E1 is ~0.9% of it), so depleting APAP to ~1% of inlet along a
   path needs ~10^3 s residence, not ~1 s. The submitted APAP/NAPQI gradients cannot come from
   physiological transit times.
2. GSH deviation from baseline is ~5e-7 mM here (the submitted 4e-4 mM excursion is not reproduced).
3. Closure: max |r| ~1e-29 mol/s against fluxes ~1e-19 mol/s (relative ~1e-10 for APAP).
   Flow balance per lobule is exact by construction (state this in the paper).
4. eps_h is ~1e-15 for APAP/GSH and ~1e-9 for NAPQI. This is NOT evidence of anything: the local
   model is autonomous and each control volume is advanced exactly, so subdividing a path cannot
   change the result (semigroup property). The refinement test only becomes informative when the
   per-control-volume input varies along the path (e.g. CYP2E1 scaling on, oxygen mid-point rule).
   The paper should say this rather than report ~0.
5. Well-mixed vs path-resolved: outlet APAP identical to 8 digits; NAPQI differs 3.6%. At
   physiological residence times the composition layer barely matters for APAP. This must be
   reported honestly in the response to Referee 1 point 6.
6. CellML re-implementation: H2O2 barely moves over seconds because the archived rates are 'per_day'
   under a time variable declared in hours (an inconsistency in the archive).

## Not yet done
- No unit tests beyond smoke_test.py; the 512-path / 16-cv / 128-lobule full run is untimed (~100 s est.)
- Multi-realisation (irregular generators) needs a Voronoi generator script (not in the released repo).
- Figures not regenerated.
