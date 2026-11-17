# Step 3 results (rebuilt code, assumed parameters; see code/DESIGN.md)

Inlet APAP 0.07 mol/m3, mean path velocity 300 um/s, radius 5 um, 1 geometry unit = 0.5 mm,
oxygen computed (25% single-pass extraction), phi = 1. CYP scaling "on" uses half-saturation
K = 1 mol/m3 (ASSUMED). All numbers regenerate with `code/run_full.py` and `code/run_seeds.py`.

## A. Released regular scaffold, 128 lobules, 512 paths
| n_cv per path | sites |
|---|---|
| 8 / 16 / 32 / 64 | 4,096 / 8,192 / 16,384 / 32,768 |

Single-pass APAP extraction: mean 0.155%, max 0.305% (identical at every refinement).
Closure: max |r| / max |J| = 7e-11 (APAP) to 1.5e-6 (GSH), set by the quadrature, not zero by construction.
Lobule flow balance: exact by construction.

### eps_h = |y_h - y_h/2| / |y_h/2| (max over 128 lobule outlets)
| refinement | CYP off: APAP / NAPQI / GSH | CYP on: APAP / NAPQI / GSH |
|---|---|---|
| 8->16  | 6e-15 / 4e-09 / 2e-15 | 2.6e-10 / 1.4e-04 / 9.8e-12 |
| 16->32 | 1e-14 / 5e-09 / 1e-15 | 6.5e-11 / 3.5e-05 / 2.5e-12 |
| 32->64 | 3e-14 / 2e-08 / 5e-15 | 1.6e-11 / 8.7e-06 / 6.2e-13 |

CYP off: values at round-off, exact because the local model is autonomous (no information in the test).
CYP on: ratios of ~4 per halving of h, i.e. second-order convergence of the oxygen mid-point rule.

## B. Well-mixed single compartment (same SBML, same total residence time and flow)
Regular scaffold, CYP off: outlet APAP differs by 3e-8 (relative, of the mean); NAPQI 0.7% of the
mean (max per lobule 7%). CYP on: NAPQI 2.5% of the mean (max per lobule 10%).

## C. Seeded irregular scaffolds (10 seeds each, n_cv = 16, regenerated geometry)
| geometry | territory volume SD | path length mm (mean±SD, in units x 0.5) | outlet APAP (mean±SD over seeds) | outlet NAPQI mean (±SD over seeds) | lobule-to-lobule CV of NAPQI | well-mixed vs resolved NAPQI (mean / max) CYP off |
|---|---|---|---|---|---|---|
| regular BCC (1 realisation) | 0.036 | 0.565±0.026 | 0.069895 | 6.26e-07 | 0.012 | 0.1% / 2% |
| jitter 0.05 | 0.047 | 0.557±0.045 | 0.069896±3e-07 | 6.20e-07±1e-09 | 0.032 | 0.3% / 6% |
| jitter 0.10 | 0.070 | 0.557±0.053 | 0.069896±5e-07 | 6.19e-07±2e-09 | 0.044 | 0.3% / 8% |
| random 128 | 0.150 | 0.582±0.077 | 0.069892±1e-06 | 6.33e-07±5e-09 | 0.061 | 0.7% / 10% |

(Path lengths are in geometry units; multiply by 0.5 mm. Territory volume is in units^3, mean 0.5.)

## Caveats
- The regenerated "regular BCC" uses my portal-selection rule, so it is not identical to the released
  network (full released scaffold mean residence time 0.99 s, regenerated 0.94 s).
- The oxygen rate constant was fixed from a 1.27 s reference (taken from an 8-lobule test), so the actual
  mean oxygen extraction over a path is about 20% (7.0 -> about 5.6 mol/m3), not the nominal 25%.
- Boundary cells are clipped by mirroring at the box planes, which makes their volumes smaller (min 0.35).
- With physiological transit times nothing in the outlets depends strongly on geometry: APAP outlets
  agree to ~1e-5 across seeds. Variability appears mainly in NAPQI (lobule-to-lobule CV 1% to 6%).
- No part of the submitted figures has been regenerated.
