# Hepatic model-composition framework: reconstructed reference implementation

Code, models, results and figures accompanying the revision of
"A framework for model composition across virtual hepatic tissue scales" (J. R. Soc. Interface, rsif-2026-1108).

**Status of this code (please read).** This is a clean re-implementation written from the
published Methods. It is NOT the code that produced the originally submitted results. Parameters not
fixed by the paper (geometry scale, vessel radius, mean velocity, inlet oxygen, oxygen consumption rate,
CYP2E1 oxygen dependence) are explicit assumptions, listed with their status in `code/DESIGN.md`.
<!-- AUTHOR: add one sentence on the relation to the original code once it has been compared (after 7 Dec 2026). -->

## Layout
- `code/hepcoupling/` library: `exchange.py` (typed contract and closure residual), `network.py`, `solver.py`, `adapters.py` (SBML and CellML), `geometry.py` (seeded Voronoi generator)
- `code/run_full.py`, `run_seeds.py`, `make_figures.py`, `smoke_test.py`
- `geometry/` released scaffold (BCC lattice, 128 lobules, 512 sinusoid paths)
- `models/` BIOMD0000000624 (CC0) and the Physiome Schmitmeier 2007 CellML (CC BY 3.0); see `models/PROVENANCE.md`
- `results/` saved outputs and `RESULTS_STEP3.md`; `figures/` Figures 2-7 (PNG 300 DPI, PDF) with `CAPTIONS.md`

## Reproduce
```bash
python -m venv .venv && .venv/Scripts/activate      # Linux/macOS: source .venv/bin/activate
pip install -r requirements.txt
python -I code/smoke_test.py            # about 10 s
python -I code/run_full.py off 8 16 32 64   # CYP2E1 off
python -I code/run_full.py on  8 16 32 64   # CYP2E1 on
python -I code/run_seeds.py 10
python -I code/make_figures.py
```
Tested with Python 3.x on Windows 11 and the package versions in `requirements.txt`.

## Known limitations
- The CellML adapter is a re-implementation; the source model's units are placeholders and the Physiome repository states that it does not reproduce published results.
- At about 1 s residence time the acetaminophen model removes about 0.16% of inlet per pass; results are verification of the coupling, not predictions.
- Oxygen enters the SBML model only through a scaling of `Vmax_2E1_APAP` (the model has no oxygen species).

## Licence and citation
Code (`code/`): MIT, see `LICENSE`. Results, figures, geometry and documentation: CC BY 4.0, see `LICENSE-DATA.md`.
<!-- AUTHOR: add Zenodo DOI. -->
Third-party models keep their own licences (see `models/PROVENANCE.md`).
