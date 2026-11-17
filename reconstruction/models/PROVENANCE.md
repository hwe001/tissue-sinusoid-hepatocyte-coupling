# Model provenance (downloaded 2026-11-17)

## 1. BIOMD0000000624 — SBML (Level 2 Version 4)
- Title: Sluka2016 - Acetaminophen metabolism (Sluka et al. 2016, PLoS ONE 11:e0162428, doi:10.1371/journal.pone.0162428)
- Source: https://www.ebi.ac.uk/biomodels/BIOMD0000000624
- Files: `biomodels_624/BIOMD0000000624_url.xml`, `biomodels_624/curation_notes.txt`
- SHA-256 (xml): 234178fd3d1f979f075bc7e31ba3c2180dc650c13de4ded25349fba53772b8ac
- Licence: CC0 (public domain dedication, stated in the SBML annotation)
- Facts checked with libRoadRunner 2.10.0:
  - species: APAP, NAPQI, GSH, NAPQIGSH, APAPconj_Glu, APAPconj_Sul
  - initial concentrations: APAP 0.1, GSH 10, others 0 (the model is a mM, 1 L compartment per the earlier draft; verify)
  - single compartment of volume 1; reactions J0-J4
  - NO oxygen species; CYP2E1 enters only via parameter `Vmax_2E1_APAP`

## 2. Schmitmeier 2007 — CellML 1.0
- Title: Improvement of metabolic performance of primary hepatocytes in hyperoxic cultures by vitamin C in a novel small-scale bioreactor (J. Membr. Sci. 298:30-40, doi:10.1016/j.memsci.2007.04.006)
- Source: https://models.physiomeproject.org/e/a2 (workspace hnielsen/schmitmeier_2007, revision b79c9b42d6844f0c29e292544e22c38f16f47f3e)
- File: `physiome_schmitmeier/schmitmeier_2007.cellml`
- SHA-256: 71e895a1c21269c127c47de82cb5dad3646a8999840b7b2fc25731ac212efde0
- Licence: Creative Commons Attribution 3.0 Unported
- Repository status note (must be disclosed in the paper): "This CellML model will run in OpenCell but does not reproduce published results. This is the 20% Oxygen version, two more versions to follow."
- Structure: components environment, O2, H2O2, V, model_parameters; two ODEs. Oxygen and H2O2 are STATE VARIABLES, so supplying oxygen from the transport layer changes model structure, not just a parameter (Referee 1, point 8).

## Environment
- Python venv: `revision/.venv` (libroadrunner 2.10.0, numpy 2.4.6, scipy 1.17.1, matplotlib)
- Downloaded files are untrusted data: they are only parsed, never executed.
