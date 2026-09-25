# Tissue–sinusoid–hepatocyte coupling

A multiscale framework linking virtual liver tissue, sinusoidal transport and
interchangeable hepatocyte SBML/CellML models.

This repository is the companion repository for the associated Interface
paper. It contains the code, generated geometry, manifests, figures and
interactive viewers required to reproduce the computational demonstrations
reported in the paper.

## Scientific aim

The project separates spatial geometry, vascular transport and local cellular
equations. A 3-D Voronoi scaffold supplies virtual lobular territories;
portal-site candidates and central-vein generators define sinusoidal paths;
local SBML/CellML models are executed at sinusoidal control volumes.

APAP metabolism and oxygen-sensitive oxidative stress are demonstrators of the
coupling architecture, not the limits of the framework.

## Repository layout

- `geometry/` — Voronoi tissue and OpenCMISS/FieldML exports
- `docs/` — coupling specification and manuscript material
- `viewers/` — interactive 3-D visualizations
- `figures/` — manuscript figure generation

## Status

Research prototype. Numerical outputs are demonstrators pending biological
calibration and structural validation.

## License and model provenance

The repository code and generated examples are released under the MIT License.
Any third-party SBML/CellML models, if used, remain subject to their original
licenses and provenance requirements.
