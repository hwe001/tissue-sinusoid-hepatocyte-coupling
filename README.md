# Tissue–sinusoid–hepatocyte coupling

A multiscale framework linking virtual liver tissue, sinusoidal transport and
interchangeable hepatocyte SBML/CellML models.

## Scientific aim

The project separates spatial geometry, vascular transport and local cellular
equations. A 3-D Voronoi scaffold supplies virtual lobular territories;
portal-site candidates and central-vein generators define sinusoidal paths;
local SBML/CellML models are executed at sinusoidal control volumes.

APAP metabolism and oxygen-sensitive oxidative stress are demonstrators of the
coupling architecture, not the limits of the framework.

## Repository layout

- `geometry/` — Voronoi tissue and OpenCMISS/FieldML exports
- `models/` — source SBML/CellML models and provenance
- `coupling/` — spatial mapping and transport code
- `viewers/` — interactive 3-D visualizations
- `figures/` — manuscript figure generation
- `docs/` — coupling specification and manuscript material

## Status

Research prototype. Numerical outputs are demonstrators pending biological
calibration and structural validation.

## License and model provenance

The framework code and generated examples will be released under an open-source
license. Third-party SBML/CellML models remain subject to their original
licenses and should be downloaded from their authoritative repositories.
