# Cross-scale coupling specification

## Purpose

This specification defines how an executable local SBML/CellML model is placed
inside a virtual hepatic tissue scaffold. The local model is replaceable; the
spatial interface is not tied to APAP or to one biochemical formalism.

## Spatial hierarchy

| Level | Entity | Primary state |
|---|---|---|
| Tissue | Voronoi scaffold | lobule topology, vascular boundaries |
| Lobule | Voronoi cell | central-vein position, portal-site set |
| Sinusoid | directed path/graph | flow, pressure, residence time |
| Control volume | axial sinusoid element | local transported concentrations |
| Cell model | SBML/CellML instance | intracellular states and fluxes |

## Required interface variables

Each local model instance receives:

| Variable | Meaning | Interface unit |
|---|---|---|
| `oxygen` | local oxygen concentration or normalized activity | explicitly declared by model |
| `flow` | local sinusoidal flow | volume/time |
| `residence_time` | time spent in control volume | time |
| `substrates` | transported extracellular species | model-declared concentration |
| `cell_volume` | represented hepatocyte volume | volume |

Each local model returns:

| Output | Meaning |
|---|---|
| `consumption_flux` | uptake or consumption of transported species |
| `production_flux` | release of metabolites or stress products |
| `functional_state` | optional scalar/vector such as Schmitmeier `V` |
| `intracellular_state` | selected states for analysis or feedback |

## Coupling sequence

For each time step and each directed sinusoidal edge:

1. Obtain flow and residence time from the vascular graph.
2. Advect extracellular concentrations into the control volume.
3. Pass local concentrations and oxygen to the cell model.
4. Integrate the local model for the residence-time interval.
5. Convert returned rates to source/sink terms.
6. Update the outgoing sinusoidal state.
7. Accumulate outputs at the central vein and tissue level.

The local model must not directly access global geometry. Geometry and
physiology communicate only through the declared interface.

## Mass-balance checks

For every transported species `q`:

`inflow(q) - outflow(q) + production(q) - consumption(q) = storage_change(q)`

The implementation should report the residual per control volume and reject a
run if the residual exceeds a user-defined tolerance.

## Spatial mapping

- A Voronoi cell is a lobular territory, not a Couinaud segment.
- Its generator is a central-vein candidate.
- Shared interior vertices are candidate portal sites.
- A portal site may serve several adjacent lobular territories.
- Sinusoid paths are directed from portal sites toward the central vein.
- Zonation is derived from local transport states, not assigned as a label.

## Demonstrator models

- `BIOMD0000000624.xml`: APAP xenobiotic metabolism, SBML.
- `schmitmeier_2007.cellml`: oxygen/oxidative-stress response, CellML.

The demonstrators use the same spatial interface and differ only in their local
equation layer. This is the central interoperability claim of the framework.

## Verification requirements

Before publication, every model adapter must report:

- source model identifier and version;
- variable mapping and units;
- initial conditions;
- solver and tolerances;
- conservation residuals;
- spatial and temporal resolution;
- whether execution is direct or equation-faithful translation.
