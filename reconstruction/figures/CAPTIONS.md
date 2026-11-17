# Draft captions (regenerated from the rebuilt code; every value depends on the assumptions in code/DESIGN.md)

Files: `fig2_scaffold`, `fig3_two_models`, `fig4_territory_fields`, `fig5_dual_inlet`, `fig6_all_paths`,
`fig7_verification` (new), each as 300 DPI PNG and vector PDF. Colours are the Okabe-Ito colour-blind-safe set.

**Figure 2. Computed virtual tissue scaffold.** (a) The 128 Voronoi territories read from
`geometry/true_3d_voronoi.json`, drawn as cell edges; red points are the 128 cell generators
(central-vein candidates) and gold points the 432 shared interior vertices (portal-site candidates).
The generators form a regular body-centred-cubic lattice (spacing 1 geometry unit); nothing is
randomised, so every interior territory is a truncated octahedron. (b) One territory (VLT 21) with its
four straight portal-to-central-vein paths. These are computational correspondences, not reconstructed
anatomical structures.

**Figure 3. Two local models on the same paths.** The shortest, a middle and the longest of the distinct
path lengths in the scaffold (legend gives length and residence time). (a, b) BIOMD0000000624 executed with
libRoadRunner: APAP removed as a percentage of the 70 umol/L inlet, and NAPQI in nmol/L. (c) Python
re-implementation of the archived Schmitmeier CellML equations (not an execution of the CellML file;
oxygen is supplied by the transport layer and its ODE removed; the archive's units carry unresolved
CONFIRM placeholders), shown as change of H2O2 from the archived initial value in ppm.

**Figure 4. Per-territory outlet fields on an irregular scaffold.** Seeded irregular scaffold (generators
displaced from the BCC lattice by N(0, 0.1^2) per axis, seed 0), because interior territories of the
regular lattice are identical. Each point is one territory in the slab -0.75 < z < 0.25 (32 of 128),
plotted at its generator in x-y projection; colour is the flow-weighted outlet value over the paths of that
territory: (a) oxygen, (b) APAP removed, (c) NAPQI. Colour limits are those of the slab; ranges over all
128 territories are O2 5.42-5.84 mol/m3, APAP removed 0.127-0.179 %, NAPQI 0.556-0.678 nmol/L. The
differences arise only from path length and flow split, not from assigned boundary conditions.

**Figure 5. Effect of hepatic-artery inflow share.** Mean over the 512 paths of the released scaffold with
the hepatic artery supplying 0, 25 and 50 % of inflow. (a) Computed oxygen (first-order consumption set for
25 % single-pass extraction; inlet 6.5 mol/m3 portal vein and 8.5 mol/m3 hepatic artery). (b-d) APAP removed,
NAPQI and GSH deviation from its 10 mmol/L baseline. The metabolite profiles are almost independent of the
oxygen supply because CYP2E1 accounts for only about 0.9 % of APAP elimination in the source model.

**Figure 6. All 512 paths.** Oxygen, APAP removed and NAPQI along every path (light) and their mean (dark).
The regular scaffold has only a few distinct path lengths, so individual traces fall into a small number of
bundles. All paths share the same inlet conditions; path length (0.28-0.57 mm) and flow split differ.

**Figure 7 (new). Verification and baselines.** (a) Maximum refinement statistic eps_h over the 128
territory outlets for NAPQI between 8->16, 16->32 and 32->64 control volumes per path. With the
oxygen-dependent CYP2E1 override on, eps_h falls by a factor of about 4 per halving (second order); with it
off, eps_h is at round-off because the local model is autonomous and each control volume is advanced exactly.
(b) NAPQI difference between the same SBML model run as one well-stirred compartment per territory (same total
flow and residence time) and the path-resolved calculation, for the worst territory and for the territory
mean. (c) Mean outlet NAPQI over 128 territories for ten seeds of each irregular-generator scheme (jitter
sigma 0.05 and 0.10, and random generators with minimum separation); dashed line is the regular lattice.
