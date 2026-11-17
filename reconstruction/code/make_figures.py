"""Regenerate manuscript figures from the rebuilt code. Output: revision/figures/*.png (300 DPI) + *.pdf (vector).
Run with the venv python: python -I code/make_figures.py"""
import json, sys
from collections import defaultdict
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from mpl_toolkits.mplot3d.art3d import Line3DCollection

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from hepcoupling import (load_network, SBMLAdapter, CellMLAdapter, Inlet, run_network, lobule_outlets, k_o2_from_extraction)

REV = ROOT.parent
GEO = REV / "geometry"
NET = GEO / "lobule_sinusoid_network.json"
SBML = REV / "models" / "biomodels_624" / "BIOMD0000000624_url.xml"
RES = REV / "results"
FIG = REV / "figures"; FIG.mkdir(exist_ok=True)
SP = ["APAP", "NAPQI", "GSH", "NAPQIGSH", "APAPconj_Glu", "APAPconj_Sul"]
# Okabe-Ito colour-blind-safe palette
C = dict(blue="#0072B2", orange="#E69F00", green="#009E73", red="#D55E00", purple="#CC79A7", sky="#56B4E9", grey="#555555")
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False, "axes.titlesize": 9,
                     "axes.labelsize": 9, "legend.frameon": False, "figure.dpi": 100, "savefig.dpi": 300})


def save(fig, name):
    fig.savefig(FIG / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)
    print("saved", name)


def cell_edges(cell):
    """Polygon outline of a triangulated convex cell: merge coplanar triangles, keep boundary edges."""
    V = np.array(cell["vertices"]); groups = defaultdict(list)
    for tri in cell["faces"]:
        a, b, c = V[tri]; n = np.cross(b - a, c - a); n /= np.linalg.norm(n)
        if n @ (a - V.mean(0)) < 0: n = -n
        groups[(tuple(np.round(n, 3)), round(float(n @ a), 3))].append(tri)
    segs = []
    for tris in groups.values():
        cnt = defaultdict(int)
        for t in tris:
            for i in range(3):
                cnt[tuple(sorted((t[i], t[(i + 1) % 3])))] += 1
        segs += [(V[i], V[j]) for (i, j), k in cnt.items() if k == 1]
    return segs


inlet = Inlet()
paths = load_network(NET, n_cv=16)
net = json.loads(NET.read_text()); vor = json.loads((GEO / "true_3d_voronoi.json").read_text())
gen = np.array([c["position"] for c in net["central_veins"]]); portal = np.array([p["position"] for p in net["portal_sites"]])
K_O2 = k_o2_from_extraction(0.25, 1.267)
ad_on = SBMLAdapter(SBML, cyp_o2_half=1.0, o2_ref=inlet.o2)
cm = CellMLAdapter(z_ref=inlet.o2)
res = run_network(paths, inlet, ad_on, K_O2, SP, second_adapter=cm)
byid = {p.pid: p for p in paths}
xs = np.linspace(0, 1, 17)

# ---------------- Figure 2: computed scaffold ----------------
fig = plt.figure(figsize=(7.2, 3.6))
ax = fig.add_subplot(1, 2, 1, projection="3d")
segs = [s for c in vor["cells"] for s in cell_edges(c)]
ax.add_collection3d(Line3DCollection(segs, colors="#9a9a9a", linewidths=0.25, alpha=0.6))
ax.scatter(*gen.T, s=7, c=C["red"], depthshade=False, label="generators = central-vein candidates (128)")
ax.scatter(*portal.T, s=3, c=C["orange"], depthshade=False, label="shared vertices = portal candidates (432)")
ax.set_xlim(-2.4, 1.8); ax.set_ylim(-2.4, 1.8); ax.set_zlim(-2.4, 1.8)
ax.set_xlabel("x (geometry units)"); ax.set_ylabel("y"); ax.set_zlabel("z"); ax.view_init(22, -55)
ax.set_title("(a) 128 Voronoi cells", loc="left")
ax2 = fig.add_subplot(1, 2, 2, projection="3d")
i0 = 21; cell = vor["cells"][i0]
ax2.add_collection3d(Line3DCollection(cell_edges(cell), colors="#555555", linewidths=0.8))
pp = [p for p in net["sinusoid_paths"] if p["lobule_id"] == i0]
for p in pp:
    poly = np.array(p["polyline"]); ax2.plot(*poly.T, c=C["blue"], lw=1.4)
    ax2.scatter(*poly[0], s=22, c=C["orange"], depthshade=False)
ax2.scatter(*gen[i0], s=40, c=C["red"], depthshade=False)
vv = np.array(cell["vertices"]); lo, hi = vv.min(0), vv.max(0)
ax2.set_xlim(lo[0], hi[0]); ax2.set_ylim(lo[1], hi[1]); ax2.set_zlim(lo[2], hi[2]); ax2.view_init(22, -55)
ax2.set_title(f"(b) one territory (VLT {i0}), 4 paths", loc="left")
for a in (ax, ax2):
    a.xaxis.pane.fill = False; a.yaxis.pane.fill = False; a.zaxis.pane.fill = False
fig.legend(*ax.get_legend_handles_labels(), loc="lower center", ncol=2, fontsize=6.5, bbox_to_anchor=(0.5, -0.07), markerscale=1.5)
save(fig, "fig2_scaffold")

# ---------------- Figure 3: two local models on the same coordinate ----------------
L = np.array([p.length for p in paths]); ul = np.unique(np.round(L, 9))   # distinct lengths (lattice paths repeat)
pick = {lab: int(np.argmin(np.abs(L - ul[i]))) for lab, i in (("shortest", 0), ("middle", len(ul) // 2), ("longest", len(ul) - 1))}
fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.5))
for (lab, k), col in zip(pick.items(), (C["blue"], C["green"], C["red"])):
    r = res[k]; p = paths[k]
    axs[0].plot(xs, (1 - r.x[:, 0] / inlet.apap) * 100, c=col, label=f"{lab} ({p.length*1e3:.2f} mm, {p.tau_total:.2f} s)")
    axs[1].plot(xs, r.x[:, 1] * 1e6, c=col)
    h = [cm.H2O2_INIT] + [y["cellml_H2O2"] for y in r.y]
    axs[2].plot(xs, (np.array(h) / cm.H2O2_INIT - 1) * 1e6, c=col)
axs[0].set_ylabel("APAP removed (% of inlet)"); axs[1].set_ylabel("NAPQI (nmol L$^{-1}$)")
axs[2].set_ylabel("H$_2$O$_2$ change (ppm)")
axs[0].set_title("(a) SBML: APAP", loc="left"); axs[1].set_title("(b) SBML: NAPQI", loc="left")
axs[2].set_title("(c) CellML re-impl.: H$_2$O$_2$", loc="left")
for a in axs: a.set_xlabel("normalised position")
fig.tight_layout()
fig.legend(*axs[0].get_legend_handles_labels(), loc="lower center", ncol=3, fontsize=6.5, bbox_to_anchor=(0.5, -0.1))
save(fig, "fig3_two_models")

# ---------------- Figure 4: per-territory outlet fields ----------------
# The regular BCC scaffold gives identical interior territories, so Figure 4 uses a seeded irregular scaffold
# (jitter sigma = 0.10, seed 0) where territories genuinely differ in volume and path length.
from hepcoupling.geometry import build
net_j, gst_j = build(0, "jitter", 0.10, unit_cv=8)
tmpj = RES / "tmp_fig4.json"; tmpj.write_text(json.dumps(net_j))
paths_j = load_network(tmpj, n_cv=16); tmpj.unlink()
res_j = run_network(paths_j, inlet, ad_on, K_O2, SP)
lo_ = lobule_outlets(res_j); ids = sorted(lo_)
gen4 = np.array([c["position"] for c in net_j["central_veins"]])
Y = np.array([lo_[i] for i in ids]); Zout = np.array([np.mean([r.z[-1] for r in res_j if r.lobule_id == i]) for i in ids])
fields = [("outlet O$_2$ (mol m$^{-3}$)", Zout), ("APAP removed (%)", (1 - Y[:, 0] / inlet.apap) * 100),
          ("outlet NAPQI (nmol L$^{-1}$)", Y[:, 1] * 1e6)]
fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.5))
ms = (gen4[:, 2] > -0.75) & (gen4[:, 2] < 0.25)                          # slab of territories shown (x-y projection)
for k, (lab, v) in enumerate(fields):
    a = axs[k]; vmin, vmax = v[ms].min(), v[ms].max()
    sc = a.scatter(gen4[ms, 0], gen4[ms, 1], c=v[ms], s=45, cmap="viridis", vmin=vmin, vmax=vmax, edgecolor="k", lw=0.3)
    a.set_aspect("equal"); a.set_xlabel("x"); a.set_ylabel("y")
    cb = fig.colorbar(sc, ax=a, shrink=0.8, pad=0.03); cb.ax.tick_params(labelsize=7)
    a.set_title(f"({'abc'[k]}) {lab}", fontsize=7, loc="left")
    print(f"fig4 {lab}: all 128 range {v.min():.4g}-{v.max():.4g}; slab (n={ms.sum()}) range {vmin:.4g}-{vmax:.4g}")
fig.tight_layout(); save(fig, "fig4_territory_fields")

# ---------------- Figure 5: dual inlet (hepatic-artery share) ----------------
fig, axs = plt.subplots(1, 4, figsize=(7.2, 2.3))
for fh, col in zip((0.0, 0.25, 0.5), (C["blue"], C["green"], C["red"])):
    inl = Inlet(frac_ha=fh)
    ad = SBMLAdapter(SBML, cyp_o2_half=1.0, o2_ref=Inlet().o2)          # same scaling reference for all
    rr = run_network(paths, inl, ad, K_O2, SP)
    mz = np.mean([r.z for r in rr], axis=0); mx = np.mean([r.x for r in rr], axis=0)
    lab = f"hepatic artery {fh*100:.0f}% of inflow"
    axs[0].plot(xs, mz, c=col, label=lab); axs[1].plot(xs, (1 - mx[:, 0] / inl.apap) * 100, c=col)
    axs[2].plot(xs, mx[:, 1] * 1e6, c=col); axs[3].plot(xs, (mx[:, 2] - inl.gsh) * 1e6, c=col)
for a, t in zip(axs, ("O$_2$ (mol m$^{-3}$)", "APAP removed (%)", "NAPQI (nmol L$^{-1}$)", "GSH $-$ baseline (nmol L$^{-1}$)")):
    a.set_ylabel(t); a.set_xlabel("normalised position")
for k, a in enumerate(axs): a.set_title("(" + "abcd"[k] + ")", loc="left")
fig.tight_layout()
fig.legend(*axs[0].get_legend_handles_labels(), loc="lower center", ncol=3, fontsize=6.5, bbox_to_anchor=(0.5, -0.1))
save(fig, "fig5_dual_inlet")

# ---------------- Figure 6: per-path traces and mean ----------------
fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.4))
for r in res:
    axs[0].plot(xs, r.z, c=C["sky"], lw=0.3, alpha=0.5); axs[1].plot(xs, (1 - r.x[:, 0] / inlet.apap) * 100, c=C["sky"], lw=0.3, alpha=0.5)
    axs[2].plot(xs, r.x[:, 1] * 1e6, c=C["sky"], lw=0.3, alpha=0.5)
axs[0].plot(xs, np.mean([r.z for r in res], 0), c=C["blue"], lw=1.6)
axs[1].plot(xs, np.mean([(1 - r.x[:, 0] / inlet.apap) * 100 for r in res], 0), c=C["blue"], lw=1.6)
axs[2].plot(xs, np.mean([r.x[:, 1] * 1e6 for r in res], 0), c=C["blue"], lw=1.6)
for a, t in zip(axs, ("O$_2$ (mol m$^{-3}$)", "APAP removed (%)", "NAPQI (nmol L$^{-1}$)")):
    a.set_ylabel(t); a.set_xlabel("normalised position")
for k, a in enumerate(axs): a.set_title("(" + "abc"[k] + f") {len(res)} paths, mean in dark", loc="left", fontsize=7.5)
fig.tight_layout(); save(fig, "fig6_all_paths")

# ---------------- Figure 7 (new): verification and baseline ----------------
d_off = json.loads((RES / "full_off.json").read_text()); d_on = json.loads((RES / "full_on.json").read_text())
rows = json.loads((RES / "seed_sweep.json").read_text())
fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.5))
for d, col, lab in ((d_on, C["blue"], "CYP2E1 oxygen scaling on"), (d_off, C["grey"], "off")):
    h = [r for r in d["runs"].values() if "eps_h_vs_previous" in r]
    xx = [r["eps_h_vs_previous"]["coarse_n_cv"] for r in h]
    axs[0].loglog(xx, [r["eps_h_vs_previous"]["max"]["NAPQI"] for r in h], "o-", c=col, ms=4,
                  label=f"{lab}" if lab != "off" else "scaling off (round-off level)")
axs[0].loglog([8, 16, 32], [1.4e-4, 3.5e-5, 8.75e-6], ":", c="k", lw=0.8); axs[0].text(10, 2.2e-5, r"$\propto h^{2}$", fontsize=7)
axs[0].set_xlabel("coarse control volumes per path"); axs[0].set_ylabel(r"max $\varepsilon_h$ over lobules")
axs[0].set_xticks([8, 16, 32]); axs[0].set_xticklabels(["8", "16", "32"]); axs[0].minorticks_off()
axs[0].legend(fontsize=6, loc="upper center", bbox_to_anchor=(0.5, -0.28))
axs[0].set_title("(a) refinement (NAPQI)", loc="left")
r64 = d_on["runs"]["16"]; r64o = d_off["runs"]["16"]
bars = [r64o["well_mixed_rel_diff_per_lobule_max"]["NAPQI"] * 100, r64["well_mixed_rel_diff_per_lobule_max"]["NAPQI"] * 100,
        r64o["well_mixed_rel_diff_of_mean"]["NAPQI"] * 100, r64["well_mixed_rel_diff_of_mean"]["NAPQI"] * 100]
axs[1].bar(range(4), bars, color=[C["grey"], C["blue"], C["grey"], C["blue"]])
axs[1].set_xticks(range(4)); axs[1].set_xticklabels(["max\nlobule\n(off)", "max\nlobule\n(on)", "mean\n(off)", "mean\n(on)"], fontsize=7)
axs[1].set_ylabel("NAPQI difference (%)"); axs[1].set_title("(b) well-mixed vs resolved", loc="left")
geoms = [("jitter", 0.05, "jitter 0.05"), ("jitter", 0.10, "jitter 0.10"), ("random", 0.0, "random")]
for k, (m, s, lab) in enumerate(geoms):
    R = [r for r in rows if r["mode"] == m and r["sigma"] == s and r["cyp"] == 1]
    axs[2].scatter(np.full(len(R), k) + np.linspace(-.12, .12, len(R)), [r["outlet_NAPQI_mean"] * 1e6 for r in R], s=10, c=C["blue"])
reg = [r for r in rows if r["sigma"] == 0.0 and r["mode"] == "jitter" and r["cyp"] == 1][0]
axs[2].axhline(reg["outlet_NAPQI_mean"] * 1e6, c=C["grey"], ls="--", lw=0.8); axs[2].text(-0.45, reg["outlet_NAPQI_mean"] * 1e6 + 0.0012, "regular BCC", fontsize=6.5, ha="left", color=C["grey"])
axs[2].set_xticks(range(3)); axs[2].set_xticklabels([g[2] for g in geoms], fontsize=7)
axs[2].set_ylabel("mean outlet NAPQI (nmol L$^{-1}$)"); axs[2].set_title("(c) 10 seeds per geometry", loc="left")
fig.tight_layout(); save(fig, "fig7_verification")
print("regular reference", reg["outlet_NAPQI_mean"] * 1e6)
