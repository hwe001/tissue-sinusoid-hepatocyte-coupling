"""Smoke test: closure, flow balance, physiological sanity, refinement. Run with .venv python -I."""
import sys, time
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from hepcoupling import (load_network, SBMLAdapter, CellMLAdapter, Inlet, run_network,
                         lobule_outlets, well_mixed_reference, k_o2_from_extraction)

REV = ROOT.parent
NET = REV / "geometry" / "lobule_sinusoid_network.json"
SBML = REV / "models" / "biomodels_624" / "BIOMD0000000624_url.xml"
SP = ["APAP", "NAPQI", "GSH", "NAPQIGSH", "APAPconj_Glu", "APAPconj_Sul"]
inlet = Inlet()

def run(n_cv, law="poiseuille", subset=None, phi=1.0):
    paths = load_network(NET, n_cv=n_cv, flow_law=law)
    if subset: paths = [p for p in paths if p.lobule_id < subset]
    tau_ref = float(np.mean([p.tau_total for p in paths]))
    k = k_o2_from_extraction(0.25, tau_ref)
    ad = SBMLAdapter(SBML, phi=phi)
    t0 = time.time()
    res = run_network(paths, inlet, ad, k, SP)
    return paths, res, ad, time.time() - t0, tau_ref

paths, res, ad, dt, tau_ref = run(16, subset=8)
q = np.array([p.q for p in paths]); L = np.array([p.length for p in paths]); tau = np.array([p.tau_total for p in paths])
print(f"[geometry] paths={len(paths)} length mm: mean={L.mean()*1e3:.3f} range=({L.min()*1e3:.3f},{L.max()*1e3:.3f})")
print(f"[flow] q per path m3/s mean={q.mean():.3e}; residence tau s: mean={tau.mean():.3f} range=({tau.min():.3f},{tau.max():.3f})")
# flow balance per lobule
from collections import defaultdict
fb = defaultdict(float)
for p in paths: fb[p.lobule_id] += p.q
Qexp = 4 * np.pi * (5e-6) ** 2 * 3e-4
print(f"[flow balance] max |sum q - Q|/Q = {max(abs(v - Qexp) / Qexp for v in fb.values()):.2e}")
# closure
R = np.concatenate([r.resid.reshape(-1, len(SP)) for r in res])
J = np.concatenate([r.J.reshape(-1, len(SP)) for r in res])
for j, s in enumerate(SP):
    scale = np.max(np.abs(J[:, j])) or 1.0
    print(f"[closure] {s:13s} max|r|={np.max(np.abs(R[:, j])):.2e} mol/s  (max|J|={np.max(np.abs(J[:, j])):.2e}; ratio {np.max(np.abs(R[:, j]))/scale:.1e})")
# physiology
out = np.array(list(lobule_outlets(res).values()))
extr = 1 - out[:, 0] / inlet.apap
print(f"[extraction] single-pass APAP extraction: mean={extr.mean():.3e} max={extr.max():.3e}  (referee estimated ~unity needed for paper's gradients)")
print(f"[O2] inlet={inlet.o2:.2f} mol/m3; outlet mean={np.mean([r.z[-1] for r in res]):.2f}")
print(f"[GSH] outlet-baseline = {np.mean(out[:,2]) - 10:.3e} mM")
wm = np.array(list(well_mixed_reference(paths, inlet, ad, SP).values()))
print(f"[well-mixed vs path-resolved] outlet APAP {wm[:,0].mean():.8f} vs {out[:,0].mean():.8f}; NAPQI {wm[:,1].mean():.3e} vs {out[:,1].mean():.3e}")
print(f"[cost] {len(paths)*16} control volumes in {dt:.1f} s")

# refinement (eps_h) on lobule-mean outlet
prev = None
for n in (8, 16, 32, 64):
    _, r, *_ = run(n, subset=8)
    o = np.array(list(lobule_outlets(r).values())).mean(axis=0)
    if prev is not None:
        eps = np.abs(prev - o) / np.maximum(np.abs(o), 1e-300)
        print(f"[eps_h] {pn}->{n}: " + " ".join(f"{s}={e:.2e}" for s, e in zip(SP[:3], eps[:3])))
    prev, pn = o, n

cm = CellMLAdapter(z_ref=inlet.o2)
from hepcoupling import ExchangeIn
rec = ExchangeIn(q_in=paths[0].q, x_in=inlet.species(), z=inlet.o2, V=paths[0].V_cv, tau=paths[0].tau_cv)
print("[cellml re-impl] ", cm.step(rec).y)
