"""Full 128-lobule run with refinement study. Usage: python -I run_full.py {off|on} [n_cv ...]"""
import json, sys, time
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from hepcoupling import (load_network, SBMLAdapter, Inlet, run_network, lobule_outlets,
                         well_mixed_reference, k_o2_from_extraction)

REV = ROOT.parent
NET = REV / "geometry" / "lobule_sinusoid_network.json"
SBML = REV / "models" / "biomodels_624" / "BIOMD0000000624_url.xml"
OUT = REV / "results"; OUT.mkdir(exist_ok=True)
SP = ["APAP", "NAPQI", "GSH", "NAPQIGSH", "APAPconj_Glu", "APAPconj_Sul"]

mode = sys.argv[1]
ncvs = [int(a) for a in sys.argv[2:]] or [8, 16, 32, 64]
inlet = Inlet()
cyp = dict(cyp_o2_half=1.0, o2_ref=inlet.o2) if mode == "on" else {}   # K=1 mol/m3 ASSUMED
out = {"mode": mode, "cyp_o2_half": cyp.get("cyp_o2_half"), "inlet": inlet.__dict__ | {"o2_mixed": inlet.o2}, "runs": {}}
prev = None
for n in ncvs:
    paths = load_network(NET, n_cv=n)
    tau_ref = float(np.mean([p.tau_total for p in paths]))
    k = k_o2_from_extraction(0.25, 1.267)   # fixed k_O2 so refinement does not change physics
    ad = SBMLAdapter(SBML, **cyp)
    t0 = time.time()
    res = run_network(paths, inlet, ad, k, SP)
    lo = lobule_outlets(res)
    ids = sorted(lo)
    Y = np.array([lo[i] for i in ids])                                   # (128, 6)
    Z = np.array([np.mean([r.z[-1] for r in res if r.lobule_id == i]) for i in ids])
    R = np.concatenate([r.resid.reshape(-1, len(SP)) for r in res]); J = np.concatenate([r.J.reshape(-1, len(SP)) for r in res])
    entry = {
        "n_cv": n, "n_sites": len(paths) * n, "seconds": time.time() - t0,
        "tau_total_mean_s": tau_ref,
        "outlet_mean": dict(zip(SP, Y.mean(0))), "outlet_min": dict(zip(SP, Y.min(0))), "outlet_max": dict(zip(SP, Y.max(0))),
        "outlet_O2_mean": float(Z.mean()),
        "apap_extraction_mean": float((1 - Y[:, 0] / inlet.apap).mean()), "apap_extraction_max": float((1 - Y[:, 0] / inlet.apap).max()),
        "closure_max_abs_resid": dict(zip(SP, np.abs(R).max(0))), "closure_max_abs_J": dict(zip(SP, np.abs(J).max(0))),
        "total_mass_residual_rel": {s: float(abs(R[:, j].sum()) / max(abs(J[:, j]).sum(), 1e-300)) for j, s in enumerate(SP)},
    }
    wm = np.array([v for _, v in sorted(well_mixed_reference(paths, inlet, ad, SP).items())])
    entry["well_mixed_outlet_mean"] = dict(zip(SP, wm.mean(0)))
    entry["well_mixed_rel_diff_per_lobule_max"] = dict(zip(SP, (np.abs(wm - Y) / np.maximum(np.abs(Y), 1e-300)).max(0)))
    entry["well_mixed_rel_diff_of_mean"] = dict(zip(SP, np.abs(wm.mean(0) - Y.mean(0)) / np.maximum(np.abs(Y.mean(0)), 1e-300)))
    if prev is not None:   # eps_h = |y_h - y_{h/2}| / |y_{h/2}|, coarse h = prev, fine = current
        e = np.abs(prev - Y) / np.maximum(np.abs(Y), 1e-300)
        entry["eps_h_vs_previous"] = {"coarse_n_cv": out["runs"][str(prev_n)]["n_cv"],
                                       "max": dict(zip(SP, e.max(0))), "mean": dict(zip(SP, e.mean(0)))}
    np.save(OUT / f"outlets_{mode}_{n}.npy", Y)
    out["runs"][str(n)] = entry
    prev, prev_n = Y, n
    (OUT / f"full_{mode}.json").write_text(json.dumps(out, indent=1, default=float))
    print(f"[{mode}] n_cv={n} done in {entry['seconds']:.0f} s", flush=True)
