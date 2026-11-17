"""Multi-realisation sweep over seeded irregular Voronoi scaffolds. Usage: python -I run_seeds.py [n_seeds]"""
import json, os, sys, tempfile
from pathlib import Path
from multiprocessing import Pool
import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
REV = ROOT.parent
SBML = REV / "models" / "biomodels_624" / "BIOMD0000000624_url.xml"
OUT = REV / "results"
SP = ["APAP", "NAPQI", "GSH", "NAPQIGSH", "APAPconj_Glu", "APAPconj_Sul"]
N_CV = 16


def one(job):
    seed, mode, sigma, cyp = job
    from hepcoupling import load_network, SBMLAdapter, Inlet, run_network, lobule_outlets, well_mixed_reference, k_o2_from_extraction
    from hepcoupling.geometry import build
    net, gst = build(seed, mode, sigma, unit_cv=8)
    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False, dir=OUT) as f:
        json.dump(net, f); tmp = f.name
    try:
        paths = load_network(tmp, n_cv=N_CV)
    finally:
        os.remove(tmp)
    inlet = Inlet()
    ad = SBMLAdapter(SBML, **(dict(cyp_o2_half=1.0, o2_ref=inlet.o2) if cyp else {}))
    k = k_o2_from_extraction(0.25, 1.267)
    res = run_network(paths, inlet, ad, k, SP)
    lo = lobule_outlets(res); ids = sorted(lo)
    Y = np.array([lo[i] for i in ids])
    wm = np.array([v for _, v in sorted(well_mixed_reference(paths, inlet, ad, SP).items())])
    tau = np.array([p.tau_total for p in paths]); q = np.array([p.q for p in paths])
    R = np.concatenate([r.resid.reshape(-1, len(SP)) for r in res]); J = np.concatenate([r.J.reshape(-1, len(SP)) for r in res])
    return {"seed": seed, "mode": mode, "sigma": sigma, "cyp": cyp, "geometry": gst,
            "tau_mean": tau.mean(), "tau_sd": tau.std(), "tau_min": tau.min(), "tau_max": tau.max(),
            "outlet_APAP_mean": Y[:, 0].mean(), "outlet_APAP_sd_over_lobules": Y[:, 0].std(),
            "outlet_NAPQI_mean": Y[:, 1].mean(), "outlet_NAPQI_sd_over_lobules": Y[:, 1].std(),
            "outlet_NAPQI_cv_over_lobules": Y[:, 1].std() / Y[:, 1].mean(),
            "apap_extraction_mean": float((1 - Y[:, 0] / inlet.apap).mean()),
            "gsh_dev_mean": float((Y[:, 2] - inlet.gsh).mean()),
            "wm_rel_diff_NAPQI_mean": float(np.abs(wm[:, 1] - Y[:, 1]).mean() / Y[:, 1].mean()),
            "wm_rel_diff_NAPQI_max": float((np.abs(wm[:, 1] - Y[:, 1]) / Y[:, 1]).max()),
            "closure_rel_APAP": float(np.abs(R[:, 0]).max() / np.abs(J[:, 0]).max())}


if __name__ == "__main__":
    n_seeds = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    jobs = [(0, "jitter", 0.0, c) for c in (0, 1)]                       # regular BCC reference
    for mode, sig in (("jitter", 0.05), ("jitter", 0.10), ("random", 0.0)):
        for s in range(n_seeds):
            for c in (0, 1):
                jobs.append((s, mode, sig, c))
    with Pool(max(1, (os.cpu_count() or 2) - 2)) as pool:
        rows = pool.map(one, jobs)
    (OUT / "seed_sweep.json").write_text(json.dumps(rows, indent=1, default=float))
    print("done", len(rows), "runs")
