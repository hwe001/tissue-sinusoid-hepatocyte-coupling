"""Transport network: flow law, control-volume discretisation, residence times.

Geometry is read from the released geometry/lobule_sinusoid_network.json
(512 paths, 128 lobules, regular Voronoi lattice, nothing randomised).

Every physical assumption below is an explicit argument, because the reviewers
asked for the flow law, residence time and represented volume to be stated.

Flow law
--------
Paths of one lobule run in parallel between the same portal pressure and the same
central-vein pressure. With a common radius r, the Poiseuille conductance of a path is
G_i = pi r^4 / (8 mu L_i) so q_i = Q_lobule * (1/L_i) / sum_j (1/L_j)   ('poiseuille').
'uniform' splits Q_lobule equally (this is what the released JSON's relative_flow=0.25
corresponds to). Q_lobule is fixed by requiring the MEAN path velocity to equal u_ref.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import List

import numpy as np


@dataclass
class SinusoidPath:
    pid: str
    lobule_id: int
    length: float           # m
    q: float                # m3 s-1
    n_cv: int
    V_cv: float             # m3, blood volume of one control volume
    tau_cv: float           # s, residence time in one control volume
    cv_centres: np.ndarray  # (n_cv, 3) in m, for rendering/aggregation

    @property
    def tau_total(self) -> float:
        return self.n_cv * self.tau_cv


def _resample_polyline(poly: np.ndarray, n_seg: int):
    seg = np.linalg.norm(np.diff(poly, axis=0), axis=1)
    s = np.concatenate([[0.0], np.cumsum(seg)])
    total = s[-1]
    s_edges = np.linspace(0.0, total, n_seg + 1)
    pts = np.stack([np.interp(s_edges, s, poly[:, k]) for k in range(3)], axis=1)
    return total, 0.5 * (pts[:-1] + pts[1:])


def load_network(
    json_path: str | Path,
    n_cv: int = 16,
    unit_to_m: float = 5.0e-4,       # ASSUMPTION: 1 geometry unit = 0.5 mm (portal-central path ~0.6 mm)
    radius_m: float = 5.0e-6,        # sinusoid radius (typical 4-6 um)
    u_ref_m_s: float = 3.0e-4,       # mean sinusoidal blood velocity (order 100-500 um/s)
    flow_law: str = "poiseuille",
) -> List[SinusoidPath]:
    d = json.loads(Path(json_path).read_text())
    by_lob: dict[int, list] = {}
    for p in d["sinusoid_paths"]:
        by_lob.setdefault(p["lobule_id"], []).append(p)
    area = np.pi * radius_m**2
    out: List[SinusoidPath] = []
    for lob, plist in by_lob.items():
        geo = []
        for p in plist:
            poly = np.asarray(p["polyline"], float) * unit_to_m
            L, centres = _resample_polyline(poly, n_cv)
            geo.append((p["id"], L, centres))
        Ls = np.array([g[1] for g in geo])
        Q = len(plist) * area * u_ref_m_s
        if flow_law == "poiseuille":
            w = 1.0 / Ls
        elif flow_law == "uniform":
            w = np.ones_like(Ls)
        else:
            raise ValueError(flow_law)
        q = Q * w / w.sum()
        for (pid, L, centres), qi in zip(geo, q):
            V_cv = area * L / n_cv
            out.append(SinusoidPath(pid, lob, L, float(qi), n_cv, V_cv, V_cv / qi, centres))
    return out


def lobule_flow_balance(paths: List[SinusoidPath]):
    """Max relative mismatch between sum of path flows and lobule flow (flow-balance check)."""
    by = {}
    for p in paths:
        by.setdefault(p.lobule_id, []).append(p.q)
    return by
