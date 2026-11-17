"""Seeded Voronoi scaffold generator (the released repository has no generator script).

Regular reference: the released 128 generators form a BCC lattice (spacing 1, all coordinates
integer or all half-integer in [-2, 1.5]); each cell is a truncated octahedron of volume 0.5.
This module reproduces that lattice (sigma = 0) and adds seeded irregular variants:
  mode 'jitter' : lattice generators displaced by N(0, sigma^2) per axis (sigma in lattice units)
  mode 'random' : 128 uniform random generators with a minimum separation (Poisson-disc-like)
The domain is the box [-2.25, 1.75]^3. Cells touching the box are closed by mirroring all
generators across the six box planes (so every original cell is bounded and clipped by the box).

Portal candidates of a cell = the 4 interior Voronoi vertices nearest to its generator
(interior = strictly inside the box). This selection rule is ASSUMED: the released geometry's
rule is not documented, and on the regular lattice all 24 vertices are equidistant.
Output follows the released network schema so load_network() reads it unchanged.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.spatial import ConvexHull, Voronoi

BOX_LO, BOX_HI = -2.25, 1.75


def bcc_generators() -> np.ndarray:
    pts = []
    for i in range(4):
        for j in range(4):
            for k in range(4):
                pts.append((-2 + i, -2 + j, -2 + k))
                pts.append((-1.5 + i, -1.5 + j, -1.5 + k))
    return np.array(pts, float)


def random_generators(n: int, rng: np.random.Generator, min_sep: float = 0.55) -> np.ndarray:
    pts: list = []
    while len(pts) < n:
        c = rng.uniform(BOX_LO, BOX_HI, 3)
        if all(np.linalg.norm(c - p) >= min_sep for p in pts):
            pts.append(c)
    return np.array(pts)


def make_generators(seed: int, mode: str, sigma: float) -> np.ndarray:
    rng = np.random.default_rng(seed)
    if mode == "random":
        return random_generators(128, rng)
    g = bcc_generators()
    if sigma > 0:
        g = np.clip(g + rng.normal(0, sigma, g.shape), BOX_LO + 0.05, BOX_HI - 0.05)
    return g + rng.normal(0, 1e-7, g.shape)          # break exact Voronoi degeneracy


def _mirror(g: np.ndarray) -> np.ndarray:
    out = [g]
    for ax in range(3):
        for plane in (BOX_LO, BOX_HI):
            m = g.copy(); m[:, ax] = 2 * plane - m[:, ax]; out.append(m)
    return np.vstack(out)


def _choose(vor, reg, gen, k):
    """k nearest vertices of a cell: strictly interior ones first, then box-boundary ones."""
    eps = 1e-6
    def strictly(v): return np.all((v > BOX_LO + eps) & (v < BOX_HI - eps))
    def within(v): return np.all((v >= BOX_LO - eps) & (v <= BOX_HI + eps))
    scored = sorted((0 if strictly(vor.vertices[vi]) else 1, float(np.linalg.norm(vor.vertices[vi] - gen)), vi)
                    for vi in reg if within(vor.vertices[vi]))
    return [vi for _, _, vi in scored[:k]]


def build(seed: int = 0, mode: str = "jitter", sigma: float = 0.0, n_portal: int = 4, unit_cv: int = 8):
    """Return (network_dict, stats_dict)."""
    g = make_generators(seed, mode, sigma)
    n = len(g)
    vor = Voronoi(_mirror(g))
    eps = 1e-6
    inside = lambda v: np.all((v > BOX_LO + eps) & (v < BOX_HI - eps))
    vol = np.zeros(n)
    portal_ids: dict = {}
    portal_pos: list = []
    incident: dict = {}
    paths, lobules, cvs = [], [], []
    for i in range(n):
        reg = vor.regions[vor.point_region[i]]
        if -1 in reg or not reg:
            raise RuntimeError(f"unbounded cell {i}")
        V = vor.vertices[reg]
        vol[i] = ConvexHull(V).volume
        chosen = _choose(vor, reg, g[i], n_portal)
        for vi in chosen:
            incident.setdefault(vi, set()).add(i)
    # portal sites (deduplicated by Voronoi vertex index), incident lobules among ORIGINAL cells
    vmap = {vi: k for k, vi in enumerate(sorted(incident))}
    true_inc: dict = {}
    for j in range(n):
        for vi in vor.regions[vor.point_region[j]]:
            if vi in vmap:
                true_inc.setdefault(vi, set()).add(j)
    for vi in sorted(incident):
        portal_pos.append({"id": vmap[vi], "position": vor.vertices[vi].tolist(), "candidate": True,
                           "incident_lobules": sorted(true_inc.get(vi, incident[vi]))})
    for i in range(n):
        reg = vor.regions[vor.point_region[i]]
        cand = [(0.0, vi) for vi in _choose(vor, reg, g[i], n_portal)]
        pids = [vmap[vi] for _, vi in cand]
        lobules.append({"id": i, "central_vein_id": i, "portal_candidate_ids": pids, "cell_model": "BIOMD0000000624", "status": "geometry_only"})
        cvs.append({"id": i, "position": g[i].tolist(), "lobule_id": i})
        for d, vi in cand:
            p0, p1 = vor.vertices[vi], g[i]
            poly = [(p0 + (p1 - p0) * t).tolist() for t in np.linspace(0, 1, unit_cv + 1)]
            paths.append({"id": f"l{i}_p{vmap[vi]}", "lobule_id": i, "portal_id": vmap[vi], "central_vein_id": i,
                          "control_volumes": unit_cv, "polyline": poly})
    net = {"schema": "lobule_sinusoid_network/regenerated", "source": f"hepcoupling.geometry seed={seed} mode={mode} sigma={sigma}",
           "lobules": lobules, "portal_sites": portal_pos, "central_vein": cvs, "central_veins": cvs, "sinusoid_paths": paths}
    L = np.array([np.linalg.norm(np.diff(np.array(p["polyline"]), axis=0), axis=1).sum() for p in paths])
    deg = np.array([len(x["incident_lobules"]) for x in portal_pos])
    stats = {"seed": seed, "mode": mode, "sigma": sigma, "cells": n, "paths": len(paths), "portal_sites": len(portal_pos),
             "volume_mean": vol.mean(), "volume_sd": vol.std(), "volume_min": vol.min(), "volume_max": vol.max(),
             "volume_sum": vol.sum(), "portal_degree_mean": deg.mean(), "portal_degree_min": int(deg.min()), "portal_degree_max": int(deg.max()),
             "path_length_mean": L.mean(), "path_length_sd": L.std(), "path_length_min": L.min(), "path_length_max": L.max()}
    return net, stats


if __name__ == "__main__":
    import sys
    net, st = build(int(sys.argv[1]) if len(sys.argv) > 1 else 0, sys.argv[2] if len(sys.argv) > 2 else "jitter",
                    float(sys.argv[3]) if len(sys.argv) > 3 else 0.0)
    print(json.dumps(st, indent=1, default=float))
