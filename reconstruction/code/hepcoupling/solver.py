"""Plug-flow marching of transported states along each path, with oxygen computed (not imposed).

Per control volume i of a path (steady plug flow, Section 2.2):
    oxygen:  dz/dt = -k_O2 z   over tau  (first-order consumption, k_O2 set from a target
             mean single-path oxygen extraction E_O2 at the reference residence time)
    local model: advances x_in -> x_out over tau (adapter), returning an independent flux J
The outgoing state of one control volume is the incoming state of the next.
A lobule outlet is the flow-weighted mixture of its paths' outlets.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List

import numpy as np

from .exchange import ExchangeIn
from .network import SinusoidPath


@dataclass
class Inlet:
    apap: float = 0.07           # mol m-3 (= 70 uM, therapeutic plasma level ~10.6 mg/L)
    o2_pv: float = 6.5           # mol m-3 total blood O2, portal vein (~73% saturation) [to be cited]
    o2_ha: float = 8.5           # mol m-3, hepatic artery (~95% saturation)           [to be cited]
    frac_ha: float = 0.25        # hepatic-artery share of inflow                       [to be cited]
    gsh: float = 10.0            # mol m-3 baseline (= GSHmax of the SBML model)

    @property
    def o2(self) -> float:
        return (1 - self.frac_ha) * self.o2_pv + self.frac_ha * self.o2_ha

    def species(self) -> Dict[str, float]:
        return {"APAP": self.apap, "NAPQI": 0.0, "GSH": self.gsh,
                "NAPQIGSH": 0.0, "APAPconj_Glu": 0.0, "APAPconj_Sul": 0.0}


def k_o2_from_extraction(extraction: float, tau_ref: float) -> float:
    """First-order rate giving single-pass oxygen extraction E over residence time tau_ref."""
    return -np.log(1.0 - extraction) / tau_ref


@dataclass
class PathResult:
    pid: str
    lobule_id: int
    q: float
    x: np.ndarray                  # (n_cv+1, n_species) concentrations at CV interfaces
    z: np.ndarray                  # (n_cv+1,) oxygen
    J: np.ndarray                  # (n_cv, n_species) sink fluxes, mol s-1
    resid: np.ndarray              # (n_cv, n_species) closure residuals, mol s-1
    y: List[dict] = field(default_factory=list)


def march_path(path: SinusoidPath, inlet: Inlet, adapter, k_o2: float, species: List[str],
               second_adapter=None) -> PathResult:
    x = np.zeros((path.n_cv + 1, len(species)))
    z = np.zeros(path.n_cv + 1)
    J = np.zeros((path.n_cv, len(species)))
    R = np.zeros_like(J)
    xin = inlet.species()
    x[0] = [xin[s] for s in species]
    z[0] = inlet.o2
    ys: List[dict] = []
    h2o2 = None
    for i in range(path.n_cv):
        z_mid = z[i] * np.exp(-0.5 * k_o2 * path.tau_cv)
        rec = ExchangeIn(q_in=path.q, x_in={s: x[i, j] for j, s in enumerate(species)},
                         z=z_mid, V=path.V_cv, tau=path.tau_cv)
        out = adapter.step(rec)
        if not out.status["ok"]:
            raise RuntimeError(f"adapter failure on {path.pid} cv {i}: {out.status}")
        x[i + 1] = [out.x_out[s] for s in species]
        z[i + 1] = z[i] * np.exp(-k_o2 * path.tau_cv)
        J[i] = [out.J[s] for s in species]
        R[i] = [out.residual[s] for s in species]
        y = dict(out.y)
        if second_adapter is not None:
            o2 = second_adapter.step(rec, h2o2_in=h2o2)
            h2o2 = o2.y["H2O2"]
            y.update({"cellml_" + k: v for k, v in o2.y.items()})
        ys.append(y)
    return PathResult(path.pid, path.lobule_id, path.q, x, z, J, R, ys)


def run_network(paths: List[SinusoidPath], inlet: Inlet, adapter, k_o2: float, species: List[str],
                second_adapter=None) -> List[PathResult]:
    return [march_path(p, inlet, adapter, k_o2, species, second_adapter) for p in paths]


def lobule_outlets(results: List[PathResult]) -> Dict[int, np.ndarray]:
    """Flow-weighted outlet mixture per lobule (central vein)."""
    acc: Dict[int, list] = {}
    for r in results:
        acc.setdefault(r.lobule_id, []).append((r.q, r.x[-1]))
    return {lob: sum(q * x for q, x in v) / sum(q for q, _ in v) for lob, v in acc.items()}


def well_mixed_reference(paths: List[SinusoidPath], inlet: Inlet, adapter, species: List[str]) -> Dict[int, np.ndarray]:
    """Same SBML model as ONE well-stirred compartment per lobule with the same total
    residence time (flow-weighted mean path residence time) and the same total flow.
    This is the baseline the reviewer asked for (cf. Dichamp et al. 2023)."""
    from .exchange import ExchangeIn
    by: Dict[int, list] = {}
    for p in paths:
        by.setdefault(p.lobule_id, []).append(p)
    xin = inlet.species()
    out: Dict[int, np.ndarray] = {}
    for lob, ps in by.items():
        Q = sum(p.q for p in ps)
        tau = sum(p.q * p.tau_total for p in ps) / Q
        rec = ExchangeIn(q_in=Q, x_in=dict(xin), z=inlet.o2, V=Q * tau, tau=tau)
        o = adapter.step(rec)
        out[lob] = np.array([o.x_out[s] for s in species])
    return out
