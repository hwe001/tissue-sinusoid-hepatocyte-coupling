"""Typed exchange record between the transport layer and a local-model adapter.

Implements the contract of Section 2.3 of the manuscript:
    I(i,k) = {q_in, x_in, z, V, tau, theta}  ->  O(i,k) = {x_out, J, y, status}
All exchanged quantities are in SI-based units declared in UNITS. The SBML model
BIOMD0000000624 is written in mmol/L (= mol m^-3) and seconds, so the conversion
factor to the contract units is exactly 1 and is nevertheless applied explicitly.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict

import numpy as np

# Declared units of every exchanged quantity (check 1 of the contract).
UNITS = {
    "q_in": "m3 s-1",
    "x_in": "mol m-3",
    "x_out": "mol m-3",
    "z": "mol m-3",      # local oxygen concentration
    "V": "m3",           # represented control-volume (blood) volume
    "tau": "s",          # residence time in the control volume
    "J": "mol s-1",      # net species SINK flux (positive = removed from the blood stream)
    "residual": "mol s-1",
}

# SBML BIOMD0000000624 declares species in mmol/L; 1 mmol/L == 1 mol m-3.
SBML_MM_TO_MOL_M3 = 1.0


@dataclass
class ExchangeIn:
    q_in: float                      # m3 s-1
    x_in: Dict[str, float]           # species -> mol m-3 (ordering fixed by the manifest)
    z: float                         # mol m-3 oxygen (mid-control-volume)
    V: float                         # m3
    tau: float                       # s
    theta: Dict[str, float] = field(default_factory=dict)  # parameter/state context


@dataclass
class ExchangeOut:
    x_out: Dict[str, float]          # mol m-3
    J: Dict[str, float]              # mol s-1, sink positive, computed INDEPENDENTLY of x_out
    y: Dict[str, float]              # observables
    status: Dict[str, object]
    residual: Dict[str, float] = field(default_factory=dict)  # r_m = q x_in - q x_out - J_m


def closure_residual(rec_in: ExchangeIn, rec_out: ExchangeOut) -> Dict[str, float]:
    """r_m = q_in x_in,m - q_out x_out,m - J_m (mol s-1). q_out = q_in (incompressible)."""
    return {
        m: rec_in.q_in * rec_in.x_in[m] - rec_in.q_in * rec_out.x_out[m] - rec_out.J.get(m, 0.0)
        for m in rec_in.x_in
    }
