"""Local-model adapters.

SBMLAdapter     executes BIOMD0000000624 (Sluka et al. 2016) with libRoadRunner.
CellMLAdapter   RE-IMPLEMENTS the archived Schmitmeier 2007 CellML equations in Python.
                It does not execute the CellML file. Name it that way in the paper.

Both expose the same step(ExchangeIn) -> ExchangeOut interface.
"""
from __future__ import annotations

from typing import Dict, Optional

import numpy as np

from .exchange import ExchangeIn, ExchangeOut, SBML_MM_TO_MOL_M3, closure_residual


class SBMLAdapter:
    """One control volume = one SBML instance run for tau_eff = phi * tau seconds.

    phi (hepatocyte-to-sinusoid volume ratio) rescales reaction time: dC/dt = phi * R(C).
    phi = 1 means the SBML compartment holds exactly the blood of the control volume.
    The CYP2E1 oxygen scaling is an EXTERNAL per-control-volume override of the SBML
    parameter Vmax_2E1_APAP (as in Sluka et al.), not an edit of the model.
    """

    model_id = "BIOMD0000000624"

    def __init__(self, sbml_path: str, phi: float = 1.0, n_points: int = 61,
                 cyp_o2_half: Optional[float] = None, o2_ref: Optional[float] = None,
                 rtol: float = 1e-10, atol: float = 1e-14):
        import roadrunner

        self.rr = roadrunner.RoadRunner(str(sbml_path))
        self.species = list(self.rr.model.getFloatingSpeciesIds())
        self.rxn = list(self.rr.model.getReactionIds())
        S = self.rr.getFullStoichiometryMatrix()
        rows = list(S.rownames)
        self.S = np.array([[S[rows.index(sp), j] for j in range(len(self.rxn))] for sp in self.species])
        self.rr.integrator.relative_tolerance = rtol
        self.rr.integrator.absolute_tolerance = atol
        self.rr.timeCourseSelections = ["time"] + [f"[{s}]" for s in self.species] + self.rxn
        self.vmax0 = self.rr.getValue("Vmax_2E1_APAP")
        self.phi, self.n_points = phi, n_points
        self.cyp_o2_half, self.o2_ref = cyp_o2_half, o2_ref
        self.rtol, self.atol = rtol, atol

    def cyp_scale(self, z: float) -> float:
        """f(O2) = [z/(K+z)] / [z_ref/(K+z_ref)]; 1.0 if scaling is switched off."""
        if self.cyp_o2_half is None:
            return 1.0
        K = self.cyp_o2_half
        return (z / (K + z)) / (self.o2_ref / (K + self.o2_ref))

    def step(self, rec: ExchangeIn) -> ExchangeOut:
        rr = self.rr
        rr.reset()
        for sp in self.species:
            rr.setValue(f"[{sp}]", rec.x_in[sp] / SBML_MM_TO_MOL_M3)
        rr.setValue("Vmax_2E1_APAP", self.vmax0 * self.cyp_scale(rec.z))
        t_end = self.phi * rec.tau
        res = rr.simulate(0.0, t_end, self.n_points)
        names = list(res.colnames)
        t = res[:, 0]
        ns = len(self.species)
        conc = res[:, 1:1 + ns]
        rates = res[:, 1 + ns:1 + ns + len(self.rxn)]
        x_out_arr = conc[-1] * SBML_MM_TO_MOL_M3
        # independent flux: integral of stoichiometry x reaction rates (Simpson), NOT x_out - x_in
        dx = np.array([_simpson(rates @ self.S[i], t) for i in range(ns)]) * SBML_MM_TO_MOL_M3
        J = {sp: -rec.q_in * dx[i] for i, sp in enumerate(self.species)}  # sink positive
        out = ExchangeOut(
            x_out={sp: float(x_out_arr[i]) for i, sp in enumerate(self.species)},
            J=J,
            y={"rate_" + r: float(rates[-1, j]) for j, r in enumerate(self.rxn)},
            status={"ok": True, "model": self.model_id, "units": "mol m-3, s",
                    "solver": f"CVODE rtol={self.rtol:g} atol={self.atol:g}",
                    "cyp_scale": self.cyp_scale(rec.z), "tau_eff": t_end},
        )
        out.residual = closure_residual(rec, out)
        return out


def _simpson(y: np.ndarray, t: np.ndarray) -> float:
    from scipy.integrate import simpson
    return float(simpson(y, x=t))


class CellMLAdapter:
    """Re-implementation of the archived Schmitmeier 2007 CellML (20% oxygen version).

    Archived equations (time in hours as declared in the file; p2 declared 'per_day'):
        dO2/dt   = p1*F1 - p2*O2
        dH2O2/dt = p2*O2 - p3*H2O2 - p4*H2O2
        V        = Vmax * ki / (ki + H2O2)
    Oxygen is a state variable in the archived model. Here oxygen is CLAMPED to the value
    supplied by the transport layer (its ODE is removed), which changes the model
    structure. Units 'CONCENTRATIONO2', 'CONFIRMp1', 'CONFIRMV', 'CONFIRMki' are
    unresolved placeholders in the archive, so the transported oxygen is mapped to model
    units by o2_model = O2_REF_MODEL * (z / z_ref); this mapping is an assumption.
    """

    model_id = "schmitmeier_2007 (re-implementation)"
    p1, F1, p2, p3, p4 = 1.43385e-5, 20.0, 0.117425199, 2.552293418, 8.379310334
    Vmax, ki = 10400000.0, 8.982e-6
    O2_INIT, H2O2_INIT = 0.00235897, 0.000027947

    def __init__(self, z_ref: float, h2o2_init: Optional[float] = None, time_to_hours: float = 1.0 / 3600.0):
        self.z_ref = z_ref
        self.o2_ref_model = self.p1 * self.F1 / self.p2   # archived steady state of O2
        self.h2o2_init = self.H2O2_INIT if h2o2_init is None else h2o2_init
        self.time_to_hours = time_to_hours

    def step(self, rec: ExchangeIn, h2o2_in: Optional[float] = None) -> ExchangeOut:
        o2 = self.o2_ref_model * rec.z / self.z_ref
        k = self.p3 + self.p4
        h_inf = self.p2 * o2 / k
        h0 = self.h2o2_init if h2o2_in is None else h2o2_in
        dt = rec.tau * self.time_to_hours
        h = h_inf + (h0 - h_inf) * np.exp(-k * dt)       # exact solution of the linear ODE
        V = self.Vmax * self.ki / (self.ki + h)
        out = ExchangeOut(
            x_out=dict(rec.x_in),                          # no transported species are altered
            J={m: 0.0 for m in rec.x_in},                  # => closure check is vacuous here
            y={"H2O2": float(h), "V": float(V), "O2_model": float(o2)},
            status={"ok": True, "model": self.model_id, "units": "archived (unresolved CONFIRM*)",
                    "note": "O2 clamped; conservation check vacuous for this adapter"},
        )
        out.residual = closure_residual(rec, out)
        return out
