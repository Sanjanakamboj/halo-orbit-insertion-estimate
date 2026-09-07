"""M4: halo-insertion phase sweep, refinement, and Δv computation.

Ties together the M3 corrected halo orbit (frozen, unmodified) with the
`arrival.py` arrival-state models to compute the local velocity-matching
insertion Δv as a function of orbital phase `tau = t/T`.

This module computes a **local insertion-Δv estimate at the corrected
halo orbit under an explicit arrival-state assumption** — it does not
compute or optimize an Earth-to-L2 transfer trajectory.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize_scalar

from .arrival import (
    ArrivalStateError,
    arrival_speed_from_jacobi,
    arrival_velocity,
    delta_v_vector,
    direction_radial_from_earth,
    direction_rotated_about_z,
    jacobi_from_state,
)
from .cr3bp import jacobi_constant
from .equilibria import find_L2
from .propagation import propagate

__all__ = [
    "InsertionCandidate",
    "halo_state_at_phase",
    "evaluate_insertion_at_phase",
    "phase_sweep",
    "refine_best_phase",
]


@dataclass
class InsertionCandidate:
    tau: float
    t_nondim: float
    t_days: float
    r_h: np.ndarray
    v_halo: np.ndarray
    v_arr: np.ndarray
    delta_v_vec: np.ndarray
    delta_v_nondim: float
    delta_v_m_s: float
    C_halo: float
    C_arr: float
    delta_C: float
    dist_from_moon_km: float
    dist_from_l2_km: float
    valid: bool
    reason: str = ""


def halo_state_at_phase(state0, period: float, mu: float, tau: float) -> np.ndarray:
    """Propagate the M3 corrected initial condition to phase `tau = t/T`
    (tau taken modulo 1.0) using the generic M2 propagator, unmodified."""
    tau_mod = tau % 1.0
    t = tau_mod * period
    if t == 0.0:
        return np.asarray(state0, dtype=float).copy()
    result = propagate(state0, (0.0, t), mu)
    if not result.success:
        raise RuntimeError(f"propagation to tau={tau!r} failed: {result.message}")
    return result.y[:, -1]


def evaluate_insertion_at_phase(
    state0,
    period: float,
    mu: float,
    tau: float,
    dC_baseline: float,
    direction_theta_rad: float = 0.0,
    direction_model: str = "radial",
) -> InsertionCandidate:
    """Evaluate the insertion candidate at phase `tau` under the
    Jacobi-consistent speed model (dC_baseline) and a chosen direction
    model ('radial': Model B rotated by direction_theta_rad about z;
    'aligned': Model A, idealized lower-bound only).
    """
    from . import constants as c

    l2 = find_L2(mu)
    state_h = halo_state_at_phase(state0, period, mu, tau)
    r_h = state_h[0:3]
    v_halo = state_h[3:6]
    C_halo = jacobi_constant(state_h, mu)

    t_nondim = (tau % 1.0) * period

    dist_from_moon_km = float(np.linalg.norm(r_h - np.array([1.0 - mu, 0.0, 0.0]))) * c.DU_KM
    dist_from_l2_km = float(np.linalg.norm(r_h - np.array([l2.x, 0.0, 0.0]))) * c.DU_KM

    try:
        C_arr = C_halo - dC_baseline
        speed = arrival_speed_from_jacobi(r_h, mu, C_arr)

        if direction_model == "aligned":
            from .arrival import direction_aligned_with_halo

            direction = direction_aligned_with_halo(v_halo)
        elif direction_model == "radial":
            base_dir = direction_radial_from_earth(r_h, mu)
            direction = direction_rotated_about_z(base_dir, direction_theta_rad)
        else:
            raise ValueError(f"unknown direction_model {direction_model!r}")

        v_arr = arrival_velocity(speed, direction)
        dv_vec = delta_v_vector(v_halo, v_arr)
        dv_nondim = float(np.linalg.norm(dv_vec))
        dv_m_s = dv_nondim * c.VSTAR_M_S

        C_arr_check = jacobi_from_state(r_h, v_arr, mu)
        delta_C = C_halo - C_arr_check

        return InsertionCandidate(
            tau=tau % 1.0,
            t_nondim=t_nondim,
            t_days=t_nondim * c.TU_DAYS,
            r_h=r_h,
            v_halo=v_halo,
            v_arr=v_arr,
            delta_v_vec=dv_vec,
            delta_v_nondim=dv_nondim,
            delta_v_m_s=dv_m_s,
            C_halo=C_halo,
            C_arr=C_arr_check,
            delta_C=delta_C,
            dist_from_moon_km=dist_from_moon_km,
            dist_from_l2_km=dist_from_l2_km,
            valid=True,
        )
    except ArrivalStateError as e:
        return InsertionCandidate(
            tau=tau % 1.0,
            t_nondim=t_nondim,
            t_days=t_nondim * c.TU_DAYS,
            r_h=r_h,
            v_halo=v_halo,
            v_arr=np.full(3, np.nan),
            delta_v_vec=np.full(3, np.nan),
            delta_v_nondim=float("nan"),
            delta_v_m_s=float("nan"),
            C_halo=C_halo,
            C_arr=float("nan"),
            delta_C=float("nan"),
            dist_from_moon_km=dist_from_moon_km,
            dist_from_l2_km=dist_from_l2_km,
            valid=False,
            reason=str(e),
        )


def phase_sweep(
    state0, period: float, mu: float, dC_baseline: float, n_points: int = 500, direction_theta_rad: float = 0.0
) -> list[InsertionCandidate]:
    """Dense phase sweep over tau in [0, 1) (endpoint excluded to avoid
    duplicating tau=0 == tau=1 in phase statistics)."""
    taus = np.linspace(0.0, 1.0, n_points, endpoint=False)
    return [
        evaluate_insertion_at_phase(state0, period, mu, tau, dC_baseline, direction_theta_rad=direction_theta_rad)
        for tau in taus
    ]


def refine_best_phase(
    state0, period: float, mu: float, dC_baseline: float, tau_bracket: tuple, direction_theta_rad: float = 0.0
) -> InsertionCandidate:
    """Bounded scalar refinement (Brent's method within `tau_bracket`)
    of the phase minimizing |Delta_v|, starting from a grid-identified
    bracket around the candidate minimum."""

    def objective(tau):
        cand = evaluate_insertion_at_phase(state0, period, mu, tau, dC_baseline, direction_theta_rad=direction_theta_rad)
        if not cand.valid:
            return 1e6  # steer optimizer away from invalid region without silently accepting it
        return cand.delta_v_nondim

    result = minimize_scalar(objective, bounds=tau_bracket, method="bounded", options={"xatol": 1e-10})
    return evaluate_insertion_at_phase(
        state0, period, mu, result.x, dC_baseline, direction_theta_rad=direction_theta_rad
    )
