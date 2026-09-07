"""M5: backward-propagated CR3BP arrival-arc validation of the M4 local
insertion estimate.

Terminology note (important): every trajectory built here is a
**backward-propagated CR3BP ballistic arc** — a real, dynamically
consistent solution of the CR3BP equations of motion, propagated
backward in time from a candidate halo-insertion pre-burn state. It is
explicitly NOT called an "Earth transfer," "optimized transfer," or
"manifold transfer" anywhere in this project, because none of those
has been demonstrated: this module does not target an Earth parking
orbit, does not optimize anything, and does not compute invariant
manifolds.

Workflow (DESIGN.md M5 section documents the full rationale):
  1. Build a candidate pre-burn arrival state s_arr = [r_h, v_arr] at a
     halo phase tau, exactly as in M4 (arrival.py's Jacobi-consistent
     speed + direction models).
  2. Propagate s_arr BACKWARD in time for a specified duration.
  3. Inspect the backward arc's closest approach to Earth and Moon,
     and its Jacobi drift.
  4. Accept/reject using explicit, pre-declared geometric criteria
     (Earthward reach, Moon-distance guard) — decided from physical
     scale reasoning, not tuned to prefer any particular Delta-v.
  5. For accepted candidates, forward-propagate the backward arc's
     far-end state for the same duration and compare against the
     original s_arr (the mandatory round-trip check).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import constants as c
from .cr3bp import jacobi_constant
from .propagation import PropagationSettings, propagate

__all__ = [
    "EARTH_DISTANCE_THRESHOLD_DU",
    "MOON_GUARD_KM",
    "ArrivalArcResult",
    "propagate_backward_arc",
    "round_trip_check",
]

# --- Acceptance criteria (decided from physical/geometric scale
# reasoning BEFORE any search was run; never adjusted afterward to
# prefer a particular result) ---

# "Earthward" reach: the backward arc's closest approach to Earth must
# fall below this nondimensional distance. This is the exact example
# threshold given in this milestone's own instructions (Section 4).
EARTH_DISTANCE_THRESHOLD_DU = 0.8

# Moon guard: reject any backward arc whose closest approach to the
# Moon falls below ~5 lunar radii (Moon mean radius ~1737.4 km), a
# conservative "several lunar radii" engineering guard against
# numerically meaningless close-encounter geometry -- NOT a collision-
# probability analysis.
MOON_RADIUS_KM = 1737.4
MOON_GUARD_KM = 5.0 * MOON_RADIUS_KM  # = 8687.0 km


@dataclass
class ArrivalArcResult:
    tau: float
    dC: float
    direction_theta_rad: float
    t_back_days: float
    t_back_nondim: float
    r_h: np.ndarray
    v_arr: np.ndarray
    s_arr: np.ndarray
    success: bool
    min_dist_earth_du: float
    min_dist_moon_du: float
    final_dist_earth_du: float
    final_state: np.ndarray
    C_arr: float
    jacobi_max_drift: float
    earthward_ok: bool
    moon_guard_ok: bool
    accepted: bool
    reason: str = ""


def propagate_backward_arc(
    r_h,
    v_arr,
    mu: float,
    t_back_days: float,
    tau: float = float("nan"),
    dC: float = float("nan"),
    direction_theta_rad: float = float("nan"),
    settings: PropagationSettings | None = None,
    n_eval: int = 500,
) -> ArrivalArcResult:
    """Backward-propagate the candidate pre-burn state [r_h, v_arr] for
    `t_back_days` days and evaluate the M5 acceptance criteria.
    """
    r_h = np.asarray(r_h, dtype=float)
    v_arr = np.asarray(v_arr, dtype=float)
    s_arr = np.concatenate([r_h, v_arr])

    t_back_nondim = t_back_days / c.TU_DAYS
    t_eval = np.linspace(0.0, -t_back_nondim, n_eval)

    result = propagate(s_arr, (0.0, -t_back_nondim), mu, settings=settings, t_eval=t_eval)

    if not result.success:
        return ArrivalArcResult(
            tau=tau, dC=dC, direction_theta_rad=direction_theta_rad,
            t_back_days=t_back_days, t_back_nondim=t_back_nondim,
            r_h=r_h, v_arr=v_arr, s_arr=s_arr, success=False,
            min_dist_earth_du=float("nan"), min_dist_moon_du=float("nan"),
            final_dist_earth_du=float("nan"), final_state=np.full(6, np.nan),
            C_arr=float("nan"), jacobi_max_drift=float("nan"),
            earthward_ok=False, moon_guard_ok=False, accepted=False,
            reason=f"backward propagation failed: {result.message}",
        )

    r = result.y[0:3]
    earth = np.array([-mu, 0.0, 0.0])
    moon = np.array([1.0 - mu, 0.0, 0.0])
    dist_earth = np.linalg.norm(r.T - earth, axis=1)
    dist_moon = np.linalg.norm(r.T - moon, axis=1)

    min_dist_earth_du = float(dist_earth.min())
    min_dist_moon_du = float(dist_moon.min())
    final_dist_earth_du = float(dist_earth[-1])
    final_state = result.y[:, -1]

    C0 = jacobi_constant(s_arr, mu)
    C_t = np.array([jacobi_constant(result.y[:, i], mu) for i in range(result.y.shape[1])])
    jacobi_max_drift = float(np.max(np.abs(C_t - C0)))

    earthward_ok = min_dist_earth_du < EARTH_DISTANCE_THRESHOLD_DU
    min_dist_moon_km = min_dist_moon_du * c.DU_KM
    moon_guard_ok = min_dist_moon_km > MOON_GUARD_KM

    accepted = earthward_ok and moon_guard_ok
    reasons = []
    if not earthward_ok:
        reasons.append(f"not earthward (min_dist_earth={min_dist_earth_du:.4f} DU >= {EARTH_DISTANCE_THRESHOLD_DU} DU)")
    if not moon_guard_ok:
        reasons.append(f"Moon guard violated (min_dist_moon={min_dist_moon_km:.0f} km <= {MOON_GUARD_KM:.0f} km)")
    reason = "; ".join(reasons)

    return ArrivalArcResult(
        tau=tau, dC=dC, direction_theta_rad=direction_theta_rad,
        t_back_days=t_back_days, t_back_nondim=t_back_nondim,
        r_h=r_h, v_arr=v_arr, s_arr=s_arr, success=True,
        min_dist_earth_du=min_dist_earth_du, min_dist_moon_du=min_dist_moon_du,
        final_dist_earth_du=final_dist_earth_du, final_state=final_state,
        C_arr=float(C0), jacobi_max_drift=jacobi_max_drift,
        earthward_ok=earthward_ok, moon_guard_ok=moon_guard_ok, accepted=accepted,
        reason=reason,
    )


def round_trip_check(s_arr, far_end_state, t_back_nondim: float, mu: float, settings: PropagationSettings | None = None):
    """Mandatory round-trip verification: forward-propagate the
    backward arc's far-end state for the same duration and compare
    against the original pre-burn state `s_arr`.

    Returns (recovered_state, position_error_nondim, velocity_error_nondim).
    """
    far_end_state = np.asarray(far_end_state, dtype=float)
    result = propagate(far_end_state, (0.0, t_back_nondim), mu, settings=settings)
    if not result.success:
        raise RuntimeError(f"round-trip forward propagation failed: {result.message}")
    recovered = result.y[:, -1]
    s_arr = np.asarray(s_arr, dtype=float)
    pos_err = float(np.linalg.norm(recovered[0:3] - s_arr[0:3]))
    vel_err = float(np.linalg.norm(recovered[3:6] - s_arr[3:6]))
    return recovered, pos_err, vel_err
