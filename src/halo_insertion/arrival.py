"""M4: transfer-arrival state models for the Earth-Moon L2 halo
insertion Δv estimate.

This module answers one question: **given a position `r_h` on the
corrected M3 halo orbit, what velocity does an incoming transfer
trajectory plausibly have there?** It does NOT propagate an actual
Earth-departure trajectory — no state here is described as "from
Earth" in the sense of an integrated trajectory. It is an explicit,
documented **arrival-state assumption**, exactly as DESIGN.md's M4
section requires.

Speed model (Jacobi-consistent, DESIGN.md M4 "baseline arrival model")
------------------------------------------------------------------------
For an assumed arrival Jacobi constant `C_arr = C_halo - dC`, the
arrival speed satisfies

    v_arr^2 = 2*Omega(r_h) - C_arr

Points where `2*Omega(r_h) - C_arr < 0` are physically invalid (no
real speed satisfies the assumed Jacobi constant there) and are
rejected explicitly, never silently clipped to zero.

Direction models
------------------
A. **Aligned** — arrival velocity is parallel to the local halo
   velocity: `v_arr_dir = v_halo / |v_halo|`. This is an
   **idealized/nonphysical lower-bound case**: combined with a free
   arrival speed this would let Delta-v collapse to a difference of
   scalars only (or to exactly zero if speed is also matched) — it is
   used ONLY as a lower-bound sanity check (DESIGN.md M4, Section 8's
   degeneracy guard), never presented as a meaningful insertion
   solution.
B. **Radial-from-Earth** (baseline) — arrival velocity points along the
   Earth-to-insertion-point radial direction:
   `v_arr_dir = (r_h - r_Earth) / |r_h - r_Earth|`, with
   `r_Earth = (-mu, 0, 0)`. This models a transfer arriving generally
   outbound from the Earth region, without claiming an actual
   propagated trajectory.
C. **Direction sweep** — model B's direction rotated by a signed angle
   `theta` about the synodic z-axis (i.e., within the local x-y plane):
   `v_arr_dir(theta) = R_z(theta) @ v_arr_dir_B`.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .cr3bp import effective_potential

__all__ = [
    "ArrivalStateError",
    "earth_position",
    "arrival_speed_from_jacobi",
    "direction_aligned_with_halo",
    "direction_radial_from_earth",
    "direction_rotated_about_z",
    "arrival_velocity",
    "delta_v_vector",
    "jacobi_from_state",
]


class ArrivalStateError(ValueError):
    """Raised when an arrival-state assumption is physically invalid
    at the requested position (e.g., negative v_arr^2)."""


def earth_position(mu: float) -> np.ndarray:
    return np.array([-mu, 0.0, 0.0])


def arrival_speed_from_jacobi(r_h, mu: float, C_arr: float) -> float:
    """v_arr = sqrt(2*Omega(r_h) - C_arr). Raises ArrivalStateError if
    the right-hand side is negative (no real arrival speed exists for
    this Jacobi assumption at this position)."""
    x, y, z = r_h
    omega = effective_potential(x, y, z, mu)
    v_arr_sq = 2.0 * omega - C_arr
    if v_arr_sq < 0:
        raise ArrivalStateError(
            f"2*Omega(r_h)-C_arr = {v_arr_sq!r} < 0 at r_h={tuple(r_h)!r}; "
            "no real arrival speed satisfies this Jacobi assumption here"
        )
    return math.sqrt(v_arr_sq)


def direction_aligned_with_halo(v_halo) -> np.ndarray:
    """Model A (idealized lower-bound only): unit vector along v_halo."""
    v_halo = np.asarray(v_halo, dtype=float)
    norm = np.linalg.norm(v_halo)
    if norm < 1e-14:
        raise ArrivalStateError("halo velocity is ~zero; aligned direction undefined")
    return v_halo / norm


def direction_radial_from_earth(r_h, mu: float) -> np.ndarray:
    """Model B (baseline): unit vector from Earth (-mu,0,0) to r_h."""
    r_h = np.asarray(r_h, dtype=float)
    rel = r_h - earth_position(mu)
    norm = np.linalg.norm(rel)
    if norm < 1e-14:
        raise ArrivalStateError("r_h coincides with Earth's position; radial direction undefined")
    return rel / norm


def direction_rotated_about_z(base_dir, theta_rad: float) -> np.ndarray:
    """Model C: rotate `base_dir` by `theta_rad` about the synodic z-axis."""
    base_dir = np.asarray(base_dir, dtype=float)
    cos_t, sin_t = math.cos(theta_rad), math.sin(theta_rad)
    Rz = np.array(
        [
            [cos_t, -sin_t, 0.0],
            [sin_t, cos_t, 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    rotated = Rz @ base_dir
    norm = np.linalg.norm(rotated)
    return rotated / norm  # Rz is orthogonal; norm should already be ~1, renormalize defensively


def arrival_velocity(speed: float, direction) -> np.ndarray:
    """v_arr = speed * unit(direction). Direction is renormalized
    defensively (callers are expected to already pass unit vectors)."""
    direction = np.asarray(direction, dtype=float)
    norm = np.linalg.norm(direction)
    if abs(norm - 1.0) > 1e-6:
        direction = direction / norm
    return speed * direction


def delta_v_vector(v_halo, v_arr) -> np.ndarray:
    """Delta_v = v_halo - v_arr (DESIGN.md Section 7's definition, unchanged since M1)."""
    return np.asarray(v_halo, dtype=float) - np.asarray(v_arr, dtype=float)


def jacobi_from_state(r, v, mu: float) -> float:
    """C = 2*Omega(r) - |v|^2, computed directly from raw position/velocity
    components (independent of cr3bp.jacobi_constant's 6-vector packing,
    used as an independent cross-check in tests/scripts)."""
    x, y, z = r
    omega = effective_potential(x, y, z, mu)
    v = np.asarray(v, dtype=float)
    return 2.0 * omega - float(np.dot(v, v))
