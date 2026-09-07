"""Normalized Earth-Moon CR3BP dynamics: effective potential, equations of
motion, and the Jacobi constant.

Convention (matches DESIGN.md Sections 4-5 exactly):
    - Earth at (-mu, 0, 0), Moon at (1-mu, 0, 0)
    - state s = [x, y, z, xdot, ydot, zdot]
    - Omega(x,y,z) = 0.5*(x^2+y^2) + (1-mu)/r1 + mu/r2
    - xddot - 2*ydot = dOmega/dx  =>  xddot = 2*ydot + dOmega/dx
    - yddot + 2*xdot = dOmega/dy  =>  yddot = -2*xdot + dOmega/dy
    - zddot          = dOmega/dz

No finite differences are used in the production right-hand side —
dOmega/dx, dOmega/dy, dOmega/dz are analytic closed-form expressions.
"""

from __future__ import annotations

import math

import numpy as np

__all__ = [
    "CR3BPStateError",
    "primary_distances",
    "effective_potential",
    "effective_potential_gradient",
    "cr3bp_rhs",
    "jacobi_constant",
]


class CR3BPStateError(ValueError):
    """Raised when a CR3BP state vector or mu value is invalid."""


def _validate_mu(mu: float) -> None:
    if not math.isfinite(mu):
        raise CR3BPStateError(f"mu must be finite, got {mu!r}")
    if not (0.0 < mu < 1.0):
        raise CR3BPStateError(f"mu must satisfy 0 < mu < 1, got {mu!r}")


def _validate_state(state) -> np.ndarray:
    s = np.asarray(state, dtype=float)
    if s.shape != (6,):
        raise CR3BPStateError(
            f"state must have shape (6,) = [x,y,z,xdot,ydot,zdot], got shape {s.shape}"
        )
    if not np.all(np.isfinite(s)):
        raise CR3BPStateError(f"state contains non-finite values: {s!r}")
    return s


def primary_distances(x: float, y: float, z: float, mu: float) -> tuple[float, float]:
    """Return (r1, r2): distances from Earth and Moon respectively.

    r1 = sqrt((x+mu)^2 + y^2 + z^2)      distance from Earth (-mu, 0, 0)
    r2 = sqrt((x-1+mu)^2 + y^2 + z^2)    distance from Moon (1-mu, 0, 0)
    """
    r1 = math.sqrt((x + mu) ** 2 + y**2 + z**2)
    r2 = math.sqrt((x - 1.0 + mu) ** 2 + y**2 + z**2)
    return r1, r2


def _check_nonsingular(r1: float, r2: float, mu: float) -> None:
    # A collision singularity makes Omega and its gradient blow up.
    # Fail clearly rather than silently returning inf/nan.
    tol = 1e-12
    if r1 < tol:
        raise CR3BPStateError(
            f"state coincides with Earth primary (r1={r1:.3e} < {tol:.0e}); "
            "effective potential is singular"
        )
    if r2 < tol:
        raise CR3BPStateError(
            f"state coincides with Moon primary (r2={r2:.3e} < {tol:.0e}); "
            "effective potential is singular"
        )


def effective_potential(x: float, y: float, z: float, mu: float) -> float:
    """Omega(x,y,z) = 0.5*(x^2+y^2) + (1-mu)/r1 + mu/r2."""
    _validate_mu(mu)
    r1, r2 = primary_distances(x, y, z, mu)
    _check_nonsingular(r1, r2, mu)
    return 0.5 * (x**2 + y**2) + (1.0 - mu) / r1 + mu / r2


def effective_potential_gradient(
    x: float, y: float, z: float, mu: float
) -> tuple[float, float, float]:
    """Analytic (dOmega/dx, dOmega/dy, dOmega/dz), closed form (no finite differences).

        dOmega/dx = x - (1-mu)*(x+mu)/r1^3   - mu*(x-1+mu)/r2^3
        dOmega/dy = y - (1-mu)*y/r1^3        - mu*y/r2^3
        dOmega/dz =   - (1-mu)*z/r1^3        - mu*z/r2^3
    """
    _validate_mu(mu)
    r1, r2 = primary_distances(x, y, z, mu)
    _check_nonsingular(r1, r2, mu)
    r1_3 = r1**3
    r2_3 = r2**3

    dOdx = x - (1.0 - mu) * (x + mu) / r1_3 - mu * (x - 1.0 + mu) / r2_3
    dOdy = y - (1.0 - mu) * y / r1_3 - mu * y / r2_3
    dOdz = -(1.0 - mu) * z / r1_3 - mu * z / r2_3
    return dOdx, dOdy, dOdz


def cr3bp_rhs(t: float, state, mu: float) -> np.ndarray:
    """First-order CR3BP right-hand side, suitable for scipy.integrate.solve_ivp.

    Parameters
    ----------
    t : float
        Time (unused — the CR3BP is autonomous — but kept for solve_ivp's
        f(t, y) signature).
    state : array-like, shape (6,)
        [x, y, z, xdot, ydot, zdot] in normalized units.
    mu : float
        Earth-Moon mass ratio, 0 < mu < 1.

    Returns
    -------
    np.ndarray, shape (6,)
        [xdot, ydot, zdot, xddot, yddot, zddot]

    Raises
    ------
    CR3BPStateError
        If `state` has the wrong shape, contains non-finite values, `mu`
        is invalid, or the state is at a collision singularity.
    """
    _validate_mu(mu)
    s = _validate_state(state)
    x, y, z, xdot, ydot, zdot = s

    dOdx, dOdy, dOdz = effective_potential_gradient(x, y, z, mu)

    xddot = 2.0 * ydot + dOdx
    yddot = -2.0 * xdot + dOdy
    zddot = dOdz

    return np.array([xdot, ydot, zdot, xddot, yddot, zddot], dtype=float)


def jacobi_constant(state, mu: float) -> float:
    """C = 2*Omega(x,y,z) - (xdot^2 + ydot^2 + zdot^2).

    Conserved under unforced CR3BP propagation; a diagnostic only — this
    function never modifies a trajectory to "enforce" conservation.
    """
    _validate_mu(mu)
    s = _validate_state(state)
    x, y, z, xdot, ydot, zdot = s
    omega = effective_potential(x, y, z, mu)
    speed_sq = xdot**2 + ydot**2 + zdot**2
    return 2.0 * omega - speed_sq
