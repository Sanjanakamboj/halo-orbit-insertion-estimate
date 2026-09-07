"""CR3BP state Jacobian (analytic second derivatives of Omega), the 6x6
variational equations, and state-transition-matrix (STM) propagation.

This module is reusable/generic — it knows nothing about halo orbits.
Halo-specific seeding and differential correction live in
`halo_seed.py` and `differential_correction.py`.

Convention (unchanged from M2, see cr3bp.py):
    state s = [x, y, z, xdot, ydot, zdot]
    A = [[ 0,   0,   0,   1, 0, 0],
         [ 0,   0,   0,   0, 1, 0],
         [ 0,   0,   0,   0, 0, 1],
         [Uxx, Uxy, Uxz,  0, 2, 0],
         [Uxy, Uyy, Uyz, -2, 0, 0],
         [Uxz, Uyz, Uzz,  0, 0, 0]]
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .cr3bp import CR3BPStateError, cr3bp_rhs, primary_distances

__all__ = [
    "omega_hessian",
    "state_jacobian_A",
    "stm_rhs",
    "augmented_rhs",
    "STMPropagationResult",
    "propagate_with_stm",
]


def omega_hessian(x: float, y: float, z: float, mu: float) -> np.ndarray:
    """Analytic 3x3 Hessian of the effective potential Omega:

        [[Uxx, Uxy, Uxz],
         [Uxy, Uyy, Uyz],
         [Uxz, Uyz, Uzz]]

    Closed-form expressions (see DESIGN.md Section "M3" for the
    derivation from Omega_x, Omega_y, Omega_z):

        Uxx = 1 - (1-mu)/r1^3 + 3(1-mu)*dx1^2/r1^5 - mu/r2^3 + 3*mu*dx2^2/r2^5
        Uyy = 1 - (1-mu)/r1^3 + 3(1-mu)*y^2 /r1^5   - mu/r2^3 + 3*mu*y^2 /r2^5
        Uzz =   - (1-mu)/r1^3 + 3(1-mu)*z^2 /r1^5   - mu/r2^3 + 3*mu*z^2 /r2^5
        Uxy = 3(1-mu)*dx1*y/r1^5 + 3*mu*dx2*y/r2^5
        Uxz = 3(1-mu)*dx1*z/r1^5 + 3*mu*dx2*z/r2^5
        Uyz = 3(1-mu)*y*z/r1^5   + 3*mu*y*z/r2^5

    with dx1 = x+mu, dx2 = x-1+mu.
    """
    r1, r2 = primary_distances(x, y, z, mu)
    dx1 = x + mu
    dx2 = x - 1.0 + mu
    r1_3, r1_5 = r1**3, r1**5
    r2_3, r2_5 = r2**3, r2**5

    Uxx = (
        1.0
        - (1.0 - mu) / r1_3
        + 3.0 * (1.0 - mu) * dx1**2 / r1_5
        - mu / r2_3
        + 3.0 * mu * dx2**2 / r2_5
    )
    Uyy = (
        1.0
        - (1.0 - mu) / r1_3
        + 3.0 * (1.0 - mu) * y**2 / r1_5
        - mu / r2_3
        + 3.0 * mu * y**2 / r2_5
    )
    Uzz = (
        -(1.0 - mu) / r1_3
        + 3.0 * (1.0 - mu) * z**2 / r1_5
        - mu / r2_3
        + 3.0 * mu * z**2 / r2_5
    )
    Uxy = 3.0 * (1.0 - mu) * dx1 * y / r1_5 + 3.0 * mu * dx2 * y / r2_5
    Uxz = 3.0 * (1.0 - mu) * dx1 * z / r1_5 + 3.0 * mu * dx2 * z / r2_5
    Uyz = 3.0 * (1.0 - mu) * y * z / r1_5 + 3.0 * mu * y * z / r2_5

    return np.array(
        [
            [Uxx, Uxy, Uxz],
            [Uxy, Uyy, Uyz],
            [Uxz, Uyz, Uzz],
        ]
    )


def state_jacobian_A(state, mu: float) -> np.ndarray:
    """The 6x6 CR3BP state Jacobian A(state) = d(sdot)/d(s)."""
    s = np.asarray(state, dtype=float)
    x, y, z = s[0], s[1], s[2]
    hess = omega_hessian(x, y, z, mu)

    A = np.zeros((6, 6))
    A[0:3, 3:6] = np.eye(3)
    A[3:6, 0:3] = hess
    A[3, 4] = 2.0
    A[4, 3] = -2.0
    return A


def stm_rhs(t: float, Phi_flat, state, mu: float) -> np.ndarray:
    """Phi_dot = A(state) @ Phi, given the (fixed, externally supplied)
    state at this instant. Provided mainly for testing; the production
    path uses `augmented_rhs`, which computes state and Phi together.
    """
    A = state_jacobian_A(state, mu)
    Phi = np.asarray(Phi_flat, dtype=float).reshape(6, 6)
    return (A @ Phi).flatten()


def augmented_rhs(t: float, Y, mu: float) -> np.ndarray:
    """RHS for the augmented state Y = [state(6), Phi.flatten()(36)],
    suitable for scipy.integrate.solve_ivp.

    sdot = cr3bp_rhs(t, state, mu)
    Phi_dot = A(state) @ Phi
    """
    Y = np.asarray(Y, dtype=float)
    if Y.shape != (42,):
        raise CR3BPStateError(f"augmented state must have shape (42,), got {Y.shape}")

    state = Y[0:6]
    Phi = Y[6:].reshape(6, 6)

    sdot = cr3bp_rhs(t, state, mu)
    A = state_jacobian_A(state, mu)
    Phi_dot = A @ Phi

    return np.concatenate([sdot, Phi_dot.flatten()])


@dataclass
class STMPropagationResult:
    """Result of propagate_with_stm: raw OdeResult plus convenience views."""

    ode_result: "object"
    t: np.ndarray
    states: np.ndarray  # shape (6, n_times)
    Phi_final: np.ndarray  # shape (6, 6), STM at the final time


def propagate_with_stm(state0, t_span, mu: float, settings=None, t_eval=None) -> STMPropagationResult:
    """Propagate state + STM together from Phi(0) = I6.

    Parameters mirror `propagation.propagate`; returns the raw
    OdeResult plus convenience `states` (6 x n) and `Phi_final` (6x6)
    views. Caller must check `ode_result.success`.
    """
    from scipy.integrate import solve_ivp

    from .propagation import DEFAULT_ATOL, DEFAULT_METHOD, DEFAULT_RTOL, PropagationSettings

    if settings is None:
        settings = PropagationSettings(method=DEFAULT_METHOD, rtol=DEFAULT_RTOL, atol=DEFAULT_ATOL)

    state0 = np.asarray(state0, dtype=float)
    Phi0 = np.eye(6)
    Y0 = np.concatenate([state0, Phi0.flatten()])

    result = solve_ivp(
        fun=augmented_rhs,
        t_span=t_span,
        y0=Y0,
        method=settings.method,
        rtol=settings.rtol,
        atol=settings.atol,
        dense_output=settings.dense_output,
        t_eval=t_eval,
        args=(mu,),
    )

    states = result.y[0:6, :]
    Phi_final = result.y[6:, -1].reshape(6, 6)

    return STMPropagationResult(ode_result=result, t=result.t, states=states, Phi_final=Phi_final)
