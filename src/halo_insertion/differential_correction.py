"""Symmetry-based Newton differential correction of an Earth-Moon L2
halo-orbit seed to a genuine numerically periodic CR3BP orbit.

Formulation (documented fully in DESIGN.md's M3 section)
-----------------------------------------------------------
Free variables:   x0, ydot0        (z0 held fixed -> selects the family member)
Target residuals: xdot(T/2) = 0,   zdot(T/2) = 0
Symmetric IC form: y0 = 0, xdot0 = 0, zdot0 = 0 (enforced by construction,
                    not corrected)

At the half-period T/2 (the first non-trivial y=0 rotating-frame plane
crossing after t=0 — see `half_period_event`), the corrector uses the
state-transition matrix Phi(T/2) to build a linear map from a small
change in (x0, ydot0) to the resulting change in (xdot(T/2), zdot(T/2)),
inverts it (least-squares if the 2x2 block is ill-conditioned), and
takes a Newton step. Iterates until both residuals are below tolerance
or `max_iter` is reached; raises on failure to converge or on a
singular/ill-conditioned correction matrix, rather than silently
returning an unconverged state.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .cr3bp import cr3bp_rhs
from .propagation import DEFAULT_ATOL, DEFAULT_METHOD, DEFAULT_RTOL
from .variational import augmented_rhs

__all__ = [
    "DifferentialCorrectionError",
    "CorrectionIteration",
    "CorrectionResult",
    "half_period_event",
    "propagate_half_period",
    "differential_correct_halo",
]


class DifferentialCorrectionError(RuntimeError):
    """Raised when the Newton corrector fails to converge or the
    correction matrix is singular/ill-conditioned."""


def half_period_event(t, Y, mu):
    """Event function: y = 0 rotating-frame crossing (Y is the 42-vector
    augmented state [state(6), Phi.flatten()(36)]).

    solve_ivp's event machinery ignores an event at t=0 by construction
    (it only reports sign changes strictly after integration begins),
    but as extra protection this project's driver always propagates
    from a small offset time before starting to look for the crossing,
    and this function's `terminal`/`direction` attributes are set by
    the caller (`propagate_half_period`) to select the correct (first,
    downward) crossing rather than the trivial t=0 start.
    """
    return Y[1]  # y


half_period_event.direction = -1.0  # y decreasing through zero
half_period_event.terminal = True


@dataclass
class HalfPeriodResult:
    t_half: float
    Y_half: np.ndarray  # 42-vector: state + Phi, at t_half
    ode_result: object


def propagate_half_period(
    state0, mu: float, t_max: float = 6.0, rtol=DEFAULT_RTOL, atol=DEFAULT_ATOL, t_warmup: float = 0.05
) -> HalfPeriodResult:
    """Propagate the augmented (state+STM) system from `state0` (with
    Phi(0)=I) until the first genuine y=0 crossing strictly after t=0.

    The initial condition itself has y(0)=0 by construction (the
    symmetric halo IC form), and ydot(0) is generally nonzero, so y(t)
    leaves zero immediately in one direction. A naive event search
    starting exactly at t=0 would immediately (mis)detect this trivial
    departure as a "crossing." To avoid that, this function propagates
    a short **warm-up** interval `[0, t_warmup]` with event detection
    disabled, then resumes from the warmed-up state (now with y != 0
    and away from the seam) with the y=0 event armed, so only the
    physically intended *return* crossing is detected. `t_warmup` is
    chosen small relative to any expected half-period in this system
    (which is O(1) nondimensional time) and is tested explicitly
    (`test_half_period_event_does_not_trigger_at_t0`).
    """
    from scipy.integrate import solve_ivp

    state0 = np.asarray(state0, dtype=float)
    Phi0 = np.eye(6)
    Y0 = np.concatenate([state0, Phi0.flatten()])

    if t_warmup <= 0 or t_warmup >= t_max:
        raise ValueError(f"t_warmup={t_warmup!r} must satisfy 0 < t_warmup < t_max={t_max!r}")

    warmup = solve_ivp(
        fun=lambda t, Y: augmented_rhs(t, Y, mu),
        t_span=(0.0, t_warmup),
        y0=Y0,
        method=DEFAULT_METHOD,
        rtol=rtol,
        atol=atol,
        dense_output=False,
    )
    if not warmup.success:
        raise DifferentialCorrectionError(f"warm-up propagation failed: {warmup.message}")

    Y_warm = warmup.y[:, -1]

    def event(t, Y):
        return half_period_event(t, Y, mu)

    event.direction = -1.0
    event.terminal = True

    result = solve_ivp(
        fun=lambda t, Y: augmented_rhs(t, Y, mu),
        t_span=(t_warmup, t_max),
        y0=Y_warm,
        method=DEFAULT_METHOD,
        rtol=rtol,
        atol=atol,
        events=event,
        dense_output=False,
    )

    if not result.success:
        raise DifferentialCorrectionError(f"half-period propagation failed: {result.message}")
    if len(result.t_events[0]) == 0:
        raise DifferentialCorrectionError(
            "no y=0 crossing found within t_max; increase t_max or check the seed"
        )

    t_half = result.t_events[0][0]
    Y_half = result.y_events[0][0]
    return HalfPeriodResult(t_half=t_half, Y_half=Y_half, ode_result=result)


@dataclass
class CorrectionIteration:
    """One row of the Newton correction convergence history."""

    iteration: int
    x0: float
    z0: float
    ydot0: float
    half_period: float
    xdot_residual: float
    zdot_residual: float
    correction_norm: float


@dataclass
class CorrectionResult:
    converged: bool
    state0: np.ndarray
    half_period: float
    iterations: list = field(default_factory=list)


def differential_correct_halo(
    state0_guess,
    mu: float,
    tol: float = 1e-12,
    max_iter: int = 25,
    t_max: float = 6.0,
    rtol=DEFAULT_RTOL,
    atol=DEFAULT_ATOL,
) -> CorrectionResult:
    """Newton-correct a symmetric halo seed to xdot(T/2)=zdot(T/2)=0.

    Free variables: x0, ydot0 (z0 held fixed). See module docstring.

    Raises DifferentialCorrectionError if the iteration does not
    converge within `max_iter` steps or the 2x2 correction block
    becomes singular/ill-conditioned (condition number > 1e12).
    """
    state = np.array(state0_guess, dtype=float)
    z0 = state[2]  # held fixed throughout
    history: list[CorrectionIteration] = []

    for it in range(max_iter):
        half = propagate_half_period(state, mu, t_max, rtol, atol)
        Y_half = half.Y_half
        s_half = Y_half[0:6]
        Phi_half = Y_half[6:].reshape(6, 6)

        xdot_res = s_half[3]
        zdot_res = s_half[5]
        residual_norm = float(np.hypot(xdot_res, zdot_res))

        history.append(
            CorrectionIteration(
                iteration=it,
                x0=float(state[0]),
                z0=float(state[2]),
                ydot0=float(state[4]),
                half_period=float(half.t_half),
                xdot_residual=float(xdot_res),
                zdot_residual=float(zdot_res),
                correction_norm=0.0,  # filled in below once computed
            )
        )

        if residual_norm < tol:
            return CorrectionResult(converged=True, state0=state.copy(), half_period=half.t_half, iterations=history)

        # Sensitivity of [xdot(T/2), zdot(T/2)] to [x0, ydot0], holding
        # z0 fixed. We must also account for the fact that T/2 itself
        # depends on x0, ydot0 (the event time moves). Use the standard
        # halo correction formula that accounts for this via the
        # instantaneous state derivative at the crossing:
        #
        #   d(target)/d(free) = Phi[rows, cols] - (sdot(T/2)[rows] * dTdfree) / ydot(T/2)
        #
        # where dTdfree solves y(T/2) held at 0 (Phi[1,cols] * d(free) + ydot(T/2)*dT = 0).
        sdot_half = cr3bp_rhs(half.t_half, s_half, mu)
        ydot_half = s_half[4]
        if abs(ydot_half) < 1e-12:
            raise DifferentialCorrectionError(
                f"ydot(T/2)={ydot_half!r} too small; cannot solve for event-time sensitivity"
            )

        free_cols = [0, 4]  # x0, ydot0
        target_rows = [3, 5]  # xdot, zdot

        # dT/d(free) from holding y(T/2) = 0: Phi[1, free] . d(free) + ydot_half * dT = 0
        Phi_y_free = Phi_half[1, free_cols]  # shape (2,)

        # d[s(row)](T/2)/d(free) = Phi[row,free] + sdot(T/2)[row] * dT/d(free)
        # (total-derivative chain rule through the implicitly-defined
        # event time T/2(free); see DESIGN.md's M3 section for the
        # full derivation).
        M = np.zeros((2, 2))
        for i, row in enumerate(target_rows):
            dT_dfree = -Phi_y_free / ydot_half  # shape (2,), for this row's linear system
            M[i, :] = Phi_half[row, free_cols] + sdot_half[row] * dT_dfree

        rhs = -np.array([xdot_res, zdot_res])

        cond = np.linalg.cond(M)
        if not np.isfinite(cond) or cond > 1e12:
            raise DifferentialCorrectionError(
                f"correction matrix is singular/ill-conditioned (cond={cond!r}) at iteration {it}"
            )

        delta, *_ = np.linalg.lstsq(M, rhs, rcond=None)
        history[-1].correction_norm = float(np.linalg.norm(delta))

        # Damped Newton step (backtracking line search): the full
        # Newton step is not guaranteed to reduce the residual far from
        # convergence (halo dynamics near L2 are locally hyperbolic in
        # the in-plane directions, so sensitivities can be large even
        # for a well-conditioned correction matrix). Halve the step
        # until the half-period residual actually improves, or give up
        # after a bounded number of halvings.
        accepted = False
        for backtrack in range(12):
            damping = 0.5**backtrack
            trial = state.copy()
            trial[0] = state[0] + damping * delta[0]
            trial[4] = state[4] + damping * delta[1]
            trial[2] = z0
            try:
                trial_half = propagate_half_period(trial, mu, t_max, rtol, atol)
            except DifferentialCorrectionError:
                continue  # this trial step leaves the region with a valid crossing; shrink further
            trial_s = trial_half.Y_half[0:6]
            trial_residual = float(np.hypot(trial_s[3], trial_s[5]))
            if trial_residual < residual_norm:
                state = trial
                accepted = True
                break

        if not accepted:
            raise DifferentialCorrectionError(
                f"Newton step did not improve the residual after backtracking at iteration {it} "
                f"(residual norm={residual_norm!r}); seed may be too far from a periodic orbit"
            )

    raise DifferentialCorrectionError(
        f"differential correction did not converge within {max_iter} iterations "
        f"(final residual norm={residual_norm!r}, tol={tol!r})"
    )
