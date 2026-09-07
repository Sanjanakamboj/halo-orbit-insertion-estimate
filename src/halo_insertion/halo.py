"""Top-level orchestration: build and independently validate a
numerically corrected Earth-Moon L2 halo orbit.

This module ties together `halo_seed`, `differential_correction`,
`variational`, and the M2 generic `propagation`/`cr3bp` machinery into
one result object. It does not introduce any new dynamics.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .cr3bp import jacobi_constant
from .differential_correction import CorrectionResult, differential_correct_halo
from .halo_seed import LinearHaloSeed, linear_halo_seed
from .normalization import nondim_time_to_days
from .propagation import PropagationSettings, propagate
from .variational import propagate_with_stm

__all__ = ["HaloOrbitResult", "build_l2_halo", "full_period_closure", "monodromy_matrix"]


@dataclass
class HaloOrbitResult:
    mu: float
    seed: LinearHaloSeed
    correction: CorrectionResult
    state0: np.ndarray
    half_period: float
    period: float
    period_days: float
    jacobi_C0: float


def build_l2_halo(x_l2: float, mu: float, Az: float, z_sign: float = 1.0, tol: float = 1e-11) -> HaloOrbitResult:
    """Construct the linear seed, differentially correct it, and return
    the corrected periodic orbit's initial condition and period.
    """
    seed = linear_halo_seed(x_l2, mu, Az, z_sign=z_sign)
    correction = differential_correct_halo(seed.state0, mu, tol=tol)
    state0 = correction.state0
    half_period = correction.half_period
    period = 2.0 * half_period
    C0 = jacobi_constant(state0, mu)

    return HaloOrbitResult(
        mu=mu,
        seed=seed,
        correction=correction,
        state0=state0,
        half_period=half_period,
        period=period,
        period_days=nondim_time_to_days(period),
        jacobi_C0=C0,
    )


def full_period_closure(state0, period: float, mu: float, settings: PropagationSettings | None = None):
    """Propagate one full period with the GENERIC M2 propagator
    (independent of the differential-correction machinery) and return
    (final_state, closure = final_state - state0, ode_result).
    """
    state0 = np.asarray(state0, dtype=float)
    result = propagate(state0, (0.0, period), mu, settings=settings)
    if not result.success:
        raise RuntimeError(f"full-period propagation failed: {result.message}")
    final_state = result.y[:, -1]
    closure = final_state - state0
    return final_state, closure, result


def monodromy_matrix(state0, period: float, mu: float, settings: PropagationSettings | None = None) -> np.ndarray:
    """Propagate the STM through one full period; M = Phi(T)."""
    stm_result = propagate_with_stm(state0, (0.0, period), mu, settings=settings)
    if not stm_result.ode_result.success:
        raise RuntimeError(f"monodromy STM propagation failed: {stm_result.ode_result.message}")
    return stm_result.Phi_final
