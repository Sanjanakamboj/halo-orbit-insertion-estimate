"""Generic CR3BP numerical propagation wrapper around scipy.integrate.solve_ivp.

This module is deliberately halo-agnostic: it knows nothing about halo
orbits, seeds, or differential correction. It exists purely to propagate
an arbitrary CR3BP state forward (or backward) in normalized time.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp

from .cr3bp import cr3bp_rhs

__all__ = ["PropagationSettings", "propagate"]

DEFAULT_METHOD = "DOP853"
DEFAULT_RTOL = 1e-11
DEFAULT_ATOL = 1e-12


@dataclass(frozen=True)
class PropagationSettings:
    """Integrator settings, kept explicit rather than buried in kwargs."""

    method: str = DEFAULT_METHOD
    rtol: float = DEFAULT_RTOL
    atol: float = DEFAULT_ATOL
    dense_output: bool = False


def propagate(
    state0,
    t_span: tuple[float, float],
    mu: float,
    settings: PropagationSettings | None = None,
    t_eval=None,
):
    """Propagate a CR3BP state from t_span[0] to t_span[1].

    Parameters
    ----------
    state0 : array-like, shape (6,)
        Initial [x, y, z, xdot, ydot, zdot] in normalized units.
    t_span : (float, float)
        Start and end nondimensional times.
    mu : float
        Earth-Moon mass ratio.
    settings : PropagationSettings, optional
        Integrator method/tolerances. Defaults to DOP853 with
        rtol=1e-11, atol=1e-12 (M2 default; see DESIGN.md Section 9).
    t_eval : array-like, optional
        Times at which to store the solution.

    Returns
    -------
    scipy.integrate.OdeResult
        The raw solve_ivp result, exposed unmodified. `result.success`
        must be checked by the caller — solver failures are not hidden
        or silently swallowed here.
    """
    if settings is None:
        settings = PropagationSettings()

    state0 = np.asarray(state0, dtype=float)

    result = solve_ivp(
        fun=cr3bp_rhs,
        t_span=t_span,
        y0=state0,
        method=settings.method,
        rtol=settings.rtol,
        atol=settings.atol,
        dense_output=settings.dense_output,
        t_eval=t_eval,
        args=(mu,),
    )
    return result
