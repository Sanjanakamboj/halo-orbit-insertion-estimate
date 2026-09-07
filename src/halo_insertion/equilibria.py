"""Collinear libration point (L1, L2, L3) solver for the Earth-Moon CR3BP.

Each point is found as a root of dOmega/dx(x, 0, 0) = 0 using
scipy.optimize.brentq on a physically motivated bracket, NOT hardcoded
to the M1 hand value. See DESIGN.md Section 6 for the derivation this
implements numerically and generalizes to L1/L3.
"""

from __future__ import annotations

from dataclasses import dataclass

from scipy.optimize import brentq

from .cr3bp import effective_potential_gradient

__all__ = ["EquilibriumPoint", "collinear_x_residual", "find_L1", "find_L2", "find_L3"]


@dataclass(frozen=True)
class EquilibriumPoint:
    """Result of a collinear equilibrium-point solve.

    Attributes
    ----------
    name : str
        "L1", "L2", or "L3".
    x : float
        Nondimensional barycentric x-coordinate (y=z=0).
    residual : float
        |dOmega/dx| evaluated at the solution — should be ~0.
    distance_from_earth_nondim : float
        |x - (-mu)|, nondimensional.
    distance_from_moon_nondim : float
        |x - (1-mu)|, nondimensional.
    """

    name: str
    x: float
    residual: float
    distance_from_earth_nondim: float
    distance_from_moon_nondim: float


def collinear_x_residual(x: float, mu: float) -> float:
    """dOmega/dx(x, 0, 0) — the scalar equation collinear equilibria satisfy."""
    dOdx, _, _ = effective_potential_gradient(x, 0.0, 0.0, mu)
    return dOdx


def _solve(name: str, mu: float, bracket: tuple[float, float]) -> EquilibriumPoint:
    a, b = bracket
    x = brentq(collinear_x_residual, a, b, args=(mu,), xtol=1e-14, rtol=1e-14, maxiter=200)
    residual = abs(collinear_x_residual(x, mu))
    return EquilibriumPoint(
        name=name,
        x=x,
        residual=residual,
        distance_from_earth_nondim=abs(x - (-mu)),
        distance_from_moon_nondim=abs(x - (1.0 - mu)),
    )


def find_L1(mu: float, eps: float = 1e-6) -> EquilibriumPoint:
    """L1: between Earth and Moon, x in (-mu + eps, 1 - mu - eps)."""
    return _solve("L1", mu, (-mu + eps, 1.0 - mu - eps))


def find_L2(mu: float, eps: float = 1e-6, x_max: float = 3.0) -> EquilibriumPoint:
    """L2: beyond the Moon, x in (1 - mu + eps, x_max)."""
    return _solve("L2", mu, (1.0 - mu + eps, x_max))


def find_L3(mu: float, eps: float = 1e-6, x_min: float = -3.0) -> EquilibriumPoint:
    """L3: beyond Earth on the far side from the Moon, x in (x_min, -mu - eps)."""
    return _solve("L3", mu, (x_min, -mu - eps))
