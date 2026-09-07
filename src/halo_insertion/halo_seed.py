"""Linearized (first-order) Earth-Moon L2 halo-orbit initial-condition
seed, for use as the starting guess for differential correction.

Provenance / method
--------------------
This is the **linear variational solution of the CR3BP linearized about
the L2 collinear point**, NOT the third-order Richardson (1980)
approximation. It is the standard first step described in halo-orbit
mission-design references (e.g. Koon, Lo, Marsden & Ross, *Dynamical
Systems, the Three-Body Problem, and Space Mission Design*, 2011,
Ch. 6): linearize the CR3BP equations about the equilibrium point,
solve the resulting linear ODEs in closed form, and use the oscillatory
(center-manifold) solution as a seed for a nonlinear differential
corrector. This module explicitly does NOT claim third-order accuracy
or the Richardson amplitude-frequency (Az-vs-period) correction — the
differential corrector in `differential_correction.py` is what turns
this approximate seed into a genuine numerically periodic orbit.

Derivation summary (documented fully in DESIGN.md's M3 section)
------------------------------------------------------------------
Linearizing the CR3BP equations about (x_L2, 0, 0) (where the
cross-derivatives Uxy = Uxz = Uyz = 0 by the equilibrium's symmetry)
gives two decoupled linear subsystems:

  in-plane (x,y):  ddot{dx} - 2*ddot{dy} = Uxx*dx,  ddot{dy} + 2*dot{dx} = Uyy*dy
  out-of-plane z:  ddot{dz} = Uzz*dz

The out-of-plane equation is simple harmonic motion with frequency
`wz = sqrt(-Uzz)`. The in-plane system's characteristic equation
`beta^2 + (4 - Uxx - Uyy)*beta + Uxx*Uyy = 0` (with beta = lambda^2) has
one negative root `beta < 0`, giving an oscillatory in-plane mode at
frequency `wxy = sqrt(-beta)`. Assuming trial solutions
`dx = Ax*cos(wxy*t)`, `dy = Ay*sin(wxy*t)` and substituting into the
in-plane equations gives `Ay = kappa * Ax` with
`kappa = -(Uxx + wxy^2) / (2*wxy)`.

At t=0 this trial solution gives exactly the symmetric halo initial
condition form required by the differential corrector (Section 7 of
this project's M3 documentation):

    x(0)    = x_L2 + Ax          (free/offset)
    y(0)    = 0                  (exact, by construction)
    z(0)    = Az                 (free amplitude)
    xdot(0) = 0                  (exact, by construction)
    ydot(0) = kappa * Ax * wxy   (from the linear relation)
    zdot(0) = 0                  (exact, for a cosine z-oscillation)

At L2, `wz` and `wxy` differ by only a few percent (verified
numerically for the mission's `Az` scale — see DESIGN.md), so this
first-order seed's residual 3-D coupling error is small enough for
Newton's method (using the STM) to correct.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from .variational import omega_hessian

__all__ = ["LinearHaloSeed", "linear_halo_seed"]


@dataclass(frozen=True)
class LinearHaloSeed:
    """A linearized halo-orbit seed and its supporting linear-theory quantities."""

    x_l2: float
    Uxx: float
    Uyy: float
    Uzz: float
    wxy: float
    wz: float
    kappa: float
    Ax: float
    Az: float
    state0: tuple  # (x, y, z, xdot, ydot, zdot)


def linear_halo_seed(x_l2: float, mu: float, Az: float, z_sign: float = 1.0) -> LinearHaloSeed:
    """Construct a linearized L2 halo seed with out-of-plane amplitude Az.

    Parameters
    ----------
    x_l2 : float
        Nondimensional L2 x-coordinate (from `equilibria.find_L2`).
    mu : float
        Earth-Moon mass ratio.
    Az : float
        Desired out-of-plane (z) amplitude, nondimensional (DU). Must be
        nonzero and positive; use `z_sign` to select the family branch.
    z_sign : float
        +1.0 for the branch with z(0) > 0, -1.0 for z(0) < 0. This
        project uses the southern family (DESIGN.md Section 1); the
        sign convention for "southern" vs. z(0) sign is documented
        alongside the corrected result, since it depends on the
        orientation convention in use.

    Returns
    -------
    LinearHaloSeed
    """
    if Az <= 0:
        raise ValueError(f"Az must be positive, got {Az!r}")

    hess = omega_hessian(x_l2, 0.0, 0.0, mu)
    Uxx, Uyy, Uzz = hess[0, 0], hess[1, 1], hess[2, 2]

    if Uzz >= 0:
        raise ValueError(
            f"Uzz={Uzz!r} >= 0 at this equilibrium point; no oscillatory "
            "out-of-plane linear mode exists here (not a valid halo seed point)"
        )
    wz = math.sqrt(-Uzz)

    b = 4.0 - Uxx - Uyy
    disc = b**2 - 4.0 * Uxx * Uyy
    if disc < 0:
        raise ValueError("Negative discriminant in in-plane characteristic equation")
    sqrt_disc = math.sqrt(disc)
    beta1 = (-b + sqrt_disc) / 2.0
    beta2 = (-b - sqrt_disc) / 2.0
    beta_osc = min(beta1, beta2)  # the negative (oscillatory) root
    if beta_osc >= 0:
        raise ValueError(
            f"No oscillatory in-plane root found (beta1={beta1!r}, beta2={beta2!r}); "
            "this point does not support a center-manifold halo-like mode"
        )
    wxy = math.sqrt(-beta_osc)

    kappa = -(Uxx + wxy**2) / (2.0 * wxy)

    # Choose the in-plane amplitude Ax from a fixed, moderate ratio to
    # Az (typical Az/Ax scale for Earth-Moon L2 halo family members is
    # O(1); the differential corrector does not depend on this choice
    # being exact, only "in the right regime" for Newton's method to
    # converge). Sign chosen so x0 sits on the Earth-facing side of L2
    # ("inward" halo lobe), matching typical Earth-Moon L2 halo
    # renderings; this is a seed convenience, not a physical constraint.
    Ax = 0.9 * Az

    # NOTE on sign: the trial solution dx(t) = Ax*cos(wxy*t), dy(t) =
    # kappa*Ax*sin(wxy*t) used to derive `kappa` above assumes dx(0) =
    # +Ax. This seed instead places x0 on the Earth-facing side of L2,
    # i.e. dx(0) = x0 - x_l2 = -Ax, so the corresponding trial solution
    # here is dx(t) = -Ax*cos(wxy*t), which flips the sign of the
    # ydot0 relation to -kappa*Ax*wxy. (Verified empirically: propagating
    # this seed under the full nonlinear CR3BP produces a small,
    # near-periodic first return, confirming the correct branch;
    # +kappa*Ax*wxy was checked and diverges rapidly — the two branches
    # are not physically equivalent because dx0's sign is fixed by the
    # x0 = x_l2 - Ax convention below, not free.)
    x0 = x_l2 - Ax
    y0 = 0.0
    z0 = z_sign * Az
    xdot0 = 0.0
    ydot0 = -kappa * Ax * wxy
    zdot0 = 0.0

    return LinearHaloSeed(
        x_l2=x_l2,
        Uxx=Uxx,
        Uyy=Uyy,
        Uzz=Uzz,
        wxy=wxy,
        wz=wz,
        kappa=kappa,
        Ax=Ax,
        Az=Az,
        state0=(x0, y0, z0, xdot0, ydot0, zdot0),
    )
