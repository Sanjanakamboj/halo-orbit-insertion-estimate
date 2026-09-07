"""Tests for the CR3BP effective potential, gradient, and equations of
motion — including an implementation-independent cross-check (B) and
structural symmetry checks.
"""

import math

import numpy as np
import pytest

from halo_insertion import constants as c
from halo_insertion.cr3bp import (
    CR3BPStateError,
    cr3bp_rhs,
    effective_potential,
    effective_potential_gradient,
    primary_distances,
)

MU = c.MU

# A handful of arbitrary nonsingular states, deliberately NOT any state
# already appearing in DESIGN.md, used for the independent cross-check.
ARBITRARY_STATES = [
    (0.5, 0.2, 0.1),
    (1.2, -0.05, 0.02),
    (-0.3, 0.4, -0.1),
    (0.85, 0.0, 0.15),
    (1.05, 0.03, 0.0),
]


def _independent_acceleration(x, y, z, xdot, ydot, zdot, mu):
    """Independently-coded (not shared code with cr3bp.py) CR3BP RHS,
    written directly from the vector form of the equations rather than
    reusing effective_potential_gradient. Used to catch:
      - wrong Coriolis signs
      - wrong Earth/Moon coordinate offsets
      - swapped mu / (1-mu)
      - missing centrifugal terms.
    """
    r1vec = np.array([x + mu, y, z])
    r2vec = np.array([x - (1 - mu), y, z])
    r1 = np.linalg.norm(r1vec)
    r2 = np.linalg.norm(r2vec)

    # Gravitational accelerations from each primary (pulls toward the
    # primary), plus centrifugal (from the x,y rotating-frame term) and
    # Coriolis (2*omega x v with omega = [0,0,1]) terms, built from
    # scratch rather than via Omega's gradient.
    grav = -(1 - mu) * r1vec / r1**3 - mu * r2vec / r2**3
    centrifugal = np.array([x, y, 0.0])
    coriolis = np.array([2 * ydot, -2 * xdot, 0.0])

    accel = grav + centrifugal + coriolis
    return accel


@pytest.mark.parametrize("x,y,z", ARBITRARY_STATES)
def test_potential_and_gradient_finite_away_from_primaries(x, y, z):
    omega = effective_potential(x, y, z, MU)
    grad = effective_potential_gradient(x, y, z, MU)
    assert math.isfinite(omega)
    assert all(math.isfinite(g) for g in grad)


@pytest.mark.parametrize("x,y,z", ARBITRARY_STATES)
def test_independent_rhs_cross_check(x, y, z):
    xdot, ydot, zdot = 0.05, -0.02, 0.01
    state = [x, y, z, xdot, ydot, zdot]

    rhs = cr3bp_rhs(0.0, state, MU)
    accel_production = rhs[3:6]

    accel_independent = _independent_acceleration(x, y, z, xdot, ydot, zdot, MU)

    assert np.allclose(accel_production, accel_independent, rtol=1e-12, atol=1e-13)


@pytest.mark.parametrize("x,y,z", ARBITRARY_STATES)
def test_analytic_gradient_matches_finite_difference(x, y, z):
    # Test-only finite-difference check of the analytic gradient (the
    # production RHS itself uses no finite differences).
    h = 1e-6
    domega_dx_fd = (
        effective_potential(x + h, y, z, MU) - effective_potential(x - h, y, z, MU)
    ) / (2 * h)
    domega_dy_fd = (
        effective_potential(x, y + h, z, MU) - effective_potential(x, y - h, z, MU)
    ) / (2 * h)
    domega_dz_fd = (
        effective_potential(x, y, z + h, MU) - effective_potential(x, y, z - h, MU)
    ) / (2 * h)

    dOdx, dOdy, dOdz = effective_potential_gradient(x, y, z, MU)
    assert dOdx == pytest.approx(domega_dx_fd, abs=1e-8)
    assert dOdy == pytest.approx(domega_dy_fd, abs=1e-8)
    assert dOdz == pytest.approx(domega_dz_fd, abs=1e-8)


def test_primary_distances_definition():
    # At the Moon's location itself, r2 should be 0 and r1 should be 1
    # (Earth-Moon separation is 1 DU by normalization).
    x_moon, y_moon, z_moon = c.MOON_POS_NONDIM
    r1, r2 = primary_distances(x_moon, y_moon, z_moon, MU)
    assert r1 == pytest.approx(1.0, abs=1e-12)
    assert r2 == pytest.approx(0.0, abs=1e-12)


def test_planar_motion_remains_planar():
    # Symmetry check: if z = zdot = 0, then zddot = dOmega/dz = 0, so a
    # state confined to the z=0 plane has no out-of-plane acceleration
    # and stays planar under the RHS.
    state = [0.8, 0.1, 0.0, 0.01, -0.02, 0.0]
    rhs = cr3bp_rhs(0.0, state, MU)
    assert rhs[2] == pytest.approx(0.0, abs=1e-14)  # zdot component of state deriv = z's own zdot=0
    assert rhs[5] == pytest.approx(0.0, abs=1e-14)  # zddot = 0


def test_z_reflection_symmetry_of_gradient():
    # Structural symmetry: Omega and dOmega/dx, dOmega/dy are even in z,
    # and dOmega/dz is odd in z (both r1, r2 depend on z^2 only).
    x, y, z = 0.9, 0.15, 0.2
    dOdx_pos, dOdy_pos, dOdz_pos = effective_potential_gradient(x, y, z, MU)
    dOdx_neg, dOdy_neg, dOdz_neg = effective_potential_gradient(x, y, -z, MU)

    assert dOdx_pos == pytest.approx(dOdx_neg, abs=1e-14)
    assert dOdy_pos == pytest.approx(dOdy_neg, abs=1e-14)
    assert dOdz_pos == pytest.approx(-dOdz_neg, abs=1e-14)


def test_y_reflection_symmetry_of_gradient():
    # Structural symmetry: Omega is even in y, dOmega/dx is even in y,
    # dOmega/dy is odd in y (r1, r2 depend on y^2 only).
    x, y, z = 0.9, 0.15, 0.05
    dOdx_pos, dOdy_pos, _ = effective_potential_gradient(x, y, z, MU)
    dOdx_neg, dOdy_neg, _ = effective_potential_gradient(x, -y, z, MU)

    assert dOdx_pos == pytest.approx(dOdx_neg, abs=1e-14)
    assert dOdy_pos == pytest.approx(-dOdy_neg, abs=1e-14)


def test_rejects_wrong_state_size():
    with pytest.raises(CR3BPStateError):
        cr3bp_rhs(0.0, [1.0, 0.0, 0.0], MU)


def test_rejects_non_finite_state():
    with pytest.raises(CR3BPStateError):
        cr3bp_rhs(0.0, [float("nan"), 0.0, 0.0, 0.0, 0.0, 0.0], MU)
    with pytest.raises(CR3BPStateError):
        cr3bp_rhs(0.0, [float("inf"), 0.0, 0.0, 0.0, 0.0, 0.0], MU)


def test_rejects_invalid_mu():
    state = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    with pytest.raises(CR3BPStateError):
        cr3bp_rhs(0.0, state, -0.1)
    with pytest.raises(CR3BPStateError):
        cr3bp_rhs(0.0, state, 1.5)
    with pytest.raises(CR3BPStateError):
        cr3bp_rhs(0.0, state, float("nan"))


def test_rejects_collision_with_earth():
    state = [-MU, 0.0, 0.0, 0.0, 0.0, 0.0]  # exactly at Earth's position
    with pytest.raises(CR3BPStateError):
        cr3bp_rhs(0.0, state, MU)


def test_rejects_collision_with_moon():
    state = [1.0 - MU, 0.0, 0.0, 0.0, 0.0, 0.0]  # exactly at Moon's position
    with pytest.raises(CR3BPStateError):
        cr3bp_rhs(0.0, state, MU)
