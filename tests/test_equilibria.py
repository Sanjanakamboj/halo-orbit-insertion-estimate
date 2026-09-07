"""Tests for the L1/L2/L3 collinear equilibrium-point solver.

Includes the M1-vs-M2 cross-check: L2 found by the general brentq
solver here must agree with the M1 hand-derived Newton's-method value
in DESIGN.md Section 6.
"""

import pytest

from halo_insertion import constants as c
from halo_insertion.cr3bp import cr3bp_rhs
from halo_insertion.equilibria import find_L1, find_L2, find_L3

MU = c.MU

# M1 hand value (DESIGN.md Section 6), reproduced with an independent
# Newton's-method solve during pre-flight verification for M2.
M1_L2_X = 1.1556821589


def test_l1_residual_near_zero():
    l1 = find_L1(MU)
    assert l1.residual < 1e-12


def test_l2_residual_near_zero():
    l2 = find_L2(MU)
    assert l2.residual < 1e-12


def test_l3_residual_near_zero():
    l3 = find_L3(MU)
    assert l3.residual < 1e-12


def test_l1_between_earth_and_moon():
    l1 = find_L1(MU)
    assert -MU < l1.x < (1.0 - MU)


def test_l2_beyond_moon():
    l2 = find_L2(MU)
    assert l2.x > (1.0 - MU)


def test_l3_beyond_earth_opposite_moon():
    l3 = find_L3(MU)
    assert l3.x < -MU
    assert l3.x < 0


def test_l2_agrees_with_m1_hand_value():
    l2 = find_L2(MU)
    assert l2.x == pytest.approx(M1_L2_X, abs=1e-9)


def test_l2_dimensional_distance_from_moon():
    l2 = find_L2(MU)
    dist_km = l2.distance_from_moon_nondim * c.DU_KM
    assert dist_km == pytest.approx(64514.906, abs=1.0)
    # Physical sanity range from DESIGN.md Section 6.
    assert 60000.0 < dist_km < 70000.0


def test_l3_farther_from_barycenter_than_l1_l2_offsets():
    # Coarse ordering sanity check: L3 lies roughly one DU on the far
    # side of Earth, well beyond L1/L2's offsets from the Moon.
    l1 = find_L1(MU)
    l2 = find_L2(MU)
    l3 = find_L3(MU)
    assert l3.x < l1.x < l2.x


@pytest.mark.parametrize("finder", [find_L1, find_L2, find_L3])
def test_equilibrium_state_has_zero_acceleration(finder):
    eq = finder(MU)
    state = [eq.x, 0.0, 0.0, 0.0, 0.0, 0.0]
    rhs = cr3bp_rhs(0.0, state, MU)
    # rhs = [xdot, ydot, zdot, xddot, yddot, zddot]; all should vanish.
    assert all(abs(v) < 1e-10 for v in rhs)


@pytest.mark.parametrize("finder", [find_L1, find_L2, find_L3])
def test_equilibrium_yz_acceleration_exactly_zero_by_symmetry(finder):
    eq = finder(MU)
    state = [eq.x, 0.0, 0.0, 0.0, 0.0, 0.0]
    rhs = cr3bp_rhs(0.0, state, MU)
    assert rhs[4] == pytest.approx(0.0, abs=1e-14)  # yddot
    assert rhs[5] == pytest.approx(0.0, abs=1e-14)  # zddot
