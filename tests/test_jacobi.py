"""Tests for the Jacobi constant diagnostic."""

import pytest

from halo_insertion import constants as c
from halo_insertion.cr3bp import effective_potential, jacobi_constant
from halo_insertion.equilibria import find_L2

MU = c.MU


def test_equilibrium_state_gives_C_equals_2_omega():
    l2 = find_L2(MU)
    state = [l2.x, 0.0, 0.0, 0.0, 0.0, 0.0]
    omega = effective_potential(l2.x, 0.0, 0.0, MU)
    C = jacobi_constant(state, MU)
    assert C == pytest.approx(2.0 * omega, abs=1e-13)


def test_increased_speed_lowers_C_at_fixed_position():
    x, y, z = 0.9, 0.1, 0.05
    state_slow = [x, y, z, 0.01, 0.0, 0.0]
    state_fast = [x, y, z, 0.05, 0.0, 0.0]

    C_slow = jacobi_constant(state_slow, MU)
    C_fast = jacobi_constant(state_fast, MU)

    assert C_fast < C_slow


def test_zero_velocity_gives_C_equals_2_omega_generic_point():
    x, y, z = 0.7, -0.2, 0.1
    state = [x, y, z, 0.0, 0.0, 0.0]
    omega = effective_potential(x, y, z, MU)
    C = jacobi_constant(state, MU)
    assert C == pytest.approx(2.0 * omega, abs=1e-13)


def test_jacobi_symmetric_in_velocity_direction():
    # C depends on speed^2, not direction, so reversing velocity leaves
    # C unchanged at fixed position — a cheap structural check.
    x, y, z = 0.8, 0.05, -0.02
    v = [0.03, -0.01, 0.02]
    C_fwd = jacobi_constant([x, y, z, *v], MU)
    C_rev = jacobi_constant([x, y, z, *(-vi for vi in v)], MU)
    assert C_fwd == pytest.approx(C_rev, abs=1e-14)
