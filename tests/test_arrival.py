"""Tests for the M4 arrival-state models (arrival.py)."""

import numpy as np
import pytest

from halo_insertion import constants as c
from halo_insertion.arrival import (
    ArrivalStateError,
    arrival_speed_from_jacobi,
    arrival_velocity,
    delta_v_vector,
    direction_aligned_with_halo,
    direction_radial_from_earth,
    direction_rotated_about_z,
    earth_position,
    jacobi_from_state,
)
from halo_insertion.cr3bp import effective_potential
from halo_insertion.equilibria import find_L2
from halo_insertion.normalization import m_s_to_nondim_velocity, nondim_velocity_to_m_s

MU = c.MU


def test_direction_vectors_are_unit_norm():
    l2 = find_L2(MU)
    r_h = np.array([l2.x - 0.05, 0.08, 0.02])
    v_halo = np.array([0.05, -0.1, 0.02])

    d_aligned = direction_aligned_with_halo(v_halo)
    d_radial = direction_radial_from_earth(r_h, MU)
    d_rotated = direction_rotated_about_z(d_radial, np.radians(15.0))

    for d in (d_aligned, d_radial, d_rotated):
        assert np.linalg.norm(d) == pytest.approx(1.0, abs=1e-12)


def test_velocity_round_trip_m_s():
    for v_m_s in [0.0, 25.0, 109.57, -50.0]:
        v_nd = m_s_to_nondim_velocity(v_m_s)
        v_back = nondim_velocity_to_m_s(v_nd)
        assert v_back == pytest.approx(v_m_s, rel=1e-13, abs=1e-9)


def test_delta_v_zero_when_arrival_equals_halo():
    v_halo = np.array([0.1, -0.05, 0.02])
    dv = delta_v_vector(v_halo, v_halo)
    assert np.linalg.norm(dv) == pytest.approx(0.0, abs=1e-14)


def test_fixed_speed_aligned_case_matches_scalar_speed_difference():
    v_halo = np.array([0.1, 0.05, -0.02])
    v_halo_speed = np.linalg.norm(v_halo)
    direction = direction_aligned_with_halo(v_halo)
    arr_speed = 0.15
    v_arr = arrival_velocity(arr_speed, direction)
    dv = delta_v_vector(v_halo, v_arr)
    assert np.linalg.norm(dv) == pytest.approx(abs(v_halo_speed - arr_speed), abs=1e-12)


def test_jacobi_from_state_matches_direct_formula():
    r = np.array([1.1, 0.05, 0.02])
    v = np.array([0.1, -0.05, 0.03])
    C_direct = 2.0 * effective_potential(*r, MU) - float(np.dot(v, v))
    C_from_func = jacobi_from_state(r, v, MU)
    assert C_from_func == pytest.approx(C_direct, abs=1e-14)


def test_arrival_speed_from_jacobi_rejects_negative_v_squared():
    l2 = find_L2(MU)
    r_h = np.array([l2.x, 0.0, 0.0])
    omega = effective_potential(*r_h, MU)
    # Choose C_arr so large that 2*Omega - C_arr < 0.
    C_arr_too_large = 2.0 * omega + 10.0
    with pytest.raises(ArrivalStateError):
        arrival_speed_from_jacobi(r_h, MU, C_arr_too_large)


def test_arrival_speed_from_jacobi_matches_definition():
    l2 = find_L2(MU)
    r_h = np.array([l2.x - 0.05, 0.05, 0.01])
    omega = effective_potential(*r_h, MU)
    C_arr = omega  # arbitrary valid choice: 2*omega - omega = omega > 0
    speed = arrival_speed_from_jacobi(r_h, MU, C_arr)
    assert speed**2 == pytest.approx(2.0 * omega - C_arr, abs=1e-12)


def test_radial_direction_points_away_from_earth():
    r_h = np.array([1.1, 0.05, 0.0])
    d = direction_radial_from_earth(r_h, MU)
    rel = r_h - earth_position(MU)
    assert np.allclose(d, rel / np.linalg.norm(rel))


def test_rotation_by_zero_is_identity():
    base = np.array([1.0, 0.0, 0.0])
    rotated = direction_rotated_about_z(base, 0.0)
    assert np.allclose(rotated, base, atol=1e-14)


def test_rotation_by_90_degrees_hand_constructed():
    base = np.array([1.0, 0.0, 0.0])
    rotated = direction_rotated_about_z(base, np.pi / 2)
    assert np.allclose(rotated, [0.0, 1.0, 0.0], atol=1e-12)
