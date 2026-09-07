"""Tests for the M5 backward-propagated CR3BP arrival-arc module."""

import numpy as np
import pytest

from halo_insertion import constants as c
from halo_insertion.arrival_arc import (
    EARTH_DISTANCE_THRESHOLD_DU,
    MOON_GUARD_KM,
    propagate_backward_arc,
    round_trip_check,
)
from halo_insertion.equilibria import find_L2

MU = c.MU


def test_backward_propagation_returns_finite_states():
    l2 = find_L2(MU)
    r_h = np.array([l2.x - 0.02, 0.05, 0.01])
    v_arr = np.array([0.1, 0.05, -0.01])
    arc = propagate_backward_arc(r_h, v_arr, MU, 10.0)
    assert arc.success
    assert np.all(np.isfinite(arc.final_state))
    assert np.isfinite(arc.min_dist_earth_du)
    assert np.isfinite(arc.min_dist_moon_du)


def test_arrival_jacobi_matches_direct_formula():
    from halo_insertion.arrival import jacobi_from_state

    l2 = find_L2(MU)
    r_h = np.array([l2.x - 0.02, 0.05, 0.01])
    v_arr = np.array([0.1, 0.05, -0.01])
    arc = propagate_backward_arc(r_h, v_arr, MU, 5.0)
    C_direct = jacobi_from_state(r_h, v_arr, MU)
    assert arc.C_arr == pytest.approx(C_direct, abs=1e-13)


def test_backward_forward_round_trip_closes():
    l2 = find_L2(MU)
    r_h = np.array([l2.x - 0.02, 0.05, 0.01])
    v_arr = np.array([0.1, 0.05, -0.01])
    arc = propagate_backward_arc(r_h, v_arr, MU, 8.0)
    assert arc.success
    recovered, pos_err, vel_err = round_trip_check(arc.s_arr, arc.final_state, arc.t_back_nondim, MU)
    assert pos_err < 1e-8
    assert vel_err < 1e-8


def test_zero_duration_propagation_is_identity():
    l2 = find_L2(MU)
    r_h = np.array([l2.x - 0.02, 0.05, 0.01])
    v_arr = np.array([0.1, 0.05, -0.01])
    s_arr = np.concatenate([r_h, v_arr])
    # A vanishingly small (not exactly zero, to avoid a degenerate
    # solve_ivp t_span) backward duration should reproduce the initial
    # state almost exactly.
    arc = propagate_backward_arc(r_h, v_arr, MU, 1e-6)
    assert np.allclose(arc.final_state, s_arr, atol=1e-6)


def test_moon_earth_distance_hand_constructed():
    # A state placed a small, non-singular offset from the Moon's
    # position: min distance to the Moon must be small (a few thousand
    # km) and distance to Earth must be close to 1 DU (Earth-Moon
    # separation), matching simple geometry.
    near_moon = np.array([1.0 - MU + 0.005, 0.0, 0.0])
    v = np.array([0.0, 0.02, 0.0])
    arc = propagate_backward_arc(near_moon, v, MU, 0.5, n_eval=200)
    assert arc.success
    assert arc.min_dist_moon_du * c.DU_KM < 5000.0
    assert arc.final_dist_earth_du == pytest.approx(1.0, abs=0.05)


def test_earthward_criterion_behaves_correctly():
    # A state whose distance from Earth never drops below the threshold
    # (e.g., an equilibrium state that stays near L2) must NOT be
    # accepted as earthward.
    l2 = find_L2(MU)
    r_h = np.array([l2.x, 0.0, 0.0])
    v_arr = np.array([0.001, 0.0, 0.0])  # tiny, stays near L2
    arc = propagate_backward_arc(r_h, v_arr, MU, 2.0)
    assert arc.min_dist_earth_du > EARTH_DISTANCE_THRESHOLD_DU
    assert not arc.earthward_ok
    assert not arc.accepted


def test_moon_guard_rejects_close_approach():
    # A trajectory that stays extremely close to the Moon (few tens of
    # km) must trip the Moon guard and be rejected.
    near_moon = np.array([1.0 - MU + 0.005, 0.0, 0.0])
    v = np.array([0.0, 0.02, 0.0])
    arc = propagate_backward_arc(near_moon, v, MU, 0.5, n_eval=200)
    assert arc.success
    assert arc.min_dist_moon_du * c.DU_KM < MOON_GUARD_KM
    assert not arc.moon_guard_ok
    assert not arc.accepted


def test_propagated_arrival_velocity_used_in_delta_v():
    from halo_insertion.arrival import delta_v_vector

    l2 = find_L2(MU)
    r_h = np.array([l2.x - 0.02, 0.05, 0.01])
    v_arr = np.array([0.1, 0.05, -0.01])
    v_halo = np.array([0.02, -0.03, 0.005])

    arc = propagate_backward_arc(r_h, v_arr, MU, 10.0)
    recovered, _, _ = round_trip_check(arc.s_arr, arc.final_state, arc.t_back_nondim, MU)
    v_arr_recovered = recovered[3:6]

    dv_recovered = delta_v_vector(v_halo, v_arr_recovered)
    dv_original = delta_v_vector(v_halo, v_arr)
    # The recovered arrival velocity should reproduce the original to
    # round-trip precision, so the two Delta_v computations must agree
    # tightly (not be independently unrelated).
    assert np.allclose(dv_recovered, dv_original, atol=1e-8)
