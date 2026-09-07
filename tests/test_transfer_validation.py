"""Tests for the M5 transfer-arrival validation: the selected arc's
consistency with M4, and no regression to M1-M4."""

import json

import numpy as np
import pytest

from halo_insertion import constants as c
from halo_insertion.arrival import delta_v_vector
from halo_insertion.arrival_arc import propagate_backward_arc, round_trip_check
from halo_insertion.cr3bp import jacobi_constant
from halo_insertion.equilibria import find_L2
from halo_insertion.insertion import evaluate_insertion_at_phase, halo_state_at_phase
from halo_insertion.propagation import PropagationSettings

M4_TAU = 0.24413043101153073
M4_DV_M_S = 109.57233088725141
DC_BASELINE = 0.01


@pytest.fixture(scope="module")
def m3_halo():
    with open("results/m3_halo_summary.json") as f:
        summary = json.load(f)
    s0 = summary["corrected_initial_state"]
    state0 = [s0["x0"], s0["y0"], s0["z0"], s0["xdot0"], s0["ydot0"], s0["zdot0"]]
    period = summary["period_nondim"]
    return state0, period, summary


def test_direct_vector_subtraction_matches_production_result(m3_halo):
    state0, period, _ = m3_halo
    cand = evaluate_insertion_at_phase(state0, period, c.MU, M4_TAU, DC_BASELINE)
    arc = propagate_backward_arc(cand.r_h, cand.v_arr, c.MU, 10.0, tau=M4_TAU, dC=DC_BASELINE)
    assert arc.accepted
    recovered, _, _ = round_trip_check(arc.s_arr, arc.final_state, arc.t_back_nondim, c.MU)
    dv_vec = delta_v_vector(cand.v_halo, recovered[3:6])
    dv_m_s = float(np.linalg.norm(dv_vec)) * c.VSTAR_M_S
    assert dv_m_s == pytest.approx(M4_DV_M_S, abs=1e-4)


def test_tighter_tolerance_reproduces_selected_arc(m3_halo):
    state0, period, _ = m3_halo
    cand = evaluate_insertion_at_phase(state0, period, c.MU, M4_TAU, DC_BASELINE)
    default_arc = propagate_backward_arc(cand.r_h, cand.v_arr, c.MU, 10.0, tau=M4_TAU, dC=DC_BASELINE)
    tight_settings = PropagationSettings(method="DOP853", rtol=1e-13, atol=1e-14)
    tight_arc = propagate_backward_arc(cand.r_h, cand.v_arr, c.MU, 10.0, tau=M4_TAU, dC=DC_BASELINE, settings=tight_settings)
    assert np.linalg.norm(default_arc.final_state - tight_arc.final_state) < 1e-6


def test_alternative_integration_method_agrees(m3_halo):
    state0, period, _ = m3_halo
    cand = evaluate_insertion_at_phase(state0, period, c.MU, M4_TAU, DC_BASELINE)
    dop853_arc = propagate_backward_arc(cand.r_h, cand.v_arr, c.MU, 10.0, tau=M4_TAU, dC=DC_BASELINE)
    rk45_settings = PropagationSettings(method="RK45", rtol=1e-12, atol=1e-13)
    rk45_arc = propagate_backward_arc(cand.r_h, cand.v_arr, c.MU, 10.0, tau=M4_TAU, dC=DC_BASELINE, settings=rk45_settings)
    assert np.linalg.norm(dop853_arc.final_state - rk45_arc.final_state) < 1e-6


def test_m4_local_result_still_reproduces_exactly(m3_halo):
    state0, period, _ = m3_halo
    cand = evaluate_insertion_at_phase(state0, period, c.MU, M4_TAU, DC_BASELINE)
    assert cand.delta_v_m_s == pytest.approx(M4_DV_M_S, abs=1e-6)


def test_m3_halo_state_at_selected_phase_independently_verified(m3_halo):
    state0, period, summary = m3_halo
    state_recompute = halo_state_at_phase(state0, period, c.MU, M4_TAU)
    C_recompute = jacobi_constant(state_recompute, c.MU)
    assert C_recompute == pytest.approx(summary["jacobi_C0"], abs=1e-9)


def test_no_m1_m4_regression():
    assert c.MU == pytest.approx(0.0121505839, abs=1e-9)
    l2 = find_L2(c.MU)
    assert l2.x == pytest.approx(1.1556821589, abs=1e-9)

    with open("results/m4_insertion_summary.json") as f:
        m4 = json.load(f)
    assert m4["selected_insertion"]["tau"] == pytest.approx(M4_TAU, abs=1e-9)
    assert m4["selected_insertion"]["delta_v_m_s"] == pytest.approx(M4_DV_M_S, abs=1e-6)
