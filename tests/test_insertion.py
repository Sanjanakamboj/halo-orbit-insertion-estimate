"""Tests for the M4 insertion phase-sweep/refinement machinery
(insertion.py), using the frozen M3 halo orbit."""

import json

import numpy as np
import pytest

from halo_insertion import constants as c
from halo_insertion.arrival import delta_v_vector, direction_rotated_about_z
from halo_insertion.cr3bp import jacobi_constant
from halo_insertion.insertion import evaluate_insertion_at_phase, halo_state_at_phase, phase_sweep, refine_best_phase


@pytest.fixture(scope="module")
def m3_halo():
    with open("results/m3_halo_summary.json") as f:
        summary = json.load(f)
    s0 = summary["corrected_initial_state"]
    state0 = [s0["x0"], s0["y0"], s0["z0"], s0["xdot0"], s0["ydot0"], s0["zdot0"]]
    period = summary["period_nondim"]
    return state0, period, summary


DC_BASELINE = 0.01


def test_phase_zero_and_one_agree(m3_halo):
    state0, period, _ = m3_halo
    s_tau0 = halo_state_at_phase(state0, period, c.MU, 0.0)
    s_tau1 = halo_state_at_phase(state0, period, c.MU, 1.0)
    assert np.allclose(s_tau0, s_tau1, atol=1e-9)


def test_dense_sweep_returns_finite_values(m3_halo):
    state0, period, _ = m3_halo
    candidates = phase_sweep(state0, period, c.MU, DC_BASELINE, n_points=100)
    for cand in candidates:
        assert cand.valid, cand.reason
        assert np.isfinite(cand.delta_v_m_s)
        assert np.isfinite(cand.delta_v_nondim)
        assert cand.delta_v_m_s > 0


def test_optimizer_does_not_worsen_grid_minimum(m3_halo):
    state0, period, _ = m3_halo
    candidates = phase_sweep(state0, period, c.MU, DC_BASELINE, n_points=200)
    dvs = np.array([cand.delta_v_m_s for cand in candidates])
    taus = np.array([cand.tau for cand in candidates])
    i = int(np.argmin(dvs))
    grid_min = dvs[i]

    refined = refine_best_phase(state0, period, c.MU, DC_BASELINE, (taus[i] - 0.01, taus[i] + 0.01))
    assert refined.delta_v_m_s <= grid_min + 1e-6


def test_direct_selected_point_recomputation_matches(m3_halo):
    state0, period, _ = m3_halo
    candidates = phase_sweep(state0, period, c.MU, DC_BASELINE, n_points=200)
    dvs = np.array([cand.delta_v_m_s for cand in candidates])
    taus = np.array([cand.tau for cand in candidates])
    i = int(np.argmin(dvs))
    refined = refine_best_phase(state0, period, c.MU, DC_BASELINE, (taus[i] - 0.01, taus[i] + 0.01))

    # Direct recomputation at the exact selected tau.
    recomputed = evaluate_insertion_at_phase(state0, period, c.MU, refined.tau, DC_BASELINE)
    assert recomputed.delta_v_m_s == pytest.approx(refined.delta_v_m_s, abs=1e-6)
    assert np.allclose(recomputed.r_h, refined.r_h, atol=1e-12)


def test_impulsive_delta_C_relation(m3_halo):
    state0, period, _ = m3_halo
    cand = evaluate_insertion_at_phase(state0, period, c.MU, 0.25, DC_BASELINE)
    # delta_C = C_halo - C_arr, and by construction C_arr = C_halo - dC_baseline,
    # so delta_C should equal dC_baseline exactly (up to numerical precision).
    assert cand.delta_C == pytest.approx(DC_BASELINE, abs=1e-10)


def test_arrival_speed_sensitivity_monotonic_in_simple_case(m3_halo):
    state0, period, _ = m3_halo
    # Larger dC (bigger Jacobi offset) at a FIXED phase should generally
    # increase arrival speed magnitude (v_arr^2 = 2*Omega - C_arr, and
    # increasing dC decreases C_arr, which increases v_arr^2).
    tau = 0.25
    speeds = []
    for dC in [0.005, 0.01, 0.02, 0.05]:
        cand = evaluate_insertion_at_phase(state0, period, c.MU, tau, dC)
        speeds.append(np.linalg.norm(cand.v_arr))
    for i in range(len(speeds) - 1):
        assert speeds[i + 1] > speeds[i]


def test_direction_perturbation_hand_constructed_case():
    # Fully hand-constructible, halo-independent check: v_halo along +x,
    # v_arr at fixed speed along the baseline (+x) direction rotated by
    # theta about z. |Delta_v| must match the law of cosines exactly.
    v_halo = np.array([0.2, 0.0, 0.0])
    speed = 0.15
    base_dir = np.array([1.0, 0.0, 0.0])
    for theta_deg in [0, 10, 30, 90]:
        theta = np.radians(theta_deg)
        direction = direction_rotated_about_z(base_dir, theta)
        v_arr = speed * direction
        dv = delta_v_vector(v_halo, v_arr)
        dv_norm = np.linalg.norm(dv)
        expected = np.sqrt(0.2**2 + speed**2 - 2 * 0.2 * speed * np.cos(theta))
        assert dv_norm == pytest.approx(expected, abs=1e-12)


def test_m3_halo_state_period_jacobi_unchanged(m3_halo):
    state0, period, summary = m3_halo
    assert summary["mu"] == pytest.approx(c.MU, abs=1e-15)
    assert period == pytest.approx(3.394396287261952, abs=1e-9)
    assert summary["period_days"] == pytest.approx(14.74009747971348, abs=1e-9)
    C0 = jacobi_constant(state0, c.MU)
    assert C0 == pytest.approx(summary["jacobi_C0"], abs=1e-12)
    assert state0[2] == pytest.approx(0.035, abs=1e-15)  # z0 amplitude untouched
