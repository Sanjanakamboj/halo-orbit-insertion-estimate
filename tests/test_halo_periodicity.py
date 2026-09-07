"""Full-orbit verification tests for the M3 corrected L2 halo orbit:
periodicity closure, Jacobi conservation, symmetry, geometry, tighter
propagation, and the monodromy matrix. Also a no-regression check
against M1/M2.
"""

import numpy as np
import pytest

from halo_insertion import constants as c
from halo_insertion.cr3bp import jacobi_constant
from halo_insertion.equilibria import find_L2
from halo_insertion.halo import build_l2_halo, full_period_closure, monodromy_matrix
from halo_insertion.propagation import PropagationSettings, propagate

MU = c.MU
AZ = 0.035  # matches scripts/build_m3_halo.py's target amplitude


@pytest.fixture(scope="module")
def halo():
    l2 = find_L2(MU)
    return build_l2_halo(l2.x, MU, AZ, z_sign=1.0, tol=1e-11)


def test_full_period_closure_small(halo):
    _, closure, _ = full_period_closure(halo.state0, halo.period, MU)
    pos_norm_km = np.linalg.norm(closure[:3]) * c.DU_KM
    vel_norm_m_s = np.linalg.norm(closure[3:]) * c.VSTAR_M_S
    assert pos_norm_km < 1e-2  # << 1 cm-scale would be excessive; this is generous but catches regressions
    assert vel_norm_m_s < 1e-2


def test_jacobi_drift_small_over_full_orbit(halo):
    t_eval = np.linspace(0.0, halo.period, 1000)
    res = propagate(halo.state0, (0.0, halo.period), MU, t_eval=t_eval)
    assert res.success
    C0 = jacobi_constant(halo.state0, MU)
    C_t = np.array([jacobi_constant(res.y[:, i], MU) for i in range(res.y.shape[1])])
    max_drift = np.max(np.abs(C_t - C0))
    assert max_drift < 1e-8


def test_orbit_is_nonplanar(halo):
    t_eval = np.linspace(0.0, halo.period, 500)
    res = propagate(halo.state0, (0.0, halo.period), MU, t_eval=t_eval)
    max_abs_z = np.max(np.abs(res.y[2]))
    assert max_abs_z > 0.01  # genuinely 3D, not a planar Lyapunov orbit (z ~ 0)


def test_orbit_centered_near_l2(halo):
    l2 = find_L2(MU)
    t_eval = np.linspace(0.0, halo.period, 500)
    res = propagate(halo.state0, (0.0, halo.period), MU, t_eval=t_eval)
    x_center = 0.5 * (res.y[0].min() + res.y[0].max())
    assert abs(x_center - l2.x) < 0.1  # near the L2 neighborhood, not another equilibrium region


def test_half_period_symmetry_residual_small(halo):
    half_state = propagate(halo.state0, (0.0, halo.half_period), MU).y[:, -1]
    assert abs(half_state[1]) < 1e-8  # y
    assert abs(half_state[3]) < 1e-8  # xdot
    assert abs(half_state[5]) < 1e-8  # zdot


def test_dimensional_conversions_self_consistent(halo):
    period_days_direct = halo.period_days
    period_days_recomputed = halo.period * c.TU_DAYS
    assert period_days_direct == pytest.approx(period_days_recomputed, rel=1e-12)

    z_amplitude_km = halo.state0[2] * c.DU_KM
    assert z_amplitude_km == pytest.approx(AZ * c.DU_KM, rel=1e-12)


def test_tighter_propagation_reproduces_corrected_orbit(halo):
    default_final, default_closure, _ = full_period_closure(halo.state0, halo.period, MU)
    tight_settings = PropagationSettings(method="DOP853", rtol=1e-13, atol=1e-14)
    tight_final, tight_closure, _ = full_period_closure(halo.state0, halo.period, MU, settings=tight_settings)

    # Tighter tolerance should not produce a wildly different final
    # state, and should not be materially worse than the default.
    assert np.linalg.norm(tight_final - default_final) < 1e-6
    assert np.linalg.norm(tight_closure[:3]) <= np.linalg.norm(default_closure[:3]) * 10  # not worse by orders of magnitude


def test_independent_integrator_cross_check(halo):
    default_final, _, _ = full_period_closure(halo.state0, halo.period, MU)
    rk45_settings = PropagationSettings(method="RK45", rtol=1e-12, atol=1e-13)
    rk45_final, _, _ = full_period_closure(halo.state0, halo.period, MU, settings=rk45_settings)
    assert np.linalg.norm(default_final - rk45_final) < 1e-6


def test_monodromy_matrix_finite_and_nonsingular(halo):
    M = monodromy_matrix(halo.state0, halo.period, MU)
    assert np.all(np.isfinite(M))
    det = np.linalg.det(M)
    assert abs(det) > 1e-6  # nonsingular
    assert det == pytest.approx(1.0, abs=1e-6)  # CR3BP flow is symplectic/volume-preserving


def test_monodromy_eigenvalues_show_reciprocal_pairing_and_instability(halo):
    M = monodromy_matrix(halo.state0, halo.period, MU)
    eigvals = np.linalg.eigvals(M)
    mags = np.sort(np.abs(eigvals))
    # Largest and smallest magnitude should be (approximately) reciprocal,
    # as expected for CR3BP Hamiltonian dynamics.
    assert mags[0] * mags[-1] == pytest.approx(1.0, abs=1e-4)
    # The orbit should be linearly unstable (a real eigenvalue >> 1).
    assert mags[-1] > 1.0 + 1e-3


def test_no_m1_m2_regression():
    from halo_insertion import constants as const
    from halo_insertion.equilibria import find_L2 as fl2

    assert const.MU == pytest.approx(0.0121505839, abs=1e-9)
    l2 = fl2(const.MU)
    assert l2.x == pytest.approx(1.1556821589, abs=1e-9)
