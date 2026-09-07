"""Tests for the generic CR3BP propagator: equilibrium-state stationarity,
Jacobi conservation on a non-equilibrium verification trajectory (NOT a
halo orbit), and a 3-tolerance numerical convergence study.
"""

import numpy as np
import pytest

from halo_insertion import constants as c
from halo_insertion.cr3bp import jacobi_constant
from halo_insertion.equilibria import find_L1, find_L2, find_L3
from halo_insertion.normalization import days_to_nondim_time
from halo_insertion.propagation import PropagationSettings, propagate

MU = c.MU


@pytest.mark.parametrize("finder", [find_L1, find_L2, find_L3])
def test_equilibrium_propagation_remains_stationary(finder):
    eq = finder(MU)
    state0 = [eq.x, 0.0, 0.0, 0.0, 0.0, 0.0]
    t_end = days_to_nondim_time(30.0)

    result = propagate(state0, (0.0, t_end), MU)
    assert result.success

    final_state = result.y[:, -1]
    drift = np.linalg.norm(final_state - np.array(state0))
    # Collinear points are dynamically unstable equilibria (this is
    # physically expected, not a bug) — round-off-level perturbations
    # from an exact fixed point can grow over 30 days. The bound here
    # is loose enough to pass on genuine round-off drift but tight
    # enough to catch a real implementation error (wrong sign, wrong
    # offset, etc.), which would show up as drift many orders larger.
    assert drift < 1e-6


def _m2_verification_trajectory_ic():
    """A benign, non-equilibrium CR3BP state near L2 — NOT a halo orbit.

    Displaced from L2 in x, with a small y/z offset and small transverse
    velocity, used purely to verify propagator behavior (Jacobi
    conservation, convergence, finiteness).
    """
    l2 = find_L2(MU)
    state0 = np.array([l2.x - 0.02, 0.01, 0.005, 0.0, 0.01, 0.005])
    t_end = days_to_nondim_time(10.0)
    return state0, t_end


def test_m2_verification_trajectory_stays_finite():
    state0, t_end = _m2_verification_trajectory_ic()
    result = propagate(state0, (0.0, t_end), MU)
    assert result.success
    assert np.all(np.isfinite(result.y))


def test_m2_verification_trajectory_conserves_jacobi():
    state0, t_end = _m2_verification_trajectory_ic()
    result = propagate(state0, (0.0, t_end), MU)
    assert result.success

    C0 = jacobi_constant(state0, MU)
    C_over_time = np.array([jacobi_constant(result.y[:, i], MU) for i in range(result.y.shape[1])])
    max_drift = np.max(np.abs(C_over_time - C0))

    # Default M2 tolerances (rtol=1e-11) should hold Jacobi drift well
    # below 1e-8 over a 10-day nondimensional propagation.
    assert max_drift < 1e-8


def test_numerical_convergence_three_tolerance_levels():
    """Loose/medium/tight integration settings on the same verification
    trajectory: terminal-state difference and Jacobi drift must shrink
    monotonically as tolerances tighten — this is what "convergence"
    means here, not a single-tolerance claim.
    """
    state0, t_end = _m2_verification_trajectory_ic()
    C0 = jacobi_constant(state0, MU)

    settings = {
        "loose": PropagationSettings(method="DOP853", rtol=1e-6, atol=1e-8),
        "medium": PropagationSettings(method="DOP853", rtol=1e-9, atol=1e-10),
        "tight": PropagationSettings(method="DOP853", rtol=1e-12, atol=1e-13),
    }

    finals = {}
    jacobi_drift = {}
    for name, s in settings.items():
        result = propagate(state0, (0.0, t_end), MU, settings=s)
        assert result.success
        finals[name] = result.y[:, -1]
        jacobi_drift[name] = abs(jacobi_constant(result.y[:, -1], MU) - C0)

    # Jacobi drift must shrink as tolerances tighten.
    assert jacobi_drift["tight"] < jacobi_drift["medium"] < jacobi_drift["loose"]

    # Terminal-state difference between medium and tight must be much
    # smaller than between loose and tight (convergence, not noise).
    diff_medium_tight = np.linalg.norm(finals["medium"] - finals["tight"])
    diff_loose_tight = np.linalg.norm(finals["loose"] - finals["tight"])
    assert diff_medium_tight < diff_loose_tight
    assert diff_medium_tight < 1e-6


def test_planar_initial_condition_stays_planar_under_propagation():
    l2 = find_L2(MU)
    state0 = [l2.x - 0.02, 0.01, 0.0, 0.0, 0.01, 0.0]  # z = zdot = 0
    t_end = days_to_nondim_time(10.0)

    result = propagate(state0, (0.0, t_end), MU)
    assert result.success

    z_values = result.y[2, :]
    zdot_values = result.y[5, :]
    assert np.max(np.abs(z_values)) < 1e-10
    assert np.max(np.abs(zdot_values)) < 1e-10


def test_propagate_exposes_raw_ode_result():
    l2 = find_L2(MU)
    state0 = [l2.x, 0.0, 0.0, 0.0, 0.0, 0.0]
    result = propagate(state0, (0.0, 1.0), MU)
    # Solver failures/success must be visible to the caller, not hidden.
    assert hasattr(result, "success")
    assert hasattr(result, "t")
    assert hasattr(result, "y")
