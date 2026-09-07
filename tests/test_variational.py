"""Tests for the CR3BP state Jacobian (A matrix) and STM propagation."""

import numpy as np
import pytest

from halo_insertion import constants as c
from halo_insertion.cr3bp import cr3bp_rhs, effective_potential_gradient
from halo_insertion.equilibria import find_L2
from halo_insertion.propagation import propagate
from halo_insertion.variational import augmented_rhs, omega_hessian, propagate_with_stm, state_jacobian_A

MU = c.MU

ARBITRARY_STATES = [
    (0.5, 0.2, 0.1),
    (1.2, -0.05, 0.02),
    (-0.3, 0.4, -0.1),
    (0.85, 0.0, 0.15),
    (1.05, 0.03, 0.07),
]


@pytest.mark.parametrize("x,y,z", ARBITRARY_STATES)
def test_hessian_matches_finite_difference_of_gradient(x, y, z):
    h = 1e-6

    def grad(x, y, z):
        return np.array(effective_potential_gradient(x, y, z, MU))

    H_fd = np.zeros((3, 3))
    for i, eps in enumerate([(h, 0, 0), (0, h, 0), (0, 0, h)]):
        gp = grad(x + eps[0], y + eps[1], z + eps[2])
        gm = grad(x - eps[0], y - eps[1], z - eps[2])
        H_fd[:, i] = (gp - gm) / (2 * h)

    H_analytic = omega_hessian(x, y, z, MU)
    assert np.allclose(H_analytic, H_fd, atol=1e-7)


@pytest.mark.parametrize("x,y,z", ARBITRARY_STATES)
def test_A_matrix_structure(x, y, z):
    state = [x, y, z, 0.01, -0.02, 0.03]
    A = state_jacobian_A(state, MU)
    assert A.shape == (6, 6)
    assert np.allclose(A[0:3, 0:3], 0.0)
    assert np.allclose(A[0:3, 3:6], np.eye(3))
    assert A[3, 4] == pytest.approx(2.0)
    assert A[4, 3] == pytest.approx(-2.0)
    assert A[3, 3] == 0.0 and A[4, 4] == 0.0 and A[5, 3] == 0.0 and A[5, 4] == 0.0 and A[5, 5] == 0.0
    hess = omega_hessian(x, y, z, MU)
    assert np.allclose(A[3:6, 0:3], hess)


def test_phi0_is_identity():
    l2 = find_L2(MU)
    state0 = [l2.x - 0.02, 0.0, 0.02, 0.0, 0.05, 0.0]
    result = propagate_with_stm(state0, (0.0, 1e-9), MU)
    Phi0_approx = result.ode_result.y[6:, 0].reshape(6, 6)
    assert np.allclose(Phi0_approx, np.eye(6), atol=1e-12)


def test_short_time_stm_approximates_I_plus_A_dt():
    l2 = find_L2(MU)
    state0 = np.array([l2.x - 0.02, 0.01, 0.02, 0.0, 0.05, 0.0])
    dt = 1e-5
    result = propagate_with_stm(state0, (0.0, dt), MU)
    Phi_dt = result.Phi_final

    A = state_jacobian_A(state0, MU)
    I_plus_Adt = np.eye(6) + A * dt

    assert np.allclose(Phi_dt, I_plus_Adt, atol=1e-8)


def test_augmented_state_matches_generic_propagator():
    l2 = find_L2(MU)
    state0 = np.array([l2.x - 0.02, 0.01, 0.02, 0.0, 0.05, 0.0])
    t_end = 0.5

    stm_result = propagate_with_stm(state0, (0.0, t_end), MU)
    state_from_stm = stm_result.states[:, -1]

    generic_result = propagate(state0, (0.0, t_end), MU)
    state_from_generic = generic_result.y[:, -1]

    assert np.allclose(state_from_stm, state_from_generic, atol=1e-9, rtol=1e-9)


def test_stm_predicts_small_perturbation_to_first_order():
    l2 = find_L2(MU)
    state0 = np.array([l2.x - 0.02, 0.01, 0.02, 0.0, 0.05, 0.0])
    dt = 0.05

    stm_result = propagate_with_stm(state0, (0.0, dt), MU)
    Phi_dt = stm_result.Phi_final
    unperturbed_final = stm_result.states[:, -1]

    eps = 1e-6
    perturbation = np.array([eps, 0, 0, 0, 0, 0])
    predicted_delta = Phi_dt @ perturbation

    perturbed_result = propagate(state0 + perturbation, (0.0, dt), MU)
    actual_delta = perturbed_result.y[:, -1] - unperturbed_final

    assert np.allclose(predicted_delta, actual_delta, atol=1e-9, rtol=1e-3)


def test_tighter_perturbations_improve_linear_prediction():
    l2 = find_L2(MU)
    state0 = np.array([l2.x - 0.02, 0.01, 0.02, 0.0, 0.05, 0.0])
    dt = 0.05

    stm_result = propagate_with_stm(state0, (0.0, dt), MU)
    Phi_dt = stm_result.Phi_final
    unperturbed_final = stm_result.states[:, -1]

    errors = []
    for eps in [1e-3, 1e-4, 1e-5, 1e-6]:
        perturbation = np.array([eps, 0, 0, 0, 0, 0])
        predicted_delta = Phi_dt @ perturbation
        perturbed_result = propagate(state0 + perturbation, (0.0, dt), MU)
        actual_delta = perturbed_result.y[:, -1] - unperturbed_final
        errors.append(np.linalg.norm(actual_delta - predicted_delta))

    # STM linearization error should shrink monotonically (quadratically)
    # as the perturbation shrinks.
    for i in range(len(errors) - 1):
        assert errors[i + 1] < errors[i]


def test_augmented_rhs_rejects_wrong_shape():
    from halo_insertion.cr3bp import CR3BPStateError

    with pytest.raises(CR3BPStateError):
        augmented_rhs(0.0, np.zeros(10), MU)
