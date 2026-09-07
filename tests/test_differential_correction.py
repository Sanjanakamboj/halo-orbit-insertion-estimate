"""Tests for the half-period event handling and Newton differential
corrector."""

import numpy as np
import pytest

from halo_insertion import constants as c
from halo_insertion.differential_correction import (
    DifferentialCorrectionError,
    differential_correct_halo,
    propagate_half_period,
)
from halo_insertion.equilibria import find_L2
from halo_insertion.halo_seed import linear_halo_seed

MU = c.MU


def test_half_period_event_does_not_trigger_at_t0():
    l2 = find_L2(MU)
    seed = linear_halo_seed(l2.x, MU, 0.03, z_sign=1.0)
    half = propagate_half_period(seed.state0, MU, t_max=6.0)
    # The detected crossing must be well after t=0, not the trivial
    # departure from the seed's own y(0)=0 initial condition.
    assert half.t_half > 0.1


def test_half_period_crossing_is_physically_the_first_return():
    l2 = find_L2(MU)
    seed = linear_halo_seed(l2.x, MU, 0.03, z_sign=1.0)
    half = propagate_half_period(seed.state0, MU, t_max=6.0)
    s_half = half.Y_half[0:6]
    assert abs(s_half[1]) < 1e-8  # y ~ 0 at the crossing


def test_corrector_converges_from_seed():
    l2 = find_L2(MU)
    seed = linear_halo_seed(l2.x, MU, 0.03, z_sign=1.0)
    result = differential_correct_halo(seed.state0, MU, tol=1e-11)
    assert result.converged


def test_correction_residual_decreases_materially():
    l2 = find_L2(MU)
    seed = linear_halo_seed(l2.x, MU, 0.03, z_sign=1.0)
    result = differential_correct_halo(seed.state0, MU, tol=1e-11)
    first_residual = np.hypot(result.iterations[0].xdot_residual, result.iterations[0].zdot_residual)
    last_residual = np.hypot(result.iterations[-1].xdot_residual, result.iterations[-1].zdot_residual)
    assert last_residual < first_residual * 1e-6


def test_corrected_half_period_xdot_residual_small():
    l2 = find_L2(MU)
    seed = linear_halo_seed(l2.x, MU, 0.03, z_sign=1.0)
    result = differential_correct_halo(seed.state0, MU, tol=1e-11)
    assert abs(result.iterations[-1].xdot_residual) < 1e-10


def test_corrected_half_period_zdot_residual_small():
    l2 = find_L2(MU)
    seed = linear_halo_seed(l2.x, MU, 0.03, z_sign=1.0)
    result = differential_correct_halo(seed.state0, MU, tol=1e-11)
    assert abs(result.iterations[-1].zdot_residual) < 1e-10


def test_z0_held_fixed_throughout():
    l2 = find_L2(MU)
    seed = linear_halo_seed(l2.x, MU, 0.03, z_sign=1.0)
    result = differential_correct_halo(seed.state0, MU, tol=1e-11)
    for it in result.iterations:
        assert it.z0 == pytest.approx(0.03)
    assert result.state0[2] == pytest.approx(0.03)


def test_corrector_raises_on_seed_too_far_from_periodic():
    # An intentionally bad "seed" (arbitrary state, not derived from
    # linear halo theory) should fail cleanly rather than silently
    # return garbage.
    l2 = find_L2(MU)
    bad_seed = [l2.x + 0.5, 0.0, 0.2, 0.0, 2.0, 0.0]
    with pytest.raises(DifferentialCorrectionError):
        differential_correct_halo(bad_seed, MU, tol=1e-11, max_iter=15)
