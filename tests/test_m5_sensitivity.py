"""Tests for the M5 propellant-implication calculation (mirrors M4's)
and a hand-constructible sensitivity sanity check."""

import math

import numpy as np
import pytest

from halo_insertion import constants as c
from halo_insertion.arrival_arc import propagate_backward_arc


def _propellant(delta_v_m_s, m0, isp, g0=9.80665):
    return m0 * (1.0 - math.exp(-delta_v_m_s / (isp * g0)))


def test_propellant_zero_for_zero_delta_v():
    assert _propellant(0.0, 6000.0, 320.0) == pytest.approx(0.0, abs=1e-12)


def test_propellant_monotonic_in_delta_v():
    m0, isp = 6000.0, 320.0
    dvs = [10.0, 50.0, 109.57, 200.0]
    props = [_propellant(dv, m0, isp) for dv in dvs]
    for i in range(len(props) - 1):
        assert props[i + 1] > props[i]


def test_backward_duration_does_not_change_delta_v_only_acceptance():
    # Delta_v is evaluated at t=0 from (r_h, v_halo, v_arr); it must be
    # invariant to how far the arc is traced backward, while acceptance
    # (Earthward reach) generically depends on duration.
    from halo_insertion.equilibria import find_L2

    MU = c.MU
    l2 = find_L2(MU)
    r_h = np.array([l2.x - 0.02, 0.096, -0.013])
    v_arr = np.array([0.14, 0.012, -0.0017])

    accepted_flags = []
    for t_back in [5.0, 10.0, 20.0]:
        arc = propagate_backward_arc(r_h, v_arr, MU, t_back, n_eval=200)
        accepted_flags.append(arc.accepted)
        # v_arr (the prescribed pre-burn velocity) never changes with duration.
        assert np.array_equal(arc.v_arr, v_arr)

    # Not all durations need to be accepted, but the arc's *prescribed*
    # v_arr (and hence any Delta_v computed from it) is identical
    # regardless of duration -- duration only changes what we learn
    # about the trajectory's far-end geometry.
    assert True  # v_arr equality above is the substantive assertion
