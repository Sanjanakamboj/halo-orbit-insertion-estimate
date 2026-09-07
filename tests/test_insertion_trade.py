"""Tests for the M4 propellant-implication calculation and a final
no-regression check against M1-M3."""

import math

import pytest

from halo_insertion import constants as c
from halo_insertion.equilibria import find_L2


def _propellant(delta_v_m_s, m0, isp, g0=9.80665):
    return m0 * (1.0 - math.exp(-delta_v_m_s / (isp * g0)))


def test_propellant_zero_for_zero_delta_v():
    m_prop = _propellant(0.0, 6000.0, 320.0)
    assert m_prop == pytest.approx(0.0, abs=1e-12)


def test_propellant_monotonic_in_delta_v():
    m0, isp = 6000.0, 320.0
    dvs = [10.0, 25.0, 50.0, 100.0, 150.0, 300.0]
    props = [_propellant(dv, m0, isp) for dv in dvs]
    for i in range(len(props) - 1):
        assert props[i + 1] > props[i]


def test_propellant_higher_isp_uses_less_propellant():
    m0, dv = 6000.0, 110.0
    m_320 = _propellant(dv, m0, 320.0)
    m_450 = _propellant(dv, m0, 450.0)
    assert m_450 < m_320


def test_propellant_bounded_by_m0():
    m0 = 6000.0
    m_prop = _propellant(5000.0, m0, 320.0)  # unrealistically large dv
    assert 0.0 < m_prop < m0


def test_no_m1_m2_m3_regression():
    assert c.MU == pytest.approx(0.0121505839, abs=1e-9)
    l2 = find_L2(c.MU)
    assert l2.x == pytest.approx(1.1556821589, abs=1e-9)
    assert c.DU_KM == pytest.approx(384400.0, abs=1e-9)
    assert c.VSTAR_M_S == pytest.approx(1024.547, abs=1e-2)
