"""Tests for M1-constant reproduction and dimensional round-trips."""

import numpy as np
import pytest

from halo_insertion import constants as c
from halo_insertion import normalization as norm


def test_m1_mu_reproduced():
    assert c.MU == pytest.approx(0.0121505839, abs=1e-9)


def test_m1_du_reproduced():
    assert c.DU_KM == pytest.approx(384400.0, abs=1e-9)


def test_m1_tu_reproduced():
    assert c.TU_DAYS == pytest.approx(4.342480, abs=1e-5)


def test_m1_vstar_reproduced():
    assert c.VSTAR_M_S == pytest.approx(1024.547, abs=1e-2)


def test_primary_positions_match_convention():
    assert c.EARTH_POS_NONDIM == (-c.MU, 0.0, 0.0)
    assert c.MOON_POS_NONDIM == (1.0 - c.MU, 0.0, 0.0)


@pytest.mark.parametrize("r_km", [0.0, 1.0, 384400.0, 64514.9, -12345.6])
def test_position_round_trip(r_km):
    r_nondim = norm.km_to_nondim_position(r_km)
    r_back = norm.nondim_position_to_km(r_nondim)
    assert r_back == pytest.approx(r_km, rel=1e-13, abs=1e-9)


@pytest.mark.parametrize("v_km_s", [0.0, 1.0245468552412853, -0.5, 3.14159])
def test_velocity_round_trip_km_s(v_km_s):
    v_nondim = norm.km_s_to_nondim_velocity(v_km_s)
    v_back = norm.nondim_velocity_to_km_s(v_nondim)
    assert v_back == pytest.approx(v_km_s, rel=1e-13, abs=1e-12)


@pytest.mark.parametrize("v_m_s", [0.0, 25.0, 100.0, -50.0])
def test_velocity_round_trip_m_s(v_m_s):
    v_nondim = norm.m_s_to_nondim_velocity(v_m_s)
    v_back = norm.nondim_velocity_to_m_s(v_nondim)
    assert v_back == pytest.approx(v_m_s, rel=1e-13, abs=1e-9)


def test_velocity_km_s_and_m_s_paths_agree():
    v_km_s = 0.05
    v_nondim_via_km = norm.km_s_to_nondim_velocity(v_km_s)
    v_nondim_via_m = norm.m_s_to_nondim_velocity(v_km_s * 1000.0)
    assert v_nondim_via_km == pytest.approx(v_nondim_via_m, rel=1e-12)


@pytest.mark.parametrize("t_s", [0.0, 375190.259, 86400.0, -1000.0])
def test_time_round_trip_seconds(t_s):
    t_nondim = norm.seconds_to_nondim_time(t_s)
    t_back = norm.nondim_time_to_seconds(t_nondim)
    assert t_back == pytest.approx(t_s, rel=1e-13, abs=1e-9)


@pytest.mark.parametrize("t_days", [0.0, 4.34248, 27.3, -1.5])
def test_time_round_trip_days(t_days):
    t_nondim = norm.days_to_nondim_time(t_days)
    t_back = norm.nondim_time_to_days(t_nondim)
    assert t_back == pytest.approx(t_days, rel=1e-13, abs=1e-9)


def test_one_tu_equals_period_over_2pi_days():
    # By construction TU = 1/n, so one full 2*pi nondimensional time unit
    # of synodic rotation equals one Earth-Moon synodic-ish period scale;
    # here we just check TU_DAYS matches n directly.
    assert c.TU_DAYS == pytest.approx(1.0 / c.N_RAD_S / 86400.0, rel=1e-14)


def test_position_conversion_is_array_safe():
    r_km = np.array([384400.0, 0.0, 64514.9])
    r_nondim = norm.km_to_nondim_position(r_km)
    r_back = norm.nondim_position_to_km(r_nondim)
    assert np.allclose(r_back, r_km, rtol=1e-13)
