"""Dimensional <-> nondimensional conversion helpers (Earth-Moon CR3BP).

All functions are pure and take/return plain floats or numpy arrays.
Units are explicit in every function name — there is no ambiguous
"convert velocity" helper. See DESIGN.md Section 3 for the conversion
rules these implement.
"""

from __future__ import annotations

from .constants import DU_KM, TU_DAYS, TU_S, VSTAR_KM_S, VSTAR_M_S

__all__ = [
    "km_to_nondim_position",
    "nondim_position_to_km",
    "km_s_to_nondim_velocity",
    "nondim_velocity_to_km_s",
    "m_s_to_nondim_velocity",
    "nondim_velocity_to_m_s",
    "seconds_to_nondim_time",
    "nondim_time_to_seconds",
    "days_to_nondim_time",
    "nondim_time_to_days",
]


def km_to_nondim_position(r_km):
    """Convert position(s) in km to nondimensional CR3BP distance units (DU)."""
    return r_km / DU_KM


def nondim_position_to_km(r_nondim):
    """Convert nondimensional CR3BP position(s) (DU) to km."""
    return r_nondim * DU_KM


def km_s_to_nondim_velocity(v_km_s):
    """Convert velocity(ies) in km/s to nondimensional CR3BP units (DU/TU)."""
    return v_km_s / VSTAR_KM_S


def nondim_velocity_to_km_s(v_nondim):
    """Convert nondimensional CR3BP velocity(ies) (DU/TU) to km/s."""
    return v_nondim * VSTAR_KM_S


def m_s_to_nondim_velocity(v_m_s):
    """Convert velocity(ies) in m/s to nondimensional CR3BP units (DU/TU)."""
    return v_m_s / VSTAR_M_S


def nondim_velocity_to_m_s(v_nondim):
    """Convert nondimensional CR3BP velocity(ies) (DU/TU) to m/s."""
    return v_nondim * VSTAR_M_S


def seconds_to_nondim_time(t_s):
    """Convert time in seconds to nondimensional CR3BP time units (TU)."""
    return t_s / TU_S


def nondim_time_to_seconds(t_nondim):
    """Convert nondimensional CR3BP time (TU) to seconds."""
    return t_nondim * TU_S


def days_to_nondim_time(t_days):
    """Convert time in days to nondimensional CR3BP time units (TU)."""
    return t_days / TU_DAYS


def nondim_time_to_days(t_nondim):
    """Convert nondimensional CR3BP time (TU) to days."""
    return t_nondim * TU_DAYS
