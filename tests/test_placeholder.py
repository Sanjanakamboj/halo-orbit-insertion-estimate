"""M1 placeholder tests.

These only verify packaging/metadata sanity and that no numerical
CR3BP solver is falsely claimed to exist yet. Real dynamics tests
(equilibrium residuals, Jacobi conservation, periodicity closure, etc.,
per DESIGN.md Section 11) arrive starting M2.
"""

import halo_insertion as hi


def test_package_imports():
    assert hi is not None


def test_version_string_present():
    assert isinstance(hi.__version__, str)
    assert hi.__version__ != ""
    parts = hi.__version__.split(".")
    assert len(parts) == 3
    assert all(p.isdigit() for p in parts)


def test_m1_milestone_recorded():
    assert "M1" in hi.IMPLEMENTED_MILESTONES


def test_no_solver_falsely_claimed_yet():
    # M1 must not expose any numerical dynamics API. This guards against
    # accidentally shipping/claiming solver functionality ahead of M2.
    for name in ("propagate", "cr3bp_rhs", "differential_correct", "solve_insertion"):
        assert not hasattr(hi, name), f"'{name}' should not exist until a later milestone"

    for name in hi.NOT_YET_IMPLEMENTED:
        assert not hasattr(hi, name), f"'{name}' is listed as not-yet-implemented but exists"
