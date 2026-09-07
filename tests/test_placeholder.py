"""Packaging/metadata sanity tests.

Originally the M1 placeholder suite. As of M2, generic CR3BP dynamics
(cr3bp_rhs, propagate, equilibria, Jacobi constant) legitimately exist,
so this module's "nothing exists yet" guard is scoped to halo-specific
machinery (seed generation, differential correction, manifolds,
insertion solving) that is still M3+ only. Dynamics-specific tests live
in test_cr3bp.py, test_equilibria.py, test_jacobi.py,
test_propagation.py, and test_normalization.py.
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


def test_m1_and_m2_milestones_recorded():
    assert "M1" in hi.IMPLEMENTED_MILESTONES
    assert "M2" in hi.IMPLEMENTED_MILESTONES


def test_no_halo_solver_falsely_claimed_yet():
    # M2 must not expose any halo-specific API (seed generation,
    # differential correction, manifolds, insertion solving). This
    # guards against accidentally shipping/claiming that functionality
    # ahead of M3+.
    for name in ("differential_correct", "solve_insertion", "richardson_halo_seed"):
        assert not hasattr(hi, name), f"'{name}' should not exist until a later milestone"

    for name in hi.NOT_YET_IMPLEMENTED:
        assert not hasattr(hi, name), f"'{name}' is listed as not-yet-implemented but exists"


def test_generic_cr3bp_dynamics_now_exist():
    # Sanity check that M2's generic (non-halo) dynamics API landed.
    for name in ("cr3bp_rhs", "propagate", "jacobi_constant", "find_L1", "find_L2", "find_L3"):
        assert hasattr(hi, name), f"'{name}' should exist as of M2"
