"""Packaging/metadata sanity tests.

Originally the M1 placeholder suite. As of M2, generic CR3BP dynamics
(cr3bp_rhs, propagate, equilibria, Jacobi constant) legitimately exist;
as of M3, the variational-equations/STM machinery and a halo seed +
differential corrector + one verified periodic orbit also exist; as of
M4, an arrival-state model and insertion phase-sweep/refinement also
exist; as of M5, dynamically propagated backward arrival arcs and
round-trip verification also exist. M6 is the final milestone (audit,
reproducibility, CI, portfolio packaging) and adds no new trajectory-
physics API. This module's "nothing exists" guard is scoped to
functionality intentionally out of scope for the whole project
(manifold design, transfer optimization, TLI, stationkeeping).
Dynamics- and halo-specific tests live in the other test_*.py modules.
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


def test_m1_through_m6_milestones_recorded():
    for m in ("M1", "M2", "M3", "M4", "M5", "M6"):
        assert m in hi.IMPLEMENTED_MILESTONES


def test_out_of_scope_solver_never_claimed():
    # This project never exposes manifold-design/transfer-optimization/
    # TLI/stationkeeping API -- that scope boundary is permanent, not
    # just "not yet". Guards against accidentally shipping/claiming
    # such functionality.
    for name in ("stable_manifold", "unstable_manifold", "transfer_trajectory_optimizer", "translunar_injection_model"):
        assert not hasattr(hi, name), f"'{name}' is out of scope for this project and must not exist"

    for name in hi.NOT_YET_IMPLEMENTED:
        assert not hasattr(hi, name), f"'{name}' is listed as not-yet-implemented but exists"


def test_generic_cr3bp_dynamics_now_exist():
    # Sanity check that M2's generic (non-halo) dynamics API landed.
    for name in ("cr3bp_rhs", "propagate", "jacobi_constant", "find_L1", "find_L2", "find_L3"):
        assert hasattr(hi, name), f"'{name}' should exist as of M2"


def test_m3_halo_machinery_now_exists():
    # Sanity check that M3's variational/seed/corrector/halo API landed.
    for name in (
        "state_jacobian_A",
        "propagate_with_stm",
        "linear_halo_seed",
        "differential_correct_halo",
        "build_l2_halo",
        "monodromy_matrix",
    ):
        assert hasattr(hi, name), f"'{name}' should exist as of M3"


def test_m4_arrival_and_insertion_machinery_now_exists():
    # Sanity check that M4's arrival-model/insertion API landed.
    for name in (
        "arrival_speed_from_jacobi",
        "direction_radial_from_earth",
        "delta_v_vector",
        "phase_sweep",
        "refine_best_phase",
        "evaluate_insertion_at_phase",
    ):
        assert hasattr(hi, name), f"'{name}' should exist as of M4"


def test_m5_arrival_arc_machinery_now_exists():
    # Sanity check that M5's backward-arc/round-trip API landed.
    for name in (
        "propagate_backward_arc",
        "round_trip_check",
        "EARTH_DISTANCE_THRESHOLD_DU",
        "MOON_GUARD_KM",
    ):
        assert hasattr(hi, name), f"'{name}' should exist as of M5"
