"""halo_insertion — Earth-Moon L2 halo-orbit insertion delta-v estimate.

Milestone status
-----------------
M1: mission definition, CR3BP governing-equation documentation, L2
location, insertion-delta-v definition, hand estimates.

M2: normalized Earth-Moon CR3BP dynamics — effective potential,
equations of motion, Jacobi constant, collinear equilibrium points
(L1/L2/L3), a generic numerical propagator, and dimensional/
nondimensional conversion helpers.

M3: the CR3BP state Jacobian and 6x6 variational equations,
state-transition-matrix (STM) propagation, a linearized L2 halo-orbit
seed, a symmetry-based Newton differential corrector, and one verified
numerically periodic Earth-Moon L2 halo orbit.

M4 (current): an explicit transfer-arrival state model (Jacobi-
consistent speed + documented direction models), a dense insertion-
phase sweep and bounded-scalar refinement over the frozen M3 halo
orbit, sensitivity studies (arrival speed, arrival direction, local
phase), and the project's first defensible local halo-insertion Delta-v
estimate. See DESIGN.md's M4 section for the full arrival-model
definition, results, and limitations.

**M4 computes a local velocity-matching insertion estimate at the
corrected halo orbit. It does not compute or optimize the complete
Earth-to-L2 transfer trajectory.** No manifold, launch/TLI, or
stationkeeping analysis is implemented. Do not import such routines
from this module — they do not exist yet.

See DESIGN.md at the repository root for full derivations, assumptions,
and the milestone roadmap.
"""

from . import constants
from .arrival import (
    ArrivalStateError,
    arrival_speed_from_jacobi,
    arrival_velocity,
    delta_v_vector,
    direction_aligned_with_halo,
    direction_radial_from_earth,
    direction_rotated_about_z,
    earth_position,
    jacobi_from_state,
)
from .cr3bp import (
    CR3BPStateError,
    cr3bp_rhs,
    effective_potential,
    effective_potential_gradient,
    jacobi_constant,
    primary_distances,
)
from .differential_correction import (
    CorrectionResult,
    DifferentialCorrectionError,
    differential_correct_halo,
    propagate_half_period,
)
from .equilibria import EquilibriumPoint, find_L1, find_L2, find_L3
from .halo import HaloOrbitResult, build_l2_halo, full_period_closure, monodromy_matrix
from .halo_seed import LinearHaloSeed, linear_halo_seed
from .insertion import InsertionCandidate, evaluate_insertion_at_phase, halo_state_at_phase, phase_sweep, refine_best_phase
from .normalization import (
    days_to_nondim_time,
    km_s_to_nondim_velocity,
    km_to_nondim_position,
    m_s_to_nondim_velocity,
    nondim_position_to_km,
    nondim_time_to_days,
    nondim_time_to_seconds,
    nondim_velocity_to_km_s,
    nondim_velocity_to_m_s,
    seconds_to_nondim_time,
)
from .propagation import PropagationSettings, propagate
from .variational import augmented_rhs, omega_hessian, propagate_with_stm, state_jacobian_A

__version__ = "0.4.0"

# Milestones implemented so far. Kept as plain data (not behavior) so
# tests can assert on project status without ambiguity.
IMPLEMENTED_MILESTONES = ("M1", "M2", "M3", "M4")

# Explicitly NOT implemented yet (M5+). Listed for clarity/testability,
# not as a promise of interface — these names are not importable.
# NOTE: `richardson_halo_seed` stays listed here deliberately: M3 uses a
# linearized (first-order) CR3BP variational seed, NOT the third-order
# Richardson approximation (see DESIGN.md's M3 section for why, and for
# the seed's actual provenance).
NOT_YET_IMPLEMENTED = (
    "richardson_halo_seed",
    "stable_manifold",
    "unstable_manifold",
    "transfer_trajectory_optimizer",
    "translunar_injection_model",
    "stationkeeping_design",
)

__all__ = [
    "constants",
    "CR3BPStateError",
    "cr3bp_rhs",
    "effective_potential",
    "effective_potential_gradient",
    "jacobi_constant",
    "primary_distances",
    "EquilibriumPoint",
    "find_L1",
    "find_L2",
    "find_L3",
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
    "PropagationSettings",
    "propagate",
    "omega_hessian",
    "state_jacobian_A",
    "augmented_rhs",
    "propagate_with_stm",
    "LinearHaloSeed",
    "linear_halo_seed",
    "DifferentialCorrectionError",
    "CorrectionResult",
    "differential_correct_halo",
    "propagate_half_period",
    "HaloOrbitResult",
    "build_l2_halo",
    "full_period_closure",
    "monodromy_matrix",
    "ArrivalStateError",
    "earth_position",
    "arrival_speed_from_jacobi",
    "direction_aligned_with_halo",
    "direction_radial_from_earth",
    "direction_rotated_about_z",
    "arrival_velocity",
    "delta_v_vector",
    "jacobi_from_state",
    "InsertionCandidate",
    "halo_state_at_phase",
    "evaluate_insertion_at_phase",
    "phase_sweep",
    "refine_best_phase",
    "IMPLEMENTED_MILESTONES",
    "NOT_YET_IMPLEMENTED",
]
