"""halo_insertion — Earth-Moon L2 halo-orbit insertion delta-v estimate.

Milestone status
-----------------
M1: mission definition, CR3BP governing-equation documentation, L2
location, insertion-delta-v definition, hand estimates.

M2 (current): normalized Earth-Moon CR3BP dynamics — effective
potential, equations of motion, Jacobi constant, collinear equilibrium
points (L1/L2/L3), a generic numerical propagator, and dimensional/
nondimensional conversion helpers. All verified against M1's hand
values and against independent cross-checks (see DESIGN.md Section 18).

No halo-orbit seed, differential correction, manifold computation, or
insertion-Delta-v solver is implemented yet. Do not import such
routines from this module — they do not exist until M3+.

See DESIGN.md at the repository root for full derivations, assumptions,
and the milestone roadmap.
"""

from . import constants
from .cr3bp import (
    CR3BPStateError,
    cr3bp_rhs,
    effective_potential,
    effective_potential_gradient,
    jacobi_constant,
    primary_distances,
)
from .equilibria import EquilibriumPoint, find_L1, find_L2, find_L3
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

__version__ = "0.2.0"

# Milestones implemented so far. Kept as plain data (not behavior) so
# tests can assert on project status without ambiguity.
IMPLEMENTED_MILESTONES = ("M1", "M2")

# Explicitly NOT implemented yet (M3+). Listed for clarity/testability,
# not as a promise of interface — these names are not importable.
NOT_YET_IMPLEMENTED = (
    "richardson_halo_seed",
    "differential_correction",
    "state_transition_matrix",
    "monodromy_matrix",
    "manifold_targeting",
    "insertion_solver",
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
    "IMPLEMENTED_MILESTONES",
    "NOT_YET_IMPLEMENTED",
]
