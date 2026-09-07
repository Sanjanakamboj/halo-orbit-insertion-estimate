"""halo_insertion — Earth-Moon L2 halo-orbit insertion delta-v estimate.

Milestone status
-----------------
M1 (current): mission definition, CR3BP governing-equation documentation,
L2 location, insertion-delta-v definition, and hand estimates only.

No CR3BP propagator, differential corrector, or insertion solver is
implemented in this package yet. Do not import numerical dynamics
routines from this module — they do not exist until M2+.

See DESIGN.md at the repository root for full derivations, assumptions,
and the milestone roadmap.
"""

__version__ = "0.1.0"

# Milestones implemented so far. Kept as plain data (not behavior) so
# tests can assert on project status without any dynamics code existing.
IMPLEMENTED_MILESTONES = ("M1",)

# Explicitly NOT implemented yet (M2+). Listed for clarity/testability,
# not as a promise of interface — these names are not importable.
NOT_YET_IMPLEMENTED = (
    "cr3bp_propagator",
    "differential_correction",
    "manifold_targeting",
    "insertion_solver",
)
