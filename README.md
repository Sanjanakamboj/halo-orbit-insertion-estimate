# Halo-Orbit Insertion Δv Estimate

**Portfolio deliverable:** Halo-orbit insertion Δv estimate for an
Earth–Moon L2 telescope-class mission, with clearly stated assumptions
and a concrete verification plan.

**Status: Milestone 5 (M5) complete — the M4 local insertion-Δv
estimate has been validated against a dynamically propagated CR3BP
arrival arc and confirmed to round-trip numerical precision (below).
The repository is NOT yet portfolio-ready: a final packaging/audit
pass (M6) remains.**

Full derivations, equations, and assumptions live in [DESIGN.md](DESIGN.md).
This README summarizes the current state.

---

## Objective

Estimate the impulsive Δv required to insert a telescope-class
spacecraft, arriving on an already-targeted transfer trajectory, onto a
representative Earth–Moon L2 halo orbit.

## What "insertion Δv" means here

At a chosen point near the target halo orbit:

```
Delta_v = v_halo - v_arrival
```

This is **only** the velocity-matching burn at the end of the
transfer — it explicitly excludes launch Δv, translunar injection,
midcourse corrections, the Earth-to-L2 transfer Δv itself, and ongoing
stationkeeping. See [DESIGN.md §7](DESIGN.md#7-definition-of-insertion-δv-for-this-project)
for the full breakdown of each term.

## M1 reference scenario

| | |
|---|---|
| System | Earth–Moon |
| Destination | L2, southern halo family |
| Spacecraft class | Telescope-class observatory |
| Wet mass | 6,000 kg |
| Maneuver model | Single impulsive burn |
| Arrival assumption | Transfer already targeted near the halo insertion region |

## Halo family / seed status

M1's target orbit was a **halo-orbit seed** (literature-scale
placeholder). **As of M3, this project has a genuine numerically
corrected periodic Earth–Moon L2 halo orbit** (see below) — the M1
placeholder is superseded for geometry/period purposes. **As of M4, a
real insertion Δv has been computed from this corrected orbit** (see
the headline result below); M1's hand estimates remain as the original
order-of-magnitude sanity bound they were always intended to be.

## Key normalization constants (Earth–Moon CR3BP)

```
mu  = 0.0121505839          (mass ratio)
DU  = 384,400 km            (distance unit)
TU  = 4.342480 days         (time unit)
V*  = 1024.547 m/s          (characteristic velocity)
```

## L2 location

```
x_L2 (nondimensional)      = 1.1556821589
distance from Moon center  ≈ 64,515 km
distance beyond Moon surface ≈ 62,778 km
```

Consistent with the expected ~60,000–70,000 km-beyond-Moon range for
Earth–Moon L2. Derivation: [DESIGN.md §6](DESIGN.md#6-l2-location).

## M2: CR3BP dynamics, equilibrium points, propagation (production code)

M2 implements and verifies the normalized Earth–Moon CR3BP: effective
potential/gradient, equations of motion, Jacobi constant, a numerical
L1/L2/L3 solver, and a generic propagator. Full details:
[DESIGN.md — Milestone 2 section](DESIGN.md#milestone-2--cr3bp-dynamics-equilibrium-points-propagation-and-jacobi-verification).

**No halo orbit exists yet** — this is generic, halo-agnostic CR3BP
infrastructure only.

Equilibrium points, solved numerically (not hardcoded) via `brentq`:

| Point | x (nondim) | Residual | Distance from Moon |
|---|---:|---:|---:|
| L1 | 0.8369151341 | 2.2e-16 | 58,019.138 km |
| L2 | 1.1556821589 | 1.9e-14 | 64,514.906 km |
| L3 | -1.0050626451 | 5.0e-16 | 766,075.396 km |

L2 reproduces the M1 hand value to 10 decimal places. Verification
highlights: Jacobi conservation to `<2e-11` over a 10-day test
trajectory (see figure below), equilibrium states stationary under
propagation to numerical round-off (L2's mild ~1.8e-8 drift over 30
days reflects its known dynamical *instability*, not an error), an
independent from-scratch acceleration cross-check agreeing to
`1.3e-15`, and a 3-tolerance convergence study showing monotonic
improvement. 85 tests pass.

![M2 verification trajectory and Jacobi conservation error](figures/m2_fig2_verification_trajectory_jacobi.png)

*A benign non-equilibrium test trajectory near L2 — explicitly NOT a
halo orbit — used only to exercise the propagator and Jacobi
diagnostic.*

## M3: corrected Earth–Moon L2 halo orbit (production code)

M3 builds a linearized CR3BP seed near L2, differentially corrects it
with a Newton/STM-based shooting method, and verifies the result is a
genuine numerically periodic 3D orbit. Full details, including the
seed derivation, correction formulation, iteration history, and two
genuine bugs found and fixed along the way:
[DESIGN.md — Milestone 3 section](DESIGN.md#milestone-3--l2-halo-orbit-seed-differential-correction-stm-verification-and-periodic-orbit-validation).

**No manifold, arrival state, or insertion Δv is computed yet** — M3
ends at a verified periodic orbit only.

| Quantity | Value |
|---|---:|
| z-amplitude (`z0` at symmetry crossing) | 0.035 DU = 13,454 km |
| max \|z\| over full orbit | 19,426.7 km |
| Period | 3.394396 DU-time = **14.74 days** |
| Full-period closure (generic propagator) | 1.92e-5 km position, 1.24e-7 m/s velocity |
| Jacobi drift over one period | 4.06e-12 (nondimensional) |
| Half-period symmetry residual | ~1e-12 to 1e-14 (round-off level) |
| Monodromy dominant eigenvalue | 999.5 → linearly **unstable** (expected; not a stationkeeping result) |

![M3 corrected Earth-Moon L2 halo orbit](figures/m3_fig1_halo_orbit_3d.png)

*Numerically differential-corrected periodic orbit (southern-convention,
`z0>0`); not an ephemeris trajectory. See DESIGN.md for the x-y/x-z/y-z
projections and the convergence/Jacobi verification figures.*

## M4: computed insertion Δv estimate (headline result)

> ## **Δv = 109.6 m/s**
>
> **Arrival-state assumption:** Jacobi-consistent arrival speed
> (`v_arr² = 2·Ω(r_h) − C_arr`, `C_arr = C_halo − 0.01` nondim) +
> velocity direction radial from Earth `(-μ,0,0)` to the insertion
> point. This is **not** a propagated Earth-departure trajectory — it
> is an explicit, documented engineering assumption. A different
> (equally defensible) assumption gives a different number — see the
> sensitivity studies below and [DESIGN.md — Milestone 4 section](DESIGN.md#milestone-4--arrival-state-definition-insertion-point-trade-and-halo-insertion-δv-estimate)
> for the full derivation.

![M4 insertion Delta-v vs. orbital phase](figures/m4_fig1_delta_v_vs_phase.png)

Selected insertion point: `tau = 0.2441` (phase along the M3 halo
period), ~3.60 days after the M3 reference crossing, 67,903 km from the
Moon, 38,303 km from L2. A secondary local minimum exists at
`tau=0.781` (224.3 m/s) — reported, not discarded. An idealized
aligned-direction lower bound (nonphysical, direction exactly matching
`v_halo`) gives 42.3 m/s — shown only as a sanity bound, never as a
meaningful insertion solution (see DESIGN.md's degeneracy-guard
discussion).

**Sensitivity:** varying the arrival-speed assumption (`dC` = 0.002 to
0.05 nondim) moves the minimum Δv from 89.5 to 199.0 m/s; a ±20°
arrival-direction perturbation moves it only ~12 m/s (109.6→100.4/112.7
m/s) — direction assumption reshapes *where* the optimum sits more than
*how deep* it is. The phase optimum itself is broad, not razor-thin:
Δv grows only ~20 m/s over ±17 hours around the selected phase.

### Comparison with M1's preliminary estimate

M1 (Section 7, hand estimate, no propagated states): **O(10¹–10²) m/s
(10–200 m/s)**. M4 (this section, computed from the actual corrected
halo orbit + an explicit arrival assumption): **109.6 m/s** — falls
within M1's range. This agreement is a genuine outcome of the chosen
assumption, not tuned to match M1 (DESIGN.md's M4 section documents
this explicitly).

**This is a local velocity-matching insertion estimate at the
corrected halo orbit under a stated arrival assumption — it is not an
optimized Earth-to-L2 transfer trajectory.**

## M5: dynamically propagated arrival-arc validation

**Did M4 survive validation? Yes.**

| | Δv |
|---|---:|
| **M4** (local arrival-state assumption) | **109.5723 m/s** |
| **M5** (dynamically propagated CR3BP arrival arc) | **109.5723 m/s** |
| Difference | 4.7×10⁻⁹ m/s (4.3×10⁻⁹ %) |

M5 takes M4's exact candidate state, propagates it **backward** as a
real CR3BP trajectory, checks that it genuinely reaches toward Earth
(`< 0.8 DU`, the threshold set *before* searching) while clearing the
Moon by a conservative 5-lunar-radii guard, then **forward**-propagates
the far end to independently recover the arrival velocity. The
recovered Δv matches M4 to round-trip numerical precision —
**classification: robust first-order estimate.**

![M4 vs M5 insertion estimate](figures/m5_fig2_m4_vs_m5.png)

**A genuine finding, not hidden:** M4's *secondary* local minimum
(`tau≈0.781`, 224.3 m/s) does **not** survive this validation — its
backward arc never gets closer than 0.959 DU to Earth at any tested
duration (10–120 days), well short of the 0.8 DU threshold. M4 alone
could not tell these two branches apart; M5's dynamical check can.

![M5 arrival arc](figures/m5_fig1_arrival_arc.png)

*The backward-propagated arrival arc (purple) swings past L1 and near
the Moon before heading toward Earth — a real CR3BP trajectory, not a
targeted Earth departure. See DESIGN.md for the full search domain,
sensitivity studies, and convergence results.*

**This validates M4's estimate under this project's CR3BP, impulsive-
burn fidelity. It does not optimize or claim a complete Earth-to-L2
transfer, and the project is not yet portfolio-ready (M6 remains).**

## Illustrative propellant cost (6,000 kg spacecraft)

M1's original order-of-magnitude table (kept for reference), plus the
M4 computed value:

| Δv (m/s) | Isp=320s | Isp=450s |
|---:|---:|---:|
| 25 | 47.6 kg | 33.9 kg |
| 50 | 94.8 kg | 67.6 kg |
| 100 | 188.2 kg | 134.4 kg |
| 150 | 280.1 kg | 200.5 kg |
| **109.6 (M4/M5, final)** | **205.9 kg** | **147.1 kg** |

Illustrative only — no propulsion architecture is finalized. Excludes
TLI, MCC, stationkeeping, and launch Δv throughout. (M4 and M5
propellant values agree to ~1e-8 kg — see M5 section above.)

## Milestone roadmap

| Milestone | Scope |
|---|---|
| **M1 ✅** | Mission definition, CR3BP derivation, L2 calculation, insertion-Δv definition, hand estimates |
| **M2 ✅** | CR3BP integrator, equilibrium-point solver, Jacobi conservation verification |
| **M3 ✅** | Linearized halo seed + STM-based differential correction to a genuine periodic L2 halo orbit |
| **M4 ✅** | Arrival-state model + insertion-point Δv calculation + sensitivity study |
| **M5 ✅** | Dynamically propagated CR3BP arrival-arc validation of the M4 estimate |
| M6 (next) | Independent validation, final figures, portfolio packaging (not yet done) |

## Limitations (see [DESIGN.md §12](DESIGN.md#12-limitations) for the full list)

CR3BP only — circular/coplanar Earth–Moon orbits, no Sun perturbation,
no lunar eccentricity, no solar radiation pressure, no ephemeris model,
no finite-burn modeling, no launch/transfer optimization, no
stationkeeping design, no navigation/dispersion or covariance analysis.
**This is an engineering-approximation estimate, not a flight
trajectory solution.**

## Running tests

```bash
pip install -e ".[dev]"
pytest -W error
```

168 tests pass as of M5: everything from M2-M4 (constants/normalization,
CR3BP potential/gradient/RHS validation, L1/L2/L3 equilibrium residuals,
Jacobi-constant identities, propagation/convergence, variational/STM,
halo seed, differential-correction, full-orbit periodicity/symmetry/
monodromy, arrival-model, phase-sweep/refinement tests), plus M5's
backward-arrival-arc tests (round-trip closure, Earthward/Moon-guard
criteria, propellant, no-regression). No manifold/transfer-
optimization/TLI/stationkeeping tests exist yet — those begin at M6, if
ever in scope.

## Repository layout

```
DESIGN.md                       full derivations, equations, assumptions
README.md                       this file
src/halo_insertion/             package: constants, normalization, CR3BP
                                 dynamics, equilibria, propagation (M2);
                                 variational/STM, halo seed, differential
                                 correction, halo orchestration (M3);
                                 arrival-state model, insertion phase
                                 sweep/refinement (M4); backward arrival-
                                 arc propagation/round-trip check (M5)
tests/                          pytest suite (168 tests as of M5)
scripts/make_m2_figures.py      generates the M2 diagnostic figures
scripts/build_m3_halo.py        builds/corrects/validates the M3 halo orbit,
                                 writes results/m3_*.{csv,json}
scripts/make_m3_figures.py      generates the M3 diagnostic figures
scripts/build_m4_insertion.py   phase sweep, refinement, sensitivity,
                                 writes results/m4_*.{csv,json}
scripts/make_m4_figures.py      generates the M4 diagnostic figures
scripts/build_m5_validation.py  backward-arc search, round-trip check,
                                 sensitivity, convergence, writes results/m5_*
scripts/make_m5_figures.py      generates the M5 diagnostic figures
figures/                        M2 + M3 + M4 + M5 diagnostic figures
results/                        m3_halo_*.{csv,json}, m4_phase_sweep.csv,
                                 m4_sensitivity.csv, m4_insertion_summary.json,
                                 m5_arrival_arc_search.csv, m5_sensitivity.csv,
                                 m5_summary.json
```
