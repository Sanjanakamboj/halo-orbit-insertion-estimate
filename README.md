# Halo-Orbit Insertion Δv Estimate

**Portfolio deliverable:** Halo-orbit insertion Δv estimate for an
Earth–Moon L2 telescope-class mission, with clearly stated assumptions
and a concrete verification plan.

**Status: Milestone 1 (M1) complete. M2 (CR3BP integrator + equilibrium
solver + Jacobi verification) not started.**

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

The target orbit is currently a **halo-orbit seed** (literature-scale
placeholder amplitude/period, southern L2 family) — **not** a
numerically corrected periodic orbit. That correction step is M3.
See [DESIGN.md §2](DESIGN.md#2-reference-halo-family-member-m1-seed).

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

## Preliminary insertion-Δv expectation (**M1, not yet computed — a regression bound only**)

Two independent M1 hand estimates (literature order-of-magnitude +
local velocity-mismatch scaling) both point to:

> **O(10¹–10²) m/s — tens to low hundreds of m/s, not km/s.**

This is **not** a final answer. It exists so that M4's actual computed
Δv (from real propagated/corrected states) has something to be sanity
checked against. See [DESIGN.md §7](DESIGN.md#7-definition-of-insertion-δv-for-this-project).

## Illustrative propellant cost (6,000 kg spacecraft)

| Δv (m/s) | Isp=320s | Isp=450s |
|---:|---:|---:|
| 25 | 47.6 kg | 33.9 kg |
| 50 | 94.8 kg | 67.6 kg |
| 100 | 188.2 kg | 134.4 kg |
| 150 | 280.1 kg | 200.5 kg |

Illustrative only — no propulsion architecture is finalized.

## Milestone roadmap

| Milestone | Scope |
|---|---|
| **M1 ✅** | Mission definition, CR3BP derivation, L2 calculation, insertion-Δv definition, hand estimates |
| M2 | CR3BP integrator, equilibrium-point solver, Jacobi conservation verification |
| M3 | Richardson analytical halo seed + differential correction to a periodic orbit |
| M4 | Arrival-state model + insertion-point Δv calculation + sensitivity study |
| M5 | Trade study: halo size, insertion point, transfer geometry, propellant implications |
| M6 | Independent validation, convergence study, final figures, portfolio packaging |

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

M1 ships one placeholder test module confirming package metadata and
that no numerical CR3BP solver is falsely claimed to exist yet.

## Repository layout

```
DESIGN.md                    full derivations, equations, assumptions
README.md                    this file
src/halo_insertion/          package (metadata only at M1)
tests/                       pytest suite
scripts/                     analysis scripts (empty at M1)
figures/                     generated figures (empty at M1)
results/                     numerical results (empty at M1)
```
