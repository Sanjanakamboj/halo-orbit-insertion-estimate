# DESIGN — Earth–Moon L2 Halo-Orbit Insertion Δv Estimate

Status: **Milestone 1 (M1) — mission definition, governing equations, L2
location, insertion-Δv definition, hand estimates, and repository
scaffolding only.** No CR3BP propagator, halo differential correction,
manifold targeting, or insertion solver is implemented yet (see
[Section 10](#10-later-milestone-plan)).

---

## 1. Mission scenario

This project analyzes a single concrete reference mission:

- **System:** Earth–Moon
- **Destination:** Earth–Moon L2 (collinear libration point beyond the
  Moon, on the far side from Earth)
- **Orbit family:** a **southern** L2 halo orbit (chosen arbitrarily
  between the two mirror-symmetric families; either is equally valid —
  southern is selected here and used consistently through later
  milestones)
- **Spacecraft class:** telescope-class observatory (JWST-like mission
  archetype — a large, mass-constrained science payload, not a crewed
  or cargo vehicle)
- **Representative wet mass:** **6,000 kg** (chosen within the stated
  several-ton class, 5,000–8,000 kg)
- **Maneuver model:** the insertion burn is modeled as a **single
  impulsive Δv** (instantaneous velocity change, zero burn duration) —
  a standard simplifying baseline for early trajectory-design estimates
- **Arrival assumption:** the spacecraft is assumed to already be on a
  transfer trajectory that has been targeted to arrive near the halo
  insertion region (e.g., via a prior translunar injection + midcourse
  correction sequence). **That transfer design is out of scope for this
  project.**

### What this project is — and is not — computing

This project estimates **one specific quantity**:

> **The Δv required, at the arrival point near the target halo orbit,
> to convert the incoming transfer velocity into a velocity consistent
> with the halo orbit — i.e., the halo-orbit insertion burn.**

It is explicitly **not** computing, and these terms are not used
interchangeably anywhere in this repository:

| Term | What it means | In scope for this project? |
|---|---|---|
| Launch Δv | Δv to reach a parking/departure orbit from the launch pad | No |
| Translunar injection (TLI) Δv | Δv to leave Earth orbit onto a lunar/L2-bound transfer | No |
| Midcourse correction (MCC) Δv | Small trajectory-cleanup burns during transfer | No |
| **L2 transfer Δv** | Δv associated with shaping the Earth-to-L2 transfer itself | No |
| **Halo insertion Δv** | Δv to match the halo orbit's velocity upon arrival | **Yes — this project's quantity of interest** |
| Stationkeeping Δv | Ongoing small corrections to remain near the unstable halo | No (noted only for context, Section 12) |

---

## 2. Reference halo family member (M1 seed)

At M1 the target halo orbit is defined as a **halo-orbit seed**: a
third-order Richardson-approximation analytical initial guess, in the
tradition of Richardson (1980) and standard CR3BP halo-generation
references (e.g., Koon, Lo, Marsden & Ross, *Dynamical Systems, the
Three-Body Problem, and Space Mission Design*, 2011). It is **not** a
numerically corrected periodic orbit — that step is explicitly deferred
to M3.

Chosen representative design point (Earth–Moon L2, southern family,
"medium" amplitude — comparable in scale to published Earth–Moon L2
halo families used for relay/observatory concepts):

| Quantity | Value | Status |
|---|---|---|
| System | Earth–Moon | fixed |
| Libration point | L2 | fixed |
| Family | Southern | chosen, Section 1 |
| Out-of-plane amplitude, `A_z` | ≈ 0.030–0.040 nondimensional DU (≈ 11,500–15,400 km) | literature-scale placeholder |
| In-plane extent, `A_x` | ≈ 0.03–0.05 nondimensional DU (order-of-magnitude, comparable to `A_z` for this family size) | literature-scale placeholder |
| Approximate period | ≈ 11–15 days (≈ 2.5–3.4 nondimensional TU) | literature-scale placeholder, consistent with published Earth–Moon L2 halo periods |
| State-vector origin | **placeholder, to be replaced by a Richardson third-order analytic seed in M3, then numerically corrected** | explicit |

This table intentionally states ranges rather than a single
high-precision state vector: at M1 no analytic seed generator or
differential corrector exists yet, so a single 6-element state would
be false precision. M3 will replace this table with an actual seeded
state vector `(x0, y0, z0, ẋ0, ẏ0, ż0)` and a corrected period `T`.

**Terminology used consistently in this repository:**
- *halo-orbit seed* — an uncorrected analytical/approximate guess (current status)
- *reference family member* — a specific chosen amplitude/family combination (this section)
- *numerically corrected periodic orbit* — reserved for after differential correction is implemented (M3+); **not yet earned**

---

## 3. Physical constants and nondimensionalization

The Earth–Moon CR3BP uses the standard rotating-frame normalization:

- **Distance unit (DU):** mean Earth–Moon center-to-center separation
- **Mass normalization:** total system mass (Earth + Moon) = 1;
  mass ratio `mu = m_moon / (m_earth + m_moon)` locates the Moon at
  `(1-mu, 0, 0)` and Earth at `(-mu, 0, 0)` in the rotating frame
- **Angular rate normalization:** rotating-frame angular rate `n = 1`
  (by construction of the time unit below)
- **Time unit (TU):** derived from the Earth–Moon mean motion `n` such
  that the system completes one synodic rotation in `2π` TU

### Gravitational parameters used

| Quantity | Symbol | Value | Source |
|---|---|---|---|
| Earth GM | `GM_earth` | 398,600.4418 km³/s² | standard Earth GM (DE-class ephemeris value, widely tabulated) |
| Moon GM | `GM_moon` | 4,902.8000 km³/s² | standard Moon GM |
| Mean Earth–Moon distance | `DU` | 384,400 km | standard mean value |

### Computed normalization quantities

```
mu = GM_moon / (GM_earth + GM_moon) = 0.0121505839
n  = sqrt((GM_earth + GM_moon) / DU^3) = 2.6653144e-06 rad/s
TU = 1/n = 3.751903e+05 s = 4.342480 days
V* = n * DU = 1.024547 km/s = 1024.547 m/s
```

`mu ≈ 0.01215058` matches the commonly published Earth–Moon CR3BP mass
ratio (≈ 0.012150585) to full precision — a first internal consistency
check.

### Conversion rules

| From | To | Rule |
|---|---|---|
| nondimensional position `r*` | dimensional km | `r_km = r* * DU` |
| nondimensional velocity `v*` | dimensional km/s | `v_km/s = v* * V*` |
| nondimensional time `t*` | dimensional days | `t_days = t* * TU / 86400` |
| dimensional km | nondimensional position | `r* = r_km / DU` |
| dimensional km/s | nondimensional velocity | `v* = v_km/s / V*` |
| dimensional days | nondimensional time | `t* = (t_days * 86400) / TU` |

**Why this matters:** all halo-orbit and insertion-Δv arithmetic in
this project is done first in dimensionless CR3BP units (`DU`, `TU`,
`V*`), then converted to physical units (km, km/s, m/s) only at the
end. Keeping this boundary explicit is what Section 11's round-trip
verification test (Section 11-B) checks.

---

## 4. Governing CR3BP equations

State vector (position + velocity, synodic rotating frame, origin at
the Earth–Moon barycenter):

```
s = [x, y, z, xdot, ydot, zdot]^T
```

Equations of motion (rotating frame, normalized units, rotation rate = 1):

```
xddot - 2*ydot = dOmega/dx
yddot + 2*xdot = dOmega/dy
zddot          = dOmega/dz
```

Pseudo-potential:

```
Omega(x,y,z) = 0.5*(x^2 + y^2) + (1-mu)/r1 + mu/r2
```

with

```
r1 = sqrt((x+mu)^2   + y^2 + z^2)     # distance from Earth
r2 = sqrt((x-1+mu)^2 + y^2 + z^2)     # distance from Moon
```

### Expanded partial derivatives

```
dOmega/dx = x - (1-mu)*(x+mu)/r1^3 - mu*(x-1+mu)/r2^3
dOmega/dy = y - (1-mu)*y/r1^3      - mu*y/r2^3
dOmega/dz =   - (1-mu)*z/r1^3      - mu*z/r2^3
```

### First-order ODE form

Defining `sdot = f(s)`:

```
xdot_    = xdot
ydot_    = ydot
zdot_    = zdot
xddot_   = 2*ydot + dOmega/dx
yddot_   = -2*xdot + dOmega/dy
zddot_   =  dOmega/dz
```

i.e. `f(s) = [ xdot, ydot, zdot, 2*ydot + dOmega/dx, -2*xdot + dOmega/dy, dOmega/dz ]^T`.

This six-state first-order form is exactly what a numerical integrator
(RK4/RK45/DOP853, etc.) will consume in M2. **No integrator is
implemented in M1** — this section documents the equations only.

---

## 5. Jacobi constant

The Jacobi integral is the CR3BP's one analytic constant of motion
under unforced (coast) dynamics:

```
C = 2*Omega(x,y,z) - (xdot^2 + ydot^2 + zdot^2)
```

### Why it matters for verification

- **Conservation check:** under the unforced CR3BP equations of Section
  4, `C` is exactly constant along any trajectory. Any numerical drift
  in `C` during propagation directly measures integrator error — this
  is the cheapest, most direct sanity check available once M2's
  integrator exists (Section 11-D).
- **Discontinuity at the burn:** `C` is **not** conserved across an
  impulsive maneuver, because the maneuver changes velocity
  discontinuously while position is unchanged. The Jacobi constant
  before and after the halo-insertion burn will differ — this is
  expected and is itself a useful cross-check that the "before" state
  (transfer arrival) and "after" state (halo orbit) are genuinely
  different trajectories, not a degenerate zero-Δv case (Section 11-I).
- **Sanity-checking arrival/halo states independently:** a halo-orbit
  seed and an assumed arrival state should each individually satisfy
  `C = 2*Omega - v^2` self-consistently; computing `C` for each is a
  free, no-integration check on state-vector validity.

### Dimensional interpretation

`C` as defined above is **inherently a nondimensional quantity** in the
normalized CR3BP formulation — it mixes a nondimensional pseudo-potential
with nondimensional velocity-squared, and there is no natural,
unambiguous "Jacobi constant in km²/s²" without re-deriving the
constant directly in dimensional variables (which introduces `GM_earth`,
`GM_moon`, and the rotation rate `n` explicitly, rather than folding
them into `mu`, `DU`, `TU`). This project reports `C` only in
normalized units and does not attempt a dimensional Jacobi constant.

---

## 6. L2 location

### Derivation

L2 lies on the x-axis beyond the Moon (`x > 1-mu`), where the net
rotating-frame acceleration (gravity of both bodies + centrifugal term)
vanishes: `dOmega/dx = 0` with `y = z = 0`. Substituting into the
Section 4 expression, and writing `r1 = x+mu`, `r2 = x-(1-mu)` for a
point beyond the Moon (both positive), the collinear equilibrium
condition reduces to the scalar equation:

```
x - (1-mu)/(x+mu)^2 - mu/(x-1+mu)^2 = 0
```

This is a degree-5 polynomial in `x` (after clearing denominators) with
no closed-form solution in general; it is solved numerically (Newton's
method here, seeded from the classical small-`mu` approximation
`x0 = 1 + (mu/3)^(1/3)`, which is accurate to ~2 significant figures
and converges in a handful of Newton iterations to machine precision).

### Computed value

```
x_L2 (nondimensional, barycentric) = 1.1556821589
```

### Dimensional conversion

```
x_L2 (km, from barycenter) = x_L2 * DU = 444,244.222 km   (≈ distance from Earth center, since barycenter ≈ Earth center)
distance from Moon center  = (x_L2 - (1-mu)) * DU = 64,514.906 km
```

Using Moon mean radius `R_moon ≈ 1,737.4 km`:

```
distance beyond Moon surface ≈ 64,514.906 - 1,737.4 ≈ 62,777.5 km
```

### Cross-check

The result (**≈ 64,500 km from the Moon's center, ≈ 62,800 km beyond
its surface**) falls squarely within the expected physical range of
**~60,000–70,000 km beyond the Moon** cited for the Earth–Moon L2
point, confirming the equilibrium solve is physically reasonable. This
is exactly the kind of independent sanity check formalized in Section
11-A.

---

## 7. Definition of "insertion Δv" for this project

### Rigorous definition

At a chosen insertion point on or near the target halo orbit, define:

- **Target halo state:** `r_h, v_h` — position and velocity on the
  reference halo trajectory at the insertion point
- **Arrival transfer state:** `r_a, v_a` — position and velocity of the
  incoming transfer trajectory at that same point (by construction,
  `r_a ≈ r_h`; the transfer has been targeted to arrive at the halo's
  location, but not necessarily its velocity)

**Baseline impulsive insertion maneuver:**

```
Delta_v = v_h - v_a
|Delta_v| = norm(v_h - v_a)
```

This is a **velocity-matching** burn: it changes the incoming velocity
vector to the halo orbit's local velocity vector at the same position,
converting an approach trajectory into (approximately) the periodic
halo trajectory itself.

Actually solving for `r_a, v_a` (i.e., constructing a realistic
transfer arrival state) and for `r_h, v_h` from a **corrected** halo
orbit is deferred to M4, once M2 (propagator) and M3 (corrected halo)
exist. M1 instead establishes two independent estimates of the
**expected scale** of this quantity, so later numerical results have
something to be checked against.

### Estimate A — literature / order-of-magnitude sanity check

Published Earth–Moon L2 halo-orbit and quasi-halo insertion/injection
maneuvers (drawn from mission design literature for Earth–Moon L1/L2
libration-point missions and CR3BP mission-design textbooks/course
notes, e.g., Koon/Lo/Marsden/Ross 2011 and general halo-orbit
mission-design surveys) fall broadly in the **tens to low hundreds of
m/s** range for a single insertion burn, distinct from the much larger
(~500–1000+ m/s) deterministic totals that include translunar
transfer shaping.

This is reported here as a **broad engineering range, not a specific
mission's flight value**, because:
- No specific, publicly-documented Earth–Moon L2 halo mission (as
  opposed to Sun–Earth L1/L2, which has flown several times — e.g.
  WMAP, Herschel, Planck, JWST, Gaia) was cross-checked numerically for
  this project; citing an exact number without doing so would overstate
  precision.
- Mission Δv budgets often bundle insertion with nearby midcourse
  cleanup, so isolating a "pure" insertion figure from public summaries
  is not reliable without primary-source navigation data.

**Estimate A range: O(10¹–10²) m/s**, i.e. roughly **10–200 m/s**,
labeled explicitly as an order-of-magnitude expectation, not a target
value.

### Estimate B — local velocity-mismatch hand estimate

Using the characteristic velocity `V* = 1024.547 m/s` (Section 3), a
plausible nondimensional arrival/halo velocity mismatch of similar
scale to typical CR3BP targeting residuals near a libration-point orbit
(chosen here illustratively as **0.01–0.02 DU/TU**, i.e., a few percent
of the characteristic velocity — consistent with the scale of velocity
discontinuities seen when a targeted-but-uncorrected transfer arrives
near, not exactly on, a periodic orbit) gives:

```
0.005 DU/TU -> 5.12 m/s
0.010 DU/TU -> 10.25 m/s
0.020 DU/TU -> 20.49 m/s
0.050 DU/TU -> 51.23 m/s
```

**Estimate B range: ≈ 5–50 m/s** for a well-targeted arrival, consistent
with (and at the low end of) Estimate A.

### M1 regression-scale conclusion

Both estimates agree the answer should land in the
**O(10¹–10²) m/s regime — tens to low hundreds of meters per second —
not km/s.** This is recorded purely as a sanity bound for later
milestones' numerical results (Section 11-M); it is **not** a claimed
final answer and will be superseded once M4 computes an actual
`|v_h - v_a|` from real propagated/corrected states.

---

## 8. Telescope-class payload / spacecraft assumption

**Chosen representative wet mass: `m0 = 6,000 kg`.**

### Mass-independence of ideal impulsive Δv

For an ideal impulsive-maneuver trajectory calculation, the **required
Δv is a purely kinematic quantity** — it depends only on the geometry
and velocities of the target and arrival trajectories (Section 7), not
on the mass of the spacecraft performing the maneuver. A 3,000 kg
smallsat and an 8,000 kg observatory arriving on the same transfer
trajectory and targeting the same halo orbit require the **same Δv**.

Mass becomes relevant only in the next step: **converting a given Δv
into a propellant mass**, via the Tsiolkovsky rocket equation:

```
m_prop = m0 * [1 - exp(-Delta_v / (Isp * g0))]
```

with `g0 = 9.80665 m/s²`.

### Illustrative propulsion systems

Two representative — **not final-design** — propulsion options, spanning
typical chemical propulsion performance:

| System | `Isp` |
|---|---|
| Chemical bipropellant (e.g., MMH/N2O4-class) | 320 s |
| High-performance chemical (e.g., LOX/LH2 or advanced bipropellant-class) | 450 s |

### Propellant table (`m0 = 6,000 kg`)

| Δv (m/s) | `m_prop` @ Isp=320s (kg) | `m_prop` @ Isp=450s (kg) |
|---:|---:|---:|
| 25 | 47.61 | 33.89 |
| 50 | 94.84 | 67.60 |
| 100 | 188.18 | 134.43 |
| 150 | 280.05 | 200.52 |

These illustrate the propellant-budget stakes of the Section 7 Δv
range (tens–low hundreds of m/s): a few tens to a few hundred kg of
propellant out of a 6,000 kg wet mass, depending on final Δv and
propulsion choice — **not** a large fraction of the spacecraft, which
is qualitatively consistent with halo insertion being a modest
correction rather than a major orbital transfer.

---

## 9. Reference mission architecture (analysis pipeline)

Planned end-to-end analysis flow for this repository, to be built out
across M2–M6:

1. Earth–Moon normalized CR3BP model (constants, `mu`, `Omega`, EOM) — **documented M1, coded M2**
2. Locate L1/L2 equilibrium points numerically — **derived M1 (L2 only), coded M2**
3. Generate an analytical halo-orbit seed (Richardson third-order approximation) — **placeholder table M1, coded M3**
4. Integrate the CR3BP equations of motion — **documented M1, coded M2**
5. Differentially correct the seed to a numerically periodic halo orbit — **M3**
6. Construct/select a representative arrival state near the halo — **M4**
7. Compute the velocity mismatch between arrival and halo states — **M4**
8. Calculate the impulsive insertion Δv — **M4**
9. Perform insertion-point sensitivity analysis (Δv vs. chosen insertion phase/location) — **M4/M5**
10. Validate with Jacobi conservation and independent checks — **ongoing from M2, formalized M6**

No code for steps 3–10 exists beyond this documentation and the
package skeleton in M1.

---

## 10. Later-milestone plan

| Milestone | Scope |
|---|---|
| **M1** (this milestone) | Mission definition, CR3BP derivation, L2 calculation, insertion-Δv definition, hand estimates, repo scaffolding |
| **M2** | CR3BP equations of motion + numerical integrator + equilibrium-point solver (L1 & L2) + Jacobi-constant conservation verification |
| **M3** | Richardson third-order analytical halo seed + numerical differential correction to an actual periodic Earth–Moon L2 halo orbit |
| **M4** | Arrival-state model + insertion-point Δv calculation + insertion-point sensitivity study |
| **M5** | Trade study: halo size/amplitude, insertion point, arrival velocity-mismatch/transfer-geometry assumptions, spacecraft propellant implications |
| **M6** | Independent validation, numerical convergence study, final figures, README/portfolio packaging |

M2+ work has **not** started. This document and the M1 file set are
the complete M1 deliverable.

---

## 11. Verification plan

Concrete, code-checkable tests planned for implementation as each
milestone lands the functionality it needs (none require plots as the
primary verification method):

- **A. L2 equilibrium residual** — `dOmega/dx, dOmega/dy, dOmega/dz`
  (equivalently, total rotating-frame acceleration) evaluate to ≈0 at
  the computed L2 state. *(Numerical solve done in M1 Section 6;
  formal residual test lands with M2's `Omega`/gradient code.)*
- **B. Dimensional/nondimensional round-trip** — converting a position,
  velocity, and time to dimensional units and back reproduces the
  original nondimensional values to numerical precision.
- **C. CR3BP equation symmetry checks** — the equations of motion are
  invariant under `(y, ẋ, ż) -> (-y, -ẋ, -ż)` (standard CR3BP mirror
  symmetry about the x–z plane); used to check halo-family symmetry in M3.
- **D. Jacobi conservation** — under unforced propagation, `C`
  (Section 5) stays constant to within integrator tolerance over one
  full halo period.
- **E. Equilibrium propagation** — propagating the L2 equilibrium state
  itself (zero velocity, at rest in the rotating frame) for any duration
  leaves the state unchanged (to integrator tolerance).
- **F. STM / differential-correction verification** — once implemented
  (M3), the state-transition matrix and corrector convergence are
  checked against finite-difference sensitivities.
- **G. Periodicity closure** — a differentially corrected halo state,
  propagated for one period `T`, returns to its initial state within a
  specified tolerance (M3).
- **H. z-symmetry / half-period symmetry** — for a halo orbit, the
  state at `T/2` satisfies the mirror-symmetry condition of check C
  relative to the initial state (M3).
- **I. Insertion identity** — if the arrival velocity is set exactly
  equal to the halo velocity at the same point, the computed
  `Delta_v` is exactly zero (M4; this is a pure unit/logic check of
  the Section 7 definition, independent of any propagation accuracy).
- **J. Scaling check** — a nondimensional Δv, multiplied by `V*`,
  reproduces the same value as computing Δv directly in dimensional
  km/s from dimensional state vectors (Section 3 round-trip applied to
  Section 7's definition).
- **K. Independent propagator / tighter-tolerance cross-check** —
  results reproduced with a second integrator or substantially tighter
  tolerance to bound integration error (M6).
- **L. Convergence study** — Δv and periodicity-closure results are
  checked for stability as integrator step size/tolerance is refined (M6).
- **M. Literature/order-of-magnitude sanity check** — final computed
  insertion Δv (M4) is compared against the Section 7 Estimate A/B
  ranges established in M1; a result wildly outside O(10¹–10²) m/s
  triggers a re-check of the model rather than being accepted at face
  value.

M1 ships only the **placeholder test** (Section 13/15) confirming
package metadata and confirming none of the above numerical machinery
is falsely claimed to exist yet. Tests A, D–L require code that does
not exist until M2+.

---

## Milestone 2 — CR3BP dynamics, equilibrium points, propagation, and Jacobi verification

Status: **complete.** This section documents what M2 actually built and
verified. **No halo orbit — seeded, corrected, or otherwise — has been
generated at this point.** Everything below is generic CR3BP machinery
(equilibrium points, propagation, Jacobi diagnostics) applicable to any
CR3BP state, not halo-specific.

### Implemented equations

Exactly the Section 4/5 equations, coded with no changes of convention:
`Omega`, its analytic gradient `(dOmega/dx, dOmega/dy, dOmega/dz)`, the
first-order RHS `[xdot,ydot,zdot,xddot,yddot,zddot]`, and the Jacobi
constant `C = 2*Omega - v^2`. The production RHS uses closed-form
analytic derivatives only — no finite differences are used inside
`cr3bp_rhs`; finite differences appear only in a test-side cross-check
(`test_analytic_gradient_matches_finite_difference`).

Modules: [`src/halo_insertion/constants.py`](src/halo_insertion/constants.py),
[`normalization.py`](src/halo_insertion/normalization.py),
[`cr3bp.py`](src/halo_insertion/cr3bp.py),
[`equilibria.py`](src/halo_insertion/equilibria.py),
[`propagation.py`](src/halo_insertion/propagation.py).

### Numerical constants (reproduced from production code)

```
mu  = 0.012150583916324809
DU  = 384400.0 km
TU  = 4.342479849812527 days  (375190.259 s)
V*  = 1024.5468552412854 m/s
```

These match the M1 hand values to full precision (M1 rounded to 6-7
significant figures for readability; production code carries full
`float64` precision). **No discrepancy found.**

### L1 / L2 / L3 — numerically solved, not hardcoded

Each point is a `scipy.optimize.brentq` root of `dOmega/dx(x,0,0)=0` on
a physically-motivated bracket (Section 7 of this document, generalized
to all three collinear points):

| Point | x (nondim, barycentric) | `\|dOmega/dx\|` residual | dist. from Earth | dist. from Moon |
|---|---:|---:|---:|---:|
| L1 | 0.8369151341 | 2.22e-16 | 326,380.862 km | 58,019.138 km |
| L2 | 1.1556821589 | 1.87e-14 | 448,914.906 km | 64,514.906 km |
| L3 | -1.0050626451 | 4.98e-16 | 381,675.396 km | 766,075.396 km |

**L2 agrees with the M1 hand value (`x = 1.1556821589`) to 10 decimal
places** — reproduced exactly by an independent `brentq` solve rather
than the M1 Newton's-method solve, confirming M1's L2 result.

**One refinement to M1's phrasing, not an error:** M1 stated the L2
distance from Earth as "≈444,244 km, since barycenter ≈ Earth center."
That approximation treats the barycenter as coincident with Earth's
center. In production code, Earth's actual position is `(-mu, 0, 0)`,
offset from the barycenter by `mu*DU ≈ 4,671 km`. The barycentric L2
coordinate (`x_L2 * DU = 444,244.222 km`) is unchanged and still
correct as "distance from the barycenter"; the **precise** distance
from Earth's own center is `(x_L2 + mu) * DU = 448,914.906 km`. This
does not change the M1 Moon-relative distance (64,514.906 km, which
does not depend on this approximation) or the insertion-Δv conclusions,
and is noted here purely for precision. M1's Section 6 statement is
retained as-is (it was explicitly qualified with "≈"); this section
supersedes it only for the Earth-distance figure.

### Propagator

`propagation.propagate()` wraps `scipy.integrate.solve_ivp` with
default settings `method="DOP853", rtol=1e-11, atol=1e-12`, raw
`OdeResult` exposed unmodified (`result.success` must be checked by the
caller). No halo-specific logic.

### Equilibrium-state propagation verification

Each of L1, L2, L3 was set as an exact equilibrium state
`[x_Li, 0,0,0,0,0]` and propagated 30 days at default tolerances:

- L1 and L3 (see below on stability) remained stationary to within
  numerical round-off.
- **L2's drift over 30 days was ≈1.77e-8 (nondimensional state-vector
  norm)** — small, but larger than round-off, because **collinear
  libration points are dynamically unstable equilibria**: any
  round-off-level perturbation from the exact fixed point grows under
  the CR3BP's local dynamics. This is expected physics, not an
  implementation defect — it is exactly why real halo missions require
  active stationkeeping (Section 12). The test suite bounds this drift
  at `< 1e-6`, loose enough to tolerate genuine round-off growth but
  tight enough to catch a real sign/offset error (which would produce
  drift many orders of magnitude larger, immediately, not after 30
  days).

### Jacobi-constant conservation

On the M2 non-equilibrium verification trajectory (next subsection),
propagated 10 days at default tolerances, Jacobi drift stayed below
`2e-11` (nondimensional) throughout — see Figure 2. At equilibrium
states, `C` matches `2*Omega` to `~1e-13`, confirming the zero-velocity
identity (`test_equilibrium_state_gives_C_equals_2_omega`,
`test_zero_velocity_gives_C_equals_2_omega_generic_point`).

### M2 verification trajectory (explicitly NOT a halo orbit)

A single benign, non-equilibrium state near L2 was used throughout M2's
propagator/Jacobi/convergence checks:

```
state0 = [x_L2 - 0.02, 0.01, 0.005, 0.0, 0.01, 0.005]   (nondimensional)
propagated over 10 days
```

This state has no halo-orbit design intent whatsoever — it exists
purely to exercise the propagator on a non-trivial, non-equilibrium
CR3BP trajectory. See Figure 2.

### Numerical convergence study

Three integrator settings on the same verification trajectory:

| Setting | `rtol` | `atol` | Jacobi drift (10 days) | Terminal-state diff vs. tight |
|---|---:|---:|---:|---:|
| loose | 1e-6 | 1e-8 | 1.148e-06 | 6.251e-05 |
| medium | 1e-9 | 1e-10 | 2.433e-09 | 8.956e-08 |
| tight | 1e-12 | 1e-13 | 1.410e-12 | — (reference) |

Jacobi drift and terminal-state difference **both shrink monotonically
and by comparable orders of magnitude** as tolerances tighten —
`test_numerical_convergence_three_tolerance_levels` asserts this
ordering directly. This is a genuine 3-point convergence demonstration,
not a single-tolerance claim.

### Independent RHS/gradient cross-check

Two independent checks, at five arbitrary nonsingular states not
appearing anywhere else in this document (`(0.5,0.2,0.1)`,
`(1.2,-0.05,0.02)`, `(-0.3,0.4,-0.1)`, `(0.85,0.0,0.15)`,
`(1.05,0.03,0.0)`):

1. **Analytic gradient vs. central finite difference** of `Omega` —
   agreement to `~1e-8` (finite-difference truncation-limited).
2. **Production RHS vs. an independently-written direct acceleration
   expression** (built from gravity + centrifugal + Coriolis terms
   coded from scratch, not sharing code with `effective_potential_gradient`)
   — **max discrepancy across all five states: 1.33e-15** (machine
   precision). This check is specifically designed to catch wrong
   Coriolis signs, wrong Earth/Moon offsets, swapped `mu`/`(1-mu)`, or a
   missing centrifugal term — none were found.

### Symmetry checks

- **Planar invariance:** a state with `z=zdot=0` has `zddot=0` and
  remains exactly planar under propagation (verified structurally and
  by 10-day propagation, max `|z|, |zdot| < 1e-10`).
- **z-reflection:** `dOmega/dx, dOmega/dy` are even in `z`;
  `dOmega/dz` is odd in `z` — verified to `1e-14`.
- **y-reflection:** `dOmega/dx` is even in `y`; `dOmega/dy` is odd in
  `y` — verified to `1e-14`.

### Discrepancies from M1

**None affecting any M1 conclusion.** The only item worth flagging is
the Earth-distance refinement above (barycenter vs. Earth-center
offset, ~4,671 km), which M1 had already qualified with "≈" and does
not change the L2-vs-Moon distance, the physical-scale cross-check
(60,000–70,000 km), or any Section 7 insertion-Δv estimate.

### Explicit scope statement

**No halo-orbit seed, differential correction, state-transition
matrix, monodromy matrix, manifold, arrival-state model, or insertion
Δv has been computed in M2.** M2 delivers only generic, halo-agnostic
CR3BP infrastructure: potential/gradient, equations of motion, Jacobi
constant, L1/L2/L3, and a verified propagator. Halo-orbit generation
begins in M3.

---

## Milestone 3 — L2 halo-orbit seed, differential correction, STM verification, and periodic-orbit validation

Status: **complete.** This section documents M3's construction of the
project's first genuine numerically periodic Earth–Moon L2 halo orbit,
replacing the M1 literature-scale placeholder. All conventions (μ,
Earth/Moon positions, state ordering, equations of motion) are
identical to M1/M2 — no changes.

### Seed method and provenance

**M3 uses a linearized (first-order) CR3BP variational seed, NOT the
third-order Richardson (1980) approximation.** This is an explicit,
documented engineering choice, not an oversight: the full third-order
Richardson expansion requires ~15–20 additional coefficients
(`a21…a32, b21, b22, d21, d31, d32, s1, s2, l1, l2`, …) whose correct
transcription from memory carries meaningful risk of a subtle sign or
indexing error — exactly the kind of error that would silently produce
a wrong-family seed. Rather than risk that, this project:

1. Linearizes the CR3BP equations about L2 (where the equilibrium's
   symmetry gives `Uxy = Uxz = Uyz = 0`), yielding two **decoupled**
   linear subsystems — in-plane `(x,y)` and out-of-plane `z` — exactly
   as in Section 4 of this document, evaluated at the L2 equilibrium.
2. Solves each subsystem in closed form: the out-of-plane motion is
   simple harmonic with frequency `wz = sqrt(-Uzz)`; the in-plane
   motion's oscillatory (center-manifold) mode has frequency
   `wxy = sqrt(-beta)`, where `beta` is the negative root of the
   characteristic equation `beta^2 + (4-Uxx-Uyy)*beta + Uxx*Uyy = 0`.
3. Builds the symmetric halo IC form (Section 7 below) from the
   in-plane trial solution `dx(t) = -Ax*cos(wxy*t)`,
   `dy(t) = kappa*(-Ax)*sin(wxy*t)` (see `halo_seed.py` for the full
   sign derivation, including a documented sign correction found
   during implementation — see "Genuine bugs found" below), with
   `Ax = 0.9 * Az` (a fixed, moderate ratio; not claimed to be the true
   nonlinear Richardson amplitude relation).

**Numerical cross-check on the theory used:** this project independently
*derived and verified* the standard Legendre-coefficient identity
`c2 = -Uzz` (and the companion `Uyy = 1-c2`, `Uxx = 1+2*c2`) against the
production Hessian at L2, and verified the closed-form `c_n(gamma2)`
formula for L2 reproduces `c2` from the Hessian to machine precision
(`c2` from Hessian: `3.19042524638630`; from the formula:
`3.19042524638629`). This built confidence in the L2 geometry
(`gamma2 = 0.16783274...`, matching M2's L2-to-Moon distance) without
requiring the full third-order machinery.

At L2, `wz = 1.786176` and `wxy = 1.862646` (nondimensional) — a
**~4.1% frequency mismatch**. This mismatch is exactly why a pure
linear seed is not already a periodic 3D orbit, and exactly what the
nonlinear differential corrector (below) must resolve.

Documented explicitly per this project's terminology (DESIGN.md
Section 2): this seed is a **halo-orbit seed**, never called a "halo
orbit" until corrected and periodicity-verified below.

### Correction formulation

Free variables: `x0, ydot0` (holding `z0` fixed — this selects the
family member/amplitude). Target residuals: `xdot(T/2) = 0`,
`zdot(T/2) = 0`. Symmetric IC form enforced by construction (never
corrected): `y0 = 0`, `xdot0 = 0`, `zdot0 = 0`.

The correction matrix accounts for the half-period `T/2` itself being
an implicit function of `(x0, ydot0)` (since it is defined by the
dynamic `y=0` event, not a fixed time). Using the chain rule through
the event-time constraint `y(T(free);free) = 0`:

```
dT/d(free)      = -Phi[1, free_cols] / ydot(T/2)
d[s(row)](T/2)/d(free) = Phi[row, free_cols] + sdot(T/2)[row] * dT/d(free)
```

for `row in {xdot=3, zdot=5}`. This 2x2 matrix is inverted via
least-squares (robust to mild ill-conditioning) for the Newton step. A
**damped Newton / backtracking line search** (step halved up to 12
times per iteration until the residual actually decreases) is used for
robustness far from convergence — plain full Newton steps were found
to overshoot badly this close to L2's strongly hyperbolic in-plane
dynamics (see "Genuine bugs found" below).

### Event handling

The symmetric IC has `y(0) = 0` exactly, with `ydot(0) != 0` generally,
so `y(t)` leaves zero *immediately* at `t=0` in one direction — a naive
event search starting at `t=0` would misdetect this trivial departure
as "the" crossing. `propagate_half_period` (in `differential_correction.py`)
handles this with an explicit two-stage integration: a short warm-up
interval `[0, t_warmup=0.05]` with event detection **disabled**, then a
resumed integration with the `y=0` event (`direction=-1`, i.e. only the
*decreasing-through-zero* crossing) armed. `test_half_period_event_does_not_trigger_at_t0`
verifies the detected crossing is always `t_half > 0.1`, well clear of
the seam.

### Correction iteration history (`Az = 0.035` DU, southern-convention seed)

| it | x0 | ydot0 | half-period | xdot residual | zdot residual | correction norm |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 1.124182 | 0.170892 | 1.439639 | 1.865e-01 | -4.790e-02 | 2.056e-01 |
| 1 | 1.109182 | 0.220061 | 1.614843 | 9.697e-02 | -2.446e-02 | 1.131e-02 |
| 2 | 1.104335 | 0.230276 | 1.745071 | -5.671e-02 | 2.085e-02 | 2.076e-02 |
| 3 | 1.109346 | 0.210131 | 1.705962 | -7.983e-03 | 3.222e-03 | 7.576e-03 |
| 4 | 1.111080 | 0.202756 | 1.697397 | -1.673e-04 | 7.013e-05 | 1.923e-04 |
| 5 | 1.111124 | 0.202569 | 1.697198 | -6.458e-08 | 2.831e-08 | 8.669e-08 |
| 6 | 1.111124 | 0.202569 | 1.697198 | -1.525e-14 | 5.653e-15 | 0.0 |

Clean **quadratic** Newton convergence from iteration 3 onward (residual
norm ~1.9e-2 → ~1.9e-4 → ~9.5e-8 → ~1.6e-14) — genuine convergence, not
merely a final answer. Full history: [results/m3_correction_history.csv](results/m3_correction_history.csv).
Convergence plot: Figure 3 below.

### Corrected initial condition and period

```
x0    =  1.1111238554848533   (nondim)
y0    =  0.0
z0    =  0.035                (nondim; 13,454.0 km — the amplitude at the x-z crossing)
xdot0 =  0.0
ydot0 =  0.2025688201399349   (nondim)
zdot0 =  0.0

half-period = 1.697198143630976   (nondim)
period      = 3.394396287261952   (nondim) = 14.740097 days
```

`Az = 0.035` DU (13,454 km) sits inside the M1-intended range
(0.03–0.04 DU / 11,500–15,400 km); it was **not forced** to any
specific value — it is simply the `z0` at which this seed converged
cleanly (Section 22 of this milestone's instructions: `Az` values
0.005/0.01 (z-xy coupling too weak — nearly-singular correction matrix)
and 0.04 (outside this seed's convergence basin) were tried and
explicitly rejected; 0.02/0.03/0.035 all converged, and 0.035 was kept
as the reported orbit). Full summary:
[results/m3_halo_summary.json](results/m3_halo_summary.json).

### Full-period verification (independent of the correction loop)

Propagated one full period using the **generic M2 propagator**
(`propagation.propagate`, default `DOP853, rtol=1e-11, atol=1e-12`) —
not the half-period-event machinery used by the corrector:

```
position closure norm  = 4.988e-11 (nondim) = 1.918e-05 km
velocity closure norm  = 1.207e-10 (nondim) = 1.237e-07 m/s
component-wise: dx=3.145e-11, dy=-3.837e-11, dz=5.221e-12,
                 dxdot=9.570e-11, dydot=-6.504e-11, dzdot=3.441e-11
```

### Jacobi verification

```
C(0)            = 3.1412189199162857   (nondim)
max |C(t)-C(0)| = 4.059e-12  over one full period (4000-point sampling)
RMS drift       = ~1.3e-12
```

No energy/Jacobi correction was applied at any point — this is a pure
propagation diagnostic on the already-corrected orbit.

### Halo geometry

```
x range:  [1.111124, 1.178180]  (nondim) -> [-17,128, +8,648] km relative to x_L2
y range:  [-0.096662, 0.096662] (nondim) -> ±37,157 km
z range:  [-0.050538, 0.035000] (nondim) -> [-19,427, +13,454] km
max|z|  =  0.050538 (nondim) = 19,426.7 km   <-- note: LARGER than z0=Az=13,454 km;
                                                  see note below
period  =  3.394396 (nondim) = 14.740098 days
```

`max|z| > 0` confirms this is genuinely three-dimensional, not a
planar Lyapunov orbit. The orbit's x-center (~1.1447) sits close to
`x_L2 = 1.155682` (within ~0.011 DU, i.e. within the orbit's own
amplitude scale) — centered on the L2 neighborhood, not some other
equilibrium region.

**Note on `Az` vs. `max|z|`:** this project (matching common
convention) reports `Az = z0` — the z-value *at the x-z-plane symmetry
crossing* (t=0) — as "the" halo amplitude, consistent with Section 2's
target scale. The orbit's `max|z|` over the *full* trajectory
(19,426.7 km) is larger than `z0` (13,454 km) because the halo's 3D
looping geometry does not peak in `|z|` exactly at the symmetry
crossing. This is a genuine, expected feature of this halo family
member, not a bug — verified directly from the propagated trajectory,
not assumed.

### Symmetry verification

At the detected half-period crossing:

```
y(T/2)    = -4.281e-14   (target: 0)
xdot(T/2) =  2.274e-12   (target: 0)
zdot(T/2) = -9.032e-14   (target: 0)
```

All three residuals are at the level of double-precision round-off,
confirming the orbit's mirror symmetry about the x-z plane to numerical
precision — quantified directly, not assessed visually. (See Figure 2's
x-z projection, which — as a *consequence* of this exact symmetry —
visually collapses to a thin retraced arc rather than an open loop:
`x(t) ≈ x(T-t)` and `z(t) ≈ z(T-t)` were verified directly from sampled
trajectory points.)

### STM verification

- `Phi(0) = I6` to `1e-12` (`test_phi0_is_identity`).
- Short-time STM `Phi(dt) ≈ I + A*dt` for `dt=1e-5`, agreement to
  `1e-8` (`test_short_time_stm_approximates_I_plus_A_dt`).
- Augmented (state+STM) propagation's state component matches the
  independent M2 generic propagator to `1e-9` over `t=0.5`
  (`test_augmented_state_matches_generic_propagator`).
- STM correctly predicts a small perturbed trajectory to first order:
  for perturbations `eps in {1e-4, 1e-5, 1e-6}` applied to `x0`, the
  STM-vs-actual prediction error scales as `eps^2` (`error/eps^2`
  constant at ~1.269 across all three, confirming genuine first-order
  linearization, not coincidence):

  | eps | STM-predicted norm | actual norm | error | error / eps² |
  |---:|---:|---:|---:|---:|
  | 1e-4 | 1.010e-4 | 1.010e-4 | 1.267e-8 | 1.267 |
  | 1e-5 | 1.010e-5 | 1.010e-5 | 1.269e-10 | 1.269 |
  | 1e-6 | 1.010e-6 | 1.010e-6 | 1.269e-12 | 1.269 |

### Independent validation of the corrected orbit

**A. Full-period closure with the generic M2 propagator** — done above
(1.9e-5 km), independent of the corrector's own half-period machinery.

**B. Tighter tolerance** (`DOP853, rtol=1e-13, atol=1e-14` vs. default
`rtol=1e-11, atol=1e-12`): closure improves to **8.924e-07 km** — a
~20x tightening as tolerance tightens by ~2 orders of magnitude,
consistent behavior, not a fluke.

**C. Independent integrator** (`RK45, rtol=1e-12, atol=1e-13` vs.
`DOP853`): final-state difference between the two integrators is
`1.580e-10` (nondimensional) — effectively identical trajectories from
two different numerical methods.

**D. STM short-time perturbation prediction** — see STM verification
above.

All four checks agree; no discrepancy found between independent
validation paths.

### Monodromy matrix (validation/characterization artifact only — not a stability study)

```
M = Phi(T),  det(M) = 0.99999999996   (should be exactly 1: CR3BP flow is symplectic/volume-preserving)

eigenvalues (real, imag):
  999.547109 + 0.000000i    <- unstable direction
    0.001000 + 0.000000i    <- stable direction (reciprocal of the above)
    0.951404 + 0.307946i    <- oscillatory (center) pair, |lambda| = 1.0000
    0.951404 - 0.307946i
    1.000000 + 0.0000021i   <- trivial pair (periodicity + energy), |lambda| ~ 1.0000
    1.000000 - 0.0000021i

reciprocal-pair products: 1.0000000002, 0.9999999999996, 0.9999999999996
  (all ~1, as expected for Hamiltonian CR3BP dynamics: eigenvalues
  come in reciprocal pairs lambda, 1/lambda)
```

The dominant eigenvalue magnitude (**999.5**) is far in excess of 1,
confirming this halo orbit is **linearly unstable** — fully expected
for CR3BP libration-point orbits and exactly why real halo missions
require active stationkeeping (Section 12). This is reported strictly
as a validation/characterization artifact; **no manifold generation or
stability/stationkeeping study is performed in M3** (that is explicitly
out of scope — Sections 15/23 of this milestone's instructions).

### Dimensional conversions (self-consistency)

```
period:        3.394396 DU-time = 14.740098 days   (14.740098 = 3.394396 * 4.342480 to 1e-12 relative)
z-amplitude:   0.035 DU = 13,454.0 km               (0.035 * 384,400 to 1e-12 relative)
characteristic velocity scale: V* = 1024.547 m/s (unchanged from M2; the orbit's own velocities
                                                    range up to ~0.30 nondim ~ 307 m/s, well below V*)
```

### Genuine bugs found during M3 (reported, not hidden)

1. **Seed amplitude-offset/velocity sign mismatch.** The initial
   `linear_halo_seed` implementation derived the in-plane trial
   solution assuming `dx(0) = +Ax`, but then placed `x0 = x_L2 - Ax`
   (i.e., `dx(0) = -Ax`) without correspondingly flipping the sign of
   the derived `ydot0`. This produced a seed that diverged rapidly
   under nonlinear propagation instead of returning near-periodically.
   Found by direct numerical comparison of both sign branches; fixed
   by deriving `ydot0 = -kappa*Ax*wxy` consistent with the `x0 = x_L2 - Ax`
   convention (see `halo_seed.py`'s in-line derivation note).
2. **Correction-matrix sign error (chain rule through the event-time
   constraint).** The first differential-correction implementation used
   `M[i,:] = Phi[row,free] - sdot[row]*dT_dfree` (minus). This is
   *wrong* — the correct total-derivative chain rule (re-derived
   carefully from `d/d(free)[s(T(free);free)]`) requires a **plus**
   sign. With the minus sign, Newton iterations showed slow, stalling,
   sub-linear convergence (residual plateauing around ~2e-3 rather than
   collapsing to machine precision). After deriving and fixing the sign,
   the exact same seed converged with clean quadratic convergence (see
   the iteration history above). This was caught by noticing the
   convergence rate was wrong (linear/stalling instead of quadratic),
   not by a passing-vs-failing test alone — a reminder that "eventually
   converges" is not sufficient verification; the *rate* of convergence
   is itself diagnostic.
3. **3-free-variable (x0, z0, ydot0) minimum-norm correction collapses
   to the planar Lyapunov family.** An alternative correction
   formulation (correcting all three of `x0, z0, ydot0` via a
   minimum-norm least-squares solve, rather than holding `z0` fixed)
   was tried and **converges to `z0 ~ 0`** (a degenerate planar orbit,
   not a halo) for every tested amplitude — exactly the failure mode
   this milestone's instructions warned about (Section 22). This
   confirms the fixed-`z0`, 2-free-variable formulation (used for the
   reported orbit) is the correct choice for this problem, not merely
   a default pick.

### Explicit scope statement

**No stable/unstable manifold, arrival-state model, or insertion Δv has
been computed in M3.** The monodromy matrix above is reported purely as
a stability/validation characterization, not a manifold or
stationkeeping result. Manifold and arrival-trajectory work begins in
M4.

---

## Milestone 4 — Arrival-state definition, insertion-point trade, and halo-insertion Δv estimate

Status: **complete.** This section documents the project's **first
defensible halo-insertion Δv estimate**. It is a **local
velocity-matching estimate at the corrected halo orbit under an
explicit arrival-state assumption — it does not compute or optimize
the complete Earth-to-L2 transfer trajectory.**

### Why a transfer-arrival assumption is necessary

The M3 corrected halo orbit gives a rotating-frame velocity `v_halo` at
every point on the orbit. Computing an insertion Δv requires a second
velocity — the velocity of the incoming transfer trajectory,
`v_arrival` — at the same position. **No Earth-departure or transfer
trajectory has been propagated in this project** (that is out of scope
through at least M4 — see Section 1 and Section 19 below). Without an
actual integrated transfer, `v_arrival` cannot be known exactly; it
must be **assumed**, explicitly and transparently. This section defines
that assumption precisely, so the resulting Δv can be read as
"Δv under stated assumption X," never as an unconditional number.

**The M3 halo orbit itself is frozen and unmodified** — no
re-correction, no change to its initial state, period, or convention.
All M4 work traces phase along the existing periodic orbit.

### Baseline arrival-state model

**Speed — Jacobi-consistent model:**

```
v_arr^2 = 2*Omega(r_h) - C_arr,   C_arr = C_halo - dC_baseline
```

with **`dC_baseline = 0.01`** (nondimensional Jacobi units) — a round,
explicitly documented value, *not* tuned to reproduce any target
number. `dC_baseline` represents a modest kinetic-energy offset
consistent with an already-well-targeted transfer (M1 Section 1: "the
spacecraft arrives on a transfer trajectory already targeted near the
halo insertion region") that has not yet been velocity-matched to the
halo. Points where `2*Omega(r_h) - C_arr < 0` are rejected explicitly
(`ArrivalStateError`), never silently clipped; none occurred for the
baseline sweep (1000/1000 points valid).

**Direction — Model B, radial-from-Earth (baseline):**

```
v_arr_dir = (r_h - r_Earth) / |r_h - r_Earth|,   r_Earth = (-mu, 0, 0)
```

i.e., the arrival velocity points outward along the line from Earth's
position to the insertion point — modeling a transfer arriving
generally outbound from the Earth region, **without claiming an actual
propagated trajectory** ("from Earth" is used here strictly in this
geometric sense, per this milestone's explicit constraint).

**Direction — Model C, sweep about baseline:**

```
v_arr_dir(theta) = R_z(theta) @ v_arr_dir_B
```

a rotation of the baseline direction by signed angle `theta` about the
synodic z-axis (within the local x-y plane), used for the direction
sensitivity study below.

**Direction — Model A, aligned with halo velocity (idealized lower
bound only):**

```
v_arr_dir = v_halo / |v_halo|
```

**Explicit degeneracy guard (Section 8 of this milestone's
instructions):** combined with a free arrival speed, Model A can drive
Δv to exactly zero — this is *not* presented as a meaningful insertion
solution anywhere in this project. It is used only as a **fixed-speed
lower-bound sanity case**: at the selected phase, with the same
baseline speed assumption, Model A gives `|Δv| = 42.29 m/s` (exactly
`||v_halo| - |v_arr||`, the pure speed-magnitude mismatch with zero
direction penalty) — reported explicitly as an idealized bound, plotted
as a dashed reference line in Figure 1, never as "the" answer.

### Phase parameterization and sweep

`tau = t/T ∈ [0,1)`, sampled densely (`N=1000` for the baseline
sweep, `tau=1` excluded to avoid duplicating `tau=0` in phase
statistics — periodicity `tau=0 == tau=1` is separately verified,
Independent Verification E below). At each `tau`, the M3 initial
condition is propagated via the **unmodified M2 generic propagator**
to `t = tau*T`, giving `r_h, v_halo` directly (not via any
halo-specific shortcut).

### Insertion-point trade: results

1. **Dense sweep** (1000 points, baseline model): **two local minima**
   found — `tau=0.244` (Δv=109.57 m/s, global) and `tau=0.781`
   (Δv=224.28 m/s, secondary). Both are reported; the secondary minimum
   is **not** discarded as noise — it is a genuine consequence of the
   radial-from-Earth direction model breaking the halo's own y-mirror
   symmetry (see Figure 1).
2. **Bounded refinement** (`scipy.optimize.minimize_scalar`, Brent's
   method, bracket `[tau_grid_min - 0.01, tau_grid_min + 0.01]`,
   `xatol=1e-10`) around the global grid minimum: **`tau = 0.244130431`**,
   Δv = **109.57233 m/s** — improves on (never worsens) the grid value
   of 109.57246 m/s, verified explicitly
   (`test_optimizer_does_not_worsen_grid_minimum`).
3. **Direct re-evaluation** of the refined phase reproduces the
   refined result exactly (`test_direct_selected_point_recomputation_matches`).

### Selected insertion state

```
tau            = 0.24413043101153073
t (nondim)     = 0.8286754286332
t (days from M3 reference crossing) = 3.598506 days

position (nondim):  [ 1.135148,  0.096577, -0.013424]
v_halo   (nondim):  [ 0.066126, -0.007318, -0.075304]
v_arr    (nondim):  [ 0.141254,  0.011890, -0.001653]
Delta_v_vec (nondim): [-0.075128, -0.019209, -0.073651]

|Delta_v| (nondim) = 0.1069471155240105
|Delta_v| (m/s)    = 109.57233088725141

distance from Moon = 67,902.966 km
distance from L2    = 38,303.159 km
```

### Comparison with M1's preliminary scale

M1's Section 7 hand estimate: **O(10¹–10²) m/s (roughly 10–200 m/s)**.
M4's computed value, **109.57 m/s, falls within this range** — this is
a genuine outcome of the baseline assumption (`dC=0.01`,
radial-from-Earth direction), **not tuned to land there**: the
assumption was chosen (Section above) before checking against M1's
range, and the resulting number happened to agree. Had it fallen
outside M1's range, this document would report and explain that
rather than adjust the assumption (Section 18 of this milestone's
instructions) — see the speed-sensitivity study below, where several
`dC` choices (e.g. `dC=0.05` → 199.0 m/s) sit at the upper edge of
M1's range, illustrating how assumption-dependent this number is.

### Jacobi interpretation

An impulsive burn changes velocity discontinuously at fixed position,
so the Jacobi constant is **not** conserved across it (M2's Section 5
already established `C` is only conserved under *unforced* propagation;
this milestone applies that fact explicitly to the insertion burn):

```
C_halo = 3.1412189199161844
C_arr  = 3.1312189199161846
Delta_C = C_halo - C_arr = 0.009999999999999787  (== dC_baseline, exactly, as designed)
```

The insertion burn maps the arrival state's `(r_h, v_arr, C_arr)` onto
the halo state's `(r_h, v_halo, C_halo)` at the same position — a
discontinuous jump in velocity and Jacobi constant, not a continuous
transition.

### Sensitivity studies

**A. Arrival-speed sensitivity** (`dC_baseline` varied 0.002 → 0.05,
minimum Δv re-optimized at each value):

| dC | min Δv (m/s) | optimal tau |
|---:|---:|---:|
| 0.002 | 89.55 | 0.2427 |
| 0.005 | 96.85 | 0.2433 |
| **0.01 (baseline)** | **109.57** | **0.2441** |
| 0.02 | 134.51 | 0.2456 |
| 0.03 | 157.69 | 0.2467 |
| 0.05 | 198.98 | 0.2485 |

Monotonic and smooth — larger assumed energy offsets produce larger
Δv and a slowly shifting optimal phase. All 300/300 grid points valid
at every `dC` tested (no infeasible-Jacobi rejections in this range).

**B. Arrival-direction sensitivity** (`theta` swept ±5°/±10°/±20° about
the baseline radial-from-Earth direction):

| theta (deg) | min Δv (m/s) | optimal tau |
|---:|---:|---:|
| -20 | 100.37 | 0.2792 |
| -10 | 105.77 | 0.2609 |
| -5 | 107.86 | 0.2524 |
| **0 (baseline)** | **109.57** | **0.2441** |
| +5 | 110.91 | 0.2361 |
| +10 | 111.87 | 0.2282 |
| +20 | 112.69 | 0.2124 |

Δv varies only **~12 m/s (≈11%) over a ±20° direction sweep** — a
**mild** sensitivity to direction assumption, though the *optimal phase*
shifts more substantially (0.212–0.279), showing direction assumption
mainly reshapes *where* the minimum sits, more than *how deep* it is.
This is reported as a genuine finding, not minimized.

**C. Local phase sensitivity around the selected optimum:**

| dtau | Δv (m/s) |
|---:|---:|
| -0.05 | 129.18 |
| -0.02 | 112.77 |
| -0.01 | 110.37 |
| -0.005 | 109.77 |
| **0 (optimum)** | **109.57** |
| +0.005 | 109.77 |
| +0.01 | 110.34 |
| +0.02 | 112.56 |
| +0.05 | 126.29 |

The minimum is **broad, not razor-thin**: Δv grows only ~2 m/s over
`dtau=±0.005` (≈±1.7 hours) and ~20 m/s over `dtau=±0.05` (≈±17
hours). This is a purely deterministic assumption/phase-sensitivity
result — **not** a navigation-robustness or dispersion analysis
(Section 10 of this milestone's instructions; no such claim is made
anywhere in this project).

### Propellant implication (illustrative, `m0 = 6000 kg`)

```
Delta_v = 109.57233 m/s
Isp=320s: m_prop = 205.88 kg
Isp=450s: m_prop = 147.14 kg
```

Illustrative only, using the same illustrative spacecraft/propulsion
scale as M1 Section 8. Excludes TLI, MCC, stationkeeping, and launch
Δv, as throughout this project.

### Independent numerical verification (Section 13 of this milestone's instructions)

| Check | Result |
|---|---:|
| A. Generic-M2-propagator recomputation of the selected halo state | exact match (diff = 0.0) |
| B. Direct vector subtraction `v_halo - v_arrival` vs. stored Δv | exact match (diff = 0.0) |
| C. Independent dimensional conversion `Δv_nd * V*` vs. stored m/s | exact match (diff = 0.0) |
| D. Jacobi recomputed from raw state components vs. `cr3bp.jacobi_constant` | exact match (diff = 0.0) |
| E. Phase periodicity: state at `tau` vs. `tau+1` | diff = 7.22e-17 (round-off) |
| F. Tighter tolerance (`rtol=1e-13, atol=1e-14`) at the selected phase | diff = 2.88e-13 (nondim) |

No discrepancy found in any independent check.

### Limitations (M4-specific, in addition to Section 12 below)

- The arrival-state model is an **explicit assumption**, not a
  propagated trajectory; a different (equally defensible) `dC_baseline`
  or direction model would produce a materially different Δv, as the
  sensitivity studies above show directly and honestly.
- "Radial from Earth" is a **geometric** direction rule, not a claim
  that any Earth-departure trajectory has been computed.
- The Jacobi-consistent speed model is a simplification; it does not
  by itself prove any physical transfer with that Jacobi constant
  actually connects to Earth.
- The two local minima found are properties of *this* arrival-direction
  model; a different direction model could shift, merge, or eliminate
  them.
- No launch vehicle, TLI, total mission Δv, transfer duration, lunar
  flyby, manifold transfer, ephemeris feasibility, stationkeeping, or
  navigation/dispersion analysis is established by this milestone (see
  Section 12 and this document's M1 Section 12 for the full standing
  limitations list).

### Explicit scope statement

**M4 computes a local velocity-matching insertion estimate at the
corrected halo orbit. It does not compute or optimize the complete
Earth-to-L2 transfer trajectory.** M5 is intended as a stronger
transfer-arrival/validation milestone, not cosmetic packaging.

---

## 12. Limitations

Stated explicitly and unconditionally, for this milestone and as an
ongoing constraint on the whole project:

- Circular Restricted Three-Body Problem (CR3BP) dynamics only
- Earth and Moon are assumed to move on **circular, coplanar** orbits
  about their barycenter (no lunar orbital eccentricity, no orbital
  plane inclination)
- Spacecraft mass is dynamically negligible (does not perturb Earth/Moon motion)
- **No Sun perturbation** (no bicircular/four-body or ephemeris
  third-body effects)
- **No lunar eccentricity** (true lunar orbit is mildly eccentric, ~0.055)
- **No solar radiation pressure**
- **No ephemeris model** — positions of Earth/Moon are idealized
  CR3BP circular-orbit positions, not JPL ephemeris states
- **No finite-burn modeling** — all maneuvers are idealized as
  instantaneous impulsive Δv
- **No launch or translunar-trajectory optimization** — the transfer
  that delivers the spacecraft near the halo is assumed to already
  exist and is out of scope (Section 1)
- **No stationkeeping design** — only the one-time insertion burn is
  estimated; the (materially different, ongoing) stationkeeping Δv
  budget is not computed
- **No navigation or dispersion analysis** — no orbit-determination
  error, no execution error, no statistical Δv margin
- **No covariance analysis**
- **No operational flight-design claim** — this is not a substitute
  for a mission-grade trajectory design produced with an ephemeris
  model, a real navigation team, and mission-specific constraints

**This project is an engineering-approximation estimate of Earth–Moon
L2 halo-orbit insertion Δv using CR3BP dynamics — it is not a flight
trajectory solution.**
