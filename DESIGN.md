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
