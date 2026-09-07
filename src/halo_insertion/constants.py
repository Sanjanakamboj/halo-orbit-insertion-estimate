"""Earth-Moon CR3BP physical constants — single source of truth.

These reproduce the M1 hand-calculated values documented in
DESIGN.md Section 3 exactly. Any code elsewhere in this package must
derive its numbers from here, not re-hardcode them.
"""

import math

# --- Primary gravitational parameters (km^3/s^2) ---
GM_EARTH_KM3_S2 = 398600.4418
GM_MOON_KM3_S2 = 4902.8000

# --- Mean Earth-Moon center-to-center distance (km) ---
DU_KM = 384400.0

# --- CR3BP mass ratio: mu = m_moon / (m_earth + m_moon) ---
MU = GM_MOON_KM3_S2 / (GM_EARTH_KM3_S2 + GM_MOON_KM3_S2)

# --- Mean motion of the Earth-Moon rotating frame (rad/s) ---
N_RAD_S = math.sqrt((GM_EARTH_KM3_S2 + GM_MOON_KM3_S2) / DU_KM**3)

# --- Time unit: 1/n, in seconds and days ---
TU_S = 1.0 / N_RAD_S
TU_DAYS = TU_S / 86400.0

# --- Characteristic velocity V* = n * DU (km/s) ---
VSTAR_KM_S = N_RAD_S * DU_KM
VSTAR_M_S = VSTAR_KM_S * 1000.0

# --- Primary positions in normalized synodic coordinates ---
# Earth at (-mu, 0, 0), Moon at (1-mu, 0, 0)
EARTH_POS_NONDIM = (-MU, 0.0, 0.0)
MOON_POS_NONDIM = (1.0 - MU, 0.0, 0.0)

# --- Standard gravity, for rocket-equation use elsewhere ---
G0_M_S2 = 9.80665
