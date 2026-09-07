"""M4 driver: sweep insertion phase over the frozen M3 halo orbit,
refine the best candidate, run sensitivity studies and independent
verification, and save results/ artifacts.

Run from the repository root (after scripts/build_m3_halo.py):
    python scripts/build_m4_insertion.py
"""

import csv
import json
import math

import numpy as np

from halo_insertion import constants as c
from halo_insertion.arrival import jacobi_from_state
from halo_insertion.cr3bp import jacobi_constant
from halo_insertion.insertion import evaluate_insertion_at_phase, halo_state_at_phase, phase_sweep, refine_best_phase
from halo_insertion.propagation import PropagationSettings, propagate

MU = c.MU
DC_BASELINE = 0.01  # nondimensional Jacobi offset defining the baseline arrival speed
N_SWEEP = 1000
M0_KG = 6000.0
ISP_VALUES = [320.0, 450.0]
G0 = 9.80665


def _load_m3():
    with open("results/m3_halo_summary.json") as f:
        summary = json.load(f)
    s0 = summary["corrected_initial_state"]
    state0 = [s0["x0"], s0["y0"], s0["z0"], s0["xdot0"], s0["ydot0"], s0["zdot0"]]
    period = summary["period_nondim"]
    return state0, period, summary


def _find_local_minima(candidates):
    n = len(candidates)
    minima = []
    for i in range(n):
        if not candidates[i].valid:
            continue
        prev_c = candidates[(i - 1) % n]
        next_c = candidates[(i + 1) % n]
        if not prev_c.valid or not next_c.valid:
            continue
        if candidates[i].delta_v_m_s < prev_c.delta_v_m_s and candidates[i].delta_v_m_s < next_c.delta_v_m_s:
            minima.append((candidates[i].tau, candidates[i].delta_v_m_s))
    return sorted(minima, key=lambda x: x[1])


def _propellant(delta_v_m_s, m0, isp, g0):
    return m0 * (1.0 - math.exp(-delta_v_m_s / (isp * g0)))


def main():
    state0, period, m3_summary = _load_m3()

    # === Baseline (theta=0) dense phase sweep ===
    candidates = phase_sweep(state0, period, MU, DC_BASELINE, n_points=N_SWEEP)
    valid = [c_ for c_ in candidates if c_.valid]
    dvs = np.array([c_.delta_v_m_s for c_ in valid])
    taus = np.array([c_.tau for c_ in valid])

    local_minima = _find_local_minima(candidates)
    grid_best_idx = int(np.argmin(dvs))
    grid_best_tau = float(taus[grid_best_idx])
    grid_best_dv = float(dvs[grid_best_idx])

    # === Refine the global minimum ===
    bracket = (grid_best_tau - 0.01, grid_best_tau + 0.01)
    refined = refine_best_phase(state0, period, MU, DC_BASELINE, bracket)

    assert refined.delta_v_m_s <= grid_best_dv + 1e-6, "refinement must not worsen the grid minimum"

    # === Save results/m4_phase_sweep.csv ===
    with open("results/m4_phase_sweep.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "tau", "t_nondim", "t_days", "x", "y", "z", "vx_halo", "vy_halo", "vz_halo",
                "vx_arr", "vy_arr", "vz_arr", "delta_v_nondim", "delta_v_m_s",
                "C_halo", "C_arr", "delta_C", "dist_from_moon_km", "dist_from_l2_km", "valid",
            ]
        )
        for cand in candidates:
            writer.writerow(
                [
                    cand.tau, cand.t_nondim, cand.t_days,
                    cand.r_h[0], cand.r_h[1], cand.r_h[2],
                    cand.v_halo[0], cand.v_halo[1], cand.v_halo[2],
                    cand.v_arr[0], cand.v_arr[1], cand.v_arr[2],
                    cand.delta_v_nondim, cand.delta_v_m_s,
                    cand.C_halo, cand.C_arr, cand.delta_C,
                    cand.dist_from_moon_km, cand.dist_from_l2_km, cand.valid,
                ]
            )

    # === Sensitivity A: arrival-speed (dC) sensitivity ===
    speed_sensitivity = []
    for dC in [0.002, 0.005, 0.01, 0.02, 0.03, 0.05]:
        cands = phase_sweep(state0, period, MU, dC, n_points=300)
        v = [c_ for c_ in cands if c_.valid]
        dv_arr = np.array([c_.delta_v_m_s for c_ in v])
        tau_arr = np.array([c_.tau for c_ in v])
        i = int(np.argmin(dv_arr))
        r = refine_best_phase(state0, period, MU, dC, (tau_arr[i] - 0.01, tau_arr[i] + 0.01))
        speed_sensitivity.append({"dC": dC, "min_delta_v_m_s": r.delta_v_m_s, "tau": r.tau, "n_valid": len(v), "n_total": len(cands)})

    # === Sensitivity B: arrival-direction sensitivity ===
    direction_sensitivity = []
    for theta_deg in [-20, -10, -5, 0, 5, 10, 20]:
        theta_rad = math.radians(theta_deg)
        cands = phase_sweep(state0, period, MU, DC_BASELINE, n_points=300, direction_theta_rad=theta_rad)
        v = [c_ for c_ in cands if c_.valid]
        dv_arr = np.array([c_.delta_v_m_s for c_ in v])
        tau_arr = np.array([c_.tau for c_ in v])
        i = int(np.argmin(dv_arr))
        r = refine_best_phase(state0, period, MU, DC_BASELINE, (tau_arr[i] - 0.01, tau_arr[i] + 0.01), direction_theta_rad=theta_rad)
        direction_sensitivity.append({"theta_deg": theta_deg, "min_delta_v_m_s": r.delta_v_m_s, "tau": r.tau})

    # === Sensitivity C: local phase sensitivity around the optimum ===
    phase_sensitivity = []
    for dtau in [-0.05, -0.02, -0.01, -0.005, 0.0, 0.005, 0.01, 0.02, 0.05]:
        cand = evaluate_insertion_at_phase(state0, period, MU, refined.tau + dtau, DC_BASELINE)
        phase_sensitivity.append({"dtau": dtau, "tau": cand.tau, "delta_v_m_s": cand.delta_v_m_s})

    # === Save results/m4_sensitivity.csv ===
    with open("results/m4_sensitivity.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["study", "parameter", "value", "min_delta_v_m_s", "tau"])
        for row in speed_sensitivity:
            writer.writerow(["arrival_speed_dC", "dC", row["dC"], row["min_delta_v_m_s"], row["tau"]])
        for row in direction_sensitivity:
            writer.writerow(["arrival_direction_theta_deg", "theta_deg", row["theta_deg"], row["min_delta_v_m_s"], row["tau"]])
        for row in phase_sensitivity:
            writer.writerow(["local_phase_dtau", "dtau", row["dtau"], row["delta_v_m_s"], row["tau"]])

    # === Degeneracy guard: idealized aligned-direction lower bound ===
    aligned_case = evaluate_insertion_at_phase(state0, period, MU, refined.tau, DC_BASELINE, direction_model="aligned")

    # === Independent verification (Section 13) ===
    # A. direct recomputation via generic M2 propagator (already what halo_state_at_phase uses)
    state_recompute = halo_state_at_phase(state0, period, MU, refined.tau)
    recompute_diff = float(np.linalg.norm(state_recompute - np.concatenate([refined.r_h, refined.v_halo])))

    # B. direct vector subtraction cross-check
    dv_direct = refined.v_halo - refined.v_arr
    dv_direct_norm = float(np.linalg.norm(dv_direct))
    b_check_diff = abs(dv_direct_norm - refined.delta_v_nondim)

    # C. independent dimensional conversion
    dv_ms_recompute = refined.delta_v_nondim * c.VSTAR_M_S
    c_check_diff = abs(dv_ms_recompute - refined.delta_v_m_s)

    # D. Jacobi recomputed from raw state components
    C_halo_direct = jacobi_from_state(refined.r_h, refined.v_halo, MU)
    C_halo_via_cr3bp = jacobi_constant(np.concatenate([refined.r_h, refined.v_halo]), MU)
    d_check_diff = abs(C_halo_direct - C_halo_via_cr3bp)

    # E. phase periodicity: tau vs tau+1
    state_tau = halo_state_at_phase(state0, period, MU, refined.tau)
    state_tau_plus1 = halo_state_at_phase(state0, period, MU, refined.tau + 1.0)
    e_check_diff = float(np.linalg.norm(state_tau - state_tau_plus1))

    # F. tighter propagation tolerances at the selected phase
    t_sel = refined.tau * period
    tight_settings = PropagationSettings(method="DOP853", rtol=1e-13, atol=1e-14)
    res_tight = propagate(state0, (0.0, t_sel), MU, settings=tight_settings)
    f_check_diff = float(np.linalg.norm(res_tight.y[:, -1] - state_tau))

    # === Propellant implication ===
    propellant = {
        f"isp_{int(isp)}_s": {
            "m_prop_kg": _propellant(refined.delta_v_m_s, M0_KG, isp, G0),
        }
        for isp in ISP_VALUES
    }

    # === Save results/m4_insertion_summary.json ===
    summary = {
        "m3_halo_reference": {
            "corrected_initial_state": m3_summary["corrected_initial_state"],
            "period_nondim": period,
            "period_days": m3_summary["period_days"],
            "jacobi_C0": m3_summary["jacobi_C0"],
        },
        "arrival_model": {
            "description": "Jacobi-consistent speed (v_arr^2 = 2*Omega(r_h) - C_arr, C_arr = C_halo - dC_baseline); "
            "baseline direction = radial from Earth (-mu,0,0) to r_h, optionally rotated by theta about the synodic z-axis.",
            "dC_baseline": DC_BASELINE,
            "direction_baseline": "radial_from_earth",
        },
        "phase_sweep": {
            "n_points": N_SWEEP,
            "n_valid": len(valid),
            "grid_min_delta_v_m_s": grid_best_dv,
            "grid_min_tau": grid_best_tau,
            "local_minima": [{"tau": t, "delta_v_m_s": dv} for t, dv in local_minima],
        },
        "selected_insertion": {
            "tau": refined.tau,
            "t_nondim": refined.t_nondim,
            "t_days": refined.t_days,
            "position_nondim": refined.r_h.tolist(),
            "v_halo_nondim": refined.v_halo.tolist(),
            "v_arr_nondim": refined.v_arr.tolist(),
            "delta_v_vec_nondim": refined.delta_v_vec.tolist(),
            "delta_v_nondim": refined.delta_v_nondim,
            "delta_v_m_s": refined.delta_v_m_s,
            "C_halo": refined.C_halo,
            "C_arr": refined.C_arr,
            "delta_C": refined.delta_C,
            "dist_from_moon_km": refined.dist_from_moon_km,
            "dist_from_l2_km": refined.dist_from_l2_km,
        },
        "degeneracy_guard_aligned_case": {
            "description": "IDEALIZED lower bound only (direction exactly aligned with v_halo): "
            "not a physically meaningful insertion solution -- see DESIGN.md Section 8.",
            "delta_v_m_s": aligned_case.delta_v_m_s,
            "tau": aligned_case.tau,
        },
        "sensitivity": {
            "arrival_speed_dC": speed_sensitivity,
            "arrival_direction_theta_deg": direction_sensitivity,
            "local_phase_dtau": phase_sensitivity,
        },
        "independent_verification": {
            "A_generic_propagator_recompute_diff": recompute_diff,
            "B_direct_vector_subtraction_diff": b_check_diff,
            "C_dimensional_conversion_diff_m_s": c_check_diff,
            "D_jacobi_raw_vs_cr3bp_diff": d_check_diff,
            "E_phase_periodicity_diff": e_check_diff,
            "F_tighter_tolerance_diff": f_check_diff,
            "F_tighter_tolerance_settings": {"method": "DOP853", "rtol": 1e-13, "atol": 1e-14},
        },
        "propellant_implication_illustrative": {
            "m0_kg": M0_KG,
            **propellant,
        },
        "comparison_with_m1_preliminary_scale": {
            "m1_range_m_s": [10, 200],
            "m4_selected_delta_v_m_s": refined.delta_v_m_s,
            "within_m1_scale": 10.0 <= refined.delta_v_m_s <= 200.0,
        },
    }
    with open("results/m4_insertion_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print("Wrote results/m4_phase_sweep.csv")
    print("Wrote results/m4_sensitivity.csv")
    print("Wrote results/m4_insertion_summary.json")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
