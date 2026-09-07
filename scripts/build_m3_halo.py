"""M3 driver: build the corrected Earth-Moon L2 halo orbit, run the
full independent validation suite, and save results/ artifacts.

Run from the repository root:
    python scripts/build_m3_halo.py
"""

import csv
import json

import numpy as np

from halo_insertion import constants as c
from halo_insertion.cr3bp import jacobi_constant
from halo_insertion.equilibria import find_L2
from halo_insertion.halo import build_l2_halo, full_period_closure, monodromy_matrix
from halo_insertion.propagation import PropagationSettings, propagate
from halo_insertion.variational import propagate_with_stm

MU = c.MU
AZ_NONDIM = 0.035  # nondimensional out-of-plane amplitude target (M1 scale: 0.03-0.04 DU)
Z_SIGN = 1.0  # this project's southern-family sign convention (see DESIGN.md M3 section)


def main():
    l2 = find_L2(MU)

    # --- Build and correct ---
    result = build_l2_halo(l2.x, MU, AZ_NONDIM, z_sign=Z_SIGN)
    state0 = result.state0
    period = result.period

    # --- Full-period closure, GENERIC M2 propagator (independent of corrector) ---
    final_state, closure, _ = full_period_closure(state0, period, MU)
    pos_closure_nondim = float(np.linalg.norm(closure[:3]))
    vel_closure_nondim = float(np.linalg.norm(closure[3:]))
    pos_closure_km = pos_closure_nondim * c.DU_KM
    vel_closure_m_s = vel_closure_nondim * c.VSTAR_M_S

    # --- Independent validation B: tighter tolerance ---
    tight_settings = PropagationSettings(method="DOP853", rtol=1e-13, atol=1e-14)
    _, closure_tight, _ = full_period_closure(state0, period, MU, settings=tight_settings)
    pos_closure_tight_km = float(np.linalg.norm(closure_tight[:3])) * c.DU_KM

    # --- Independent validation C: different integrator (RK45) ---
    rk45_settings = PropagationSettings(method="RK45", rtol=1e-12, atol=1e-13)
    final_rk45, closure_rk45, _ = full_period_closure(state0, period, MU, settings=rk45_settings)
    dop853_vs_rk45_diff = float(np.linalg.norm(final_state - final_rk45))

    # --- Independent validation D: STM short-time perturbation prediction ---
    dt_pred = 0.01
    stm_short = propagate_with_stm(state0, (0.0, dt_pred), MU)
    Phi_dt = stm_short.Phi_final
    eps_values = [1e-4, 1e-5, 1e-6]
    stm_prediction_errors = []
    unperturbed = propagate(state0, (0.0, dt_pred), MU).y[:, -1]
    for eps in eps_values:
        perturbation = np.array([eps, 0, 0, 0, 0, 0])
        predicted_delta = Phi_dt @ perturbation
        perturbed_final = propagate(state0 + perturbation, (0.0, dt_pred), MU).y[:, -1]
        actual_delta = perturbed_final - unperturbed
        err = float(np.linalg.norm(actual_delta - predicted_delta))
        stm_prediction_errors.append({"eps": eps, "error": err, "error_over_eps_sq": err / eps**2})

    # --- Jacobi conservation over full orbit ---
    n_eval = 4000
    t_eval = np.linspace(0.0, period, n_eval)
    res_eval = propagate(state0, (0.0, period), MU, t_eval=t_eval)
    C0 = jacobi_constant(state0, MU)
    C_t = np.array([jacobi_constant(res_eval.y[:, i], MU) for i in range(res_eval.y.shape[1])])
    jacobi_max_drift = float(np.max(np.abs(C_t - C0)))
    jacobi_rms_drift = float(np.sqrt(np.mean((C_t - C0) ** 2)))

    # --- Symmetry residual at half period ---
    half_state = propagate(state0, (0.0, result.half_period), MU).y[:, -1]
    symmetry_residual = {
        "y_at_half_period": float(half_state[1]),
        "xdot_at_half_period": float(half_state[3]),
        "zdot_at_half_period": float(half_state[5]),
    }

    # --- Geometry ---
    x, y, z = res_eval.y[0], res_eval.y[1], res_eval.y[2]
    x_min, x_max = float(x.min()), float(x.max())
    y_min, y_max = float(y.min()), float(y.max())
    z_min, z_max = float(z.min()), float(z.max())
    max_abs_z = float(np.max(np.abs(z)))

    # --- Monodromy matrix ---
    M = monodromy_matrix(state0, period, MU)
    eigvals = np.linalg.eigvals(M)
    eig_mags = np.abs(eigvals)
    det_M = complex(np.linalg.det(M))
    # pair up reciprocal magnitudes (sorted descending) for reporting
    sorted_idx = np.argsort(-eig_mags)
    sorted_mags = eig_mags[sorted_idx]
    reciprocal_products = [
        float(sorted_mags[i] * sorted_mags[-(i + 1)]) for i in range(len(sorted_mags) // 2)
    ]

    # === Save results/m3_halo_initial_state.csv ===
    with open("results/m3_halo_initial_state.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["field", "value"])
        writer.writerow(["mu", MU])
        writer.writerow(["x_L2", l2.x])
        writer.writerow(["x0", state0[0]])
        writer.writerow(["y0", state0[1]])
        writer.writerow(["z0", state0[2]])
        writer.writerow(["xdot0", state0[3]])
        writer.writerow(["ydot0", state0[4]])
        writer.writerow(["zdot0", state0[5]])
        writer.writerow(["half_period_nondim", result.half_period])
        writer.writerow(["period_nondim", period])
        writer.writerow(["period_days", result.period_days])

    # === Save results/m3_correction_history.csv ===
    with open("results/m3_correction_history.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            ["iteration", "x0", "z0", "ydot0", "half_period", "xdot_residual", "zdot_residual", "correction_norm"]
        )
        for it in result.correction.iterations:
            writer.writerow(
                [it.iteration, it.x0, it.z0, it.ydot0, it.half_period, it.xdot_residual, it.zdot_residual, it.correction_norm]
            )

    # === Save results/m3_halo_summary.json ===
    summary = {
        "mu": MU,
        "x_L2": l2.x,
        "seed_method": "linearized (first-order) CR3BP variational solution about L2, NOT third-order Richardson",
        "Az_target_nondim": AZ_NONDIM,
        "z_sign_convention": Z_SIGN,
        "corrected_initial_state": {
            "x0": float(state0[0]),
            "y0": float(state0[1]),
            "z0": float(state0[2]),
            "xdot0": float(state0[3]),
            "ydot0": float(state0[4]),
            "zdot0": float(state0[5]),
        },
        "half_period_nondim": float(result.half_period),
        "period_nondim": float(period),
        "period_days": float(result.period_days),
        "correction_iterations": len(result.correction.iterations),
        "correction_converged": bool(result.correction.converged),
        "jacobi_C0": float(C0),
        "jacobi_max_drift": jacobi_max_drift,
        "jacobi_rms_drift": jacobi_rms_drift,
        "geometry": {
            "x_min_nondim": x_min,
            "x_max_nondim": x_max,
            "y_min_nondim": y_min,
            "y_max_nondim": y_max,
            "z_min_nondim": z_min,
            "z_max_nondim": z_max,
            "max_abs_z_nondim": max_abs_z,
            "max_abs_z_km": max_abs_z * c.DU_KM,
            "z0_amplitude_km": float(state0[2]) * c.DU_KM,
            "x_extent_from_L2_km": [
                (x_min - l2.x) * c.DU_KM,
                (x_max - l2.x) * c.DU_KM,
            ],
            "y_extent_km": [y_min * c.DU_KM, y_max * c.DU_KM],
        },
        "full_period_closure_generic_propagator": {
            "position_norm_nondim": pos_closure_nondim,
            "position_norm_km": pos_closure_km,
            "velocity_norm_nondim": vel_closure_nondim,
            "velocity_norm_m_s": vel_closure_m_s,
            "component_wise": {
                "dx": float(closure[0]),
                "dy": float(closure[1]),
                "dz": float(closure[2]),
                "dxdot": float(closure[3]),
                "dydot": float(closure[4]),
                "dzdot": float(closure[5]),
            },
        },
        "independent_validation": {
            "tighter_tolerance_closure_km": pos_closure_tight_km,
            "tighter_tolerance_settings": {"method": "DOP853", "rtol": 1e-13, "atol": 1e-14},
            "rk45_vs_dop853_final_state_diff_nondim": dop853_vs_rk45_diff,
            "stm_short_time_prediction": stm_prediction_errors,
        },
        "symmetry_residual": symmetry_residual,
        "monodromy": {
            "determinant_real": det_M.real,
            "determinant_imag": det_M.imag,
            "eigenvalues_real": [float(e.real) for e in eigvals],
            "eigenvalues_imag": [float(e.imag) for e in eigvals],
            "eigenvalue_magnitudes": [float(m) for m in eig_mags],
            "reciprocal_pair_products": reciprocal_products,
            "max_eigenvalue_magnitude": float(np.max(eig_mags)),
            "linearly_unstable": bool(np.max(eig_mags) > 1.0 + 1e-6),
        },
    }
    with open("results/m3_halo_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print("Wrote results/m3_halo_initial_state.csv")
    print("Wrote results/m3_correction_history.csv")
    print("Wrote results/m3_halo_summary.json")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
