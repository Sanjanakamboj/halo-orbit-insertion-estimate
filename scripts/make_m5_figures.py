"""Generate the M5 diagnostic figures.

Figure 1: dynamically propagated arrival arc + halo.
Figure 2: M4 (local) vs M5 (dynamically propagated) insertion estimate.
Figure 3: sensitivity (direction offset) + convergence (round-trip/Jacobi vs tolerance).

Run from the repository root (after scripts/build_m5_validation.py):
    python scripts/make_m5_figures.py
"""

import csv
import json

import matplotlib.pyplot as plt
import numpy as np

from halo_insertion import constants as c
from halo_insertion.equilibria import find_L1, find_L2
from halo_insertion.propagation import propagate

MU = c.MU


def _load():
    with open("results/m5_summary.json") as f:
        m5 = json.load(f)
    with open("results/m4_insertion_summary.json") as f:
        m4 = json.load(f)
    with open("results/m3_halo_summary.json") as f:
        m3 = json.load(f)
    return m3, m4, m5


def make_figure_1(m3, m5):
    l1 = find_L1(MU)
    l2 = find_L2(MU)
    s0 = m3["corrected_initial_state"]
    state0 = np.array([s0["x0"], s0["y0"], s0["z0"], s0["xdot0"], s0["ydot0"], s0["zdot0"]])
    period = m3["period_nondim"]

    t_eval = np.linspace(0, period, 2000)
    res_halo = propagate(state0, (0, period), MU, t_eval=t_eval)

    sel = m5["selected_arc"]
    s_arr = np.array(sel["position_nondim"] + sel["v_arr_prescribed_nondim"])
    t_back_nondim = sel["t_back_days"] / c.TU_DAYS
    t_eval_arc = np.linspace(0, -t_back_nondim, 1000)
    res_arc = propagate(s_arr, (0, -t_back_nondim), MU, t_eval=t_eval_arc)

    fig, ax = plt.subplots(figsize=(8, 7))
    ax.plot(res_halo.y[0], res_halo.y[1], lw=1.0, color="0.7", label="M3 halo orbit")
    ax.plot(res_arc.y[0], res_arc.y[1], lw=1.3, color="tab:purple", label="M5 backward-propagated arrival arc")

    ax.scatter([-MU], [0], color="tab:blue", s=70, zorder=4, label="Earth")
    ax.scatter([1 - MU], [0], color="0.4", s=50, zorder=4, label="Moon")
    ax.scatter([l1.x], [0], color="tab:green", marker="x", s=50, zorder=4, label="L1")
    ax.scatter([l2.x], [0], color="tab:red", marker="x", s=50, zorder=4, label="L2")
    ax.scatter([sel["position_nondim"][0]], [sel["position_nondim"][1]], color="black", s=50, zorder=5, label="insertion point")

    ax.set_xlabel("x (DU)")
    ax.set_ylabel("y (DU)")
    ax.set_aspect("equal", adjustable="datalim")
    ax.legend(fontsize=8, loc="best")
    ax.grid(True, alpha=0.3)

    fig.suptitle(
        "M5 CR3BP arrival-arc validation of L2 halo insertion\n"
        "Backward-propagated ballistic arrival arc; not an optimized Earth-to-halo transfer",
        fontsize=10.5,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.90])
    fig.savefig("figures/m5_fig1_arrival_arc.png", dpi=150)
    plt.close(fig)


def make_figure_2(m4, m5):
    taus, dvs = [], []
    with open("results/m4_phase_sweep.csv") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row["valid"] == "True":
                taus.append(float(row["tau"]))
                dvs.append(float(row["delta_v_m_s"]))

    stage_a_taus, stage_a_dvs, stage_a_accepted = [], [], []
    with open("results/m5_arrival_arc_search.csv") as f:
        reader = csv.DictReader(f)
        for row in reader:
            stage_a_taus.append(float(row["tau"]))
            stage_a_dvs.append(float(row["delta_v_m_s"]))
            stage_a_accepted.append(row["accepted"] == "True")

    stage_a_taus = np.array(stage_a_taus)
    stage_a_dvs = np.array(stage_a_dvs)
    stage_a_accepted = np.array(stage_a_accepted)

    fig, ax = plt.subplots(figsize=(10, 5.5))
    ax.plot(taus, dvs, lw=1.2, color="0.6", label="M4 local arrival-state estimate (Δv vs phase)", zorder=1)

    ax.scatter(
        stage_a_taus[stage_a_accepted], stage_a_dvs[stage_a_accepted],
        color="tab:green", marker="o", s=35, zorder=3,
        label="M5 candidates: dynamically valid (Earthward + Moon-safe)",
    )
    ax.scatter(
        stage_a_taus[~stage_a_accepted], stage_a_dvs[~stage_a_accepted],
        color="tab:red", marker="x", s=30, zorder=2,
        label="M5 candidates: rejected (fails acceptance criteria)",
    )

    m4_sel = m4["selected_insertion"]
    ax.scatter([m4_sel["tau"]], [m4_sel["delta_v_m_s"]], color="black", marker="*", s=180, zorder=5,
               label=f"M4 selected minimum: {m4_sel['delta_v_m_s']:.2f} m/s")

    m5_sel = m5["selected_arc"]
    ax.scatter([m5_sel["tau"]], [m5_sel["delta_v_m_s"]], color="tab:blue", marker="D", s=70, zorder=6,
               label=f"M5 selected (propagated CR3BP arrival arc): {m5_sel['delta_v_m_s']:.2f} m/s")

    ax.set_xlabel("orbital phase  tau = t / T")
    ax.set_ylabel("insertion Delta-v [m/s]")
    ax.legend(fontsize=7.5, loc="upper center")
    ax.grid(True, alpha=0.3)

    fig.suptitle(
        "M4 local arrival-state estimate vs. M5 propagated CR3BP arrival arcs\n"
        "M5 candidates are only shown where dynamically evaluated (Stage-A phase scan, T_back=10 days)",
        fontsize=10,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.88])
    fig.savefig("figures/m5_fig2_m4_vs_m5.png", dpi=150)
    plt.close(fig)


def make_figure_3(m5):
    dir_data = [row for row in _read_sensitivity_csv() if row["study"] == "direction_offset_theta_deg"]
    conv_data = m5["convergence_study"]

    fig, (ax_dir, ax_conv) = plt.subplots(1, 2, figsize=(11, 4.5))

    thetas = [float(row["value"]) for row in dir_data]
    dvs = [float(row["delta_v_m_s"]) for row in dir_data]
    accepted = [row["accepted"] == "True" for row in dir_data]
    colors = ["tab:green" if a else "tab:red" for a in accepted]
    ax_dir.plot(thetas, dvs, color="0.6", lw=1, zorder=1)
    ax_dir.scatter(thetas, dvs, c=colors, zorder=2, s=50)
    ax_dir.axvline(0, color="0.5", ls="--", lw=1, label="baseline direction")
    ax_dir.set_xlabel("direction offset from baseline [deg]")
    ax_dir.set_ylabel("insertion Delta-v [m/s]")
    ax_dir.set_title("Direction sensitivity (green=accepted, red=rejected)", fontsize=9.5)
    ax_dir.legend(fontsize=8)
    ax_dir.grid(True, alpha=0.3)

    tol_names = [row["tolerance"] for row in conv_data]
    round_trip_pos = [row["round_trip_pos_err_nondim"] * c.DU_KM for row in conv_data]
    jacobi_drift = [row["jacobi_max_drift"] for row in conv_data]
    x = np.arange(len(tol_names))
    ax_conv.semilogy(x, round_trip_pos, marker="o", color="tab:blue", label="round-trip position error [km]")
    ax_conv.semilogy(x, jacobi_drift, marker="s", color="tab:orange", label="Jacobi max drift (nondim)")
    ax_conv.set_xticks(x)
    ax_conv.set_xticklabels(tol_names)
    ax_conv.set_xlabel("integration tolerance")
    ax_conv.set_ylabel("error (log scale)")
    ax_conv.set_title("Convergence: round-trip & Jacobi drift vs. tolerance", fontsize=9.5)
    ax_conv.legend(fontsize=8)
    ax_conv.grid(True, which="both", alpha=0.3)

    fig.suptitle("M5 sensitivity and numerical convergence", fontsize=10.5)
    fig.tight_layout(rect=[0, 0, 1, 0.90])
    fig.savefig("figures/m5_fig3_sensitivity_convergence.png", dpi=150)
    plt.close(fig)


def _read_sensitivity_csv():
    with open("results/m5_sensitivity.csv") as f:
        return list(csv.DictReader(f))


if __name__ == "__main__":
    m3, m4, m5 = _load()
    make_figure_1(m3, m5)
    make_figure_2(m4, m5)
    make_figure_3(m5)
    print("Wrote figures/m5_fig1_arrival_arc.png")
    print("Wrote figures/m5_fig2_m4_vs_m5.png")
    print("Wrote figures/m5_fig3_sensitivity_convergence.png")
