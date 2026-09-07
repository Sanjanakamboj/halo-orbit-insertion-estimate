"""Generate the M4 diagnostic figures.

Figure 1: insertion Delta-v vs. halo phase (headline M4 figure).
Figure 2: halo orbit with insertion point and velocity vectors.
Figure 3: arrival-speed and arrival-direction sensitivity.

Run from the repository root (after scripts/build_m4_insertion.py):
    python scripts/make_m4_figures.py
"""

import csv
import json

import matplotlib.pyplot as plt
import numpy as np

from halo_insertion import constants as c
from halo_insertion.equilibria import find_L2
from halo_insertion.propagation import propagate

MU = c.MU


def _load():
    with open("results/m4_insertion_summary.json") as f:
        summary = json.load(f)

    taus, dvs, valid = [], [], []
    with open("results/m4_phase_sweep.csv") as f:
        reader = csv.DictReader(f)
        for row in reader:
            taus.append(float(row["tau"]))
            dvs.append(float(row["delta_v_m_s"]))
            valid.append(row["valid"] == "True")

    return summary, np.array(taus), np.array(dvs), np.array(valid)


def make_figure_1(summary, taus, dvs, valid):
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(taus[valid], dvs[valid], lw=1.3, color="tab:blue")

    sel = summary["selected_insertion"]
    ax.scatter([sel["tau"]], [sel["delta_v_m_s"]], color="tab:red", zorder=5, s=60,
               label=f"selected minimum: {sel['delta_v_m_s']:.1f} m/s @ tau={sel['tau']:.3f}")

    for lm in summary["phase_sweep"]["local_minima"]:
        ax.scatter([lm["tau"]], [lm["delta_v_m_s"]], color="0.3", marker="x", zorder=4, s=40)

    ax.axhline(summary["degeneracy_guard_aligned_case"]["delta_v_m_s"], color="tab:gray", ls="--", lw=1,
                label=f"idealized aligned-direction lower bound: {summary['degeneracy_guard_aligned_case']['delta_v_m_s']:.1f} m/s (nonphysical)")

    ax.set_xlabel("orbital phase  tau = t / T")
    ax.set_ylabel("insertion Delta-v [m/s]")
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(True, alpha=0.3)

    fig.suptitle(
        "M4 halo-insertion velocity mismatch vs. orbital phase\n"
        "Engineering arrival-state assumption; not an optimized Earth-to-halo transfer",
        fontsize=10.5,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.90])
    fig.savefig("figures/m4_fig1_delta_v_vs_phase.png", dpi=150)
    plt.close(fig)


def make_figure_2(summary):
    l2 = find_L2(MU)
    s0 = summary["m3_halo_reference"]["corrected_initial_state"]
    state0 = np.array([s0["x0"], s0["y0"], s0["z0"], s0["xdot0"], s0["ydot0"], s0["zdot0"]])
    period = summary["m3_halo_reference"]["period_nondim"]

    t_eval = np.linspace(0, period, 2000)
    res = propagate(state0, (0, period), MU, t_eval=t_eval)
    x, y, z = res.y[0], res.y[1], res.y[2]

    sel = summary["selected_insertion"]
    r_h = np.array(sel["position_nondim"])
    v_halo = np.array(sel["v_halo_nondim"])
    v_arr = np.array(sel["v_arr_nondim"])
    dv_vec = np.array(sel["delta_v_vec_nondim"])

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.2))

    # Visual velocity-vector scale (position is in DU; velocities are DU/TU --
    # scaled for visibility, explicitly labeled as such).
    vec_scale = 0.15

    for ax, (i, j), xlabel, ylabel in [
        (axes[0], (0, 1), "x (DU)", "y (DU)"),
        (axes[1], (0, 2), "x (DU)", "z (DU)"),
    ]:
        a, b = [x, y, z][i], [x, y, z][j]
        ax.plot(a, b, lw=1.0, color="0.7", label="halo orbit", zorder=1)
        ax.scatter([1 - MU], [0], color="0.4", s=50, zorder=2, label="Moon" if ax is axes[0] else None)
        ax.scatter([l2.x], [0], color="tab:red", marker="x", s=50, zorder=2, label="L2" if ax is axes[0] else None)
        ax.scatter([r_h[i]], [r_h[j]], color="black", s=40, zorder=3, label="insertion point" if ax is axes[0] else None)

        ax.quiver(r_h[i], r_h[j], v_halo[i] * vec_scale, v_halo[j] * vec_scale,
                   angles="xy", scale_units="xy", scale=1, color="tab:green", width=0.008,
                   label="v_halo" if ax is axes[0] else None)
        ax.quiver(r_h[i], r_h[j], v_arr[i] * vec_scale, v_arr[j] * vec_scale,
                   angles="xy", scale_units="xy", scale=1, color="tab:orange", width=0.008,
                   label="v_arrival" if ax is axes[0] else None)
        ax.quiver(r_h[i] + v_arr[i] * vec_scale, r_h[j] + v_arr[j] * vec_scale,
                   dv_vec[i] * vec_scale, dv_vec[j] * vec_scale,
                   angles="xy", scale_units="xy", scale=1, color="tab:red", width=0.008,
                   label="Delta-v" if ax is axes[0] else None)

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_aspect("equal", adjustable="datalim")

    axes[0].set_title("x-y projection", fontsize=10)
    axes[1].set_title("x-z projection", fontsize=10)
    axes[0].legend(fontsize=7, loc="upper left")

    fig.suptitle(
        f"M3 halo orbit with selected insertion point and velocity vectors\n"
        f"Velocity vectors scaled x{1/vec_scale:.1f} for visibility -- NOT to position scale; "
        f"Delta-v = {sel['delta_v_m_s']:.1f} m/s",
        fontsize=10,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.88])
    fig.savefig("figures/m4_fig2_insertion_point_vectors.png", dpi=150)
    plt.close(fig)


def make_figure_3(summary):
    speed_data = summary["sensitivity"]["arrival_speed_dC"]
    dir_data = summary["sensitivity"]["arrival_direction_theta_deg"]

    fig, (ax_speed, ax_dir) = plt.subplots(1, 2, figsize=(11, 4.5))

    dCs = [row["dC"] for row in speed_data]
    dvs = [row["min_delta_v_m_s"] for row in speed_data]
    ax_speed.plot(dCs, dvs, marker="o", color="tab:blue")
    baseline_dC = summary["arrival_model"]["dC_baseline"]
    ax_speed.axvline(baseline_dC, color="0.5", ls="--", lw=1, label=f"baseline dC={baseline_dC}")
    ax_speed.set_xlabel("arrival-speed assumption  dC  (nondim Jacobi offset)")
    ax_speed.set_ylabel("minimum insertion Delta-v [m/s]")
    ax_speed.set_title("Arrival-speed sensitivity", fontsize=10)
    ax_speed.legend(fontsize=8)
    ax_speed.grid(True, alpha=0.3)

    thetas = [row["theta_deg"] for row in dir_data]
    dvs2 = [row["min_delta_v_m_s"] for row in dir_data]
    ax_dir.plot(thetas, dvs2, marker="o", color="tab:orange")
    ax_dir.axvline(0, color="0.5", ls="--", lw=1, label="baseline direction (radial from Earth)")
    ax_dir.set_xlabel("direction offset from baseline [deg, about z-axis]")
    ax_dir.set_ylabel("minimum insertion Delta-v [m/s]")
    ax_dir.set_title("Arrival-direction sensitivity", fontsize=10)
    ax_dir.legend(fontsize=8)
    ax_dir.grid(True, alpha=0.3)

    fig.suptitle(
        "M4 sensitivity to arrival-state assumptions\n"
        "Jacobi-consistent speed model; radial-from-Earth baseline direction",
        fontsize=10.5,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.88])
    fig.savefig("figures/m4_fig3_sensitivity.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    summary, taus, dvs, valid = _load()
    make_figure_1(summary, taus, dvs, valid)
    make_figure_2(summary)
    make_figure_3(summary)
    print("Wrote figures/m4_fig1_delta_v_vs_phase.png")
    print("Wrote figures/m4_fig2_insertion_point_vectors.png")
    print("Wrote figures/m4_fig3_sensitivity.png")
