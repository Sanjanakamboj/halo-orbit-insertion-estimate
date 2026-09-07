"""Generate the M3 diagnostic figures.

Figure 1: corrected L2 halo orbit, 3D.
Figure 2: rotating-frame x-y / x-z / y-z projections.
Figure 3: differential-correction convergence + periodicity/Jacobi verification.

Run from the repository root (after scripts/build_m3_halo.py):
    python scripts/make_m3_figures.py
"""

import json

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

from halo_insertion import constants as c
from halo_insertion.cr3bp import jacobi_constant
from halo_insertion.equilibria import find_L2
from halo_insertion.propagation import propagate

MU = c.MU


def _load_summary():
    with open("results/m3_halo_summary.json") as f:
        return json.load(f)


def _state0_and_period(summary):
    s0 = summary["corrected_initial_state"]
    state0 = np.array([s0["x0"], s0["y0"], s0["z0"], s0["xdot0"], s0["ydot0"], s0["zdot0"]])
    period = summary["period_nondim"]
    return state0, period


def make_figure_1(summary):
    state0, period = _state0_and_period(summary)
    l2 = find_L2(MU)
    t_eval = np.linspace(0, period, 3000)
    res = propagate(state0, (0, period), MU, t_eval=t_eval)
    x, y, z = res.y[0], res.y[1], res.y[2]

    fig = plt.figure(figsize=(8, 7))
    ax = fig.add_subplot(111, projection="3d")
    ax.plot(x, y, z, lw=1.3, color="tab:purple", label="corrected halo orbit")
    ax.scatter([1.0 - MU], [0], [0], color="0.4", s=60, label="Moon")
    ax.scatter([l2.x], [0], [0], color="tab:red", marker="x", s=60, label="L2")
    ax.scatter([state0[0]], [state0[1]], [state0[2]], color="tab:green", s=50, label="initial point (x-z crossing)")

    ax.set_xlabel("x (DU)")
    ax.set_ylabel("y (DU)")
    ax.set_zlabel("z (DU)")

    # Equal-ish aspect ratio so geometry isn't visually distorted.
    max_range = np.array([x.max() - x.min(), y.max() - y.min(), z.max() - z.min()]).max() / 2.0
    mid_x, mid_y, mid_z = (x.max() + x.min()) / 2, (y.max() + y.min()) / 2, (z.max() + z.min()) / 2
    ax.set_xlim(mid_x - max_range, mid_x + max_range)
    ax.set_ylim(mid_y - max_range, mid_y + max_range)
    ax.set_zlim(mid_z - max_range, mid_z + max_range)

    ax.legend(fontsize=8, loc="upper left")
    fig.suptitle(
        "M3 corrected Earth-Moon L2 halo orbit — CR3BP\n"
        "Numerically differential-corrected periodic orbit; not an ephemeris trajectory",
        fontsize=10.5,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig("figures/m3_fig1_halo_orbit_3d.png", dpi=150)
    plt.close(fig)


def make_figure_2(summary):
    state0, period = _state0_and_period(summary)
    l2 = find_L2(MU)
    t_eval = np.linspace(0, period, 3000)
    res = propagate(state0, (0, period), MU, t_eval=t_eval)
    x, y, z = res.y[0], res.y[1], res.y[2]

    fig, axes = plt.subplots(1, 3, figsize=(13, 4.3))

    panels = [
        (axes[0], x, y, "x (DU)", "y (DU)", "x-y projection"),
        (axes[1], x, z, "x (DU)", "z (DU)", "x-z projection"),
        (axes[2], y, z, "y (DU)", "z (DU)", "y-z projection"),
    ]
    for ax, a, b, xlabel, ylabel, title in panels:
        ax.plot(a, b, lw=1.2, color="tab:purple")
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title, fontsize=10)
        ax.set_aspect("equal", adjustable="datalim")

    axes[0].scatter([l2.x], [0], color="tab:red", marker="x", zorder=3, label="L2")
    axes[0].scatter([1 - MU], [0], color="0.4", zorder=3, label="Moon")
    axes[0].scatter([state0[0]], [state0[1]], color="tab:green", zorder=3, label="start")
    axes[0].legend(fontsize=7, loc="best")

    axes[1].scatter([l2.x], [0], color="tab:red", marker="x", zorder=3)
    axes[1].scatter([1 - MU], [0], color="0.4", zorder=3)
    axes[1].scatter([state0[0]], [state0[2]], color="tab:green", zorder=3)

    axes[2].scatter([0], [0], color="tab:red", marker="x", zorder=3)
    axes[2].scatter([state0[1]], [state0[2]], color="tab:green", zorder=3)

    fig.suptitle(
        "M3 corrected L2 halo orbit — rotating-frame projections\n"
        "Numerically differential-corrected periodic orbit; not an ephemeris trajectory",
        fontsize=10.5,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.88])
    fig.savefig("figures/m3_fig2_projections.png", dpi=150)
    plt.close(fig)


def make_figure_3(summary):
    import csv

    iterations, residuals = [], []
    with open("results/m3_correction_history.csv") as f:
        reader = csv.DictReader(f)
        for row in reader:
            iterations.append(int(row["iteration"]))
            xdot_r = float(row["xdot_residual"])
            zdot_r = float(row["zdot_residual"])
            residuals.append(np.hypot(xdot_r, zdot_r))
    residuals = np.array(residuals)
    residuals_floor = np.maximum(residuals, 1e-16)

    state0, period = _state0_and_period(summary)
    t_eval = np.linspace(0, period, 3000)
    res = propagate(state0, (0, period), MU, t_eval=t_eval)
    C0 = jacobi_constant(state0, MU)
    C_t = np.array([jacobi_constant(res.y[:, i], MU) for i in range(res.y.shape[1])])
    jacobi_err = np.maximum(np.abs(C_t - C0), 1e-16)

    fig, (ax_conv, ax_jacobi) = plt.subplots(1, 2, figsize=(11, 4.3))

    ax_conv.semilogy(iterations, residuals_floor, marker="o", color="tab:blue")
    ax_conv.set_xlabel("Newton iteration")
    ax_conv.set_ylabel("half-period residual norm  |[xdot, zdot](T/2)|")
    ax_conv.set_title("Differential-correction convergence (solver)", fontsize=10)
    ax_conv.grid(True, which="both", alpha=0.3)

    t_days = t_eval * c.TU_DAYS
    ax_jacobi.semilogy(t_days, jacobi_err, lw=1.2, color="tab:orange")
    ax_jacobi.set_xlabel("time (days)")
    ax_jacobi.set_ylabel("|C(t) - C(0)|")
    ax_jacobi.set_title("Jacobi conservation over one period (propagation)", fontsize=10)
    ax_jacobi.grid(True, which="both", alpha=0.3)

    fig.suptitle(
        "M3 solver convergence vs. propagation conservation — distinct checks\n"
        f"Full-period closure (generic propagator): "
        f"{summary['full_period_closure_generic_propagator']['position_norm_km']:.2e} km",
        fontsize=10,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.88])
    fig.savefig("figures/m3_fig3_convergence_and_verification.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    summary = _load_summary()
    make_figure_1(summary)
    make_figure_2(summary)
    make_figure_3(summary)
    print("Wrote figures/m3_fig1_halo_orbit_3d.png")
    print("Wrote figures/m3_fig2_projections.png")
    print("Wrote figures/m3_fig3_convergence_and_verification.png")
