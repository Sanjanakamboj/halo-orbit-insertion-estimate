"""Generate the M2 diagnostic figures.

Figure 1: Earth-Moon rotating-frame equilibrium geometry (L1/L2/L3).
Figure 2: M2 verification trajectory (NOT a halo orbit) + Jacobi error.

Run from the repository root:
    python scripts/make_m2_figures.py
"""

import numpy as np
import matplotlib.pyplot as plt

from halo_insertion import constants as c
from halo_insertion.cr3bp import jacobi_constant
from halo_insertion.equilibria import find_L1, find_L2, find_L3
from halo_insertion.normalization import days_to_nondim_time
from halo_insertion.propagation import propagate

MU = c.MU


def make_figure_1():
    l1, l2, l3 = find_L1(MU), find_L2(MU), find_L3(MU)
    earth_x, moon_x = -MU, 1.0 - MU

    fig, (ax_full, ax_zoom) = plt.subplots(1, 2, figsize=(11, 4.2))

    points = [
        ("Earth", earth_x, "tab:blue", 90),
        ("Moon", moon_x, "0.4", 50),
        ("L1", l1.x, "tab:red", 35),
        ("L2", l2.x, "tab:red", 35),
        ("L3", l3.x, "tab:red", 35),
    ]

    for ax, title in [(ax_full, "Full system"), (ax_zoom, "Earth-Moon zoom")]:
        for name, x, color, size in points:
            ax.scatter([x], [0], s=size, color=color, zorder=3)
            ax.annotate(
                name,
                (x, 0),
                textcoords="offset points",
                xytext=(0, 10),
                ha="center",
                fontsize=9,
            )
        ax.axhline(0, color="0.8", lw=0.8, zorder=1)
        ax.set_xlabel("x (DU, synodic rotating frame)")
        ax.set_yticks([])
        ax.set_title(title, fontsize=10)

    ax_full.set_xlim(l3.x - 0.3, l2.x + 0.3)
    ax_zoom.set_xlim(earth_x - 0.15, l2.x + 0.15)

    fig.suptitle(
        "CR3BP equilibrium geometry — diagnostic, not a halo orbit\n"
        "(collinear points L1/L2/L3 are equilibria; L1/L2/L3 are NOT dynamically stable)",
        fontsize=10.5,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.88])
    fig.savefig("figures/m2_fig1_equilibrium_geometry.png", dpi=150)
    plt.close(fig)


def make_figure_2():
    l2 = find_L2(MU)
    state0 = np.array([l2.x - 0.02, 0.01, 0.005, 0.0, 0.01, 0.005])
    t_end = days_to_nondim_time(10.0)

    n_eval = 2000
    t_eval = np.linspace(0.0, t_end, n_eval)
    result = propagate(state0, (0.0, t_end), MU, t_eval=t_eval)
    assert result.success

    C0 = jacobi_constant(state0, MU)
    C_t = np.array([jacobi_constant(result.y[:, i], MU) for i in range(result.y.shape[1])])
    jacobi_err = np.abs(C_t - C0)
    jacobi_err_floor = np.maximum(jacobi_err, 1e-16)  # avoid log(0)

    fig, (ax_traj, ax_jacobi) = plt.subplots(1, 2, figsize=(11, 4.2))

    ax_traj.plot(result.y[0], result.y[1], lw=1.2, color="tab:purple")
    ax_traj.scatter([state0[0]], [state0[1]], color="tab:green", zorder=3, label="start")
    ax_traj.scatter([l2.x], [0], color="tab:red", marker="x", zorder=3, label="L2")
    ax_traj.set_xlabel("x (DU)")
    ax_traj.set_ylabel("y (DU)")
    ax_traj.set_title("x-y rotating-frame trajectory", fontsize=10)
    ax_traj.legend(fontsize=8, loc="best")
    ax_traj.set_aspect("equal", adjustable="datalim")

    t_days = result.t * c.TU_DAYS
    ax_jacobi.semilogy(t_days, jacobi_err_floor, lw=1.2, color="tab:orange")
    ax_jacobi.set_xlabel("time (days)")
    ax_jacobi.set_ylabel("|C(t) - C(0)|")
    ax_jacobi.set_title("Jacobi conservation error (log scale)", fontsize=10)
    ax_jacobi.grid(True, which="both", alpha=0.3)

    fig.suptitle(
        "M2 CR3BP verification trajectory — NOT a halo orbit\n"
        "(benign near-L2 test state; propagator/Jacobi diagnostic only)",
        fontsize=10.5,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.88])
    fig.savefig("figures/m2_fig2_verification_trajectory_jacobi.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    make_figure_1()
    make_figure_2()
    print("Wrote figures/m2_fig1_equilibrium_geometry.png")
    print("Wrote figures/m2_fig2_verification_trajectory_jacobi.png")
