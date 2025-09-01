#!/usr/bin/env python3
# Apache-2.0
"""
Push–Pull (Futile) Cycle driven by ATP/ADP using log-linear reaction quotient dynamics.

We model a two-state substrate S <-> S~P where phosphorylation is driven by the
ATP/ADP ratio (chemostatted) and dephosphorylation is passive (water/Pi chemostatted).
Following the repo's framework:
    d/dt ln Q = -k * ln(Q/Keq) + u(t)
Here we take Q to be the ratio r = [S~P]/[S] (i.e., Q ≡ r). Chemostats (ATP/ADP)
enter through the control input u(t) = ln([ATP]/[ADP]) + beta, where 'beta'
captures standard-state free-energy bias and any fixed chemostats (e.g., Pi).

At steady state: ln(Q*/Keq) = u/k  ⇒  Q* = Keq * exp(u/k).
With Keq=1 and u = ln(ATP/ADP) + ln(r0), we obtain:
    r* = r0 * (ATP/ADP)
and the phosphorylated fraction is f* = r*/(1 + r*), a logistic curve in log-space.

This script:
  1) Plots f* vs ATP/ADP over several decades (log x-axis).
  2) Simulates a step in ATP/ADP (1 → 10) and plots f(t).

Requires: numpy, scipy, matplotlib (as in the repo README).
"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# ---------- Parameters (edit here to explore) ----------
k = 1.0            # relaxation rate (1/s) for ln Q
Keq = 1.0          # equilibrium constant for Q=r=[S~P]/[S] absent drive
r0 = 0.1           # baseline ratio r at ATP/ADP = 1 (encodes standard-state bias)
                   # i.e., when ATP/ADP = 1, r* = r0 and f* = r0/(1+r0)

# Range of ATP/ADP to explore for steady-state curve
ATP_ADP_grid = np.logspace(-1, 2, 200)  # 0.1 to 100

# Step response settings
t_end = 10.0       # seconds
t_step = 2.0       # time of step
ATP_ADP_lo = 1.0   # before step
ATP_ADP_hi = 10.0  # after step
# ------------------------------------------------------


def u_from_atp_adp(atp_adp: float) -> float:
    """
    Control input u(t) entering d/dt ln Q = -k ln(Q/Keq) + u.
    We fold fixed standard-state and Pi effects into ln(r0).
    """
    return np.log(atp_adp) + np.log(r0)  # = ln(ATP/ADP) + beta, with beta = ln r0


def f_from_r(r: np.ndarray) -> np.ndarray:
    """Phosphorylated fraction f = r/(1+r)."""
    return r / (1.0 + r)


def steady_state_fraction_curve():
    """Compute steady-state f* vs ATP/ADP grid and save a figure."""
    # At steady state with Keq=1: r* = exp(u/k) = r0 * (ATP/ADP)^(1/k)
    # With k=1, slope in ln-space is unity. If you change k, the slope scales by 1/k.
    r_ss = np.exp(u_from_atp_adp(ATP_ADP_grid) / k) * (Keq ** (1.0 / k))
    f_ss = f_from_r(r_ss)

    plt.figure(figsize=(6, 4))
    plt.semilogx(ATP_ADP_grid, f_ss, linewidth=2)
    plt.xlabel("ATP/ADP ratio")
    plt.ylabel("Phosphorylated fraction  f*")
    plt.title("Push–Pull Cycle: steady-state f* vs ATP/ADP")
    plt.grid(True, which="both", linestyle=":")
    plt.tight_layout()
    plt.savefig("push_pull_fraction_vs_ATP.png", dpi=180)
    plt.close()


def simulate_step_response():
    """
    Simulate a step in ATP/ADP:
      ATP/ADP = 1 for t < t_step, then 10 for t >= t_step.
    State x = ln(Q/Keq) with Q ≡ r; dynamics x' = -k x + u(t).
    """
    def u_of_t(t):
        return u_from_atp_adp(ATP_ADP_lo if t < t_step else ATP_ADP_hi)

    def rhs(t, x):
        return -k * x + u_of_t(t)

    # Initial condition: start at steady state for ATP/ADP_lo
    x0 = u_from_atp_adp(ATP_ADP_lo) / k  # since x* = u/k for Keq=1
    sol = solve_ivp(rhs, (0.0, t_end), [x0], max_step=0.02, dense_output=True)

    tt = np.linspace(0.0, t_end, 500)
    x = sol.sol(tt)[0]
    r = np.exp(x) * Keq
    f = f_from_r(r)

    # Plot fraction vs time
    plt.figure(figsize=(6, 4))
    plt.plot(tt, f, linewidth=2)
    plt.axvline(t_step, linestyle="--")
    plt.xlabel("Time (s)")
    plt.ylabel("Phosphorylated fraction  f(t)")
    plt.title(f"Push–Pull Cycle: step ATP/ADP {ATP_ADP_lo} → {ATP_ADP_hi}")
    plt.grid(True, linestyle=":")
    plt.tight_layout()
    plt.savefig("push_pull_step_response.png", dpi=180)
    plt.close()


if __name__ == "__main__":
    steady_state_fraction_curve()
    simulate_step_response()
    print("Saved figures:\n  - push_pull_fraction_vs_ATP.png\n  - push_pull_step_response.png")

