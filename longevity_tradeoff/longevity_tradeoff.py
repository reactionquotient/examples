"""
Longevity tradeoff (AMPK/NAD+ vs mTOR/nutrients)
------------------------------------------------
Two log-linear reaction-quotient modes:
  Q1 ~ cellular maintenance / autophagy (longevity-leaning)
  Q2 ~ anabolic growth / mTOR activity (growth-leaning)

Dynamics in log-space:
  d/dt ln Q = -K ln(Q/Keq) + u
with constant control inputs u determined by "energy state" (NAD+/NADH)
and "nutrient signal" (e.g., insulin/IGF/mTOR drive).

We compare two steady drive conditions:
  1) Ad libitum (high nutrients, lower NAD+/NADH)
  2) Caloric restriction (lower nutrients, higher NAD+/NADH)

This is a didactic toy model—parameters are arbitrary and chosen to
illustrate the qualitative tradeoff, not to recapitulate specific biology.

Requires: numpy, scipy, matplotlib
Run:     python longevity_tradeoff/longevity_tradeoff.py
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import expm

# ----- Model setup -----

# Rate/coupling matrix K (must be positive-stable for nice behavior)
# Diagonals are "relaxation rates"; small symmetric off-diagonal coupling.
k1 = 1.0   # maintenance/autophagy mode
k2 = 1.2   # growth/mTOR mode
c  = 0.20  # weak coupling
K = np.array([[k1, c],
              [c,  k2]])

# Equilibrium reaction quotients (dimensionless)
Keq = np.array([1.0, 1.0])

# Initial condition in log-space: at equilibrium (ln(Q/Keq) = 0)
x0 = np.zeros(2)

# Helper to turn “state” dials into control input u
def control_u(nad_ratio, nutrient, k_nad=0.6, k_nutr=0.8, k_ampk=0.3):
    """
    u1 increases with NAD+/NADH (maintenance drive via AMPK/Sirtuins).
    u2 increases with nutrients (growth drive via mTOR) and is opposed by NAD+/NADH.
    All signals are log-transformed to match the framework.
    """
    u1 = k_nad * np.log(nad_ratio)
    u2 = k_nutr * np.log(nutrient) - k_ampk * np.log(nad_ratio)
    return np.array([u1, u2])

# Two simple operating conditions
conditions = {
    "Ad libitum": {
        "nad_ratio": 1.2,  # lower NAD+/NADH (fed)
        "nutrient":  5.0   # strong nutrient signal
    },
    "Caloric restriction": {
        "nad_ratio": 3.0,  # higher NAD+/NADH (fasted/CR)
        "nutrient":  1.2   # weak nutrient signal
    }
}

# ----- Closed-form trajectory for constant u -----
# x(t) = expm(-K t) x0 + K^{-1} (I - expm(-K t)) u
# where x = ln(Q/Keq). We only need exp(x) = Q/Keq for plotting.
K_inv = np.linalg.inv(K)

def simulate(K, x0, u, t_grid):
    X = np.zeros((len(t_grid), len(x0)))
    I = np.eye(K.shape[0])
    for i, t in enumerate(t_grid):
        A = expm(-K * t)
        X[i] = A @ x0 + K_inv @ (I - A) @ u
    return X  # log-space

# ----- Run + plot -----
t = np.linspace(0, 10, 200)
fig, ax = plt.subplots(figsize=(7.5, 4.8))

colors = {'maintenance': '#1f77b4', 'growth': '#ff7f0e'}  # blue for maintenance, orange for growth
linestyles = {'Ad libitum': '-', 'Caloric restriction': '--'}
for label, pars in conditions.items():
    u = control_u(pars["nad_ratio"], pars["nutrient"])
    X = simulate(K, x0, u, t)
    R = np.exp(X)  # R = Q/Keq (dimensionless)
    ax.plot(t, R[:, 0], lw=2, color=colors['maintenance'], ls=linestyles[label], label=f"{label} — maintenance (Q1/Keq)")
    ax.plot(t, R[:, 1], lw=2, color=colors['growth'], ls=linestyles[label], label=f"{label} — growth (Q2/Keq)")

ax.set_xlabel("time (a.u.)")
ax.set_ylabel("relative reaction quotient (Q/Keq)")
ax.set_title("Longevity tradeoff: maintenance vs growth (log-linear RQ dynamics)")
ax.legend(ncol=2, fontsize=9)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("longevity_tradeoff.png", dpi=160)
plt.show()

