#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
General single reaction under log-linear reaction-quotient dynamics.

Reaction:
    sum_i alpha[i] * X_i   <->   sum_i beta[i] * X_i
with stoichiometric net coefficients:
    nu[i] = beta[i] - alpha[i]   (products positive, reactants negative)

Core model (same as other examples):
    d/dt ln Q = -k * ln(Q / Keq)
=>  Q(t) = Keq * (Q0/Keq) ** exp(-k t)

Closed system with one extent of reaction ξ:
    c_i(t) = c_i0 + nu[i] * ξ(t)

Given Q(ξ) = prod_i c_i(ξ) ** nu[i], we solve for ξ(t) by root-finding:
    f(ξ) := sum_i nu[i] * ln(c_i0 + nu[i] * ξ) - ln Q_target(t) = 0
f is strictly increasing on the feasible interval (all concentrations > 0),
so a unique solution exists if Q_target is reachable; otherwise we clamp to
the closest feasible boundary.

Optionally, we compare to general mass-action kinetics:
    v_f = kf * Π_i c_i ** alpha[i]
    v_r = kr * Π_i c_i ** beta[i]
    dc_i/dt = nu[i] * (v_f - v_r)
with Keq = kf / kr (ideal activities, elementary orders = stoich).

Requires: numpy, scipy, matplotlib
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq
from scipy.integrate import solve_ivp

# --------------------------
# User parameters
# --------------------------
# Species (names only used for labels)
species = ["A", "B", "C", "D"]

# Stoichiometry (alpha: reactants, beta: products)
# Example: 2A + B <-> C + 2D
alpha = {"A": 2, "B": 1, "C": 0, "D": 0}
beta  = {"A": 0, "B": 0, "C": 1, "D": 2}

# Initial concentrations (must be >=0; positive for any species with |nu|>0)
c0 = {"A": 3.5, "B": 1.5, "C": 0.6, "D": 0.3}

# Log-linear dynamics parameters
k = 1.0       # relaxation rate (1/time)
Keq = 4.0     # equilibrium constant for Q = Π c_i^{nu_i}

# Time grid
t_end = 10.0
N = 400
t = np.linspace(0.0, t_end, N)

# Mass-action comparison toggle & rates
COMPARE_MASS_ACTION = True
kf = 1.0
kr = kf / Keq   # choose so Keq = kf/kr at equilibrium

# --------------------------
# Derived arrays
# --------------------------
# Build aligned vectors for faster math
sp_to_idx = {s: i for i, s in enumerate(species)}
m = len(species)

nu = np.zeros(m, dtype=float)
a = np.zeros(m, dtype=float)  # alpha
b = np.zeros(m, dtype=float)  # beta
c0_vec = np.zeros(m, dtype=float)

for s in species:
    a[sp_to_idx[s]] = float(alpha.get(s, 0))
    b[sp_to_idx[s]] = float(beta.get(s, 0))
    nu[sp_to_idx[s]] = b[sp_to_idx[s]] - a[sp_to_idx[s]]
    c0_vec[sp_to_idx[s]] = float(c0.get(s, 0.0))

# Sanity checks
if np.any((np.abs(nu) > 0) & (c0_vec <= 0)):
    raise ValueError("Initial concentrations must be >0 for species participating in the reaction.")

# --------------------------
# Utilities
# --------------------------
def reaction_str():
    lhs = " + ".join([f"{int(a[i]) if a[i].is_integer() else a[i]} {species[i]}".replace("1 ", "") 
                      for i in range(m) if a[i] > 0])
    rhs = " + ".join([f"{int(b[i]) if b[i].is_integer() else b[i]} {species[i]}".replace("1 ", "") 
                      for i in range(m) if b[i] > 0])
    return f"{lhs} \u21C4 {rhs}"  # ⇄

def feasible_xi_interval(c0v, nuv, eps=1e-12):
    """
    For each i, require c_i = c0_i + nu_i * ξ >= 0.
    If nu_i > 0  => ξ >= -c0_i/nu_i  (lower bound)
    If nu_i < 0  => ξ <=  c0_i/(-nu_i) (upper bound)
    """
    xi_min = -np.inf
    xi_max = +np.inf
    for c0_i, nu_i in zip(c0v, nuv):
        if nu_i > 0:
            xi_min = max(xi_min, -c0_i / nu_i)
        elif nu_i < 0:
            xi_max = min(xi_max, c0_i / (-nu_i))
    # shrink slightly to keep logs finite
    if not np.isneginf(xi_min):
        xi_min += eps
    else:
        xi_min = -1e12
    if not np.isposinf(xi_max):
        xi_max -= eps
    else:
        xi_max = 1e12
    if xi_min > xi_max:
        raise ValueError("Infeasible initial state / stoichiometry (no valid extent interval).")
    return xi_min, xi_max

def S_of_xi(xi, c0v, nuv):
    """S(ξ) = sum_i nu_i * ln(c0_i + nu_i ξ) (requires all terms > 0)."""
    vals = c0v + nuv * xi
    if np.any(vals <= 0):
        return np.nan
    return float(np.sum(nuv * np.log(vals)))

def Q_of_c(c):
    """Q(c) = prod_i c_i ** nu_i, robust to floating errors."""
    return float(np.exp(np.sum(nu * np.log(np.maximum(c, 1e-300)))) if np.all(c > 0) else np.nan)

def xi_from_Q(Q_target, c0v, nuv):
    """
    Solve S(ξ) = ln Q_target for ξ in [xi_min, xi_max].
    If ln Q_target lies outside [S(xi_min), S(xi_max)], clamp to the nearest boundary.
    """
    xi_lo, xi_hi = feasible_xi_interval(c0v, nuv)
    S_lo = S_of_xi(xi_lo, c0v, nuv)
    S_hi = S_of_xi(xi_hi, c0v, nuv)
    lnQ = np.log(Q_target)

    # Monotone: S_lo < S_hi unless degenerate
    if lnQ <= S_lo:
        return xi_lo
    if lnQ >= S_hi:
        return xi_hi

    # Root-finding
    f = lambda x: S_of_xi(x, c0v, nuv) - lnQ
    try:
        return brentq(f, xi_lo, xi_hi, xtol=1e-12, rtol=1e-10, maxiter=200)
    except ValueError:
        # Numerical safety: fall back to midpoint if bracketing fails (should be rare)
        return 0.5 * (xi_lo + xi_hi)

# --------------------------
# Log-linear Q(t) and concentrations
# --------------------------
Q0 = np.prod(c0_vec ** nu)  # requires positive c0 for participating species
Q_t = Keq * (Q0 / Keq) ** np.exp(-k * t)

xi = np.array([xi_from_Q(Qt, c0_vec, nu) for Qt in Q_t])
C = (c0_vec[:, None] + nu[:, None] * xi[None, :])  # shape (m, N)

# --------------------------
# (Optional) Mass-action comparison
# --------------------------
if COMPARE_MASS_ACTION:
    def rhs(_t, y):
        y = np.maximum(y, 0.0)  # prevent tiny negatives
        vf = kf * np.prod(y ** a)
        vr = kr * np.prod(y ** b)
        dy = nu * (vf - vr)
        return dy

    sol = solve_ivp(rhs, (0.0, t_end), c0_vec, t_eval=t, rtol=1e-8, atol=1e-10)
    C_ma = sol.y  # shape (m, N)

# --------------------------
# Plots
# --------------------------
# Define colors for consistent plotting
colors = plt.cm.tab10(np.linspace(0, 1, len(species)))

title = f"{reaction_str()}: log-linear dynamics"
plt.figure(figsize=(7.4, 5.0))
for i, s in enumerate(species):
    plt.plot(t, C[i], label=f"{s} (log-linear)", color=colors[i])
if COMPARE_MASS_ACTION:
    for i, s in enumerate(species):
        plt.plot(t, C_ma[i], "--", label=f"{s} (mass-action)", color=colors[i])
plt.xlabel("time")
plt.ylabel("concentration")
plt.title(title)
plt.legend(ncol=2, fontsize=8)
plt.tight_layout()
plt.savefig("general_single_reaction_conc.png", dpi=150)
print("Saved figure: general_single_reaction_conc.png")

plt.figure(figsize=(6.4, 4.4))
plt.plot(t, Q_t, label="Q(t)")
plt.axhline(Keq, linestyle=":", label="Keq")
plt.xlabel("time")
plt.ylabel("Q = Π c_i^{ν_i}")
plt.title("Reaction quotient relaxation")
plt.legend()
plt.tight_layout()
plt.savefig("general_single_reaction_Q.png", dpi=150)
print("Saved figure: general_single_reaction_Q.png")

