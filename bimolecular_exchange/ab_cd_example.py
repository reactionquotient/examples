#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A + B <-> C + D  (bimolecular exchange) under log-linear reaction-quotient dynamics.

We evolve the reaction quotient Q(t) toward equilibrium with:
    d/dt ln Q = -k * ln(Q / Keq)
which has the closed-form solution:
    Q(t) = Keq * (Q0/Keq) ** exp(-k t)
(see repo README)                                         

We then recover species concentrations using the closed-system extent of reaction ξ:
    A(t) = A0 - ξ(t)
    B(t) = B0 - ξ(t)
    C(t) = C0 + ξ(t)
    D(t) = D0 + ξ(t)

Given Q(t) = (C D)/(A B) and the relations above, ξ(t) satisfies a quadratic
for each time t; we choose the physically valid, continuous root.

Optionally, we compare to classical mass-action kinetics:
    v_f = kf * A * B,  v_r = kr * C * D
    dA/dt = -v_f + v_r, etc., with Keq = kf/kr.

Requires: numpy, scipy, matplotlib
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# --------------------------
# Parameters (feel free to edit)
# --------------------------
# Log-linear dynamics
k = 1.0               # relaxation rate (1/time)
Keq = 10.0            # equilibrium constant for Q = (C D)/(A B)

# Initial concentrations (must be nonnegative)
A0, B0, C0, D0 = 5.0, 4.0, 1.0, 0.5

# Time grid
t_end = 10.0
N = 400
t = np.linspace(0.0, t_end, N)

# Mass-action comparison (set to True to overlay)
COMPARE_MASS_ACTION = True
kf = 1.0              # forward rate constant (for A+B -> C+D)
kr = kf / Keq         # choose kr so that Keq = kf/kr

# --------------------------
# Helper: compute ξ from Q
# --------------------------
def extent_from_Q_series(Q_series, A0, B0, C0, D0, eps=1e-12):
    """
    Solve (C0+ξ)(D0+ξ) = Q (A0-ξ)(B0-ξ) for ξ at each Q.
    Quadratic: (Q-1) ξ^2 + [ -Q(A0+B0) - (C0+D0) ] ξ + (Q A0 B0 - C0 D0) = 0

    We pick the physically valid root (all concentrations >= 0) and enforce
    continuity by staying closest to the previous ξ.
    """
    xi = np.zeros_like(Q_series, dtype=float)
    xi_prev = 0.0

    for i, Q in enumerate(Q_series):
        a = Q - 1.0
        b = -Q * (A0 + B0) - (C0 + D0)
        c = Q * A0 * B0 - C0 * D0

        if abs(a) < eps:
            # Degenerates to linear: b ξ + c = 0
            xi_candidate = -c / b if abs(b) > eps else 0.0
            roots = [xi_candidate]
        else:
            disc = b*b - 4.0*a*c
            # Guard against tiny negative due to FP errors
            disc = max(disc, 0.0)
            s = np.sqrt(disc)
            roots = [(-b + s) / (2.0*a), (-b - s) / (2.0*a)]

        # Filter roots by nonnegativity of concentrations
        valid = []
        for r in roots:
            A = A0 - r
            B = B0 - r
            C = C0 + r
            D = D0 + r
            if (A >= -1e-10) and (B >= -1e-10) and (C >= -1e-10) and (D >= -1e-10):
                valid.append(r)

        if len(valid) == 0:
            # Fallback: pick root that least violates nonnegativity
            # (should be rare and only due to numerical edge cases)
            penalties = []
            for r in roots:
                A = A0 - r
                B = B0 - r
                C = C0 + r
                D = D0 + r
                penalties.append(
                    sum(max(-A, 0),) + sum(max(-B, 0),) + sum(max(-C, 0),) + sum(max(-D, 0),)
                )
            r = roots[int(np.argmin(penalties))]
        else:
            # Enforce continuity
            r = min(valid, key=lambda z: abs(z - xi_prev))

        xi[i] = r
        xi_prev = r

    return xi

# --------------------------
# Log-linear Q(t) and concentrations
# --------------------------
Q0 = (C0 * D0) / (A0 * B0)
Q = Keq * (Q0 / Keq) ** np.exp(-k * t)  # closed-form solution

xi = extent_from_Q_series(Q, A0, B0, C0, D0)
A = A0 - xi
B = B0 - xi
C = C0 + xi
D = D0 + xi

# --------------------------
# (Optional) Mass-action comparison
# --------------------------
if COMPARE_MASS_ACTION:
    def rhs(_t, y):
        A, B, C, D = y
        vf = kf * A * B
        vr = kr * C * D
        dA = -vf + vr
        dB = -vf + vr
        dC =  vf - vr
        dD =  vf - vr
        return [dA, dB, dC, dD]

    y0 = [A0, B0, C0, D0]
    sol = solve_ivp(rhs, (0.0, t_end), y0, t_eval=t, rtol=1e-8, atol=1e-10)
    A_ma, B_ma, C_ma, D_ma = sol.y

# --------------------------
# Plot
# --------------------------
plt.figure(figsize=(7, 4.8))
colors = ['C0', 'C1', 'C2', 'C3']
plt.plot(t, A, color=colors[0], label='A (log-linear)')
plt.plot(t, B, color=colors[1], label='B (log-linear)')
plt.plot(t, C, color=colors[2], label='C (log-linear)')
plt.plot(t, D, color=colors[3], label='D (log-linear)')

if COMPARE_MASS_ACTION:
    plt.plot(t, A_ma, '--', color=colors[0], label='A (mass-action)')
    plt.plot(t, B_ma, '--', color=colors[1], label='B (mass-action)')
    plt.plot(t, C_ma, '--', color=colors[2], label='C (mass-action)')
    plt.plot(t, D_ma, '--', color=colors[3], label='D (mass-action)')

plt.xlabel('time')
plt.ylabel('concentration')
plt.title('A + B ⇌ C + D: log-linear dynamics (± mass-action)')
plt.legend(ncol=2, fontsize=9)
plt.tight_layout()
plt.savefig('ab_cd_reaction.png', dpi=150)
print("Saved figure: ab_cd_reaction.png")

# Also visualize Q(t) itself
plt.figure(figsize=(6.2, 4.2))
plt.plot(t, Q, label='Q(t)')
plt.axhline(Keq, linestyle=':', label='Keq')
plt.xlabel('time')
plt.ylabel('Q = (C·D)/(A·B)')
plt.title('Reaction quotient relaxation')
plt.legend()
plt.tight_layout()
plt.savefig('ab_cd_Q_relaxation.png', dpi=150)
print("Saved figure: ab_cd_Q_relaxation.png")

