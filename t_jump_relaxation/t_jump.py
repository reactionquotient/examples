# t_jump_relaxation/t_jump.py
import numpy as np
import matplotlib.pyplot as plt

R = 8.314  # J/mol/K

def ln_Keq(T, dH=-50_000.0, dS=-100.0):
    # van ’t Hoff: ln Keq = -(ΔH)/(RT) + ΔS/R
    return (-dH / (R * T)) + (dS / R)

# time and temperature profile: 298 K → 330 K at t=2 s → back at t=7 s
t_end, dt = 12.0, 0.002
t = np.arange(0, t_end + dt, dt)
T = np.where(t < 2.0, 298.0, np.where(t < 7.0, 330.0, 298.0))

lnKeq = ln_Keq(T)
Keq = np.exp(lnKeq)

k = 1.0  # s^-1 relaxation rate (single fit parameter in classic T-jump)
Q = np.empty_like(t)
lnQ = np.empty_like(t)

# start away from equilibrium at initial temperature
Q0 = 0.2 * np.exp(ln_Keq(298.0))
lnQ[0] = np.log(Q0)

# Discrete integration of: d/dt ln Q = -k (ln Q - ln Keq(t))
for i in range(len(t) - 1):
    lnQ[i+1] = lnQ[i] + dt * (-k * (lnQ[i] - lnKeq[i]))

Q = np.exp(lnQ)
deviation = lnQ - lnKeq  # = ln(Q/Keq), relaxes exponentially with time constant 1/k

# Plot: Q vs Keq (semilogy) and log deviation
plt.figure(figsize=(6.0, 4.0), dpi=150)
plt.semilogy(t, Keq, label="Keq(T)")
plt.semilogy(t, Q, label="Q(t)")
plt.xlabel("Time (s)"); plt.ylabel("Value"); plt.legend(); plt.tight_layout()
plt.savefig("t_jump.png")

# Optional: print the expected half-time from k
print("Half-time (s) ~ ln(2)/k =", np.log(2)/k)

