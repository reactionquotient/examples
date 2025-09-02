# conservation_decoupling/conservation_decoupling.py
import numpy as np
import matplotlib.pyplot as plt

# --- log-linear RQ dynamics params (single reaction) ---
k = 0.8           # relaxation rate (1/s)
Keq = 0.6         # equilibrium reaction quotient
Q0 = 3.0          # initial reaction quotient (>0)

# time grid
t = np.linspace(0, 10, 801)

# analytic solution: Q(t) = Keq * exp( ln(Q0/Keq) * exp(-k t) )
Q = Keq * np.exp(np.log(Q0/Keq) * np.exp(-k * t))

# --- map Q(t) to concentrations for different conserved totals ---
totals = [1.0, 3.0]     # two pool sizes to illustrate decoupling
AB = []
for Ctot in totals:
    A = Ctot / (1.0 + Q)
    B = Ctot * Q / (1.0 + Q)
    AB.append((A, B))

# --- plots ---
plt.figure(figsize=(6,4), dpi=140)

# Q(t): identical across cases
plt.plot(t, Q, color='#2E8B57', linewidth=2, label="Q(t)")  # Sea Green
plt.axhline(Keq, color='#B22222', linestyle="--", linewidth=2, label="Keq")  # Fire Brick

plt.xlabel("time")
plt.ylabel("reaction quotient Q")
plt.title("Log-linear RQ dynamics (same Q for all conserved totals)")
plt.legend()
plt.tight_layout()
plt.savefig('reaction_quotient_dynamics.png', dpi=300, bbox_inches='tight')
plt.close()

# concentrations: differ with Ctot though Q(t) is the same
plt.figure(figsize=(6,4), dpi=140)

# Define distinct colors for different total concentrations
colors = ['#1f77b4', '#ff7f0e']  # Blue and Orange
line_styles = ['-', '--']  # Solid and dashed

for i, (Ctot, (A, B)) in enumerate(zip(totals, AB)):
    plt.plot(t, A, color=colors[i], linestyle='-', linewidth=2, 
             label=f"[A](t), Ctot={Ctot}")
    plt.plot(t, B, color=colors[i], linestyle=':', linewidth=2, 
             label=f"[B](t), Ctot={Ctot}")
plt.xlabel("time")
plt.ylabel("concentration")
plt.title("Concentrations depend on C_tot; Q(t) does not")
plt.legend()
plt.tight_layout()
plt.savefig('concentration_dynamics.png', dpi=300, bbox_inches='tight')
plt.close()

