import numpy as np
import matplotlib.pyplot as plt

# Parameters
k = 1.0
K_eq = 2.0
C_total = 10.0
u = 3.0  # External drive

fig, axes = plt.subplots(1, 2, figsize=(10, 4))

# Panel A: [B] concentration over time
t = np.linspace(0, 5, 100)
alpha_values = [0, 0.5, 1.0, 2.0]
Q0 = 0.5  # Start below equilibrium

ax = axes[0]
for alpha in alpha_values:
    if alpha == 0:
        # Without feedback - use small alpha for numerical stability
        alpha_eff = 0 #0.01
        label = 'No feedback (α=0)'
    else:
        alpha_eff = alpha
        label = f'α={alpha}'
    
    # Analytical solution: Q(t) = K_eq * exp(u/(k+α)) * (Q0/K_eq * exp(-u/(k+α)))^exp(-(k+α)t)
    # Or more simply: Q(t) = Q_ss + (Q0 - Q_ss) * exp(-(k+α)t) in log space
    Q_ss = K_eq * np.exp(u / (k + alpha_eff))
    Q = Q_ss * (Q0 / Q_ss) ** np.exp(-(k + alpha_eff) * t)
    
    # Convert to [B] using conservation law
    B = Q * C_total / (1 + Q)
    
    ax.plot(t, B, label=label, linewidth=2)

ax.axhline(C_total, color='k', linestyle=':', alpha=0.5, label='C_total')
ax.set_xlabel('Time (s)')
ax.set_ylabel('[B] (mM)')
ax.set_title('(A) Feedback controls product accumulation')
ax.legend()
ax.grid(True, alpha=0.3)
ax.set_ylim([0, C_total*1.05])

# Panel B: Input-output relationship
ax = axes[1]
u_range = np.linspace(0, 5, 100)
for alpha in [0.5, 1.0, 2.0]:
    Q_ss = K_eq * np.exp(u_range / (k + alpha))
    ax.semilogy(u_range, Q_ss, label=f'α={alpha}', linewidth=2)

ax.set_xlabel('Drive strength u')
ax.set_ylabel('Steady-state Q')
ax.set_title('(B) Feedback reduces sensitivity to drive')
ax.legend()
ax.grid(True, alpha=0.3)
ax.set_ylim([1, 100])

plt.tight_layout()
plt.savefig('feedback_simple.png', dpi=300, bbox_inches='tight')
plt.show()
