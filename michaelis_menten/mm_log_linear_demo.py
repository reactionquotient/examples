
"""
mm_log_linear_demo.py

Simple demonstration of log-linear reaction-quotient (RQ) dynamics for
Michaelis-Menten kinetics using analytical solutions.

Model: E + S <-> ES -> E + P (binding equilibrium + irreversible catalysis)

Key insight: In log-linear RQ dynamics, we evolve x = ln(Q/K) where:
- Q(t) = K * exp(-k*t + integral(u)) is the analytical solution
- When u=0 (no drive), Q(t) = K * exp(-k*t) relaxes to equilibrium
- At equilibrium (t→∞), [ES] = E_tot * K*S / (1 + K*S) gives MM kinetics
"""

import numpy as np
import matplotlib.pyplot as plt

def mm_rate(S, Vmax, Km):
    """Classic Michaelis-Menten rate law."""
    return Vmax * S / (Km + S)

def rq_analytical(t, Q0, K, k, u=0):
    """Analytical solution for log-linear RQ dynamics.
    
    From dx/dt = -k*x + u where x = ln(Q/K):
    Solution: x(t) = x₀*exp(-kt) + (u/k)(1-exp(-kt))
    Therefore: Q(t) = K * exp(x(t))
    """
    x0 = np.log(Q0 / K)
    if u == 0:
        # No drive: x(t) = x₀*exp(-kt)
        x_t = x0 * np.exp(-k * t)
    else:
        # With drive: x(t) = x₀*exp(-kt) + (u/k)(1-exp(-kt))
        x_t = x0 * np.exp(-k * t) + (u / k) * (1 - np.exp(-k * t))
    
    return K * np.exp(x_t)

def es_concentration(E_tot, S, Q):
    """Enzyme-substrate concentration from binding quotient."""
    return E_tot * Q * S / (1 + Q * S)

def demo_rq_dynamics():
    """Demonstrate log-linear RQ dynamics and its connection to MM kinetics."""
    # Parameters
    E_tot = 1.0    # enzyme total (μM) 
    K = 0.1        # binding equilibrium constant (1/μM)
    k = 10.0       # RQ relaxation rate (1/s) 
    k_cat = 1.0    # catalytic turnover (1/s)
    
    # Derived MM parameters (fast binding limit where Q₁ ≈ K)
    K_D = 1.0 / K  # dissociation constant (μM)
    V_max = k_cat * E_tot
    K_M = K_D      # in rapid equilibrium approximation
    
    # Substrate range for steady-state comparison
    S_vals = np.logspace(-1, 2, 100)  # 0.1 to 100 μM
    
    # Time points for relaxation dynamics  
    t = np.linspace(0, 0.5, 200)  # 0 to 0.5 seconds
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # LEFT: Steady-state v(S) comparison
    # MM formula
    v_mm = mm_rate(S_vals, V_max, K_M)
    
    # RQ at equilibrium (Q = K, no dynamics)
    ES_eq = es_concentration(E_tot, S_vals, K)
    v_rq_eq = k_cat * ES_eq
    
    ax1.loglog(S_vals, v_mm, 'b-', label='Michaelis-Menten', linewidth=2)
    ax1.loglog(S_vals, v_rq_eq, 'r:', label='RQ (fast binding limit)', linewidth=2)
    ax1.set_xlabel('[S] (μM)')
    ax1.set_ylabel('Rate v (μM/s)')
    ax1.set_title('Steady-State: MM vs RQ Framework')
    ax1.legend(loc="lower right")
    ax1.grid(True, alpha=0.3)
    
    # Add annotation about the equivalence
    ax1.text(0.02, 0.95, 'Identical when Q₁ ≈ K (fast binding)', 
             transform=ax1.transAxes, fontsize=9, 
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    # RIGHT: RQ dynamics with different control inputs
    S_fixed = 10.0  # fixed substrate (μM)
    Q_0 = 3.0 * K   # start 3x above equilibrium
    
    # Three cases: no drive, positive drive, negative drive
    u_values = [0.0, 2.0, -2.0]
    colors = ['green', 'red', 'blue']
    labels = ['u = 0 (thermal equilibrium)', 'u = +2 (energy input)', 'u = -2 (energy dissipation)']
    
    # Equilibrium reference (MM prediction when u=0)
    v_eq_mm = mm_rate(S_fixed, V_max, K_M)
    
    for u, color, label in zip(u_values, colors, labels):
        Q_t = rq_analytical(t, Q_0, K, k, u=u)
        ES_t = es_concentration(E_tot, S_fixed, Q_t)
        v_t = k_cat * ES_t
        ax2.plot(t, v_t, color=color, linewidth=2, label=label)
    
    ax2.axhline(y=v_eq_mm, color='black', linestyle='--', linewidth=1, alpha=0.7,
                label='MM equilibrium (u=0)')
    ax2.set_xlabel('Time (s)')
    ax2.set_ylabel('Rate v (μM/s)')
    ax2.set_title(f'RQ Dynamics with Control Input u\n[S] = {S_fixed} μM')
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.3)
    
    # Add annotation about control input
    ax2.text(0.02, 0.98, 'Control input u shifts steady-state\naway from thermal equilibrium', 
             transform=ax2.transAxes, fontsize=9, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig('mm_rq_comparison.png', dpi=150, bbox_inches='tight')
    
    # Print key insights
    print(f"Parameters:")
    print(f"  E_tot = {E_tot} μM")
    print(f"  K = {K} μM⁻¹ (binding equilibrium constant)")  
    print(f"  K_D = 1/K = {K_D} μM (dissociation constant)")
    print(f"  k = {k} s⁻¹ (RQ relaxation rate)")
    print(f"  k_cat = {k_cat} s⁻¹ (catalytic rate)")
    print(f"\nDerived MM parameters (fast binding limit):")
    print(f"  V_max = k_cat × E_tot = {V_max} μM/s")
    print(f"  K_M ≈ K_D = {K_M} μM")
    print(f"\nRelaxation timescale: τ = 1/k = {1/k:.1f} s")


if __name__ == "__main__":
    demo_rq_dynamics()
    print("\nSaved: mm_rq_comparison.png")
    print("\nKey insights demonstrated:")
    print("• Left panel: MM and RQ frameworks are identical when Q₁ ≈ K (fast binding)")
    print("• Right panel: Control input u drives the system away from thermal equilibrium") 
    print("• u = 0: Relaxes to MM equilibrium (green curve)")
    print("• u > 0: Energy input pushes reaction quotient above equilibrium (red)")
    print("• u < 0: Energy dissipation pulls reaction quotient below equilibrium (blue)")
