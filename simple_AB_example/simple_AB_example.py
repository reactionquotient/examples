import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# Set style for better plots
try:
    plt.style.use('seaborn-v0_8-darkgrid')
except OSError:
    try:
        plt.style.use('seaborn-darkgrid')
    except OSError:
        plt.style.use('default')

class SimpleABReaction:
    """
    Simple A ⇌ B reaction comparing log-linear dynamics with mass action kinetics
    """
    
    def __init__(self, k_f=1.0, K_eq=2.0):
        """
        Initialize with forward rate k_f and equilibrium constant K_eq
        """
        self.k_f = k_f  # Forward rate constant
        self.K_eq = K_eq  # Equilibrium constant [B]/[A]
        self.k_r = k_f / K_eq  # Reverse rate from detailed balance
        
        # Log-linear relaxation rate (matched to mass action near equilibrium)
        # From the paper: k = k_r(1 + K_eq)
        self.k_log = self.k_r * (1 + self.K_eq)
    
    def mass_action_Q_dynamics(self, Q, t):
        """
        Mass action dynamics for Q = [B]/[A]
        From paper equation (2): dQ/dt = k_r(1 + Q)(K_eq - Q)
        """
        return self.k_r * (1 + Q) * (self.K_eq - Q)
    
    def log_linear_dynamics(self, ln_Q, t):
        """
        Log-linear dynamics: d(ln Q)/dt = -k ln(Q/K_eq)
        From paper equation (1)
        """
        return -self.k_log * (ln_Q - np.log(self.K_eq))
    
    def simulate(self, Q0_values, t_max=10):
        """
        Simulate both dynamics for different initial conditions
        """
        t = np.linspace(0, t_max, 500)
        results = []
        
        for Q0 in Q0_values:
            # Mass action simulation
            Q_mass = odeint(self.mass_action_Q_dynamics, Q0, t).flatten()
            
            # Log-linear simulation
            ln_Q0 = np.log(Q0)
            ln_Q = odeint(self.log_linear_dynamics, ln_Q0, t).flatten()
            Q_log = np.exp(ln_Q)
            
            # Analytical solution for log-linear
            Q_analytical = self.K_eq * (Q0/self.K_eq)**np.exp(-self.k_log * t)
            
            results.append({
                'Q0': Q0,
                't': t,
                'Q_mass': Q_mass,
                'Q_log': Q_log,
                'Q_analytical': Q_analytical
            })
        
        return results
    
    def plot_comparison(self):
        """
        Create a simple, clear comparison figure with just the two best panels
        """
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        
        # Test with different initial conditions
        Q0_values = [0.5, 1.0, 4.0, 8.0]  # Below, at, and above equilibrium
        colors = ['blue', 'green', 'orange', 'red']
        
        # Get results
        results = self.simulate(Q0_values, t_max=5)
        
        # Panel A: Q vs time for all initial conditions
        ax = axes[0]
        for res, color in zip(results, colors):
            ax.plot(res['t'], res['Q_mass'], '-', color=color, linewidth=2,
                   label=f"Q₀={res['Q0']:.1f} (MA)", alpha=0.8)
            ax.plot(res['t'], res['Q_log'], '--', color=color, linewidth=2,
                   label=f"Q₀={res['Q0']:.1f} (LL)", alpha=0.8)
        
        ax.axhline(y=self.K_eq, color='black', linestyle=':', linewidth=1.5, 
                  label=f'K_eq={self.K_eq}')
        ax.set_xlabel('Time (s)', fontsize=12)
        ax.set_ylabel('Q = [B]/[A]', fontsize=12)
        ax.set_title('(A) Reaction Quotient Evolution', fontsize=14)#, fontweight='bold')
        ax.legend(loc='best', fontsize=8, ncol=2)
        ax.grid(True, alpha=0.3)
        ax.set_xlim([0, 5])
        
        # Panel B: Near equilibrium comparison
        ax = axes[1]
        Q0_near = [1.8, 2.0, 2.2]  # Very close to K_eq = 2
        results_near = self.simulate(Q0_near, t_max=5)
        
        for res in results_near:
            ax.plot(res['t'], res['Q_mass'], '-', linewidth=2,
                   label=f"Q₀={res['Q0']:.1f} (Mass Action)")
            ax.plot(res['t'], res['Q_log'], '--', linewidth=2,
                   label=f"Q₀={res['Q0']:.1f} (Log-Linear)")
        
        ax.axhline(y=self.K_eq, color='black', linestyle=':', linewidth=1.5,
                  label=f'K_eq={self.K_eq}')
        ax.set_xlabel('Time (s)', fontsize=12)
        ax.set_ylabel('Q = [B]/[A]', fontsize=12)
        ax.set_title('(B) Near-Equilibrium Behavior', fontsize=14)#, fontweight='bold')
        ax.legend(loc='best', fontsize=9)
        ax.grid(True, alpha=0.3)
        ax.set_xlim([0, 5])
        ax.set_ylim([1.7, 2.3])
        
        #plt.suptitle('A ↔ B Reaction: Mass Action vs Log-Linear Dynamics', 
        #            fontsize=14, fontweight='bold', y=1.02)
        plt.tight_layout()
        
        return fig
    
    def print_summary(self):
        """Print key parameters"""
        print("=== A ⇌ B Reaction Parameters ===")
        print(f"Forward rate k_f: {self.k_f:.3f} /s")
        print(f"Reverse rate k_r: {self.k_r:.3f} /s")
        print(f"Equilibrium constant K_eq: {self.K_eq:.2f}")
        print(f"Log-linear rate k: {self.k_log:.3f} /s")
        print(f"Time constant τ = 1/k: {1/self.k_log:.3f} s")
        print(f"\nAt equilibrium: [B]/[A] = {self.K_eq:.2f}")

# Run the example
if __name__ == "__main__":
    # Create model with simple parameters
    model = SimpleABReaction(k_f=1.0, K_eq=2.0)
    
    # Print parameters
    model.print_summary()
    
    # Create comparison plots
    fig = model.plot_comparison()
    
    # Save figure
    fig.savefig('AB_reaction_comparison.png', dpi=300, bbox_inches='tight')
    print("\nFigure saved as 'AB_reaction_comparison.png'")
    
    plt.show()
