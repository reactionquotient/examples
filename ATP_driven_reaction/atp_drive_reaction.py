import numpy as np
import matplotlib.pyplot as plt

class ATPDrivenReaction:
    """
    Hexokinase reaction: Glucose + ATP → G6P + ADP
    Modeled as Q = [G6P]/[Glucose] with ATP/ADP ratio as external drive
    """
    
    def __init__(self):
        # Reaction parameters
        self.K_eq = 0.5  # Equilibrium constant favors glucose
        self.k = 1.0  # Relaxation rate (1/s)
        
        # ATP coupling
        self.k_ATP = 2.0  # Coupling strength to ATP/ADP ratio
        self.ATP_ADP_ratio = 10.0  # Typical cellular ratio
        
        # Conservation
        self.C_total = 10.0  # Total glucose + G6P (mM)
        
    def calculate_drive(self, ATP_ADP_ratio=None):
        """Calculate driving force from ATP/ADP ratio"""
        if ATP_ADP_ratio is None:
            ATP_ADP_ratio = self.ATP_ADP_ratio
        return self.k_ATP * np.log(ATP_ADP_ratio)
    
    def steady_state_Q(self, ATP_ADP_ratio=None):
        """Calculate steady-state reaction quotient"""
        u = self.calculate_drive(ATP_ADP_ratio)
        return self.K_eq * np.exp(u / self.k)
    
    def dynamics(self, t, Q0, ATP_ADP_ratio=None):
        """
        Analytical solution for Q(t)
        d ln Q/dt = -k ln(Q/K_eq) + u
        """
        u = self.calculate_drive(ATP_ADP_ratio)
        Q_ss = self.steady_state_Q(ATP_ADP_ratio)
        
        # Analytical solution
        Q = Q_ss * (Q0 / Q_ss) ** np.exp(-self.k * t)
        return Q
    
    def get_concentrations(self, Q):
        """Convert Q to actual concentrations using conservation"""
        glucose = self.C_total / (1 + Q)
        G6P = Q * self.C_total / (1 + Q)
        return glucose, G6P
    
    def plot_analysis(self):
        """Create focused two-panel figure"""
        t = np.linspace(0, 5, 100)
        Q0 = 10.0  # Start with some glucose already phosphorylated
        
        fig, axes = plt.subplots(1, 2, figsize=(10, 4))
        
        # Panel A: Q dynamics for different ATP/ADP ratios
        ax = axes[0]
        ATP_ADP_ratios = [0.1, 1.0, 10.0, 100.0]
        colors = ['red', 'orange', 'green', 'blue']
        
        for ratio, color in zip(ATP_ADP_ratios, colors):
            Q = self.dynamics(t, Q0, ratio)
            ax.semilogy(t, Q, color=color, linewidth=2, 
                       label=f'ATP/ADP = {ratio}')
        
        ax.axhline(self.K_eq, color='black', linestyle=':', 
                  label='K_eq (no ATP)', alpha=0.5)
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Q = [G6P]/[Glucose]')
        ax.set_title('(A) ATP drives reaction quotient')
        ax.legend(loc='best')
        ax.grid(True, alpha=0.3)
        ax.set_ylim([1e-4, 1e4])
        
        # Panel B: Glucose trapping efficiency vs ATP/ADP
        ax = axes[1]
        ratio_range = np.logspace(-1, 2, 100)
        Q_ss = [self.steady_state_Q(r) for r in ratio_range]
        
        # Convert to fraction phosphorylated
        fraction_G6P = np.array(Q_ss) / (1 + np.array(Q_ss))
        
        ax.semilogx(ratio_range, fraction_G6P, 'b-', linewidth=2.5)
        ax.axvline(self.ATP_ADP_ratio, color='red', linestyle='--', 
                  label='Cellular (ATP/ADP=10)', alpha=0.7, linewidth=2)
        ax.axhline(0.5, color='gray', linestyle=':', alpha=0.5)
        ax.set_xlabel('ATP/ADP ratio')
        ax.set_ylabel('Fraction as G6P')
        ax.set_title('(B) Glucose trapping efficiency')
        ax.grid(True, alpha=0.3)
        ax.legend(loc='lower right')
        ax.set_ylim([0, 1.05])
        
        #plt.suptitle('ATP-Driven Glucose Phosphorylation', 
        #            fontsize=13, fontweight='bold')
        plt.tight_layout()
        
        return fig
    
    def print_analysis(self):
        """Print key insights"""
        Q_no_ATP = self.K_eq
        Q_with_ATP = self.steady_state_Q()
        
        glucose_no_ATP, G6P_no_ATP = self.get_concentrations(Q_no_ATP)
        glucose_with_ATP, G6P_with_ATP = self.get_concentrations(Q_with_ATP)
        
        print("=== ATP-Driven Hexokinase ===\n")
        
        print("Without ATP drive:")
        print(f"  Q = {Q_no_ATP:.1f} (equilibrium favors glucose)")
        print(f"  Only {G6P_no_ATP/self.C_total:.1%} phosphorylated\n")
        
        print(f"With cellular ATP/ADP = {self.ATP_ADP_ratio}:")
        print(f"  Q = {Q_with_ATP:.1f} (driven far from equilibrium)")
        print(f"  {G6P_with_ATP/self.C_total:.1%} phosphorylated (glucose trapped!)\n")
        
        print("Key insight:")
        print(f"  ATP drives {(Q_with_ATP/Q_no_ATP):.0f}-fold increase in Q")
        print(f"  This exemplifies how cells use energy to control metabolism")

# Run the example
if __name__ == "__main__":
    model = ATPDrivenReaction()
    model.print_analysis()
    fig = model.plot_analysis()
    plt.savefig('hexokinase_ATP_drive.png', dpi=300, bbox_inches='tight')
    plt.show()
