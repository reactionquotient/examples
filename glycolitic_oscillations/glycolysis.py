import numpy as np
import matplotlib.pyplot as plt

class GlycolyticOscillations:
    """
    Analytical solution for glycolytic oscillations using eigenvalue analysis.
    
    Q1 = [FBP]/[F6P] (PFK reaction)
    Q2 = [Products]/[FBP] (downstream reactions)
    """
    
    def __init__(self):
        # Parameters from K matrix eigenvalues
        self.damping = 0.5  # Damping rate (real part)
        self.frequency = 2.0  # Oscillation frequency (imaginary part)
        self.period = 2 * np.pi / self.frequency
        self.Keq = np.array([2.0, 1.0])  # Equilibrium constants
        
        # Initial condition coefficients
        self.A1, self.B1 = 0.5, 0.2  # Coefficients for Q1
        self.A2, self.B2 = -0.3, 0.4  # Coefficients for Q2
        
        # Steady-state amplitude for sustained oscillations
        self.amplitude = 1.6  # Amplitude in log-space
        
        print(f"Eigenvalues: {self.damping} ± {self.frequency}i")
        print(f"Oscillation period: {self.period:.2f} seconds")
        
    def solve_dynamics(self, t_span=20):
        """Generate time series for both damped and sustained oscillations."""
        t = np.linspace(0, t_span, 1000)
        
        # Damped oscillations (no ATP drive)
        envelope = np.exp(-self.damping * t)
        x1_damped = envelope * (self.A1 * np.cos(self.frequency * t) + 
                               self.B1 * np.sin(self.frequency * t))
        Q1_damped = self.Keq[0] * np.exp(x1_damped)
        
        # Sustained oscillations (with ATP drive)
        x1_sustained = self.amplitude * np.sin(self.frequency * t + 0.3)
        Q1_sustained = self.Keq[0] * np.exp(x1_sustained)
        
        # Phase space coordinates
        x2_damped = envelope * (self.A2 * np.cos(self.frequency * t) + 
                               self.B2 * np.sin(self.frequency * t))
        x2_sustained = self.amplitude * np.cos(self.frequency * t + 0.3)
        
        return t, Q1_damped, Q1_sustained, x1_damped, x1_sustained, x2_damped, x2_sustained
    
    def create_figure(self):
        """Create two-panel figure showing oscillations and phase portrait."""
        t, Q1_damped, Q1_sustained, x1_damped, x1_sustained, x2_damped, x2_sustained = self.solve_dynamics()
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

        # Panel A: Time series
        ax1.plot(t, Q1_damped, 'b--', linewidth=2, alpha=0.6, 
                label='No ATP drive')
        ax1.plot(t, Q1_sustained, 'b-', linewidth=2.5, 
                label='With ATP drive')
        
        # Add equilibrium line
        ax1.axhline(y=self.Keq[0], color='b', linestyle='--', alpha=0.3, 
                   label='K_eq')
        
        # Panel formatting
        ax1.set_xlabel('Time (s)', fontsize=12)
        ax1.set_ylabel('Reaction Quotient Q₁', fontsize=12)
        ax1.set_title('(A) Glycolytic Oscillations', fontsize=14)
        ax1.legend(loc='upper right')
        ax1.grid(True, alpha=0.3)
        ax1.set_xlim([0, 20])
        
        # Add period annotation - position above the sustained oscillation
        t_mark = 15
        # Find max value of sustained oscillation in this time range
        t_indices = (t >= t_mark) & (t <= t_mark + self.period)
        y_max = np.max(Q1_sustained[t_indices])
        y_pos = y_max + 0.5
        
        #ax1.annotate('', xy=(t_mark + self.period, y_pos), xytext=(t_mark, y_pos),
        #            arrowprops=dict(arrowstyle='<->', color='gray', lw=1.5))
        #ax1.text(t_mark + self.period/2, y_pos + 0.2, f'τ = {self.period:.1f} s', 
        #        ha='center', color='gray', fontsize=10)
        
        # Panel B: Phase portrait)
        # Damped spiral (first part only for clarity)
        t_spiral = t[t < 8]
        x1_spiral = np.exp(-self.damping * t_spiral) * (
            self.A1 * np.cos(self.frequency * t_spiral) + 
            self.B1 * np.sin(self.frequency * t_spiral))
        x2_spiral = np.exp(-self.damping * t_spiral) * (
            self.A2 * np.cos(self.frequency * t_spiral) + 
            self.B2 * np.sin(self.frequency * t_spiral))
        
        ax2.plot(x1_spiral, x2_spiral, 'gray', linewidth=1.5, alpha=0.5, 
                label='Damped (no drive)')
        
        # Sustained limit cycle
        theta = np.linspace(0, 2*np.pi, 200)
        x1_cycle = self.amplitude * np.sin(theta + 0.3)
        x2_cycle = self.amplitude * np.cos(theta + 0.3)
        ax2.plot(x1_cycle, x2_cycle, 'b-', linewidth=2.5, 
                label='Limit cycle (ATP drive)')
        
        # Mark key points
        ax2.plot(0, 0, 'r*', markersize=15, label='Equilibrium', zorder=5)
        ax2.plot(x1_spiral[0], x2_spiral[0], 'go', markersize=10, 
                label='Initial state', zorder=5)
        
        # Add direction arrow on limit cycle
        idx = 50
        dx = x1_cycle[idx+1] - x1_cycle[idx]
        dy = x2_cycle[idx+1] - x2_cycle[idx]
        ax2.arrow(x1_cycle[idx], x2_cycle[idx], 5*dx, 5*dy, 
                 head_width=0.1, head_length=0.08, fc='b', ec='b')
        
        # Panel formatting
        ax2.set_xlabel('ln(Q₁/K₁)', fontsize=12)
        ax2.set_ylabel('ln(Q₂/K₂)', fontsize=12)
        ax2.set_title('(B) Phase Portrait in Log-Space', fontsize=14)
        ax2.legend(loc='upper left')
        ax2.grid(True, alpha=0.3)
        ax2.set_aspect('equal')
        ax2.set_xlim([-2, 2])
        ax2.set_ylim([-2, 2])
        
        # Add eigenvalue annotation
        #ax2.text(0.02, 0.98, f'K eigenvalues: {self.damping:.1f} ± {self.frequency:.0f}i', 
        #        transform=ax2.transAxes, fontsize=10,
        #        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7),
        #        verticalalignment='top')
        
        #plt.suptitle('Log-Linear Framework: Analytical Solution', 
        #            fontsize=16, fontweight='bold')
        plt.tight_layout()
        return fig


if __name__ == "__main__":
    print("="*60)
    print("GLYCOLYTIC OSCILLATIONS - ANALYTICAL")
    print("="*60)
    
    oscillator = GlycolyticOscillations()
    fig = oscillator.create_figure()
    
    print("\nKey features:")
    print("- Analytical solution from eigenvalue analysis")
    print("- Damped oscillations decay without ATP drive")
    print("- Sustained oscillations maintained by ATP forcing")
    
    #plt.savefig('glycolytic_oscillations.pdf', dpi=300, bbox_inches='tight')
    plt.savefig('glycolytic_oscillations.png', dpi=150, bbox_inches='tight')
    plt.show()
