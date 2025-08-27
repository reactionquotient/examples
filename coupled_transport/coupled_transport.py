import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import expm

class CoupledTransport:
    """
    Two coupled membrane transporters with clean eigenvalues
    """
    
    def __init__(self):
        # Equilibrium constants
        
        # K matrix designed for eigenvalues λ = 0.5 and 2.0
        # K = V @ diag(λ) @ V^(-1)
        # Using symmetric matrix for real, clean eigenvalues
        
        # More realistic parameters
        self.K_eq = np.array([1.0, 2.0])  # Different equilibrium constants
        self.k1 = 1.0    
        self.k2 = 2.0    # H⁺ transport faster
        self.alpha = 0.5  # Moderate coupling

        self.K_matrix = np.array([[self.k1,     self.alpha],
				  [self.alpha,  self.k2]])
        # Check eigenvalues
        eigenvals = np.linalg.eigvals(self.K_matrix)
        print(f"K matrix eigenvalues: {sorted(eigenvals)}")

        
    def dynamics(self, t, Q0, membrane_potential=1.0):
        """Solve dynamics with membrane potential driving both transporters"""
        Q = np.zeros((len(t), 2))
        Q[0] = Q0
        
        # Different drives
        u = np.array([np.log(membrane_potential),  # Na⁺ has additional drive
        	      np.log(membrane_potential)])        # H⁺ only membrane potential
        
        # Steady state
        x_ss = np.linalg.solve(self.K_matrix, u)
        
        for i in range(1, len(t)):
            dt = t[i] - t[i-1]
            x = np.log(Q[i-1] / self.K_eq)
            exp_Kt = expm(-self.K_matrix * dt)
            x_new = exp_Kt @ x + (np.eye(2) - exp_Kt) @ x_ss
            Q[i] = self.K_eq * np.exp(x_new)
        
        return Q
    
    def plot_analysis(self):
        """Two panel figure showing coupling and eigenmodes"""
        t = np.linspace(0, 8, 200)
        
        fig, axes = plt.subplots(1, 2, figsize=(10, 4))
        
        # Panel A: Coupled vs uncoupled dynamics
        ax = axes[0]
        Q0 = np.array([0.2, 0.2])
        
        # With coupling
        Q_coupled = self.dynamics(t, Q0, membrane_potential=3.0)
        
        # Without coupling
        K_uncoupled = np.array([[self.k1, 0], [0, self.k2]])
        self.K_matrix = K_uncoupled
        Q_uncoupled = self.dynamics(t, Q0, membrane_potential=3.0)
        self.K_matrix = np.array([[self.k1, self.alpha],
                                  [self.alpha, self.k2]])  # Restore
        
        ax.plot(t, Q_coupled[:, 0], 'b-', linewidth=2.5, 
                label='Q₁ (Na⁺) coupled')
        ax.plot(t, Q_coupled[:, 1], 'r-', linewidth=2.5, 
                label='Q₂ (H⁺) coupled')
        ax.plot(t, Q_uncoupled[:, 0], 'b--', linewidth=2, 
                label='Q₁ uncoupled', alpha=0.6)
        ax.plot(t, Q_uncoupled[:, 1], 'r--', linewidth=2, 
                label='Q₂ uncoupled', alpha=0.6)
        
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Transport quotient')
        ax.set_title('(A) Coupling synchronizes transport')
        ax.legend()
        ax.grid(True, alpha=0.3)
        #ax.set_ylim([0, 3.5])
        
        # Panel B: Eigenmode decomposition
        ax = axes[1]
        
        # Get eigenvectors and eigenvalues
        eigenvals, eigenvecs = np.linalg.eig(self.K_matrix)
        idx = eigenvals.argsort()
        eigenvals = eigenvals[idx]
        eigenvecs = eigenvecs[:, idx]
        
        # Initial condition with both modes excited
        Q0 = np.array([0.2, 3.0])
        Q = self.dynamics(t, Q0, membrane_potential=2.0)
        
        # Transform to eigenmode coordinates
        x = np.log(Q / self.K_eq)
        z = np.zeros((len(t), 2))
        for i in range(len(t)):
            z[i] = eigenvecs.T @ x[i]
        
        ax.plot(t, z[:, 0], 'purple', linewidth=2.5, 
                label=f'Slow mode (λ = {eigenvals[0]:.1f})')
        ax.plot(t, z[:, 1], 'teal', linewidth=2.5, 
                label=f'Fast mode (λ = {eigenvals[1]:.1f})')
        
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Mode amplitude')
        ax.set_title('(B) Eigenmodes reveal time scales')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        return fig

# Run the example
if __name__ == "__main__":
    model = CoupledTransport()
    fig = model.plot_analysis()
    plt.savefig('coupled_transport.png', dpi=300, bbox_inches='tight')
    plt.show()
