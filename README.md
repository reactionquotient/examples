# Log-Linear Reaction Quotient Dynamics: Implementation Examples

This repository contains Python implementations of the examples presented in the paper **"Log-Linear Reaction Quotient Dynamics"** by Steven Diamond ([arXiv:2508.18523](https://arxiv.org/pdf/2508.18523)).

## 📋 Overview

This framework introduces a novel approach to modeling chemical reaction networks where reaction quotients (Q) evolve exponentially toward equilibrium when viewed on a logarithmic scale. Unlike traditional mass action kinetics, this yields analytically tractable linear dynamics in log-space.

### Key Equation

For a single reaction, the dynamics follow:
```
d ln Q/dt = -k ln(Q/Keq)
```

Where:
- **Q**: Reaction quotient measuring distance from equilibrium
- **Keq**: Equilibrium constant
- **k**: Relaxation rate constant

### Framework Advantages

1. **Analytical Solutions**: Exact solutions exist for arbitrary network topologies
2. **Thermodynamic Integration**: Automatic incorporation of constraints via ΔG = RT ln(Q/Keq)
3. **Decoupled Dynamics**: Conservation laws separate from reaction quotient evolution
4. **Linear Control**: External energy sources (e.g., ATP) couple linearly to dynamics
5. **Tractable Analysis**: Decades of linear systems theory become applicable

## 🧪 Examples

### 1. Simple A ⇌ B Reaction (`simple_AB_example/`)
**File**: `simple_AB_example.py`  
**Figure**: `AB_reaction_comparison.png`  
**Section**: 4 (Mass action comparison)

Validates the log-linear framework against traditional mass action kinetics for the simplest reversible reaction. The parameter matching k = kr(1 + Keq) ensures agreement near equilibrium, while far from equilibrium the log-linear model shows faster relaxation due to exponential decay in log-space.

**Key Insights**:
- Near equilibrium: virtually indistinguishable from mass action
- Far from equilibrium: faster convergence due to log-space dynamics
- Analytical solution: Q(t) = Keq(Q0/Keq)^exp(-kt)

### 2. Feedback Inhibition (`feedback_inhibition/`)
**File**: `feedback_simple.py`  
**Figure**: `feedback_simple.png`  
**Section**: 4 (Feedback inhibition)

Demonstrates the minimal feedback motif where product B inhibits its own formation from A. Feedback effectively increases the relaxation rate from k to k + α, implementing a natural robustness-responsiveness tradeoff.

**Key Features**:
- Steady state: Qss = Keq exp(u/(k + α))
- Feedback prevents unlimited accumulation despite constant drive
- Reduced sensitivity to drive fluctuations with stronger feedback
- System remains stable for all α > 0

### 3. ATP-Driven Hexokinase Reaction (`ATP_driven_reaction/`)
**File**: `atp_drive_reaction.py`  
**Figure**: `hexokinase_ATP_drive.png`  
**Section**: 4 (ATP-driven reaction)

Models the first step of glycolysis where glucose is phosphorylated using ATP energy. Demonstrates how cellular energy state can drive reactions far from chemical equilibrium through thermodynamic coupling.

**Highlights**:
- Control input: u = kATP ln([ATP]/[ADP])
- Sharp sigmoid response between ATP/ADP ratios of 1-10
- At cellular conditions: 98% glucose trapped as G6P despite unfavorable chemistry
- Biological switching behavior emerges from thermodynamic coupling

### 4. Coupled Membrane Transport (`coupled_transport/`)
**File**: `coupled_transport.py`  
**Figure**: `coupled_transport.png`  
**Section**: 4 (Coupled transport)

Illustrates two membrane transporters (Na+ and H+) sharing membrane potential through off-diagonal coupling in the rate matrix K. Demonstrates emergence of complex non-monotonic dynamics.

**Phenomena**:
- Overshoot behavior in H+ transport
- Eigenmode decomposition explains dynamics
- Fast mode (λ = 2.2): transporters oppose
- Slow mode (λ = 0.8): transporters move together
- Interference between modes creates transient overshoot

### 5. Glycolytic Oscillations (`glycolitic_oscillations/`)
**File**: `glycolysis.py`  
**Figure**: `glycolytic_oscillations.png`  
**Section**: 4 (Glycolytic oscillations)

Models oscillatory behavior in a simplified glycolytic pathway with FBP feedback activation. Shows how metabolic rhythms arise from linear dynamics in log-space without requiring Hill functions.

**Analysis**:
- Eigenvalues λ = 0.5 ± 2i determine oscillation properties
- Period: τ = 2π/ω ≈ 3.1 seconds (matches yeast experiments)
- Damping without ATP drive
- Sustained limit cycle with resonant ATP forcing
- Analytical solution available for transient and steady-state behavior

### 6. Bimolecular Exchange Reaction (`bimolecular_exchange/`)
**File**: `ab_cd_example.py`  
**Figures**: `ab_cd_reaction.png`, `ab_cd_Q_relaxation.png`  
**Section**: Extension of framework

Models the A + B ⇌ C + D reaction under log-linear dynamics. Demonstrates quadratic root-finding to recover species concentrations from reaction quotient Q(t), with physically valid solutions ensuring all concentrations remain nonnegative.

**Key Features**:
- Closed-form Q(t) = Keq(Q0/Keq)^exp(-kt) solution
- Quadratic solver for extent of reaction ξ(t) from Q(t)
- Continuity enforcement for smooth concentration trajectories
- Mass-action comparison

### 7. Push–Pull Cycle (`push_pull_cycle/`)

File: `push_pull.py`  
Figures: `push_pull_fraction_vs_ATP.png`, `push_pull_step_response.png`  
Section: 4 (Linear control / energy-driven switching)

Shows a thermodynamic push–pull (futile) phosphorylation cycle driven by ATP/ADP.
In log-space the dynamics are linear: $\frac{d}{dt}\ln Q = -k\ln(Q/K_{eq}) + u$.
The steady-state phosphorylated fraction is a logistic function of $\ln(\mathrm{ATP}/\mathrm{ADP})$,
and step changes in drive yield exponential responses in $\ln Q$.

### 8. General Single Reaction (`single_reaction_general/`)
**File**: `single_reaction_general.py`  
**Figures**: `general_single_reaction_conc.png`, `general_single_reaction_Q.png`  
**Section**: Framework extension

Demonstrates the log-linear framework for arbitrary single reactions with general stoichiometry. Uses the example reaction 2A + B ⇌ C + 2D, showing how to solve for species concentrations from reaction quotient Q(t) via root-finding on the extent of reaction ξ.

**Key Features**:
- General stoichiometric coefficients: α (reactants), β (products)
- Root-finding solution: S(ξ) = ln Q(t) with feasibility constraints
- Robust numerical methods with boundary clamping
- Direct comparison to mass-action kinetics
- Demonstrates crossing dynamics as system evolves toward equilibrium

## 🔧 Installation & Usage

### Requirements
```bash
pip install numpy scipy matplotlib
```

### Running Examples
Each example can be run independently:
```python
python simple_AB_example/simple_AB_example.py
python feedback_inhibition/feedback_simple.py
python ATP_driven_reaction/atp_drive_reaction.py
python coupled_transport/coupled_transport.py
python glycolitic_oscillations/glycolysis.py
python bimolecular_exchange/ab_cd_example.py
python single_reaction_general/single_reaction_general.py
```

### Parameter Exploration
All examples include adjustable parameters at the top of each file. Modify these to explore different dynamical regimes:
- Relaxation rates (k)
- Equilibrium constants (Keq)
- Coupling strengths (α, kATP)
- Initial conditions (Q0)

## 📚 Mathematical Framework

### Single Reaction Dynamics
Solution: `Q(t) = Keq * (Q0/Keq)^exp(-kt)`

With control input u:
`Q(t) = Keq * exp[ln(Q0/Keq - u/k) * exp(-kt) + u/k]`

### Multiple Reactions
Vector form: `d/dt ln Q = -K ln(Q/Keq) + u`

Where K is the relaxation rate matrix determining:
- Diagonal elements: individual reaction rates
- Off-diagonal elements: reaction coupling
- Eigenvalues: stability and oscillation properties

### Thermodynamic Interpretation
The framework naturally incorporates free energy:
- ΔG = RT ln(Q/Keq)
- Dynamics become: dΔG/dt = -kΔG
- Control inputs represent external energy gradients

## 🎯 Applications

This framework may enable:
- **Metabolic Engineering**: Optimize pathway design using K as design variable
- **Drug Discovery**: Predict drug effects throughout metabolic networks
- **Systems Medicine**: Classify metabolic disorders via eigenvalue analysis
- **Control Theory**: Apply optimal control to cellular metabolism

## 📖 Citation

If you use this code in your research, please cite:

```bibtex
@article{diamond2025loglinear,
  title={Log-Linear Reaction Quotient Dynamics},
  author={Diamond, Steven},
  journal={arXiv preprint arXiv:2508.18523},
  year={2025}
}
```

## 📄 License

This code is licensed under the Apache License 2.0. See the `LICENSE` file for details.

## 🔗 Resources

- **Paper**: [arXiv:2508.18523](https://arxiv.org/pdf/2508.18523)
- **Repository**: [github.com/reactionquotient/examples](https://github.com/reactionquotient/examples)
- **Contact**: steven@gridmatic.com

