# Examples: Log-Linear Dynamics and the Reaction Quotient

This directory contains Python implementations of the examples presented in the paper **"Log-Linear Reaction Quotient Dynamics"** ([arXiv:2508.18523](https://arxiv.org/pdf/2508.18523)).

## Overview

These examples demonstrate log-linear dynamics for biochemical reaction networks, where the dynamics are expressed in terms of reaction quotients Q rather than individual concentrations. The log-linear approach provides:

- Analytical solutions
- Natural incorporation of thermodynamic constraints
- Direct connection to free energy landscapes

## Examples

### 1. Simple A ⇌ B Reaction (`simple_AB_example/`)
**File**: `simple_AB_example.py`  
**Figure**: `AB_reaction_comparison.png`

Compares log-linear dynamics with traditional mass action kinetics for a simple reversible reaction A ⇌ B. Shows how log-linear dynamics provide accurate approximations, particularly near equilibrium.

### 2. Feedback Inhibition (`feedback_inhibition/`)
**File**: `feedback_simple.py`  
**Figure**: `feedback_simple.png`

Shows the stabilizing effect of feedback on reaction rates and steady states within the log-linear framework..

### 3. ATP-Driven Reaction (`ATP_driven_reaction/`)
**File**: `atp_drive_reaction.py`  
**Figure**: `hexokinase_ATP_drive.png`

Models the hexokinase reaction (Glucose + ATP → G6P + ADP) as an example of how ATP hydrolysis can drive reactions away from equilibrium. Demonstrates how control inputs factor into the log-linear framework.

### 4. Coupled Transport (`coupled_transport/`)
**File**: `coupled_transport.py`  
**Figure**: `coupled_transport.png`

Illustrates coupled transport processes where the flow of one species drives the transport of another against its concentration gradient. Shows how overshoots and other non-monotonic behavior can emerge from the log-linear framework.

### 5. Glycolytic Oscillations (`glycolitic_oscillations/`)
**File**: `glycolysis.py`  
**Figure**: `glycolytic_oscillations.png`

Models oscillatory behavior in glycolysis using eigenvalue analysis of the log-linear system. Compares damped oscillations (no ATP drive) with sustained oscillations (with ATP drive).

## License

This code is licensed under the Apache License 2.0. See the `LICENSE` file for details.
