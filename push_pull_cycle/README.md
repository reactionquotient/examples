# Push–Pull (Futile) Cycle — ATP/ADP-driven control

This example models a classic **push–pull phosphorylation cycle** (S ⇌ S~P) using the
**log-linear reaction quotient dynamics** from the paper “Log-Linear Reaction Quotient Dynamics.”
We treat the reaction quotient as the ratio \(Q \equiv r = [\mathrm{S\sim P}]/[\mathrm{S}]\) and include
the chemostatted ATP/ADP drive as an **additive control input** \(u(t)\) in log-space,

\[
\frac{d}{dt}\ln Q = -k\ln(Q/K_{eq}) + u(t), \quad
u(t) = \ln\!\big(\mathrm{ATP}/\mathrm{ADP}\big) + \beta,
\]

where \(\beta\) lumps standard-state and fixed-chemostat offsets (e.g., \( \ln r_0\)).
At steady state \( \ln(Q^*/K_{eq}) = u/k \), so for \(K_{eq}=1\) and \(k=1\)
we get \( r^* = r_0 \cdot (\mathrm{ATP}/\mathrm{ADP}) \) and the phosphorylated fraction
\( f^* = r^*/(1+r^*) \), i.e. a **logistic curve** in \(\ln(\mathrm{ATP}/\mathrm{ADP})\).

## What it shows
- **Steady-state switching:** \(f^*\) vs ATP/ADP is sigmoidal (logistic in log-space),
  illustrating energy-driven switching.
- **Time response:** step change in ATP/ADP produces a clean exponential trajectory in
  \(\ln Q\), hence a smooth transition in the phosphorylated fraction.

## Run
```bash
pip install numpy scipy matplotlib
python push_pull.py

