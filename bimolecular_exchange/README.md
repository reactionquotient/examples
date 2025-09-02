# A + B ⇌ C + D (Bimolecular Exchange)

This example extends the repo’s A ⇌ B demo to a 2→2 reaction. We evolve the reaction quotient

> Q(t) = Keq · (Q0/Keq)^{exp(-k t)}

which is the closed-form solution of the log-linear dynamics

> d/dt ln Q = -k · ln(Q/Keq)

and then reconstruct species in a closed system via the extent of reaction ξ(t):

- A(t) = A0 − ξ, B(t) = B0 − ξ  
- C(t) = C0 + ξ, D(t) = D0 + ξ  

with ξ found at each time by solving a quadratic implied by Q(t) = (C·D)/(A·B).

Run:

```bash
python bimolecular_exchange/ab_cd_example.py

