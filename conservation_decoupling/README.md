# Conservation–Decoupling (Log‑Linear Reaction Quotient)

This tiny example demonstrates a core property of **log‑linear reaction quotient (RQ) dynamics** on the simplest network \(A \rightleftharpoons B\):  
**the trajectory of the reaction quotient \(Q(t)\)** (and thus its relaxation to equilibrium) **is independent of the conserved pool size** \(C_{\text{tot}}=[A]+[B]\). Concentrations are reconstructed *after* the RQ dynamics via conservation, so different totals produce different concentration time-courses while sharing the **same** \(Q(t)\).

---

## What this shows (and why it’s distinct)

- **State variable:** \(Q = \dfrac{[B]}{[A]}\).  
- **Dynamics (log‑linear form):**  
  \[
  \frac{d}{dt}\ln Q \;=\; -k\,\ln\!\Big(\frac{Q}{K_{\mathrm{eq}}}\Big)
  \]
- **Closed‑form solution:**  
  \[
  Q(t)\;=\;K_{\mathrm{eq}}\exp\!\Big(\ln\!\frac{Q_0}{K_{\mathrm{eq}}}\;e^{-kt}\Big).
  \]
- **Concentration reconstruction (via conservation for \(A\rightleftharpoons B\))** with \(C_{\text{tot}}=[A]+[B]\):  
  \[
  [A](t)=\frac{C_{\text{tot}}}{1+Q(t)} ,\qquad
  [B](t)=\frac{C_{\text{tot}}\,Q(t)}{1+Q(t)} .
  \]
- **Key takeaway:** Changing \(C_{\text{tot}}\) alters \([A](t),[B](t)\) but **does not change** \(Q(t)\).  
- **Non‑redundant focus:** No temperature jumps, no feedback, no transport, no multi‑reaction coupling—just the RQ law plus conservation mapping.

---

## Files

```
conservation_decoupling/
├─ conservation_decoupling.py   # script to generate the plots
└─ README.md                    # this file
```

---

## How to run

```bash
# from the repository root (or this folder), with Python 3.9+
pip install numpy matplotlib

python conservation_decoupling.py
```

The script will produce two figures:

1. **RQ dynamics:** \(Q(t)\) relaxing toward \(K_{\mathrm{eq}}\) (identical for all chosen \(C_{\text{tot}}\)).  
2. **Concentrations:** \([A](t)\) and \([B](t)\) for two different totals (e.g., \(C_{\text{tot}}=1\) vs \(3\)), illustrating that concentrations differ even though \(Q(t)\) is the same.

---

## Parameters (edit in the script)

- `k` — relaxation rate (s\(^{-1}\)) in the log‑linear law.  
- `Keq` — equilibrium reaction quotient.  
- `Q0` — initial reaction quotient.  
- `totals` — list of pool sizes (e.g., `[1.0, 3.0]`) used to reconstruct concentrations from the same \(Q(t)\).

---

## Minimal theory recap

The **log‑linear RQ dynamics** specify the evolution of \(Q\) in **log space**, not in concentration space. After solving for \(Q(t)\), any species levels are computed via algebraic constraints (here, a single conservation law for \(A\) and \(B\)). This cleanly separates thermodynamic/kinetic relaxation (first plot) from bookkeeping of conserved matter (second plot).

---

## Optional tiny extension (still single‑reaction)

Add a constant drive \(u\) to the log‑linear ODE,
\(\frac{d}{dt}\ln Q = -k\,\ln(Q/K_{\mathrm{eq}})+u\),
to create a non‑equilibrium offset with steady value \(Q_\infty = K_{\mathrm{eq}}\,e^{u/k}\). Keep this extension in the same folder to preserve the **single‑reaction, no‑feedback** focus.
