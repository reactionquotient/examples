# What is *u*? (Mass-Action vs. Log-Linear)

## TL;DR
**\(u\)** is a **dimensionless driving force**—typically a log ratio of chemostatted species (e.g., \(u=\ln([\mathrm{ATP}]/[\mathrm{ADP}])\)) or, equivalently, a chemical-potential difference divided by \(RT\). In **mass‑action**, it acts by **rescaling effective rate constants** when you coarse‑grain out the cofactor. In the **log‑linear** quotient model, it appears **additively** as the natural thermodynamic drive.

---

## Definition
- **Dimensionless:** \(u \equiv \Delta\mu/RT\), or concretely a log‑ratio like \(u=\ln\!\frac{[\mathrm{ATP}]}{[\mathrm{ADP}]}\).
- **Physical meaning:** a free‑energy bias supplied by an external bath (a chemostat).

---

## How *u* fits **mass‑action**

Start from a cofactor‑using step:
\[
A + \mathrm{ATP} \;\rightleftharpoons\; B + \mathrm{ADP}
\]
with mass‑action rates
\[
v_f = k_f [A][\mathrm{ATP}],\qquad
v_r = k_r [B][\mathrm{ADP}].
\]

If \([\mathrm{ATP}]\) and \([\mathrm{ADP}]\) are **chemostatted** (held by an external reservoir), define
\[
u \equiv \ln\!\frac{[\mathrm{ATP}]}{[\mathrm{ADP}]}\,.
\]
When we **coarse‑grain** this into an effective two‑state step \(A \rightleftharpoons B\), the cofactor concentrations fold into the *effective* rate constants:
\[
k_f^{\mathrm{eff}} = k_f\,[\mathrm{ATP}] = k_f^{0}\,e^{u},\qquad
k_r^{\mathrm{eff}} = k_r\,[\mathrm{ADP}] = k_r^{0}\,e^{-u}.
\]

So in mass‑action it is legitimate to model a time‑varying drive by
\[
k_f(t) = k_f^{0} e^{u(t)},
\]
i.e., **\(u\)** is a control input that rescales the forward rate constant in a way consistent with transition‑state theory \(\big(k \propto e^{-\Delta G^\ddagger/RT}\big)\).

---

## Why *u* is natural in the **log‑linear** framework

Let the reaction quotient be \(Q\) (e.g., \(Q = B/A\) for \(A \rightleftharpoons B\)), and write
\[
x \equiv \ln\!\frac{Q}{K_{eq}}\,.
\]
With chemostatted drive folded into the kinetics, the **near‑equilibrium** dynamics become
\[
\frac{d}{dt}\ln Q \;=\; -k\,\ln\!\frac{Q}{K_{eq}'} 
\;\approx\; -k\,\ln\!\frac{Q}{K_{eq}} \;+\; \underbrace{\Big(\frac{\partial \ln K_{eq}'}{\partial u}\Big)}_{\text{constant}}\,u,
\]
which is the **log‑linear form**:
\[
\dot{x} \;=\; -K\,x \;+\; b\,u(t).
\]
Here \(u\) enters **additively** with a clear thermodynamic interpretation (log‑ratio / chemical‑potential bias).

---

## Bottom line
- **Mass‑action:** \(u\) is a principled way to encode how external cofactors (kept at fixed levels) **rescale effective rate constants** after coarse‑graining. Writing \(k_f=k_f^{0}e^{u}\) follows directly from eliminating the explicit cofactor variables.
- **Log‑linear:** \(u\) is the **natural drive**—a log‑ratio/free‑energy term that couples **linearly** into \(\dot{x}\).

If energy units are preferred, multiply by \(RT\): \(RT\,u = \Delta\mu\).
