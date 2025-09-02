# Ligand–Receptor Competitive Binding (dose–response shifts)

**What this shows**  
A minimal competitive-binding model where ligand **L** and antagonist/competitor **C** vie for the same receptor **R**. The figure plots the percentage of receptors occupied by **L** as a function of **[L]**, for several fixed **[C]** values. The hallmark of competitive antagonism appears: **parallel rightward shifts** of the dose–response curve as [C] increases, with the **same top asymptote**.

**Model (normalized units)**  
Let
- `Lk = [L]/Kd_L`
- `Ck = [C]/Kd_C`

Under ligand excess (standard binding assay assumption), the steady‑state fraction of receptors occupied by **L** is
\[
f_L = \frac{Lk}{1 + Lk + Ck}.
\]
This is the classic competitive‑binding form. In these units, the **dose ratio** follows
\[
\mathrm{EC}_{50}(C)/\mathrm{EC}_{50}(0) = 1 + Ck.
\]

**Files**
- `receptor_competition.py` — generates the plot using normalized units.
- `receptor_competition_dose_response.png` — output figure.

**How to run**
```bash
python receptor_competition/receptor_competition.py
```

**What to look for**
- Increasing `Ck` shifts the curve to the right without lowering the maximal occupancy.
- At `Ck = 0`, the half‑max occurs at `Lk ≈ 1` (i.e., `[L] ≈ Kd_L`). With competitor, the EC50 increases by `1 + Ck`.

**Play with it**
- Change the list `C_over_Kd_levels` to explore different antagonist concentrations (e.g., `0, 0.1, 0.3, 1, 3`).
- Extend the `L_over_Kd` range to confirm the plateau (e.g., up to `1e5`).
- (Optional) Produce a quick Schild‑style check by measuring `EC50` at different `Ck` and verifying that `EC50(C)/EC50(0) - 1` equals `Ck` (slope ≈ 1 on a log–log plot).

**Notes**
- This example is about **occupancy**, not downstream efficacy. Partial agonism, cooperativity (Hill slopes ≠ 1), or non‑competitive mechanisms are intentionally omitted to keep the example focused and minimal.
