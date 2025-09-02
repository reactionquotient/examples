# Longevity Tradeoff

Two coupled reaction-quotient “modes” model a simple maintenance-vs-growth tradeoff in cells:

- **Q₁:** maintenance / autophagy (AMPK/Sirtuin-leaning)  
- **Q₂:** anabolic growth (mTOR-leaning)

The dynamics follow the same log-linear RQ form used elsewhere in this repo:

\[
\frac{d}{dt}\ln Q = -K\,\ln(Q/K_{eq}) + u
\]

where \(Q=[Q_1,Q_2]\), \(K\) is a positive-stable rate/coupling matrix, and \(u\) is a constant input representing “drives.”

---

## What it shows

We compare two operating conditions:

1) **Ad libitum** (fed): lower NAD⁺/NADH, high nutrient/mTOR drive  
2) **Caloric restriction** (CR): higher NAD⁺/NADH, lower nutrient/mTOR drive

Under fed conditions, the growth mode rises while maintenance lags; under CR, maintenance dominates and growth subsides.

---

## Files

- `longevity_tradeoff.py` — runs the simulation and produces a single figure
- `longevity_tradeoff.png` — saved plot (created on run)
- `README.md` — this file

---

## Quickstart

```bash
# from the repo root or this folder
pip install numpy scipy matplotlib

python longevity_tradeoff/longevity_tradeoff.py
# or, if you're already in this folder:
# python longevity_tradeoff.py
```

You’ll get a figure comparing \(Q_1/K_{eq,1}\) (maintenance) and \(Q_2/K_{eq,2}\) (growth) over time for both conditions.

---

## Model details

**State:** \(x=\ln(Q/K_{eq})\)  
**ODE:** \(\dot{x} = -Kx + u\) with constant \(u\)  
**Closed form:**  
\[
x(t) = e^{-Kt}\,x_0 + K^{-1}\left(I - e^{-Kt}\right)u
\]
The script uses `scipy.linalg.expm` to evaluate \(e^{-Kt}\).

**Inputs (toy mapping):**
```python
u1 = k_nad  * ln(NAD+/NADH)             # maintenance drive (AMPK/Sirtuins)
u2 = k_nutr * ln(nutrient) - k_ampk*ln(NAD+/NADH)  # growth drive (mTOR) opposed by NAD+
```

**Default parameters (illustrative only):**
```text
K = [[1.0, 0.20],
     [0.20, 1.2]]     # positive-stable with weak coupling
Keq = [1.0, 1.0]
x0  = [0.0, 0.0]      # start at equilibrium in log-space
```

**Conditions:**
```text
Ad libitum:           NAD+/NADH = 1.2, nutrient = 5.0
Caloric restriction:  NAD+/NADH = 3.0, nutrient = 1.2
```

---

## How to tweak

- **Tradeoff strength:** increase/decrease the off-diagonal \(K_{12}=K_{21}=c\).  
- **Timescales:** change the diagonals of \(K\) (higher = faster relaxation).  
- **Energy state:** adjust `nad_ratio` to push maintenance up/down.  
- **Nutrient drive:** adjust `nutrient` to push growth up/down.  
- **Opposition strength:** tune `k_ampk` to make NAD⁺ more strongly suppress growth.

---

## Notes

- This is a **didactic toy model** to illustrate the log-linear RQ framework, not a calibrated biological model.  
- Parameter choices are arbitrary but chosen for clear qualitative behavior and numerical stability.
