# receptor_competition/receptor_competition.py
import numpy as np
import matplotlib.pyplot as plt

# Parameters (units arbitrary; think nM)
Kd_L = 10.0     # agonist Kd
Kd_C = 100.0    # antagonist Kd
Keq_L = 1.0 / Kd_L
Keq_C = 1.0 / Kd_C

# Normalized grids (dimensionless)
L_over_Kd = np.logspace(-2, 4, 200)    # [L]/Kd_L from 0.01 to 10,000
C_over_Kd_levels = [0.0, 0.3, 1.0, 3.0]

def frac_L_bound_norm(Lk, Ck):
    # In normalized units, f = Lk / (1 + Lk + Ck)
    return Lk / (1.0 + Lk + Ck)

plt.figure(figsize=(6,4.5))
for Ck in C_over_Kd_levels:
    f = frac_L_bound_norm(L_over_Kd, Ck)
    plt.semilogx(L_over_Kd, 100*f, label=f"C/Kd_C={Ck:g}")

plt.xlabel("[L]/Kd_L (dimensionless)")
plt.ylabel("% receptor bound by L")
plt.title("Competitive binding:\nDose–response shifts with antagonist (normalized units)")
plt.legend(frameon=False)
plt.tight_layout(pad=1.2)
plt.savefig("receptor_competition_dose_response.png", dpi=200, bbox_inches="tight")
print("Wrote receptor_competition_dose_response.png")

