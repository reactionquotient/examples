# Complex example: fit a *coupled* log-linear K matrix from data
# ---------------------------------------------------------------
# We generate data from a 3-species reversible *cycle* (A <-> B <-> C <-> A).
# The forward rate on A<->B is modulated by a stepped input u(t) (e.g., ATP/ADP effect).
# We then *fit* the log-linear model:
#         d x / dt = -K x + b * u(t)             where x = ln(Q / Keq),  Q = [B/A, C/B, A/C]
# from the simulated "measurements" (with noise). Finally we validate by simulating
# the fitted model forward and overlaying predictions vs. data.
#
# If Tellurium is available, we use an Antimony/RoadRunner model with EVENTS for u-steps.
# Otherwise we fall back to a pure-Python RK4 mass-action simulator that accepts u(t).
#
# Artifacts saved under ./output/: CSVs + figures + fitted_K.json

import numpy as np
import math
import json
import os
import matplotlib.pyplot as plt

# ---------- Try Tellurium ----------
try:
    import tellurium as te
    TELLURIUM_AVAILABLE = True
except Exception as e:
    TELLURIUM_AVAILABLE = False
    tellurium_err = str(e)

# ---------- Utility: RK4 ODE ----------
def rk4(f, y0, t, dt):
    y = np.zeros((len(t), len(y0)), dtype=float)
    y[0] = y0
    for i in range(len(t)-1):
        ti = t[i]
        yi = y[i]
        k1 = f(ti, yi)
        k2 = f(ti + 0.5*dt, yi + 0.5*dt*k1)
        k3 = f(ti + 0.5*dt, yi + 0.5*dt*k2)
        k4 = f(ti + dt, yi + dt*k3)
        y[i+1] = yi + (dt/6.0)*(k1 + 2*k2 + 2*k3 + k4)
    return y

def safe_ratio(a,b,eps=1e-12):
    return a / (b + eps)

# ---------- Experiment design ----------
np.random.seed(7)

T = 60.0
dt = 0.01
t = np.arange(0.0, T+1e-12, dt)

# Input u(t): step, return, negative step (drive/inhibition)
def u_of_t(tt):
    if tt < 10.0: return 0.0
    if tt < 30.0: return 0.8     # strong drive
    if tt < 45.0: return 0.0
    return -0.6                  # inhibition
u = np.array([u_of_t(tt) for tt in t])

# Base equilibrium constants for the three reversible edges (A<->B, B<->C, C<->A)
Keq1, Keq2, Keq3 = 2.0, 1.5, 0.8

# Base relaxation scales (pick reverse rates, compute forward by Keq)
s1, s2, s3 = 0.8, 1.1, 0.6
k1r, k2r, k3r = s1, s2, s3
k1f_base, k2f, k3f = Keq1*k1r, Keq2*k2r, Keq3*k3r

# The input multiplies the *forward* rate of A->B: k1f = k1f_base * exp(u)
# Initial condition:
A0, B0, C0 = 0.9, 0.05, 0.05

# ---------- Simulate mass-action to produce "data" ----------
if TELLURIUM_AVAILABLE:
    ant = f"""
    model ABC_cycle()
      // species
      A = {A0}; B = {B0}; C = {C0};
      // base kinetic params
      k1f_base = {k1f_base}; k1r = {k1r};
      k2f = {k2f}; k2r = {k2r};
      k3f = {k3f}; k3r = {k3r};
      // input (dimensionless), used as assignment into k1f
      u = 0;
      k1f := k1f_base * exp(u);
      // reactions
      R1f: A -> B; k1f*A;
      R1r: B -> A; k1r*B;
      R2f: B -> C; k2f*B;
      R2r: C -> B; k2r*C;
      R3f: C -> A; k3f*C;
      R3r: A -> C; k3r*A;

      // events for u-steps
      at (time >= 10): u = 0.8;
      at (time >= 30): u = 0.0;
      at (time >= 45): u = -0.6;
    end
    """
    rr = te.loada(ant)
    res = rr.simulate(0, T, int(T/dt)+1)
    t_data = res[:,0]
    # Order: time, A, B, C, ... (RR may include fluxes). We'll map by species IDs to be safe.
    species = rr.model.getFloatingSpeciesIds()
    A_idx = list(species).index("A")+1
    B_idx = list(species).index("B")+1
    C_idx = list(species).index("C")+1
    A = res[:,A_idx]; B = res[:,B_idx]; C = res[:,C_idx]
    # Reconstruct u(t) from events for completeness
    u_data = np.zeros_like(t_data)
    u_data[t_data>=10] = 0.8
    u_data[t_data>=30] = 0.0
    u_data[t_data>=45] = -0.6
else:
    # Pure-Python mass-action with time-varying input on k1f
    def f(tt, y):
        A, B, C = y
        k1f = k1f_base * math.exp(u_of_t(tt))
        dA = -k1f*A + k1r*B + k3f*C - k3r*A
        dB =  k1f*A - k1r*B - k2f*B + k2r*C
        dC =  k2f*B - k2r*C - k3f*C + k3r*A
        return np.array([dA, dB, dC], dtype=float)
    Y = rk4(f, np.array([A0,B0,C0], dtype=float), t, dt)
    A, B, C = Y[:,0], Y[:,1], Y[:,2]
    t_data = t
    u_data = u

# Add measurement noise (e.g., 2% CV)
noise_scale = 0.02
A_meas = A * (1 + noise_scale*np.random.randn(len(A)))
B_meas = B * (1 + noise_scale*np.random.randn(len(B)))
C_meas = C * (1 + noise_scale*np.random.randn(len(C)))

# Compute quotients and x = ln(Q/Keq)
Q1 = safe_ratio(B_meas, A_meas);  Q2 = safe_ratio(C_meas, B_meas);  Q3 = safe_ratio(A_meas, C_meas)
Q = np.vstack([Q1,Q2,Q3]).T
Keq_vec = np.array([Keq1, Keq2, Keq3], dtype=float)
x = np.log(Q / Keq_vec)

# Numerical derivative x_dot (central diff, skip endpoints)
x_dot = np.zeros_like(x)
x_dot[1:-1] = (x[2:] - x[:-2]) / (2*dt)
# Trim ends to avoid edge effects
mask = np.ones(len(t_data), dtype=bool)
mask[[0,-1]] = False

# Build regression: x_dot ≈ -K x + b*u
X_feat = np.hstack([ -x[mask],  u_data[mask,None] ])   # shape (N, 3+1)
Y_resp = x_dot[mask]                                    # shape (N, 3)

# Ridge regression (small λ) for stability
lam = 1e-3
Z = X_feat
W = np.linalg.solve(Z.T@Z + lam*np.eye(Z.shape[1]), Z.T@Y_resp)  # shape (4,3)

# Unpack K_hat and b_hat
K_hat = W[:3,:].T   # (3x3)
b_hat = W[3,:]      # (3,)

# Forward-simulate fitted log-linear model with the same u(t)
def sim_loglinear(K, b, x0, t, dt, u_of_t):
    xsim = np.zeros((len(t), len(x0)))
    xsim[0] = x0
    for i in range(len(t)-1):
        u_val = u_of_t(t[i])
        dx = -K @ xsim[i] + b * u_val
        xsim[i+1] = xsim[i] + dt*dx
    return xsim

# We'll reconstruct a continuous u(t) function from the data grid
def u_func(tt):
    # identical to generation
    if tt < 10.0: return 0.0
    if tt < 30.0: return 0.8
    if tt < 45.0: return 0.0
    return -0.6

x0 = x[0]
x_pred = sim_loglinear(K_hat, b_hat, x0, t_data, dt, u_func)
Q_pred = np.exp(x_pred) * Keq_vec

# Evaluation metrics (RMSE for each quotient)
rmse = np.sqrt(np.mean((Q_pred - Q)**2, axis=0))

# Create output dir
outdir = "output"
os.makedirs(outdir, exist_ok=True)

# Save datasets
import pandas as pd
df_conc = pd.DataFrame({"time": t_data, "A": A_meas, "B": B_meas, "C": C_meas, "u": u_data})
df_quot = pd.DataFrame({"time": t_data, "Q1": Q[:,0], "Q2": Q[:,1], "Q3": Q[:,2]})
df_conc.to_csv(os.path.join(outdir, "concentrations.csv"), index=False)
df_quot.to_csv(os.path.join(outdir, "quotients.csv"), index=False)

# Save parameters
with open(os.path.join(outdir, "fitted_K.json"), "w") as f:
    json.dump({
        "K_hat": K_hat.tolist(),
        "b_hat": b_hat.tolist(),
        "ridge_lambda": lam,
        "notes": "Fitted from x_dot ≈ -K x + b u with central-diff derivative; 3-quotient cycle."
    }, f, indent=2)

# ---------- Figures ----------

# 1) Time courses: measured Q vs model Q_pred
plt.figure()
plt.plot(t_data, Q[:,0], label="Q1 data")
plt.plot(t_data, Q_pred[:,0], label="Q1 model", linestyle="--")
plt.title("Q1=B/A: data vs fitted log-linear model")
plt.xlabel("time")
plt.ylabel("Q1")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(outdir, "Q1_fit.png"))
plt.show()

plt.figure()
plt.plot(t_data, Q[:,1], label="Q2 data")
plt.plot(t_data, Q_pred[:,1], label="Q2 model", linestyle="--")
plt.title("Q2=C/B: data vs fitted log-linear model")
plt.xlabel("time")
plt.ylabel("Q2")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(outdir, "Q2_fit.png"))
plt.show()

plt.figure()
plt.plot(t_data, Q[:,2], label="Q3 data")
plt.plot(t_data, Q_pred[:,2], label="Q3 model", linestyle="--")
plt.title("Q3=A/C: data vs fitted log-linear model")
plt.xlabel("time")
plt.ylabel("Q3")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(outdir, "Q3_fit.png"))
plt.show()

# 2) Derivative fit quality: x_dot vs (-Kx + b u)
lhs = x_dot[mask]                      # true derivative (from data)
rhs = (- (K_hat @ x[mask].T).T + b_hat[None,:] * u_data[mask,None])  # model RHS
# R^2 for each component
def r2(y_true, y_pred):
    ss_res = np.sum((y_true - y_pred)**2, axis=0)
    ss_tot = np.sum((y_true - np.mean(y_true, axis=0))**2, axis=0)
    return 1 - ss_res/ss_tot

r2_vec = r2(lhs, rhs)

plt.figure()
plt.plot(t_data[mask], lhs[:,0], label="x1_dot data")
plt.plot(t_data[mask], rhs[:,0], label="x1_dot model", linestyle="--")
plt.title("Derivative fit (component 1)")
plt.xlabel("time")
plt.ylabel("dx1/dt")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(outdir, "derivative_fit_x1.png"))
plt.show()

plt.figure()
plt.plot(t_data[mask], lhs[:,1], label="x2_dot data")
plt.plot(t_data[mask], rhs[:,1], label="x2_dot model", linestyle="--")
plt.title("Derivative fit (component 2)")
plt.xlabel("time")
plt.ylabel("dx2/dt")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(outdir, "derivative_fit_x2.png"))
plt.show()

plt.figure()
plt.plot(t_data[mask], lhs[:,2], label="x3_dot data")
plt.plot(t_data[mask], rhs[:,2], label="x3_dot model", linestyle="--")
plt.title("Derivative fit (component 3)")
plt.xlabel("time")
plt.ylabel("dx3/dt")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(outdir, "derivative_fit_x3.png"))
plt.show()

# 3) Show fitted K and eigenvalues
eigvals, eigvecs = np.linalg.eig(K_hat)
plt.figure()
plt.plot(np.arange(1,4), np.diag(K_hat), marker="o", linestyle="")
plt.title("Fitted K: diagonal entries")
plt.xlabel("index")
plt.ylabel("K_ii")
plt.tight_layout()
plt.savefig(os.path.join(outdir, "K_diagonal.png"))
plt.show()

plt.figure()
plt.plot(np.arange(1,4), np.sort(eigvals.real), marker="o", linestyle="")
plt.title("Fitted K eigenvalues (real parts)")
plt.xlabel("index")
plt.ylabel("Re(λ)")
plt.tight_layout()
plt.savefig(os.path.join(outdir, "K_eigs.png"))
plt.show()

# Print summary
print("Tellurium available:", TELLURIUM_AVAILABLE)
if not TELLURIUM_AVAILABLE:
    print("  (Used pure-Python mass-action simulator.)")

print("\nFITTED PARAMETERS")
print("K_hat =")
print(K_hat)
print("b_hat =", b_hat)
print("\nRMSE on Q components:", rmse)
print("R^2 for derivative fit components:", r2_vec)

print("\nConcentrations (noisy) used for fitting:")
print(df_conc.head(10))
print("\nQuotients computed for fitting:")
print(df_quot.head(10))

print("\nSaved artifacts in:", outdir)
for name in sorted(os.listdir(outdir)):
    print(" -", os.path.join(outdir, name))

