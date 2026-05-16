"""Step 3: Diagnostics, posterior predictive, R0, residuals plots
4-parameter model: beta, gamma, rho, phi
"""
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt, json, os

# Get the directory of this script for relative paths
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(SCRIPT_DIR)

np.random.seed(42)
PLOT_DIR = os.path.join(PROJECT_DIR, "outputs", "plots")
d = np.load(os.path.join(PROJECT_DIR, "outputs", "step1_data.npz"))
weeks, cases, n_weeks, N = d['weeks'], d['cases'], int(d['n_weeks']), int(d['N'])
m = np.load(os.path.join(PROJECT_DIR, "outputs", "mcmc_results.npz"))
chains_nat = m['chains_nat']  # (4, 3000, 4) => beta, gamma, rho, phi
n_chains = chains_nat.shape[0]
plt.rcParams.update({'font.family': 'serif', 'font.size': 11, 'axes.grid': True,
                     'grid.alpha': 0.3, 'grid.color': 'gray', 'figure.facecolor': 'white'})

param_names = ['beta', 'gamma', 'rho', 'phi']

def compute_rhat(cp):
    m2_chains = []
    for c in range(cp.shape[0]):
        h = cp.shape[1] // 2
        m2_chains.append(cp[c, :h]); m2_chains.append(cp[c, h:2*h])
    sc = np.array(m2_chains); mc = sc.shape[0]; nc = sc.shape[1]
    cm = sc.mean(axis=1); gm = cm.mean()
    B = nc * np.sum((cm - gm) ** 2) / (mc - 1)
    W = np.mean(np.var(sc, axis=1, ddof=1))
    return np.sqrt(((1 - 1/nc) * W + B/nc) / W) if W > 0 else 999

def compute_ess(cp):
    x = cp.flatten(); n = len(x); mu = x.mean(); v = x.var()
    if v == 0: return n
    ml = min(n // 2, 300)
    ac = np.correlate(x - mu, x - mu, mode='full')
    ac = ac[n-1:n-1+ml] / (v * n)
    tau = 1.0
    for k in range(1, ml - 1, 2):
        rs = ac[k] + ac[k+1]
        if rs < 0: break
        tau += 2 * rs
    return n / tau

# Compute diagnostics for all 4 parameters
rhat, ess, summary = {}, {}, {}
all_post = np.concatenate([chains_nat[c] for c in range(n_chains)], axis=0)
for j, nm in enumerate(param_names):
    rhat[nm] = float(compute_rhat(chains_nat[:, :, j]))
    ess[nm] = float(compute_ess(chains_nat[:, :, j]))
    v = all_post[:, j]
    summary[nm] = {'mean': float(np.mean(v)), 'std': float(np.std(v)),
        'q2.5': float(np.percentile(v, 2.5)), 'q25': float(np.percentile(v, 25)),
        'q50': float(np.percentile(v, 50)), 'q75': float(np.percentile(v, 75)),
        'q97.5': float(np.percentile(v, 97.5)), 'rhat': rhat[nm], 'ess': ess[nm]}
    print(f"{nm}: mean={summary[nm]['mean']:.4f} median={summary[nm]['q50']:.4f} "
          f"95%CI=[{summary[nm]['q2.5']:.4f},{summary[nm]['q97.5']:.4f}] "
          f"R-hat={rhat[nm]:.4f} ESS={ess[nm]:.0f}")

# R0 = beta / gamma
R0 = all_post[:, 0] / all_post[:, 1]
R0s = {'mean': float(np.mean(R0)), 'std': float(np.std(R0)),
       'q2.5': float(np.percentile(R0, 2.5)), 'q50': float(np.percentile(R0, 50)),
       'q97.5': float(np.percentile(R0, 97.5))}
print(f"R0: mean={R0s['mean']:.3f} median={R0s['q50']:.3f} 95%CI=[{R0s['q2.5']:.3f},{R0s['q97.5']:.3f}]")

# Plot 4: Trace plots (4 panels)
cc = ['black', '#444444', '#777777', '#AAAAAA']
fig, axes = plt.subplots(4, 1, figsize=(10, 10), sharex=True)
for j, nm in enumerate(param_names):
    for c in range(n_chains):
        axes[j].plot(chains_nat[c, :, j], color=cc[c], alpha=0.7, lw=0.4,
                     label=f'Chain {c+1}' if j == 0 else None)
    axes[j].set_ylabel(nm)
    axes[j].set_title(f"Trace: {nm} (R-hat = {rhat[nm]:.4f})")
axes[3].set_xlabel("Iteration")
axes[0].legend(loc='upper right', fontsize=8, frameon=True, edgecolor='black')
plt.tight_layout()
fig.savefig(f"{PLOT_DIR}/04_trace_plots.png", dpi=200, bbox_inches='tight')
plt.close()

# Plot 5: Posterior histograms (5 panels: beta, gamma, rho, phi, R0)
fig, axes = plt.subplots(1, 5, figsize=(16, 3.5))
for j, nm in enumerate(param_names):
    axes[j].hist(all_post[:, j], bins=50, color='gray', edgecolor='black', lw=0.5, density=True)
    axes[j].axvline(summary[nm]['q50'], color='black', lw=2, label='Median')
    axes[j].axvline(summary[nm]['q2.5'], color='black', lw=1, ls='--', label='95% CI')
    axes[j].axvline(summary[nm]['q97.5'], color='black', lw=1, ls='--')
    axes[j].set_xlabel(nm); axes[j].set_ylabel("Density")
    axes[j].set_title(f"Posterior: {nm}")
    axes[j].legend(fontsize=6, frameon=True, edgecolor='black')
axes[4].hist(R0, bins=50, color='gray', edgecolor='black', lw=0.5, density=True)
axes[4].axvline(R0s['q50'], color='black', lw=2, label='Median')
axes[4].axvline(R0s['q2.5'], color='black', lw=1, ls='--', label='95% CI')
axes[4].axvline(R0s['q97.5'], color='black', lw=1, ls='--')
axes[4].set_xlabel("R0"); axes[4].set_ylabel("Density"); axes[4].set_title("Posterior: R0")
axes[4].legend(fontsize=6, frameon=True, edgecolor='black')
plt.tight_layout()
fig.savefig(f"{PLOT_DIR}/05_posterior_histograms.png", dpi=200, bbox_inches='tight')
plt.close()

# Plot 6: Pair plots (4x4)
fig, axes = plt.subplots(4, 4, figsize=(12, 12))
for i in range(4):
    for j in range(4):
        if i == j:
            axes[i][j].hist(all_post[:, i], bins=40, color='gray', edgecolor='black', lw=0.3)
            axes[i][j].set_xlabel(param_names[i])
        elif i > j:
            idx = np.random.choice(len(all_post), min(1500, len(all_post)), replace=False)
            axes[i][j].scatter(all_post[idx, j], all_post[idx, i], s=1, color='black', alpha=0.3)
            axes[i][j].set_xlabel(param_names[j]); axes[i][j].set_ylabel(param_names[i])
        else:
            axes[i][j].axis('off')
plt.suptitle("Posterior Pair Plots", fontsize=14, y=1.01)
plt.tight_layout()
fig.savefig(f"{PLOT_DIR}/06_pair_plots.png", dpi=200, bbox_inches='tight')
plt.close()

# Posterior predictive with reporting rate
def solve_sir_fast(beta, gamma, N, nw):
    dt = 0.25; S, I = float(N - 1), 1.0; Io = np.empty(nw)
    for w in range(nw):
        for _ in range(4):
            inf = beta * S * I / N * dt; rec = gamma * I * dt
            if inf > S: inf = S
            if rec > I: rec = I
            S -= inf; I += inf - rec
        Io[w] = I
    return Io

n_pp = 200
pp_c = np.zeros((n_pp, n_weeks))   # rho * I(t) curves
pp_s = np.zeros((n_pp, n_weeks))   # simulated observations
idx_pp = np.random.choice(len(all_post), n_pp, replace=False)
for i, ix in enumerate(idx_pp):
    b, g, rho, phi = all_post[ix]
    Ip = solve_sir_fast(b, g, N, n_weeks)
    mu_curve = rho * Ip  # expected cases = rho * I(t)
    pp_c[i] = mu_curve
    for t in range(n_weeks):
        mu = max(mu_curve[t], 0.01)
        prob = phi / (phi + mu)
        try:
            pp_s[i, t] = np.random.negative_binomial(phi, prob)
        except:
            pp_s[i, t] = mu

med = np.median(pp_c, axis=0)
q5, q25, q75, q95 = [np.percentile(pp_c, p, axis=0) for p in [5, 25, 75, 95]]

# Plot 7: Posterior predictive intervals
fig, ax = plt.subplots(figsize=(10, 5))
ax.fill_between(weeks, q5, q95, color='#CCCCCC', alpha=0.7, label='90% CI')
ax.fill_between(weeks, q25, q75, color='#888888', alpha=0.7, label='50% CI')
ax.plot(weeks, med, color='black', lw=2, label='Median prediction')
ax.scatter(weeks, cases, color='black', s=30, zorder=5, label='Observed data')
ax.set_xlabel("Week Index"); ax.set_ylabel("Number of Cases")
ax.set_title("Posterior Predictive Check: Model Fit vs Observed Data")
ax.legend(frameon=True, edgecolor='black'); ax.set_ylim(bottom=0)
plt.tight_layout()
fig.savefig(f"{PLOT_DIR}/07_posterior_predictive.png", dpi=200, bbox_inches='tight')
plt.close()

# Plot 8: Simulated observations
sm = np.median(pp_s, axis=0)
s5, s25, s75, s95 = [np.percentile(pp_s, p, axis=0) for p in [5, 25, 75, 95]]
fig, ax = plt.subplots(figsize=(10, 5))
ax.fill_between(weeks, s5, s95, color='#CCCCCC', alpha=0.7, label='90% CI (simulated)')
ax.fill_between(weeks, s25, s75, color='#888888', alpha=0.7, label='50% CI (simulated)')
ax.plot(weeks, sm, color='black', lw=2, label='Median (simulated)')
ax.scatter(weeks, cases, color='black', s=30, zorder=5, label='Observed data')
ax.set_xlabel("Week Index"); ax.set_ylabel("Number of Cases")
ax.set_title("Posterior Predictive Check: Simulated Observations vs Actual Data")
ax.legend(frameon=True, edgecolor='black'); ax.set_ylim(bottom=0)
plt.tight_layout()
fig.savefig(f"{PLOT_DIR}/08_posterior_predictive_sim.png", dpi=200, bbox_inches='tight')
plt.close()

# Plot 9: R0 histogram
fig, ax = plt.subplots(figsize=(8, 5))
ax.hist(R0, bins=60, color='gray', edgecolor='black', lw=0.5, density=True)
ax.axvline(R0s['q50'], color='black', lw=2, label=f"Median = {R0s['q50']:.2f}")
ax.axvline(R0s['q2.5'], color='black', lw=1.5, ls='--',
           label=f"95% CI: [{R0s['q2.5']:.2f}, {R0s['q97.5']:.2f}]")
ax.axvline(R0s['q97.5'], color='black', lw=1.5, ls='--')
ax.axvline(1.0, color='black', lw=1, ls=':', label='R0 = 1 (epidemic threshold)')
ax.set_xlabel("R0 (Basic Reproduction Number)"); ax.set_ylabel("Density")
ax.set_title("Posterior Distribution of R0")
ax.legend(frameon=True, edgecolor='black')
plt.tight_layout()
fig.savefig(f"{PLOT_DIR}/09_R0_distribution.png", dpi=200, bbox_inches='tight')
plt.close()

# Plot 10: Residuals
res = cases - med
fig, axes = plt.subplots(2, 1, figsize=(10, 7))
axes[0].bar(weeks, res, color='gray', edgecolor='black', lw=0.5)
axes[0].axhline(0, color='black', lw=1)
axes[0].set_xlabel("Week Index"); axes[0].set_ylabel("Residual (Observed - Predicted)")
axes[0].set_title("Residuals Over Time")
axes[1].scatter(med, res, color='black', s=25)
axes[1].axhline(0, color='black', lw=1)
axes[1].set_xlabel("Predicted Cases"); axes[1].set_ylabel("Residual")
axes[1].set_title("Residuals vs Predicted Values")
plt.tight_layout()
fig.savefig(f"{PLOT_DIR}/10_residuals.png", dpi=200, bbox_inches='tight')
plt.close()

# Posterior predictive fit check
ss_res = np.sum((cases - med) ** 2)
ss_tot = np.sum((cases - cases.mean()) ** 2)
r2 = 1 - ss_res / ss_tot
rmse = np.sqrt(np.mean((cases - med) ** 2))
print(f"\nFIT CHECK: R2={r2:.4f}, RMSE={rmse:.2f}")
print(f"Peak predicted: week {np.argmax(med)}, Peak observed: week {cases.argmax()}")
print(f"STEP 3 FIT: {'PASS' if r2 > 0.7 else 'FAIL'}")

# Save all results as JSON
results = {
    'N': N, 'n_weeks': int(n_weeks), 'total_cases': int(cases.sum()),
    'peak_cases': int(cases.max()), 'peak_week': int(cases.argmax()),
    'summary': summary, 'R0': R0s,
    'rhat': {k: float(v) for k, v in rhat.items()},
    'ess': {k: float(v) for k, v in ess.items()},
    'fit': {'R2': r2, 'RMSE': rmse}
}
with open(os.path.join(PROJECT_DIR, "outputs", "results.json"), 'w') as f:
    json.dump(results, f, indent=2)

print("Step 3 done: All plots 04-10 and results saved.")
