"""Step 2: MCMC for 4-parameter SIR model
Params: R0, gamma, rho (reporting rate), phi (NegBin overdispersion)
All sampled in log-space: (log_R0, log_gamma, log_rho, log_phi)
Observation model: cases(t) ~ NegBin(rho * I(t), phi)
"""
import numpy as np
from math import lgamma

# Get the directory of this script for relative paths
import os
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(SCRIPT_DIR)

np.random.seed(42)
d = np.load(os.path.join(PROJECT_DIR, "outputs", "step1_data.npz"))
weeks, cases, n_weeks, N = d['weeks'], d['cases'], int(d['n_weeks']), int(d['N'])

def solve_sir(beta, gamma, N, nw):
    """Euler with dt=0.25, 4 sub-steps per week (fast for MCMC)"""
    dt = 0.25; S, I = float(N - 1), 1.0; Io = np.empty(nw)
    for w in range(nw):
        for _ in range(4):
            inf = beta * S * I / N * dt
            rec = gamma * I * dt
            if inf > S: inf = S
            if rec > I: rec = I
            S -= inf; I += inf - rec
        Io[w] = I
    return Io

def neg_binom_ll(k_arr, mu_arr, phi):
    r = phi; ll = 0.0
    for i in range(len(k_arr)):
        mu = max(mu_arr[i], 0.01); k = int(k_arr[i])
        p = r / (r + mu)
        if p <= 0 or p >= 1: return -1e10
        ll += lgamma(k + r) - lgamma(k + 1) - lgamma(r) + r * np.log(p) + k * np.log(1 - p)
    return ll

def log_post(t):
    """t = [log_R0, log_gamma, log_rho, log_phi]"""
    lR0, lg, lrho, lp = t
    R0 = np.exp(lR0); gamma = np.exp(lg); rho = np.exp(lrho); phi = np.exp(lp)
    beta = R0 * gamma
    # Bounds
    if beta > 10 or gamma > 5 or phi > 100 or phi < 0.01 or R0 < 0.5 or R0 > 10:
        return -1e10
    if rho < 0.001 or rho > 1.0:
        return -1e10
    # Priors (all log-normal, specified on log-scale)
    lpr = -0.5 * ((lR0 - np.log(1.7)) / 0.25) ** 2    # R0 ~ LogNormal(log(1.7), 0.25)
    lpr += -0.5 * ((lg - np.log(0.44)) / 0.3) ** 2     # gamma ~ LogNormal(log(0.44), 0.3)
    lpr += -0.5 * ((lrho - np.log(0.05)) / 0.5) ** 2   # rho ~ LogNormal(log(0.05), 0.5)
    lpr += -phi + lp                                     # phi ~ Exp(1), + Jacobian for log_phi
    # Likelihood
    Ip = solve_sir(beta, gamma, N, n_weeks)
    mu_arr = rho * Ip  # expected cases = rho * I(t)
    return lpr + neg_binom_ll(cases, mu_arr, phi)

# Grid search for initialization
print("Grid search for 4 parameters...")
best_lp, best_t = -1e20, None
for r0 in np.linspace(1.2, 2.2, 12):
    for g in np.linspace(0.3, 0.8, 10):
        for rho in [0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.10]:
            for p in [0.5, 1.0, 2.0, 5.0]:
                t = [np.log(r0), np.log(g), np.log(rho), np.log(p)]
                lp = log_post(t)
                if lp > best_lp:
                    best_lp, best_t = lp, t[:]
print(f"Grid best: R0={np.exp(best_t[0]):.3f} gamma={np.exp(best_t[1]):.3f} rho={np.exp(best_t[2]):.4f} phi={np.exp(best_t[3]):.3f}")

# Phase 1: Tuning chain to learn proposal covariance
print("Tuning chain (6000 iterations)...")
ndim = 4
sc = np.array([0.04, 0.04, 0.06, 0.06])
theta = np.array(best_t)
np.random.seed(999)
lp_cur = log_post(theta)
rm = theta.copy(); rc = np.diag(sc ** 2)
tune = np.empty((6000, ndim))
for i in range(6000):
    if i > 300:
        cov = rc * (2.38 ** 2 / ndim) + 1e-8 * np.eye(ndim)
        try:
            L = np.linalg.cholesky(cov)
            prop = theta + L @ np.random.randn(ndim)
        except:
            prop = theta + sc * np.random.randn(ndim)
    else:
        prop = theta + sc * np.random.randn(ndim)
    lp_p = log_post(prop)
    if np.log(np.random.rand()) < (lp_p - lp_cur):
        theta = prop; lp_cur = lp_p
    tune[i] = theta
    delta = theta - rm
    rm += delta / (i + 1)
    if i > 100:
        rc = rc * (1 - 1 / (i + 1)) + np.outer(delta, delta) / (i + 1)

learned_cov = np.cov(tune[3000:].T) * (2.38 ** 2 / ndim) * 0.35
print(f"Learned proposal sd: {np.sqrt(np.diag(learned_cov))}")

# Phase 2: Production chains
def run_chain(seed, init, cov, n_iter):
    np.random.seed(seed)
    theta = np.array(init)
    L = np.linalg.cholesky(cov + 1e-8 * np.eye(ndim))
    samp = np.empty((n_iter, ndim)); lp_cur = log_post(theta); acc = 0
    for i in range(n_iter):
        prop = theta + L @ np.random.randn(ndim)
        lp_p = log_post(prop)
        if np.log(np.random.rand()) < (lp_p - lp_cur):
            theta = prop; lp_cur = lp_p; acc += 1
        samp[i] = theta
    return samp, acc / n_iter

n_chains, n_sample = 4, 3000
all_chains, all_acc = [], []
for ci in range(n_chains):
    init = [best_t[j] + np.random.RandomState(42 + ci * 7).normal(0, 0.08) for j in range(ndim)]
    samp, ar = run_chain(200 + ci * 19, init, learned_cov, n_sample)
    all_chains.append(samp)
    all_acc.append(ar)
    print(f"Chain {ci + 1}: accept={ar:.3f}")

chains = np.array(all_chains)  # (4, 3000, 4) in log-space

# Transform to natural scale: [beta, gamma, rho, phi]
chains_nat = np.zeros_like(chains)
for c in range(n_chains):
    chains_nat[c, :, 0] = np.exp(chains[c, :, 0]) * np.exp(chains[c, :, 1])  # beta = R0 * gamma
    chains_nat[c, :, 1] = np.exp(chains[c, :, 1])   # gamma
    chains_nat[c, :, 2] = np.exp(chains[c, :, 2])   # rho
    chains_nat[c, :, 3] = np.exp(chains[c, :, 3])   # phi

# Also store R0 separately
chains_R0 = np.exp(chains[:, :, 0:1])  # (4, 3000, 1)

np.savez(os.path.join(PROJECT_DIR, "outputs", "mcmc_results.npz"),
    chains=chains, chains_nat=chains_nat, chains_R0=chains_R0, acc=np.array(all_acc))
print("Step 2 done.")
