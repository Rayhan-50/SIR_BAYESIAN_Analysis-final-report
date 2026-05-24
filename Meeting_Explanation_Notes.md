# Meeting Explanation Notes — Bayesian SIR Model
**Raihan Research Group | For next meeting with Prof. Abu Muhammad Hossain**

---

## 1. What is the SIR Model?

The SIR model divides the total population **N** into three groups at each point in time:

| Compartment | Symbol | Meaning |
|-------------|--------|---------|
| Susceptible | S(t)   | People who can still catch the disease |
| Infectious  | I(t)   | People currently infected and spreading it |
| Recovered   | R(t)   | People who recovered (assumed immune) |

The model says:

```
dS/dt = -β · S(t) · I(t) / N      ← susceptible people becoming infected
dI/dt =  β · S(t) · I(t) / N  - γ · I(t)   ← new infections minus recoveries
dR/dt =  γ · I(t)              ← people recovering
```

---

## 2. What are the Parameters and Where Do They Come From?

> **KEY POINT for the professor:**
> We did NOT manually pick or "change" these numbers.
> They are **estimated from the data** using Bayesian statistics.

| Parameter | Meaning | How it was set |
|-----------|---------|----------------|
| **β (beta)** | Transmission rate — how fast the disease spreads | **Estimated by MCMC from data** |
| **γ (gamma)** | Recovery rate — how fast people recover | **Estimated by MCMC from data** |
| **R₀** | Basic reproduction number = β / γ (average people one case infects) | **Derived: R₀ = β / γ** |
| **ρ (rho)** | Reporting rate — fraction of true cases actually detected/reported | **Estimated by MCMC from data** |
| **φ (phi)** | Overdispersion — controls noise/variability in counts | **Estimated by MCMC from data** |

---

## 3. How Does Bayesian Estimation Work? (Step by Step)

### Step 1 — Prior Beliefs (before seeing data)
We start by saying: *"Based on flu biology, what are plausible ranges for each parameter?"*

```
R₀    ~ LogNormal(log(1.7), 0.25)   → we expect R₀ around 1.7, but allow 1.2–2.5
γ     ~ LogNormal(log(0.44), 0.3)   → recovery takes ~2 weeks (1/0.44 ≈ 2.3 weeks)
ρ     ~ LogNormal(log(0.05), 0.5)   → we expect ~5% of cases are reported
φ     ~ Exponential(1)              → weak prior on overdispersion
```

### Step 2 — Likelihood (comparing model to data)
For each candidate set of parameters, we compute:
- "If β and γ were these values, what case curve would the SIR model predict?"
- "How well does that predicted curve match the actual observed weekly cases?"

```
Observed_cases(t) ~ NegativeBinomial( ρ · I(t),  φ )
```
*The Negative Binomial distribution accounts for week-to-week random variation.*

### Step 3 — MCMC Sampling (Stan does this automatically)
Stan's **Hamiltonian Monte Carlo** sampler:
- Explores thousands of candidate parameter combinations
- Keeps combinations that fit the data well
- Discards combinations that fit poorly
- The final collection of kept samples = **posterior distribution**

### Step 4 — Posterior = Updated Beliefs
After seeing the data, we get **distributions** for each parameter:

| Parameter | Prior expectation | Posterior median (aggregated site) | 95% Credible Interval |
|-----------|-------------------|------------------------------------|----------------------|
| R₀        | ~1.7              | **1.44**                           | [1.31, 1.65]         |
| β         | ~0.75             | **0.95**                           | [0.80, 1.14]         |
| γ         | ~0.44             | **0.66**                           | [0.49, 0.86]         |
| ρ         | ~0.05             | **0.079**                          | [0.049, 0.131]       |

---

## 4. How to Explain "How Did You Change the Transmission Rate?"

**Say this to the professor:**

> *"We did not manually change beta. Instead, we gave the model a prior belief
> that beta is probably around 0.75 (based on flu literature), and then Stan's
> MCMC algorithm adjusted it by fitting the observed weekly case counts.
> The data told us beta is actually closer to 0.95. The Bayesian framework
> automatically finds the value of beta—and all other parameters—that best
> explains what we observed, while staying consistent with biological plausibility."*

---

## 5. What the Model Output Means

### R₀ = 1.44
- Each infected person infects **1.44 others** on average
- Since R₀ > 1 → epidemic can grow (confirmed by the rising case curve)
- Since R₀ < 2 → a moderate, controllable epidemic

### ρ = 0.079 (8%)
- Only **~8% of true infections** were officially reported
- This is realistic for influenza surveillance (under-reporting is common)

### R² = 0.84
- The model explains **84% of the variance** in observed weekly case counts
- This is a good fit for an epidemic model

---

## 6. What Still Needs to Be Done (Before Final Report)

- [x] Aggregated analysis done (one combined dataset)
- [ ] **Per-site analysis: ILI sentinel site**
- [ ] **Per-site analysis: Severe/Hospital sentinel site**
- [ ] **Per-site analysis: SARI sentinel site**
- [ ] Compare R₀ across all three sites
- [ ] Professor review and approval of all results
- [ ] Final report writing

---

## 7. Useful Terms to Know for the Meeting

| Term | Plain English |
|------|---------------|
| MCMC | A computer algorithm that explores parameter space randomly but intelligently |
| Posterior distribution | Our updated beliefs about a parameter AFTER seeing the data |
| Credible interval | The range where the true parameter value probably lies (95% CI = 95% probability) |
| R-hat (R̂) | Convergence check — should be < 1.05 (means all 4 chains agree) |
| n_eff | Effective sample size — should be > 400 (means we have enough samples) |
| Negative Binomial | A probability distribution for count data with extra variability |
| rho (ρ) | Under-reporting factor — not all true cases are diagnosed and reported |
