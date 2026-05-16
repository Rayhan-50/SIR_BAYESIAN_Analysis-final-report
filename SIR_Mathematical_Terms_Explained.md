# Mathematical Terms Explained — Bayesian SIR Model
## Plain-Language Guide for the 2025 Epidemic Analysis

---

> [!NOTE]
> Every symbol below is explained in simple words first, then with its formal meaning. No advanced math background is needed to follow this guide.

---

## Part 1 — The SIR Model: What Are S, I, R?

The **SIR model** divides the entire population into three groups (called **compartments**). At any point in time, every person belongs to exactly one group.

```
  S  ──────────►  I  ──────────►  R
(Susceptible)   (Infectious)   (Recovered)
```

| Symbol | Full name | Plain meaning |
|---|---|---|
| **S** | Susceptible | People who have **not yet caught** the disease. They are at risk of being infected. |
| **I** | Infectious | People who are **currently sick** and can spread the disease to others. |
| **R** | Recovered | People who have **recovered** (or died). They are immune and no longer spread the disease. |
| **N** | Total population | The total number of people: N = S + I + R. In this study, **N = 11,750**. |

> **Example in plain words:** On Day 1 of the epidemic, almost everyone is Susceptible (S ≈ 11,749), only 1 person is Infectious (I = 1), and nobody has Recovered yet (R = 0).

---

## Part 2 — The Greek Letters (Model Parameters)

These are the "knobs" of the model — numbers the model learns from the data.

---

### β (Beta) — Transmission Rate

> **Plain meaning:** How fast does the disease spread from sick people to healthy people?

- **Technical definition:** The average number of **new infections** one infectious person causes per susceptible person per week.
- **In this study:** β median ≈ **0.95** per week
- **Intuition:** A higher β → disease spreads faster → epidemic grows more steeply. A β of ~0.95 means each infectious person "contacts" susceptibles at rate 0.95/week.
- **In the ODE equation:** New infections per week = β × S × I / N

> **Analogy:** Think of β as how "contagious" a person is — how hard they sneeze and how close they stand to others.

---

### γ (Gamma) — Recovery Rate

> **Plain meaning:** How quickly do sick people get better (or leave the infectious pool)?

- **Technical definition:** The fraction of infectious people who recover each week. It is the **rate of leaving** the I compartment.
- **In this study:** γ median ≈ **0.66** per week
- **Intuition:** If γ = 0.66, then about 66% of sick people recover each week.
- **Mean infectious period** = 1/γ = 1/0.66 ≈ **1.5 weeks** — an average sick person is contagious for about 1.5 weeks.
- **In the ODE equation:** Recoveries per week = γ × I

> **Analogy:** Think of γ as how quickly sick people "leave the stage." A high γ means people recover fast; a low γ means they stay sick a long time.

---

### ρ (Rho) — Reporting Fraction

> **Plain meaning:** What fraction of true cases actually get reported/detected?

- **Technical definition:** The probability that a truly infectious person is counted in the official case data.
- **In this study:** ρ median ≈ **0.079** → only about **8% of true cases** are reported!
- **Why this matters:** If the model says I(t) = 1,000 infectious people at some week, but ρ = 0.08, then the expected reported cases = ρ × I(t) = 80. The remaining 920 cases are real but undetected (mild, asymptomatic, or untested).
- **In the model:** Expected reported cases = ρ × I(t) → this is called **μ (mu)**

> **Analogy:** Think of a net with holes. The true epidemic is the fish in the sea, ρ is how much of the net actually catches fish. A small ρ (8%) means the net is very leaky — most cases slip through uncounted.

---

### φ (Phi) — Overdispersion Parameter

> **Plain meaning:** How much extra randomness/variability is there in the case counts, beyond what a simple model would predict?

- **Technical definition:** The "size" parameter of the **Negative Binomial distribution** (the statistical model for case counts). A higher φ means counts are **closer to a Poisson distribution** (less extra noise). A lower φ means counts are more **over-dispersed** (more bursty and variable).
- **In this study:** φ median ≈ **4.4** → moderate overdispersion
- **Why this matters:** Real epidemic data has more variability than a simple Poisson model allows (some weeks are unusually high due to outbreaks, testing spikes, etc.). The Negative Binomial with φ accounts for this.
- **Variance formula:** Var(cases) = μ + μ²/φ — as φ → ∞, this approaches Poisson; as φ → small, variance is much larger than the mean.

> **Analogy:** Think of weekly case counts as rainfall measurements. Some weeks it rains way more than average (a storm), some much less. φ controls how "stormy" the data is allowed to be.

---

### μ (Mu) — Expected Cases

> **Plain meaning:** The model's prediction of how many cases we expect to see in a given week.

- **Formula:** μ(t) = ρ × I(t)
- This is the **mean** of the Negative Binomial distribution at each time step.
- It combines the epidemic dynamics (how many people are truly infectious) with the observation process (what fraction gets reported).

---

## Part 3 — R₀ (R-naught): The Most Important Number

> **Plain meaning:** On average, how many new people does ONE sick person infect in a completely susceptible population?

- **Formula:** R₀ = β / γ
- **In this study:** R₀ median = **1.44** (95% credible interval: 1.31 to 1.65)

| R₀ value | What it means |
|---|---|
| R₀ < 1 | Disease dies out on its own — not enough transmission to sustain an epidemic |
| R₀ = 1 | Disease stays at a constant level (endemic equilibrium) |
| R₀ > 1 | Disease **spreads exponentially** → epidemic! |

- In this study, R₀ ≈ 1.44 means each infected person infects ~1.44 others on average → confirmed epidemic.
- **Why R₀ = β/γ?** Because β is "how fast you infect others" and γ is "how fast you stop being infectious." The ratio tells you the net number of infections per person over their entire infectious period (1/γ weeks).

> **Famous examples:** COVID-19 original strain R₀ ≈ 2–3; Measles R₀ ≈ 12–18; Seasonal flu R₀ ≈ 1.2–1.4. This epidemic's R₀ ≈ 1.44 is close to **seasonal influenza**.

---

## Part 4 — The ODE System (Differential Equations)

The SIR model is defined by three **ordinary differential equations (ODEs)** — these are equations that describe the **rate of change** of each compartment over time.

```
dS/dt = - β·S·I / N          (susceptibles leaving → becoming infected)

dI/dt = + β·S·I / N - γ·I    (new infections arriving − recoveries leaving)

dR/dt = + γ·I                (recoveries accumulating)
```

| Term | Meaning |
|---|---|
| dS/dt | Rate of change of Susceptible people per week |
| dI/dt | Rate of change of Infectious people per week |
| dR/dt | Rate of change of Recovered people per week |
| β·S·I/N | New infections per week (mass-action contact) |
| γ·I | Recoveries per week |

> **Plain reading of dI/dt:** The number of infectious people changes each week by gaining new infections (β·S·I/N) and losing recoveries (γ·I). When β·S·I/N > γ·I, I grows (epidemic rising). When β·S·I/N < γ·I, I falls (epidemic declining).

---

## Part 5 — Bayesian Statistical Terms

---

### Prior (Prior Distribution)

> **Plain meaning:** What do we believe about a parameter **before** seeing the data?

- A prior is a probability distribution that captures existing knowledge or assumptions.
- **Example:** We believe R₀ is probably around 1.7 (informed by similar diseases), so we set R₀ ~ LogNormal(log(1.7), 0.25). This says "R₀ is most likely near 1.7, but could plausibly range from ~1.1 to ~2.8."
- A **weakly informative prior** allows a wide range of values — it guides the model but doesn't force it.
- A **non-informative prior** says "any value is equally plausible" — rarely used in practice.

---

### Likelihood

> **Plain meaning:** Given a specific set of parameter values, how probable is the observed data?

- **Formula used here:** cases(t) ~ Negative-Binomial(μ(t), φ)
- This says: at week t, the observed case count comes from a Negative Binomial distribution with mean μ(t) and dispersion φ.
- The likelihood "scores" each possible set of (β, γ, ρ, φ) values by how well they explain the actual observed case counts.

---

### Posterior (Posterior Distribution)

> **Plain meaning:** What do we believe about a parameter **after** seeing the data? It is the combination of prior beliefs + evidence from the data.

- **Bayes' Theorem:** Posterior ∝ Likelihood × Prior
- The posterior is the final answer — it tells us the full probability distribution of each parameter given the observed epidemic data.
- In this study, all our estimates (mean, median, credible intervals) come from the posterior.

---

### Credible Interval (CrI) / Credible Region

> **Plain meaning:** A range that contains the true parameter value with a stated probability.

- A **95% credible interval** means: "We are 95% certain the true parameter lies within this range."
- **Example:** R₀ 95% CrI = [1.31, 1.65] means we are 95% sure the true R₀ is between 1.31 and 1.65.
- ⚠️ This is **different** from a frequentist confidence interval! A CrI is a direct probability statement about the parameter.

---

## Part 6 — MCMC Diagnostics

### MCMC (Markov Chain Monte Carlo)

> **Plain meaning:** A computer algorithm that **explores** the posterior distribution by taking millions of random steps through the parameter space.

- Because the posterior has no simple formula, we use MCMC to draw thousands of samples from it.
- These samples approximate the posterior — from them we compute means, medians, and credible intervals.
- **Stan** (used in this study) uses a specific MCMC method called HMC (Hamiltonian Monte Carlo) which is particularly efficient.

---

### Chain

> **Plain meaning:** One independent "run" of the MCMC algorithm.

- In this study, 4 chains were run in parallel (each starting from a different random point).
- If all 4 chains converge to the same distribution → the algorithm has successfully found the posterior.

---

### R̂ (R-hat) — Gelman-Rubin Convergence Statistic

> **Plain meaning:** Did all the MCMC chains agree with each other? (Did the algorithm converge?)

- **Formula:** R̂ compares the variance **within** each chain to the variance **between** chains.
- **R̂ = 1.0** → perfect convergence (chains are identical in distribution)
- **R̂ < 1.05** → acceptable convergence ✅
- **R̂ ≥ 1.05** → chains have not converged ❌ — results are unreliable

> **In this study:** All R̂ values are between 1.000 and 1.005 — **excellent convergence.**

---

### ESS (Effective Sample Size)

> **Plain meaning:** How many truly independent samples does our MCMC chain give us?

- Because consecutive MCMC samples are correlated (each step is close to the previous), the actual number of **independent** samples is less than the total number of iterations.
- ESS tells us the equivalent number of independent samples.
- **ESS > 400** is generally considered acceptable for reliable inference.

> **In this study:** ESS ranges from 946 to 1510 for all parameters — **well above the threshold.**

---

### Trace Plot

> **Plain meaning:** A time-series plot showing the value of a parameter at each MCMC iteration. Used to visually check convergence.

- A **good trace plot** looks like a "hairy caterpillar" — all chains mixed together, stationary around a stable value, no drifting or getting stuck.
- A **bad trace plot** shows chains drifting apart, getting stuck in one region, or showing trends.

---

## Part 7 — Model Fit Statistics

### R² (R-squared) — Coefficient of Determination

> **Plain meaning:** What fraction of the variability in the data is explained by the model?

- **Formula:** R² = 1 − (Sum of squared residuals) / (Total sum of squares)
- **Range:** 0 to 1 (higher is better)
- R² = 1.0 → perfect fit | R² = 0 → model explains nothing

| R² | Interpretation |
|---|---|
| > 0.90 | Excellent fit |
| 0.75 – 0.90 | Good fit |
| 0.50 – 0.75 | Moderate fit |
| < 0.50 | Poor fit |

> **In this study:** R² = **0.84** → **Good fit** ✅

---

### RMSE (Root Mean Squared Error)

> **Plain meaning:** On average, how far off is the model's prediction from the observed case count (in the same units as the data)?

- **Formula:** RMSE = √( mean of (observed − predicted)² )
- Lower RMSE = better fit
- RMSE is in **the same units as the cases** (cases per week)

> **In this study:** RMSE = **8.63 cases/week** — the model's median prediction is off by ~8–9 cases per week on average. Given that the peak is 69, this is an error of ~12–13%, which is acceptable.

---

### Residual

> **Plain meaning:** The difference between what was actually observed and what the model predicted.

- **Formula:** Residual = Observed cases − Model predicted (median)
- A **positive residual** → model underpredicted (reality was higher than forecast)
- A **negative residual** → model overpredicted (reality was lower than forecast)
- **Ideal residuals** are small, centered on zero, and have no pattern over time.
- **Systematic patterns** in residuals indicate the model is missing something structurally.

---

### Negative Binomial Distribution

> **Plain meaning:** A flexible probability distribution for **count data** (like weekly case counts) that allows for extra variability beyond what Poisson allows.

- Used instead of Poisson because real epidemic data is "clumpy" — some weeks have unusually many or few cases due to testing variability, reporting delays, cluster outbreaks, etc.
- **Parameters:** mean (μ) and dispersion (φ)
- As φ → ∞, Negative Binomial → Poisson (no overdispersion)
- As φ → small, variance ≫ mean (lots of overdispersion)

---

### LogNormal Distribution

> **Plain meaning:** A distribution for **positive** numbers where the logarithm of the number follows a Normal (bell-curve) distribution.

- Used as priors for β, γ, ρ because these are positive rates/fractions.
- **Notation:** X ~ LogNormal(log(m), σ) means the median of X is m.
- **Example:** R₀ ~ LogNormal(log(1.7), 0.25) → median R₀ ≈ 1.7, with 95% prior range roughly [1.05, 2.75]

---

## Quick Reference Card

| Symbol | Name | Plain meaning | Value in this study |
|---|---|---|---|
| **N** | Population size | Total people | 11,750 |
| **S** | Susceptible | Never infected, at risk | Starts at 11,749 |
| **I** | Infectious | Currently sick & spreading | Peaks at ~week 27 |
| **R** | Recovered | No longer infectious | Accumulates over time |
| **β** | Transmission rate | Spreading speed | 0.947/week |
| **γ** | Recovery rate | Healing speed | 0.657/week |
| **1/γ** | Mean infectious period | How long people are sick | ~1.5 weeks |
| **ρ** | Reporting fraction | % of cases detected | ~8% |
| **φ** | Overdispersion | Extra randomness in counts | 4.4 |
| **μ** | Expected cases | Model prediction each week | ρ × I(t) |
| **R₀** | Basic reproduction number | How many people one person infects | **1.44** |
| **R̂** | R-hat | Did MCMC chains converge? | 1.000–1.005 ✅ |
| **ESS** | Effective sample size | How many independent samples | 946–1510 ✅ |
| **R²** | R-squared | Overall model fit quality | **0.84** ✅ |
| **RMSE** | Root mean squared error | Average prediction error | 8.63 cases/week |
| **CrI** | Credible interval | Range containing true value with X% probability | e.g., R₀ ∈ [1.31, 1.65] |
| **Prior** | Prior distribution | Belief before seeing data | e.g., R₀ ~ LogNormal |
| **Posterior** | Posterior distribution | Updated belief after seeing data | Final estimates |
| **MCMC** | Markov Chain Monte Carlo | Algorithm to explore the posterior | Stan/HMC used |
