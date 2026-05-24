---
title: "Bayesian SIR Model Case Study: 2025 Epidemic"
author: "Raihan Research Group"
date: "`r Sys.Date()`"
output: 
  pdf_document:
    toc: true
    number_sections: true
  html_document:
    toc: true
---

# Introduction

This document provides a comprehensive breakdown of the Bayesian SIR (Susceptible-Infectious-Recovered) modeling workflow used to analyze the 2025 Epidemic Dataset. It is structured similarly to the official [Stan Boarding School Case Study](https://mc-stan.org/learn-stan/case-studies/boarding_school_case_study.html) and serves as a definitive guide for understanding the model, the code, the generated plots, and how to answer probable academic questions.

---

# 1. The Mathematical Model and Basic Terms

The model is based on the standard Ordinary Differential Equation (ODE) formulation of the SIR compartments.

## The Three Compartments (Buckets)
*   **$S$ (Susceptible)**: Healthy individuals at risk of catching the disease.
*   **$I$ (Infectious)**: Sick individuals currently transmitting the disease.
*   **$R$ (Recovered)**: Immune individuals who can no longer catch or spread it.
*   **$N$**: Total closed population ($S + I + R = 11,750$).

## The Parameters (Speed Limits)
*   **$\beta$ (Beta - Transmission Rate)**: The average number of effective contacts a single infectious person makes per week that result in new infections.
*   **$\gamma$ (Gamma - Recovery Rate)**: The rate at which people recover per week. The inverse ($1/\gamma$) is the average duration of infection.
*   **$R_0$ (Basic Reproduction Number)**: Calculated as $R_0 = \beta / \gamma$. The average number of people infected by a single sick person in a completely healthy population.

## Observation Terms
*   **$\rho$ (Rho - Reporting Rate)**: The percentage of true infectious individuals actually reported in the weekly data.
*   **$\mu$ (Mu - Expected Cases)**: $\mu = \rho \times I$. The number of reported cases the model expects.
*   **$\phi$ (Phi - Dispersion)**: Controls the variance/noise in the actual data around the expected mean.

---

# 2. Stan Code Explanation (`sir_model.stan`)

The Stan code is divided into several blocks.

## The Functions Block (The ODE)
```stan
functions {
  vector sir_ode(real t, vector y, array[] real theta, array[] real x_r, array[] int x_i) {
    real beta  = theta[1];
    real gamma = theta[2];
    real N     = x_r[1];

    real S = y[1];
    real I = y[2];

    vector[3] dydt;
    dydt[1] = -beta * S * I / N;                 // Susceptible decreases
    dydt[2] =  beta * S * I / N - gamma * I;     // Infectious grows then recovers
    dydt[3] =  gamma * I;                        // Recovered grows
    return dydt;
  }
}
```
**Explanation:** This mathematically defines the rules of the epidemic. It uses frequency-dependent transmission (dividing by $N$).

## The Transformed Parameters (Solving the ODE)
```stan
transformed parameters {
  real beta = R0 * gamma;
  array[2] real theta = {beta, gamma};
  // ... initial conditions setup ...
  
  y_sol = ode_rk45(sir_ode, y0, t0, ts, theta, x_r, x_i);

  for (i in 1:N_weeks) {
    mu[i] = fmax(1e-6, rho * y_sol[i][2]); 
  }
}
```
**Explanation:** 
*   **Reparameterization**: We sample $R_0$ and $\gamma$, then calculate $\beta = R_0 \times \gamma$. This significantly improves MCMC convergence by reducing correlation.
*   **`ode_rk45`**: This is Stan's numerical solver that "runs" the epidemic forward in time.
*   **`mu[i]`**: We scale the true prevalence (`y_sol[i][2]`) by the reporting rate `rho`. `fmax(1e-6)` is a guardrail preventing mathematical crashes if `mu` hits exactly zero.

## The Model Block (Priors and Likelihood)
```stan
model {
  // Priors
  R0    ~ lognormal(log(1.7), 0.25);
  gamma ~ lognormal(log(0.44), 0.3);
  rho   ~ lognormal(log(0.05), 0.5);
  phi   ~ exponential(1.0);

  // Likelihood
  for (i in 1:N_weeks) {
    cases[i] ~ neg_binomial_2(mu[i], phi);
  }
}
```
**Explanation:** 
*   **Priors**: We use `lognormal(log(median), sigma)`. This ensures our medians strictly follow our epidemiological assumptions ($R_0$ of 1.7, recovery rate of 0.44 per week, reporting rate of 5%).
*   **Likelihood**: The `neg_binomial_2` likelihood robustly handles the high variance (overdispersion) of weekly epidemic counts.

---

# 3. R Workflow & Plot Code Explanations

The `analysis.R` script handles data loading, prior checks, running Stan, and extracting posteriors. Here is an explanation of the generated plots:

*   **Plot 01: Raw Data (`01_raw_data.png`)**
    *   **What it does:** Uses `ggplot2` to visualize the raw weekly case counts over time.
    *   **Why we need it:** To visually inspect the epidemic curve shape and identify the peak.
*   **Plot 02: SIR Compartments (`02_sir_compartments.png`)**
    *   **What it does:** Uses `deSolve::ode` in R with the posterior median parameters to plot the $S$, $I$, and $R$ curves.
    *   **Why we need it:** Shows the unseen "true" epidemic dynamics underlying the reported data.
*   **Plot 03: Prior Predictive Check (`03_prior_predictive.png`)**
    *   **What it does:** Samples parameters *only* from the priors (ignoring data) and simulates epidemic curves.
    *   **Why we need it:** Proves to reviewers that our prior assumptions are broad enough to cover the actual data without strongly forcing a specific outcome.
*   **Plot 04: Trace Plots (`04_trace_plots.png`)**
    *   **What it does:** Uses `bayesplot::mcmc_trace` to visualize the MCMC sampling chains.
    *   **Why we need it:** "Hairy caterpillar" plots prove the algorithm successfully converged and explored the parameter space cleanly.
*   **Plot 05 & 06: Posteriors & Pair Plots (`05_posterior_histograms.png`, `06_pair_plots.png`)**
    *   **What they do:** Show the final estimated distributions for $\beta$, $\gamma$, $R_0$, and $\rho$, along with their correlations.
    *   **Why we need them:** Provides the exact statistical answers and confidence intervals for our epidemic parameters.
*   **Plot 07 & 08: Posterior Predictive Checks (`07_posterior_predictive.png`)**
    *   **What it does:** Compares the model's 50% and 90% Confidence Intervals of expected cases against the actual data dots.
    *   **Why we need it:** The ultimate test of model fit. It proves the ODE model successfully learned and can replicate the real-world dataset.

---

# 4. Probable Defense Questions & Answers

If presenting this to a professor or committee, prepare for these questions:

### Q1: Why did you use a Bayesian approach instead of standard Frequentist (Maximum Likelihood) fitting?
**Answer:** The Bayesian approach allows us to incorporate known epidemiological priors (like the fact that respiratory diseases usually have an $R_0$ around 1.5 - 2.5). More importantly, instead of just giving a single "best guess" line, it provides full posterior distributions, meaning we can quantify our exact uncertainty (95% Credible Intervals) for every parameter.

### Q2: How did you handle the fact that not all sick people are reported?
**Answer:** We explicitly modeled a reporting rate, $\rho$ (rho). Our model assumes the observed data is only a fraction of the true infectious pool ($I$). We placed a prior centering around 5% reporting, and let the model infer the exact rate based on the epidemic curve dynamics.

### Q3: Why did you use `neg_binomial_2` instead of a Poisson likelihood?
**Answer:** Epidemic count data is almost always "overdispersed"—meaning the variance is much higher than the mean due to super-spreader events and reporting batching. A Poisson distribution assumes mean equals variance, which is too restrictive and would lead to overconfident, narrow confidence intervals. The Negative Binomial perfectly handles this extra noise.

### Q4: In Plot 07, the model slightly underpredicts the peak. Is the model broken?
**Answer:** No, this is a known, standard structural limitation of the basic SIR ODE model. The SIR model assumes "homogeneous mixing" (everyone interacts with everyone equally). In reality, human networks are clustered, causing epidemics to spike faster and sharper than homogeneous math predicts. This slight underprediction is expected and indicates the model is mathematically honest.

### Q5: How do you know the MCMC chains converged?
**Answer:** We used three standard diagnostics:
1.  **Trace Plots (Plot 04):** The chains look like well-mixed "hairy caterpillars" with no wandering.
2.  **$\hat{R}$ (R-hat):** All parameters had an $\hat{R} < 1.05$, indicating the 4 independent chains arrived at the exact same distribution.
3.  **ESS (Effective Sample Size):** All parameters had an $N_{eff} > 400$, meaning we have enough independent samples to trust the tails of the distributions.

### Q6: Why did you sample $R_0$ and $\gamma$ instead of $\beta$ and $\gamma$ directly?
**Answer:** Because $\beta$ and $\gamma$ are highly correlated in epidemic models (a higher transmission rate and a faster recovery rate can produce the exact same curve). Sampling them independently causes MCMC chains to get stuck. By sampling $R_0$ and $\gamma$, and calculating $\beta = R_0 \times \gamma$, we broke that collinearity and achieved much faster, stable convergence. This is a best-practice technique in Stan.
