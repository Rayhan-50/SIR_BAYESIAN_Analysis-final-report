# Bayesian SIR Model Analysis — 2025 Epidemic Dataset

**Raihan Research Group | May 2026**

This repository contains a fully Bayesian analysis of 52 weeks of epidemic case data using a stochastic SIR compartmental model implemented in Stan.

## Overview

We fit a four-parameter SIR model to weekly case counts from a 2025 epidemic (N = 11,750, total cases = 833). The model uses Hamiltonian Monte Carlo sampling via RStan and a negative-binomial observation model with an explicit reporting rate to account for under-ascertainment.

**Key results:**

| Parameter | Median | 95% CI |
|-----------|--------|--------|
| R0        | 1.44   | [1.31, 1.65] |












| beta      | 0.95   | [0.80, 1.14] |
| gamma     | 0.66   | [0.49, 0.86] |
| rho       | 0.079  | [0.049, 0.131] |
| Model R²  | 0.84   | — |

## Repository Structure

```
├── README.md
├── data/
│   └── clean_epidemic_dataset_2025.csv   # 52-week case counts
├── code/
│   ├── sir_model.stan                    # Stan model with ODE solver
│   ├── analysis.R                        # Full Bayesian workflow in R
│   └── build_report.R                    # Word document report builder
├── plots/
│   ├── 01_raw_data.png
│   ├── 02_sir_compartments.png
│   ├── 03_prior_predictive.png
│   ├── 04_trace_plots.png
│   ├── 05_posterior_histograms.png
│   ├── 06_pair_plots.png
│   ├── 07_posterior_predictive.png
│   ├── 08_posterior_predictive_sim.png
│   ├── 09_R0_distribution.png
│   └── 10_residuals.png
├── outputs/
│   └── results_R.json                    # Posterior summary statistics
└── SIR_Bayesian_Analysis_Report.docx     # Final report
```

## How to Run

### Prerequisites

- R >= 4.4
- RStan >= 2.32 (`install.packages("rstan", repos = "https://stan-dev.r-universe.dev")`)
- Required R packages: `deSolve`, `ggplot2`, `bayesplot`, `gridExtra`, `jsonlite`, `officer`, `dplyr`, `tidyr`, `GGally`

### Steps

```r
# Install packages (first time only)
install.packages(c("deSolve", "bayesplot", "officer"),
                 repos = "https://cran.rstudio.com/")
install.packages("rstan", repos = "https://stan-dev.r-universe.dev")

# Run full analysis (MCMC sampling, all 10 plots, results JSON)
Rscript code/analysis.R

# Build Word report
Rscript code/build_report.R
```

## Model Specification

The SIR ODE system:

```
dS/dt = -beta * S * I / N
dI/dt =  beta * S * I / N - gamma * I
dR/dt =  gamma * I
```

with `S(0) = N-1`, `I(0) = 1`, `R(0) = 0`, and `N = 11,750`.

Observation model: `cases(t) ~ NegBin(rho * I(t), phi)`

Stan's `ode_rk45` integrator is used to solve the ODE system within the model block.

## References

- Carpenter et al. (2017). Stan: A probabilistic programming language. *Journal of Statistical Software*, 76(1).
- Gabry et al. (2019). Visualization in Bayesian workflow. *JRSS-A*, 182(2), 389–402.
- Kermack & McKendrick (1927). A contribution to the mathematical theory of epidemics. *Proc. Royal Soc. A*, 115, 700–721.
# SIR_Bayesian_Analysis_Reporting
