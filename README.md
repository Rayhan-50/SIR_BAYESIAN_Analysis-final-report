# Bayesian SIR Model Analysis — 2025 Epidemic Dataset
**Raihan Research Group | May 2026**

This repository contains a fully Bayesian analysis of 52 weeks of epidemic case data using a stochastic SIR compartmental model implemented in Stan — run **separately for each sentinel surveillance site**.

## Sentinel Sites

| Site | Description | Population (N) |
|------|-------------|---------------|
| **ILI** | Influenza-Like Illness (Outpatient Sentinel) | 11,750 |
| **Severe** | Severe/Hospital Sentinel | 11,750 |
| **SARI** | Severe Acute Respiratory Infection | 500,000 |

## Aggregated Site Results (Baseline)

| Parameter | Median | 95% CI |
|-----------|--------|--------|
| R₀        | 1.44   | [1.31, 1.65] |
| beta      | 0.95   | [0.80, 1.14] |
| gamma     | 0.66   | [0.49, 0.86] |
| rho       | 0.079  | [0.049, 0.131] |
| Model R²  | 0.84   | — |

> Per-site results are in `outputs/results_ILI.json`, `outputs/results_Severe.json`, `outputs/results_SARI.json` after running the analysis.

## Repository Structure

```
├── README.md
├── Meeting_Explanation_Notes.md         # Notes for professor meeting
├── data/
│   ├── clean_epidemic_dataset_2025.csv  # Aggregated 52-week case counts
│   ├── unified_weekly_dataset.xlsx      # Raw multi-site sentinel data
│   ├── site_ILI.csv                     # ILI sentinel weekly cases
│   ├── site_Severe.csv                  # Severe/Hospital sentinel weekly cases
│   └── site_SARI.csv                    # SARI sentinel weekly cases
├── code/
│   ├── sir_model.stan                   # Stan model with ODE solver
│   ├── analysis.R                       # Aggregated Bayesian workflow
│   └── analysis_per_site.R              # Per-site loop (ILI, Severe, SARI)
├── plots/
│   ├── 01_raw_data.png ... 10_residuals.png  (aggregated)
│   ├── ILI/                             # 10 plots for ILI site
│   ├── Severe/                          # 10 plots for Severe site
│   └── SARI/                            # 10 plots for SARI site
└── outputs/
    ├── results_R.json                   # Aggregated posterior summary
    ├── results_ILI.json                 # ILI site posterior summary
    ├── results_Severe.json              # Severe site posterior summary
    └── results_SARI.json                # SARI site posterior summary
```

## How to Run

### Prerequisites

- R >= 4.4
- RStan >= 2.32
- Required packages: `deSolve`, `ggplot2`, `bayesplot`, `gridExtra`, `jsonlite`, `dplyr`, `tidyr`, `GGally`

```r
install.packages(c("deSolve","bayesplot","gridExtra","GGally"),
                 repos = "https://cran.rstudio.com/")
install.packages("rstan", repos = "https://stan-dev.r-universe.dev")
```

### Steps

```r
# Aggregated analysis (original single-site)
Rscript code/analysis.R

# Per-site analysis — ILI, Severe, SARI (professor's request)
Rscript code/analysis_per_site.R
```

## Model Specification

### SIR ODE System

```
dS/dt = -β · S(t) · I(t) / N
dI/dt =  β · S(t) · I(t) / N  - γ · I(t)
dR/dt =  γ · I(t)
```

Initial conditions: `S(0) = N-1`, `I(0) = 1`, `R(0) = 0`

### Parameters (all estimated by Bayesian MCMC — not manually set)

| Parameter | Prior | Interpretation |
|-----------|-------|---------------|
| R₀        | LogNormal(log(1.7), 0.25) | Basic reproduction number |
| γ (gamma) | LogNormal(log(0.44), 0.3) | Recovery rate (~2 week infectious period) |
| ρ (rho)   | LogNormal(log(0.05), 0.5) | Case reporting rate |
| φ (phi)   | Exponential(1)            | Negative-binomial overdispersion |
| β (beta)  | Derived: β = R₀ × γ      | Transmission rate |

### Observation Model

```
cases(t) ~ NegBinomial( ρ · I(t),  φ )
```

Stan's `ode_rk45` integrator solves the ODE within each MCMC iteration.

## References

- Carpenter et al. (2017). Stan: A probabilistic programming language. *Journal of Statistical Software*, 76(1).
- Gabry et al. (2019). Visualization in Bayesian workflow. *JRSS-A*, 182(2), 389–402.
- Kermack & McKendrick (1927). A contribution to the mathematical theory of epidemics. *Proc. Royal Soc. A*, 115, 700–721.
