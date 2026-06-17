# Bayesian SIR Model Analysis — 2025 Epidemic Dataset
**Raihan Research Group | May 2026**

This repository contains a fully Bayesian analysis of 52 weeks of epidemic surveillance data using a stochastic SIR compartmental model implemented in Stan — run independently for each surveillance stream.

## Surveillance Streams

| Stream | Description | Population (N) | Model |
|--------|-------------|----------------|-------|
| **Community** | Primary community surveillance (independent stream) | 11,750 | Standard SIR |
| **ILI** | Influenza-Like Illness (Outpatient Sentinel) | 11,750 | Standard SIR |
| **Severe** | Severe/Hospital Sentinel | 11,750 | Standard SIR |

> **Important Data Clarification:**
> The Community surveillance dataset (`clean_epidemic_dataset_2025.csv`) is an **independent** reporting system. It is **NOT** the arithmetic sum of the three sentinel sites. ILI and Severe share the same catchment (N=11,750), so summing their cases would double-count. Each surveillance stream is analysed independently.

## Community Surveillance Results (Primary)

| Parameter | Median | 95% CrI |
|-----------|--------|---------|
| R₀ | 1.44 | [1.31, 1.65] |
| beta | 0.95 | [0.80, 1.14] |
| gamma | 0.66 | [0.49, 0.86] |
| rho | 0.079 | [0.049, 0.131] |
| Model R² | 0.84 | — |

> Per-site results are in `outputs/results_ILI.json`, `outputs/results_Severe.json` after running the analysis. Data provenance is in `outputs/data_provenance.json`.

## Known Data Quality Issues (Manuscript Limitations)

1. **Stair-step artifact (ILI & Severe):** Identical case counts appear for 3–4 consecutive weeks (e.g., exactly 27 cases × 4 weeks). This likely reflects monthly data interpolated to weekly by dividing by 4. The ODE model cannot perfectly fit step-function data; structured residuals will appear. This is disclosed in the manuscript Limitations section.

## Repository Structure

```
├── README.md
├── data/
│   ├── clean_epidemic_dataset_2025.csv  # Community surveillance (PRIMARY — independent stream)
│   ├── unified_weekly_dataset.xlsx      # Raw multi-site sentinel data
│   ├── site_ILI.csv                     # ILI sentinel weekly cases
│   ├── site_Severe.csv                  # Severe/Hospital sentinel weekly cases
├── code/
│   ├── prepare_data.R                   # [RUN FIRST] Data validation & provenance
│   ├── sir_model.stan                   # Standard SIR Stan model (community, ILI, Severe)
│   ├── analysis.R                       # Community surveillance Bayesian workflow
│   ├── analysis_per_site.R              # Per-site loop (ILI, Severe)
│   ├── loo_sensitivity.R                # LOO-CV + prior sensitivity analysis
│   └── plot_fitted.R                    # Publication-quality fitted plots (N auto-read from JSON)
├── plots/
│   ├── 01_raw_data.png ... 10_residuals.png  (community surveillance)
│   ├── fitted_model_presentation.png
│   ├── fitted_model_with_residuals.png
│   ├── ILI/                             # 10 plots for ILI site
│   └── Severe/                          # 10 plots for Severe site
└── outputs/
    ├── data_provenance.json             # Dataset validation report
    ├── results_R.json                   # Community surveillance posterior summary
    ├── results_ILI.json                 # ILI site posterior summary
    └── results_Severe.json              # Severe site posterior summary
```

## How to Run (Full Reproducible Pipeline)

### Prerequisites

- R >= 4.4
- RStan >= 2.32
- Required packages: `deSolve`, `ggplot2`, `bayesplot`, `gridExtra`, `jsonlite`, `dplyr`, `tidyr`, `GGally`, `loo`

```r
install.packages(c("deSolve","bayesplot","gridExtra","GGally","loo"),
                 repos = "https://cran.rstudio.com/")
install.packages("rstan", repos = "https://stan-dev.r-universe.dev")
```

### Steps (in order)

```r
# Step 0 — Validate data & generate provenance report
Rscript code/prepare_data.R

# Step 1 — Community surveillance analysis (primary Bayesian SIR)
Rscript code/analysis.R

# Step 2 — Generate publication-quality fitted plots (reads N from JSON automatically)
Rscript code/plot_fitted.R

# Step 3 — Per-site analysis: ILI, Severe (standard SIR)
Rscript code/analysis_per_site.R

# Step 4 — LOO cross-validation + prior sensitivity analysis
Rscript code/loo_sensitivity.R
```

## Model Specification

### Standard SIR ODE System

```
dS/dt = -β · S(t) · I(t) / N
dI/dt =  β · S(t) · I(t) / N  - γ · I(t)
dR/dt =  γ · I(t)
```

Initial conditions: `S(0) = N-1`, `I(0) = 1`, `R(0) = 0`


### Parameters (all estimated by Bayesian MCMC — not manually set)

| Parameter | Prior | Interpretation |
|-----------|-------|---------------|
| R₀ | LogNormal(log(2.0), 0.3) | Basic reproduction number; median 2.0, 95% in [1.1, 3.5] |
| γ (gamma) | LogNormal(log(0.5), 0.3) | Recovery rate; median ~2-week infectious period |
| ρ (rho) | Beta(2, 10) | Case reporting rate; mean ~17%, constrained to (0, 1) |
| φ (phi) | Exponential(0.5) | Negative-binomial overdispersion; mean 2 |
| I₀ | LogNormal(log(1), 2.0) | Initial infectious seed; wide prior |
| λ (lambda) | Exponential(0.5) | Background sporadic case rate |
| β (beta) | Derived: β = R₀ × γ | Transmission rate |

### Observation Model

```
incidence(t)  =  S(t-1) - S(t)                   [new infections per week — a flow]
mu(t)         =  ρ · incidence(t)  +  λ           [expected reported cases]
cases(t)      ~  NegBinomial_2( mu(t),  φ )
```

## References

- Carpenter et al. (2017). Stan: A probabilistic programming language. *Journal of Statistical Software*, 76(1).
- Gabry et al. (2019). Visualization in Bayesian workflow. *JRSS-A*, 182(2), 389–402.
- Gelman et al. (2013). *Bayesian Data Analysis* (3rd ed.). CRC Press.
- Kermack & McKendrick (1927). A contribution to the mathematical theory of epidemics. *Proc. Royal Soc. A*, 115, 700–721.
- Vehtari et al. (2017). Practical Bayesian model evaluation using LOO-CV and WAIC. *Statistics and Computing*, 27, 1413–1432.
- Vehtari et al. (2021). Rank-normalization, folding, and localization: An improved R-hat for assessing convergence of MCMC. *Bayesian Analysis*, 16(2), 667–718.
