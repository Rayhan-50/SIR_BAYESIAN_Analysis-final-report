# Guide to Bayesian SIR Model Plots (2025 Epidemic Analysis)
This document provides a rigorous, academic breakdown of all 12 plots generated in this study, designed for direct presentation to your professor.

---

## Executive Summary for Your Professor
We fit a mechanistic **Susceptible-Infectious-Recovered (SIR)** model parameterized by Ordinary Differential Equations (ODEs) to a 52-week case dataset ($N = 11,750$) from 2025. The parameters were estimated using **Bayesian Markov Chain Monte Carlo (MCMC)** via Stan. 

The model achieves **excellent MCMC convergence** ($\hat{R} \le 1.005$) and explains **84% of the variance** in case counts ($R^2 = 0.84$). However, it displays a characteristic **systematic underprediction of the epidemic peak** (predicting a median peak of ~47 cases vs. the observed 69). This is a known structural limitation of homogeneous-mixing compartmental models, likely caused by spatial clustering or super-spreading events in the real-world population.

---

## 1. Raw Data Analysis
### [Plot 01: Weekly Case Counts (Raw Data)](file:///d:/research/research/SIR_Bayesian_Analysis_Report/plots/01_raw_data.png)

* **What it represents:** The time series of weekly reported case counts across 52 weeks of the year 2025.
* **Key Visual Features:** A textbook single-wave, unimodal (bell-shaped) outbreak. Slow growth from weeks 0 to 20, acceleration to a sharp peak of **69 cases in Week 30 (early August)**, followed by a symmetric, exponential decay back to near-zero by Week 48.
* **Interpretation for your Professor:** 
  * The unimodal, symmetric shape justifies using a standard deterministic SIR model structure.
  * There are no secondary waves or complex endemic dynamics in this year-long span.
  * The cumulative reported cases equal $833$ out of $N = 11,750$, indicating a raw observed attack rate of $\approx 7.1\%$.

---

## 2. Compartmental Dynamics
### [Plot 02: SIR Compartmental Trajectories](file:///d:/research/research/SIR_Bayesian_Analysis_Report/plots/02_sir_compartments.png)

* **What it represents:** The underlying, unobserved state variables ($S(t)$, $I(t)$, $R(t)$) of the SIR ODE system, overlaying the scaled infectious population ($\rho \cdot I(t)$) with the observed case counts.
* **Key Visual Features:**
  * **Susceptible ($S(t)$):** Drops monotonically from $11,749$ to $\approx 5,400$ by Week 52.
  * **Infectious ($I(t)$):** Rises to peak around Weeks 27–28 before returning to near-zero.
  * **Recovered ($R(t)$):** Grows in an S-shaped curve (sigmoidal) to $\approx 6,350$.
  * **Expected Reported Cases ($\rho \cdot I(t)$):** Dotted line tracking closely with the observed data points (circles) on the secondary axis.
* **Interpretation for your Professor:**
  * While the reported attack rate was only $7.1\%$ ($833$ cases), the model estimates that the **true attack rate** was **$54\%$** ($6,350$ infections).
  * This discrepancy is resolved by the reporting fraction ($\rho \approx 0.08$), meaning only $8\%$ of actual infections were clinically recorded. The remaining $92\%$ represents asymptomatic or sub-clinical spread.

---

## 3. Prior Predictive Validation
### [Plot 03: Prior Predictive Check](file:///d:/research/research/SIR_Bayesian_Analysis_Report/plots/03_prior_predictive.png)

* **What it represents:** Simulating epidemic trajectories using parameter values drawn entirely from the prior distributions *before* exposing the model to the data.
* **Key Visual Features:** A wide cloud of 200 grey curves. Peak heights range from $0$ to $>300$ cases/week, and peak timings span from Week 10 to Week 45.
* **Interpretation for your Professor:**
  * **Verdict: PASS.** The priors are weakly informative and appropriate.
  * The observed dataset (peak = 69, peaking in Week 30) lies comfortably within the center of the prior predictive cloud.
  * This proves that our priors are broad enough not to bias the posterior calculations, but tight enough to prevent the model from generating mathematically impossible dynamics (e.g., negative cases or infinite growth).

---

## 4. MCMC Performance & Convergence
### [Plot 04: Trace Plots](file:///d:/research/research/SIR_Bayesian_Analysis_Report/plots/04_trace_plots.png)

* **What it represents:** The values of parameters ($\beta$, $\gamma$, $\rho$, $\phi$) plotted against the iteration step for four independent, parallel MCMC chains.
* **Key Visual Features:** Overlapping, highly dense, stationary "hairy caterpillars" for all parameters. No chain drifts, diverges, or gets stuck in local minima.
* **Interpretation for your Professor:**
  * **Verdict: CONVERGED.** The trace plots show excellent mixing.
  * The Gelman-Rubin convergence diagnostic **$\hat{R} \le 1.005$** (ideal is $1.000$) and the Effective Sample Sizes **$ESS \ge 946$** (ideal is $>400$) quantitatively confirm that the MCMC chains have fully converged to the joint posterior distribution.

---

## 5. Posterior Estimations
### [Plot 05: Posterior Histograms](file:///d:/research/research/SIR_Bayesian_Analysis_Report/plots/05_posterior_histograms.png)

* **What it represents:** The marginal posterior probability density distributions for the model parameters after fitting.
* **Key Visual Features:** Clean, unimodal, approximately symmetric (bell-shaped) distributions.
* **Parameter Breakdown for your Professor:**
  * **$\beta$ (transmission rate):** Median $0.947$ (95% CrI: $[0.804, 1.137]$). Indicates an infectious individual contacts susceptible individuals at a rate of ~0.95/week.
  * **$\gamma$ (recovery rate):** Median $0.657$ (95% CrI: $[0.488, 0.865]$). This yields a **mean infectious period** ($1/\gamma$) of **$1.52$ weeks** (~10.6 days).
  * **$\rho$ (reporting fraction):** Median $0.079$ (95% CrI: $[0.048, 0.131]$). Quantifies the severe underreporting (~8% detected).
  * **$\phi$ (overdispersion):** Median $4.44$ (95% CrI: $[2.41, 7.65]$). Confirms overdispersion relative to a Poisson distribution (which requires $\phi \to \infty$).

### [Plot 09: $R_0$ Posterior Distribution](file:///d:/research/research/SIR_Bayesian_Analysis_Report/plots/09_R0_distribution.png)

* **What it represents:** The probability density of the Basic Reproduction Number, $R_0$, calculated as $R_0 = \beta / \gamma$.
* **Key Visual Features:** Narrow, bell-shaped distribution centered at **$1.44$** (95% CrI: $[1.31, 1.65]$).
* **Interpretation for your Professor:**
  * The entire posterior distribution lies **firmly above the epidemic threshold of $1.0$**. The probability that $R_0 \le 1.0$ is $0.000$.
  * An $R_0 \approx 1.44$ indicates a moderately transmissible pathogen, highly consistent with seasonal influenza.

---

## 6. Joint Posteriors & Identifiability
### [Plot 06: Pair Plots (Joint Densities)](file:///d:/research/research/SIR_Bayesian_Analysis_Report/plots/06_pair_plots.png)

* **What it represents:** Scatter plots and contour diagrams showing correlations between parameters in their joint posterior distribution.
* **Key Visual Features:**
  * Extreme linear collinearity between **$\beta$ and $\gamma$ ($r = 0.997$)**.
  * Strong collinearity between **$\beta$ and $\rho$ ($r = 0.934$)** and **$\gamma$ and $\rho$ ($r = 0.933$)**.
* **Interpretation for your Professor:**
  * **Critical Scientific Insight:** This correlation highlights a well-known *structural non-identifiability* in compartmental models. The clinical data (weekly counts) only constraints the ratio $R_0 = \beta / \gamma$, but cannot easily separate how *infectious* the disease is ($\beta$) from how *long* patients remain infectious ($1/\gamma$).
  * Similarly, a higher transmission rate ($\beta$) combined with lower reporting ($\rho$) can produce identical case data to a lower transmission rate with higher reporting.
  * Despite this collinearity, the ratio $R_0$ is extremely well-identified (CrI width of only $0.34$).

---

## 7. Model Fit & Predictive Capability
### [Plot 07: Posterior Predictive Check (Expected Cases)](file:///d:/research/research/SIR_Bayesian_Analysis_Report/plots/07_posterior_predictive.png)

* **What it represents:** The model's expected reported case counts (posterior median trajectory) with 50% and 90% credible intervals, compared to the actual observed cases.
* **Key Visual Features:** 
  * The overall fit is excellent ($R^2 = 0.84$, $RMSE = 8.63$).
  * The model captures the growth phase (Weeks 1–22) and the decay phase (Weeks 33–52) with high fidelity.
  * **Systematic Mismatch:** The model underpredicts the peak. Observed cases reached 69 in Week 30, but the model's median prediction peaks at **47 cases**, and the 90% CI upper bound (~57) does not cover the observed peak.
* **Interpretation for your Professor:**
  * Compartmental models with homogeneous mixing assume uniform contact rates, which naturally "smooths out" peaks. 
  * In reality, the peak was likely driven by temporary contact heterogeneity (e.g., school terms, superspreader events, or localized outbreaks), which a basic SIR model cannot capture.

### [Plot 08: Posterior Predictive Simulations](file:///d:/research/research/SIR_Bayesian_Analysis_Report/plots/08_posterior_predictive_sim.png)

* **What it represents:** 200 individual simulation runs where parameter uncertainty and observation noise (Negative Binomial distribution) are both simulated.
* **Key Visual Features:** A cloud of grey trajectories. The observed data points (black dots) are entirely enveloped by the cloud.
* **Interpretation for your Professor:**
  * **Verdict: PASS.** While the expected mean underpredicts the peak (Plot 07), the simulated data shows that an observed peak of 69 is highly plausible under the Negative Binomial observation model. 
  * It validates that the error model ($\phi \approx 4.4$) is correctly sized to account for real-world stochasticity.

---

## 8. Residual Analysis
### [Plot 10: Residual Diagnostics](file:///d:/research/research/SIR_Bayesian_Analysis_Report/plots/10_residuals.png)

* **What it represents:** Plot of residuals (Observed $-$ Fitted) over time (left) and residuals vs. fitted value (right).
* **Key Visual Features:**
  * **Left Panel:** A clear U-shaped pattern. Residuals are slightly positive in weeks 1–18, become negative in weeks 19–24, spike highly positive (+15 to +29) during the peak (weeks 25–31), and return to near-zero.
  * **Right Panel:** Residual variance increases as the fitted value increases (heteroscedasticity).
* **Interpretation for your Professor:**
  * The U-shaped residual trend confirms a systematic structural misfit around the peak (the model is too flat).
  * The heteroscedasticity is mathematically expected since we used a Negative Binomial likelihood where variance scales quadratically with the mean ($Var = \mu + \mu^2/\phi$).

---

## 9. Presentation Figures
### [fitted_model_presentation.png](file:///d:/research/research/SIR_Bayesian_Analysis_Report/plots/fitted_model_presentation.png)
* A high-resolution, publication-ready visualization with calendar dates on the X-axis, explicit annotations of the median parameters ($\beta$, $\gamma$, $\rho$), $R_0$ with its 95% CrI, $R^2$, and RMSE. It is styled in a clean serif typography suitable for presentation slides.

### [fitted_model_with_residuals.png](file:///d:/research/research/SIR_Bayesian_Analysis_Report/plots/fitted_model_with_residuals.png)
* A vertically stacked publication figure combining the fitted trajectory (top panel) and the weekly residuals (bottom panel). It is ideal for inclusion in the methods section of a manuscript or thesis report, as it displays both fit quality and systematic errors transparently.
