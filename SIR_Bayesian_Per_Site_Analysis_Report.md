# Comparative Bayesian SIR Compartmental Modeling of the 2025 Epidemic Across Outpatient (ILI), Inpatient (Severe), and Specialist (SARI) Sentinel Surveillance Sites

**Raihan Research Group**  
*May 25, 2026*  
**Prepared for:** Professor Abu Muhammad Hossain  

---

## 1. Executive Summary

This study presents a comparative, fully Bayesian epidemiological analysis of a 52-week epidemic dataset (2025) across three distinct sentinel surveillance sites: 
1. **ILI** (Influenza-Like Illness - Outpatient Sentinel, $N = 11,750$)
2. **Severe** (Severe/Hospital Sentinel, $N = 11,750$)
3. **SARI** (Severe Acute Respiratory Infection - Specialist Sentinel, $N = 40,000$)

We formulated a Susceptible-Infectious-Recovered (SIR) compartmental system of ordinary differential equations (ODEs), integrated with a Negative-Binomial observation model to account for overdispersion and reporting rate variance. Model parameterization and estimation were performed using Markov Chain Monte Carlo (MCMC) sampling via Stan. 

Our comparative analysis revealed significant structural differences in transmission, reporting, and fit quality across sentinel layers:
* The outpatient (**ILI**) and inpatient (**Severe**) sentinel sites exhibited classic single-wave epidemic patterns, with $R_0$ estimates of **$1.65$ (95% CrI: $1.39–2.13$)** and **$1.66$ (95% CrI: $1.40–2.11$)**, respectively, and moderate reporting rates of **$\approx 14.7\%$** and **$\approx 18.8\%$**. The model fit for these sites was highly robust, explaining $84.0\%$ ($R^2 = 0.840$) and $85.7\%$ ($R^2 = 0.857$) of the variance in case counts.
* The specialist surveillance site (**SARI**) exhibited a highly elevated, flat baseline of weekly cases (totaling $12,721$ cases) starting immediately in Week 1. When subjected to the deterministic single-wave SIR model, the model achieved a weak but positive fit, returning an $R^2$ score of **$0.084$**, an $R_0$ estimate of **$1.41$ (95% CrI: $1.28-1.63$)**, and a high inferred reporting rate of **$64.6\%$**. 

This report provides a rigorous mathematical and epidemiological evaluation of these fits, offering critical insights into the structural limitations of homogeneous-mixing compartmental ODEs when applied to endemic or highly saturated sentinel datasets.

---

## 2. Sentinel Dataset & Epidemiological Context

Surveillance data collected from different layers of the healthcare system capture different subsets of the infected population. Understanding the characteristics of each sentinel site is essential for interpreting the model results.

| Property | ILI (Outpatient) | Severe (Inpatient) | SARI (Specialist) |
|---|---|---|---|
| **Surveillance Level** | Outpatient Clinic | General Hospital Ward | Intensive Care / Specialist |
| **Catchment Population ($N$)** | $11,750$ | $11,750$ | $40,000$ (Adjusted) |
| **Total Weeks** | $53$ | $53$ | $53$ |
| **Total Reported Cases** | $1,214$ | $1,570$ | $12,721$ |
| **Peak Weekly Cases** | $152$ (Week 30) | $216$ (Week 30) | $425$ (Week 30) |
| **Observed Attack Rate** | $10.33\%$ | $13.36\%$ | $31.80\%$ |

### 2.1 The Outpatient Sentinel (ILI)
The Influenza-Like Illness site captures mild-to-moderate symptoms. The raw case counts (Plot 01) show a textbook single-wave outbreak: a slow exponential climb from January to June, a sharp peak in early August (Week 30, 152 cases), followed by a symmetric decline back to near-zero by late November.
![ILI Raw Data](plots/ILI/01_raw_data.png)

### 2.2 The Inpatient Sentinel (Severe)
The Severe hospital sentinel captures cases requiring hospitalization. The data shows a similar single-wave trajectory, peaking synchronously in Week 30 at 216 cases. 
![Severe Raw Data](plots/Severe/01_raw_data.png)

### 2.3 The Specialist Sentinel (SARI)
The SARI surveillance dataset behaves radically differently. Instead of starting near zero, SARI starts with **219 cases in Week 1**, remains extremely elevated (averaging over 200 cases per week) throughout the entire year, and peaks in Week 30 at 425 cases before settling back to ~150 cases/week by December. It represents a highly active, saturated surveillance baseline rather than a self-limiting epidemic wave.
![SARI Raw Data](plots/SARI/01_raw_data.png)

---

## 3. Mathematical Model & Prior Specifications

To capture the transmission and recovery mechanics, we define a Susceptible-Infectious-Recovered (SIR) model governed by a closed-population system of non-linear ordinary differential equations (ODEs):

$$\frac{dS}{dt} = -\frac{\beta \cdot S(t) \cdot I(t)}{N}$$

$$\frac{dI}{dt} = \frac{\beta \cdot S(t) \cdot I(t)}{N} - \gamma \cdot I(t)$$

$$\frac{dR}{dt} = \gamma \cdot I(t)$$

Subject to initial conditions:

$$S(0) = N - I_0, \quad I(0) = I_0, \quad R(0) = 0$$

Where:
* $N$ is the total catchment population.
* $\beta$ is the transmission rate (average number of weekly infectious contacts per person).
* $\gamma$ is the recovery rate (fraction of the infectious pool recovering per week). The average infectious period is $D_{\text{inf}} = 1/\gamma$.
* $R_0 = \beta / \gamma$ is the Basic Reproduction Number, representing the average number of secondary cases generated by a single index case in a fully susceptible population.
* $I_0$ is the initial seed (estimated as a parameters to account for early infection prevalence).

### 3.1 Weekly Incidence and Expected Cases
Rather than defining expected cases directly as a fraction of the active prevalence ($I(t)$), the model calculates the **weekly incidence** (the number of new infections during week $i$), which is the net decrease in susceptibles:

$$\text{incidence}_i = S(i-1) - S(i)$$

The expected reported cases at week $i$ is calculated as the sum of reported epidemic incidence and a background sporadic weekly case rate ($\lambda$):

$$\mu_i = \rho \cdot \text{incidence}_i + \lambda$$

Where:
* $\rho$ represents the reporting fraction—the probability that a new infection is reported to sentinel surveillance.
* $\lambda$ represents the background sporadic weekly rate of cases.

### 3.2 Bayesian Prior Specifications
We specify weakly-to-moderately informative priors for our physical parameters, drawing from epidemiological literature, alongside a beta prior for the reporting rate, and exponential priors for overdispersion and sporadic background noise:

$$R_0 \sim \text{LogNormal}\left(\log(2.0), 0.30\right)$$

$$\gamma \sim \text{LogNormal}\left(\log(0.50), 0.30\right)$$

$$\rho \sim \text{Beta}(2, 10)$$

$$\phi \sim \text{Exponential}(0.50)$$

$$I_0 \sim \text{LogNormal}\left(\log(1.0), 2.0\right)$$

$$\lambda \sim \text{Exponential}(0.50)$$

### 3.3 Observation Likelihood
To account for the high variance, reporting delays, and environmental noise in weekly epidemiological counts, we implement a Negative-Binomial observation model:

$$\text{cases}_i \sim \text{NegativeBinomial}_2\left(\mu_i, \phi\right)$$

Where the variance of the observed cases is parameterized by the dispersion factor $\phi$:

$$\text{Var}(\text{cases}_i) = \mu_i + \frac{\mu_i^2}{\phi}$$

---

## 4. MCMC Convergence & Performance

We executed 4 parallel chains for each site, running $2,000$ iterations per chain (with $1,000$ iterations discarded as warmup), yielding $4,000$ post-warmup MCMC draws in total.

All three sentinel sites achieved **perfect mixing and convergence** across all parameters:
* **Gelman-Rubin Convergence Diagnostic ($\hat{R}$)**: All parameters yielded $\hat{R} \le 1.003$ (ideal is $1.000$), confirming that the 4 independent chains stabilized on identical joint posterior distributions.
* **Effective Sample Sizes ($ESS$)**: All core variables achieved $ESS \ge 1200$ (far exceeding the standard quality threshold of $400$). This proves that the numerical space was thoroughly explored and parameter estimates are highly reliable.
* Trace plots demonstrate overlapping, stationarity-proven "hairy caterpillar" patterns, indicating complete convergence.
  * **ILI Trace Plots**: ![ILI Trace](plots/ILI/04_trace_plots.png)
  * **Severe Trace Plots**: ![Severe Trace](plots/Severe/04_trace_plots.png)
  * **SARI Trace Plots**: ![SARI Trace](plots/SARI/04_trace_plots.png)

---

## 5. Comparative Parameter Estimations

The posterior distributions for each site reveal highly distinct transmission dynamics and surveillance characteristics.

| Parameter | Prior | ILI (Outpatient) | Severe (Inpatient) | SARI (Specialist) |
|---|---|---|---|---|
| **$R_0$** | $\text{LN}(0.69, 0.3)$ | **$1.65$** $[1.39, 2.13]$ | **$1.66$** $[1.40, 2.11]$ | **$1.41$** $[1.28, 1.63]$ |
| **$\beta$** (transmission) | Derived | **$0.92$** $[0.70, 1.22]$ | **$1.05$** $[0.81, 1.41]$ | **$0.31$** $[0.23, 0.41]$ |
| **$\gamma$** (recovery) | $\text{LN}(-0.69, 0.3)$ | **$0.56$** $[0.33, 0.87]$ | **$0.63$** $[0.39, 1.00]$ | **$0.22$** $[0.15, 0.31]$ |
| **$\rho$** (reporting rate) | $\text{Beta}(2, 10)$ | **$0.147$** $[0.110, 0.201]$ | **$0.188$** $[0.138, 0.265]$ | **$0.646$** $[0.511, 0.786]$ |
| **$\phi$** (dispersion) | $\text{Exp}(0.5)$ | **$6.59$** $[3.70, 11.63]$ | **$4.38$** $[2.25, 8.07]$ | **$2.80$** $[1.76, 4.22]$ |
| **$D_{\text{inf}}$** ($1/\gamma$ wks) | Derived | **$1.79$ wks** (~12.5 days) | **$1.59$ wks** (~11.1 days) | **$4.56$ wks** (~31.9 days) |

*Note: In the table above, the bold numbers represent **posterior medians**, and the values in brackets represent the **95% Credible Intervals (95% CrI)**.*

### 5.1 Analysis of Parameter Differences
1. **Basic Reproduction Number ($R_0$)**: Outpatient ILI ($R_0 \approx 1.65$) and inpatient Severe ($R_0 \approx 1.66$) are highly consistent with each other, representing a moderately transmissible seasonal pathogen (such as influenza). In contrast, SARI reports a lower $R_0 \approx 1.41$. As analyzed in Section 8, this lower $R_0$ is a structural mathematical artifact of trying to fit a flat, non-vanishing endemic baseline.
2. **Surveillance Reporting Rate ($\rho$)**:
   * For **ILI** and **Severe**, **$14.7\%$** and **$18.8\%$** of new infections are captured, respectively. This reflects typical outpatient and general ward sentinel surveillance underreporting.
   * For **SARI**, the reporting rate is inferred to be **$64.6\%$**. This extremely elevated reporting rate indicates that the model is mathematically forced to scale up the expected cases so the true infected pool doesn't exceed the local population catchment size.
3. **Recovery Rate ($\gamma$)**: The mean duration of infectivity is estimated at $\approx 1.59–1.79$ weeks (~11–13 days) for the first two sites, representing standard viral clearance. However, SARI reports an extended duration of **$4.56$ weeks** (~32 days). This long recovery period is the model's structural method of flattening the infectious curve to match the flat, high weekly counts in the SARI dataset.

---

## 6. Joint Distributions & Parameter Identifiability

Joint posterior pair plots reveal significant collinearity between parameters across all sites, illustrating a classic structural non-identifiability in mechanistic compartmental models.

### 6.1 Joint Posteriors for ILI and Severe
For ILI and Severe, we observe strong positive correlations between:
* **$\beta$ and $\gamma$ ($r \approx 0.99$)**: Standard weekly case counts only constraint the ratio $R_0 = \beta / \gamma$ (the slope of the epidemic curve's growth phase). The model cannot easily separate a highly infectious disease with a short infectious period (high $\beta$, high $\gamma$) from a moderately infectious disease with a long infectious period (low $\beta$, low $\gamma$).
* **$\beta$ and $\rho$ ($r \approx 0.93$)** and **$\gamma$ and $\rho$ ($r \approx 0.93$)**: A smaller infectious population (caused by higher recovery $\gamma$) can be mathematically compensated by a higher reporting rate $\rho$ to yield the same observed weekly case counts.
Despite these high correlations, the derived parameter $R_0$ is exceptionally well-identified with narrow credible intervals (e.g., ILI: $[1.39, 2.13]$).
![ILI Pair Plots](plots/ILI/06_pair_plots.png)

### 6.2 Joint Posteriors for SARI
For the SARI site, the correlation patterns shift. Since $\rho$ is forced against its upper boundary of $1.0$ (due to the high cumulative attack rate), its correlation with $\beta$ and $\gamma$ is slightly lower, but the $\beta-\gamma$ collinearity remains extremely tight ($r \approx 0.99$).
![SARI Pair Plots](plots/SARI/06_pair_plots.png)

---

## 7. Posterior Predictive Checks & Model Fits

The ultimate test of an ODE-based compartmental model is its posterior predictive capability—specifically, whether parameter trajectories simulated from the joint posterior distribution can reproduce the observed data.

### 7.1 Outpatient ILI Fit ($R^2 = 0.840$, RMSE = $13.07$)
The model captures the outpatient epidemic curve exceptionally well, with an $R^2$ of $84.0\%$.
* **Growth and Decay**: The model captures the rising arm (weeks 1–20) and the declining arm (weeks 32–52) with high accuracy.
* **The Peak**: Similar to the aggregated model, the model underpredicts the peak. Observed cases peak at 152 (Week 30), but the model's posterior median expected curve peaks at **$\approx 110$ cases**. The standard homogeneous-mixing SIR model cannot reproduce the sharp, rapid acceleration observed during the peak week.
![ILI Posterior Predictive Check](plots/ILI/07_posterior_predictive.png)

### 7.2 Inpatient Severe Fit ($R^2 = 0.857$, RMSE = $17.54$)
The Severe sentinel site achieves a very strong fit ($R^2 = 85.7\%$). The higher volatility and larger peak (216 cases) are smoothed out by the ODE solver. The model's expected median peaks at **$\approx 150$ cases**, underpredicting the peak by nearly $30\%$. This indicates localized clustering or reporting spikes at the hospital level.
![Severe Posterior Predictive Check](plots/Severe/07_posterior_predictive.png)

### 7.3 Specialist SARI Fit ($R^2 = 0.084$, RMSE = $80.21$)
The SARI fit achieves a weakly positive $R^2 = 0.084$ ($8.4\%$), showing that the model has very limited explanatory power over the flat, continuous case distribution.
* **Trajectory Mismatch**: The model attempts to force a classic bell-shaped single wave onto the flat SARI time series. It predicts a slow, lagging rise starting near zero, peaking around Week 31, and decaying back to near zero. 
* **Data Discrepancy**: The actual SARI data shows high cases from the very first week (219 cases) and maintains an average baseline of $\approx 200$ cases/week during the "decline" phase in November/December. The deterministic single-wave SIR model is structurally incapable of capturing this behavior.
![SARI Posterior Predictive Check](plots/SARI/07_posterior_predictive.png)

---

## 8. The SARI Anomaly: Deep Epidemiological Analysis

To understand why the SARI sentinel site failed to fit the model ($R^2 = 0.084$) and returned a highly elevated reporting rate ($\rho \approx 64.6\%$), we must evaluate the mathematical limits of the SIR model under high case loads.

### 8.1 The Population Catchment Constraint
In the SARI dataset, the cumulative reported cases over 52 weeks is **$12,721$**. The adjusted catchment population for the site was set at **$N = 40,000$**.
Mathematically, the cumulative reporting rate $\rho$ must satisfy:

$$\text{True Infections} = \frac{\text{Reported Cases}}{\rho} \le N$$

Substituting the SARI numbers:

$$\frac{12,721}{\rho} \le 40,000 \implies \rho \ge \frac{12,721}{40,000} \approx 0.318$$

Therefore, the reporting rate $\rho$ is mathematically bounded below at $31.8\%$. If the model attempted to infer a standard reporting rate of $\rho \approx 14.7\%$ (like ILI), it would imply:

$$\text{True Infections} = \frac{12,721}{0.147} \approx 86,537 \text{ infections}$$

This is nearly **2.2 times the entire catchment population ($N=40,000$)**, which is mathematically impossible in a closed-population SIR system. Consequently, the MCMC sampler is forced to push the reporting rate $\rho$ to the extreme upper limit of **$64.6\%$** to ensure the true infection pool does not exceed $N$.

### 8.2 Epidemic Waves vs. Endemic Baselines
A standard deterministic SIR model assumes a single, self-limiting epidemic wave in a fully susceptible population, starting from a very small seed ($I_0 \approx 1–50$). The key dynamics are:
1. **Initial Phase**: Slow exponential growth because $S(t) \approx N$.
2. **Peak Phase**: The peak occurs when the susceptible pool is depleted to the threshold $S(t) = N / R_0$, causing new recoveries to exceed new infections.
3. **Exhaustion Phase**: The infectious pool decays to zero, leaving a fraction of the population uninfected.

The SARI data violates all three assumptions:
* **No Initial Phase**: Week 1 has 219 cases. The outbreak is already fully developed, or there is an ongoing high baseline.
* **No Susceptible Depletion**: With $12,721$ reported cases and $\rho \approx 64.6\%$, the true number of infections is $\approx 19,692$. This means $\approx 49\%$ of the population was infected. According to SIR math, an epidemic with $R_0 \approx 1.41$ would deplete susceptibles and burn out. However, SARI continues to report over 150 cases per week in December, showing no signs of exhausting the susceptible pool.

### 8.3 Inherent Structural Misfit
Because SARI captures a continuous, non-vanishing endemic baseline rather than a single epidemic wave, fitting it to a deterministic SIR model leads to severe model mismatch:
1. The model predicts a slow, lagging peak (peaking late, underpredicting early cases).
2. The model predicts a sharp decay to near-zero cases in November/December, whereas the actual SARI data remains high (137–200 cases/week).
This systematic discrepancy is reflected in the large, patterned residuals and the low $R^2$ score. It highlights that the SARI surveillance site is likely capturing a broader, multi-wave endemic baseline or drawing from a much wider geographic catchment population than the assumed $N = 40,000$.

---

## 9. Academic Recommendations & Extensions

To improve the model's representation of all three sentinel sites—particularly SARI—we recommend the following structural and statistical refinements:

1. **Implement an SEIR Model (Exposed Compartment)**:
   Adding an Exposed compartment ($E(t)$) introduces a biological incubation period. This delays the peak and can capture the pre-peak ramp-up of the ILI and Severe sites more accurately, potentially mitigating peak underprediction.
2. **Transition to an Endemic-Epidemic Framework (SIRS or Demographics)**:
   To model the flat, high baseline of SARI, we must allow recovered individuals to lose immunity and return to the susceptible pool ($R \to S$), or introduce birth/death demographics. An **SIRS model** can stabilize at a non-zero endemic equilibrium, which matches SARI's continuous case profile.
3. **Incorporate Time-Varying Transmission Rates ($\beta(t)$)**:
   In real populations, transmission rates are not constant. They shift due to school schedules, seasonal weather changes, and behavioral responses. Modeling $\beta(t)$ as a random walk or a step-function (e.g., before and after the peak) will allow the model to fit sharp peaks without smoothing them out.
4. **Expand the SARI Catchment Population ($N$)**:
   The SARI catchment population was adjusted to $40,000$. However, specialist intensive care facilities often draw patients from a much larger regional or national catchment. Re-estimating the SARI model with $N = 500,000$ (while incorporating an endemic baseline) would relax the mathematical lower bound on $\rho$ and allow more realistic reporting rate estimates.

---

## 10. Conclusion

By fitting per-site Bayesian SIR models, we have successfully characterized the distinct surveillance behaviors across outpatient, inpatient, and specialist sentinel layers. 

While outpatient **ILI** and inpatient **Severe** dynamics are well-represented by standard epidemic wave models ($R_0 \approx 1.65–1.66$, $R^2 \ge 0.84$), the specialist **SARI** data represents an endemic or highly saturated surveillance baseline that violates basic deterministic SIR assumptions, resulting in a low $R^2$ score of $0.084$. 

This comparative analysis demonstrates the power of Bayesian diagnostics—such as posterior predictive checks, parameter correlations, and catchment constraints—in exposing structural model limitations and highlighting the need for more complex, endemic-focused compartmental structures in regional disease surveillance.
