# Line-by-Line Code Explanation for Professor Meeting

This guide provides a detailed line-by-line breakdown of the two core scripts driving our Bayesian SIR project: `code/sir_model.stan` (the mathematical framework) and `code/analysis.R` (the execution, fitting, and plotting script). 

Use this document to confidently answer your professor when they ask: *"What exactly is this code doing here?"*

---

## Part 1: The Bayesian Model (`code/sir_model.stan`)

Stan is a probabilistic programming language used for Bayesian inference. This script is compiled into C++ and runs the Hamiltonian Monte Carlo (HMC) algorithm to find our parameter estimates.

### 1. Functions Block (The Mathematics)
```stan
functions {
  vector sir_ode(real t, vector y, array[] real theta, array[] real x_r, array[] int x_i) {
```
* **Lines 1-7:** Defines the ordinary differential equation (ODE) function. `t` is time, `y` holds the compartments (S, I, R), and `theta` holds our parameters ($\beta$ and $\gamma$).

```stan
    real beta  = theta[1];
    real gamma = theta[2];
    real N     = x_r[1];
    real S = y[1];
    real I = y[2];
```
* **Lines 8-14:** Extracts the specific values from the input arrays so we can use readable names like `beta`, `gamma`, and `S` (Susceptible), `I` (Infectious).

```stan
    vector[3] dydt;
    dydt[1] = -beta * S * I / N;
    dydt[2] =  beta * S * I / N - gamma * I;
    dydt[3] =  gamma * I;
    return dydt;
  }
}
```
* **Lines 15-21:** The core SIR mathematical equations. 
  * `dydt[1]` represents $dS/dt$: Susceptibles decrease as they get infected.
  * `dydt[2]` represents $dI/dt$: Infectious people increase by new infections and decrease as they recover.
  * `dydt[3]` represents $dR/dt$: Recovered people increase based on the recovery rate ($\gamma$).

### 2. Data Block (What we feed the model)
```stan
data {
  int<lower=1> N_weeks;
  array[N_weeks] int<lower=0> cases;
  real<lower=0> pop;
}
```
* **Lines 23-27:** Declares the data we provide from R. We give it the total number of weeks (`N_weeks`), the weekly observed reported `cases`, and the total population (`pop`).

### 3. Parameters Block (What we want to learn)
```stan
parameters {
  real<lower=0>          R0;
  real<lower=0>          gamma;
  real<lower=0, upper=1> rho;
  real<lower=0>          phi;
  real<lower=0.001>      I0;      
  real<lower=0>          lambda;  
}
```
* **Lines 37-44:** These are the unknowns the MCMC algorithm is trying to estimate. 
  * `R0`: Basic reproduction number.
  * `gamma`: Recovery rate.
  * `rho`: The reporting fraction (how many true cases actually get diagnosed, constrained between 0 and 1).
  * `phi`: Overdispersion (noise in the data).
  * `I0`: Initial infected seed at week 0.
  * `lambda`: Background (sporadic) cases not caused by the main epidemic wave.

### 4. Transformed Parameters (Connecting math to data)
```stan
transformed parameters {
  real beta = R0 * gamma;
```
* **Line 47:** We calculate $\beta$ (transmission rate) directly from $R_0$ and $\gamma$. We estimate $R_0$ rather than $\beta$ directly because priors on $R_0$ are easier to justify biologically.

```stan
  array[N_weeks] vector[3] y_sol;
  y_sol = ode_rk45(sir_ode, y0, t0, ts, theta, x_r, x_i);
```
* **Lines 55-56:** This tells Stan to solve the ODE system over the specified weeks using the Runge-Kutta 4th order (`rk45`) numerical solver.

```stan
  array[N_weeks] real mu;
  for (i in 1:N_weeks) {
    mu[i] = fmax(1e-6, rho * incidence[i] + lambda);
  }
}
```
* **Lines 65-69:** Calculates the *expected* reported cases (`mu`) for each week. It equals the true new infections (`incidence`) multiplied by the reporting fraction (`rho`), plus background noise (`lambda`).

### 5. Model Block (Priors and Likelihood)
```stan
model {
  // Priors
  R0     ~ lognormal(log(2.0),  0.3);   
  gamma  ~ lognormal(log(0.5),  0.3);   
  rho    ~ beta(2, 10);                  
  phi    ~ exponential(0.5);
  I0     ~ lognormal(log(1),    2.0);   
  lambda ~ exponential(0.5);             
```
* **Lines 72-79:** **PRIORS.** This is the "Bayesian" part. Before looking at the data, we tell the model what values are biologically realistic (e.g., $R_0$ is probably around 2.0 based on flu literature).

```stan
  for (i in 1:N_weeks) {
    cases[i] ~ neg_binomial_2(mu[i], phi);
  }
}
```
* **Lines 81-84:** **LIKELIHOOD.** This connects our expected cases (`mu`) to the real observed `cases` using a Negative Binomial probability distribution. The model asks: *"How likely is the real data if our expected 'mu' is correct?"*

---

## Part 2: The R Analysis Script (`code/analysis.R`)

This script cleans data, runs the Stan model, and generates all the plots.

### 1. Setup and Data Loading
```R
suppressPackageStartupMessages({
  library(rstan)
  library(deSolve)
  # ... other libraries ...
})
```
* **Lines 6-16:** Loads required R packages. `rstan` interfaces with Stan, and `deSolve` is used to solve ODEs in R for our predictive checks.

```R
dat <- read.csv(file.path(proj_dir, "data/clean_epidemic_dataset_2025.csv"), stringsAsFactors = FALSE)
N      <- 11750L
```
* **Lines 46-50:** Loads the real 2025 case data from our CSV file. It sets the total population `N` to 11,750 (the size of the study site).

### 2. Prior Predictive Check
```R
prior_curves <- lapply(seq_len(n_prior), function(i) {
  R0_s     <- rlnorm(1, log(2.0), 0.3)
  # ... sample parameters ...
  sol_s    <- solve_sir(beta_s, gamma_s, I0 = I0_s)
```
* **Lines 89-104:** Generates **Plot 03**. It randomly draws parameters from our prior beliefs (ignoring real data) and simulates what epidemics *could* look like. This proves to the professor our priors are reasonable and not forcing impossible outcomes.

### 3. Running the MCMC Algorithm
```R
fit <- stan(
  file    = stan_file,
  data    = list(N_weeks = n_wks, cases = cases, pop = as.double(N)),
  chains  = 4,
  iter    = 2000,
  warmup  = 1000,
  seed    = 12345
)
```
* **Lines 123-132:** The core execution step. It compiles the `sir_model.stan` file and runs 4 independent MCMC chains. Each chain takes 2000 steps (discarding the first 1000 as "warmup/burn-in" to find the optimal zone). 

### 4. Convergence Diagnostics
```R
rhat_vals <- diag_df[,"Rhat"]
neff_vals <- diag_df[,"n_eff"]
```
* **Lines 140-156:** Checks if the MCMC worked correctly. $\hat{R}$ (R-hat) must be < 1.05 (proving all 4 chains agree on the results), and $n_{eff}$ (Effective Sample Size) must be > 400 (proving we have enough independent samples).

### 5. Plotting Results
```R
sol <- solve_sir(beta_med, gamma_med, I0 = max(0.001, I0_med))
```
* **Lines 178:** Takes the final estimated median values for $\beta$ and $\gamma$ and runs the SIR model one last time to generate **Plot 02** (the compartmental curves for S, I, and R).

```R
p04_trace <- mcmc_trace(posterior_arr, pars = c("R0","gamma","rho","phi","I0","lambda"))
```
* **Lines 221-225:** Generates **Plot 04** (Trace plots). This visualizes the MCMC algorithm exploring the parameter space. It looks like a "fuzzy caterpillar", which is the visual proof that the chains mixed well.

```R
p07 <- ggplot(pp_df, aes(x = week)) +
  geom_ribbon(aes(ymin = lo90, ymax = hi90), fill = "#CCCCCC", alpha = 0.7) +
  geom_line(aes(y = med), linewidth = 0.9, color = "black") +
  geom_point(aes(y = cases), shape = 19, size = 2.5, color = "black")
```
* **Lines 287-296:** Generates **Plot 07** (Posterior Predictive Check). This overlays the model's final prediction (line and shaded confidence intervals) on top of the real observed data points (black dots) to show the final model fit.

---

### What about the Images (Plots)?
To explain what each generated plot actually means mathematically and biologically to your professor, please refer strictly to the **`Plot_Explanations_Guide.md`** file, which breaks down every single image in the `plots/` folder perfectly.
