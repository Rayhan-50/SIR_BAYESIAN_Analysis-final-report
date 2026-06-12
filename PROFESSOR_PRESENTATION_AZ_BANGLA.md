# সম্পূর্ণ উপস্থাপনা গাইড — A থেকে Z
# Complete A-to-Z Professor Presentation Guide
### Bayesian SIR Epidemic Model | Raihan Research Group | 2025
**ভাষা: বাংলা + English (Bilingual)**

---

> **কিভাবে ব্যবহার করবেন:**  
> প্রতিটি সেকশনে **বাংলায়** বলার স্ক্রিপ্ট আছে, এবং সাথে **ইংরেজিতে** একই কথা আছে।  
> প্রফেসরকে দেখানোর সময় — প্রতিটি ধাপে সঠিক ফাইল বা প্লট খুলুন।

---

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# ধাপ A: প্রজেক্টের পরিচয়
# STEP A: Project Introduction
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

**📂 ফাইল খুলুন:** `README.md`

### 🇧🇩 বাংলায় বলুন:
> "স্যার, এই প্রজেক্টে আমরা ২০২৫ সালের একটি মহামারি (epidemic) ডেটাসেটকে গাণিতিকভাবে মডেল করেছি।  
> আমাদের মূল প্রশ্ন ছিল — এই রোগটি কত দ্রুত ছড়িয়েছে, কতজন আসলে আক্রান্ত হয়েছে, এবং রোগটি কতটা বিপজ্জনক।  
> আমরা তিনটি আলাদা সার্ভেইলেন্স সাইট থেকে ডেটা নিয়েছি: ILI (বাইরের রোগী), Severe (হাসপাতালে ভর্তি), এবং SARI (বিশেষজ্ঞ হাসপাতাল)।"

### 🇬🇧 English:
> "Professor, in this project we mathematically modeled a 2025 epidemic dataset. Our core questions were: how fast did this disease spread, how many people were truly infected, and what are the key transmission parameters. We collected data from three distinct sentinel surveillance sites: ILI (outpatient), Severe (hospital inpatient), and SARI (specialist)."

---

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# ধাপ B: SIR মডেলের তত্ত্ব
# STEP B: SIR Model Theory
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

**📂 কোথাও কিছু খুলতে হবে না — মুখে বলুন ও বোর্ডে লিখুন**

### 🇧🇩 বাংলায় বলুন:
> "স্যার, SIR মডেলে মোট জনগোষ্ঠীকে তিনটি ভাগে ভাগ করা হয়:  
> **S (Susceptible)** = যারা এখনো আক্রান্ত হয়নি কিন্তু হতে পারে  
> **I (Infectious)** = যারা এখন আক্রান্ত এবং অন্যকে ছড়াচ্ছে  
> **R (Recovered)** = যারা সেরে গেছে এবং রোগপ্রতিরোধ ক্ষমতা অর্জন করেছে  
>
> এই তিনটি গ্রুপ সময়ের সাথে পরিবর্তন হয়, এবং সেটা তিনটি ডিফারেনশিয়াল সমীকরণ দিয়ে লেখা হয়:"

```
dS/dt = -β × S(t) × I(t) / N     ← প্রতি সপ্তাহে কতজন নতুন আক্রান্ত হয়
dI/dt =  β × S(t) × I(t) / N  - γ × I(t)   ← নতুন আক্রান্ত বিয়োগ যারা সেরে গেছে
dR/dt =  γ × I(t)              ← প্রতি সপ্তাহে কতজন সেরে উঠছে
```

> **β (beta)** = ট্রান্সমিশন রেট — কত দ্রুত রোগ ছড়ায়  
> **γ (gamma)** = রিকভারি রেট — কত দ্রুত মানুষ সেরে ওঠে  
> **R₀ = β / γ** = মূল প্রজনন সংখ্যা — একজন আক্রান্ত ব্যক্তি গড়ে কতজনকে আক্রান্ত করে"

### 🇬🇧 English:
> "The SIR model divides the total population N into three groups. S: people who can catch the disease. I: people currently spreading it. R: people who recovered and are immune. The three differential equations describe how these groups change over time. β is the transmission rate, γ is the recovery rate, and R₀ = β/γ is the basic reproduction number."

---

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# ধাপ C: বেয়েসিয়ান পদ্ধতি কেন?
# STEP C: Why Bayesian Approach?
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

### 🇧🇩 বাংলায় বলুন:
> "স্যার, সাধারণত মানুষ β এবং γ এর মান নিজে থেকে বেছে নেয় বা গেস করে।  
> কিন্তু আমরা **Bayesian পদ্ধতি** ব্যবহার করেছি, যেখানে:  
>
> **১. Prior (আগের বিশ্বাস):** ফ্লু গবেষণার সাহিত্য থেকে আমরা জানি R₀ সাধারণত ১.৫ থেকে ২.৫ এর মধ্যে হয়। আমরা এটা মডেলকে বলে দিই।  
>
> **২. Likelihood (ডেটার সাথে মিলানো):** মডেলটি হাজার হাজার সম্ভাব্য parameter combination চেষ্টা করে এবং দেখে কোনগুলো আসল ডেটার সাথে ভালোভাবে মেলে।  
>
> **৩. Posterior (আপডেট করা বিশ্বাস):** ডেটা দেখার পর আমরা জানতে পারি β, γ, R₀ এর সঠিক মান এবং সেগুলোর সম্ভাবনার বিতরণ।  
>
> এটি শুধু একটি সংখ্যা না — এটি **পুরো সম্ভাবনার বিতরণ** দেয়, যা থেকে আমরা ৯৫% Credible Interval বলতে পারি।"

### 🇬🇧 English:
> "Instead of manually guessing parameters, we used the Bayesian framework. Step 1: Set a prior belief (e.g., R₀ is probably around 2.0 based on flu literature). Step 2: The MCMC algorithm tests thousands of parameter combinations against the real data. Step 3: The posterior distribution gives us the updated probability distribution for each parameter — not just one number, but the full uncertainty range."

---

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# ধাপ D: সফটওয়্যার সেটআপ
# STEP D: Software Setup
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

**📂 ফাইল খুলুন:** `code/install_packages.R`

### 🇧🇩 বাংলায় বলুন:
> "স্যার, আমরা দুটি প্রোগ্রামিং ল্যাঙ্গুয়েজ ব্যবহার করেছি:  
> ১. **R** — ডেটা লোড, প্লট তৈরি, এবং Stan কে নিয়ন্ত্রণ করার জন্য  
> ২. **Stan** — MCMC অ্যালগরিদম চালানো এবং গাণিতিক মডেল সংজ্ঞায়িত করার জন্য  
>
> প্রথমে এই ফাইলটি রান করলে সব প্যাকেজ ইনস্টল হয়ে যায়।"

### কোড ব্যাখ্যা (line by line):
```r
# Line 9-11: কোন কোন প্যাকেজ লাগবে তার তালিকা
pkgs <- c("deSolve", "ggplot2", "bayesplot", "gridExtra",
          "jsonlite", "dplyr", "tidyr", "GGally", ...)

# Line 14-21: প্রতিটি প্যাকেজ চেক করে, না থাকলে ইনস্টল করে
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, quiet = TRUE)    # ← ইনস্টল করো
  }
}

# Line 24-31: Stan (MCMC এর জন্য) আলাদাভাবে ইনস্টল হয়
install.packages("rstan", repos = c("https://stan-dev.r-universe.dev", ...))
```

> **প্যাকেজগুলো কী করে:**  
> - `rstan` → Stan মডেল চালায়  
> - `deSolve` → ODE সমীকরণ সমাধান করে  
> - `ggplot2` → গ্রাফ তৈরি করে  
> - `bayesplot` → MCMC trace plot তৈরি করে  
> - `dplyr`, `tidyr` → ডেটা পরিষ্কার করে

---

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# ধাপ E: Stan মডেল (গণিতের কোড)
# STEP E: The Stan Model (Mathematical Code)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

**📂 ফাইল খুলুন:** `code/sir_model.stan`

### 🇧🇩 বাংলায় বলুন:
> "স্যার, এই ফাইলটি হলো আমাদের গাণিতিক মডেলের সংজ্ঞা। Stan একটি বিশেষ ভাষা যা স্বয়ংক্রিয়ভাবে C++ কোডে রূপান্তরিত হয় এবং MCMC চালায়।  
> ফাইলটি পাঁচটি অংশে বিভক্ত:"

---

### ব্লক ১: functions (SIR সমীকরণ)
```stan
functions {
  vector sir_ode(real t, vector y, array[] real theta, ...) {
    real beta  = theta[1];   // ← β বের করা হচ্ছে
    real gamma = theta[2];   // ← γ বের করা হচ্ছে
    real N     = x_r[1];     // ← জনসংখ্যা

    real S = y[1];           // ← Susceptible
    real I = y[2];           // ← Infectious

    vector[3] dydt;
    dydt[1] = -beta * S * I / N;            // dS/dt
    dydt[2] =  beta * S * I / N - gamma*I;  // dI/dt
    dydt[3] =  gamma * I;                   // dR/dt
    return dydt;
  }
}
```
> **বাংলা:** এই অংশে আমরা SIR এর তিনটি ডিফারেনশিয়াল সমীকরণ লিখেছি।  
> dydt[1] মানে S কমছে (মানুষ আক্রান্ত হচ্ছে)  
> dydt[2] মানে I বাড়ছে (নতুন আক্রান্ত) এবং কমছে (যারা সেরে গেছে)  
> dydt[3] মানে R বাড়ছে (সুস্থ হওয়া মানুষ)

---

### ব্লক ২: data (আমরা কী দিচ্ছি)
```stan
data {
  int<lower=1> N_weeks;              // মোট কত সপ্তাহের ডেটা
  array[N_weeks] int<lower=0> cases; // প্রতি সপ্তাহের কেস সংখ্যা
  real<lower=0> pop;                 // মোট জনসংখ্যা
}
```
> **বাংলা:** আমরা R থেকে Stan কে তিনটি জিনিস দিই: সপ্তাহের সংখ্যা, কেস ডেটা, এবং জনসংখ্যা।

---

### ব্লক ৩: parameters (যা শিখতে হবে)
```stan
parameters {
  real<lower=0>          R0;       // মূল প্রজনন সংখ্যা (>0)
  real<lower=0>          gamma;    // রিকভারি রেট (>0)
  real<lower=0, upper=1> rho;      // রিপোর্টিং রেট (০ থেকে ১)
  real<lower=0>          phi;      // ওভারডিসপার্শন
  real<lower=0.001>      I0;       // শুরুতে কতজন আক্রান্ত ছিল
  real<lower=0>          lambda;   // ব্যাকগ্রাউন্ড কেস রেট
}
```
> **বাংলা:** Stan এই ছয়টি parameter এর সেরা মান খুঁজে বের করবে। আমরা ম্যানুয়ালি কিছু ঠিক করিনি।

---

### ব্লক ৪: transformed parameters (গণনা)
```stan
transformed parameters {
  real beta = R0 * gamma;      // ← β = R₀ × γ সূত্র থেকে বের করা হচ্ছে

  // ODE সমাধান করা হচ্ছে (Runge-Kutta পদ্ধতিতে)
  y_sol = ode_rk45(sir_ode, y0, t0, ts, theta, x_r, x_i);

  // প্রতি সপ্তাহে নতুন কেস = S কমে যাওয়া
  incidence[i] = y_sol[i-1][1] - y_sol[i][1];

  // প্রত্যাশিত রিপোর্টেড কেস = ρ × incidence + background
  mu[i] = rho * incidence[i] + lambda;
}
```
> **বাংলা:** এখানে আমরা β হিসাব করি (β = R₀ × γ)।  
> তারপর ODE সমীকরণ সমাধান করে প্রতি সপ্তাহে S, I, R এর মান বের করি।  
> μ হলো আমাদের প্রত্যাশিত রিপোর্টেড কেস।

---

### ব্লক ৫: model (Prior + Likelihood)
```stan
model {
  // PRIOR: আগের বিশ্বাস (ফ্লু সাহিত্য থেকে)
  R0     ~ lognormal(log(2.0),  0.3);  // R₀ সম্ভবত ~২ এর কাছাকাছি
  gamma  ~ lognormal(log(0.5),  0.3);  // রিকভারি সম্ভবত ~২ সপ্তাহ
  rho    ~ beta(2, 10);                // রিপোর্টিং ~১৭% মানে
  phi    ~ exponential(0.5);           // নয়েজের উপর দুর্বল prior
  I0     ~ lognormal(log(1), 2.0);     // শুরুতে ১ জনের কাছাকাছি
  lambda ~ exponential(0.5);           // ব্যাকগ্রাউন্ড কম থাকবে

  // LIKELIHOOD: আসল ডেটার সাথে মিলানো
  for (i in 1:N_weeks) {
    cases[i] ~ neg_binomial_2(mu[i], phi);  // Negative Binomial বিতরণ
  }
}
```
> **বাংলা:**  
> Prior অংশে আমরা বলছি — "ফ্লুর জন্য R₀ সাধারণত ২ এর কাছে হয়, সেটা আমাদের শুরুর অনুমান।"  
> Likelihood অংশে Stan জিজ্ঞেস করছে — "যদি μ সঠিক হয়, তাহলে আসল কেস ডেটা কতটা সম্ভব?"  
> Negative Binomial ব্যবহার করা হয়েছে কারণ সাপ্তাহিক কেস ডেটায় স্বাভাবিকের চেয়ে বেশি variability থাকে।

---

### ব্লক ৬: generated quantities (চেক ও সিমুলেশন)
```stan
generated quantities {
  array[N_weeks] real cases_rep;   // মডেল থেকে সিমুলেটেড কেস
  array[N_weeks] real log_lik;     // প্রতি সপ্তাহের log-likelihood

  for (i in 1:N_weeks) {
    cases_rep[i] = neg_binomial_2_rng(mu[i], phi);     // নতুন ডেটা সিমুলেট করো
    log_lik[i]   = neg_binomial_2_lpmf(cases[i] | mu[i], phi); // fit কতটা ভালো
  }
}
```
> **বাংলা:** এই অংশটি Plot 08 (সিমুলেটেড কেস) তৈরি করে এবং মডেলের ফিটনেস পরিমাপ করে।

---

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# ধাপ F: R কোড — প্রধান বিশ্লেষণ
# STEP F: R Code — Main Analysis
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

**📂 ফাইল খুলুন:** `code/analysis_per_site.R`

### 🇧🇩 বাংলায় বলুন:
> "স্যার, এই R ফাইলটি পুরো বিশ্লেষণ পরিচালনা করে। এটি তিনটি সাইটের জন্য একই কোড চালায়।"

---

### অংশ ১: প্যাকেজ লোড (Line 10-20)
```r
suppressPackageStartupMessages({
  library(rstan)      # Stan চালানোর জন্য
  library(deSolve)    # ODE সমাধান করার জন্য (prior predictive check এ)
  library(ggplot2)    # গ্রাফ আঁকার জন্য
  library(bayesplot)  # MCMC trace plot এর জন্য
  library(gridExtra)  # একসাথে অনেক plot দেখানোর জন্য
  library(dplyr)      # ডেটা পরিষ্কার ও ফিল্টার করার জন্য
})
```
> **বাংলা:** suppressPackageStartupMessages মানে হলো প্যাকেজ লোডের সময় অপ্রয়োজনীয় বার্তা লুকানো।

---

### অংশ ২: সাইট সংজ্ঞায়িত করা (Line 40-59)
```r
sites <- list(
  list(label="ILI",    csv="data/site_ILI.csv",    pop=11750L, desc="Outpatient"),
  list(label="Severe", csv="data/site_Severe.csv", pop=11750L, desc="Hospital"),
  list(label="SARI",   csv="data/site_SARI.csv",   pop=40000L, desc="Specialist")
)
```
> **বাংলা:** এখানে তিনটি সাইটের তথ্য একটি তালিকায় রাখা হয়েছে।  
> ILI এবং Severe এর জনসংখ্যা 11,750 (একই ক্যাচমেন্ট এলাকা)।  
> SARI এর জনসংখ্যা 40,000 কারণ এটি বড় বিশেষজ্ঞ হাসপাতাল থেকে ডেটা নিয়েছে।

---

### অংশ ৩: Stan মডেল কম্পাইল (Line 98-101)
```r
stan_model_obj <- stan_model(file = stan_file)
```
> **বাংলা:** Stan ফাইলটি একবার C++ এ কম্পাইল হয়, তারপর তিনটি সাইটের জন্য পুনরায় ব্যবহার হয়। এতে সময় বাঁচে।

---

### অংশ ৪: মূল লুপ (Line 108-433)
```r
for (site in sites) {   # ← প্রতিটি সাইটের জন্য
  # ডেটা লোড করো
  dat   <- read.csv(site$csv)
  cases <- dat$cases
  
  # Stan মডেল ফিট করো
  fit <- sampling(
    stan_model_obj,
    data    = list(N_weeks = n_wks, cases = cases, pop = as.double(N)),
    chains  = 4,       # ← ৪টি স্বাধীন MCMC চেইন চালাও
    iter    = 2000,    # ← প্রতি চেইনে ২০০০ ধাপ
    warmup  = 1000,    # ← প্রথম ১০০০ ধাপ বাতিল (burn-in)
    seed    = 12345    # ← পুনরুৎপাদনযোগ্যতার জন্য (reproducibility)
  )
}
```
> **বাংলা:**  
> `chains = 4` মানে আমরা একই মডেল ৪টি ভিন্ন শুরু থেকে চালাই।  
> যদি সব চেইন একই জায়গায় এসে পৌঁছায়, তাহলে মডেল কনভার্জ করেছে।  
> `warmup = 1000` মানে প্রথম ১০০০ ধাপ হলো "অন্বেষণ" সময়, এগুলো ফলাফলে গণনা হয় না।  
> বাকি ১০০০ × ৪ = **৪০০০ posterior sample** ব্যবহার করা হয়।

---

### অংশ ৫: কনভার্জেন্স চেক (Line 225-249)
```r
rhat_vals <- diag_df[,"Rhat"]    # R-hat মান বের করা হচ্ছে
neff_vals <- diag_df[,"n_eff"]   # Effective Sample Size বের করা হচ্ছে

if (all(rhat_vals < 1.05)) {
  cat("PASS: All R-hat < 1.05")  # ← চেইনগুলো একমত
}
if (all(neff_vals > 400)) {
  cat("PASS: All n_eff > 400")   # ← যথেষ্ট sample আছে
}
```
> **বাংলা:**  
> **R-hat (R̂):** যদি < ১.০৫ হয়, তাহলে সব ৪টি চেইন একই posterior এ পৌঁছেছে — **PASS ✅**  
> **n_eff:** যদি > ৪০০ হয়, তাহলে আমাদের কাছে যথেষ্ট স্বাধীন sample আছে — **PASS ✅**  
> আমাদের সব সাইটে R̂ ≤ ১.০০৩ এবং n_eff ≥ ১২০০ ছিল — অত্যন্ত ভালো ফলাফল।

---

### অংশ ৬: Posterior Sample বের করা (Line 252-257)
```r
post       <- rstan::extract(fit)   # সব ৪০০০ sample বের করা হচ্ছে
beta_post  <- post$beta             # β এর ৪০০০ মান
gamma_post <- post$gamma            # γ এর ৪০০০ মান
rho_post   <- post$rho             # ρ এর ৪০০০ মান
R0_post    <- post$R0              # R₀ এর ৪০০০ মান
```
> **বাংলা:** এই ৪০০০ টি মান হলো আমাদের **posterior distribution**। এগুলো দিয়ে আমরা মেডিয়ান, ৯৫% CI সব হিসাব করি।

---

### অংশ ৭: R² হিসাব (Line 313-316)
```r
ss_res <- sum((cases - mu_med)^2)    # Residual Sum of Squares
ss_tot <- sum((cases - mean(cases))^2) # Total Sum of Squares
R2     <- 1 - ss_res / ss_tot        # R² = মডেলের ব্যাখ্যা করার ক্ষমতা
```
> **বাংলা:** R² মান বলে দেয় আমাদের মডেল ডেটার কত শতাংশ ব্যাখ্যা করতে পারে।  
> ILI এর জন্য R² = ০.৮৪ → মডেল ৮৪% ব্যাখ্যা করতে পারে।

---

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# ধাপ G: ১০টি প্লট বিস্তারিত ব্যাখ্যা
# STEP G: All 10 Plots Explained
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

**📂 ফোল্ডার খুলুন:** `plots/ILI/` বা `plots/Severe/` বা `plots/SARI/`

---

## 🖼️ Plot 01: Raw Data (কাঁচা ডেটা)
**📂 ফাইল:** `plots/ILI/01_raw_data.png`

```r
# R কোড (line 140-148):
p01 <- ggplot(dat, aes(x = date, y = cases)) +
  geom_point(shape = 19, size = 2.5, color = "black") +  # ← কালো বিন্দু
  geom_line(linewidth = 0.4, color = "gray40") +          # ← সংযোগকারী রেখা
  labs(title = "[ILI] Weekly Case Counts", x = "Date", y = "Reported cases")
```

> **🇧🇩 বাংলা:** "স্যার, এই প্লটে ২০২৫ সালের প্রতিটি সপ্তাহে কতটি কেস রিপোর্ট হয়েছে তা দেখানো হয়েছে। ILI তে জানুয়ারি থেকে শুরু হয়ে আগস্টে সর্বোচ্চ ১৫২ কেসে পৌঁছে এবং ডিসেম্বরে প্রায় শূন্যে নেমে আসে। এই ক্লাসিক ঘণ্টা আকৃতির curve দেখে আমরা সিদ্ধান্ত নিয়েছি যে একটি standard SIR মডেল এখানে উপযুক্ত।"

---

## 🖼️ Plot 02: SIR Compartments (S, I, R ট্র্যাজেক্টরি)
**📂 ফাইল:** `plots/ILI/02_sir_compartments.png`

```r
# ODE সমাধান করে গ্রাফ আঁকা হচ্ছে
sol <- solve_sir(beta_med, gamma_med, n_wks, N)
comp_long <- sol_df %>% pivot_longer(cols = c(Susceptible, Infectious, Recovered))

p02 <- ggplot() +
  geom_line(data = comp_long, aes(x = week, y = Count, linetype = Compartment)) +
  geom_point(...)  # ← আসল কেস ডেটা ওপরে দেখানো হচ্ছে
```

> **🇧🇩 বাংলা:** "এই প্লটে তিনটি আলাদা রেখা দেখাচ্ছে — S (যারা আক্রান্ত হতে পারে), I (এখন আক্রান্ত), R (সেরে উঠেছে)। S রেখাটি নিচে নামছে কারণ মানুষ আক্রান্ত হচ্ছে। I রেখাটি উপরে উঠে নামছে — এটিই epidemic wave। R রেখাটি ক্রমাগত বাড়ছে। গুরুত্বপূর্ণ বিষয়: মডেল বলছে মোট জনগোষ্ঠীর ~৫৪% আসলে আক্রান্ত হয়েছে, কিন্তু মাত্র ~১৫% রিপোর্ট হয়েছে।"

---

## 🖼️ Plot 03: Prior Predictive Check (Prior বৈধতা)
**📂 ফাইল:** `plots/ILI/03_prior_predictive.png`

```r
# ২০০ বার random parameter দিয়ে epidemic simulate করা হচ্ছে
prior_curves <- lapply(seq_len(200), function(i) {
  R0_s    <- rlnorm(1, log(1.7),  0.25)    # ← random R₀ নেওয়া হচ্ছে
  gamma_s <- rlnorm(1, log(0.44), 0.3)     # ← random γ নেওয়া হচ্ছে
  beta_s  <- R0_s * gamma_s                # ← β হিসাব করা হচ্ছে
  sol_s   <- solve_sir(beta_s, gamma_s, ...)  # ← ODE সমাধান
  mu_s    <- rho_s * sol_s[, "I"]           # ← প্রত্যাশিত কেস
})
```

> **🇧🇩 বাংলা:** "এই প্লটটি প্রমাণ করে যে ডেটা দেখার আগে আমাদের Prior অনুমান সঠিক ছিল। ধূসর রেখাগুলো হলো ২০০টি random simulation, যেখানে আমরা শুধু Prior থেকে parameter নিয়েছি। কালো বিন্দু (আসল ডেটা) এই মেঘের মাঝখানে আছে — প্রমাণ করে যে আমাদের Prior অসম্ভব কোনো epidemic forced করছে না।"

---

## 🖼️ Plot 04: Trace Plots (MCMC কনভার্জেন্স)
**📂 ফাইল:** `plots/ILI/04_trace_plots.png`

```r
posterior_arr <- as.array(fit, pars = c("R0", "gamma", "rho", "phi"))
p04 <- mcmc_trace(posterior_arr,   # ← bayesplot প্যাকেজ থেকে
                  facet_args = list(nrow = 4))
```

> **🇧🇩 বাংলা:** "এই প্লটটি আমাদের সবচেয়ে গুরুত্বপূর্ণ algorithmic proof। প্রতিটি রঙিন রেখা হলো একটি আলাদা MCMC চেইন। যদি সব চেইন একসাথে মিশে 'লোমশ শুঁয়োপোকার' মতো দেখায়, মানে MCMC সফলভাবে কাজ করেছে। আমাদের সব parameter এ R̂ ≤ ১.০০৩ — যা আদর্শ মানের খুব কাছে।"

---

## 🖼️ Plot 05: Posterior Histograms (parameter অনুমান)
**📂 ফাইল:** `plots/ILI/05_posterior_histograms.png`

```r
make_hist <- function(x, lab) {
  med <- median(x); lo <- quantile(x, 0.025); hi <- quantile(x, 0.975)
  ggplot(data.frame(v = x), aes(x = v)) +
    geom_histogram(fill = "gray70", color = "black", bins = 40) +
    geom_vline(xintercept = med, linetype = "solid") +   # ← মেডিয়ান রেখা
    geom_vline(xintercept = lo,  linetype = "dashed") +  # ← ২.৫% রেখা
    geom_vline(xintercept = hi,  linetype = "dashed")    # ← ৯৭.৫% রেখা
}
```

> **🇧🇩 বাংলা:** "প্রতিটি histogram হলো একটি parameter এর ৪০০০ MCMC sample এর বিতরণ। শক্ত রেখাটি মেডিয়ান (সেরা অনুমান), আর ড্যাশ রেখা দুটি ৯৫% Credible Interval এর সীমা। উদাহরণ: ILI তে β মেডিয়ান = ০.৯২, ৯৫% CrI = [০.৭০, ১.২২]।"

---

## 🖼️ Plot 06: Pair Plots (parameter এর মধ্যে সম্পর্ক)
**📂 ফাইল:** `plots/ILI/06_pair_plots.png`

```r
pair_df <- data.frame(beta = beta_post, gamma = gamma_post,
                      rho  = rho_post,  phi   = phi_post)
p06 <- ggpairs(pair_df,
               upper = list(continuous = wrap("cor", ...)),  # ← correlation
               lower = list(continuous = wrap("points", ...))) # ← scatter plot
```

> **🇧🇩 বাংলা:** "এই প্লট দেখায় parameter গুলো একে অপরের সাথে কতটা সম্পর্কযুক্ত। β এবং γ এর মধ্যে r ≈ ০.৯৯ — প্রায় পূর্ণ correlation। এর কারণ হলো সাপ্তাহিক কেস ডেটা শুধু R₀ = β/γ অনুপাতটি ভালোভাবে নির্ধারণ করতে পারে, কিন্তু β এবং γ আলাদাভাবে নির্ধারণ করতে পারে না। এই ঘটনাটিকে 'structural non-identifiability' বলে।"

---

## 🖼️ Plot 07: Posterior Predictive Check (চূড়ান্ত ফিট)
**📂 ফাইল:** `plots/ILI/07_posterior_predictive.png`

```r
mu_med  <- apply(mu_mat, 2, median)     # ← প্রতি সপ্তাহে মেডিয়ান প্রত্যাশা
mu_lo90 <- apply(mu_mat, 2, quantile, 0.05)  # ← ৯০% CI এর নিচের সীমা
mu_hi90 <- apply(mu_mat, 2, quantile, 0.95)  # ← ৯০% CI এর উপরের সীমা

p07 <- ggplot(pp_df, aes(x = week)) +
  geom_ribbon(aes(ymin = lo90, ymax = hi90), fill = "#CCCCCC") + # ← ৯০% CI
  geom_ribbon(aes(ymin = lo50, ymax = hi50), fill = "#888888") + # ← ৫০% CI
  geom_line(aes(y = med), linewidth = 0.9) +    # ← মেডিয়ান রেখা
  geom_point(aes(y = cases), shape = 19, size = 2.5) # ← আসল ডেটা
```

> **🇧🇩 বাংলা:** "এটি সবচেয়ে গুরুত্বপূর্ণ প্লট। কালো রেখা হলো মডেলের সেরা অনুমান। ছায়া এলাকা হলো অনিশ্চয়তার পরিসীমা। কালো বিন্দু হলো আসল রিপোর্টেড কেস। ILI এর জন্য R² = ০.৮৪ — মডেল ৮৪% ব্যাখ্যা করতে পারে। একটি সীমাবদ্ধতা: মডেল peak এর সময় একটু কম অনুমান করে কারণ SIR মডেল সবসময় একটি smooth curve তৈরি করে, কিন্তু আসল দুনিয়ায় super-spreading ঘটনায় তীক্ষ্ণ spike হয়।"

---

## 🖼️ Plot 08: Simulated Observations (সিমুলেটেড কেস)
**📂 ফাইল:** `plots/ILI/08_posterior_predictive_sim.png`

```r
yrep_mat <- post$cases_rep        # ← Stan থেকে সিমুলেটেড কেস বের করা
draw_idx <- sample(nrow(yrep_mat), 200)  # ← ২০০টি random draw নেওয়া

p08 <- ggplot() +
  geom_line(data = yrep_long, aes(x = week, y = y, group = draw),
            color = "gray70", alpha = 0.2) +   # ← ২০০ সিমুলেশন রেখা
  geom_point(...)   # ← আসল ডেটা উপরে
```

> **🇧🇩 বাংলা:** "Plot 07 এ আমরা শুধু মেডিয়ান দেখেছি। এই প্লটে ২০০টি সম্পূর্ণ সিমুলেশন দেখানো হয়েছে, যেখানে parameter uncertainty এবং Negative Binomial noise উভয়ই অন্তর্ভুক্ত। আসল ডেটার বিন্দুগুলো সম্পূর্ণভাবে ধূসর মেঘের মধ্যে আছে — প্রমাণ করে যে মডেলের noise model (φ ≈ ৪.৪) সঠিক।"

---

## 🖼️ Plot 09: R₀ Distribution (মূল প্রজনন সংখ্যা)
**📂 ফাইল:** `plots/ILI/09_R0_distribution.png`

```r
p09 <- ggplot(data.frame(R0 = R0_post), aes(x = R0)) +
  geom_histogram(fill = "gray70", bins = 50) +
  geom_vline(xintercept = 1.0, linetype = "dashed") +       # ← epidemic threshold
  geom_vline(xintercept = median(R0_post), linetype = "solid") # ← মেডিয়ান
```

> **🇧🇩 বাংলা:** "এই প্লটটি সবচেয়ে গুরুত্বপূর্ণ epidemiological ফলাফল। ড্যাশ রেখাটি হলো ১.০ — epidemic threshold। যদি R₀ > ১ হয়, epidemic বাড়বে; < ১ হলে কমবে। আমাদের ILI সাইটে R₀ এর মেডিয়ান = ১.৬৫, এবং পুরো posterior distribution ১.০ এর উপরে। এর মানে আমরা ১০০% নিশ্চিত যে এটি একটি সক্রিয় epidemic ছিল।"

---

## 🖼️ Plot 10: Residuals (মডেলের ত্রুটি বিশ্লেষণ)
**📂 ফাইল:** `plots/ILI/10_residuals.png`

```r
resid_df <- data.frame(
  week     = 1:n_wks,
  residual = cases - mu_med,   # ← আসল - প্রত্যাশিত
  fitted   = mu_med
)

r10a <- ggplot(resid_df, aes(x = week, y = residual)) +
  geom_hline(yintercept = 0, linetype = "dashed") +  # ← শূন্য রেখা
  geom_point(shape = 19, size = 2)

r10b <- ggplot(resid_df, aes(x = fitted, y = residual)) + ...
```

> **🇧🇩 বাংলা:** "Residual = আসল কেস বিয়োগ মডেলের অনুমান। একজন ভালো গবেষক সবসময় এই প্লট দেখান কারণ এটি মডেলের দুর্বলতা প্রকাশ করে। Peak এর সময়ে (week 25-31) residual সবচেয়ে বড় — মানে মডেল সেই সময়ে একটু কম অনুমান করেছে। এটি SIR মডেলের একটি পরিচিত সীমাবদ্ধতা।"

---

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# ধাপ H: তিন সাইটের তুলনামূলক বিশ্লেষণ
# STEP H: Three-Site Comparative Analysis
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

### 🇧🇩 বাংলায় বলুন:
> "স্যার, এখন আমি তিনটি সাইটের ফলাফল তুলনা করব।"

| পরিমাপ | ILI (বাইরের রোগী) | Severe (হাসপাতাল) | SARI (বিশেষজ্ঞ) |
|--------|-------------------|-------------------|-----------------|
| **ডেটার আকৃতি** | ঘণ্টা আকৃতি ✅ | ঘণ্টা আকৃতি ✅ | সমতল ও উঁচু ❌ |
| **মোট কেস** | ১,২১৪ | ১,৫৭০ | ১২,৭২১ |
| **মডেল ফিট R²** | **০.৮৪** ✅ | **০.৮৬** ✅ | **০.০৮** ❌ |
| **R₀** | ১.৬৫ | ১.৬৬ | ১.৪১ |
| **রিপোর্টিং রেট ρ** | ১৪.৭% | ১৮.৮% | ৬৪.৬% |
| **সংক্রমণকাল** | ~১২.৫ দিন | ~১১.১ দিন | ~৩১.৯ দিন |

### ILI ও Severe একই R₀ কেন?
> "স্যার, ILI (১.৬৫) এবং Severe (১.৬৬) এর R₀ প্রায় একই — এটাই প্রমাণ করে যে দুটি সাইট একই epidemic wave ধারণ করছে। শুধু রোগীর তীব্রতা আলাদা।"

### SARI তে R² মাত্র ০.০৮ কেন?
> "স্যার, SARI সাইটটি জানুয়ারি থেকেই ২১৯ কেস নিয়ে শুরু হয়েছে। একটি standard SIR মডেল ধরে নেয় শুরুতে প্রায় শূন্য কেস থাকে এবং ধীরে বাড়ে। কিন্তু SARI তে এটি সত্য না — এটি একটি endemic baseline, epidemic wave না। তাই মডেল SARI তে ভালো fit করতে পারেনি। এটি মডেলের ব্যর্থতা না — এটি Bayesian diagnostics এর সাফল্য যে সে এই সমস্যা চিহ্নিত করতে পেরেছে।"

---

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# ধাপ I: চূড়ান্ত ফলাফল সারসংক্ষেপ
# STEP I: Final Results Summary
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

### 🇧🇩 মূল আবিষ্কারগুলো বলুন:

> "স্যার, আমাদের গবেষণার মূল চারটি আবিষ্কার:

> **১. R₀ = ১.৬৫ (ILI) ও ১.৬৬ (Severe):**  
> এই রোগটি মাঝারি মাত্রার সংক্রামক, seasonal flu এর মতো। প্রতিটি আক্রান্ত ব্যক্তি গড়ে ১.৬৫ জনকে আক্রান্ত করেছে।

> **২. Under-reporting আবিষ্কার:**  
> ILI তে মাত্র ১৪.৭% এবং Severe তে ১৮.৮% আসল কেস রিপোর্ট হয়েছে। মানে আসল সংক্রমণ রিপোর্টের চেয়ে ৫-৭ গুণ বেশি ছিল।

> **৩. SARI সাইটের বিশেষ সমস্যা:**  
> SARI তে ৬৪.৬% রিপোর্টিং রেট এবং R² মাত্র ০.০৮ — এটি প্রমাণ করে যে SARI একটি epidemic নয়, এটি endemic baseline। ভবিষ্যতে SIRS বা SEIR মডেল ব্যবহার করা উচিত।

> **৪. উৎকৃষ্ট MCMC convergence:**  
> সব সাইটে R̂ ≤ ১.০০৩ এবং n_eff ≥ ১২০০ — প্রমাণ করে আমাদের statistical estimation নির্ভরযোগ্য।"

---

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# ধাপ J: প্রফেসরের সম্ভাব্য প্রশ্ন ও উত্তর
# STEP J: Expected Professor Questions & Answers
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

---

### ❓ প্রশ্ন ১: "তুমি β এর মান কোথা থেকে পেলে?"
> **🇧🇩 উত্তর:** "স্যার, আমরা β ম্যানুয়ালি বেছে নিইনি। আমরা Stan MCMC অ্যালগরিদম ব্যবহার করেছি যা ৪০০০ বিভিন্ন parameter combination চেষ্টা করেছে এবং যেগুলো আসল ডেটার সাথে সবচেয়ে ভালো মেলে সেগুলো রেখে দিয়েছে। ডেটা নিজেই আমাদের বলেছে β = ০.৯২।"

---

### ❓ প্রশ্ন ২: "Negative Binomial কেন? Normal distribution কেন না?"
> **🇧🇩 উত্তর:** "স্যার, সাপ্তাহিক কেস ডেটায় সাধারণ Normal distribution ব্যবহার করলে সমস্যা হয় কারণ কেস কখনো negative হতে পারে না। আরও গুরুত্বপূর্ণভাবে, epidemic ডেটায় variance সাধারণত mean এর চেয়ে অনেক বেশি (overdispersion)। Negative Binomial তে এই extra variability parameter φ দিয়ে নিয়ন্ত্রণ করা হয়।"

---

### ❓ প্রশ্ন ৩: "৪টি chain কেন? একটি দিলে হতো না?"
> **🇧🇩 উত্তর:** "স্যার, যদি আমরা একটি মাত্র chain চালাই, আমরা কখনো নিশ্চিত হতে পারব না যে সেটি সত্যিকারের posterior distribution খুঁজে পেয়েছে নাকি শুধু একটি local mode তে আটকে গেছে। ৪টি আলাদা শুরু থেকে ৪টি chain চালালে এবং সব মিলে গেলে আমরা নিশ্চিত হতে পারি।"

---

### ❓ প্রশ্ন ৪: "SARI তে R² এত কম কেন? মডেল ব্যর্থ হয়েছে?"
> **🇧🇩 উত্তর:** "স্যার, না — মডেল ব্যর্থ হয়নি, মডেল সঠিকভাবে দেখিয়েছে যে SARI ডেটা একটি single-wave epidemic এর সাথে মেলে না। SARI তে Week 1 থেকেই ২১৯ কেস — এটি endemic baseline। একটি standard SIR মডেল এই ধরনের ডেটার জন্য প্রযোজ্য নয়। এটি আমাদের recommendation: SARI এর জন্য SIRS বা SEIR মডেল ব্যবহার করা উচিত।"

---

### ❓ প্রশ্ন ৫: "R-hat কী? কেন ১.০৫ এর নিচে থাকতে হবে?"
> **🇧🇩 উত্তর:** "স্যার, R-hat হলো Gelman-Rubin convergence diagnostic। এটি মাপে ৪টি chain এর মধ্যে কতটা পার্থক্য আছে। যদি R̂ = ১.০০০ হয় মানে সব chain পুরোপুরি একমত। যদি ১.০৫ এর উপরে যায় মানে chain গুলো এখনো সঠিক posterior খুঁজে পায়নি। আমাদের সব R̂ ≤ ১.০০৩ ছিল — প্রায় আদর্শ।"

---

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# ধাপ K: উপসংহার বলার স্ক্রিপ্ট
# STEP K: Closing Statement
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

### 🇧🇩 বাংলায় শেষ করুন:
> "স্যার, সংক্ষেপে বলতে গেলে:  
>
> আমরা একটি Bayesian SIR মডেল তৈরি করেছি এবং তিনটি sentinel সাইটের ডেটায় প্রয়োগ করেছি।  
>
> ILI এবং Severe সাইটে মডেল চমৎকারভাবে কাজ করেছে (R² = ০.৮৪ এবং ০.৮৬)। আমরা প্রমাণ করেছি যে R₀ ≈ ১.৬৫, অর্থাৎ প্রতিটি আক্রান্ত ব্যক্তি গড়ে আরও ১.৬৫ জনকে আক্রান্ত করেছে। আসল সংক্রমণ রিপোর্টের ৫-৭ গুণ বেশি ছিল।  
>
> SARI সাইটে মডেল fit কম হয়েছে (R² = ০.০৮), কিন্তু এটি আমাদের একটি গুরুত্বপূর্ণ খোঁজ — SARI একটি endemic surveillance baseline, একটি isolated epidemic নয়।  
>
> সব মিলিয়ে এই গবেষণা দেখায় যে Bayesian পদ্ধতি শুধু parameters অনুমান করে না, বরং মডেলের সীমাবদ্ধতাও স্বয়ংক্রিয়ভাবে চিহ্নিত করে। ধন্যবাদ।"

### 🇬🇧 English closing:
> "In summary, Professor, we successfully fit Bayesian SIR models to three sentinel sites. ILI and Severe achieved excellent fits (R² ≥ 0.84) with R₀ ≈ 1.65, confirming a moderately transmissible seasonal pathogen. The SARI site's low fit (R² = 0.08) is not a failure — it is a diagnostic success revealing that SARI represents an endemic baseline requiring a more complex SIRS or SEIR framework. All MCMC diagnostics were excellent (R̂ ≤ 1.003, n_eff ≥ 1200). Thank you."

---

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# দ্রুত রেফারেন্স কার্ড
# Quick Reference Card
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

| শব্দ | বাংলা অর্থ |
|------|------------|
| SIR Model | তিনটি compartment এ জনগোষ্ঠীকে ভাগ করা মহামারি মডেল |
| β (beta) | সংক্রমণ রেট — কত দ্রুত রোগ ছড়ায় |
| γ (gamma) | রিকভারি রেট — কত দ্রুত সেরে ওঠে |
| R₀ | মূল প্রজনন সংখ্যা = β/γ |
| ρ (rho) | রিপোর্টিং হার — কত ভাগ কেস ধরা পড়ে |
| φ (phi) | Overdispersion — ডেটার extra variability |
| Prior | ডেটা দেখার আগের অনুমান |
| Posterior | ডেটা দেখার পরের আপডেট করা অনুমান |
| MCMC | Markov Chain Monte Carlo — parameter খোঁজার অ্যালগরিদম |
| Credible Interval | সম্ভাবনার পরিসীমা (Bayesian version of confidence interval) |
| R̂ (R-hat) | MCMC convergence এর প্রমাণ (< ১.০৫ হওয়া উচিত) |
| n_eff | কার্যকর sample সংখ্যা (> ৪০০ হওয়া উচিত) |
| Negative Binomial | Count ডেটার জন্য overdispersion সহ probability distribution |
| ILI | Influenza-Like Illness — বাইরের ক্লিনিকের রোগী |
| Severe | হাসপাতালে ভর্তি গুরুতর রোগী |
| SARI | Severe Acute Respiratory Infection — বিশেষজ্ঞ হাসপাতালের রোগী |
| R² | মডেলের ব্যাখ্যা ক্ষমতা (১.০ = নিখুঁত, ০.৮৪ = ভালো) |
| Endemic | রোগ সবসময় বিদ্যমান (কমে না) |
| Epidemic | রোগ হঠাৎ বেড়ে আবার কমে |
