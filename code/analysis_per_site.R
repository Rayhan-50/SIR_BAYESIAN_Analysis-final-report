# ============================================================
# Bayesian SIR Model — Per Sentinel Site Analysis
# Raihan Research Group | 2026
# Sites: ILI (outpatient), Severe (hospital), SARI
# ============================================================
# Run: Rscript code/analysis_per_site.R
# Outputs: plots/<site>/ and outputs/results_<site>.json
# ============================================================

suppressPackageStartupMessages({
  library(rstan)
  library(loo)
  library(deSolve)
  library(ggplot2)
  library(bayesplot)
  library(gridExtra)
  library(jsonlite)
  library(dplyr)
  library(tidyr)
  library(GGally)
})

rstan_options(auto_write = TRUE)
options(mc.cores = min(parallel::detectCores(), 4))

# ── Project paths ─────────────────────────────────────────────────────────────
# Use script location if available, else fall back to working directory
proj_dir <- tryCatch(
  normalizePath(file.path(dirname(sys.frame(1)$ofile), ".."), mustWork = TRUE),
  error = function(e) getwd()
)
if (!dir.exists(proj_dir)) proj_dir <- getwd()
cat(sprintf("Project dir: %s\n", proj_dir))

stan_file <- file.path(proj_dir, "code", "sir_model.stan")
out_dir   <- file.path(proj_dir, "outputs")
dir.create(out_dir, showWarnings = FALSE)

# ── Site definitions ──────────────────────────────────────────────────────────
# Each site: label, CSV file, total catchment population (N)
sites <- list(
  list(
    label    = "ILI",
    csv      = file.path(proj_dir, "data", "site_ILI.csv"),
    pop      = 11750L,
    desc     = "Influenza-Like Illness (Outpatient Sentinel)"
  ),
  list(
    label    = "Severe",
    csv      = file.path(proj_dir, "data", "site_Severe.csv"),
    pop      = 11750L,
    desc     = "Severe/Hospital Sentinel"
  )
)

# ── Black-and-white theme ─────────────────────────────────────────────────────
bw_theme <- theme_minimal(base_family = "serif") +
  theme(
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    text             = element_text(color = "black"),
    axis.text        = element_text(color = "black"),
    axis.title       = element_text(color = "black"),
    plot.title       = element_text(color = "black", face = "bold"),
    plot.subtitle    = element_text(color = "gray30", size = 9),
    strip.text       = element_text(color = "black")
  )

# ── SIR ODE ───────────────────────────────────────────────────────────────────
sir_ode <- function(t, y, parms) {
  S <- y[1]; I <- y[2]; R <- y[3]
  beta  <- parms["beta"]
  gamma <- parms["gamma"]
  N     <- parms["N"]
  list(c(-beta * S * I / N,
          beta * S * I / N - gamma * I,
          gamma * I))
}

solve_sir <- function(beta, gamma, n_wks, pop) {
  parms <- c(beta = beta, gamma = gamma, N = pop)
  y0    <- c(S = pop - 1, I = 1, R = 0)
  times <- seq(0, n_wks, by = 1)
  ode(y = y0, times = times, func = sir_ode, parms = parms, method = "rk4")
}

# ── Helper: save plot ─────────────────────────────────────────────────────────
save_plot <- function(p, path, w = 7, h = 5) {
  ggsave(path, plot = p, width = w, height = h, dpi = 200, bg = "white")
  invisible(p)
}

# ── Compile Stan model ONCE ───────────────────────────────────────────────────
cat("\n=== Compiling Stan model ===\n")
stan_model_obj <- stan_model(file = stan_file)
cat("Stan model compiled.\n\n")

# ══════════════════════════════════════════════════════════════════════════════
# MAIN LOOP — one iteration per sentinel site
# ══════════════════════════════════════════════════════════════════════════════
all_results <- list()

for (site in sites) {

  site_label <- site$label
  site_desc  <- site$desc
  N          <- site$pop

  cat(sprintf("\n\n══════════════════════════════════════════\n"))
  cat(sprintf("  SITE: %s — %s\n", site_label, site_desc))
  cat(sprintf("══════════════════════════════════════════\n\n"))

  # Create site-specific plots directory
  plots_dir <- file.path(proj_dir, "plots", site_label)
  dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

  # ── 1. Load data ────────────────────────────────────────────────────────────
  dat   <- read.csv(site$csv, stringsAsFactors = FALSE)
  dat$date <- as.Date(dat$date)
  n_wks <- nrow(dat)
  cases <- dat$cases

  # Skip site if all zeros
  if (sum(cases) == 0) {
    cat(sprintf("SKIP: %s has no cases.\n", site_label))
    next
  }

  # ── Data quality checks ─────────────────────────────────────────────────────
  # 1. Stair-step artifact check (monthly data divided by 4)
  rle_obj <- rle(cases)
  max_run <- max(rle_obj$lengths[rle_obj$values > 0])
  if (max_run >= 3) {
    cat(sprintf("*** WARNING [%s]: Stair-step artifact detected.\n", site_label))
    cat(sprintf("    Longest run of identical nonzero values: %d weeks.\n", max_run))
    cat("    This may indicate monthly data interpolated to weekly.\n")
    cat("    Acknowledge in manuscript Limitations section.\n\n")
  }

  cat(sprintf("Weeks       : %d\n",      n_wks))
  cat(sprintf("Total cases : %d\n",      sum(cases)))
  cat(sprintf("Peak cases  : %d (week %d)\n", max(cases), which.max(cases) - 1))
  cat(sprintf("Population N: %d\n\n",    N))

  # ── Plot 01: Raw case counts ─────────────────────────────────────────────────
  p01 <- ggplot(dat, aes(x = date, y = cases)) +
    geom_point(shape = 19, size = 2.5, color = "black") +
    geom_line(linewidth = 0.4, color = "gray40") +
    labs(
      title    = sprintf("[%s] Weekly Case Counts", site_label),
      subtitle = site_desc,
      x = "Date", y = "Reported cases"
    ) + bw_theme
  save_plot(p01, file.path(plots_dir, "01_raw_data.png"))
  cat("Plot 01 saved.\n")

  # Plot 02 is generated after Stan fitting using posterior medians (see below)

  # ── Plot 03: Prior predictive check ──────────────────────────────────────────
  # Priors MUST match sir_model.stan exactly:
  #   R0    ~ lognormal(log(2.0), 0.3)
  #   gamma ~ lognormal(log(0.5), 0.3)
  #   rho   ~ beta(2, 10)
  #   lambda ~ exponential(0.5)
  # Observation model: mu(t) = rho * incidence(t) + lambda
  #   where incidence(t) = S(t-1) - S(t)  [new infections per week]
  set.seed(42)
  n_prior <- 200
  prior_curves <- lapply(seq_len(n_prior), function(i) {
    R0_s     <- rlnorm(1, log(2.0), 0.3)
    gamma_s  <- rlnorm(1, log(0.5), 0.3)
    rho_s    <- rbeta(1, 2, 10)
    lambda_s <- rexp(1, 0.5)
    beta_s   <- R0_s * gamma_s
    tryCatch({
      sol_s  <- solve_sir(beta_s, gamma_s, n_wks, N)
      sol_df_s <- as.data.frame(sol_s)
      S_vec  <- sol_df_s$S
      inc    <- pmax(c(S_vec[1] - S_vec[2], -diff(S_vec[-1])), 1e-6)
      mu_s   <- pmax(rho_s * inc + lambda_s, 1e-6)
      data.frame(week = 1:n_wks, mu = mu_s, draw = i)
    }, error = function(e) NULL)
  })
  prior_df <- bind_rows(Filter(Negate(is.null), prior_curves))

  p03 <- ggplot() +
    geom_line(data = prior_df, aes(x = week, y = mu, group = draw),
              color = "gray70", alpha = 0.3, linewidth = 0.3) +
    geom_point(data = data.frame(week = 1:n_wks, cases = cases),
               aes(x = week, y = cases), shape = 19, size = 2, color = "black") +
    coord_cartesian(ylim = c(0, max(cases) * 5)) +
    labs(
      title    = sprintf("[%s] Prior Predictive Check", site_label),
      subtitle = "mu(t) = rho * incidence(t) + lambda  [incidence = S(t-1) - S(t)]",
      x = "Week", y = "Expected cases"
    ) + bw_theme
  save_plot(p03, file.path(plots_dir, "03_prior_predictive.png"))
  cat("Plot 03 saved.\n")

  # ── 4. Fit Stan model ─────────────────────────────────────────────────────────
  cat(sprintf("\n=== Fitting Stan model for %s ===\n", site_label))
  baseline_priors <- list(
    R0_mu       = log(2.0), R0_sigma     = 0.3,
    gamma_mu    = log(0.5), gamma_sigma  = 0.3,
    rho_alpha   = 2,        rho_beta     = 10,
    phi_rate    = 0.5,
    I0_mu       = log(1),   I0_sigma     = 2.0,
    lambda_rate = 0.5
  )

  stan_data <- c(
    list(N_weeks = n_wks, cases = cases, pop = as.double(N)),
    baseline_priors
  )

  fit <- sampling(
    stan_model_obj,
    data    = stan_data,
    chains  = 4,
    iter    = 3000,
    warmup  = 1500,
    seed    = 12345,
    control = list(adapt_delta = 0.98, max_treedepth = 15),
    refresh = 200
  )

  # ── Convergence diagnostics ───────────────────────────────────────────────────
  cat(sprintf("\n=== Convergence Diagnostics [%s] ===\n", site_label))
  sum_fit <- summary(fit)$summary
  params  <- c("R0", "gamma", "rho", "phi", "beta", "I0", "lambda")
  diag_df <- sum_fit[params, c("mean","sd","2.5%","50%","97.5%","n_eff","Rhat")]
  print(round(diag_df, 4))

  rhat_vals <- diag_df[,"Rhat"]
  neff_vals <- diag_df[,"n_eff"]

  cat("\n--- R-hat check (threshold 1.01, Vehtari et al. 2021) ---\n")
  if (all(rhat_vals < 1.01, na.rm = TRUE)) {
    cat("PASS: All R-hat < 1.01\n")
  } else if (all(rhat_vals < 1.05, na.rm = TRUE)) {
    warning(sprintf("[%s] WARN: R-hat in [1.01, 1.05) — consider more iterations", site_label))
    print(rhat_vals[rhat_vals >= 1.01])
  } else {
    warning(sprintf("[%s] FAIL: Some R-hat >= 1.05", site_label))
    print(rhat_vals[rhat_vals >= 1.05])
  }

  cat("--- ESS check (min 1000 recommended) ---\n")
  if (all(neff_vals > 1000, na.rm = TRUE)) {
    cat("PASS: All n_eff > 1000\n")
  } else if (all(neff_vals > 400, na.rm = TRUE)) {
    cat("WARN: Some n_eff in [400, 1000):\n")
    print(neff_vals[neff_vals <= 1000])
  } else {
    cat("FAIL: Some n_eff <= 400:\n")
    print(neff_vals[neff_vals <= 400])
  }

  n_div <- get_num_divergent(fit)
  if (n_div > 0) warning(sprintf("[%s] DIVERGENCES: %d divergent transitions!", site_label, n_div))

  # ── 5. Extract posterior samples ──────────────────────────────────────────────
  post        <- rstan::extract(fit)
  beta_post   <- post$beta
  gamma_post  <- post$gamma
  rho_post    <- post$rho
  phi_post    <- post$phi
  R0_post     <- post$R0
  I0_post     <- post$I0
  lambda_post <- post$lambda

  beta_med   <- median(beta_post)
  gamma_med  <- median(gamma_post)
  rho_med    <- median(rho_post)
  I0_med     <- median(I0_post)
  lambda_med <- median(lambda_post)

  # ── LOO-CV ────────────────────────────────────────────────────────────────────
  cat(sprintf("\n=== LOO Cross-Validation [%s] ===\n", site_label))
  log_lik_mat <- extract_log_lik(fit, parameter_name = "log_lik", merge_chains = FALSE)
  loo_res     <- loo(log_lik_mat, r_eff = relative_eff(exp(log_lik_mat)))
  print(loo_res)
  loo_elpd    <- loo_res$estimates["elpd_loo", "Estimate"]
  loo_se      <- loo_res$estimates["elpd_loo", "SE"]
  k_vals      <- loo_res$pointwise[, "influence_pareto_k"]
  n_bad_k     <- sum(k_vals > 0.7)
  if (n_bad_k > 0) cat(sprintf("WARNING: %d observations with Pareto k > 0.7\n", n_bad_k))

  # ── Plot 02: SIR compartments with posterior medians (incidence-based) ────────
  sol_p <- solve_sir(beta_med, gamma_med, n_wks, N)
  sol_pdf <- as.data.frame(sol_p) %>% filter(time >= 0)
  S_full <- sol_pdf$S
  incidence_fit <- pmax(c(S_full[1] - S_full[2], -diff(S_full[-1])), 1e-6)
  fitted_cases  <- rho_med * incidence_fit + lambda_med

  comp_long2 <- as.data.frame(sol_p) %>%
    filter(time >= 1) %>%
    mutate(week = time) %>%
    rename(Susceptible = S, Infectious = I, Recovered = R) %>%
    select(week, Susceptible, Infectious, Recovered) %>%
    pivot_longer(cols = c(Susceptible, Infectious, Recovered),
                 names_to = "Compartment", values_to = "Count")

  p02 <- ggplot() +
    geom_line(data = comp_long2,
              aes(x = week, y = Count, linetype = Compartment), linewidth = 0.8) +
    geom_line(data = data.frame(week = 1:n_wks, fitted = fitted_cases),
              aes(x = week, y = fitted * (N / max(fitted_cases))),
              linetype = "dotted", linewidth = 1.0, color = "black") +
    geom_point(data = data.frame(week = 1:n_wks,
                                 obs = cases * (N / max(fitted_cases))),
               aes(x = week, y = obs), shape = 1, size = 2, color = "black") +
    scale_y_continuous(
      name = "Compartment count",
      sec.axis = sec_axis(~ . / (N / max(fitted_cases)),
                          name = "Cases / fitted (rho * incidence + lambda)")
    ) +
    scale_linetype_manual(values = c(Susceptible = "solid",
                                     Infectious   = "dashed",
                                     Recovered    = "dotdash")) +
    labs(
      title    = sprintf("[%s] SIR Compartments (posterior medians)", site_label),
      subtitle = sprintf("R0=%.2f | beta=%.3f | gamma=%.3f | rho=%.3f | lambda=%.2f",
                         median(R0_post), beta_med, gamma_med, rho_med, lambda_med),
      x = "Week", linetype = ""
    ) + bw_theme + theme(legend.position = "bottom")
  save_plot(p02, file.path(plots_dir, "02_sir_compartments.png"), w = 8, h = 5)
  cat("Plot 02 saved.\n")

  # ── Plot 04: Trace plots ────────────────────────────────────────────────────
  posterior_arr <- as.array(fit, pars = c("R0","gamma","rho","phi"))
  color_scheme_set("gray")
  p04 <- mcmc_trace(posterior_arr,
                    pars = c("R0","gamma","rho","phi"),
                    facet_args = list(nrow = 4, labeller = label_parsed)) +
    labs(title = sprintf("[%s] MCMC Trace Plots", site_label)) + bw_theme
  save_plot(p04, file.path(plots_dir, "04_trace_plots.png"), w = 9, h = 7)
  cat("Plot 04 saved.\n")

  # ── Plot 05: Posterior histograms ────────────────────────────────────────────
  make_hist <- function(x, lab, xlab = lab) {
    med <- median(x); lo <- quantile(x, 0.025); hi <- quantile(x, 0.975)
    ggplot(data.frame(v = x), aes(x = v)) +
      geom_histogram(fill = "gray70", color = "black", bins = 40) +
      geom_vline(xintercept = med, linetype = "solid",  linewidth = 0.8) +
      geom_vline(xintercept = lo,  linetype = "dashed", linewidth = 0.6) +
      geom_vline(xintercept = hi,  linetype = "dashed", linewidth = 0.6) +
      labs(title = lab, x = xlab, y = "Count") + bw_theme
  }
  p05 <- grid.arrange(
    make_hist(beta_post,  expression(beta),  "beta"),
    make_hist(gamma_post, expression(gamma), "gamma"),
    make_hist(rho_post,   expression(rho),   "rho"),
    make_hist(phi_post,   expression(phi),   "phi"),
    make_hist(R0_post,    expression(R[0]),  "R0"),
    nrow = 2, ncol = 3,
    top = sprintf("[%s] Posterior Distributions (median + 95%% CI)", site_label)
  )
  ggsave(file.path(plots_dir, "05_posterior_histograms.png"),
         plot = p05, width = 12, height = 7, dpi = 200, bg = "white")
  cat("Plot 05 saved.\n")

  # ── Plot 06: Pair plots ──────────────────────────────────────────────────────
  pair_df  <- data.frame(beta = beta_post, gamma = gamma_post,
                         rho  = rho_post,  phi   = phi_post)
  idx      <- sample(nrow(pair_df), min(1000, nrow(pair_df)))
  p06 <- ggpairs(pair_df[idx, ],
                 upper = list(continuous = wrap("cor", size = 3, color = "black")),
                 lower = list(continuous = wrap("points", alpha = 0.15, size = 0.5, color = "black")),
                 diag  = list(continuous = wrap("densityDiag", fill = "gray70", color = "black"))) +
    labs(title = sprintf("[%s] Joint Posterior Distributions", site_label)) + bw_theme
  ggsave(file.path(plots_dir, "06_pair_plots.png"),
         plot = p06, width = 9, height = 8, dpi = 200, bg = "white")
  cat("Plot 06 saved.\n")

  # ── Plot 07: Posterior predictive check ─────────────────────────────────────
  mu_mat  <- post$mu
  mu_med  <- apply(mu_mat, 2, median)
  mu_lo50 <- apply(mu_mat, 2, quantile, 0.25)
  mu_hi50 <- apply(mu_mat, 2, quantile, 0.75)
  mu_lo90 <- apply(mu_mat, 2, quantile, 0.05)
  mu_hi90 <- apply(mu_mat, 2, quantile, 0.95)

  ss_res <- sum((cases - mu_med)^2)
  ss_tot <- sum((cases - mean(cases))^2)
  R2     <- 1 - ss_res / ss_tot
  cat(sprintf("\nPosterior predictive R² [%s] = %.3f\n", site_label, R2))

  pp_df <- data.frame(week = 1:n_wks, med = mu_med,
                      lo50 = mu_lo50, hi50 = mu_hi50,
                      lo90 = mu_lo90, hi90 = mu_hi90,
                      cases = cases)

  p07 <- ggplot(pp_df, aes(x = week)) +
    geom_ribbon(aes(ymin = lo90, ymax = hi90), fill = "#CCCCCC", alpha = 0.7) +
    geom_ribbon(aes(ymin = lo50, ymax = hi50), fill = "#888888", alpha = 0.7) +
    geom_line(aes(y = med), linewidth = 0.9, color = "black") +
    geom_point(aes(y = cases), shape = 19, size = 2.5, color = "black") +
    annotate("text", x = 5, y = max(cases) * 0.9,
             label = sprintf("R² = %.2f", R2), size = 4, family = "serif") +
    labs(
      title    = sprintf("[%s] Posterior Predictive Check", site_label),
      subtitle = "Dark gray: 50% CI | Light gray: 90% CI | Line: median",
      x = "Week", y = "Expected cases"
    ) + bw_theme
  save_plot(p07, file.path(plots_dir, "07_posterior_predictive.png"))
  cat("Plot 07 saved.\n")

  # ── Plot 08: Simulated observations ─────────────────────────────────────────
  yrep_mat  <- post$cases_rep
  draw_idx  <- sample(nrow(yrep_mat), min(200, nrow(yrep_mat)))
  yrep_df   <- as.data.frame(t(yrep_mat[draw_idx, ]))
  yrep_df$week <- 1:n_wks
  yrep_long <- pivot_longer(yrep_df, -week, names_to = "draw", values_to = "y")

  p08 <- ggplot() +
    geom_line(data = yrep_long, aes(x = week, y = y, group = draw),
              color = "gray70", alpha = 0.2, linewidth = 0.3) +
    geom_point(data = data.frame(week = 1:n_wks, cases = cases),
               aes(x = week, y = cases), shape = 19, size = 2.5, color = "black") +
    labs(
      title = sprintf("[%s] Posterior Predictive — Simulated Observations", site_label),
      x = "Week", y = "Cases"
    ) + bw_theme
  save_plot(p08, file.path(plots_dir, "08_posterior_predictive_sim.png"))
  cat("Plot 08 saved.\n")

  # ── Plot 09: R0 distribution ─────────────────────────────────────────────────
  p09 <- ggplot(data.frame(R0 = R0_post), aes(x = R0)) +
    geom_histogram(fill = "gray70", color = "black", bins = 50) +
    geom_vline(xintercept = 1.0,             linetype = "dashed", linewidth = 1.0) +
    geom_vline(xintercept = median(R0_post), linetype = "solid",  linewidth = 0.9) +
    annotate("text", x = 1.02, y = Inf, label = "Epidemic\nthreshold",
             hjust = 0, vjust = 1.5, size = 3.5, family = "serif") +
    annotate("text", x = median(R0_post) + 0.02, y = Inf,
             label = sprintf("Median = %.2f", median(R0_post)),
             hjust = 0, vjust = 3.0, size = 3.5, family = "serif") +
    labs(
      title = sprintf("[%s] Posterior Distribution of R0", site_label),
      x = expression(R[0]), y = "Count"
    ) + bw_theme
  save_plot(p09, file.path(plots_dir, "09_R0_distribution.png"))
  cat("Plot 09 saved.\n")

  # ── Plot 10: Residuals ────────────────────────────────────────────────────────
  resid_df <- data.frame(week = 1:n_wks,
                         residual = cases - mu_med,
                         fitted   = mu_med)
  r10a <- ggplot(resid_df, aes(x = week, y = residual)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_point(shape = 19, size = 2) +
    geom_line(linewidth = 0.4, color = "gray50") +
    labs(title = sprintf("[%s] Residuals Over Time", site_label),
         x = "Week", y = "Residual") + bw_theme
  r10b <- ggplot(resid_df, aes(x = fitted, y = residual)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_point(shape = 19, size = 2) +
    labs(title = "Residuals vs. Fitted", x = "Fitted (median)", y = "Residual") + bw_theme
  p10 <- grid.arrange(r10a, r10b, nrow = 1)
  ggsave(file.path(plots_dir, "10_residuals.png"),
         plot = p10, width = 10, height = 4.5, dpi = 200, bg = "white")
  cat("Plot 10 saved.\n")

  # ── Save results JSON ─────────────────────────────────────────────────────────
  q <- function(x, p) quantile(x, p, names = FALSE)
  results_list <- list(
    site        = site_label,
    description = site_desc,
    N           = N,
    n_weeks     = n_wks,
    total_cases = sum(cases),
    peak_cases  = max(cases),
    peak_week   = which.max(cases) - 1L,
    summary = list(
      beta  = list(mean = mean(beta_post),  sd = sd(beta_post),
                   q2.5 = q(beta_post,0.025), q50 = q(beta_post,0.5),
                   q97.5 = q(beta_post,0.975),
                   rhat = diag_df["beta","Rhat"], ess = diag_df["beta","n_eff"]),
      gamma = list(mean = mean(gamma_post), sd = sd(gamma_post),
                   q2.5 = q(gamma_post,0.025), q50 = q(gamma_post,0.5),
                   q97.5 = q(gamma_post,0.975),
                   rhat = diag_df["gamma","Rhat"], ess = diag_df["gamma","n_eff"]),
      rho   = list(mean = mean(rho_post),   sd = sd(rho_post),
                   q2.5 = q(rho_post,0.025), q50 = q(rho_post,0.5),
                   q97.5 = q(rho_post,0.975),
                   rhat = diag_df["rho","Rhat"], ess = diag_df["rho","n_eff"]),
      phi    = list(mean = mean(phi_post),    sd = sd(phi_post),
                    q2.5 = q(phi_post,0.025), q50 = q(phi_post,0.5),
                    q97.5 = q(phi_post,0.975),
                    rhat = diag_df["phi","Rhat"], ess = diag_df["phi","n_eff"]),
      I0     = list(mean = mean(I0_post),     sd = sd(I0_post),
                    q2.5 = q(I0_post,0.025),  q50 = q(I0_post,0.5),
                    q97.5 = q(I0_post,0.975),
                    rhat = diag_df["I0","Rhat"], ess = diag_df["I0","n_eff"]),
      lambda = list(mean = mean(lambda_post), sd = sd(lambda_post),
                    q2.5 = q(lambda_post,0.025), q50 = q(lambda_post,0.5),
                    q97.5 = q(lambda_post,0.975),
                    rhat = diag_df["lambda","Rhat"], ess = diag_df["lambda","n_eff"])
    ),
    R0  = list(mean = mean(R0_post), sd = sd(R0_post),
               q2.5 = q(R0_post,0.025), q50 = q(R0_post,0.5),
               q97.5 = q(R0_post,0.975)),
    fit = list(
      R2           = R2,
      RMSE         = sqrt(mean(resid_df$residual^2)),
      elpd_loo     = loo_elpd,
      elpd_loo_se  = loo_se,
      n_pareto_k_bad = n_bad_k
    )
  )

  json_path <- file.path(out_dir, sprintf("results_%s.json", site_label))
  write(toJSON(results_list, auto_unbox = TRUE, pretty = TRUE), json_path)
  cat(sprintf("Results saved to outputs/results_%s.json\n", site_label))

  all_results[[site_label]] <- results_list
  cat(sprintf("\n=== Site %s complete ===\n", site_label))
}

# ── Cross-site summary table ───────────────────────────────────────────────────
cat("\n\n══════════════════════════════════════════\n")
cat("  CROSS-SITE SUMMARY\n")
cat("══════════════════════════════════════════\n")
cat(sprintf("%-10s %8s %8s %8s %8s %6s\n",
            "Site", "R0_med", "beta_med", "gamma_med", "rho_med", "R2"))
cat(strrep("-", 58), "\n")
for (nm in names(all_results)) {
  r <- all_results[[nm]]
  cat(sprintf("%-10s %8.3f %8.3f %8.3f %8.4f %6.3f\n",
              nm,
              r$R0$q50,
              r$summary$beta$q50,
              r$summary$gamma$q50,
              r$summary$rho$q50,
              r$fit$R2))
}
cat("\n=== All sites complete ===\n")
