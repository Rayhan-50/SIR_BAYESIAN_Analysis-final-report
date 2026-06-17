# =============================================================================
# LOO Cross-Validation + Full Prior Sensitivity Analysis
# Requires analysis.R to have been run first (Stan model compiled).
# Run: Rscript code/loo_sensitivity.R
# Outputs:
#   plots/11_loo_pareto_k.png
#   plots/12_sensitivity_R0.png
#   plots/13_sensitivity_rho.png
#   plots/14_ppc_coverage.png
#   outputs/loo_sensitivity_results.json
# =============================================================================

suppressPackageStartupMessages({
  library(rstan)
  library(loo)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(gridExtra)
  library(jsonlite)
})

rstan_options(auto_write = TRUE)
options(mc.cores = min(parallel::detectCores(), 4))

proj_dir <- tryCatch(
  normalizePath(file.path(dirname(sys.frame(1)$ofile), ".."), mustWork = TRUE),
  error = function(e) getwd()
)
plots_dir <- file.path(proj_dir, "plots")
out_dir   <- file.path(proj_dir, "outputs")
dir.create(plots_dir, showWarnings = FALSE)
dir.create(out_dir,   showWarnings = FALSE)

bw_theme <- theme_minimal(base_family = "serif") +
  theme(
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    text       = element_text(color = "black"),
    axis.text  = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    plot.title = element_text(color = "black", face = "bold"),
    strip.text = element_text(color = "black")
  )

# ── Load data ──────────────────────────────────────────────────────────────────
dat   <- read.csv(file.path(proj_dir, "data/clean_epidemic_dataset_2025.csv"),
                  stringsAsFactors = FALSE)
n_wks <- nrow(dat)
cases <- dat$cases
N     <- 11750L

stan_file  <- file.path(proj_dir, "code", "sir_model.stan")
stan_model_obj <- stan_model(file = stan_file)
cat("Stan model compiled.\n")

# ── Helper: fit one scenario ───────────────────────────────────────────────────
fit_scenario <- function(priors, seed = 12345) {
  stan_data <- c(list(N_weeks = n_wks, cases = cases, pop = as.double(N)), priors)
  sampling(
    stan_model_obj,
    data    = stan_data,
    chains  = 4,
    iter    = 3000,
    warmup  = 1500,
    seed    = seed,
    control = list(adapt_delta = 0.98, max_treedepth = 15),
    refresh = 0
  )
}

# ── Sensitivity scenarios ──────────────────────────────────────────────────────
# Each scenario changes ONE prior at a time relative to baseline.
# Baseline matches sir_model.stan documentation exactly.
scenarios <- list(

  baseline = list(
    label       = "Baseline\nR0~LN(log2, 0.3), rho~Beta(2,10)",
    R0_mu       = log(2.0), R0_sigma     = 0.3,
    gamma_mu    = log(0.5), gamma_sigma  = 0.3,
    rho_alpha   = 2,        rho_beta     = 10,
    phi_rate    = 0.5,
    I0_mu       = log(1),   I0_sigma     = 2.0,
    lambda_rate = 0.5
  ),

  R0_narrow = list(
    label       = "R0 tighter\nR0~LN(log2, 0.15)",
    R0_mu       = log(2.0), R0_sigma     = 0.15,
    gamma_mu    = log(0.5), gamma_sigma  = 0.3,
    rho_alpha   = 2,        rho_beta     = 10,
    phi_rate    = 0.5,
    I0_mu       = log(1),   I0_sigma     = 2.0,
    lambda_rate = 0.5
  ),

  R0_wide = list(
    label       = "R0 wider\nR0~LN(log2, 0.5)",
    R0_mu       = log(2.0), R0_sigma     = 0.5,
    gamma_mu    = log(0.5), gamma_sigma  = 0.3,
    rho_alpha   = 2,        rho_beta     = 10,
    phi_rate    = 0.5,
    I0_mu       = log(1),   I0_sigma     = 2.0,
    lambda_rate = 0.5
  ),

  rho_diffuse = list(
    label       = "rho diffuse\nrho~Beta(1,5)",
    R0_mu       = log(2.0), R0_sigma     = 0.3,
    gamma_mu    = log(0.5), gamma_sigma  = 0.3,
    rho_alpha   = 1,        rho_beta     = 5,
    phi_rate    = 0.5,
    I0_mu       = log(1),   I0_sigma     = 2.0,
    lambda_rate = 0.5
  ),

  rho_tight = list(
    label       = "rho tight\nrho~Beta(3,20)",
    R0_mu       = log(2.0), R0_sigma     = 0.3,
    gamma_mu    = log(0.5), gamma_sigma  = 0.3,
    rho_alpha   = 3,        rho_beta     = 20,
    phi_rate    = 0.5,
    I0_mu       = log(1),   I0_sigma     = 2.0,
    lambda_rate = 0.5
  )
)

# ── Run all scenarios ──────────────────────────────────────────────────────────
cat("\n=== Running", length(scenarios), "sensitivity scenarios ===\n")
results <- list()

for (nm in names(scenarios)) {
  s <- scenarios[[nm]]
  cat(sprintf("\n--- Scenario: %s ---\n", nm))

  fit_s <- tryCatch(fit_scenario(s[setdiff(names(s), "label")]),
                    error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL })
  if (is.null(fit_s)) next

  post_s      <- rstan::extract(fit_s)
  R0_post_s   <- post_s$R0
  rho_post_s  <- post_s$rho
  gamma_post_s <- post_s$gamma

  # LOO
  ll_mat <- extract_log_lik(fit_s, parameter_name = "log_lik", merge_chains = FALSE)
  loo_s  <- loo(ll_mat, r_eff = relative_eff(exp(ll_mat)))

  # R-hat
  sum_s   <- summary(fit_s)$summary
  rhat_ok <- all(sum_s[c("R0","gamma","rho","phi"), "Rhat"] < 1.01, na.rm = TRUE)

  results[[nm]] <- list(
    label    = s$label,
    R0_med   = median(R0_post_s),
    R0_lo    = quantile(R0_post_s, 0.025),
    R0_hi    = quantile(R0_post_s, 0.975),
    rho_med  = median(rho_post_s),
    rho_lo   = quantile(rho_post_s, 0.025),
    rho_hi   = quantile(rho_post_s, 0.975),
    gamma_med = median(gamma_post_s),
    elpd_loo  = loo_s$estimates["elpd_loo", "Estimate"],
    elpd_se   = loo_s$estimates["elpd_loo", "SE"],
    n_bad_k   = sum(loo_s$pointwise[, "influence_pareto_k"] > 0.7),
    rhat_ok   = rhat_ok
  )

  cat(sprintf("  R0 = %.3f [%.3f, %.3f]  rho = %.3f  ELPD = %.1f (SE %.1f)  Rhat OK: %s\n",
              results[[nm]]$R0_med, results[[nm]]$R0_lo, results[[nm]]$R0_hi,
              results[[nm]]$rho_med,
              results[[nm]]$elpd_loo, results[[nm]]$elpd_se,
              results[[nm]]$rhat_ok))

  # Save Pareto k plot for baseline only
  if (nm == "baseline") {
    k_df <- data.frame(
      week = seq_along(loo_s$pointwise[, "influence_pareto_k"]),
      k    = loo_s$pointwise[, "influence_pareto_k"]
    )
    p_k <- ggplot(k_df, aes(x = week, y = k)) +
      geom_point(aes(shape = ifelse(k > 0.7, "triangle", "circle")),
                 size = 2.5, fill = "gray50", color = "black") +
      geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50") +
      geom_hline(yintercept = 0.7, linetype = "solid",  color = "black") +
      scale_shape_manual(values = c(circle = 21, triangle = 24),
                         labels = c("k < 0.7", "k > 0.7 (unreliable)"),
                         name = "Quality") +
      annotate("text", x = max(k_df$week) * 0.95, y = 0.73,
               label = "Unreliable threshold (k = 0.7)",
               hjust = 1, size = 3.2, family = "serif") +
      labs(title = "Pareto-k Diagnostics — LOO-CV (Baseline Model)",
           subtitle = "k > 0.7 flags highly influential observations requiring moment-matching",
           x = "Week index", y = expression(hat(k))) +
      bw_theme
    ggsave(file.path(plots_dir, "11_loo_pareto_k.png"),
           plot = p_k, width = 9, height = 4.5, dpi = 200, bg = "white")
    cat("  Plot 11 saved.\n")
  }
}

# ── Plot 12: R0 forest plot (sensitivity) ─────────────────────────────────────
sens_df <- bind_rows(lapply(names(results), function(nm) {
  r <- results[[nm]]
  data.frame(scenario = r$label, R0_med = r$R0_med, R0_lo = r$R0_lo, R0_hi = r$R0_hi,
             is_baseline = nm == "baseline", stringsAsFactors = FALSE)
})) %>% mutate(scenario = factor(scenario, levels = rev(scenario)))

p12 <- ggplot(sens_df, aes(x = R0_med, y = scenario)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  geom_errorbarh(aes(xmin = R0_lo, xmax = R0_hi, linetype = is_baseline),
                 height = 0.25, linewidth = 0.8) +
  geom_point(aes(shape = is_baseline), size = 3.5, fill = "gray30", color = "black") +
  scale_shape_manual(values = c(`TRUE` = 22, `FALSE` = 21),
                     labels = c(`TRUE` = "Baseline", `FALSE` = "Sensitivity"),
                     name = "") +
  scale_linetype_manual(values = c(`TRUE` = "solid", `FALSE` = "dashed"),
                        guide = "none") +
  labs(title = expression("Prior Sensitivity Analysis — Posterior " * R[0]),
       subtitle = "Squares = baseline; circles = alternative priors. Error bars = 95% CrI.",
       x = expression(R[0] ~ "(posterior median + 95% credible interval)"),
       y = "Prior scenario") +
  bw_theme + theme(legend.position = "bottom")
ggsave(file.path(plots_dir, "12_sensitivity_R0.png"),
       plot = p12, width = 9, height = 5, dpi = 200, bg = "white")
cat("Plot 12 saved.\n")

# ── Plot 13: rho forest plot ───────────────────────────────────────────────────
sens_df2 <- bind_rows(lapply(names(results), function(nm) {
  r <- results[[nm]]
  data.frame(scenario = r$label, rho_med = r$rho_med,
             rho_lo = r$rho_lo, rho_hi = r$rho_hi,
             is_baseline = nm == "baseline", stringsAsFactors = FALSE)
})) %>% mutate(scenario = factor(scenario, levels = rev(scenario)))

p13 <- ggplot(sens_df2, aes(x = rho_med, y = scenario)) +
  geom_errorbarh(aes(xmin = rho_lo, xmax = rho_hi, linetype = is_baseline),
                 height = 0.25, linewidth = 0.8) +
  geom_point(aes(shape = is_baseline), size = 3.5, fill = "gray30", color = "black") +
  scale_shape_manual(values = c(`TRUE` = 22, `FALSE` = 21), guide = "none") +
  scale_linetype_manual(values = c(`TRUE` = "solid", `FALSE` = "dashed"), guide = "none") +
  labs(title = expression("Prior Sensitivity Analysis — Reporting Rate " * rho),
       subtitle = "Squares = baseline; circles = alternative priors. Error bars = 95% CrI.",
       x = expression(rho ~ "(case reporting fraction)"),
       y = "Prior scenario") +
  bw_theme
ggsave(file.path(plots_dir, "13_sensitivity_rho.png"),
       plot = p13, width = 9, height = 5, dpi = 200, bg = "white")
cat("Plot 13 saved.\n")

# ── Plot 14: ELPD comparison ───────────────────────────────────────────────────
elpd_df <- bind_rows(lapply(names(results), function(nm) {
  r <- results[[nm]]
  data.frame(scenario = r$label, elpd = r$elpd_loo, se = r$elpd_se,
             stringsAsFactors = FALSE)
})) %>% mutate(
  scenario = factor(scenario, levels = rev(scenario)),
  lo = elpd - 2 * se,
  hi = elpd + 2 * se
)

p14 <- ggplot(elpd_df, aes(x = elpd, y = scenario)) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.25, linewidth = 0.7) +
  geom_point(size = 3.5, shape = 21, fill = "gray40", color = "black") +
  labs(title = "LOO-ELPD by Prior Scenario",
       subtitle = "Error bars = ±2 SE. Higher ELPD = better out-of-sample predictive accuracy.",
       x = "ELPD (LOO cross-validation)", y = "Prior scenario") +
  bw_theme
ggsave(file.path(plots_dir, "14_elpd_sensitivity.png"),
       plot = p14, width = 9, height = 5, dpi = 200, bg = "white")
cat("Plot 14 saved.\n")

# ── Posterior predictive coverage (baseline) ───────────────────────────────────
# Re-fit baseline if needed (reuses fit from first scenario)
cat("\n=== Posterior Predictive Coverage Check ===\n")
baseline_fit <- tryCatch(fit_scenario(scenarios$baseline[setdiff(names(scenarios$baseline), "label")]),
                         error = function(e) NULL)
if (!is.null(baseline_fit)) {
  post_b   <- rstan::extract(baseline_fit)
  yrep_mat <- post_b$cases_rep

  cov90 <- mean(sapply(seq_along(cases), function(t) {
    lo <- quantile(yrep_mat[, t], 0.05); hi <- quantile(yrep_mat[, t], 0.95)
    cases[t] >= lo & cases[t] <= hi
  }))
  cov50 <- mean(sapply(seq_along(cases), function(t) {
    lo <- quantile(yrep_mat[, t], 0.25); hi <- quantile(yrep_mat[, t], 0.75)
    cases[t] >= lo & cases[t] <= hi
  }))
  ppp_mean <- mean(apply(yrep_mat, 1, mean) >= mean(cases))
  ppp_max  <- mean(apply(yrep_mat, 1, max)  >= max(cases))

  cat(sprintf("  90%% PI coverage: %.1f%%  (nominal: 90%%)\n", cov90 * 100))
  cat(sprintf("  50%% PI coverage: %.1f%%  (nominal: 50%%)\n", cov50 * 100))
  cat(sprintf("  Bayesian p-value (mean): %.3f\n", ppp_mean))
  cat(sprintf("  Bayesian p-value (max):  %.3f\n", ppp_max))

  # Plot 15: coverage visualisation
  mu_mat  <- post_b$mu
  cov_df  <- data.frame(
    week  = seq_along(cases),
    obs   = cases,
    med   = apply(mu_mat, 2, median),
    lo50  = apply(mu_mat, 2, quantile, 0.25),
    hi50  = apply(mu_mat, 2, quantile, 0.75),
    lo90  = apply(mu_mat, 2, quantile, 0.05),
    hi90  = apply(mu_mat, 2, quantile, 0.95)
  )
  p15 <- ggplot(cov_df, aes(x = week)) +
    geom_ribbon(aes(ymin = lo90, ymax = hi90), fill = "#CCCCCC", alpha = 0.7) +
    geom_ribbon(aes(ymin = lo50, ymax = hi50), fill = "#888888", alpha = 0.7) +
    geom_line(aes(y = med), linewidth = 0.9) +
    geom_point(aes(y = obs), shape = 19, size = 2.5) +
    annotate("text", x = 2, y = max(cases) * 0.95,
             label = sprintf("90%% PI: %.0f%% coverage\n50%% PI: %.0f%% coverage\nBayes p (mean)=%.2f",
                             cov90 * 100, cov50 * 100, ppp_mean),
             hjust = 0, vjust = 1, size = 3.5, family = "serif") +
    labs(title = "Posterior Predictive Coverage — Baseline Model",
         subtitle = "Dark: 50% PI  |  Light: 90% PI  |  Dots: observed",
         x = "Week", y = "Weekly cases") + bw_theme
  ggsave(file.path(plots_dir, "15_ppc_coverage.png"),
         plot = p15, width = 9, height = 5, dpi = 200, bg = "white")
  cat("Plot 15 saved.\n")
} else {
  cov90 <- cov50 <- ppp_mean <- ppp_max <- NA
}

# ── Summary table ──────────────────────────────────────────────────────────────
cat("\n================================================================\n")
cat("  SENSITIVITY SUMMARY\n")
cat("================================================================\n")
cat(sprintf("%-30s %7s %12s %8s %7s %6s\n",
            "Scenario", "R0_med", "R0_95CrI", "rho_med", "ELPD", "BadK"))
cat(strrep("-", 78), "\n")
for (nm in names(results)) {
  r <- results[[nm]]
  cat(sprintf("%-30s %7.3f [%5.3f,%5.3f] %8.4f %7.1f %6d\n",
              gsub("\n", " ", r$label),
              r$R0_med, r$R0_lo, r$R0_hi,
              r$rho_med, r$elpd_loo, r$n_bad_k))
}

# ── Save all results ───────────────────────────────────────────────────────────
out <- list(
  scenarios       = lapply(results, function(r) r[names(r) != "label"]),
  ppc_coverage_90 = cov90,
  ppc_coverage_50 = cov50,
  ppp_mean        = ppp_mean,
  ppp_max         = ppp_max
)
write(toJSON(out, auto_unbox = TRUE, pretty = TRUE),
      file.path(out_dir, "loo_sensitivity_results.json"))
cat("\nResults saved to outputs/loo_sensitivity_results.json\n")
cat("=== LOO and sensitivity analysis complete ===\n")
