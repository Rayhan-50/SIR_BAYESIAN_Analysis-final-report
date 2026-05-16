# =============================================================================
# Bayesian SIR Model — Publication-Quality Fitted Plot
# Shows: Observed data vs. Posterior Predictive (fitted model)
# Run this AFTER running analysis.R (which saves the Stan fit and posteriors)
# =============================================================================

suppressPackageStartupMessages({
  library(rstan)
  library(deSolve)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(jsonlite)
  library(gridExtra)
  library(grid)
})

# ── Paths (update proj_dir to your local path) ────────────────────────────────
proj_dir  <- "d:/research/research/SIR_Bayesian_Analysis_Report"
plots_dir <- file.path(proj_dir, "plots")
out_dir   <- file.path(proj_dir, "outputs")
dir.create(plots_dir, showWarnings = FALSE)

# ── Load data ─────────────────────────────────────────────────────────────────
dat   <- read.csv(file.path(proj_dir, "data/clean_epidemic_dataset_2025.csv"),
                  stringsAsFactors = FALSE)
dat$date <- as.Date(dat$date)
N     <- 11750L
cases <- dat$cases
n_wks <- nrow(dat)

# ── Load posterior results (from results_R.json) ──────────────────────────────
res   <- fromJSON(file.path(out_dir, "results_R.json"))

# Posterior median parameters
beta_med  <- res$summary$beta$q50
gamma_med <- res$summary$gamma$q50
rho_med   <- res$summary$rho$q50
R2_val    <- res$fit$R2
RMSE_val  <- res$fit$RMSE
R0_med    <- res$R0$q50
R0_lo     <- res$R0$q2.5
R0_hi     <- res$R0$q97.5

cat(sprintf("Using: beta=%.4f  gamma=%.4f  rho=%.4f\n",
            beta_med, gamma_med, rho_med))
cat(sprintf("R0 = %.3f [%.3f, %.3f]\n", R0_med, R0_lo, R0_hi))

# ── SIR ODE solver ────────────────────────────────────────────────────────────
sir_ode <- function(t, y, parms) {
  S <- y[1]; I <- y[2]
  beta  <- parms["beta"]
  gamma <- parms["gamma"]
  list(c(-beta * S * I / N,
          beta * S * I / N - gamma * I,
          gamma * I))
}

solve_sir <- function(beta, gamma, pop = N, n = n_wks) {
  y0    <- c(S = pop - 1, I = 1, R = 0)
  times <- seq(0, n, by = 1)
  as.data.frame(ode(y = y0, times = times, func = sir_ode,
                    parms = c(beta = beta, gamma = gamma), method = "rk4"))
}

# ── Compute fitted trajectory from posterior medians ──────────────────────────
sol    <- solve_sir(beta_med, gamma_med)
sol_wk <- sol %>% filter(time >= 1) %>%
  mutate(week = as.integer(time),
         date = dat$date,
         fitted = rho_med * I)

# ── Re-run Stan to extract posterior predictive intervals ─────────────────────
# (If the Stan fit object is not in memory, we reconstruct CI from the median
#  solution and known overdispersion — a fast approximation)

# We simulate uncertainty bands analytically using negative-binomial quantiles
phi_med <- res$summary$phi$q50
phi_lo  <- res$summary$phi$q2.5
phi_hi  <- res$summary$phi$q97.5

mu_fit  <- sol_wk$fitted  # posterior median expected cases

# 50% and 90% CI from the Negative-Binomial predictive distribution
# Using nbinom parameterisation: mean = mu, size = phi
# lower.tail CI
ci_lo90 <- qnbinom(0.05,  size = phi_med, mu = pmax(mu_fit, 1e-6))
ci_hi90 <- qnbinom(0.95,  size = phi_med, mu = pmax(mu_fit, 1e-6))
ci_lo50 <- qnbinom(0.25,  size = phi_med, mu = pmax(mu_fit, 1e-6))
ci_hi50 <- qnbinom(0.75,  size = phi_med, mu = pmax(mu_fit, 1e-6))

plot_df <- data.frame(
  week    = sol_wk$week,
  date    = sol_wk$date,
  fitted  = mu_fit,
  lo90    = ci_lo90,
  hi90    = ci_hi90,
  lo50    = ci_lo50,
  hi50    = ci_hi50,
  observed = cases
)

# ── Annotation labels ─────────────────────────────────────────────────────────
ann_text <- sprintf(
  "R\u2080 = %.2f  [95%% CrI: %.2f \u2013 %.2f]\nR\u00b2 = %.3f  |  RMSE = %.1f cases/wk\n\u03b2 = %.3f  |  \u03b3 = %.3f  |  \u03c1 = %.3f",
  R0_med, R0_lo, R0_hi, R2_val, RMSE_val,
  beta_med, gamma_med, rho_med
)

# ── PLOT A: Date on X-axis (for presentation) ─────────────────────────────────
pA <- ggplot(plot_df, aes(x = date)) +
  # 90% predictive band
  geom_ribbon(aes(ymin = lo90, ymax = hi90),
              fill = "#4E79A7", alpha = 0.18) +
  # 50% predictive band
  geom_ribbon(aes(ymin = lo50, ymax = hi50),
              fill = "#4E79A7", alpha = 0.35) +
  # Fitted median curve
  geom_line(aes(y = fitted), color = "#2C5F8A",
            linewidth = 1.1, linetype = "solid") +
  # Observed data points
  geom_point(aes(y = observed), shape = 21, size = 3,
             color = "black", fill = "#E15759", stroke = 0.8) +
  # Annotation box
  annotate("label", x = as.Date("2025-01-20"), y = max(cases) * 0.90,
           label = ann_text, hjust = 0, vjust = 1,
           size = 3.5, family = "serif", color = "black",
           fill = "white", label.size = 0.3,
           label.padding = unit(0.4, "lines")) +
  # Vertical line at peak
  geom_vline(xintercept = dat$date[which.max(cases)],
             linetype = "dotted", color = "gray40", linewidth = 0.6) +
  annotate("text",
           x = dat$date[which.max(cases)] + 4,
           y = max(cases),
           label = sprintf("Peak: %d cases\n(%s)", max(cases),
                           format(dat$date[which.max(cases)], "%b %d")),
           hjust = 0, size = 3.2, family = "serif", color = "gray30") +
  scale_x_date(date_labels = "%b %Y", date_breaks = "2 months",
               expand = expansion(mult = 0.02)) +
  scale_y_continuous(breaks = pretty_breaks(6),
                     limits = c(0, max(cases) * 1.12),
                     expand = expansion(mult = 0.01)) +
  labs(
    title    = "Bayesian SIR Model — Fitted vs. Observed Cases",
    subtitle = "Blue line: posterior median fit  |  Shaded: 50% and 90% predictive intervals  |  Red dots: observed weekly cases",
    x        = "Date (2025)",
    y        = "Weekly reported cases",
    caption  = "Model: SIR-ODE fitted via Stan MCMC (4 chains × 1000 post-warmup samples)\nLikelihood: Negative-Binomial | Population N = 11,750"
  ) +
  theme_classic(base_family = "serif", base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 15, color = "black"),
    plot.subtitle = element_text(size = 10, color = "gray30"),
    plot.caption  = element_text(size = 8.5, color = "gray40", hjust = 0),
    axis.title    = element_text(face = "bold"),
    axis.text.x   = element_text(angle = 30, hjust = 1),
    panel.grid.major.y = element_line(color = "gray92", linewidth = 0.4),
    panel.grid.major.x = element_blank(),
    panel.border  = element_rect(color = "gray70", fill = NA),
    plot.margin   = margin(12, 16, 10, 10)
  )

ggsave(file.path(plots_dir, "fitted_model_presentation.png"),
       plot = pA, width = 10, height = 6, dpi = 300, bg = "white")
cat("Saved: fitted_model_presentation.png\n")

# ── PLOT B: Week on X-axis + residuals panel (for professor/report) ───────────
pB_top <- ggplot(plot_df, aes(x = week)) +
  geom_ribbon(aes(ymin = lo90, ymax = hi90),
              fill = "#4E79A7", alpha = 0.18) +
  geom_ribbon(aes(ymin = lo50, ymax = hi50),
              fill = "#4E79A7", alpha = 0.35) +
  geom_line(aes(y = fitted), color = "#2C5F8A",
            linewidth = 1.1) +
  geom_point(aes(y = observed), shape = 21, size = 2.8,
             color = "black", fill = "#E15759", stroke = 0.7) +
  annotate("label",
           x = 1, y = max(cases) * 0.95,
           label = sprintf("R\u2080 = %.2f [%.2f, %.2f]\nR\u00b2 = %.3f | RMSE = %.1f",
                           R0_med, R0_lo, R0_hi, R2_val, RMSE_val),
           hjust = 0, vjust = 1, size = 3.5, family = "serif",
           fill = "white", label.size = 0.3,
           label.padding = unit(0.35, "lines")) +
  scale_y_continuous(limits = c(0, max(cases) * 1.12),
                     breaks = pretty_breaks(5)) +
  labs(
    title = "Bayesian SIR Model — Fit to Observed Data",
    x = NULL, y = "Weekly cases"
  ) +
  theme_classic(base_family = "serif", base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    panel.grid.major.y = element_line(color = "gray92", linewidth = 0.4),
    panel.border = element_rect(color = "gray70", fill = NA)
  )

# Residuals sub-panel
plot_df$residual <- plot_df$observed - plot_df$fitted
pB_bot <- ggplot(plot_df, aes(x = week, y = residual)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50",
             linewidth = 0.7) +
  geom_segment(aes(xend = week, yend = 0),
               color = "gray70", linewidth = 0.4) +
  geom_point(shape = 21, size = 2.5, fill = ifelse(plot_df$residual >= 0,
             "#E15759", "#4E79A7"), color = "black", stroke = 0.6) +
  scale_y_continuous(breaks = pretty_breaks(4)) +
  labs(
    x = "Week (2025)",
    y = "Residual\n(Obs − Fitted)",
    caption = "Red dots: positive residual (underprediction) | Blue dots: negative residual (overprediction)"
  ) +
  theme_classic(base_family = "serif", base_size = 11) +
  theme(
    axis.title = element_text(face = "bold"),
    plot.caption = element_text(size = 8, color = "gray40", hjust = 0),
    panel.grid.major.y = element_line(color = "gray92", linewidth = 0.4),
    panel.border = element_rect(color = "gray70", fill = NA)
  )

pB_combined <- grid.arrange(pB_top, pB_bot,
                             nrow = 2, heights = c(2.8, 1))
ggsave(file.path(plots_dir, "fitted_model_with_residuals.png"),
       plot = pB_combined, width = 10, height = 7.5, dpi = 300, bg = "white")
cat("Saved: fitted_model_with_residuals.png\n")

cat("\n=== Done! Show your professor: ===\n")
cat("  1. fitted_model_presentation.png  (clean single plot)\n")
cat("  2. fitted_model_with_residuals.png (fit + residuals, for reports)\n")
