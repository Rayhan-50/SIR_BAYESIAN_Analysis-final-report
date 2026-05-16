# Bayesian SIR Model Analysis — Raihan Research Group
# Epidemic dataset 2025: weekly case counts

suppressPackageStartupMessages({
  library(rstan)
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
options(mc.cores = 4)

# ── Directories ───────────────────────────────────────────────────────────────
proj_dir  <- "/Users/yusuf/Downloads/Raihan Rsearch"
plots_dir <- file.path(proj_dir, "plots")
out_dir   <- file.path(proj_dir, "outputs")
dir.create(plots_dir, showWarnings = FALSE)
dir.create(out_dir,   showWarnings = FALSE)

save_plot <- function(p, name, w = 7, h = 5) {
  ggsave(file.path(plots_dir, name), plot = p,
         width = w, height = h, dpi = 200, bg = "white")
  invisible(p)
}

bw_theme <- theme_minimal(base_family = "serif") +
  theme(
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    text             = element_text(color = "black"),
    axis.text        = element_text(color = "black"),
    axis.title       = element_text(color = "black"),
    plot.title       = element_text(color = "black", face = "bold"),
    strip.text       = element_text(color = "black")
  )

# ── 1. Load data ──────────────────────────────────────────────────────────────
dat <- read.csv(file.path(proj_dir, "data/clean_epidemic_dataset_2025.csv"),
                stringsAsFactors = FALSE)
dat$date <- as.Date(dat$date)

N      <- 11750L
n_wks  <- nrow(dat)
cases  <- dat$cases

cat("=== Dataset Summary ===\n")
cat(sprintf("Weeks       : %d\n",      n_wks))
cat(sprintf("Total cases : %d\n",      sum(cases)))
cat(sprintf("Peak cases  : %d (week %d)\n", max(cases), which.max(cases) - 1))
cat(sprintf("Population N: %d\n\n",    N))

# ── Plot 01: Raw data ─────────────────────────────────────────────────────────
p01 <- ggplot(dat, aes(x = date, y = cases)) +
  geom_point(shape = 19, size = 2.5, color = "black") +
  labs(title = "Weekly Case Counts — 2025 Epidemic",
       x = "Date", y = "Reported cases") +
  bw_theme
save_plot(p01, "01_raw_data.png")
cat("Plot 01 saved.\n")

# ── 2. SIR ODE helper ────────────────────────────────────────────────────────
sir_ode <- function(t, y, parms) {
  S <- y[1]; I <- y[2]; R <- y[3]
  beta  <- parms["beta"]
  gamma <- parms["gamma"]
  list(c(-beta * S * I / N,
          beta * S * I / N - gamma * I,
          gamma * I))
}

solve_sir <- function(beta, gamma, n = n_wks, pop = N) {
  parms  <- c(beta = beta, gamma = gamma)
  y0     <- c(S = pop - 1, I = 1, R = 0)
  times  <- seq(0, n, by = 1)
  ode(y = y0, times = times, func = sir_ode, parms = parms,
      method = "rk4")
}

# ── Plot 02: SIR compartments with posterior medians ─────────────────────────
beta_med  <- 0.955; gamma_med <- 0.659; rho_med <- 0.080
sol <- solve_sir(beta_med, gamma_med)
sol_df <- as.data.frame(sol) %>%
  filter(time >= 1) %>%
  mutate(week = time, Observed_x_rho = rho_med * I) %>%
  rename(Susceptible = S, Infectious = I, Recovered = R) %>%
  select(week, Susceptible, Infectious, Recovered, Observed_x_rho)

comp_long <- sol_df %>%
  pivot_longer(cols = c(Susceptible, Infectious, Recovered),
               names_to = "Compartment", values_to = "Count")

p02 <- ggplot() +
  geom_line(data = comp_long,
            aes(x = week, y = Count, linetype = Compartment), linewidth = 0.8) +
  geom_line(data = sol_df,
            aes(x = week, y = Observed_x_rho * 10),
            linetype = "dotted", linewidth = 1.1, color = "black") +
  geom_point(data = data.frame(week = 1:n_wks, cases = cases),
             aes(x = week, y = cases * 10),
             shape = 1, size = 2, color = "black") +
  scale_y_continuous(
    name = "Compartment count",
    sec.axis = sec_axis(~ . / 10, name = "Cases / rho*I(t)")
  ) +
  scale_linetype_manual(values = c(Susceptible = "solid",
                                   Infectious   = "dashed",
                                   Recovered    = "dotdash")) +
  labs(title = "SIR Compartments (posterior medians)",
       x = "Week", linetype = "") +
  bw_theme +
  theme(legend.position = "bottom")
save_plot(p02, "02_sir_compartments.png", w = 8, h = 5)
cat("Plot 02 saved.\n")

# ── 3. Prior predictive check ─────────────────────────────────────────────────
set.seed(42)
n_prior <- 200
prior_curves <- lapply(seq_len(n_prior), function(i) {
  R0_s    <- rlnorm(1, log(1.7),  0.25)
  gamma_s <- rlnorm(1, log(0.44), 0.3)
  rho_s   <- rlnorm(1, log(0.05), 0.5)
  phi_s   <- rexp(1, 1)
  beta_s  <- R0_s * gamma_s
  tryCatch({
    sol_s <- solve_sir(beta_s, gamma_s)
    mu_s  <- rho_s * sol_s[-1, "I"]
    mu_s  <- pmax(mu_s, 1e-6)
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
  labs(title = "Prior Predictive Check",
       x = "Week", y = "Expected cases (rho * I(t))") +
  bw_theme
save_plot(p03, "03_prior_predictive.png")
cat("Plot 03 saved.\n")

ok_prior <- any(prior_df$mu > max(cases) * 0.5)
if (!ok_prior) warning("Prior predictive curves may not bracket the data well.")

# ── 4. Compile and fit Stan model ─────────────────────────────────────────────
cat("\n=== Compiling Stan model ===\n")
stan_file <- file.path(proj_dir, "code/sir_model.stan")

fit <- stan(
  file    = stan_file,
  data    = list(N_weeks = n_wks, cases = cases, pop = as.double(N)),
  chains  = 4,
  iter    = 2000,
  warmup  = 1000,
  seed    = 12345,
  control = list(adapt_delta = 0.9, max_treedepth = 12),
  refresh = 200
)

cat("\n=== Convergence Diagnostics ===\n")
sum_fit  <- summary(fit)$summary
params   <- c("R0", "gamma", "rho", "phi", "beta")
diag_df  <- sum_fit[params, c("mean","sd","2.5%","50%","97.5%","n_eff","Rhat")]
print(round(diag_df, 4))

rhat_vals <- diag_df[,"Rhat"]
neff_vals <- diag_df[,"n_eff"]

cat("\n--- R-hat check ---\n")
if (all(rhat_vals < 1.05, na.rm = TRUE)) {
  cat("PASS: All R-hat < 1.05\n")
} else {
  warning("FAIL: Some R-hat >= 1.05 — chains may not have converged!")
  print(rhat_vals[rhat_vals >= 1.05])
}

cat("--- ESS check ---\n")
if (all(neff_vals > 400, na.rm = TRUE)) {
  cat("PASS: All n_eff > 400\n")
} else {
  cat("WARNING: Some n_eff <= 400:\n")
  print(neff_vals[neff_vals <= 400])
}

# ── 5. Extract posterior samples ──────────────────────────────────────────────
post <- rstan::extract(fit)
beta_post  <- post$beta
gamma_post <- post$gamma
rho_post   <- post$rho
phi_post   <- post$phi
R0_post    <- post$R0

# ── Plot 04: Trace plots ──────────────────────────────────────────────────────
posterior_arr <- as.array(fit, pars = c("R0","gamma","rho","phi"))
color_scheme_set("gray")

p04_trace <- mcmc_trace(posterior_arr,
                         pars = c("R0","gamma","rho","phi"),
                         facet_args = list(nrow = 4, labeller = label_parsed)) +
  labs(title = "MCMC Trace Plots") +
  bw_theme
save_plot(p04_trace, "04_trace_plots.png", w = 9, h = 7)
cat("Plot 04 saved.\n")

# ── Plot 05: Posterior histograms ─────────────────────────────────────────────
make_hist <- function(x, lab, xlab = lab) {
  med   <- median(x)
  lo    <- quantile(x, 0.025)
  hi    <- quantile(x, 0.975)
  df    <- data.frame(v = x)
  ggplot(df, aes(x = v)) +
    geom_histogram(fill = "gray70", color = "black", bins = 40) +
    geom_vline(xintercept = med, linetype = "solid",  linewidth = 0.8) +
    geom_vline(xintercept = lo,  linetype = "dashed", linewidth = 0.6) +
    geom_vline(xintercept = hi,  linetype = "dashed", linewidth = 0.6) +
    labs(title = lab, x = xlab, y = "Count") +
    bw_theme
}

h_beta  <- make_hist(beta_post,  expression(beta),  "beta")
h_gamma <- make_hist(gamma_post, expression(gamma), "gamma")
h_rho   <- make_hist(rho_post,   expression(rho),   "rho")
h_phi   <- make_hist(phi_post,   expression(phi),   "phi")
h_R0    <- make_hist(R0_post,    expression(R[0]),  "R0")

p05 <- grid.arrange(h_beta, h_gamma, h_rho, h_phi, h_R0,
                    nrow = 2, ncol = 3,
                    top = "Posterior Distributions (median + 95% CI)")
ggsave(file.path(plots_dir, "05_posterior_histograms.png"),
       plot = p05, width = 12, height = 7, dpi = 200, bg = "white")
cat("Plot 05 saved.\n")

# ── Plot 06: Pair plots ───────────────────────────────────────────────────────
pair_df <- data.frame(beta  = beta_post,
                      gamma = gamma_post,
                      rho   = rho_post,
                      phi   = phi_post)
# Subsample for speed
idx <- sample(nrow(pair_df), min(1000, nrow(pair_df)))
pair_sub <- pair_df[idx, ]

p06 <- ggpairs(pair_sub,
               upper = list(continuous = wrap("cor", size = 3,
                                               color = "black")),
               lower = list(continuous = wrap("points", alpha = 0.15,
                                               size = 0.5, color = "black")),
               diag  = list(continuous = wrap("densityDiag",
                                               fill = "gray70", color = "black"))) +
  labs(title = "Joint Posterior Distributions") +
  bw_theme
ggsave(file.path(plots_dir, "06_pair_plots.png"),
       plot = p06, width = 9, height = 8, dpi = 200, bg = "white")
cat("Plot 06 saved.\n")

# ── 6. Posterior predictive check ─────────────────────────────────────────────
mu_mat <- post$mu   # [draws x weeks]

mu_med <- apply(mu_mat, 2, median)
mu_lo50 <- apply(mu_mat, 2, quantile, 0.25)
mu_hi50 <- apply(mu_mat, 2, quantile, 0.75)
mu_lo90 <- apply(mu_mat, 2, quantile, 0.05)
mu_hi90 <- apply(mu_mat, 2, quantile, 0.95)

pp_df <- data.frame(
  week   = 1:n_wks,
  med    = mu_med,
  lo50   = mu_lo50, hi50 = mu_hi50,
  lo90   = mu_lo90, hi90 = mu_hi90,
  cases  = cases
)

# R²
ss_res <- sum((cases - mu_med)^2)
ss_tot <- sum((cases - mean(cases))^2)
R2 <- 1 - ss_res / ss_tot
cat(sprintf("\nPosterior predictive R² = %.3f\n", R2))
if (R2 < 0.70) warning("R² < 0.70 — model fit may be poor!")

p07 <- ggplot(pp_df, aes(x = week)) +
  geom_ribbon(aes(ymin = lo90, ymax = hi90), fill = "#CCCCCC", alpha = 0.7) +
  geom_ribbon(aes(ymin = lo50, ymax = hi50), fill = "#888888", alpha = 0.7) +
  geom_line(aes(y = med), linewidth = 0.9, color = "black") +
  geom_point(aes(y = cases), shape = 19, size = 2.5, color = "black") +
  annotate("text", x = 5, y = max(cases) * 0.9,
           label = sprintf("R² = %.2f", R2), size = 4, family = "serif") +
  labs(title = "Posterior Predictive Check — Expected Cases",
       subtitle = "Dark gray: 50% CI | Light gray: 90% CI | Line: median",
       x = "Week", y = "Expected cases (rho * I(t))") +
  bw_theme
save_plot(p07, "07_posterior_predictive.png")
cat("Plot 07 saved.\n")

# ── Plot 08: Posterior predictive simulated observations ──────────────────────
yrep_mat <- post$cases_rep
# draw 200 replicated trajectories
n_draws <- min(200, nrow(yrep_mat))
draw_idx <- sample(nrow(yrep_mat), n_draws)
yrep_sub <- yrep_mat[draw_idx, ]

yrep_df <- as.data.frame(t(yrep_sub))
yrep_df$week <- 1:n_wks
yrep_long <- pivot_longer(yrep_df, -week, names_to = "draw", values_to = "y")

p08 <- ggplot() +
  geom_line(data = yrep_long, aes(x = week, y = y, group = draw),
            color = "gray70", alpha = 0.2, linewidth = 0.3) +
  geom_point(data = data.frame(week = 1:n_wks, cases = cases),
             aes(x = week, y = cases), shape = 19, size = 2.5, color = "black") +
  labs(title = "Posterior Predictive — Simulated Observations",
       x = "Week", y = "Cases") +
  bw_theme
save_plot(p08, "08_posterior_predictive_sim.png")
cat("Plot 08 saved.\n")

# ── Plot 09: R0 posterior ─────────────────────────────────────────────────────
R0_df <- data.frame(R0 = R0_post)
p09 <- ggplot(R0_df, aes(x = R0)) +
  geom_histogram(fill = "gray70", color = "black", bins = 50) +
  geom_vline(xintercept = 1.0, linetype = "dashed", linewidth = 1.0) +
  geom_vline(xintercept = median(R0_post), linetype = "solid", linewidth = 0.9) +
  annotate("text", x = 1.02, y = Inf, label = "Epidemic\nthreshold",
           hjust = 0, vjust = 1.5, size = 3.5, family = "serif") +
  annotate("text", x = median(R0_post) + 0.02, y = Inf,
           label = sprintf("Median = %.2f", median(R0_post)),
           hjust = 0, vjust = 3.0, size = 3.5, family = "serif") +
  labs(title = expression("Posterior Distribution of " * R[0]),
       x = expression(R[0]), y = "Count") +
  bw_theme
save_plot(p09, "09_R0_distribution.png")
cat("Plot 09 saved.\n")

# ── Plot 10: Residuals ────────────────────────────────────────────────────────
resid_df <- data.frame(
  week     = 1:n_wks,
  residual = cases - mu_med,
  fitted   = mu_med
)

r10a <- ggplot(resid_df, aes(x = week, y = residual)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(shape = 19, size = 2) +
  geom_line(linewidth = 0.4, color = "gray50") +
  labs(title = "Residuals Over Time", x = "Week", y = "Residual") +
  bw_theme

r10b <- ggplot(resid_df, aes(x = fitted, y = residual)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(shape = 19, size = 2) +
  labs(title = "Residuals vs. Fitted", x = "Fitted (median)", y = "Residual") +
  bw_theme

p10 <- grid.arrange(r10a, r10b, nrow = 1)
ggsave(file.path(plots_dir, "10_residuals.png"),
       plot = p10, width = 10, height = 4.5, dpi = 200, bg = "white")
cat("Plot 10 saved.\n")

# ── 7. Save results summary ───────────────────────────────────────────────────
q <- function(x, p) quantile(x, p, names = FALSE)

results_list <- list(
  N        = N,
  n_weeks  = n_wks,
  total_cases = sum(cases),
  peak_cases  = max(cases),
  peak_week   = which.max(cases) - 1L,
  summary  = list(
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
    phi   = list(mean = mean(phi_post),   sd = sd(phi_post),
                 q2.5 = q(phi_post,0.025), q50 = q(phi_post,0.5),
                 q97.5 = q(phi_post,0.975),
                 rhat = diag_df["phi","Rhat"], ess = diag_df["phi","n_eff"])
  ),
  R0 = list(mean = mean(R0_post), sd = sd(R0_post),
            q2.5 = q(R0_post,0.025), q50 = q(R0_post,0.5),
            q97.5 = q(R0_post,0.975)),
  fit = list(R2 = R2, RMSE = sqrt(mean(resid_df$residual^2)))
)

write(toJSON(results_list, auto_unbox = TRUE, pretty = TRUE),
      file.path(out_dir, "results_R.json"))
cat("\nResults saved to outputs/results_R.json\n")
cat("\n=== Analysis complete ===\n")
