# ============================================================
# Package installer for Bayesian SIR per-site analysis
# Run once: Rscript code/install_packages.R
# ============================================================

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

pkgs <- c("deSolve", "ggplot2", "bayesplot", "gridExtra",
          "jsonlite", "dplyr", "tidyr", "GGally", "officer",
          "readxl", "scales", "patchwork", "loo")

cat("=== Installing CRAN packages ===\n")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    cat(sprintf("Installing %s...\n", p))
    install.packages(p, quiet = TRUE)
  } else {
    cat(sprintf("  %s already installed.\n", p))
  }
}

# Install RStan (from Stan universe)
if (!requireNamespace("rstan", quietly = TRUE)) {
  cat("Installing rstan from stan-dev r-universe...\n")
  install.packages("rstan",
    repos = c("https://stan-dev.r-universe.dev",
              "https://cloud.r-project.org"))
} else {
  cat("  rstan already installed.\n")
}

cat("\n=== Verifying installations ===\n")
all_pkgs <- c(pkgs, "rstan")
ok <- sapply(all_pkgs, requireNamespace, quietly = TRUE)
for (p in all_pkgs) {
  cat(sprintf("  %-15s %s\n", p, if(ok[p]) "OK" else "FAILED"))
}

if (all(ok)) {
  cat("\nAll packages installed successfully!\n")
} else {
  cat("\nWARNING: Some packages failed. Re-run this script.\n")
}
