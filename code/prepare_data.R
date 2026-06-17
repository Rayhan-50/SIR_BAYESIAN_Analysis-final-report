# =============================================================================
# Data Preparation & Provenance Script
# Bayesian SIR Model Analysis — Raihan Research Group
# =============================================================================
# PURPOSE:
#   This script documents the origin of all datasets, validates their
#   structure, and produces a clear provenance log saved to outputs/.
#
# DATASET CLARIFICATION (Important for reproducibility):
#   - clean_epidemic_dataset_2025.csv  : PRIMARY community surveillance stream
#     (N = 11,750 catchment). This is an INDEPENDENT reporting system from
#     the sentinel sites. It is NOT the sum of the three sentinel sites.
#
#   - site_ILI.csv    : Outpatient ILI sentinel (same N=11,750 catchment)
#   - site_Severe.csv : Inpatient hospital sentinel (same N=11,750 catchment)
#   - site_SARI.csv   : SARI specialist sentinel (N=40,000 hospital ward)
#
#   NOTE: ILI and Severe cover overlapping catchment populations (N=11,750
#   each), so summing their cases would represent DOUBLE-COUNTING. Each site
#   represents a distinct clinical reporting stream, not a sub-population.
#
# Run: Rscript code/prepare_data.R
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(jsonlite)
})

proj_dir <- tryCatch(
  normalizePath(file.path(dirname(sys.frame(1)$ofile), ".."), mustWork = TRUE),
  error = function(e) getwd()
)
out_dir  <- file.path(proj_dir, "outputs")
dir.create(out_dir, showWarnings = FALSE)

cat("=============================================================\n")
cat(" DATA PREPARATION & VALIDATION REPORT\n")
cat(" Bayesian SIR Model — Raihan Research Group\n")
cat("=============================================================\n\n")

# ── Load all datasets ──────────────────────────────────────────────────────────
files <- list(
  community = list(
    path = file.path(proj_dir, "data", "clean_epidemic_dataset_2025.csv"),
    N    = 11750L,
    label = "Community Surveillance (Primary)",
    type  = "epidemic"
  ),
  ILI = list(
    path = file.path(proj_dir, "data", "site_ILI.csv"),
    N    = 11750L,
    label = "ILI Outpatient Sentinel",
    type  = "epidemic"
  ),
  Severe = list(
    path = file.path(proj_dir, "data", "site_Severe.csv"),
    N    = 11750L,
    label = "Severe/Hospital Inpatient Sentinel",
    type  = "epidemic"
  )
)

provenance <- list()

for (nm in names(files)) {
  spec <- files[[nm]]
  if (!file.exists(spec$path)) {
    cat(sprintf("[MISSING] %s — file not found: %s\n", nm, spec$path))
    next
  }

  dat <- read.csv(spec$path, stringsAsFactors = FALSE)
  dat$date <- as.Date(dat$date)

  n_weeks   <- nrow(dat)
  total     <- sum(dat$cases)
  peak_val  <- max(dat$cases)
  peak_wk   <- which.max(dat$cases)
  week1_val <- dat$cases[1]
  last_val  <- dat$cases[n_weeks]

  # ── SIR assumption check: epidemic should start near 0 ────────────────────
  sir_ok <- week1_val <= 5
  stair  <- any(rle(dat$cases)$lengths >= 3 & rle(dat$cases)$values > 0)

  cat(sprintf("── %s (%s) ──\n", nm, spec$label))
  cat(sprintf("   File       : %s\n", basename(spec$path)))
  cat(sprintf("   N          : %s\n", formatC(spec$N, format="d", big.mark=",")))
  cat(sprintf("   Weeks      : %d\n", n_weeks))
  cat(sprintf("   Total cases: %d\n", total))
  cat(sprintf("   Peak       : %d (week %d)\n", peak_val, peak_wk))
  cat(sprintf("   Week 1     : %d cases\n", week1_val))
  cat(sprintf("   SIR start OK (week1<=5) : %s\n", ifelse(sir_ok, "YES", "NO — possible endemic baseline")))
  cat(sprintf("   Stair-step artifact    : %s\n", ifelse(stair, "YES — check monthly reporting", "no")))
  cat(sprintf("   Model type : %s\n\n", spec$type))

  provenance[[nm]] <- list(
    label       = spec$label,
    file        = basename(spec$path),
    N           = spec$N,
    n_weeks     = n_weeks,
    total_cases = total,
    peak_cases  = peak_val,
    peak_week   = peak_wk,
    week1_cases = week1_val,
    sir_start_valid = sir_ok,
    stair_artifact  = stair,
    model_type  = spec$type
  )
}

# ── Dataset relationship note ──────────────────────────────────────────────────
cat("=============================================================\n")
cat(" DATASET RELATIONSHIPS\n")
cat("=============================================================\n")
cat(" Community surveillance (clean_epidemic_dataset_2025.csv) is an\n")
cat(" INDEPENDENT reporting stream from the three sentinel sites.\n")
cat(" ILI and Severe share the same catchment (N=11,750), so their\n")
cat(" case counts CANNOT be summed without double-counting.\n\n")
cat(" Correct analytical approach:\n")
cat("   [1] Community surveillance: independent primary Bayesian SIR fit\n")
cat("   [2] ILI & Severe: separate per-site Bayesian SIR fits\n\n")

# ── Save provenance JSON ───────────────────────────────────────────────────────
prov_path <- file.path(out_dir, "data_provenance.json")
write(toJSON(provenance, auto_unbox = TRUE, pretty = TRUE), prov_path)
cat(sprintf("Provenance saved: %s\n", prov_path))
cat("=== Data preparation complete ===\n")
