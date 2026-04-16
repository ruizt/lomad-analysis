#!/usr/bin/env Rscript
# sim_job.R ---------------------------------------------------------------
# One simulation replicate for the lomad power / type-I-error study.
#
# Pipeline:
#   make_trends_dist()  →  add_noise()  →  lomad()
#
# All parameters are read from environment variables so the same Docker image
# can cover any cell in a parameter grid (job array pattern).
#
# Environment variables (with defaults):
#   SIM_N            integer  series length             (500)
#   SIM_D            numeric  trend separation distance  (1)
#   SIM_LAMBDA       numeric  target SNR                 (4)
#   SIM_ARMA_P       integer  AR order                   (2)
#   SIM_ARMA_Q       integer  MA order                   (1)
#   SIM_SCALE        numeric  trend scale factor         (5)
#   SIM_B            integer  bootstrap replicates       (500)
#   SIM_METHOD       string   "boot" | "mc" | "analytic" ("boot")
#   SIM_NCORES       integer  parallel cores             (1)
#   SIM_SEED         integer  RNG seed                   (42)
#   SIM_OUT_DIR      string   output directory           (".")
#   SIM_JOB_ID       string   label for output filename  ("job")
# -------------------------------------------------------------------------

env_int <- function(key, default) {
  v <- Sys.getenv(key, unset = NA_character_)
  if (is.na(v) || !nzchar(v)) return(default)
  as.integer(v)
}

env_dbl <- function(key, default) {
  v <- Sys.getenv(key, unset = NA_character_)
  if (is.na(v) || !nzchar(v)) return(default)
  as.numeric(v)
}

env_chr <- function(key, default) {
  v <- Sys.getenv(key, unset = NA_character_)
  if (is.na(v) || !nzchar(v)) return(default)
  v
}

# --- Read parameters -------------------------------------------------------

params <- list(
  n       = env_int("SIM_N",       500L),
  d       = env_dbl("SIM_D",       1),
  lambda  = env_dbl("SIM_LAMBDA",  4),
  arma_p  = env_int("SIM_ARMA_P",  2L),
  arma_q  = env_int("SIM_ARMA_Q",  1L),
  scale   = env_dbl("SIM_SCALE",   5),
  B       = env_int("SIM_B",       500L),
  method  = env_chr("SIM_METHOD",  "boot"),
  ncores  = env_int("SIM_NCORES",  1L),
  seed    = env_int("SIM_SEED",    42L),
  out_dir = env_chr("SIM_OUT_DIR", "."),
  job_id  = env_chr("SIM_JOB_ID",  "job")
)

cat("=== lomad simulation job ===\n")
cat(sprintf(
  "  n=%d  d=%.2f  lambda=%.2f  ARMA(%d,%d)  scale=%.1f\n",
  params$n, params$d, params$lambda, params$arma_p, params$arma_q, params$scale
))
cat(sprintf(
  "  B=%d  method=%s  ncores=%d  seed=%d\n",
  params$B, params$method, params$ncores, params$seed
))

# --- Load package ----------------------------------------------------------

suppressPackageStartupMessages(library(lomad))

# --- Step 1: Generate trend pair ------------------------------------------

cat("\n[1/3] Generating trends ...\n")
trends <- make_trends_dist(
  n    = params$n,
  d    = params$d,
  seed = params$seed
)

# --- Step 2: Add calibrated ARMA noise ------------------------------------

cat("[2/3] Adding ARMA noise ...\n")
sim <- add_noise(
  trends        = trends,
  h             = 30L,                          # MA window for SNR calibration
  lambda_target = params$lambda,
  scale         = params$scale,
  order         = c(params$arma_p, params$arma_q),
  seed          = params$seed + 1L
)

cat(sprintf(
  "      series1 mean_snr=%.3f  series2 mean_snr=%.3f\n",
  sim$noise$series1$mean_snr,
  sim$noise$series2$mean_snr
))

# --- Step 3: Fit and test --------------------------------------------------

cat(sprintf("[3/3] Running lomad (method=%s, B=%d) ...\n", params$method, params$B))
t0  <- proc.time()

result <- lomad(
  x1     = sim$y1,
  x2     = sim$y2,
  method = params$method,
  B      = params$B,
  seed   = params$seed + 2L,
  ncores = params$ncores,
  verbose = FALSE
)

elapsed <- (proc.time() - t0)[["elapsed"]]
cat(sprintf("      done in %.1f s\n", elapsed))

# --- Collect output -------------------------------------------------------

output <- list(
  params   = params,
  p_values = result$p_values,
  observed = result$observed,
  expected = result$expected,
  elapsed  = elapsed,
  noise_diagnostics = list(
    series1 = sim$noise$series1,
    series2 = sim$noise$series2
  )
)

cat("\n--- p-values ---\n")
for (nm in names(result$p_values)) {
  cat(sprintf("  %-20s %.4f\n", nm, result$p_values[[nm]]))
}

# --- Save results ---------------------------------------------------------

dir.create(params$out_dir, showWarnings = FALSE, recursive = TRUE)
out_file <- file.path(
  params$out_dir,
  sprintf("result_%s_n%d_d%.1f_lam%.1f_arma%d%d_seed%d.rds",
          params$job_id, params$n, params$d, params$lambda,
          params$arma_p, params$arma_q, params$seed)
)

saveRDS(output, out_file)
cat(sprintf("\nSaved → %s\n", out_file))
