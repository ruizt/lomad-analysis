## Power study — Kubernetes container entrypoint
##
## One container = all S replicates for one (structure, d, n, snr) combination.
## Parameters are passed as environment variables by the job spec.
##
## Environment variables:
##   SIM_D         — L² separation (default: 0)
##   SIM_STRUCTURE — trend structure: smooth, cross, rate (default: smooth)
##   SIM_N         — series length (default: 500)
##   SIM_SNR       — signal-to-noise ratio (default: 1.5)
##   SIM_PHI       — AR(1) coefficient, fixed (default: 0.5)
##   SIM_S         — number of replicates (default: 200)
##   SIM_SEED      — base seed (default: 2847)
##   SIM_OUT_DIR   — output directory (default: /jobs/output)
##
## Test locally:
##   SIM_D=0 SIM_STRUCTURE=smooth SIM_N=500 SIM_SNR=1.5 SIM_S=5 \
##     SIM_OUT_DIR=dev/sims/power/results/raw \
##     Rscript dev/sims/power/tide/sim.R

library(lomad)
library(dplyr)

# ---- Parameters from environment -------------------------------------------

d         <- as.numeric(Sys.getenv("SIM_D",         "0"))
structure <- Sys.getenv("SIM_STRUCTURE", "smooth")
n         <- as.integer(Sys.getenv("SIM_N",         "500"))
snr       <- as.numeric(Sys.getenv("SIM_SNR",       "1.5"))
phi       <- as.numeric(Sys.getenv("SIM_PHI",       "0.5"))
S         <- as.integer(Sys.getenv("SIM_S",         "200"))
seed0     <- as.integer(Sys.getenv("SIM_SEED",      "2847"))
out_dir   <- Sys.getenv("SIM_OUT_DIR", "/jobs/output")

# ---- Fixed parameters (must match template.R) --------------------------------

alpha <- 0.05

h_win <- max(5L, floor(n / 200L))
s_win <- min(60L * h_win, floor(n / 4L))

struct_params <- list(
  dist   = list(),
  smooth = list(bw = 50),
  cross  = list(bw = 50),
  rate   = list(rate = 0.01)
)

# ---- Copy run_rep() from template.R -----------------------------------------

run_rep <- function(d, structure, seed) {
  set.seed(seed)

  trends <- do.call(sim_trends,
    c(list(n = n, d = d, method = structure, seed = seed),
      struct_params[[structure]]))

  sim <- sim_noise_pair(trends, h = h_win, lambda_target = snr,
                        ar.coefs = phi, seed = seed + 1L)

  fit <- tryCatch(
    lomad_fit(sim$y1, sim$y2, h = h_win, s = s_win),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(data.frame(d = d, structure = structure, seed = seed,
                      detected = NA, sensitivity = NA_real_, fdr = NA_real_,
                      n_flagged = NA_integer_))
  }

  tst <- lomad_test(fit, alpha = alpha)
  vi  <- fit$valid_idx
  rejected <- tst$rejected[vi]

  w <- if (!is.null(trends$w)) trends$w[vi] else NULL
  w_thresh <- 0.1

  detected    <- any(rejected, na.rm = TRUE)
  sensitivity <- if (!is.null(w) && any(abs(w - 1) > w_thresh))
    mean(rejected[abs(w - 1) > w_thresh], na.rm = TRUE) else NA_real_
  fdr_val     <- if (!is.null(w) && any(rejected, na.rm = TRUE))
    mean(abs(w[rejected] - 1) <= w_thresh, na.rm = TRUE) else NA_real_

  data.frame(d = d, structure = structure, seed = seed,
             detected = detected, sensitivity = sensitivity,
             fdr = fdr_val, n_flagged = sum(rejected, na.rm = TRUE))
}

# ---- Simulation loop ---------------------------------------------------------

set.seed(seed0 + as.integer(d * 100))
seeds <- sample.int(1e6, S)

results <- bind_rows(lapply(seq_len(S), function(i) {
  run_rep(d, structure, seeds[i])
}))

# ---- Summary -----------------------------------------------------------------

results_summary <- results |>
  summarise(
    S           = n(),
    detection   = mean(detected, na.rm = TRUE),
    sensitivity = mean(sensitivity, na.rm = TRUE),
    fdr         = mean(fdr, na.rm = TRUE)
  )

# ---- Save results ------------------------------------------------------------

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
filename <- sprintf("%s_d%s_n%d_snr%s_phi%s.rds",
                    structure,
                    gsub("\\.", "-", format(d,   nsmall = 1)),
                    n,
                    gsub("\\.", "-", format(snr, nsmall = 1)),
                    gsub("\\.", "-", format(phi, nsmall = 1)))
saveRDS(list(d = d, structure = structure, n = n, snr = snr, phi = phi,
             S = S, seed0 = seed0,
             results = results, results_summary = results_summary),
        file.path(out_dir, filename))
