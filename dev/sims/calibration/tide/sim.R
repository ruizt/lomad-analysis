## Calibration simulation — Kubernetes container entrypoint
##
## One container = all S replicates for one value of d.
## Parameters are passed as environment variables by the job spec.
##
## Environment variables:
##   SIM_D       — L² separation (default: 0)
##   SIM_S       — number of replicates (default: 200)
##   SIM_SEED    — base seed (default: 4853)
##   SIM_OUT_DIR — output directory (default: /jobs/output)
##
## Test locally before submitting to Tide:
##   SIM_D=0 SIM_S=5 SIM_SEED=4853 SIM_OUT_DIR=dev/sims/calibration/results/raw \
##     Rscript dev/sims/calibration/tide/sim.R

library(lomad)
library(dplyr)

# ---- Parameters from environment -------------------------------------------

d       <- as.numeric(Sys.getenv("SIM_D",       "0"))
S       <- as.integer(Sys.getenv("SIM_S",       "200"))
seed0   <- as.integer(Sys.getenv("SIM_SEED",    "4853"))
out_dir <- Sys.getenv("SIM_OUT_DIR", "/jobs/output")

# ---- Fixed parameters (must match template.R) -----------------------------

n     <- 1000
phi   <- 0.5
snr   <- 1
alpha <- 0.05

h_win <- 10
s_win <- 50

# ---- Single-replicate function (copied from template.R) -------------------

run_rep <- function(d, seed) {
  trends <- sim_trends(n = n, d = d, method = "dist", seed = seed)
  sim    <- suppressMessages(
    sim_noise_pair(trends, h = h_win, s = s_win, lambda_target = snr,
                   ar.coefs = phi, seed = seed + 1L)
  )

  # 1. CLT test: estimated pipeline
  fit <- suppressMessages(
    lomad_fit(sim$y1, sim$y2, h = h_win, s = s_win)
  )
  tst <- suppressMessages(lomad_test(fit, alpha = alpha))
  clt_rejected <- if (length(fit$valid_idx) > 0) {
    any(tst$rejected[fit$valid_idx], na.rm = TRUE)
  } else NA

  # 2. CLT test: oracle (true per-series noise parameters)
  fit_orc <- suppressMessages(
    lomad_fit(sim$y1, sim$y2, h = h_win, s = s_win,
              noise_override = list(
                list(ar = phi, sigma2 = sim$noise$series1$sigma^2),
                list(ar = phi, sigma2 = sim$noise$series2$sigma^2)
              ))
  )
  tst_orc <- lomad_test(fit_orc, alpha = alpha)
  oracle_rejected <- if (length(fit_orc$valid_idx) > 0) {
    any(tst_orc$rejected[fit_orc$valid_idx], na.rm = TRUE)
  } else NA

  # 3. Identity test (oracle, global Bonferroni p-value)
  sigma2_innov <- mean(c(sim$noise$series1$sigma,
                         sim$noise$series2$sigma)^2)
  ident <- suppressMessages(
    lomad_test_identity(sim$y1, sim$y2, q = h_win, alpha = alpha,
                        noise_override = list(ar = phi, sigma2 = sigma2_innov))
  )
  identity_rejected <- ident$global_p < alpha

  data.frame(d = d, seed = seed,
             clt_rejected      = clt_rejected,
             oracle_rejected   = oracle_rejected,
             identity_rejected = identity_rejected)
}

# ---- Simulation loop -------------------------------------------------------

set.seed(seed0 + as.integer(d * 100))
seeds <- sample.int(1e6, S)

results <- bind_rows(lapply(seq_len(S), function(i) {
  run_rep(d, seeds[i])
}))

# ---- Save results ----------------------------------------------------------

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
filename <- sprintf("d%s.rds", gsub("\\.", "-", format(d, nsmall = 1)))
saveRDS(list(d = d, S = S, seed0 = seed0, results = results),
        file.path(out_dir, filename))

cat(sprintf("Saved: %s\n", file.path(out_dir, filename)))
