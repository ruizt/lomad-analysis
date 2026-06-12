## Power study — Kubernetes container entrypoint
##
## One container = all S replicates for one (structure, d, n, snr) combination.
## Parameters are passed as environment variables by the job spec.
##
## Outputs per combo:
##   {combo}.rds         — list with metadata + per-rep summary (detection rate)
##   {combo}-series.rds  — list keyed by seed with w, vi, p_raw, p_adj, rejected
##
## Environment variables:
##   SIM_D         — L² separation (default: 0)
##   SIM_STRUCTURE — trend structure: smooth, cross, rate (default: smooth)
##   SIM_N         — series length (default: 500)
##   SIM_SNR       — signal-to-noise ratio (default: 1.5)
##   SIM_PHI       — AR(1) coefficient, fixed (default: 0.5)
##   SIM_S         — number of replicates (default: 200)
##   SIM_SEED      — base seed (default: 2847)
##   SIM_ORACLE    — use true noise params, bypassing estimation (default: FALSE)
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
struct    <- Sys.getenv("SIM_STRUCTURE", "smooth")
n         <- as.integer(Sys.getenv("SIM_N",         "500"))
snr       <- as.numeric(Sys.getenv("SIM_SNR",       "1.5"))
phi       <- as.numeric(Sys.getenv("SIM_PHI",       "0.5"))
S         <- as.integer(Sys.getenv("SIM_S",         "200"))
seed0     <- as.integer(Sys.getenv("SIM_SEED",      "2847"))
oracle    <- as.logical(Sys.getenv("SIM_ORACLE",    "FALSE"))
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

run_rep <- function(d, struct, seed) {
  set.seed(seed)

  trends <- do.call(sim_trends,
    c(list(n = n, d = d, method = struct, seed = seed),
      struct_params[[struct]]))

  w <- if (!is.null(trends$w)) trends$w else NULL

  sim <- sim_noise_pair(trends, h = h_win, lambda_target = snr,
                        ar.coefs = phi, seed = seed + 1L)

  # Oracle: bypass noise estimation with true AR params
  noise_ov <- NULL
  if (oracle) {
    z1 <- sim$y1 - sim$x1
    innov1 <- z1[-1] - phi * z1[-length(z1)]
    noise_ov <- list(ar = phi, sigma2 = var(innov1))
  }

  fit <- tryCatch(
    lomad_fit(sim$y1, sim$y2, h = h_win, s = s_win, noise_override = noise_ov),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(list(
      summary = data.frame(d = d, struct = struct, n = n,
                           phi = phi, snr = snr, seed = seed,
                           detected = NA),
      series = NULL
    ))
  }

  tst <- lomad_test(fit, alpha = alpha)

  list(
    summary = data.frame(d = d, n = n, phi = phi, snr = snr,
                         struct = struct, seed = seed,
                         detected = any(tst$rejected, na.rm = TRUE)),
    series = list(w = w,
                  sep = abs(trends$x1 - trends$x2),
                  vi = fit$valid_idx,
                  p_raw = tst$p_values,
                  p_adj = tst$p_adj,
                  rejected = tst$rejected)
  )
}

# ---- Simulation loop ---------------------------------------------------------

set.seed(seed0 + as.integer(d * 100))
seeds <- sample.int(1e6, S)

reps <- lapply(seq_len(S), function(i) run_rep(d, struct, seeds[i]))

results <- bind_rows(lapply(reps, `[[`, "summary"))
series  <- setNames(lapply(reps, `[[`, "series"), seeds)

# ---- Summary -----------------------------------------------------------------

results_summary <- results |>
  summarise(
    S         = n(),
    detection = mean(detected, na.rm = TRUE)
  )

# ---- Save results ------------------------------------------------------------

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
basename <- sprintf("%s_d%s_n%d_snr%s_phi%s%s",
                    struct,
                    gsub("\\.", "-", format(d,   nsmall = 1)),
                    n,
                    gsub("\\.", "-", format(snr, nsmall = 1)),
                    gsub("\\.", "-", format(phi, nsmall = 1)),
                    if (oracle) "-oracle" else "")

saveRDS(list(d = d, struct = struct, n = n, snr = snr, phi = phi,
             S = S, seed0 = seed0,
             results = results, results_summary = results_summary),
        file.path(out_dir, paste0(basename, ".rds")))

saveRDS(series, file.path(out_dir, paste0(basename, "-series.rds")))
