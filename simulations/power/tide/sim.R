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
##     SIM_OUT_DIR=simulations/power/results/_raw \
##     Rscript simulations/power/tide/sim.R

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

# ---- Fixed parameters (must match simulation-template.R) --------------------------------

alpha <- 0.05

h_win <- max(5L, floor(n / 200L))
s_win <- min(60L * h_win, floor(n / 4L))

struct_params <- list(
  dist   = list(),
  smooth = list(bw = 50),
  cross  = list(bw = 50),
  # bump = "gaussian": the shape-2 gamma default puts a corner at each event
  # onset, which difference-based noise estimation cannot cancel (Hall and Van
  # Keilegom 2003 require a bounded derivative). The leftover biases the
  # residual autocovariance upward, and near the unit root that bias is
  # amplified by ~2/(1-phi)^2, so at phi = 0.8 the fixed-rate structure loses
  # almost all detection. The gaussian pulse is matched on width and smooth at
  # onset; it still trails the stochastic structures but is no longer anomalous.
  rate   = list(rate = 0.01, bump = "gaussian")
)

# ---- run_rep(): the authoritative simulation logic ---------------------------
# This is what runs on the cluster and produces the archived results.
# simulation-template.R mirrors it for local inspection; change this first, then mirror.

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

  # Realized affine effect size on the windows the test uses, from the true
  # noise-free trends smoothed exactly as the observed series are. d is a
  # design knob; delta_t = sqrt(1 - r_t^2) is what the test has power against,
  # so it is recorded rather than assumed.
  kern <- rep(1 / h_win, h_win)
  t1s  <- as.numeric(stats::filter(trends$x1, kern, sides = 1))
  t2s  <- as.numeric(stats::filter(trends$x2, kern, sides = 1))
  r_t  <- rep(NA_real_, n)
  for (tt in s_win:n) {
    ww <- (tt - s_win + 1L):tt
    a <- t1s[ww]; b <- t2s[ww]
    if (anyNA(a) || anyNA(b) || sd(a) == 0 || sd(b) == 0) next
    r_t[tt] <- suppressWarnings(stats::cor(a, b))
  }
  delta_t <- sqrt(pmax(0, 1 - r_t^2))

  # Realized per-series SNR on Proposition 1's definition: window signal
  # variance over smoothed noise variance, via the package's own Var_W.
  # Measured on s_win rather than sim_noise_pair()'s 2h calibration window, so
  # the level sits above `snr`; the ratio is the quantity of interest.
  eta1 <- as.numeric(stats::filter(sim$y1 - sim$x1, kern, sides = 1))
  eta2 <- as.numeric(stats::filter(sim$y2 - sim$x2, kern, sides = 1))
  tau_w <- function(z) mean(compute_tau_sq(z, s_win), na.rm = TRUE)
  lam1 <- tau_w(t1s) / var(eta1, na.rm = TRUE)
  lam2 <- tau_w(t2s) / var(eta2, na.rm = TRUE)

  fit <- tryCatch(
    lomad_fit(sim$y1, sim$y2, h = h_win, s = s_win, noise_override = noise_ov),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(list(
      summary = data.frame(d = d, struct = struct, n = n,
                           phi = phi, snr = snr, seed = seed,
                           detected = NA,
                           delta_sup = NA_real_, delta_bar = NA_real_,
                           lambda1 = NA_real_, lambda2 = NA_real_),
      series = NULL
    ))
  }

  tst <- lomad_test(fit, alpha = alpha)

  list(
    summary = data.frame(d = d, n = n, phi = phi, snr = snr,
                         struct = struct, seed = seed,
                         detected = any(tst$rejected, na.rm = TRUE),
                         delta_sup = suppressWarnings(max(delta_t, na.rm = TRUE)),
                         delta_bar = mean(delta_t, na.rm = TRUE),
                         lambda1 = lam1, lambda2 = lam2),
    series = list(w = w,
                  sep = abs(trends$x1 - trends$x2),
                  delta_t = delta_t,
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
