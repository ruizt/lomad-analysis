## Power study — Kubernetes container entrypoint
##
## One container = all S replicates for one (structure, d, s, snr, phi) cell.
## Parameters are passed as environment variables by the job spec.
##
## The study reports rejection as a function of *realized* local separation
## delta_t, not of d. d is a generator knob whose mapping to delta_t depends on
## the trend structure, so it is scaled per structure (see
## _notes/delta-calibration.md) and never appears in a figure.
##
## Outputs per cell:
##   {cell}.rds         — metadata + per-replicate summary
##   {cell}-windows.rds — per-window (delta_t, lambda, p_raw, rejected)
##
## Environment variables:
##   SIM_D         — separation knob, already scaled for the structure
##   SIM_STRUCTURE — trend structure: rs, rm, fr (default: rs)
##   SIM_S         — correlation window s (default: 100); n follows as N_OVER_S * s
##   SIM_SNR       — target signal-to-noise ratio (default: 1.5)
##   SIM_PHI       — AR(1) coefficient (default: 0.5)
##   SIM_REPS      — number of replicates (default: 200)
##   SIM_SEED      — base seed (default: 2847)
##   SIM_OUT_DIR   — output directory (default: /jobs/output)
##
## Test locally:
##   SIM_D=0.6 SIM_STRUCTURE=rs SIM_S=100 SIM_SNR=1.5 SIM_REPS=3 \
##     SIM_OUT_DIR=simulations/power/results/_raw \
##     Rscript simulations/power/tide/sim.R

library(lomad)
library(dplyr)

# ---- Parameters from environment --------------------------------------------

d       <- as.numeric(Sys.getenv("SIM_D",         "0"))
struct  <- Sys.getenv("SIM_STRUCTURE", "rs")
s_win   <- as.integer(Sys.getenv("SIM_S",         "100"))
snr     <- as.numeric(Sys.getenv("SIM_SNR",       "1.5"))
phi     <- as.numeric(Sys.getenv("SIM_PHI",       "0.5"))
S       <- as.integer(Sys.getenv("SIM_REPS",      "200"))
seed0   <- as.integer(Sys.getenv("SIM_SEED",      "2847"))
out_dir <- Sys.getenv("SIM_OUT_DIR", "/jobs/output")

# ---- Fixed parameters (must match simulation-template.R) --------------------

alpha <- 0.05

# The window s is the design factor; n follows it. Accumulated affine drift is
# a function of the number of windows, n/s, so holding that ratio fixed keeps
# the drift identical across window sizes instead of confounding the two. At
# n/s = 25 the coefficient b_t swings about 17% end to end at every s.
N_OVER_S <- 25L
n     <- N_OVER_S * s_win
h_win <- 5L

# Basis count that holds the shortest basis period at s/3 for any s: with
# n = 25 s, K = n / (s/3) = 75, so nb = 2K + 1 is the same at every cell. This
# is what keeps trend smoothness fixed relative to the window.
nb <- 151L

# Affine layer. The cap is per *window*, which is what keeps the pair locally
# affine similar while the coefficients accumulate across the series. 1.5% is
# the most that holds the null in the worst cell (high SNR); see
# _notes/delta-calibration.md.
affine_cap <- 0.015
affine_bw  <- 0.5

struct_params <- list(
  dist = list(),
  rs   = list(bw = 50),
  rm   = list(bw = 50),
  # bump = "gaussian": the shape-2 gamma default puts a corner at each event
  # onset, which difference-based noise estimation cannot cancel (Hall and Van
  # Keilegom 2003 require a bounded derivative). The leftover biases the
  # residual autocovariance upward, and near the unit root that bias is
  # amplified by ~2/(1-phi)^2, so at phi = 0.8 the fixed-rate structure loses
  # almost all detection. The gaussian pulse is matched on width and smooth at
  # onset; it still trails the stochastic structures but is no longer anomalous.
  fr   = list(rate = 0.01, bump = "gaussian")
)

# ---- run_rep(): the authoritative simulation logic --------------------------
# This is what runs on the cluster and produces the archived results.
# simulation-template.R mirrors it for local inspection; change this first,
# then mirror.

run_rep <- function(d, struct, seed) {
  set.seed(seed)

  trends <- do.call(sim_trends,
    c(list(n = n, d = d, method = struct, nb = nb, seed = seed,
           affine_s = s_win, affine_cap = affine_cap, affine_bw = affine_bw),
      struct_params[[struct]]))

  sim <- sim_noise_pair(trends, h = h_win, lambda_target = snr,
                        ar.coefs = phi, seed = seed + 1L)

  # Ground truth. delta_t is measured on the noise-free trends smoothed exactly
  # as the observed series are, by least squares within each window -- the same
  # sqrt(1 - r_t^2) the paper defines. The generating coefficients are NOT used
  # to remove the affine map: fixing the slope at its window average charges the
  # base separation twice once d > 0, and delta_t then exceeds 1.
  kern <- rep(1 / h_win, h_win)
  t1s  <- as.numeric(stats::filter(trends$x1, kern, sides = 1))
  t2s  <- as.numeric(stats::filter(trends$x2, kern, sides = 1))
  delta_t <- lomad:::.compute_delta(t1s, t2s, s_win)

  # Realized per-series SNR on Proposition 1's definition: window signal
  # variance over smoothed noise variance, via the package's own Var_W.
  # Measured on s_win rather than sim_noise_pair()'s 2h calibration window, so
  # the level sits above `snr`; the ratio is the quantity of interest.
  eta1 <- as.numeric(stats::filter(sim$y1 - sim$x1, kern, sides = 1))
  eta2 <- as.numeric(stats::filter(sim$y2 - sim$x2, kern, sides = 1))
  tau_w <- function(z) mean(compute_tau_sq(z, s_win), na.rm = TRUE)
  lam1 <- tau_w(t1s) / var(eta1, na.rm = TRUE)
  lam2 <- tau_w(t2s) / var(eta2, na.rm = TRUE)

  # Realized affine drift, so the layer's magnitude is recoverable per cell
  b_range <- diff(range(trends$b))

  fit <- tryCatch(
    lomad_fit(sim$y1, sim$y2, h = h_win, s = s_win),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(list(
      summary = data.frame(d = d, struct = struct, s = s_win, n = n,
                           phi = phi, snr = snr, seed = seed,
                           detected = NA, n_win = NA_integer_,
                           delta_med = NA_real_, lambda1 = NA_real_,
                           lambda2 = NA_real_, b_range = b_range),
      windows = NULL
    ))
  }

  tst <- lomad_test(fit, alpha = alpha)

  keep <- fit$valid_idx

  list(
    summary = data.frame(d = d, struct = struct, s = s_win, n = n,
                         phi = phi, snr = snr, seed = seed,
                         detected = any(tst$rejected, na.rm = TRUE),
                         n_win = length(keep),
                         delta_med = median(delta_t[keep], na.rm = TRUE),
                         lambda1 = lam1, lambda2 = lam2, b_range = b_range),
    # p_raw is kept so the whole study can be rethresholded at another alpha
    # without regenerating anything; rejected is what BY gave at `alpha`.
    windows = data.frame(
      seed     = seed,
      t        = keep,
      delta    = delta_t[keep],
      lambda   = pmin(fit$lambda1[keep], fit$lambda2[keep]),
      p_raw    = tst$p_values[keep],
      rejected = tst$rejected[keep]
    )
  )
}

# ---- Simulation loop --------------------------------------------------------

set.seed(seed0 + as.integer(d * 100))
seeds <- sample.int(1e6, S)

reps <- lapply(seq_len(S), function(i) run_rep(d, struct, seeds[i]))

results <- bind_rows(lapply(reps, `[[`, "summary"))
windows <- bind_rows(lapply(reps, `[[`, "windows"))

# ---- Save results -----------------------------------------------------------

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
basename <- sprintf("%s_d%s_s%d_snr%s_phi%s",
                    struct,
                    gsub("\\.", "-", format(d,   nsmall = 2)),
                    s_win,
                    gsub("\\.", "-", format(snr, nsmall = 1)),
                    gsub("\\.", "-", format(phi, nsmall = 1)))

saveRDS(list(d = d, struct = struct, s = s_win, n = n, snr = snr, phi = phi,
             S = S, seed0 = seed0,
             affine_cap = affine_cap, affine_bw = affine_bw,
             results = results),
        file.path(out_dir, paste0(basename, ".rds")))

saveRDS(windows, file.path(out_dir, paste0(basename, "-windows.rds")))
