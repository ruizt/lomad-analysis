## simulation-template.R — local illustration of the power study
##
## Runs S replicates for each (structure, d, s, phi, snr) cell so the simulation
## can be inspected and reasoned about locally, at a scale that runs in seconds
## rather than on the cluster.
##
## NOT the source of truth. tide/sim.R is what actually runs on Tide and
## produces the archived results; run_rep() below deliberately mirrors it so the
## illustration is faithful. If you change the simulation logic, change
## tide/sim.R first, then mirror it here. The two differ only in how they take
## parameters -- arguments here, environment variables there -- and in the loop
## that drives them.
##
## Writes nothing: only tide/sim.R and collect-results.R touch results/.
##
## Usage (from repo root):
##   source("simulations/power/simulation-template.R")

library(lomad)
library(dplyr)

# ---- Parameters -------------------------------------------------------------

alpha <- 0.05
S     <- 20

# The window s is the design factor and n follows it. Accumulated affine drift
# depends on the number of windows n/s, so holding that ratio fixed keeps the
# drift identical across window sizes rather than confounding the two.
N_OVER_S <- 25L
h_win    <- 5L
nb       <- 151L        # holds the shortest basis period at s/3 for any s

affine_cap <- 0.015     # per-window; the most that holds H_0 in the worst cell
affine_bw  <- 0.5

s_vals   <- c(50L, 100L, 150L)
phi_vals <- c(0.3, 0.5, 0.8)
snr_vals <- c(0.5, 1.5)
structs  <- c("rs", "rm", "fr")

# d is a generator knob, not an effect size: the same d gives ~2x different
# local separation across structures because they distribute it differently in
# time. It is therefore scaled per structure onto a common realized delta_t,
# and never appears in a figure. Base grid gives median delta_t of roughly
# 0, 0.20 and 0.45; see _notes/delta-calibration.md.
d_base   <- c(0, 0.60, 1.55)
d_factor <- c(rs = 1.00, rm = 0.60, fr = 0.45)
d_vals   <- lapply(d_factor, function(f) round(f * d_base, 2))

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

THIN <- 5L              # windows overlap by s - 1, so neighbours are redundant

# ---- Single-replicate function ----------------------------------------------

run_rep <- function(d, struct, s_win, phi, snr, seed, oracle = FALSE) {
  set.seed(seed)
  n <- N_OVER_S * s_win

  trends <- do.call(sim_trends,
    c(list(n = n, d = d, method = struct, nb = nb, seed = seed,
           affine_s = s_win, affine_cap = affine_cap, affine_bw = affine_bw),
      struct_params[[struct]]))

  sim <- sim_noise_pair(trends, h = h_win, lambda_target = snr,
                        ar.coefs = phi, seed = seed + 1L)

  # Oracle: bypass noise estimation with the true AR params
  noise_ov <- NULL
  if (oracle) {
    z1 <- sim$y1 - sim$x1
    innov1 <- z1[-1] - phi * z1[-length(z1)]
    noise_ov <- list(ar = phi, sigma2 = var(innov1))
  }

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

  b_range <- diff(range(trends$b))

  fit <- tryCatch(
    lomad_fit(sim$y1, sim$y2, h = h_win, s = s_win, noise_override = noise_ov),
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

  vi   <- fit$valid_idx
  keep <- vi[seq(1L, length(vi), by = THIN)]

  # Schema matches tide/sim.R exactly: collect-results.R assumes it. p_raw is
  # kept so the study can be rethresholded at another alpha without regenerating
  # anything; rejected is what BY gave at `alpha`.
  list(
    summary = data.frame(d = d, struct = struct, s = s_win, n = n,
                         phi = phi, snr = snr, seed = seed,
                         detected = any(tst$rejected, na.rm = TRUE),
                         n_win = length(vi),
                         delta_med = median(delta_t[vi], na.rm = TRUE),
                         lambda1 = lam1, lambda2 = lam2, b_range = b_range),
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

# summary output
run_rep(d = 0.6, struct = "rs", s_win = 100L, phi = 0.5, snr = 1.5,
        seed = 123)$summary

# per-window output: this is what the local power curves are built from
run_rep(d = 0.6, struct = "rs", s_win = 100L, phi = 0.5, snr = 1.5,
        seed = 123)$windows |> head()

# ---- Main loop --------------------------------------------------------------

set.seed(2847)
all_seeds <- sample.int(1e6, S)

results <- lapply(structs, function(struct) {
  lapply(d_vals[[struct]], function(d) {
    lapply(s_vals, function(s_win) {
      lapply(phi_vals, function(phi) {
        lapply(snr_vals, function(snr) {
          reps <- lapply(all_seeds,
                         function(sd) run_rep(d, struct, s_win, phi, snr, sd))
          bind_rows(lapply(reps, `[[`, "summary"))
        }) |> bind_rows()
      }) |> bind_rows()
    }) |> bind_rows()
  }) |> bind_rows()
}) |> bind_rows()
