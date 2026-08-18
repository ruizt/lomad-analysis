## Validation study — Kubernetes container entrypoint
##
## One container = all S replicates for one experiment.
## Parameters are passed as environment variables by the job spec.
##
## To test locally (requires lomad installed):
##   SIM_EXPERIMENT=clt-s80 SIM_S=5 SIM_SEED=7291 \
##     SIM_OUT_DIR=simulations/validation/results/_raw \
##     Rscript simulations/validation/tide/sim.R

# In the container, lomad is pre-installed. For local testing, install it
# first or use simulation-template.R with library(lomad) for development.
library(lomad)

# ---- Parameters from environment ---------------------------------------------

experiment <- Sys.getenv("SIM_EXPERIMENT", "clt-s80")
S          <- as.integer(Sys.getenv("SIM_S",       "200"))
seed0      <- as.integer(Sys.getenv("SIM_SEED",    "7291"))
out_dir    <- Sys.getenv("SIM_OUT_DIR", "/jobs/output")

# Parsed here rather than at dispatch: the oracle quantities below are built at
# the window size, so s_val has to be known before them.
parts    <- strsplit(experiment, "-")[[1]]
exp_type <- parts[1]
s_val    <- as.integer(sub("^s", "", parts[2]))

cat(sprintf("experiment = %s (type %s, s = %d), S = %d, seed0 = %d\n",
            experiment, exp_type, s_val, S, seed0))

# ---- Fixed parameters --------------------------------------------------------

n_obs <- 2000
h_win <- 5

# Shared trend via sim_trends (common trend, d = 0)
# Two trends that are equal under affine transformation: nu_2 = B_SCALE * nu_1
# with B_SCALE > 0, so H_0 holds by construction while tau_2^2 = B_SCALE^2
# tau_1^2. The equal-amplitude case B_SCALE = 1 is a special case in which the
# revised Proposition 1 collapses to its predecessor -- and in which a
# cross-pairing error in V is undetectable, since tau_2^2 sigma_1^4 L_1 and
# tau_1^2 sigma_1^4 L_1 coincide. At B_SCALE = 2 the correct and swapped forms
# differ by 53%, so the study can actually see the error it is most at risk of.
#
# Only the scale matters. Correlation is location invariant, so an offset
# nu_2 = a + b nu_1 changes R, tau, rho and V not at all; the location half of
# the affine null needs no simulation.
B_SCALE <- 2

tr     <- sim_trends(n = n_obs, d = 0, nb = 25, sd0 = 50, p = 1.5, seed = 5381)
trend  <- tr$x1
trend2 <- B_SCALE * trend

# Noiseless MA-smoothed trend — used for oracle tau_sq everywhere.
# Because d = 0 (shared trend), this equals filter(trend, ...) directly;
# no bias correction is needed.
ma_trend  <- as.numeric(stats::filter(trend,  rep(1 / h_win, h_win), sides = 1))
ma_trend2 <- as.numeric(stats::filter(trend2, rep(1 / h_win, h_win), sides = 1))

# AR(1) noise
e2e_ar1 <- 0.5;  e2e_sd1 <- 0.8
e2e_ar2 <- 0.3;  e2e_sd2 <- 0.8

.ma_filter_acov <- lomad:::.ma_filter_acov

# Filtered autocovariances and oracle quantities
e2e_acov_raw1  <- arma_acov(e2e_ar1, numeric(0), e2e_sd1^2, lag_max = h_win + 100)
e2e_acov_raw2  <- arma_acov(e2e_ar2, numeric(0), e2e_sd2^2, lag_max = h_win + 100)
e2e_acov_filt1 <- .ma_filter_acov(e2e_acov_raw1, h_win, lag_max = 100)
e2e_acov_filt2 <- .ma_filter_acov(e2e_acov_raw2, h_win, lag_max = 100)
e2e_acov_filt1 <- e2e_acov_filt1[!is.na(e2e_acov_filt1)]
e2e_acov_filt2 <- e2e_acov_filt2[!is.na(e2e_acov_filt2)]
e2e_ml         <- min(length(e2e_acov_filt1), length(e2e_acov_filt2))
e2e_acov_filt1 <- e2e_acov_filt1[1:e2e_ml]
e2e_acov_filt2 <- e2e_acov_filt2[1:e2e_ml]

e2e_sigma1     <- e2e_acov_filt1[1]
e2e_sigma2     <- e2e_acov_filt2[1]
e2e_sums       <- acov_sums(e2e_acov_filt1, e2e_acov_filt2)
e2e_s          <- s_val
# Use noiseless ma_trend directly — no bias correction needed in oracle setting
e2e_tau_sq     <- compute_tau_sq(ma_trend,  e2e_s)
e2e_tau2_sq    <- compute_tau_sq(ma_trend2, e2e_s)
e2e_rho_oracle <- compute_rho(e2e_tau_sq, e2e_tau2_sq, e2e_sigma1, e2e_sigma2)
e2e_V_oracle   <- compute_V(e2e_tau_sq, e2e_tau2_sq, e2e_sigma1, e2e_sigma2,
                              e2e_sums$L1, e2e_sums$L2,
                              e2e_sums$Q1, e2e_sums$Q2, e2e_sums$Q12)
# Four evaluation points at even 500-unit spacing, chosen by position rather
# than by their local quantities so the set carries no selection. They span
# lambda_1 from 0.024 to 1.402 and rho from 0.059 to 0.730. t = 1900 sits in a
# near-flat window where the test has little power; a dense sweep over every
# 5th t shows coverage there is closest to nominal, so including it costs
# nothing in calibration terms and widens the range on display.
e2e_eval_pts <- c(400, 900, 1400, 1900)
e2e_band_pts <- seq(h_win + e2e_s, n_obs, by = 5L)
E2E_RHO_REPS <- 100L

# ---- run_rep functions -------------------------------------------------------

run_rep_e2e <- function(s, seed, eval_pts, alpha = 0.05) {
  set.seed(seed)
  z1 <- arima.sim(model = list(ar = e2e_ar1), n = n_obs, sd = e2e_sd1)
  z2 <- arima.sim(model = list(ar = e2e_ar2), n = n_obs, sd = e2e_sd2)
  y1 <- trend + z1
  y2 <- trend2 + z2

  fit <- suppressMessages(lomad_fit(y1, y2, h = h_win, s = s))
  tst <- suppressMessages(lomad_test(fit, alpha = alpha))

  # One fit serves both: the four evaluation points and the dense band grid.
  list(
    eval = data.frame(
      t        = eval_pts,
      R        = fit$R[eval_pts],
      rho_hat  = fit$rho[eval_pts],
      V_hat    = fit$V[eval_pts],
      Z_pipe   = tst$Z[eval_pts],
      rejected = tst$rejected[eval_pts]
    ),
    band = data.frame(
      R      = fit$R[e2e_band_pts],
      Z_pipe = tst$Z[e2e_band_pts]
    )
  )
}

# ---- Dispatch ----------------------------------------------------------------

set.seed(seed0 + s_val)
seeds <- sample.int(1e6, S)

# ---- Run simulation ----------------------------------------------------------

if (exp_type == "e2e") {
  e2e_results <- vector("list", S)
  for (i in seq_along(seeds)) {
    e2e_results[[i]] <- run_rep_e2e(s_val, seeds[i], e2e_eval_pts)
    if (i %% 100 == 0) cat(sprintf("  %d / %d\n", i, S))
  }

  n_pts   <- length(e2e_eval_pts)
  R_mat   <- matrix(NA_real_, S, n_pts)
  rho_mat <- matrix(NA_real_, S, n_pts)
  V_mat   <- matrix(NA_real_, S, n_pts)
  Z_mat   <- matrix(NA_real_, S, n_pts)

  for (i in seq_along(e2e_results)) {
    R_mat[i, ]   <- e2e_results[[i]]$eval$R
    rho_mat[i, ] <- e2e_results[[i]]$eval$rho_hat
    V_mat[i, ]   <- e2e_results[[i]]$eval$V_hat
    Z_mat[i, ]   <- e2e_results[[i]]$eval$Z_pipe
  }
  e2e_results_band <- lapply(e2e_results, `[[`, "band")

  # Per-window rejection rates over the dense grid. The test forms
  # p = pnorm(Z) and rejects when p <= alpha, so only the lower tail fires;
  # two-sided coverage would let an inflated lower tail cancel a deflated
  # upper one.
  bR <- matrix(NA_real_, S, length(e2e_band_pts))
  bZ <- bR
  for (i in seq_along(e2e_results_band)) {
    bR[i, ] <- e2e_results_band[[i]]$R
    bZ[i, ] <- e2e_results_band[[i]]$Z_pipe
  }
  e2e_band <- do.call(rbind, lapply(seq_along(e2e_band_pts), function(j) {
    tt <- e2e_band_pts[j]
    Zo <- sqrt(s_val) * (bR[, j] - e2e_rho_oracle[tt]) / sqrt(e2e_V_oracle[tt])
    Zo <- Zo[is.finite(Zo)]
    Zp <- bZ[, j]; Zp <- Zp[is.finite(Zp)]
    f <- function(Z, ty) data.frame(
      t = tt, lambda1 = e2e_tau_sq[tt] / e2e_sigma1, rho = e2e_rho_oracle[tt],
      # Moments of R over replicates, so rho and V accuracy come from this
      # experiment rather than separate oracle runs. The mean uses the first
      # E2E_RHO_REPS only: over all S its Monte Carlo error is thinner than the
      # plotted line and the empirical curve vanishes under the theoretical.
      R_mean = mean(bR[seq_len(min(E2E_RHO_REPS, nrow(bR))), j], na.rm = TRUE),
      R_mean_n = min(E2E_RHO_REPS, nrow(bR)),
      R_var  = stats::var(bR[, j], na.rm = TRUE),
      V_th   = e2e_V_oracle[tt],
      cov = mean(abs(Z) < stats::qnorm(0.975)),
      t05 = mean(Z < stats::qnorm(0.05)),
      t01 = mean(Z < stats::qnorm(0.01)),
      n = length(Z), type = ty)
    rbind(f(Zo, "Oracle"), f(Zp, "End-to-end"))
  }))

  result <- list(
    experiment  = experiment, s = s_val, S = S, seed0 = seed0, h = h_win,
    sigma2_eta = c(e2e_sigma1, e2e_sigma2),
    eval_pts    = e2e_eval_pts,
    R_mat       = R_mat,
    rho_est_mat = rho_mat,
    V_est_mat   = V_mat,
    Z_est_mat   = Z_mat,
    b_scale     = B_SCALE,
    rho_oracle  = e2e_rho_oracle[e2e_eval_pts],
    V_oracle    = e2e_V_oracle[e2e_eval_pts],
    # Full-length, so the figure can draw local SNR over t without repeating
    # the noise constants. lambda_2 is a fixed multiple of lambda_1 here --
    # b^2 sigma_1^2 / sigma_2^2, since both trends are one curve up to scale.
    lambda1     = e2e_tau_sq  / e2e_sigma1,
    lambda2     = e2e_tau2_sq / e2e_sigma2,
    tau2_1      = e2e_tau_sq,
    tau2_2      = e2e_tau2_sq,
    band        = e2e_band
  )

} else {
  stop("Unknown experiment type: ", exp_type)
}

# ---- Save results ------------------------------------------------------------

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
filename <- paste0(experiment, ".rds")
saveRDS(result, file.path(out_dir, filename))
cat(sprintf("Saved: %s\n", file.path(out_dir, filename)))
