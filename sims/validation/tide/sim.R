## Validation study — Kubernetes container entrypoint
##
## One container = all S replicates for one experiment.
## Parameters are passed as environment variables by the job spec.
##
## To test locally (requires lomad installed):
##   SIM_EXPERIMENT=clt-s80 SIM_S=5 SIM_SEED=7291 \
##     SIM_OUT_DIR=sims/validation/results/raw \
##     Rscript sims/validation/tide/sim.R

# In the container, lomad is pre-installed. For local testing, install it
# first or use template.R with library(lomad) for development.
library(lomad)

# ---- Parameters from environment ---------------------------------------------

experiment <- Sys.getenv("SIM_EXPERIMENT", "clt-s80")
S          <- as.integer(Sys.getenv("SIM_S",       "200"))
seed0      <- as.integer(Sys.getenv("SIM_SEED",    "7291"))
out_dir    <- Sys.getenv("SIM_OUT_DIR", "/jobs/output")

cat(sprintf("experiment = %s, S = %d, seed0 = %d\n", experiment, S, seed0))

# ---- Fixed parameters --------------------------------------------------------

n_obs <- 2000
h_win <- 20

# Shared trend via sim_trends (common trend, d = 0)
tr     <- sim_trends(n = n_obs, d = 0, nb = 25, sd0 = 50, p = 1.5, seed = 5381)
trend  <- tr$x1

# Noiseless MA-smoothed trend — used for oracle tau_sq everywhere.
# Because d = 0 (shared trend), this equals filter(trend, ...) directly;
# no bias correction is needed.
ma_trend <- as.numeric(stats::filter(trend, rep(1 / h_win, h_win), sides = 1))

# Oracle ARMA(1,1) noise (Figures 1–2)
oracle_ar1 <- 0.6;  oracle_ma1 <- 0.3;   oracle_sd1 <- 0.8
oracle_ar2 <- 0.4;  oracle_ma2 <- -0.2;  oracle_sd2 <- 1.0

# End-to-end AR(1) noise (Figure 3)
e2e_ar1 <- 0.5;  e2e_sd1 <- 0.8
e2e_ar2 <- 0.3;  e2e_sd2 <- 0.8

# Precompute filtered autocovariances for oracle experiments
.ma_filter_acov <- lomad:::.ma_filter_acov

acov_raw1  <- arma_acov(oracle_ar1, oracle_ma1, oracle_sd1^2, lag_max = 200)
acov_raw2  <- arma_acov(oracle_ar2, oracle_ma2, oracle_sd2^2, lag_max = 200)
acov_eta1  <- .ma_filter_acov(acov_raw1, h_win, lag_max = 200)
acov_eta2  <- .ma_filter_acov(acov_raw2, h_win, lag_max = 200)
sigma1_sq  <- acov_eta1[1]
sigma2_sq  <- acov_eta2[1]
acov_eta1  <- acov_eta1[!is.na(acov_eta1)]
acov_eta2  <- acov_eta2[!is.na(acov_eta2)]
ml         <- min(length(acov_eta1), length(acov_eta2))
acov_eta1  <- acov_eta1[1:ml]
acov_eta2  <- acov_eta2[1:ml]
cov_sums   <- acov_sums(acov_eta1, acov_eta2)

# Precompute e2e oracle quantities
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
e2e_s          <- 150
# Use noiseless ma_trend directly — no bias correction needed in oracle setting
e2e_tau_sq     <- compute_tau_sq(ma_trend, e2e_s)
e2e_rho_oracle <- compute_rho(e2e_tau_sq, e2e_sigma1, e2e_sigma2)
e2e_V_oracle   <- compute_V(e2e_tau_sq, e2e_sigma1, e2e_sigma2,
                              e2e_sums$L1, e2e_sums$L2,
                              e2e_sums$Q1, e2e_sums$Q2, e2e_sums$Q12)
e2e_eval_pts   <- c(500, 850, 1000, 1400, 1600)

# ---- run_rep functions -------------------------------------------------------

run_rep_clt <- function(s, seed, eval_t = 600) {
  set.seed(seed)
  z1 <- arima.sim(model = list(ar = oracle_ar1, ma = oracle_ma1),
                  n = n_obs, sd = oracle_sd1)
  z2 <- arima.sim(model = list(ar = oracle_ar2, ma = oracle_ma2),
                  n = n_obs, sd = oracle_sd2)
  y1 <- trend + z1
  y2 <- trend + z2
  m1 <- as.numeric(stats::filter(y1, rep(1 / h_win, h_win), sides = 1))
  m2 <- as.numeric(stats::filter(y2, rep(1 / h_win, h_win), sides = 1))
  idx <- (eval_t - s + 1):eval_t
  if (any(is.na(m1[idx])) || any(is.na(m2[idx]))) return(NA_real_)
  cor(m1[idx], m2[idx])
}

run_rep_rho <- function(s, seed) {
  set.seed(seed)
  z1 <- arima.sim(model = list(ar = oracle_ar1, ma = oracle_ma1),
                  n = n_obs, sd = oracle_sd1)
  z2 <- arima.sim(model = list(ar = oracle_ar2, ma = oracle_ma2),
                  n = n_obs, sd = oracle_sd2)
  y1 <- trend + z1
  y2 <- trend + z2
  m1 <- as.numeric(stats::filter(y1, rep(1 / h_win, h_win), sides = 1))
  m2 <- as.numeric(stats::filter(y2, rep(1 / h_win, h_win), sides = 1))

  R <- rep(NA_real_, n_obs)
  for (t in s:n_obs) {
    idx <- (t - s + 1):t
    a <- m1[idx]; b <- m2[idx]
    if (any(is.na(a)) || any(is.na(b))) next
    R[t] <- cor(a, b)
  }
  R
}

run_rep_var <- function(s, seed, eval_pts) {
  set.seed(seed)
  z1 <- arima.sim(model = list(ar = oracle_ar1, ma = oracle_ma1),
                  n = n_obs, sd = oracle_sd1)
  z2 <- arima.sim(model = list(ar = oracle_ar2, ma = oracle_ma2),
                  n = n_obs, sd = oracle_sd2)
  y1 <- trend + z1
  y2 <- trend + z2
  m1 <- as.numeric(stats::filter(y1, rep(1 / h_win, h_win), sides = 1))
  m2 <- as.numeric(stats::filter(y2, rep(1 / h_win, h_win), sides = 1))

  R <- numeric(length(eval_pts))
  for (j in seq_along(eval_pts)) {
    t0  <- eval_pts[j]
    idx <- (t0 - s + 1):t0
    a <- m1[idx]; b <- m2[idx]
    R[j] <- if (any(is.na(a)) || any(is.na(b))) NA_real_ else cor(a, b)
  }
  R
}

run_rep_e2e <- function(s, seed, eval_pts, alpha = 0.05) {
  set.seed(seed)
  z1 <- arima.sim(model = list(ar = e2e_ar1), n = n_obs, sd = e2e_sd1)
  z2 <- arima.sim(model = list(ar = e2e_ar2), n = n_obs, sd = e2e_sd2)
  y1 <- trend + z1
  y2 <- trend + z2

  fit <- suppressMessages(lomad_fit(y1, y2, method = "clt", h = h_win, s = s))
  tst <- suppressMessages(lomad_test(fit, alpha = alpha))

  data.frame(
    t        = eval_pts,
    R        = fit$R[eval_pts],
    rho_hat  = fit$rho[eval_pts],
    V_hat    = fit$V[eval_pts],
    Z_pipe   = tst$Z[eval_pts],
    rejected = tst$rejected[eval_pts]
  )
}

# ---- Parse experiment ID and dispatch ----------------------------------------

parts <- strsplit(experiment, "-")[[1]]
exp_type <- parts[1]
s_val    <- as.integer(sub("^s", "", parts[2]))

cat(sprintf("Experiment type: %s, s = %d, S = %d reps\n", exp_type, s_val, S))

set.seed(seed0 + s_val)
seeds <- sample.int(1e6, S)

# ---- Run simulation ----------------------------------------------------------

if (exp_type == "clt") {
  R_vec <- vapply(seeds, function(sd) {
    run_rep_clt(s_val, sd, eval_t = 600)
  }, numeric(1))

  tau_sq <- compute_tau_sq(ma_trend, s_val)
  rho_t  <- compute_rho(tau_sq, sigma1_sq, sigma2_sq)
  V_t    <- compute_V(tau_sq, sigma1_sq, sigma2_sq,
                       cov_sums$L1, cov_sums$L2,
                       cov_sums$Q1, cov_sums$Q2, cov_sums$Q12)

  result <- list(
    experiment = experiment, s = s_val, S = S, seed0 = seed0,
    eval_t     = 600,
    R_vec      = R_vec,
    rho_oracle = rho_t[600],
    V_oracle   = V_t[600]
  )

} else if (exp_type == "rho") {
  R_accum <- rep(0, n_obs)
  R_count <- rep(0L, n_obs)

  for (i in seq_along(seeds)) {
    R_rep <- run_rep_rho(s_val, seeds[i])
    ok <- !is.na(R_rep)
    R_accum[ok] <- R_accum[ok] + R_rep[ok]
    R_count[ok] <- R_count[ok] + 1L
    if (i %% 100 == 0) cat(sprintf("  %d / %d\n", i, S))
  }

  R_mean <- ifelse(R_count > 0, R_accum / R_count, NA_real_)
  tau_sq <- compute_tau_sq(ma_trend, s_val)
  rho_th <- compute_rho(tau_sq, sigma1_sq, sigma2_sq)

  result <- list(
    experiment = experiment, s = s_val, S = S, seed0 = seed0,
    R_mean     = R_mean,
    rho_th     = rho_th
  )

} else if (exp_type == "var") {
  eval_pts <- seq(s_val + h_win, n_obs, by = 20)

  R_mat <- matrix(NA_real_, S, length(eval_pts))
  for (i in seq_along(seeds)) {
    R_mat[i, ] <- run_rep_var(s_val, seeds[i], eval_pts)
    if (i %% 100 == 0) cat(sprintf("  %d / %d\n", i, S))
  }

  tau_sq <- compute_tau_sq(ma_trend, s_val)
  V_th   <- compute_V(tau_sq, sigma1_sq, sigma2_sq,
                       cov_sums$L1, cov_sums$L2,
                       cov_sums$Q1, cov_sums$Q2, cov_sums$Q12)

  result <- list(
    experiment = experiment, s = s_val, S = S, seed0 = seed0,
    eval_pts   = eval_pts,
    R_mat      = R_mat,
    V_theory   = V_th[eval_pts]
  )

} else if (exp_type == "e2e") {
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
    R_mat[i, ]   <- e2e_results[[i]]$R
    rho_mat[i, ] <- e2e_results[[i]]$rho_hat
    V_mat[i, ]   <- e2e_results[[i]]$V_hat
    Z_mat[i, ]   <- e2e_results[[i]]$Z_pipe
  }

  result <- list(
    experiment  = experiment, s = s_val, S = S, seed0 = seed0,
    eval_pts    = e2e_eval_pts,
    R_mat       = R_mat,
    rho_est_mat = rho_mat,
    V_est_mat   = V_mat,
    Z_est_mat   = Z_mat,
    rho_oracle  = e2e_rho_oracle[e2e_eval_pts],
    V_oracle    = e2e_V_oracle[e2e_eval_pts]
  )

} else {
  stop("Unknown experiment type: ", exp_type)
}

# ---- Save results ------------------------------------------------------------

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
filename <- paste0(experiment, ".rds")
saveRDS(result, file.path(out_dir, filename))
cat(sprintf("Saved: %s\n", file.path(out_dir, filename)))
