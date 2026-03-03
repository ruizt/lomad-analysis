## fit_null_model.R
##
## fit_null_model():
##   Estimate null model and observed low-correlation statistics for two
##   *pre-smoothed* series, using an effective BY threshold (non-adaptive),
##   and compute:
##     - HAC-based approximate n_eff, tau^2, and
##     - asymptotic expectations under the null for:
##         entry rate, mean run length, fraction of time in state,
##       with respect to a user-specified null correlation rho0.
##
## bootstrap_local_corr():
##   Given the fitted null model, run a parametric bootstrap to obtain
##   null expectations, p-values, and refined n_eff and tau^2.
##
## -------------------------------------------------------------------
## fit_null_model.R
##
## INPUTS:
##   x1, x2 : numeric vectors (same length) = pre-smoothed series
##   q      : (optional) integer, moving-average window for trend smoothing;
##            if NULL, chosen automatically from series length
##   h      : (optional) integer, rolling correlation window length;
##            if NULL, chosen automatically to target n_eff ~ 30
##   alpha  : overall FDR level for BY threshold (default 0.05)
##   rho0   : null correlation value to test against (default 0)
##            (H0: rho = rho0 vs H1: rho < rho0)
##   max_pq : maximum AR and MA orders for ARMA order selection (default 2)

fit_null_model <- function(x1,
                           x2,
                           q = NULL,
                           h = NULL,
                           alpha = 0.05,
                           rho0 = 0,
                           max_pq = 2) {
  
  print("Starting null model estimation...")
  
  if (rho0 <= -0.999 || rho0 >= 0.999) {
    stop("rho0 must be strictly between -1 and 1.")
  }
  
  stopifnot(length(x1) == length(x2))
  x1 <- as.numeric(x1)
  x2 <- as.numeric(x2)
  Tn <- length(x1)
  
  ## Helper: choose default q and h based on theory if not supplied
  choose_q_h <- function(Tn) {
    ## q ≈ T/200 (at least 5), h ≈ 60*q (n_eff ~ 30), capped at T/4
    q_def <- max(5L, floor(Tn / 200L))
    h_raw <- 60L * q_def
    h_def <- min(max(h_raw, 10L * q_def), floor(Tn / 4L))
    list(q = q_def, h = h_def)
  }
  
  if (is.null(q) || is.null(h)) {
    print("Choosing q and h automatically...")
    defaults <- choose_q_h(Tn)
    if (is.null(q)) q <- defaults$q
    if (is.null(h)) h <- defaults$h
  }
  q <- as.integer(q)
  h <- as.integer(h)
  if (h <= 3L) stop("h must be > 3 (for Fisher z approximation).")
  if (q < 1L) stop("q must be >= 1.")
  
  ## Helper: ARMA(p,q) order selection by AIC
  select_arma_order <- function(resid, max_pq = 2) {
    best_aic <- Inf
    best_fit <- NULL
    best_order <- c(0, 0, 0)
    
    for (p in 0:max_pq) {
      for (q in 0:max_pq) {
        if (p == 0 && q == 0) next
        fit <- try(stats::arima(resid, order = c(p, 0, q)), silent = TRUE)
        if (inherits(fit, "try-error")) next
        aic <- fit$aic
        if (!is.na(aic) && aic < best_aic) {
          best_aic <- aic
          best_fit <- fit
          best_order <- c(p, 0, q)
        }
      }
    }
    
    if (is.null(best_fit)) {
      best_fit <- stats::arima(resid, order = c(0, 0, 0))
      best_order <- c(0, 0, 0)
    }
    list(order = best_order, fit = best_fit)
  }
  
  ## Step 1: Trend via moving average of width q
  print("Estimating shared trend via moving average...")
  ma_kernel <- rep(1 / q, q)
  ma1 <- as.numeric(stats::filter(x1, ma_kernel, sides = 2))
  ma2 <- as.numeric(stats::filter(x2, ma_kernel, sides = 2))
  trend_hat <- (ma1 + ma2) / 2
  
  ## Residuals for ARMA fitting and lambda
  print("Computing residuals and signal-to-noise ratios (lambda)...")
  idx_trend <- which(!is.na(trend_hat))
  if (length(idx_trend) < (5 * (max_pq + 1))) {
    warning("Very few non-NA points for ARMA fitting; results may be unstable.")
  }
  resid1 <- x1[idx_trend] - trend_hat[idx_trend]
  resid2 <- x2[idx_trend] - trend_hat[idx_trend]
  
  var_trend  <- stats::var(trend_hat[idx_trend], na.rm = TRUE)
  var_resid1 <- stats::var(resid1, na.rm = TRUE)
  var_resid2 <- stats::var(resid2, na.rm = TRUE)
  
  lambda1 <- if (is.finite(var_trend) && is.finite(var_resid1) && var_resid1 > 0)
    var_trend / var_resid1 else NA_real_
  lambda2 <- if (is.finite(var_trend) && is.finite(var_resid2) && var_resid2 > 0)
    var_trend / var_resid2 else NA_real_
  
  ## Step 2: Fit separate ARMA models to residuals of x1 and x2
  print("Fitting ARMA models for series 1...")
  arma_sel1 <- select_arma_order(resid1, max_pq = max_pq)
  arma_fit1 <- arma_sel1$fit
  arma_order1 <- arma_sel1$order
  arma_coef1 <- arma_fit1$coef
  arma_sigma2_1 <- arma_fit1$sigma2
  ar1_coefs <- arma_coef1[grepl("^ar", names(arma_coef1))]
  ma1_coefs <- arma_coef1[grepl("^ma", names(arma_coef1))]
  
  print("Fitting ARMA models for series 2...")
  arma_sel2 <- select_arma_order(resid2, max_pq = max_pq)
  arma_fit2 <- arma_sel2$fit
  arma_order2 <- arma_sel2$order
  arma_coef2 <- arma_fit2$coef
  arma_sigma2_2 <- arma_fit2$sigma2
  ar2_coefs <- arma_coef2[grepl("^ar", names(arma_coef2))]
  ma2_coefs <- arma_coef2[grepl("^ma", names(arma_coef2))]
  
  ## Step 3: Rolling correlations on MA-smoothed series
  print("Computing rolling correlations and pointwise p-values...")
  Y1 <- ma1
  Y2 <- ma2
  R <- rep(NA_real_, Tn)
  p_series <- rep(NA_real_, Tn)
  
  z0 <- atanh(rho0)  # Fisher z of null correlation
  
  for (t in seq_len(Tn)) {
    start <- t - h + 1L
    if (start < 1L) next
    idx <- start:t
    if (any(is.na(Y1[idx])) || any(is.na(Y2[idx]))) next
    y1w <- Y1[idx]
    y2w <- Y2[idx]
    if (stats::sd(y1w) == 0 || stats::sd(y2w) == 0) next
    r <- suppressWarnings(stats::cor(y1w, y2w))
    if (!is.finite(r) || abs(r) >= 1) next
    R[t] <- r
    
    ## Fisher z test for H0: rho = rho0 vs H1: rho < rho0
    z_r <- 0.5 * log((1 + r) / (1 - r))
    z_stat <- (z_r - z0) * sqrt(h - 3)
    p_series[t] <- stats::pnorm(z_stat, lower.tail = TRUE)
  }
  
  valid_idx <- which(!is.na(R) & !is.na(p_series))
  m <- length(valid_idx)
  if (m < 5L) stop("Too few valid rolling correlation estimates; check q and h.")
  
  ## Step 4: Effective BY threshold and observed state statistics
  print("Computing effective BY threshold and observed state statistics...")
  c_m <- sum(1 / seq_len(m))
  alpha_eff <- alpha / c_m
  p_thresh <- alpha_eff
  
  ## Invert Fisher z test to get x_eff (threshold on R_t)
  z_crit <- stats::qnorm(p_thresh, lower.tail = TRUE)
  z_x <- z0 + z_crit / sqrt(h - 3)
  x_eff <- tanh(z_x)
  
  I <- rep(NA_integer_, Tn)
  I[valid_idx] <- as.integer(R[valid_idx] < x_eff)
  
  I_obs <- I[valid_idx]
  n_I <- length(I_obs)
  
  frac_state_obs <- mean(I_obs)
  
  entry_rate_obs <- NA_real_
  if (n_I > 1L) {
    I_prev <- I_obs[-n_I]
    I_curr <- I_obs[-1L]
    steps01 <- sum(I_prev == 0L)
    entries <- sum(I_curr == 1L & I_prev == 0L)
    entry_rate_obs <- if (steps01 > 0) entries / steps01 else NA_real_
  }
  
  run_lengths <- integer(0)
  if (any(I_obs == 1L)) {
    rle_I <- rle(I_obs)
    run_lengths <- rle_I$lengths[rle_I$values == 1L]
  }
  mean_run_length_obs <- if (length(run_lengths) > 0L) mean(run_lengths) else NA_real_
  
  ## Step 5: HAC-based approximate n_eff and tau^2
  print("Estimating HAC-based approximate n_eff and tau^2...")
  ## Score-like process: W_t = Y1_t * Y2_t on valid_idx
  W <- Y1[valid_idx] * Y2[valid_idx]
  W <- W[is.finite(W)]
  Nw <- length(W)
  
  neff_approx <- NA_real_
  tau2_approx <- NA_real_
  
  if (Nw > 10L) {
    Wc <- W - mean(W)
    gamma0 <- stats::var(Wc)  # sample variance
    
    ## bandwidth K for HAC (Bartlett): simple rule K ~ N^(1/3), cap at 50
    K <- min(floor(Nw^(1/3)), 50L)
    if (K >= 1L && is.finite(gamma0) && gamma0 > 0) {
      gamma_k <- numeric(K)
      for (k in 1:K) {
        gamma_k[k] <- mean(Wc[1:(Nw - k)] * Wc[(1 + k):Nw])
      }
      w_k <- 1 - (1:K) / (K + 1)  # Bartlett taper
      sigma2_LR <- gamma0 + 2 * sum(w_k * gamma_k, na.rm = TRUE)
      if (is.finite(sigma2_LR) && sigma2_LR > 0) {
        tau2_approx <- sigma2_LR / gamma0
        neff_approx <- h / tau2_approx
      }
    }
  }
  
  print(paste("Approximate n_eff (HAC) ~", signif(neff_approx, 3),
              ", tau^2 ~", signif(tau2_approx, 3)))
  
  ## Step 6: Asymptotic expectations under the null (centered at rho0)
  print("Computing asymptotic expectations under the null...")
  
  entry_rate_th      <- NA_real_
  mean_run_length_th <- NA_real_
  frac_state_th      <- NA_real_
  
  if (is.finite(neff_approx) && neff_approx > 0) {
    R_valid <- R[valid_idx]
    if (length(R_valid) > 2L) {
      phi_R_hat <- suppressWarnings(
        stats::cor(R_valid[-length(R_valid)], R_valid[-1L],
                   use = "complete.obs")
      )
      if (!is.finite(phi_R_hat)) phi_R_hat <- 0
    } else {
      phi_R_hat <- 0
    }
    
    z_x_asym <- (x_eff - rho0) * sqrt(neff_approx)
    p_x <- stats::pnorm(z_x_asym)
    frac_state_th <- p_x
    
    if (p_x > 0 && p_x < 1 && requireNamespace("mvtnorm", quietly = TRUE)) {
      sigma_mat <- matrix(c(1, phi_R_hat,
                            phi_R_hat, 1), 2, 2)
      joint_prob <- as.numeric(
        mvtnorm::pmvnorm(
          upper = c(z_x_asym, z_x_asym),
          mean  = c(0, 0),
          sigma = sigma_mat
        )
      )
      pi11_th <- joint_prob / p_x
      pi01_th <- (p_x - joint_prob) / (1 - p_x)
      
      entry_rate_th <- pi01_th
      mean_run_length_th <- if (pi11_th < 1) 1 / (1 - pi11_th) else NA_real_
    }
  }
  
  print(paste("Asymptotic expected frac_state ~", signif(frac_state_th, 3),
              ", entry_rate ~", signif(entry_rate_th, 3),
              ", mean_run_length ~", signif(mean_run_length_th, 3)))
  
  print("Null model estimation complete.")
  
  null_model <- list(
    series1 = list(
      arma_order = arma_order1,
      ar = as.numeric(ar1_coefs),
      ma = as.numeric(ma1_coefs),
      sigma2 = arma_sigma2_1,
      lambda = lambda1
    ),
    series2 = list(
      arma_order = arma_order2,
      ar = as.numeric(ar2_coefs),
      ma = as.numeric(ma2_coefs),
      sigma2 = arma_sigma2_2,
      lambda = lambda2
    ),
    q_ma = q,
    h_corr = h,
    alpha = alpha,
    rho0 = rho0,
    neff_approx = neff_approx,
    tau2_approx = tau2_approx
  )
  
  observed <- list(
    entry_rate = entry_rate_obs,
    mean_run_length = mean_run_length_obs,
    frac_state = frac_state_obs
  )
  
  expected_asymptotic <- list(
    entry_rate      = entry_rate_th,
    mean_run_length = mean_run_length_th,
    frac_state      = frac_state_th
  )
  
  inputs <- list(
    q = q,
    h = h,
    alpha = alpha,
    rho0 = rho0,
    max_pq = max_pq,
    T = Tn
  )
  
  thresholds <- list(
    p_thresh = p_thresh,
    x_eff = x_eff,
    alpha_eff = alpha_eff
  )
  
  list(
    null_model = null_model,
    observed = observed,
    expected_asymptotic = expected_asymptotic,
    inputs = inputs,
    trend_hat = trend_hat,
    ma1 = ma1,
    ma2 = ma2,
    R = R,
    p_series = p_series,
    thresholds = thresholds,
    I = I,
    valid_idx = valid_idx
  )
}
