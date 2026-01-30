## bootstrap_local_corr.R
##
## Parametric bootstrap for locally low correlation statistics, given
## a fitted null model from fit_null_model_local_corr().
##
## INPUTS:
##   fit_obj   : list returned by fit_null_model_local_corr()
##   B         : number of bootstrap replicates (default 500)
##   boot_seed : (optional) integer RNG seed for reproducible bootstrap
##
## OUTPUT:
##   A list with:
##     $null_model : original null_model augmented with:
##                     $neff : bootstrap-estimated effective sample size
##                     $tau2 : corresponding variance inflation factor
##                     $boot_seed
##     $expected   : bootstrap means under the null:
##                   $entry_rate, $mean_run_length, $frac_state
##     $observed   : observed statistics (copied from fit_obj$observed)
##     $expected_asymptotic :
##                   asymptotic expectations (copied from fit_obj)
##     $p_values   : one-sided bootstrap p-values for each statistic
##     $inputs     : original inputs plus B and boot_seed
##     $trend_hat, $ma1, $ma2, $R, $p_series, $thresholds :
##                   copied from fit_obj for convenience

bootstrap_local_corr <- function(fit_obj,
                                 B = 500,
                                 boot_seed = NULL) {
  
  print("Starting parametric bootstrap...")
  
  if (!is.null(boot_seed)) {
    set.seed(boot_seed)
  }
  
  null_model <- fit_obj$null_model
  inputs <- fit_obj$inputs
  thresholds <- fit_obj$thresholds
  trend_hat <- fit_obj$trend_hat
  
  q <- null_model$q_ma
  h <- null_model$h_corr
  alpha <- null_model$alpha
  
  Tn <- inputs$T
  stopifnot(length(trend_hat) == Tn)
  
  series1 <- null_model$series1
  series2 <- null_model$series2
  
  ar1_coefs <- series1$ar
  ma1_coefs <- series1$ma
  sigma2_1  <- series1$sigma2
  
  ar2_coefs <- series2$ar
  ma2_coefs <- series2$ma
  sigma2_2  <- series2$sigma2
  
  x_eff <- thresholds$x_eff
  p_thresh <- thresholds$p_thresh
  
  print(paste("Using q =", q, "and h =", h, "with alpha =", alpha))
  print(paste("Effective BY p-threshold:", signif(p_thresh, 3),
              " --> R_t threshold x_eff =", signif(x_eff, 3)))
  
  ## Moving-average kernel for within-bootstrap trend smoothing
  ma_kernel <- rep(1 / q, q)
  df <- h - 2L
  
  ## Helper: simulate one pair under the null
  simulate_pair <- function() {
    model_list1 <- list(ar = as.numeric(ar1_coefs),
                        ma = as.numeric(ma1_coefs))
    model_list2 <- list(ar = as.numeric(ar2_coefs),
                        ma = as.numeric(ma2_coefs))
    
    e1 <- as.numeric(stats::arima.sim(model = model_list1,
                                      n = Tn,
                                      sd = sqrt(sigma2_1)))
    e2 <- as.numeric(stats::arima.sim(model = model_list2,
                                      n = Tn,
                                      sd = sqrt(sigma2_2)))
    
    list(x1 = trend_hat + e1,
         x2 = trend_hat + e2)
  }
  
  ## Helper: compute stats for a simulated pair
  compute_stats_for_pair <- function(x1_star, x2_star) {
    ma1_star <- as.numeric(stats::filter(x1_star, ma_kernel, sides = 2))
    ma2_star <- as.numeric(stats::filter(x2_star, ma_kernel, sides = 2))
    Y1_star <- ma1_star
    Y2_star <- ma2_star
    
    R_star <- rep(NA_real_, Tn)
    p_star <- rep(NA_real_, Tn)
    
    for (t in seq_len(Tn)) {
      start <- t - h + 1L
      if (start < 1L) next
      idx <- start:t
      if (any(is.na(Y1_star[idx])) || any(is.na(Y2_star[idx]))) next
      y1w <- Y1_star[idx]
      y2w <- Y2_star[idx]
      if (stats::sd(y1w) == 0 || stats::sd(y2w) == 0) next
      r <- suppressWarnings(stats::cor(y1w, y2w))
      if (!is.finite(r)) next
      R_star[t] <- r
      t_stat <- r * sqrt(df / (1 - r^2))
      p_star[t] <- stats::pt(t_stat, df = df, lower.tail = TRUE)
    }
    
    valid_star <- which(!is.na(R_star) & !is.na(p_star))
    if (length(valid_star) < 5L) {
      return(list(entry_rate = NA_real_,
                  mean_run_length = NA_real_,
                  frac_state = NA_real_,
                  varR = NA_real_))
    }
    
    varR <- stats::var(R_star[valid_star], na.rm = TRUE)
    
    I_star <- rep(NA_integer_, Tn)
    I_star[valid_star] <- as.integer(R_star[valid_star] < x_eff)
    I_s <- I_star[valid_star]
    n_s <- length(I_s)
    
    frac_state <- mean(I_s)
    
    entry_rate <- NA_real_
    mean_run_length <- NA_real_
    
    if (n_s > 1L) {
      I_prev_s <- I_s[-n_s]
      I_curr_s <- I_s[-1L]
      steps01_s <- sum(I_prev_s == 0L)
      entries_s <- sum(I_curr_s == 1L & I_prev_s == 0L)
      entry_rate <- if (steps01_s > 0) entries_s / steps01_s else NA_real_
    }
    
    if (any(I_s == 1L)) {
      rle_s <- rle(I_s)
      run_lengths_s <- rle_s$lengths[rle_s$values == 1L]
      mean_run_length <- if (length(run_lengths_s) > 0L) mean(run_lengths_s) else NA_real_
    }
    
    list(entry_rate = entry_rate,
         mean_run_length = mean_run_length,
         frac_state = frac_state,
         varR = varR)
  }
  
  entry_boot  <- numeric(B)
  runlen_boot <- numeric(B)
  frac_boot   <- numeric(B)
  varR_boot   <- numeric(B)
  
  cat("Running parametric bootstrap with", B, "replicates:\n")
  for (b in seq_len(B)) {
    if (b == 1) {
      cat("  Iteration:")
    }
    if (b %% 10 == 1) {
      cat("\n   ")
    }
    cat(b, "")
    
    pair <- simulate_pair()
    stats_b <- compute_stats_for_pair(pair$x1, pair$x2)
    entry_boot[b]  <- stats_b$entry_rate
    runlen_boot[b] <- stats_b$mean_run_length
    frac_boot[b]   <- stats_b$frac_state
    varR_boot[b]   <- stats_b$varR
  }
  cat("\nBootstrap complete.\n")
  
  entry_boot  <- entry_boot[is.finite(entry_boot)]
  runlen_boot <- runlen_boot[is.finite(runlen_boot)]
  frac_boot   <- frac_boot[is.finite(frac_boot)]
  varR_boot   <- varR_boot[is.finite(varR_boot)]
  
  ## Expected values and p-values
  observed <- fit_obj$observed
  
  expected_entry <- if (length(entry_boot) > 0L) mean(entry_boot) else NA_real_
  expected_run   <- if (length(runlen_boot) > 0L) mean(runlen_boot) else NA_real_
  expected_frac  <- if (length(frac_boot) > 0L) mean(frac_boot) else NA_real_
  
  p_entry <- if (!is.na(observed$entry_rate) && length(entry_boot) > 0L) {
    mean(entry_boot >= observed$entry_rate)
  } else NA_real_
  
  p_run <- if (!is.na(observed$mean_run_length) && length(runlen_boot) > 0L) {
    mean(runlen_boot >= observed$mean_run_length)
  } else NA_real_
  
  p_frac <- if (!is.na(observed$frac_state) && length(frac_boot) > 0L) {
    mean(frac_boot >= observed$frac_state)
  } else NA_real_
  
  ## Estimate Var(R), n_eff, tau^2
  mean_varR <- if (length(varR_boot) > 0L) mean(varR_boot) else NA_real_
  neff_est  <- if (is.finite(mean_varR) && mean_varR > 0) 1 / mean_varR else NA_real_
  tau2_est  <- if (is.finite(neff_est)) h / neff_est else NA_real_
  
  print(paste("Estimated n_eff (bootstrap) ~", signif(neff_est, 3),
              "and tau^2 ~", signif(tau2_est, 3)))
  
  ## Augment null_model with neff, tau2, and boot_seed
  null_model$neff <- neff_est
  null_model$tau2 <- tau2_est
  null_model$boot_seed <- boot_seed
  
  expected <- list(
    entry_rate = expected_entry,
    mean_run_length = expected_run,
    frac_state = expected_frac
  )
  
  p_values <- list(
    entry_rate = p_entry,
    mean_run_length = p_run,
    frac_state = p_frac
  )
  
  inputs_out <- c(
    inputs,
    list(
      B = B,
      boot_seed = boot_seed
    )
  )
  
  print("Bootstrap procedure complete.")
  
  list(
    null_model = null_model,
    expected = expected,
    observed = observed,
    expected_asymptotic = fit_obj$expected_asymptotic,
    p_values = p_values,
    inputs = inputs_out,
    trend_hat = fit_obj$trend_hat,
    ma1 = fit_obj$ma1,
    ma2 = fit_obj$ma2,
    R = fit_obj$R,
    p_series = fit_obj$p_series,
    thresholds = fit_obj$thresholds
  )
}