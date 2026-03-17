#' Parametric bootstrap for local correlation statistics
#'
#' Runs a parametric bootstrap under the fitted null model from
#' [lomad_fit()] to obtain null distributions and p-values for the
#' entry rate, mean run length, and fraction of time in the low-correlation
#' state.
#'
#' @param fit List returned by [lomad_fit()].
#' @param B Integer. Number of bootstrap replicates (default 500).
#' @param seed Integer or NULL. RNG seed for reproducibility.
#'
#' @return A named list:
#'   \describe{
#'     \item{null_model}{Null model from `fit`, augmented with
#'       bootstrap-estimated `neff` and `tau2`.}
#'     \item{observed}{Observed statistics from `fit` (including `n_entries`).}
#'     \item{expected}{Bootstrap means under the null.}
#'     \item{expected_asymptotic}{Asymptotic expectations from `fit`.}
#'     \item{p_values}{One-sided bootstrap p-values (P(stat >= observed)).}
#'     \item{inputs}{`n`, `max_pq`, `B`, `seed`.}
#'     \item{trend_hat, ma1, ma2, R, thresholds, I, valid_idx}{Passed through
#'       from `fit` for plotting and downstream use.}
#'   }
#'
#' @export
lomad_test <- function(fit,
                                 B    = 500,
                                 seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  null_model <- fit$null_model
  n          <- fit$inputs$n
  trend_hat  <- fit$trend_hat
  x_eff      <- fit$thresholds$x_eff
  z0         <- atanh(null_model$rho0)

  q  <- null_model$q
  h  <- null_model$h

  ma_kernel <- rep(1 / q, q)

  # --- Helper: simulate one series pair under the null ---

  sim_pair <- function() {
    e1 <- as.numeric(stats::arima.sim(
      model = list(ar = null_model$series1$ar, ma = null_model$series1$ma),
      n = n, sd = sqrt(null_model$series1$sigma2)
    ))
    e2 <- as.numeric(stats::arima.sim(
      model = list(ar = null_model$series2$ar, ma = null_model$series2$ma),
      n = n, sd = sqrt(null_model$series2$sigma2)
    ))
    list(x1 = trend_hat + e1, x2 = trend_hat + e2)
  }

  # --- Helper: state statistics from one simulated pair ---

  sim_stats <- function(x1, x2) {
    s1 <- as.numeric(stats::filter(x1, ma_kernel, sides = 2))
    s2 <- as.numeric(stats::filter(x2, ma_kernel, sides = 2))
    R  <- rep(NA_real_, n)

    for (t in seq_len(n)) {
      start <- t - h + 1L
      if (start < 1L) next
      w <- start:t
      y1w <- s1[w]; y2w <- s2[w]
      if (any(is.na(y1w)) || any(is.na(y2w))) next
      if (stats::sd(y1w) == 0 || stats::sd(y2w) == 0) next
      r <- suppressWarnings(stats::cor(y1w, y2w))
      if (is.finite(r)) R[t] <- r
    }

    v <- which(!is.na(R))
    if (length(v) < 5L)
      return(list(entry_rate = NA_real_, mean_run_length = NA_real_,
                  frac_state = NA_real_, var_R = NA_real_))

    I_v        <- as.integer(R[v] < x_eff)
    m          <- length(I_v)
    frac_state <- mean(I_v)

    steps01    <- sum(I_v[-m] == 0L)
    entry_rate <- if (steps01 > 0L)
      sum(I_v[-1L] == 1L & I_v[-m] == 0L) / steps01 else NA_real_

    runs           <- rle(I_v)
    run_lengths    <- runs$lengths[runs$values == 1L]
    n_entries_sim  <- length(run_lengths)
    mean_run_length <- if (n_entries_sim > 0L) mean(run_lengths) else NA_real_

    list(entry_rate      = entry_rate,
         mean_run_length = mean_run_length,
         frac_state      = frac_state,
         n_entries       = as.numeric(n_entries_sim),
         var_R           = stats::var(R[v], na.rm = TRUE))
  }

  # --- Bootstrap loop ---

  boot <- lapply(seq_len(B), function(b) {
    pair <- sim_pair()
    sim_stats(pair$x1, pair$x2)
  })

  entry_v  <- Filter(is.finite, vapply(boot, `[[`, numeric(1), "entry_rate"))
  run_v    <- Filter(is.finite, vapply(boot, `[[`, numeric(1), "mean_run_length"))
  frac_v   <- Filter(is.finite, vapply(boot, `[[`, numeric(1), "frac_state"))
  nent_v   <- Filter(is.finite, vapply(boot, `[[`, numeric(1), "n_entries"))
  var_R_v  <- Filter(is.finite, vapply(boot, `[[`, numeric(1), "var_R"))

  # --- Summaries ---

  obs <- fit$observed

  expected <- list(
    entry_rate      = if (length(entry_v)) mean(entry_v) else NA_real_,
    mean_run_length = if (length(run_v))   mean(run_v)   else NA_real_,
    frac_state      = if (length(frac_v))  mean(frac_v)  else NA_real_,
    n_entries       = if (length(nent_v))  mean(nent_v)  else NA_real_
  )

  p_values <- list(
    entry_rate      = if (!is.na(obs$entry_rate)      && length(entry_v))
      mean(entry_v >= obs$entry_rate)      else NA_real_,
    mean_run_length = if (!is.na(obs$mean_run_length) && length(run_v))
      mean(run_v   >= obs$mean_run_length) else NA_real_,
    frac_state      = if (!is.na(obs$frac_state)      && length(frac_v))
      mean(frac_v  >= obs$frac_state)      else NA_real_,
    n_entries       = if (!is.na(obs$n_entries)       && length(nent_v))
      mean(nent_v  >= obs$n_entries)       else NA_real_
  )

  mean_var_R  <- if (length(var_R_v)) mean(var_R_v) else NA_real_
  neff_boot   <- if (is.finite(mean_var_R) && mean_var_R > 0) 1 / mean_var_R else NA_real_
  tau2_boot   <- if (is.finite(neff_boot)) h / neff_boot else NA_real_

  message(sprintf(
    "B = %d | p: entry_rate = %s, run_length = %s, frac_state = %s, n_entries = %s",
    B,
    signif(p_values$entry_rate, 3),
    signif(p_values$mean_run_length, 3),
    signif(p_values$frac_state, 3),
    signif(p_values$n_entries, 3)
  ))

  null_model$neff <- neff_boot
  null_model$tau2 <- tau2_boot
  null_model$seed <- seed

  list(
    null_model          = null_model,
    observed            = obs,
    expected            = expected,
    expected_asymptotic = fit$expected_asymptotic,
    p_values            = p_values,
    inputs              = list(n = fit$inputs$n, max_pq = fit$inputs$max_pq,
                               B = B, seed = seed),
    trend_hat  = fit$trend_hat,
    ma1        = fit$ma1,
    ma2        = fit$ma2,
    R          = fit$R,
    thresholds = fit$thresholds,
    I          = fit$I,
    valid_idx  = fit$valid_idx
  )
}
