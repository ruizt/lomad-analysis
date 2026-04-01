#' Analytic CLT test for local correlation decoupling
#'
#' Closed-form inference using a 2-state Markov chain CLT approximation derived
#' from [lomad_fit()]. No simulation is required. Currently provides an
#' analytic p-value only for `frac_state`; p-values for the other statistics
#' are `NA`.
#'
#' @param fit List returned by [lomad_fit()] or [lomad_fit_blocks()].
#' @param T_eff Integer or NULL. Effective chain length for the CLT. If `NULL`,
#'   defaults to `length(fit$valid_idx)` for single-series fits, or
#'   `median(fit$inputs$block_lengths)` for multi-block fits.
#'
#' @return A named list with the same structure as [lomad_test_boot()]:
#'   \describe{
#'     \item{null_model}{Null model from `fit`, augmented with Markov
#'       transition parameters (`transitions`) and `T_eff`.}
#'     \item{observed}{Observed statistics from `fit`.}
#'     \item{expected}{Analytic Markov chain expectations (same as
#'       `expected_asymptotic`).}
#'     \item{expected_asymptotic}{Asymptotic expectations from `fit`.}
#'     \item{p_values}{One-sided p-values: `frac_state` via CLT z-test;
#'       `entry_rate`, `mean_run_length`, `n_entries` are `NA` (no analytic
#'       form currently available).}
#'     \item{inputs}{`n`, `max_pq`, `T_eff`.}
#'     \item{trend_hat, ma1, ma2, R, thresholds, I, valid_idx}{Passed through
#'       from `fit` for plotting and downstream use.}
#'   }
#'
#' @export
lomad_test_analytic <- function(fit,
                                 T_eff = NULL) {

  observed  <- fit$observed
  expectedA <- fit$expected_asymptotic

  # --- Markov transition parameters from asymptotic expectations ---

  pi1_th  <- expectedA$frac_state
  L_th    <- expectedA$mean_run_length
  pi11_th <- if (!is.na(L_th) && L_th > 1) 1 - 1 / L_th else NA_real_
  pi0_th  <- if (!is.na(pi1_th)) 1 - pi1_th else NA_real_
  pi01_th <- if (!is.na(pi1_th) && !is.na(pi11_th) &&
                 pi1_th > 0 && pi1_th < 1)
    pi1_th * (1 - pi11_th) / (1 - pi1_th) else NA_real_
  pi00_th <- if (!is.na(pi01_th)) 1 - pi01_th else NA_real_

  if (is.null(T_eff)) {
    T_eff <- if (!is.null(fit$inputs$block_lengths)) {
      as.integer(stats::median(fit$inputs$block_lengths))
    } else {
      length(fit$valid_idx)
    }
  }

  # --- CLT for frac_state under 2-state Markov chain ---
  # Var(pi1_hat) ~ Var(I_t) * (1 + theta) / (1 - theta) / T_eff
  # where theta = pi11 + pi00 - 1 is the lag-1 autocorrelation of I_t.

  p_frac <- NA_real_

  if (is.finite(pi1_th) && pi1_th > 0 && pi1_th < 1 &&
      is.finite(pi11_th) && is.finite(pi00_th) && T_eff > 0) {

    theta <- pi11_th + pi00_th - 1

    if (is.finite(theta) && abs(theta) < 1) {
      var_CLT <- pi1_th * (1 - pi1_th) * (1 + theta) / (1 - theta) / T_eff
      if (var_CLT > 0)
        p_frac <- 1 - stats::pnorm((observed$frac_state - pi1_th) / sqrt(var_CLT))
    }
  }

  message(sprintf("p (frac_state, analytic CLT) = %s", signif(p_frac, 3)))

  null_model             <- fit$null_model
  null_model$transitions <- list(pi1  = pi1_th,  pi0  = pi0_th,
                                  pi11 = pi11_th, pi01 = pi01_th, pi00 = pi00_th)
  null_model$T_eff       <- T_eff

  list(
    null_model          = null_model,
    observed            = observed,
    expected            = expectedA,
    expected_asymptotic = expectedA,
    p_values            = list(entry_rate      = NA_real_,
                               mean_run_length = NA_real_,
                               frac_state      = p_frac,
                               n_entries       = NA_real_),
    inputs              = list(n = fit$inputs$n, max_pq = fit$inputs$max_pq,
                               T_eff = T_eff),
    trend_hat  = fit$trend_hat,
    ma1        = fit$ma1,
    ma2        = fit$ma2,
    R          = fit$R,
    thresholds = fit$thresholds,
    I          = fit$I,
    valid_idx  = fit$valid_idx
  )
}
