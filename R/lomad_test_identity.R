#' Naive pointwise test of exact trend identity
#'
#' Tests \eqn{H_0: \nu_{1t} = \nu_{2t}} at each time point using the
#' pointwise difference in smoothed moving averages. Included as a comparison
#' method: because it tests exact identity rather than similarity, it rejects
#' whenever the trends differ by any nonzero amount.
#'
#' @param x1 Numeric vector. First observed series.
#' @param x2 Numeric vector. Second observed series.
#' @param q Integer. Moving-average window width for trend smoothing.
#' @param max_pq Integer. Maximum AR and MA order for AIC selection (default 3).
#' @param alpha Numeric. Significance level for BY correction (default 0.05).
#' @param q_se Integer. Pilot MA window for detrending the difference series
#'   before ARMA fitting. Must satisfy `q_se > q`. Defaults to `10L * q`.
#' @param noise_override Optional list with elements `ar`, `ma` (optional),
#'   `sigma2`. If supplied, these parameters are used to compute the variance
#'   of the smoothed difference directly, bypassing ARMA estimation. The
#'   parameters describe one noise series; the difference of two independent
#'   copies has twice the variance, which is handled internally. Useful for
#'   oracle experiments where the true noise process is known.
#'
#' @return A list with:
#'   \describe{
#'     \item{I}{Integer vector. 1 where adjusted p-value < alpha, 0 elsewhere,
#'       NA at edges.}
#'     \item{p_raw}{Numeric vector. Raw two-sided p-values.}
#'     \item{p_adj}{Numeric vector. BY-adjusted p-values.}
#'     \item{D}{Numeric vector. Pointwise difference of smoothed series.}
#'     \item{se_D}{Numeric scalar. Standard error of D_t.}
#'     \item{m_eff}{Integer. Effective number of independent tests.}
#'     \item{fit_diff}{List. ARMA fit for the pilot-detrended difference
#'       (NULL when `noise_override` is used).}
#'   }
#'
#' @export
lomad_test_identity <- function(x1, x2, q, max_pq = 3L, alpha = 0.05,
                                q_se = 10L * q, noise_override = NULL) {
  n    <- length(x1)
  q    <- as.integer(q)
  q_se <- as.integer(q_se)

  # Centered MA trend estimate
  ma_kernel <- rep(1 / q, q)
  ma1       <- as.numeric(stats::filter(x1, ma_kernel, sides = 2))
  ma2       <- as.numeric(stats::filter(x2, ma_kernel, sides = 2))

  valid_t <- which(!is.na(ma1) & !is.na(ma2))

  if (!is.null(noise_override)) {
    # Oracle: build a fit_diff for the difference series (2x variance)
    fit_diff <- list(
      ar     = noise_override$ar,
      ma     = if (!is.null(noise_override$ma)) noise_override$ma else numeric(0),
      sigma2 = 2 * noise_override$sigma2
    )
  } else {
    # Estimate: pilot-detrend then select ARMA
    diff_full    <- x1 - x2
    pilot_kernel <- rep(1 / q_se, q_se)
    d_trend_est  <- as.numeric(stats::filter(diff_full, pilot_kernel, sides = 2))
    d_detrended  <- diff_full - d_trend_est
    valid_arma   <- which(!is.na(d_detrended))
    fit_diff     <- .select_arma(d_detrended[valid_arma], max_pq)
  }

  # Var(MA_q(Delta))
  se_D <- sqrt(.smooth_noise_var(fit_diff, q))

  # Pointwise z-test
  D     <- ma1 - ma2
  z_t   <- D / se_D
  p_raw <- 2 * stats::pnorm(-abs(z_t))

  # BY correction using effective number of independent tests
  n_valid <- length(valid_t)
  m_eff   <- max(1L, floor(n_valid / q))

  alpha_eff <- alpha / sum(1 / seq_len(m_eff))
  p_adj     <- rep(NA_real_, n)
  I         <- rep(NA_integer_, n)

  p_adj[valid_t] <- pmin(p_raw[valid_t] * sum(1 / seq_len(m_eff)), 1)
  I[valid_t]     <- as.integer(p_raw[valid_t] <= alpha_eff)

  n_rejected <- sum(I, na.rm = TRUE)
  message(sprintf("Identity test: %d / %d rejections at BY-FDR = %.2f",
                  n_rejected, n_valid, alpha))

  list(
    I        = I,
    p_raw    = p_raw,
    p_adj    = p_adj,
    D        = D,
    se_D     = se_D,
    m_eff    = m_eff,
    fit_diff = if (is.null(noise_override)) fit_diff else NULL
  )
}
