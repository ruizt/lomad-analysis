#' Naive pointwise test of exact trend identity
#'
#' Tests $H_0: \nu_{1t} = \nu_{2t}$ at each time point using the pointwise
#' difference in smoothed moving averages. The standard error is derived from
#' the lag-zero autocovariance of the MA-smoothed noise process, computed
#' analytically from ARMA parameter estimates.
#'
#' This test is intended as a naive comparator for the local-similarity
#' method: because it tests exact identity rather than similarity, it will
#' reject whenever the trends differ by any nonzero amount, regardless of
#' whether that difference is materially meaningful.
#'
#' @details
#' The procedure is:
#' \enumerate{
#'   \item Compute centered moving averages \eqn{Y_{it}^{(q)}} for each series.
#'   \item Pilot-detrend the difference \eqn{\Delta_t = X_{1t} - X_{2t}} using
#'     a wider MA of width \code{q_se} to remove the trend component
#'     \eqn{\nu_{1t} - \nu_{2t}}, yielding residuals \eqn{\tilde\Delta_t \approx
#'     Z_{1t} - Z_{2t}}. This prevents a non-zero trend difference from
#'     inflating \eqn{\hat\sigma^2} and making the test conservative under
#'     \eqn{H_1}.
#'   \item Fit a single ARMA model to \eqn{\tilde\Delta_t} via AIC-based order
#'     selection up to order \code{max_pq}.
#'   \item Compute the lag-zero autocovariance of the MA-smoothed noise:
#'     \deqn{\mathrm{Var}(D_t) =
#'       \frac{1}{q^2}\left[q\gamma_\Delta(0) +
#'       2\sum_{l=1}^{q-1}(q-l)\gamma_\Delta(l)\right]}
#'     where \eqn{\gamma_\Delta(l)} is the autocovariance of the fitted
#'     noise model at lag \eqn{l}.
#'   \item Form the pointwise test statistic
#'     \eqn{z_t = D_t / \mathrm{SE}(D_t)} where
#'     \eqn{D_t = Y_{1t}^{(q)} - Y_{2t}^{(q)}} and
#'     \eqn{\mathrm{SE}(D_t) = \sqrt{\mathrm{Var}(D_t)}}.
#'   \item Apply a Benjamini-Yekutieli correction using the effective number of
#'     independent tests \eqn{m_\text{eff} = \lfloor n_\text{valid} / q \rfloor}.
#'     Adjacent \eqn{D_t} values share \eqn{q - 1} of \eqn{q} observations, so
#'     the \eqn{n_\text{valid}} test statistics contain only \eqn{\approx
#'     n_\text{valid}/q} independent pieces of information. Using
#'     \eqn{m_\text{eff}} in place of \eqn{n_\text{valid}} in the BY constant
#'     corrects for this without discarding any test statistics.
#'     Flag \eqn{I_t = 1} where the adjusted p-value is below \code{alpha}.
#' }
#'
#' @param x1 Numeric vector. First observed series.
#' @param x2 Numeric vector. Second observed series.
#' @param q Integer. Moving-average window width for trend smoothing.
#' @param max_pq Integer. Maximum AR and MA order considered in AIC selection
#'   (default 3).
#' @param alpha Numeric. Significance level for Benjamini-Yekutieli correction
#'   (default 0.05).
#' @param q_se Integer. Pilot MA window width used to remove the trend component
#'   from \eqn{X_{1t} - X_{2t}} before ARMA fitting. Must satisfy
#'   \code{q_se > q}. Defaults to \code{10L * q}.
#'
#' @return A list with:
#' \describe{
#'   \item{I}{Integer vector. Binary detection indicator: 1 where the adjusted
#'     p-value is below \code{alpha}, 0 elsewhere, NA at edge time points.}
#'   \item{p_raw}{Numeric vector. Raw two-sided pointwise p-values.}
#'   \item{p_adj}{Numeric vector. BY-adjusted p-values.}
#'   \item{D}{Numeric vector. Pointwise difference \eqn{D_t = Y_{1t}^{(q)} -
#'     Y_{2t}^{(q)}}.}
#'   \item{se_D}{Numeric scalar. Standard error of \eqn{D_t} (constant across
#'     \eqn{t} under stationarity).}
#'   \item{m_eff}{Integer. Effective number of independent tests used in BY correction.}
#'   \item{fit_diff}{List. ARMA fit for the pilot-detrended difference series.}
#' }
#'
#' @export
test_pointwise_identity <- function(x1, x2, q, max_pq = 3L, alpha = 0.05,
                                    q_se = 10L * q) {

  n    <- length(x1)
  q    <- as.integer(q)
  q_se <- as.integer(q_se)

  # --- Step 1: centered MA trend estimate (test window) ---

  ma_kernel <- rep(1 / q, q)
  ma1       <- as.numeric(stats::filter(x1, ma_kernel, sides = 2))
  ma2       <- as.numeric(stats::filter(x2, ma_kernel, sides = 2))

  # --- Step 2: pilot-detrended difference series ---
  # Apply a wider MA (q_se >> q) to x1 - x2 to remove the trend component
  # nu1 - nu2, leaving approximately Z1 - Z2 for ARMA fitting. Without this
  # step, a non-zero trend difference under H1 inflates sigma2 and makes
  # se_D too large, causing the test to become more conservative as d grows.

  diff_full    <- x1 - x2
  pilot_kernel <- rep(1 / q_se, q_se)
  d_trend_est  <- as.numeric(stats::filter(diff_full, pilot_kernel, sides = 2))
  d_detrended  <- diff_full - d_trend_est
  valid_arma   <- which(!is.na(d_detrended))

  valid_t      <- which(!is.na(ma1) & !is.na(ma2))

  # --- Step 3: AIC-based ARMA selection on detrended difference ---

  select_arma <- function(resid) {
    best_aic <- Inf
    best_fit <- NULL
    for (p in 0:max_pq) {
      for (qi in 0:max_pq) {
        if (p == 0L && qi == 0L) next
        converged <- TRUE
        fit <- withCallingHandlers(
          try(stats::arima(resid, order = c(p, 0L, qi),
                           optim.control = list(maxit = 500L)),
              silent = TRUE),
          warning = function(w) {
            if (grepl("convergence problem", conditionMessage(w))) {
              converged <<- FALSE
              invokeRestart("muffleWarning")
            }
          }
        )
        if (inherits(fit, "try-error") || is.na(fit$aic) || !converged) next
        if (fit$aic < best_aic) {
          best_aic <- fit$aic
          best_fit <- fit
        }
      }
    }
    if (is.null(best_fit))
      best_fit <- stats::arima(resid, order = c(0L, 0L, 0L))
    coef <- best_fit$coef
    list(ar     = as.numeric(coef[grepl("^ar", names(coef))]),
         ma     = as.numeric(coef[grepl("^ma", names(coef))]),
         sigma2 = best_fit$sigma2)
  }

  fit_diff <- select_arma(d_detrended[valid_arma])

  # --- Step 4: Var(MA_q(Delta)) = se_D^2 ---

  smooth_noise_var <- function(fit, q) {
    ar     <- fit$ar
    ma     <- fit$ma
    sigma2 <- fit$sigma2
    max_lag <- q - 1L

    if (length(ar) == 0L && length(ma) == 0L) {
      return(sigma2 / q)
    }

    psi    <- c(1, stats::ARMAtoMA(ar = ar, ma = ma,
                                   lag.max = max(200L, max_lag)))
    gamma0 <- sigma2 * sum(psi^2)
    acf    <- stats::ARMAacf(ar = ar, ma = ma, lag.max = max_lag)
    gamma  <- gamma0 * as.numeric(acf)

    weights <- q - seq_len(max_lag)
    (q * gamma[1L] + 2 * sum(weights * gamma[-1L])) / q^2
  }

  se_D <- sqrt(smooth_noise_var(fit_diff, q))

  # --- Step 5: pointwise z-test ---

  D     <- ma1 - ma2
  z_t   <- D / se_D
  p_raw <- 2 * stats::pnorm(-abs(z_t))

  # --- Step 6: BY correction using effective number of independent tests ---

  valid   <- !is.na(p_raw)
  n_valid <- sum(valid)
  m_eff   <- max(1L, floor(n_valid / q))
  H_meff  <- sum(1 / seq_len(m_eff))
  const   <- m_eff * H_meff

  p_sub   <- p_raw[valid]
  ord     <- order(p_sub)
  raw_adj <- p_sub[ord] * const / seq_len(n_valid)
  step_up <- rev(cummin(rev(raw_adj)))

  p_adj <- rep(NA_real_, n)
  p_adj[valid] <- pmin(1, step_up)[order(ord)]

  I <- rep(NA_integer_, n)
  I[valid] <- as.integer(p_adj[valid] < alpha)

  list(
    I        = I,
    p_raw    = p_raw,
    p_adj    = p_adj,
    D        = D,
    se_D     = se_D,
    m_eff    = m_eff,
    fit_diff = fit_diff
  )
}
