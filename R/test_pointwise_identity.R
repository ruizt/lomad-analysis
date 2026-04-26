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
#'   \item Estimate the common trend as \eqn{\hat\nu_t = (Y_{1t}^{(q)} +
#'     Y_{2t}^{(q)})/2} and compute noise residuals \eqn{\hat Z_{it} = X_{it}
#'     - \hat\nu_t}.
#'   \item Fit ARMA models to \eqn{\hat Z_{1t}} and \eqn{\hat Z_{2t}} via
#'     AIC-based order selection up to order \code{max_pq}.
#'   \item Compute the lag-zero autocovariance of the MA-smoothed noise:
#'     \deqn{\mathrm{Var}(\eta_{it}^{(q)}) =
#'       \frac{1}{q^2}\left[q\gamma_i(0) +
#'       2\sum_{l=1}^{q-1}(q-l)\gamma_i(l)\right]}
#'     where \eqn{\gamma_i(l)} is the process autocovariance at lag \eqn{l},
#'     derived from the ARMA parameters.
#'   \item Form the pointwise test statistic
#'     \eqn{z_t = D_t / \mathrm{SE}(D_t)} where
#'     \eqn{D_t = Y_{1t}^{(q)} - Y_{2t}^{(q)}} and
#'     \eqn{\mathrm{SE}(D_t) = \sqrt{\mathrm{Var}(\eta_{1t}^{(q)}) +
#'     \mathrm{Var}(\eta_{2t}^{(q)})}}.
#'   \item Apply the Benjamini-Yekutieli correction across time points and
#'     flag \eqn{I_t = 1} where the adjusted p-value is below \code{alpha}.
#' }
#'
#' @param x1 Numeric vector. First observed series.
#' @param x2 Numeric vector. Second observed series.
#' @param q Integer. Moving-average window width for trend smoothing.
#' @param max_pq Integer. Maximum AR and MA order considered in AIC selection
#'   (default 3).
#' @param alpha Numeric. Significance level for Benjamini-Yekutieli correction
#'   (default 0.05).
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
#'   \item{fit1}{List. ARMA fit for series 1 noise residuals.}
#'   \item{fit2}{List. ARMA fit for series 2 noise residuals.}
#' }
#'
#' @export
test_pointwise_identity <- function(x1, x2, q, max_pq = 3L, alpha = 0.05) {

  n <- length(x1)
  q <- as.integer(q)

  # --- Step 1: centered MA trend estimate ---

  ma_kernel <- rep(1 / q, q)
  ma1       <- as.numeric(stats::filter(x1, ma_kernel, sides = 2))
  ma2       <- as.numeric(stats::filter(x2, ma_kernel, sides = 2))
  trend_hat <- (ma1 + ma2) / 2

  # --- Step 2: noise residuals ---

  valid_t <- which(!is.na(trend_hat))
  resid1  <- x1[valid_t] - trend_hat[valid_t]
  resid2  <- x2[valid_t] - trend_hat[valid_t]

  # --- Step 3: AIC-based ARMA selection ---

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

  fit1 <- select_arma(resid1)
  fit2 <- select_arma(resid2)

  # --- Step 4: lag-zero autocovariance of MA_q-smoothed noise ---
  # Var(MA_q(Z)) = (1/q^2) * [q*gamma(0) + 2*sum_{l=1}^{q-1}(q-l)*gamma(l)]

  smooth_noise_var <- function(fit, q) {
    ar     <- fit$ar
    ma     <- fit$ma
    sigma2 <- fit$sigma2
    max_lag <- q - 1L

    if (length(ar) == 0L && length(ma) == 0L) {
      # white noise: Var(MA_q(Z)) = sigma2 / q
      return(sigma2 / q)
    }

    # process autocovariance gamma(l) = sigma2 * sum(psi_k * psi_{k+l})
    psi    <- c(1, stats::ARMAtoMA(ar = ar, ma = ma,
                                   lag.max = max(200L, max_lag)))
    gamma0 <- sigma2 * sum(psi^2)
    acf    <- stats::ARMAacf(ar = ar, ma = ma, lag.max = max_lag)
    gamma  <- gamma0 * as.numeric(acf)

    weights <- q - seq_len(max_lag)
    (q * gamma[1L] + 2 * sum(weights * gamma[-1L])) / q^2
  }

  v1   <- smooth_noise_var(fit1, q)
  v2   <- smooth_noise_var(fit2, q)
  se_D <- sqrt(v1 + v2)

  # --- Step 5: pointwise z-test ---

  D     <- ma1 - ma2
  z_t   <- D / se_D
  p_raw <- 2 * stats::pnorm(-abs(z_t))

  # --- Step 6: BY correction over valid time points ---

  valid <- !is.na(p_raw)
  p_adj <- rep(NA_real_, n)
  p_adj[valid] <- stats::p.adjust(p_raw[valid], method = "BY")

  I <- rep(NA_integer_, n)
  I[valid] <- as.integer(p_adj[valid] < alpha)

  list(
    I     = I,
    p_raw = p_raw,
    p_adj = p_adj,
    D     = D,
    se_D  = se_D,
    fit1  = fit1,
    fit2  = fit2
  )
}
