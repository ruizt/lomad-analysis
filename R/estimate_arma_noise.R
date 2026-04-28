#' Estimate ARMA noise models for two de-trended series
#'
#' Subtracts the estimated shared trend from each observed series and fits
#' an ARMA(p, q) model to each set of residuals via AIC-based order selection.
#' The ARMA fit is global (using all valid, non-NA residuals), consistent with
#' the locally-stationary assumption that the noise covariance structure does
#' not vary with time.
#'
#' @param y1 Numeric vector. First observed series.
#' @param y2 Numeric vector. Second observed series, same length as `y1`.
#' @param trend Numeric vector (length `n`). Estimated shared trend (as
#'   returned by [estimate_trends()]`$trend`). `NA` values are excluded
#'   from residual fitting.
#' @param max_pq Non-negative integer. Maximum AR and MA order to consider
#'   in AIC search. Default 2. A white-noise model (0, 0) is always included
#'   as a fallback.
#'
#' @return A named list with elements `series1` and `series2`, each containing:
#'   \describe{
#'     \item{ar}{Numeric vector of AR coefficients (length p).}
#'     \item{ma}{Numeric vector of MA coefficients (length q).}
#'     \item{sigma2}{Positive numeric. Innovation variance \eqn{\hat\sigma^2}.}
#'     \item{order}{Integer vector \eqn{(p, 0, q)}.}
#'     \item{resid}{Numeric vector of residuals over valid (non-NA) indices.}
#'   }
#'
#' @export
estimate_arma_noise <- function(y1, y2, trend, max_pq = 2L) {
  if (length(y1) != length(y2) || length(y1) != length(trend))
    stop("`y1`, `y2`, and `trend` must all have the same length.")
  max_pq <- as.integer(max_pq)

  valid_t <- which(!is.na(trend))
  if (length(valid_t) < 5L * (max_pq + 1L))
    warning("Very few trend-valid points; ARMA fits may be unstable.")

  resid1 <- y1[valid_t] - trend[valid_t]
  resid2 <- y2[valid_t] - trend[valid_t]

  list(
    series1 = .fit_arma(resid1, max_pq),
    series2 = .fit_arma(resid2, max_pq)
  )
}


# Internal: AIC-based ARMA order selection and fitting.
.fit_arma <- function(resid, max_pq) {
  best_aic   <- Inf
  best_fit   <- NULL
  best_order <- c(0L, 0L, 0L)

  for (p in 0:max_pq) {
    for (q in 0:max_pq) {
      if (p == 0L && q == 0L) next
      converged <- TRUE
      fit <- withCallingHandlers(
        try(stats::arima(resid, order = c(p, 0L, q),
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
        best_aic   <- fit$aic
        best_fit   <- fit
        best_order <- c(p, 0L, q)
      }
    }
  }

  if (is.null(best_fit)) {
    best_fit   <- stats::arima(resid, order = c(0L, 0L, 0L))
    best_order <- c(0L, 0L, 0L)
  }

  coef <- best_fit$coef
  list(
    ar     = as.numeric(coef[grepl("^ar", names(coef))]),
    ma     = as.numeric(coef[grepl("^ma", names(coef))]),
    sigma2 = best_fit$sigma2,
    order  = best_order,
    resid  = as.numeric(stats::residuals(best_fit))
  )
}
