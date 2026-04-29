#' Estimate AR(1) noise parameters via the variogram
#'
#' Estimates the AR(1) coefficient and innovation variance for each series
#' using the lag-1 and lag-2 sample variograms of shared-trend residuals.
#' The variogram at short lags is dominated by the noise (not the trend),
#' making the estimates robust to trend differences when \eqn{d > 0}.
#'
#' This is a special case of the difference-based approach of
#' Hall and Van Keilegom (2003) for the AR(1) model. For general AR(p)
#' noise, see [estimate_arma_noise()].
#'
#' Residuals are computed as \eqn{y_k - \hat s_t}, where \eqn{\hat s_t} is
#' the shared trend estimate from [estimate_trends()].  The variogram
#' \deqn{V(l) = \frac{1}{2(n - l)} \sum_{t=1}^{n-l} (r_{t+l} - r_t)^2}
#' satisfies \eqn{V(l) = \gamma(0) - \gamma(l)} for a stationary process.
#' For AR(1) with coefficient \eqn{\phi} and innovation variance
#' \eqn{\sigma^2}: \eqn{V(1) = \sigma^2 / (1 + \phi)} and
#' \eqn{V(2) = \sigma^2}, giving \eqn{\phi = V(2)/V(1) - 1}.
#'
#' @param y1 Numeric vector. First observed series.
#' @param y2 Numeric vector. Second observed series, same length as `y1`.
#' @param trend Numeric vector (length `n`). Shared trend estimate
#'   \eqn{\hat s_t} as returned by [estimate_trends()]`$trend`.
#'
#' @return A named list with elements `series1` and `series2`, each containing:
#'   \describe{
#'     \item{ar}{Numeric scalar. Estimated AR(1) coefficient \eqn{\hat\phi}.}
#'     \item{ma}{`numeric(0)` (no MA component).}
#'     \item{sigma2}{Positive numeric. Innovation variance \eqn{\hat\sigma^2}.}
#'     \item{order}{Integer vector `c(1, 0, 0)`.}
#'   }
#'
#' @references
#' Hall, P. and Van Keilegom, I. (2003). Using difference-based methods for
#' inference in nonparametric regression with time series errors. \emph{Journal
#' of the Royal Statistical Society Series B}, 65(2), 443--456.
#'
#' @seealso [estimate_arma_noise()] for general AR(p) noise estimation.
#'
#' @export
estimate_ar1_noise <- function(y1, y2, trend) {
  if (length(y1) != length(y2) || length(y1) != length(trend))
    stop("`y1`, `y2`, and `trend` must all have the same length.")

  valid_t <- which(!is.na(trend))
  if (length(valid_t) < 10L)
    stop("Fewer than 10 valid trend points; cannot estimate noise.")

  resid1 <- y1[valid_t] - trend[valid_t]
  resid2 <- y2[valid_t] - trend[valid_t]

  list(
    series1 = .variogram_ar1(resid1),
    series2 = .variogram_ar1(resid2)
  )
}


#' Estimate AR(p) noise parameters via the multi-lag variogram
#'
#' Generalises [estimate_ar1_noise()] to AR(p) noise using multi-lag
#' variograms and Yule-Walker estimation, following the difference-based
#' framework of Hall and Van Keilegom (2003).
#'
#' The variogram \eqn{V(l) = \gamma(0) - \gamma(l)} is computed at lags
#' \eqn{1, \ldots, L}.  At sufficiently large \eqn{L}, \eqn{V(L) \approx
#' \gamma(0)}, giving an estimate of the process variance.  The recovered
#' autocovariances \eqn{\hat\gamma(l) = \hat\gamma(0) - V(l)} are then
#' plugged into the Yule-Walker equations to obtain AR coefficients and the
#' innovation variance.
#'
#' If `p` is supplied, the AR order is fixed.
#' If `p = NULL` (default), BIC selects among orders \eqn{1, \ldots, p_{max}}.
#'
#' @param y1 Numeric vector. First observed series.
#' @param y2 Numeric vector. Second observed series, same length as `y1`.
#' @param trend Numeric vector (length `n`). Shared trend estimate
#'   \eqn{\hat s_t} as returned by [estimate_trends()]`$trend`.
#' @param p Integer or `NULL`. If specified, fit AR(p) directly (no order
#'   selection). If `NULL`, select order via BIC.
#' @param p_max Positive integer. Maximum AR order for BIC selection.
#'   Ignored when `p` is specified. Default 5.
#' @param L Positive integer. Maximum variogram lag used to estimate
#'   \eqn{\gamma(0)}. Default 20.
#'
#' @return A named list with elements `series1` and `series2`, each containing:
#'   \describe{
#'     \item{ar}{Numeric vector of AR coefficients.}
#'     \item{ma}{`numeric(0)` (no MA component in the AR approximation).}
#'     \item{sigma2}{Positive numeric. Innovation variance \eqn{\hat\sigma^2}.}
#'     \item{order}{Integer vector `c(p, 0, 0)`.}
#'   }
#'
#' @references
#' Hall, P. and Van Keilegom, I. (2003). Using difference-based methods for
#' inference in nonparametric regression with time series errors. \emph{Journal
#' of the Royal Statistical Society Series B}, 65(2), 443--456.
#'
#' @seealso [estimate_ar1_noise()] for the AR(1) special case.
#'
#' @export
estimate_arma_noise <- function(y1, y2, trend,
                                p     = NULL,
                                p_max = 5L,
                                L     = 20L) {
  if (length(y1) != length(y2) || length(y1) != length(trend))
    stop("`y1`, `y2`, and `trend` must all have the same length.")

  p_max <- as.integer(p_max)
  L     <- as.integer(L)
  if (!is.null(p)) p <- as.integer(p)

  valid_t <- which(!is.na(trend))
  if (length(valid_t) < 10L)
    stop("Fewer than 10 valid trend points; cannot estimate noise.")

  resid1 <- y1[valid_t] - trend[valid_t]
  resid2 <- y2[valid_t] - trend[valid_t]

  list(
    series1 = .variogram_ar(resid1, p = p, p_max = p_max, L = L),
    series2 = .variogram_ar(resid2, p = p, p_max = p_max, L = L)
  )
}
