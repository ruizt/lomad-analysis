#' Fit the LOMAD local correlation model
#'
#' Estimates quantities needed for inference on the local correlation between
#' two time series. Runs the CLT pipeline: MA trend estimation, noise
#' fitting, and rolling computation of \eqn{R_t}, \eqn{\rho_t},
#' and \eqn{V_t}.
#'
#' @param x1 Numeric vector. First time series.
#' @param x2 Numeric vector. Second time series, same length as `x1`.
#' @param h Integer or NULL. Smoothing window half-bandwidth. Auto-selected
#'   if `NULL`.
#' @param s Integer or NULL. Rolling correlation window length.
#'   Auto-selected if `NULL`.
#' @param noise_method Character. Noise estimation method: `"ar1"` (default)
#'   uses variogram-based AR(1); `"arma"` uses BIC-selected AR(p) via
#'   [estimate_arma_noise()]. Ignored when `noise_override` is supplied.
#' @param noise_override Optional list with elements `ar`, `ma` (optional),
#'   `sigma2`. If supplied, these noise parameters are used directly instead
#'   of being estimated from the data. Useful for oracle experiments or when
#'   noise is estimated externally (e.g. via [estimate_arma_noise()]). May
#'   also be a list of two such lists (one per series).
#' @param lag_max Integer. Lag truncation for autocovariance sums.
#'   Default 100.
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{method}{Character `"clt"`.}
#'     \item{trend}{Shared trend estimate.}
#'     \item{ma1, ma2}{MA-smoothed series.}
#'     \item{noise}{Per-series noise parameter estimates.}
#'     \item{tau_sq, rho, V, R}{CLT quantities for inference.}
#'     \item{valid_idx}{Integer indices of time points with valid test stats.}
#'     \item{inputs}{List of input parameters (n, h, s, lag_max).}
#'   }
#'
#' @seealso [lomad_test()], [lomad()], [lomad_plot()]
#'
#' @export
lomad_fit <- function(x1             = NULL,
                      x2             = NULL,
                      h              = NULL,
                      s              = NULL,
                      noise_method   = c("ar1", "arma"),
                      noise_override = NULL,
                      lag_max        = 100L) {

  if (is.null(x1) || is.null(x2))
    stop("`x1` and `x2` are required.")
  if (length(x1) != length(x2))
    stop("`x1` and `x2` must have the same length.")

  noise_method <- match.arg(noise_method)

  .lomad_fit_clt(y1 = x1, y2 = x2, h = h, s = s, lag_max = lag_max,
                 noise_method = noise_method, noise_override = noise_override)
}
