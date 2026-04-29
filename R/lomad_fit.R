#' Fit the LOMAD local correlation model
#'
#' Estimates quantities needed for inference on the local correlation between
#' two time series. The default `method = "clt"` runs the paper's pipeline:
#' MA trend estimation, variogram-based AR(1) noise fitting, and rolling
#' computation of \eqn{R_t}, \eqn{\rho_t}, and \eqn{V_t}. The legacy
#' `method = "state"` fits a state-process model with ARMA noise and
#' Fisher-z thresholding.
#'
#' @param x1 Numeric vector. First time series.
#' @param x2 Numeric vector. Second time series, same length as `x1`.
#' @param method Character. `"clt"` (default, paper method) or `"state"`
#'   (legacy state-process pipeline).
#' @param h Integer or NULL. Smoothing window width. Auto-selected if `NULL`.
#' @param s Integer or NULL. Rolling correlation window length (CLT only).
#'   Auto-selected if `NULL`.
#' @param noise_override Optional list with elements `ar`, `ma`, `sigma2`
#'   (CLT only). If supplied, these noise parameters are used directly instead
#'   of being estimated from the data. Useful for oracle experiments where the
#'   true noise process is known. May also be a list of two such lists (one per
#'   series).
#' @param lag_max Integer. Lag truncation for autocovariance sums (CLT only).
#'   Default 100.
#' @param q Integer or NULL. MA trend window (state only). Auto-selected if
#'   `NULL`.
#' @param alpha Numeric. FDR level (state only, default 0.05).
#' @param rho0 Numeric. Null correlation (state only, default 0).
#' @param max_pq Integer. Maximum ARMA order (state only, default 2).
#'
#' @return A named list whose contents depend on `method`. All fit objects
#'   include a `$method` field (`"clt"` or `"state"`) used by [lomad_test()]
#'   for automatic dispatch. See the internal implementations for full details.
#'
#' @seealso [lomad_test()], [lomad()], [lomad_plot()]
#'
#' @export
lomad_fit <- function(x1             = NULL,
                      x2             = NULL,
                      method         = c("clt", "state"),
                      noise_override = NULL,
                      h              = NULL,
                      s              = NULL,
                      lag_max        = 100L,
                      q              = NULL,
                      alpha          = 0.05,
                      rho0           = 0,
                      max_pq         = 2) {

  method <- match.arg(method)

  if (is.null(x1) || is.null(x2))
    stop("`x1` and `x2` are required.")
  if (length(x1) != length(x2))
    stop("`x1` and `x2` must have the same length.")

  switch(method,
    clt   = .lomad_fit_clt(y1 = x1, y2 = x2, h = h, s = s, lag_max = lag_max,
                           noise_override = noise_override),
    state = .lomad_fit_state(x1 = x1, x2 = x2, q = q, h = h,
                             alpha = alpha, rho0 = rho0, max_pq = max_pq)
  )
}
