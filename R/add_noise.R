#' Add calibrated ARMA noise to a pair of trend series
#'
#' Convenience wrapper around [make_noise()] for use with output from
#' [make_trends_dist()]. Applies noise models to both series and returns the
#' observed series alongside high-level diagnostics.
#'
#' `lambda_target` and `order` can be supplied as a single value (shared across
#' both series) or as a length-2 vector / list (one per series). Seeds for the
#' two calls are derived independently from `seed` so that the two noise
#' processes are not identical.
#'
#' @param trends List returned by [make_trends_dist()], with elements `x1`
#'   and `x2`.
#' @param h Integer. MA window for SNR calibration (passed to [make_noise()]).
#' @param lambda_target Numeric scalar or length-2 vector. Target SNR for each
#'   series. If scalar, the same target is used for both (default 1).
#' @param scale Numeric. Multiplicative scale applied to both trend series
#'   before adding noise (default 1).
#' @param ar.coefs Numeric vector. AR coefficients (see [make_noise()]).
#' @param ma.coefs Numeric vector. MA coefficients (see [make_noise()]).
#' @param order `c(p, q)` integer vector, or a list of two such vectors. If
#'   supplied, auto-generates ARMA coefficients independently for each series
#'   (see [make_noise()]). Cannot be combined with `ar.coefs` / `ma.coefs`.
#' @param s Integer or NULL. Window for local variance estimation (see
#'   [make_noise()]).
#' @param n_start Integer. Burn-in observations for `arima.sim` (default 100).
#' @param seed Integer or NULL. RNG seed for reproducibility. Seeds for the two
#'   series are derived as `seed` and `seed + 1` so their noise processes are
#'   independent.
#'
#' @return A list:
#'   \describe{
#'     \item{y1, y2}{Observed series (signal + noise), as plain numeric
#'       vectors.}
#'     \item{x1, x2}{Scaled trend series (i.e. `scale * trends$x1/x2`).}
#'     \item{noise}{List with `series1` and `series2`, each containing:
#'       `ar`, `ma` (coefficients used), `sigma` (innovation SD), `mean_snr`
#'       (mean local SNR), and `noise_var_ratio` (empirical / theoretical
#'       smooth noise variance; should be near 1).}
#'   }
#'
#' @seealso [make_trends_dist()], [make_noise()]
#'
#' @export
add_noise <- function(trends,
                      h,
                      lambda_target = 1,
                      scale         = 1,
                      ar.coefs      = numeric(0),
                      ma.coefs      = numeric(0),
                      order         = NULL,
                      s             = NULL,
                      n_start       = 100,
                      seed          = NULL) {

  # --- Expand per-series arguments ---

  if (length(lambda_target) == 1L) {
    lambda_target <- c(lambda_target, lambda_target)
  } else if (length(lambda_target) != 2L) {
    stop("`lambda_target` must be a scalar or length-2 vector.")
  }

  if (is.null(order)) {
    order1 <- order2 <- NULL
  } else if (is.list(order)) {
    if (length(order) != 2L)
      stop("`order` as a list must have exactly 2 elements.")
    order1 <- order[[1]]
    order2 <- order[[2]]
  } else {
    # single c(p, q) shared across both
    order1 <- order2 <- order
  }

  seed1 <- if (!is.null(seed)) seed      else NULL
  seed2 <- if (!is.null(seed)) seed + 1L else NULL

  x1 <- scale * trends$x1
  x2 <- scale * trends$x2

  n1 <- make_noise(x.state = x1, h = h, lambda_target = lambda_target[1],
                   ar.coefs = ar.coefs, ma.coefs = ma.coefs, order = order1,
                   s = s, n_start = n_start, seed = seed1)
  n2 <- make_noise(x.state = x2, h = h, lambda_target = lambda_target[2],
                   ar.coefs = ar.coefs, ma.coefs = ma.coefs, order = order2,
                   s = s, n_start = n_start, seed = seed2)

  diag_slim <- function(nm) {
    list(
      ar              = nm$inputs$ar,
      ma              = nm$inputs$ma,
      sigma           = nm$outputs$sigma,
      mean_snr        = nm$diagnostics$mean_snr,
      noise_var_ratio = nm$diagnostics$noise_var_ratio
    )
  }

  list(
    y1    = n1$outputs$y,
    y2    = n2$outputs$y,
    x1    = x1,
    x2    = x2,
    noise = list(
      series1 = diag_slim(n1),
      series2 = diag_slim(n2)
    )
  )
}
