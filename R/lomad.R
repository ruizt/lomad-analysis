#' Detect and test local correlation decoupling
#'
#' Convenience wrapper around [lomad_fit()] followed by
#' [lomad_test()]. Fits the null model, then runs a parametric
#' bootstrap to test whether the observed decoupling statistics exceed what
#' is expected under the null.
#'
#' @param x1 Numeric vector. First time series.
#' @param x2 Numeric vector. Second time series, same length as `x1`.
#' @param q Integer or NULL. MA window for trend smoothing (see
#'   [lomad_fit()]).
#' @param h Integer or NULL. Rolling correlation window length (see
#'   [lomad_fit()]).
#' @param alpha Numeric. FDR level for BY threshold (default 0.05).
#' @param rho0 Numeric. Null correlation value (default 0).
#' @param max_pq Integer. Maximum ARMA order (default 2).
#' @param B Integer. Number of bootstrap replicates (default 500).
#' @param seed Integer or NULL. RNG seed passed to [lomad_test()].
#'
#' @return The list returned by [lomad_test()], which contains
#'   observed statistics, bootstrap p-values, expected values under the null,
#'   and the full fitted null model.
#'
#' @seealso [lomad_fit()], [lomad_test()], [plot_lomad_fit()]
#'
#' @export
lomad <- function(x1,
                  x2,
                  q      = NULL,
                  h      = NULL,
                  alpha  = 0.05,
                  rho0   = 0,
                  max_pq = 2,
                  B      = 500,
                  seed   = NULL) {

  fit <- lomad_fit(x1, x2,
                   q      = q,
                   h      = h,
                   alpha  = alpha,
                   rho0   = rho0,
                   max_pq = max_pq)

  lomad_test(fit, B = B, seed = seed)
}
