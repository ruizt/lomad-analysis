#' Detect and test local correlation decoupling
#'
#' Convenience wrapper around [lomad_fit()] followed by one of the test
#' functions. Fits the null model, then tests whether the observed decoupling
#' statistics exceed what is expected under the null.
#'
#' @param x1 Numeric vector. First time series.
#' @param x2 Numeric vector. Second time series, same length as `x1`.
#' @param q Integer or NULL. MA window for trend smoothing (see [lomad_fit()]).
#' @param h Integer or NULL. Rolling correlation window length (see
#'   [lomad_fit()]).
#' @param alpha Numeric. FDR level for BY threshold (default 0.05).
#' @param rho0 Numeric. Null correlation value (default 0).
#' @param max_pq Integer. Maximum ARMA order (default 2).
#' @param method Character. Inference method: `"boot"` (full parametric
#'   bootstrap, default), `"mc"` (Markov-chain bootstrap on state process),
#'   or `"analytic"` (closed-form CLT for `frac_state`).
#' @param B Integer. Number of bootstrap replicates for `"boot"` and `"mc"`
#'   (default 500).
#' @param T_sim Integer or NULL. Simulated chain length for `"mc"`. Defaults
#'   to `length(fit$valid_idx)`.
#' @param T_eff Integer or NULL. Effective chain length for `"analytic"`.
#'   Defaults to `length(fit$valid_idx)`.
#' @param seed Integer or NULL. RNG seed for `"boot"` and `"mc"`.
#' @param ncores Integer. Number of cores for parallel bootstrap when
#'   `method = "boot"` (default 1). Passed to [lomad_test_boot()]; ignored for
#'   other methods and on Windows.
#'
#' @return The list returned by the selected test function
#'   ([lomad_test_boot()], [lomad_test_mc()], or [lomad_test_analytic()]),
#'   which contains observed statistics, p-values, expected values under the
#'   null, and the full fitted null model.
#'
#' @seealso [lomad_fit()], [lomad_test_boot()], [lomad_test_mc()],
#'   [lomad_test_analytic()], [plot_lomad_fit()]
#'
#' @export
lomad <- function(x1,
                  x2,
                  q      = NULL,
                  h      = NULL,
                  alpha  = 0.05,
                  rho0   = 0,
                  max_pq = 2,
                  method = c("boot", "mc", "analytic"),
                  B      = 500,
                  T_sim  = NULL,
                  T_eff  = NULL,
                  seed   = NULL,
                  ncores = 1L) {

  method <- match.arg(method)

  fit <- lomad_fit(x1, x2,
                   q      = q,
                   h      = h,
                   alpha  = alpha,
                   rho0   = rho0,
                   max_pq = max_pq)

  switch(method,
    boot     = lomad_test_boot(fit,     B = B, seed = seed, ncores = ncores),
    mc       = lomad_test_mc(fit,       B = B, T_sim = T_sim, seed = seed),
    analytic = lomad_test_analytic(fit, T_eff = T_eff)
  )
}
