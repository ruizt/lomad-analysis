#' Detect and test local correlation decoupling
#'
#' Convenience wrapper that calls [lomad_fit()] followed by [lomad_test()].
#' Returns the combined fit and test results.
#'
#' @param x1 Numeric vector. First time series.
#' @param x2 Numeric vector. Second time series, same length as `x1`.
#' @param method Character. Inference pipeline: `"clt"` (default, paper
#'   method) or `"state"` (legacy). Passed to [lomad_fit()].
#' @param test_method Character or NULL. Test method within the pipeline.
#'   If `NULL`, auto-selected: CLT fits use the pointwise Z-test; state fits
#'   use `"boot"`. Passed to [lomad_test()].
#' @param alpha Numeric. FDR level (default 0.05).
#' @param ... Additional arguments passed to [lomad_fit()] and/or
#'   [lomad_test()] (e.g. `h`, `s`, `lag_max`, `B`, `seed`, `ncores`).
#'
#' @return A list combining the output of [lomad_fit()] and [lomad_test()].
#'   The fit is stored as `$fit` and the test as `$test`.
#'
#' @seealso [lomad_fit()], [lomad_test()], [lomad_plot()]
#'
#' @export
lomad <- function(x1,
                  x2,
                  method      = c("clt", "state"),
                  test_method = NULL,
                  alpha       = 0.05,
                  ...) {

  method <- match.arg(method)
  dots   <- list(...)

  # Split ... into fit args and test args
  fit_args  <- c(list(x1 = x1, x2 = x2, method = method),
                 dots[names(dots) %in% c("h", "s", "lag_max", "q", "rho0",
                                         "max_pq", "blocks",
                                         "noise_override")])
  test_args <- c(list(method = test_method, alpha = alpha),
                 dots[names(dots) %in% c("B", "seed", "ncores", "verbose",
                                         "T_sim", "T_eff")])

  fit <- do.call(lomad_fit, fit_args)
  tst <- do.call(lomad_test, c(list(fit = fit), test_args))

  list(fit = fit, test = tst)
}
