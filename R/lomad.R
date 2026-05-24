#' Detect and test local correlation decoupling
#'
#' Convenience wrapper that calls [lomad_fit()] followed by [lomad_test()].
#' Returns the combined fit and test results.
#'
#' @param x1 Numeric vector. First time series.
#' @param x2 Numeric vector. Second time series, same length as `x1`.
#' @param alpha Numeric. FDR level (default 0.05).
#' @param ... Additional arguments passed to [lomad_fit()] (e.g. `h`, `s`,
#'   `lag_max`, `noise_method`, `noise_override`).
#'
#' @return A list with elements `$fit` ([lomad_fit()] result) and `$test`
#'   ([lomad_test()] result).
#'
#' @seealso [lomad_fit()], [lomad_test()], [lomad_plot()]
#'
#' @export
lomad <- function(x1, x2, alpha = 0.05, ...) {
  fit <- lomad_fit(x1, x2, ...)
  tst <- lomad_test(fit, alpha = alpha)
  list(fit = fit, test = tst)
}
