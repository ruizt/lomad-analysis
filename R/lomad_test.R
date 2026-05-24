#' Test for local correlation decoupling
#'
#' Given the output of [lomad_fit()], tests whether the local sample
#' correlation is significantly below the estimated population correlation
#' using the BY-corrected pointwise Z-test.
#'
#' @param fit List returned by [lomad_fit()].
#' @param alpha Numeric. FDR level for BY correction (default 0.05).
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{Z}{Numeric vector of Z-statistics (NA outside valid indices).}
#'     \item{p_values}{Raw p-values.}
#'     \item{p_adj}{BY-adjusted p-values.}
#'     \item{rejected}{Logical vector of rejections.}
#'     \item{alpha_eff}{Effective significance level after BY correction.}
#'   }
#'
#' @seealso [lomad_fit()], [lomad()]
#'
#' @export
lomad_test <- function(fit, alpha = 0.05) {
  .lomad_test_clt(fit, alpha = alpha)
}
