#' Test for local correlation decoupling
#'
#' Given the output of [lomad_fit()], tests whether the local sample
#' correlation is significantly below the estimated population correlation.
#' For CLT fits (the paper method), runs the BY-corrected pointwise Z-test
#' automatically. For legacy state fits, dispatches to `"boot"`, `"mc"`, or
#' `"analytic"` test methods.
#'
#' @param fit List returned by [lomad_fit()].
#' @param method Character or NULL. Test method to use. If `NULL` (default),
#'   auto-dispatches based on `fit$method`: CLT fits use the pointwise Z-test;
#'   state fits default to `"boot"`. For state fits, can be `"boot"`, `"mc"`,
#'   or `"analytic"`.
#' @param alpha Numeric. FDR level for BY correction (CLT, default 0.05).
#' @param B Integer. Bootstrap replicates (boot/mc, default 500).
#' @param seed Integer or NULL. RNG seed (boot/mc).
#' @param ncores Integer. Parallel cores for boot (default 1).
#' @param verbose Logical. Progress bar for boot/mc (default FALSE).
#' @param T_sim Integer or NULL. Simulated chain length (mc only).
#' @param T_eff Integer or NULL. Effective chain length (analytic only).
#'
#' @return A named list whose contents depend on the test method. See
#'   `.lomad_test_clt`, `.lomad_test_boot`, `.lomad_test_mc`, or
#'   `.lomad_test_analytic` for details.
#'
#' @seealso [lomad_fit()], [lomad()]
#'
#' @export
lomad_test <- function(fit,
                       method  = NULL,
                       alpha   = 0.05,
                       B       = 500,
                       seed    = NULL,
                       ncores  = 1L,
                       verbose = FALSE,
                       T_sim   = NULL,
                       T_eff   = NULL) {

  fit_method <- fit$method
  if (is.null(fit_method))
    stop("fit object has no `$method` field. Was it created by `lomad_fit()`?")

  if (fit_method == "clt") {
    if (!is.null(method) && method != "clt")
      warning("CLT fit only supports the CLT test; ignoring `method`.")
    return(.lomad_test_clt(fit, alpha = alpha))
  }

  # State pipeline: dispatch to legacy test methods
  if (is.null(method)) method <- "boot"
  method <- match.arg(method, c("boot", "mc", "analytic"))

  switch(method,
    boot     = .lomad_test_boot(fit, B = B, seed = seed, ncores = ncores,
                                verbose = verbose),
    mc       = .lomad_test_mc(fit, B = B, T_sim = T_sim, seed = seed,
                              verbose = verbose),
    analytic = .lomad_test_analytic(fit, T_eff = T_eff)
  )
}
