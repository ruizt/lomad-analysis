#' Simulate a pair of Fourier-basis trend series with controlled separation
#'
#' Generates two time series from a shared Fourier basis with a target L2
#' distance `d`. The coupling weight `w` controls the time-varying mixing
#' between the two underlying trends: `w = 0` means fully decoupled (maximum
#' separation), `w = 1` means fully coupled (both series track their shared
#' mean). The output is rescaled so that `||x1 - x2||_2 = d` exactly.
#'
#' The coupling weight can be specified in three ways:
#' \enumerate{
#'   \item **Named method** (`method`): one of `"dist"` (static separation,
#'     default), `"smooth"` (stochastic repulsion), `"cross"` (stochastic
#'     crossing), or `"rate"` (periodic event-based decoupling). Method-specific
#'     parameters are passed via `...`.
#'   \item **Numeric vector** (`w`): a precomputed coupling weight of length
#'     `n`, overriding `method`.
#'   \item **Function** (`w`): a function `f(n)` returning a numeric vector of
#'     length `n`, called at generation time.
#' }
#'
#' @param n Integer. Length of the output series (default 500).
#' @param d Numeric. Target L2 distance between the two output series
#'   (default 1).
#' @param method Character. Coupling method when `w` is not supplied. One of
#'   `"dist"` (default), `"smooth"`, `"cross"`, or `"rate"`.
#' @param w Coupling weight: a numeric vector of length `n`, a function
#'   `f(n) -> numeric(n)`, or `NULL` (default, uses `method`).
#' @param nb Integer. Number of Fourier basis functions (must be odd,
#'   default 25).
#' @param sd0 Numeric. Standard deviation of the lowest-frequency Fourier
#'   coefficient (default 2). Higher-frequency coefficients decay as
#'   `sd0 / k^p`.
#' @param p Numeric. Spectral decay exponent (default 2.5).
#' @param k_min Integer. Minimum harmonic index to include (default 1). Setting
#'   `k_min > 1` excludes low-frequency components, concentrating signal power
#'   at shorter periods.
#' @param seed Integer or NULL. RNG seed for reproducibility.
#' @param ... Additional arguments passed to the coupling weight generator when
#'   using a named `method`:
#'   \describe{
#'     \item{`"smooth"`}{`bw` (bandwidth, default 50), `coupling` (fraction of
#'       time in coupled state, default 0.8).}
#'     \item{`"cross"`}{`bw` (bandwidth, default 50), `coupling` (fraction of
#'       time near coupled state, default 0.8).}
#'     \item{`"rate"`}{`rate` (decoupling events per unit time, default 0.01).}
#'   }
#'
#' @return A list with:
#'   \describe{
#'     \item{x1}{Numeric vector of length `n`. First output series.}
#'     \item{x2}{Numeric vector of length `n`. Second output series.}
#'     \item{x_mean}{Numeric vector of length `n`. Shared mean trend.}
#'     \item{w}{Numeric vector of length `n`. Coupling weight used.}
#'   }
#'
#' @examples
#' # Unstructured (static separation)
#' tr <- sim_trends(500, d = 2)
#'
#' # Named method with tuning arguments
#' tr <- sim_trends(500, d = 2, method = "smooth", bw = 50, coupling = 0.8)
#'
#' # Custom coupling weight vector
#' w_custom <- rep(c(0, 1), each = 250)
#' tr <- sim_trends(500, d = 2, w = w_custom)
#'
#' @export
sim_trends <- function(n      = 500,
                       d      = 1,
                       method = c("dist", "smooth", "cross", "rate"),
                       w      = NULL,
                       nb     = 25,
                       sd0    = 2,
                       p      = 2.5,
                       k_min  = 1L,
                       seed   = NULL,
                       ...) {
  if ((nb %% 2) == 0) stop("`nb` must be odd.")

  # Resolve coupling weight
  if (is.null(w)) {
    method <- match.arg(method)
    w <- switch(method,
      dist   = rep(0, n),
      smooth = .make_w_smooth(n, ...),
      cross  = .make_w_cross(n, ...),
      rate   = .make_w_rate(n, ...)
    )
  } else if (is.function(w)) {
    w <- w(n)
  }

  if (length(w) != n) {
    stop("`w` must have length `n` (", n, "), got ", length(w), ".")
  }

  # Build base trends from Fourier basis
  basis <- .make_basis_trends(n = n, nb = nb, sd0 = sd0, d = 1,
                              p = p, k_min = k_min, seed = seed)

  # Mix and rescale
  .apply_w(basis$mu1, basis$mu2, basis$x_mean, w, d)
}
