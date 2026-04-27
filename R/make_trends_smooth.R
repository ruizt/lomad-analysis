#' Generate a pair of Fourier-basis trend series with stochastic repulsion
#'
#' Generates two time series that alternate organically between coupled and
#' decoupled states. Coupling is governed by a smooth latent process `lambda`
#' derived from kernel-smoothed Gaussian noise, giving irregular timing and
#' duration of decoupling events rather than a periodic structure.
#'
#' During coupling (lambda near 1), both series track their shared mean.
#' During decoupling (lambda near 0), they repel symmetrically away from it.
#' Lambda is restricted to (0, 1), so the series always stay on their
#' respective sides of the mean and never cross. See \code{\link{make_trends_cross}}
#' for a variant that allows crossings.
#'
#' The output distance `d` is the Euclidean distance between the two resulting
#' series. `lambda` is returned as ground truth for the coupling state at each
#' time point.
#'
#' @param n Integer. Length of the output series.
#' @param nb Integer. Number of Fourier basis functions (must be odd).
#' @param d Numeric. Target Euclidean distance between the two output series
#'   (default 1).
#' @param bw Numeric. Bandwidth of the Gaussian kernel used to smooth the
#'   latent coupling process, in units of time (default 50). Larger values
#'   produce longer, more sustained coupling and decoupling stretches.
#' @param coupling Numeric in (0, 1). Proportion of time the series spend in
#'   a coupled state (lambda > 0.5), approximately (default 0.8).
#' @param sd0 Numeric. Standard deviation of the lowest-frequency Fourier
#'   coefficient (default 2). Higher-frequency coefficients decay as
#'   \code{sd0 / k^p}.
#' @param p Numeric. Spectral decay exponent (default 2.5).
#' @param k_min Integer. Minimum harmonic index to include (default 1). Setting
#'   `k_min > 1` excludes low-frequency components; see [make_trends_dist()] for
#'   details.
#' @param seed Integer or NULL. RNG seed for reproducibility.
#'
#' @return A list with:
#'   \describe{
#'     \item{y1}{Numeric vector of length \code{n}. First output series.}
#'     \item{y2}{Numeric vector of length \code{n}. Second output series.}
#'     \item{x_mean}{Numeric vector of length \code{n}. Shared mean trend.}
#'     \item{lambda}{Numeric vector of length \code{n}. Latent coupling
#'       weights in \eqn{(0, 1)}; values near 1 indicate coupling, near 0
#'       indicate decoupling. This is the ground truth coupling state.}
#'   }
#'
#' @examples
#' sim <- make_trends_smooth(n = 500, d = 2, bw = 50, coupling = 0.8)
#' plot(sim$y1, type = "l")
#' lines(sim$y2, col = "blue")
#' lines(sim$lambda, col = "gray")
#'
#' @export
make_trends_smooth <- function(n        = 500,
                               nb       = 25,
                               d        = 1,
                               bw       = 50,
                               coupling = 0.8,
                               sd0      = 2,
                               p        = 2.5,
                               k_min    = 1L,
                               seed     = NULL) {
  if ((nb %% 2) == 0) stop("`nb` must be odd.")
  if (coupling <= 0 || coupling >= 1) stop("`coupling` must be in (0, 1).")

  if (!is.null(seed)) set.seed(seed)

  # smooth latent process: kernel-smoothed Gaussian noise -> lambda in (0, 1)
  z_raw    <- stats::rnorm(n)
  z_smooth <- stats::ksmooth(seq_len(n), z_raw, kernel = "normal",
                             bandwidth = bw, x.points = seq_len(n))$y
  z        <- (z_smooth - mean(z_smooth)) / stats::sd(z_smooth)
  # shift threshold so ~`coupling` fraction of timepoints have lambda > 0.5
  lambda   <- stats::pnorm(stats::qnorm(coupling) - z)

  # underlying Fourier series at unit coefficient distance
  coefs <- generate_coef_pair(nb = nb, sd0 = sd0, d = 1, p = p,
                              k_min = k_min, seed = NULL)
  fb    <- fda::create.fourier.basis(rangeval = c(0, n), nbasis = nb, period = n)
  Phi   <- fda::eval.basis(seq_len(n), fb)[, -1]
  x1    <- as.numeric(Phi %*% coefs$coef1)
  x2    <- as.numeric(Phi %*% coefs$coef2)
  x_mean <- (x1 + x2) / 2

  # symmetric repulsion at unit scale
  y1_unit <- x_mean + (1 - lambda) * (x1 - x_mean)
  y2_unit <- x_mean + (1 - lambda) * (x2 - x_mean)

  # rescale so that ||y1 - y2|| == d
  r_unit <- sqrt(sum((y1_unit - y2_unit)^2))
  s      <- d / r_unit

  list(
    y1     = x_mean + s * (y1_unit - x_mean),
    y2     = x_mean + s * (y2_unit - x_mean),
    x_mean = x_mean,
    lambda = lambda
  )
}
