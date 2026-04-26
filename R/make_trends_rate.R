#' Generate a pair of Fourier-basis trend series with repulsion at a controlled rate
#'
#' Generates two time series that are coupled most of the time but undergo
#' symmetric repulsion events at a controlled rate. During each event, both
#' series deviate away from their shared mean in opposite directions; between
#' events they track together. The output distance `d` is the Euclidean
#' distance between the two resulting series.
#'
#' Decoupling events are spaced evenly across the series, with each event
#' producing a gamma-shaped repulsion spike. The event width narrows
#' automatically as `rate` increases. The coefficient amplitude is rescaled so
#' that `||y1 - y2|| == d` exactly.
#'
#' @param n Integer. Length of the output series.
#' @param nb Integer. Number of Fourier basis functions (must be odd).
#' @param d Numeric. Target Euclidean distance between the two output series
#'   (default 1).
#' @param rate Numeric. Decoupling events per unit time (default 0.01, i.e. one
#'   event per 100 time points). Determines \code{n_events = round(rate * n)}.
#' @param sd0 Numeric. Standard deviation of the lowest-frequency Fourier
#'   coefficient (default 2). Higher-frequency coefficients decay as
#'   \code{sd0 / k^p}.
#' @param p Numeric. Spectral decay exponent (default 2.5).
#' @param seed Integer or NULL. RNG seed for reproducibility.
#'
#' @return A list with:
#'   \describe{
#'     \item{y1}{Numeric vector of length \code{n}. First output series.}
#'     \item{y2}{Numeric vector of length \code{n}. Second output series.}
#'     \item{x_mean}{Numeric vector of length \code{n}. Shared mean trend.}
#'     \item{lambda}{Numeric vector of length \code{n}. Mixing weights in
#'       \eqn{[0, 1]}; values near 1 indicate coupling, near 0 indicate
#'       decoupling.}
#'   }
#'
#' @examples
#' sim <- make_trends_rate(n = 500, d = 2, rate = 0.01)
#' plot(sim$y1, type = "l")
#' lines(sim$y2, col = "blue")
#'
#' @export
make_trends_rate <- function(n    = 500,
                             nb   = 25,
                             d    = 1,
                             rate = 0.01,
                             sd0  = 2,
                             p    = 2.5,
                             seed = NULL) {
  if ((nb %% 2) == 0) stop("`nb` must be odd.")

  n_events <- round(rate * n)
  if (n_events < 1) stop("`rate * n` must be at least 1; increase `rate` or `n`.")

  # internal decoupling parameters: full-depth dips, width derived from n_events
  delta            <- 1 / n_events
  gap              <- n / n_events
  decouple_length  <- round(gap * delta)   # = round(n / n_events^2)
  decouple_strength <- n_events            # delta * n_events = 1; dips reach 0
  sigma_scale      <- 0.5
  shape            <- 2

  # generate lambda: evenly spaced gamma-shaped dips
  events <- round(seq(from = gap / 2, by = gap, length.out = n_events))
  t      <- seq_len(n)
  lambda <- rep(1, n)
  for (i in events) {
    bump   <- stats::dgamma(t - i, shape = shape, scale = decouple_length * sigma_scale)
    bump   <- bump / max(bump)
    lambda <- lambda - delta * decouple_strength * bump
  }
  lambda <- pmax(0, pmin(1, lambda))

  # generate underlying Fourier series at unit coefficient distance
  coefs <- generate_coef_pair(nb = nb, sd0 = sd0, d = 1, p = p, seed = seed)
  fb    <- fda::create.fourier.basis(rangeval = c(0, n), nbasis = nb, period = n)
  Phi   <- fda::eval.basis(t, fb)[, -1]
  x1    <- as.numeric(Phi %*% coefs$coef1)
  x2    <- as.numeric(Phi %*% coefs$coef2)
  x_mean <- (x1 + x2) / 2

  # repulsion at unit scale
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
