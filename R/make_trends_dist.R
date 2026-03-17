# Internal helper: draw Fourier coefficients with spectral decay
generate_fourier_coef <- function(nb, sd0 = 2, p = 2.5) {
  K   <- (nb - 1) / 2
  k   <- rep(1:K, each = 2)
  sd_k <- sd0 / (k^p)
  stats::rnorm(2 * K, mean = 0, sd = sd_k)
}

# Internal helper: generate a pair of coefficient vectors separated by
# distance d using the sphere-to-ellipse transform method
generate_coef_pair <- function(nb, sd0 = 2, d = 1, p = 2.5, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  coef1 <- generate_fourier_coef(nb, sd0, p)

  z <- stats::rnorm(length(coef1))
  u <- z / sqrt(sum(z^2))

  # Ellipsoid axes: decay relaxes as d increases
  axes <- 1 / (seq_along(coef1))^(p - 0.1 * d)
  dir  <- u * axes
  dir  <- dir / sqrt(sum(dir^2))

  coef2 <- coef1 + d * dir

  list(coef1 = coef1, coef2 = coef2)
}


#' Generate a pair of Fourier-basis trend series at a controlled distance
#'
#' Generates two time series from a shared Fourier basis with a controlled
#' separation distance `d` between their coefficient vectors. Larger `d`
#' produces series that are further apart in coefficient space and tend toward
#' lower correlation.
#'
#' Coefficient vectors are drawn with spectral decay (sd proportional to
#' `k^{-p}`), and the second vector is displaced from the first by distance
#' `d` along a direction drawn from a sphere-to-ellipse transform.
#'
#' @param n Integer. Length of the output series.
#' @param nb Integer. Number of Fourier basis functions (must be odd).
#' @param d Numeric. Euclidean distance between the two coefficient vectors
#'   (default 1). Controls how distinct the two series are.
#' @param sd0 Numeric. Standard deviation of the lowest-frequency coefficient
#'   (default 2). Higher-frequency coefficients decay as `sd0 / k^p`.
#' @param p Numeric. Spectral decay exponent (default 2.5).
#' @param seed Integer or NULL. RNG seed for reproducibility.
#'
#' @return A list with:
#'   \describe{
#'     \item{x1}{Numeric vector of length `n`. First evaluated series.}
#'     \item{x2}{Numeric vector of length `n`. Second evaluated series.}
#'     \item{coef1}{Fourier coefficients for `x1`.}
#'     \item{coef2}{Fourier coefficients for `x2`.}
#'   }
#'
#' @examples
#' sim <- make_trends_dist(n = 500, nb = 25, d = 2)
#' plot(sim$x1, type = "l")
#' lines(sim$x2, col = "blue")
#'
#' @export
make_trends_dist <- function(n   = 500,
                                 nb  = 25,
                                 d   = 1,
                                 sd0 = 2,
                                 p   = 2.5,
                                 seed = NULL) {
  if ((nb %% 2) == 0) stop("`nb` must be odd.")

  coefs <- generate_coef_pair(nb = nb, sd0 = sd0, d = d, p = p, seed = seed)

  fb  <- fda::create.fourier.basis(rangeval = c(0, n), nbasis = nb, period = n)
  Phi <- fda::eval.basis(seq_len(n), fb)[, -1]

  x1 <- as.numeric(Phi %*% coefs$coef1)
  x2 <- as.numeric(Phi %*% coefs$coef2)

  list(x1    = x1,
       x2    = x2,
       coef1 = coefs$coef1,
       coef2 = coefs$coef2)
}
