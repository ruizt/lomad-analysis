# Internal utilities for sim_trends() and sim_noise()
#
# Fourier basis:
#   .generate_fourier_coef()  — draw coefficients with spectral decay
#   .generate_coef_pair()     — pair of coefficient vectors at distance d
#   .make_basis_trends()      — evaluate Fourier basis to get mu1, mu2, x_mean
#
# Coupling weights:
#   .make_w_smooth()          — stochastic repulsion
#   .make_w_cross()           — stochastic crossing
#   .make_w_rate()            — periodic event-based decoupling
#
# Mixing:
#   .apply_w()                — mix mu1, mu2 via w and rescale to target d
#
# Noise:
#   .pacf_to_arma_coefs()     — PACF parameterisation → ARMA coefficients


# ---- Fourier basis -------------------------------------------------------

# Draw Fourier coefficients with spectral decay: sd_k = sd0 / k^p
.generate_fourier_coef <- function(nb, sd0 = 2, p = 2.5, k_min = 1L) {
  K    <- (nb - 1L) / 2L
  k    <- rep(seq_len(K), each = 2L)
  sd_k <- ifelse(k >= k_min, sd0 / (k^p), 0)
  stats::rnorm(2L * K, mean = 0, sd = sd_k)
}

# Generate a pair of coefficient vectors separated by distance d using the
# sphere-to-ellipse transform method.
.generate_coef_pair <- function(nb, sd0 = 2, d = 1, p = 2.5, k_min = 1L,
                                seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  coef1 <- .generate_fourier_coef(nb, sd0, p, k_min)

  z <- stats::rnorm(length(coef1))
  u <- z / sqrt(sum(z^2))

  k_idx <- rep(seq_len((nb - 1L) / 2L), each = 2L)
  axes  <- ifelse(k_idx >= k_min, 1 / (k_idx^max(0, p - 0.1 * d)), 0)
  dir   <- u * axes
  dir   <- dir / sqrt(sum(dir^2))

  coef2 <- coef1 + d * dir

  list(coef1 = coef1, coef2 = coef2)
}

# Evaluate Fourier basis pair: returns mu1, mu2, x_mean
.make_basis_trends <- function(n, nb, sd0, d, p, k_min, seed) {
  coefs <- .generate_coef_pair(nb = nb, sd0 = sd0, d = d, p = p,
                               k_min = k_min, seed = seed)

  fb  <- fda::create.fourier.basis(rangeval = c(0, n), nbasis = nb, period = n)
  Phi <- fda::eval.basis(seq_len(n), fb)[, -1]

  mu1    <- as.numeric(Phi %*% coefs$coef1)
  mu2    <- as.numeric(Phi %*% coefs$coef2)
  x_mean <- (mu1 + mu2) / 2

  list(mu1 = mu1, mu2 = mu2, x_mean = x_mean)
}


# ---- Coupling weight generators -----------------------------------------

# Stochastic repulsion: w in (0, 1), no crossing.
# coupling controls the fraction of time w > 0.5.
.make_w_smooth <- function(n, bw = 50, coupling = 0.8) {
  if (coupling <= 0 || coupling >= 1) stop("`coupling` must be in (0, 1).")
  z_raw    <- stats::rnorm(n)
  z_smooth <- stats::ksmooth(seq_len(n), z_raw, kernel = "normal",
                             bandwidth = bw, x.points = seq_len(n))$y
  z <- (z_smooth - mean(z_smooth)) / stats::sd(z_smooth)
  stats::pnorm(stats::qnorm(coupling) - z)
}

# Stochastic crossing: w centered at 1, can go below 0 or above 1.
# coupling controls the fraction of time |w - 1| < 0.5.
.make_w_cross <- function(n, bw = 50, coupling = 0.8) {
  if (coupling <= 0 || coupling >= 1) stop("`coupling` must be in (0, 1).")
  z_raw    <- stats::rnorm(n)
  z_smooth <- stats::ksmooth(seq_len(n), z_raw, kernel = "normal",
                             bandwidth = bw, x.points = seq_len(n))$y
  z        <- (z_smooth - mean(z_smooth)) / stats::sd(z_smooth)
  sigma_w  <- 0.5 / stats::qnorm((1 + coupling) / 2)
  1 - sigma_w * z
}

# Periodic event-based decoupling: evenly spaced gamma-shaped dips.
.make_w_rate <- function(n, rate = 0.01) {
  n_events <- round(rate * n)
  if (n_events < 1) stop("`rate * n` must be at least 1; increase `rate` or `n`.")

  delta     <- 1 / n_events
  gap       <- n / n_events
  dl        <- round(gap * delta)
  strength  <- n_events

  events <- round(seq(from = gap / 2, by = gap, length.out = n_events))
  t      <- seq_len(n)
  w      <- rep(1, n)
  for (i in events) {
    bump <- stats::dgamma(t - i, shape = 2, scale = dl * 0.5)
    bump <- bump / max(bump)
    w    <- w - delta * strength * bump
  }
  pmax(0, pmin(1, w))
}


# ---- Mixing and rescaling -----------------------------------------------

# Mix mu1, mu2 via coupling weight w and rescale to target distance d.
#
# Mixing formula: x_i = x_mean + (1 - w) * (mu_i - x_mean)
# Then rescale so that ||x1 - x2|| = d.
.apply_w <- function(mu1, mu2, x_mean, w, d) {
  x1_unit <- x_mean + (1 - w) * (mu1 - x_mean)
  x2_unit <- x_mean + (1 - w) * (mu2 - x_mean)

  r_unit <- sqrt(sum((x1_unit - x2_unit)^2))
  s      <- d / r_unit

  list(
    x1     = x_mean + s * (x1_unit - x_mean),
    x2     = x_mean + s * (x2_unit - x_mean),
    x_mean = x_mean,
    w      = w
  )
}


# ---- Noise utilities -----------------------------------------------------

# Convert p PACF values (drawn uniformly from (-0.95, 0.95)) to ARMA
# coefficients via the Durbin-Levinson recursion. All resulting processes
# are causal/invertible by construction.
.pacf_to_arma_coefs <- function(p) {
  if (p == 0) return(numeric(0))
  pacf <- stats::runif(p, min = -0.95, max = 0.95)
  coefs <- pacf[1]
  if (p == 1) return(coefs)
  for (k in 2:p) {
    coefs <- c(coefs - pacf[k] * rev(coefs), pacf[k])
  }
  coefs
}
