# Internal fit implementation for lomad_fit()
#
# .lomad_fit_clt() — CLT pipeline (paper method)


# CLT pipeline: the paper method.
# Estimates all quantities needed for pointwise inference on R_t.
.lomad_fit_clt <- function(y1, y2, h, s, lag_max,
                          noise_method = "ar1", noise_override = NULL) {

  y1 <- as.numeric(y1)
  y2 <- as.numeric(y2)
  n  <- length(y1)

  if (is.null(h)) h <- max(5L, floor(n / 200L))
  if (is.null(s)) {
    s <- min(60L * as.integer(h), floor(n / 4L))
    message(sprintf("h = %d, s = %d (auto)", h, s))
  }
  h       <- as.integer(h)
  s       <- as.integer(s)
  lag_max <- as.integer(lag_max)
  if (s <= 3L) stop("`s` must be > 3.")
  if (h < 1L)  stop("`h` must be >= 1.")

  tr <- estimate_trends(y1, y2, h)

  if (!is.null(noise_override)) {
    # Accept a single spec (shared) or a list of two (per-series)
    if (!is.null(noise_override$ar)) {
      spec1 <- spec2 <- noise_override
    } else {
      spec1 <- noise_override[[1]]
      spec2 <- noise_override[[2]]
    }
    noise <- list(
      series1 = list(ar = spec1$ar,
                     ma = if (!is.null(spec1$ma)) spec1$ma else numeric(0),
                     sigma2 = spec1$sigma2),
      series2 = list(ar = spec2$ar,
                     ma = if (!is.null(spec2$ma)) spec2$ma else numeric(0),
                     sigma2 = spec2$sigma2)
    )
    message("Using noise_override (oracle parameters)")
  } else if (noise_method == "arma") {
    noise <- estimate_arma_noise(y1, y2, tr$trend)
    message(sprintf("AR(%d) / AR(%d): sigma2 = %.4f / %.4f",
                    length(noise$series1$ar), length(noise$series2$ar),
                    noise$series1$sigma2, noise$series2$sigma2))
  } else {
    noise <- estimate_ar1_noise(y1, y2, tr$trend)
    message(sprintf("AR(1): phi = %.3f / %.3f",
                    noise$series1$ar, noise$series2$ar))
  }

  raw_lag_max <- lag_max + h - 1L
  acov1_raw <- arma_acov(ar = noise$series1$ar, ma = noise$series1$ma,
                          sigma2 = noise$series1$sigma2, lag_max = raw_lag_max)
  acov2_raw <- arma_acov(ar = noise$series2$ar, ma = noise$series2$ma,
                          sigma2 = noise$series2$sigma2, lag_max = raw_lag_max)

  acov1 <- .ma_filter_acov(acov1_raw, h, lag_max)
  acov2 <- .ma_filter_acov(acov2_raw, h, lag_max)
  sums  <- acov_sums(acov1, acov2)

  sigma1_sq <- acov1[1L]
  sigma2_sq <- acov2[1L]

  noise_bias <- (sigma1_sq + sigma2_sq) / 4L

  tau_sq <- pmax(0, compute_tau_sq(tr$trend, s) - noise_bias)
  rho    <- compute_rho(tau_sq, sigma1_sq, sigma2_sq)
  V      <- compute_V(tau_sq, sigma1_sq, sigma2_sq,
                      sums$L1, sums$L2, sums$Q1, sums$Q2, sums$Q12)

  R <- rep(NA_real_, n)
  for (t in s:n) {
    w   <- (t - s + 1L):t
    y1w <- tr$ma1[w]
    y2w <- tr$ma2[w]
    if (any(is.na(y1w)) || any(is.na(y2w))) next
    if (stats::sd(y1w) == 0 || stats::sd(y2w) == 0) next
    r <- suppressWarnings(stats::cor(y1w, y2w))
    if (is.finite(r) && abs(r) < 1) R[t] <- r
  }

  valid_idx <- which(is.finite(R) & is.finite(rho) & is.finite(V) & V > 0)

  list(
    method    = "clt",
    trend     = tr$trend,
    ma1       = tr$ma1,
    ma2       = tr$ma2,
    noise     = noise,
    acov_sums = sums,
    tau_sq    = tau_sq,
    rho       = rho,
    V         = V,
    R         = R,
    valid_idx = valid_idx,
    inputs    = list(n = n, h = h, s = s, lag_max = lag_max)
  )
}
