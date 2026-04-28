#' Rolling signal variance (tau-squared) in each window
#'
#' For each time point \eqn{t}, computes the variance of the trend over the
#' left-aligned window \eqn{W_t = \{t - s + 1, \ldots, t\}}:
#' \deqn{\hat\tau^2_t = \frac{1}{s} \sum_{u \in W_t} (s_u - \bar s)^2.}
#' Returns `NA` for positions where fewer than `s` trend values are available.
#'
#' @param trend Numeric vector. Estimated (or true) shared trend.
#' @param s Positive integer. Rolling window length (paper notation: \eqn{s}).
#'
#' @return Numeric vector of length `length(trend)` with `NA` for the first
#'   `s - 1` positions.
#'
#' @export
compute_tau_sq <- function(trend, s) {
  s <- as.integer(s)
  stopifnot(s >= 2L, is.numeric(trend))
  n      <- length(trend)
  tau_sq <- rep(NA_real_, n)
  for (t in s:n) {
    w <- trend[(t - s + 1L):t]
    if (any(is.na(w))) next
    tau_sq[t] <- mean((w - mean(w))^2)
  }
  tau_sq
}


#' Local population correlation under the common-trend approximation
#'
#' Computes the approximation for the local population correlation \eqn{\rho}
#' from Proposition 1:
#' \deqn{\rho \approx
#'   \frac{\tau^2}{\sqrt{(\tau^2 + \sigma_1^2)(\tau^2 + \sigma_2^2)}}.}
#'
#' @param tau_sq Non-negative numeric (scalar or vector). Signal variance
#'   \eqn{\tau^2} over the local window.
#' @param sigma1_sq Positive numeric (scalar). Marginal noise variance for
#'   series 1, \eqn{\sigma_1^2 = \gamma_1(0)}.
#' @param sigma2_sq Positive numeric (scalar). Marginal noise variance for
#'   series 2, \eqn{\sigma_2^2 = \gamma_2(0)}.
#'
#' @return Numeric (same length as `tau_sq`) in \eqn{[0, 1)}.
#'
#' @export
compute_rho <- function(tau_sq, sigma1_sq, sigma2_sq) {
  stopifnot(sigma1_sq > 0, sigma2_sq > 0, all(tau_sq >= 0, na.rm = TRUE))
  B <- tau_sq + sigma1_sq
  D <- tau_sq + sigma2_sq
  tau_sq / sqrt(B * D)
}


#' Asymptotic variance of the local sample correlation (V)
#'
#' Computes the approximation for \eqn{V} from Proposition 1. The pointwise
#' CLT gives \eqn{\sqrt{s}(R_t - \rho_t) \xrightarrow{d} N(0, V_t)}, so the
#' standard error of \eqn{R_t} is \eqn{\sqrt{V_t / s}}.
#'
#' @param tau_sq Non-negative numeric (scalar or vector). Signal variance.
#' @param sigma1_sq Positive numeric. Marginal noise variance, series 1.
#' @param sigma2_sq Positive numeric. Marginal noise variance, series 2.
#' @param L1 Numeric. \eqn{L_1 = \sum_{l} \gamma_1(l)}, long-run variance of
#'   series 1 noise.
#' @param L2 Numeric. \eqn{L_2 = \sum_{l} \gamma_2(l)}.
#' @param Q1 Numeric. \eqn{Q_1 = \sum_{l} \gamma_1(l)^2}.
#' @param Q2 Numeric. \eqn{Q_2 = \sum_{l} \gamma_2(l)^2}.
#' @param Q12 Numeric. \eqn{Q_{12} = \sum_{l} \gamma_1(l)\gamma_2(l)}.
#'
#' @return Numeric (same length as `tau_sq`), non-negative.
#'
#' @export
compute_V <- function(tau_sq, sigma1_sq, sigma2_sq,
                      L1, L2, Q1, Q2, Q12) {
  stopifnot(sigma1_sq > 0, sigma2_sq > 0)
  B  <- tau_sq + sigma1_sq
  D  <- tau_sq + sigma2_sq
  Q12 / (B * D) +
    tau_sq * sigma1_sq^2 * L1 / (B^3 * D) +
    tau_sq * sigma2_sq^2 * L2 / (B * D^3) +
    tau_sq^2 * Q1 / (2 * B^3 * D) +
    tau_sq^2 * Q2 / (2 * B * D^3)
}
