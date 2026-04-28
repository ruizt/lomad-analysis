#' Theoretical autocovariance sequence for an ARMA process
#'
#' Computes the autocovariance sequence \eqn{\gamma(0), \gamma(1), \ldots,
#' \gamma(\text{lag\_max})} for a stationary ARMA(p, q) process with
#' innovation variance \eqn{\sigma^2}.
#'
#' @param ar Numeric vector of AR coefficients (length p; empty vector or
#'   `numeric(0)` for a pure MA process).
#' @param ma Numeric vector of MA coefficients (length q; empty vector or
#'   `numeric(0)` for a pure AR process).
#' @param sigma2 Positive numeric. Innovation variance.
#' @param lag_max Non-negative integer. Maximum lag to compute. Default 100.
#'
#' @return Numeric vector of length `lag_max + 1` containing
#'   \eqn{\gamma(0), \gamma(1), \ldots, \gamma(\text{lag\_max})}.
#'
#' @export
arma_acov <- function(ar = numeric(0), ma = numeric(0),
                      sigma2 = 1, lag_max = 100L) {
  lag_max <- as.integer(lag_max)
  stopifnot(sigma2 > 0, lag_max >= 0L)

  # White noise special case: ARMAacf / ARMAtoMA fail on empty models
  if (length(ar) == 0L && length(ma) == 0L) {
    return(c(sigma2, rep(0, lag_max)))
  }

  # gamma(0) = sigma2 * sum(psi_k^2), where psi are MA-infinity coefficients
  psi    <- c(1, stats::ARMAtoMA(ar = ar, ma = ma, lag.max = lag_max))
  gamma0 <- sigma2 * sum(psi^2)

  # ARMAacf returns the normalised ACF rho(l) = gamma(l)/gamma(0)
  acf_vals <- stats::ARMAacf(ar = ar, ma = ma, lag.max = lag_max,
                              pacf = FALSE)
  gamma0 * acf_vals
}


# Internal: ACVF of the MA(h)-smoothed process.
#
# If Z_t is a stationary process with ACVF gamma_Z(l) and we form
# eta_t = (1/h) * sum_{j=0}^{h-1} Z_{t-j}, then eta_t has ACVF
#
#   gamma_eta(l) = (1/h^2) * sum_{m=-(h-1)}^{h-1} (h - |m|) * gamma_Z(|l + m|)
#
# acov_raw  : vector gamma_Z(0), gamma_Z(1), ... of length >= lag_max + h
# h         : MA smoothing window width (paper notation: h)
# lag_max   : maximum lag to compute for eta
#
# Returns gamma_eta(0), gamma_eta(1), ..., gamma_eta(lag_max).
.ma_filter_acov <- function(acov_raw, h, lag_max) {
  m_vals <- -(h - 1L):(h - 1L)
  wts    <- (h - abs(m_vals)) / h^2
  result <- numeric(lag_max + 1L)
  for (l in 0L:lag_max) {
    lag_vals       <- abs(l + m_vals)          # symmetric: Gamma_Z(-k)=Gamma_Z(k)
    result[l + 1L] <- sum(wts * acov_raw[lag_vals + 1L])
  }
  result
}


# Internal: variance of MA(h)-smoothed process (scalar shortcut, lag 0 only).
.ma_filtered_var <- function(acov, h) {
  .ma_filter_acov(acov, h, lag_max = 0L)[1L]
}


#' Autocovariance sums for the asymptotic variance approximation
#'
#' Computes the three covariance sums that appear in Proposition 1 of the
#' asymptotic variance formula:
#' \deqn{L_k = \sum_{l \in \mathbb{Z}} \gamma_k(l), \quad
#'       Q_k = \sum_{l \in \mathbb{Z}} \gamma_k(l)^2, \quad
#'       Q_{12} = \sum_{l \in \mathbb{Z}} \gamma_1(l)\,\gamma_2(l).}
#' All sums are approximated by truncating at the supplied lag vectors.
#'
#' @param acov1 Numeric vector. Autocovariance sequence for series 1,
#'   \eqn{\gamma_1(0), \gamma_1(1), \ldots}, as returned by [arma_acov()].
#' @param acov2 Numeric vector. Autocovariance sequence for series 2 (same
#'   length as `acov1`). If `NULL`, series 2 is treated identically to
#'   series 1 (i.e. \eqn{Q_{12} = Q_1}).
#'
#' @return A named list with elements `L1`, `Q1`, `L2`, `Q2`, `Q12`.
#'
#' @export
acov_sums <- function(acov1, acov2 = NULL) {
  if (is.null(acov2)) acov2 <- acov1

  if (length(acov1) != length(acov2))
    stop("`acov1` and `acov2` must have the same length.")

  # Each sequence is gamma(0), gamma(1), ..., gamma(M).
  # Full two-sided sum: gamma(0) + 2 * sum_{l=1}^{M} gamma(l)
  two_sided_sum  <- function(g) g[1L] + 2 * sum(g[-1L])
  two_sided_sum2 <- function(g) g[1L]^2 + 2 * sum(g[-1L]^2)
  two_sided_cross <- function(g1, g2) g1[1L] * g2[1L] + 2 * sum(g1[-1L] * g2[-1L])

  list(
    L1  = two_sided_sum(acov1),
    Q1  = two_sided_sum2(acov1),
    L2  = two_sided_sum(acov2),
    Q2  = two_sided_sum2(acov2),
    Q12 = two_sided_cross(acov1, acov2)
  )
}
