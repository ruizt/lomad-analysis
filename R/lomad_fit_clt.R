#' Fit the CLT-based local correlation model
#'
#' Estimates all quantities needed for pointwise inference on the local
#' correlation between two time series under the signal-plus-noise model.
#' Runs the full pipeline: MA trend estimation, global ARMA noise fitting,
#' autocovariance sums, and rolling computation of \eqn{R_t}, \eqn{\rho_t},
#' and \eqn{V_t} at each valid time point.
#'
#' @param y1 Numeric vector. First observed series.
#' @param y2 Numeric vector. Second observed series, same length as `y1`.
#' @param h Positive integer. MA smoothing window width (paper notation:
#'   \eqn{h}). If `NULL`, chosen automatically as `max(5, floor(n / 200))`.
#' @param s Positive integer. Rolling correlation window length (paper
#'   notation: \eqn{s}). If `NULL`, chosen automatically as
#'   `min(60 * h, floor(n / 4))`.
#' @param max_pq Non-negative integer. Maximum AR/MA order for AIC-based noise
#'   model selection. Default 2.
#' @param lag_max Positive integer. Lag truncation for autocovariance sums
#'   \eqn{L_k}, \eqn{Q_k}, \eqn{Q_{12}}. Default 100, sufficient for AR(1)
#'   up to \eqn{\phi \approx 0.95}.
#'
#' @return A named list:
#'   \describe{
#'     \item{trend}{Estimated shared trend (length `n`).}
#'     \item{ma1, ma2}{MA-smoothed series (length `n`).}
#'     \item{noise}{Output of [estimate_arma_noise()]: ARMA fits for each
#'       series.}
#'     \item{acov_sums}{Output of [acov_sums()]: \eqn{L_1, Q_1, L_2, Q_2,
#'       Q_{12}}.}
#'     \item{tau_sq}{Numeric vector (length `n`). Rolling signal variance
#'       \eqn{\hat\tau^2_t}.}
#'     \item{rho}{Numeric vector (length `n`). Estimated local population
#'       correlation \eqn{\hat\rho_t}.}
#'     \item{V}{Numeric vector (length `n`). Estimated asymptotic variance
#'       \eqn{\hat V_t}.}
#'     \item{R}{Numeric vector (length `n`). Rolling sample correlation
#'       \eqn{R_t}; `NA` where the window is incomplete or degenerate.}
#'     \item{valid_idx}{Integer vector. Indices where \eqn{R_t}, \eqn{\rho_t},
#'       and \eqn{V_t} are all finite.}
#'     \item{inputs}{List of `n`, `h`, `s`, `max_pq`, `lag_max`.}
#'   }
#'
#' @seealso [lomad_test_clt()], [estimate_trends()], [estimate_arma_noise()],
#'   [compute_rho()], [compute_V()], [acov_sums()]
#'
#' @export
lomad_fit_clt <- function(y1, y2,
                           h       = NULL,
                           s       = NULL,
                           max_pq  = 2L,
                           lag_max = 100L) {

  if (length(y1) != length(y2))
    stop("`y1` and `y2` must have the same length.")

  y1 <- as.numeric(y1)
  y2 <- as.numeric(y2)
  n  <- length(y1)

  # --- Window selection ---
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

  # --- Layer 3: estimate trends and ARMA noise ---
  tr    <- estimate_trends(y1, y2, h)
  noise <- estimate_arma_noise(y1, y2, tr$trend, max_pq)

  message(sprintf("ARMA(%d,%d) / ARMA(%d,%d)",
                  noise$series1$order[1], noise$series1$order[3],
                  noise$series2$order[1], noise$series2$order[3]))

  # --- Layer 1: autocovariance sums ---
  # The CLT applies to the MA(h)-smoothed series, so all quantities in
  # Proposition 1 (sigma_k^2, L_k, Q_k, Q_12) must use the ACVF of the
  # *smoothed* noise eta_t = (1/h)*sum_{j=0}^{h-1} Z_{t-j}, not the raw
  # ARMA ACVF.  Compute raw ACVF with extended lag_max so the filter has
  # enough support (needs lags up to lag_max + h - 1).
  raw_lag_max <- lag_max + h - 1L
  acov1_raw <- arma_acov(ar      = noise$series1$ar,
                          ma      = noise$series1$ma,
                          sigma2  = noise$series1$sigma2,
                          lag_max = raw_lag_max)
  acov2_raw <- arma_acov(ar      = noise$series2$ar,
                          ma      = noise$series2$ma,
                          sigma2  = noise$series2$sigma2,
                          lag_max = raw_lag_max)

  # Smooth-filter: gamma_eta(l) = (1/h^2) sum_m (h-|m|) Gamma_Z(|l+m|)
  acov1 <- .ma_filter_acov(acov1_raw, h, lag_max)
  acov2 <- .ma_filter_acov(acov2_raw, h, lag_max)
  sums  <- acov_sums(acov1, acov2)

  # gamma_eta(0) is the smoothed-noise variance; use this as sigma_k^2
  sigma1_sq <- acov1[1L]
  sigma2_sq <- acov2[1L]

  # --- Debias tau^2 ---
  # trend_hat = (MA_h(y1) + MA_h(y2)) / 2 absorbs noise of variance
  # (sigma1_sq + sigma2_sq) / 4 from the smoothed residuals.
  noise_bias <- (sigma1_sq + sigma2_sq) / 4L

  # --- Layer 2: rolling tau^2, rho, V ---
  tau_sq <- pmax(0, compute_tau_sq(tr$trend, s) - noise_bias)
  rho    <- compute_rho(tau_sq, sigma1_sq, sigma2_sq)
  V      <- compute_V(tau_sq, sigma1_sq, sigma2_sq,
                      sums$L1, sums$L2, sums$Q1, sums$Q2, sums$Q12)

  # --- Rolling correlation R_t ---
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
    inputs    = list(n = n, h = h, s = s, max_pq = max_pq, lag_max = lag_max)
  )
}
