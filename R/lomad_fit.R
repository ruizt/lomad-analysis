#' Fit the null model for local correlation decoupling
#'
#' Estimates a null model for two time series, smooths out the shared trend,
#' fits ARMA residual models, computes a rolling correlation series, and
#' classifies each time point as inside or outside a low-correlation state
#' using a Benjamini-Yekutieli-corrected Fisher z threshold. Returns observed
#' state statistics and asymptotic expectations under the null.
#'
#' @param x1 Numeric vector. First time series.
#' @param x2 Numeric vector. Second time series, same length as `x1`.
#' @param q Integer or NULL. MA window for shared trend smoothing. If `NULL`,
#'   chosen automatically as `max(5, floor(n/200))`.
#' @param h Integer or NULL. Rolling correlation window length. If `NULL`,
#'   chosen automatically to target `n_eff ~ 30`. Must be > 3.
#' @param alpha Numeric. FDR level for the BY multiple testing correction
#'   (default 0.05).
#' @param rho0 Numeric. Null correlation to test against (default 0). Tests
#'   H0: rho = rho0 vs H1: rho < rho0. Must be in (-1, 1).
#' @param max_pq Integer. Maximum AR/MA order for AIC-based model selection
#'   (default 2).
#'
#' @return A named list:
#'   \describe{
#'     \item{null_model}{ARMA fits, window sizes, HAC n_eff/tau2, and
#'       lambda (signal-to-noise) for each series.}
#'     \item{observed}{`entry_rate`, `mean_run_length`, `frac_state`, and
#'       `n_entries` (number of distinct decoupling episodes) from the data.}
#'     \item{expected_asymptotic}{Asymptotic null expectations for the same
#'       statistics.}
#'     \item{inputs}{`n` and `max_pq`.}
#'     \item{trend_hat}{Estimated shared trend (length `n`).}
#'     \item{ma1, ma2}{MA-smoothed series (length `n`).}
#'     \item{R}{Rolling correlation series (length `n`, NA before window fills).}
#'     \item{thresholds}{`alpha_eff` (BY-corrected p-threshold) and `x_eff`
#'       (corresponding correlation threshold).}
#'     \item{I}{Integer vector (length `n`): 1 = low-correlation state, 0 =
#'       normal, NA = no valid estimate.}
#'     \item{valid_idx}{Integer indices where R is defined.}
#'   }
#'
#' @export
lomad_fit <- function(x1,
                           x2,
                           q      = NULL,
                           h      = NULL,
                           alpha  = 0.05,
                           rho0   = 0,
                           max_pq = 2) {

  if (length(x1) != length(x2))
    stop("x1 and x2 must have the same length.")
  if (!is.numeric(rho0) || length(rho0) != 1 || rho0 <= -1 || rho0 >= 1)
    stop("rho0 must be a single numeric value in (-1, 1).")
  if (alpha <= 0 || alpha >= 1)
    stop("alpha must be in (0, 1).")

  x1 <- as.numeric(x1)
  x2 <- as.numeric(x2)
  n  <- length(x1)

  # --- Window selection ---

  if (is.null(q)) q <- max(5L, floor(n / 200L))
  if (is.null(h)) {
    q_tmp <- as.integer(q)
    h_raw <- 60L * q_tmp
    h <- min(max(h_raw, 10L * q_tmp), floor(n / 4L))
    message(sprintf("q = %d, h = %d (auto)", q, h))
  }
  q <- as.integer(q)
  h <- as.integer(h)
  if (h <= 3L) stop("h must be > 3.")
  if (q < 1L)  stop("q must be >= 1.")

  # --- Helper: AIC-based ARMA order selection ---

  select_arma <- function(resid) {
    best_aic <- Inf
    best_fit <- NULL
    best_order <- c(0L, 0L, 0L)
    for (p in 0:max_pq) {
      for (qi in 0:max_pq) {
        if (p == 0L && qi == 0L) next
        fit <- try(stats::arima(resid, order = c(p, 0L, qi)), silent = TRUE)
        if (inherits(fit, "try-error") || is.na(fit$aic)) next
        if (fit$aic < best_aic) {
          best_aic   <- fit$aic
          best_fit   <- fit
          best_order <- c(p, 0L, qi)
        }
      }
    }
    if (is.null(best_fit)) {
      best_fit   <- stats::arima(resid, order = c(0L, 0L, 0L))
      best_order <- c(0L, 0L, 0L)
    }
    coef   <- best_fit$coef
    list(order  = best_order,
         ar     = as.numeric(coef[grepl("^ar", names(coef))]),
         ma     = as.numeric(coef[grepl("^ma", names(coef))]),
         sigma2 = best_fit$sigma2)
  }

  # --- Step 1: Shared trend and ARMA residual models ---

  ma_kernel <- rep(1 / q, q)
  ma1       <- as.numeric(stats::filter(x1, ma_kernel, sides = 2))
  ma2       <- as.numeric(stats::filter(x2, ma_kernel, sides = 2))
  trend_hat <- (ma1 + ma2) / 2

  valid_t <- which(!is.na(trend_hat))
  if (length(valid_t) < 5L * (max_pq + 1L))
    warning("Very few trend-valid points; ARMA fits may be unstable.")

  resid1 <- x1[valid_t] - trend_hat[valid_t]
  resid2 <- x2[valid_t] - trend_hat[valid_t]

  fit1 <- select_arma(resid1)
  fit2 <- select_arma(resid2)

  var_trend <- stats::var(trend_hat[valid_t], na.rm = TRUE)
  lambda1   <- var_trend / stats::var(resid1, na.rm = TRUE)
  lambda2   <- var_trend / stats::var(resid2, na.rm = TRUE)

  # --- Step 2: Rolling correlations and Fisher z p-values ---

  z0       <- atanh(rho0)
  R        <- rep(NA_real_, n)
  p_vals   <- rep(NA_real_, n)

  for (t in seq_len(n)) {
    start <- t - h + 1L
    if (start < 1L) next
    w <- start:t
    y1w <- ma1[w]; y2w <- ma2[w]
    if (any(is.na(y1w)) || any(is.na(y2w))) next
    if (stats::sd(y1w) == 0 || stats::sd(y2w) == 0) next
    r <- suppressWarnings(stats::cor(y1w, y2w))
    if (!is.finite(r) || abs(r) >= 1) next
    R[t]      <- r
    z_r       <- 0.5 * log((1 + r) / (1 - r))
    p_vals[t] <- stats::pnorm((z_r - z0) * sqrt(h - 3L), lower.tail = TRUE)
  }

  valid_idx <- which(!is.na(R) & !is.na(p_vals))
  if (length(valid_idx) < 5L)
    stop("Too few valid rolling correlation windows. Check q and h.")

  # --- Step 3: BY threshold and state classification ---

  m         <- length(valid_idx)
  alpha_eff <- alpha / sum(1 / seq_len(m))
  x_eff     <- tanh(z0 + stats::qnorm(alpha_eff) / sqrt(h - 3L))

  I              <- rep(NA_integer_, n)
  I[valid_idx]   <- as.integer(R[valid_idx] < x_eff)
  I_v            <- I[valid_idx]

  frac_state_obs <- mean(I_v)

  # Entry rate: P(I_t = 1 | I_{t-1} = 0)
  steps01        <- sum(I_v[-m] == 0L)
  entries        <- sum(I_v[-1L] == 1L & I_v[-m] == 0L)
  entry_rate_obs <- if (steps01 > 0L) entries / steps01 else NA_real_

  # Run lengths of decoupling episodes
  runs               <- rle(I_v)
  run_lengths        <- runs$lengths[runs$values == 1L]
  n_entries          <- length(run_lengths)
  mean_run_length_obs <- if (n_entries > 0L) mean(run_lengths) else NA_real_

  # --- Step 4: HAC-based n_eff ---

  W  <- ma1[valid_idx] * ma2[valid_idx]
  W  <- W[is.finite(W)]
  Nw <- length(W)

  neff_approx <- NA_real_
  tau2_approx <- NA_real_

  if (Nw > 10L) {
    Wc     <- W - mean(W)
    gamma0 <- stats::var(Wc)
    K      <- min(floor(Nw^(1/3)), 50L)
    if (K >= 1L && is.finite(gamma0) && gamma0 > 0) {
      gamma_k  <- vapply(1:K, function(k) mean(Wc[1:(Nw-k)] * Wc[(k+1):Nw]),
                         numeric(1L))
      w_k      <- 1 - (1:K) / (K + 1L)
      sigma2_LR <- gamma0 + 2 * sum(w_k * gamma_k)
      if (is.finite(sigma2_LR) && sigma2_LR > 0) {
        tau2_approx <- sigma2_LR / gamma0
        neff_approx <- h / tau2_approx
      }
    }
  }

  # --- Step 5: Asymptotic expectations under the null ---

  entry_rate_th      <- NA_real_
  mean_run_length_th <- NA_real_
  frac_state_th      <- NA_real_

  if (is.finite(neff_approx) && neff_approx > 0) {
    phi_R <- suppressWarnings(
      stats::cor(R[valid_idx][-m], R[valid_idx][-1L], use = "complete.obs")
    )
    if (!is.finite(phi_R)) phi_R <- 0

    z_x_asym  <- (x_eff - rho0) * sqrt(neff_approx)
    p_x       <- stats::pnorm(z_x_asym)
    frac_state_th <- p_x

    if (p_x > 0 && p_x < 1 && requireNamespace("mvtnorm", quietly = TRUE)) {
      sigma_mat  <- matrix(c(1, phi_R, phi_R, 1), 2, 2)
      joint_prob <- as.numeric(
        mvtnorm::pmvnorm(upper = c(z_x_asym, z_x_asym),
                         mean  = c(0, 0),
                         sigma = sigma_mat)
      )
      pi11_th        <- joint_prob / p_x
      pi01_th        <- (p_x - joint_prob) / (1 - p_x)
      entry_rate_th      <- pi01_th
      mean_run_length_th <- if (pi11_th < 1) 1 / (1 - pi11_th) else NA_real_
    }
  }

  message(sprintf(
    "ARMA(%d,%d)/ARMA(%d,%d) | n_eff ~ %s | x_eff = %s | n_entries = %d/%d",
    fit1$order[1], fit1$order[3],
    fit2$order[1], fit2$order[3],
    signif(neff_approx, 3),
    signif(x_eff, 3),
    n_entries, m
  ))

  # --- Output ---

  list(
    null_model = list(
      series1     = list(order = fit1$order, ar = fit1$ar, ma = fit1$ma,
                         sigma2 = fit1$sigma2, lambda = lambda1),
      series2     = list(order = fit2$order, ar = fit2$ar, ma = fit2$ma,
                         sigma2 = fit2$sigma2, lambda = lambda2),
      q           = q,
      h           = h,
      alpha       = alpha,
      rho0        = rho0,
      neff_approx = neff_approx,
      tau2_approx = tau2_approx
    ),
    observed = list(
      entry_rate      = entry_rate_obs,
      mean_run_length = mean_run_length_obs,
      frac_state      = frac_state_obs,
      n_entries       = n_entries
    ),
    expected_asymptotic = list(
      entry_rate      = entry_rate_th,
      mean_run_length = mean_run_length_th,
      frac_state      = frac_state_th
    ),
    inputs    = list(n = n, max_pq = max_pq),
    trend_hat = trend_hat,
    ma1       = ma1,
    ma2       = ma2,
    R         = R,
    thresholds = list(alpha_eff = alpha_eff, x_eff = x_eff),
    I         = I,
    valid_idx = valid_idx
  )
}
