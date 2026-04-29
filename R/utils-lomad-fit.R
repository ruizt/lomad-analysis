# [LEGACY] Internal fit implementations for lomad_fit() dispatch
#
# .lomad_fit_clt()    — CLT pipeline (paper method)
# .lomad_fit_state()  — [LEGACY] state-process pipeline
# .lomad_fit_blocks() — [LEGACY] multi-block state-process pipeline


# CLT pipeline: the paper method.
# Estimates all quantities needed for pointwise inference on R_t.
.lomad_fit_clt <- function(y1, y2, h, s, lag_max, noise_override = NULL) {

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


# [LEGACY] State-process pipeline.
.lomad_fit_state <- function(x1, x2, q, h, alpha, rho0, max_pq) {

  if (length(x1) != length(x2))
    stop("x1 and x2 must have the same length.")
  if (!is.numeric(rho0) || length(rho0) != 1 || rho0 <= -1 || rho0 >= 1)
    stop("rho0 must be a single numeric value in (-1, 1).")
  if (alpha <= 0 || alpha >= 1)
    stop("alpha must be in (0, 1).")

  x1 <- as.numeric(x1)
  x2 <- as.numeric(x2)
  n  <- length(x1)

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

  ma_kernel <- rep(1 / q, q)
  ma1       <- as.numeric(stats::filter(x1, ma_kernel, sides = 2))
  ma2       <- as.numeric(stats::filter(x2, ma_kernel, sides = 2))
  trend_hat <- (ma1 + ma2) / 2

  valid_t <- which(!is.na(trend_hat))
  if (length(valid_t) < 5L * (max_pq + 1L))
    warning("Very few trend-valid points; ARMA fits may be unstable.")

  resid1 <- x1[valid_t] - trend_hat[valid_t]
  resid2 <- x2[valid_t] - trend_hat[valid_t]

  fit1 <- .select_arma(resid1, max_pq)
  fit2 <- .select_arma(resid2, max_pq)

  var_trend <- stats::var(trend_hat[valid_t], na.rm = TRUE)
  lambda1   <- var_trend / stats::var(resid1, na.rm = TRUE)
  lambda2   <- var_trend / stats::var(resid2, na.rm = TRUE)

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

  m         <- length(valid_idx)
  alpha_eff <- alpha / sum(1 / seq_len(m))
  x_eff     <- tanh(z0 + stats::qnorm(alpha_eff) / sqrt(h - 3L))

  I              <- rep(NA_integer_, n)
  I[valid_idx]   <- as.integer(R[valid_idx] < x_eff)
  I_v            <- I[valid_idx]

  frac_state_obs <- mean(I_v)

  steps01        <- sum(I_v[-m] == 0L)
  entries        <- sum(I_v[-1L] == 1L & I_v[-m] == 0L)
  entry_rate_obs <- if (steps01 > 0L) entries / steps01 else NA_real_

  runs               <- rle(I_v)
  run_lengths        <- runs$lengths[runs$values == 1L]
  n_entries          <- length(run_lengths)
  mean_run_length_obs <- if (n_entries > 0L) mean(run_lengths) else NA_real_

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

  list(
    method  = "state",
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


# [LEGACY] Multi-block state-process pipeline.
.lomad_fit_blocks <- function(blocks, q, h, alpha, rho0, max_pq) {

  if (!is.list(blocks) || length(blocks) == 0)
    stop("blocks must be a non-empty list.")
  if (!is.numeric(rho0) || length(rho0) != 1 || rho0 <= -1 || rho0 >= 1)
    stop("rho0 must be a single numeric value in (-1, 1).")
  if (alpha <= 0 || alpha >= 1)
    stop("alpha must be in (0, 1).")

  for (k in seq_along(blocks)) {
    b <- blocks[[k]]
    if (is.null(b$x1) || is.null(b$x2))
      stop(sprintf("Block %d must have x1 and x2 components.", k))
    if (length(b$x1) != length(b$x2))
      stop(sprintf("Block %d: x1 and x2 must have the same length.", k))
  }

  block_lengths_orig <- vapply(blocks, function(b) length(b$x1), integer(1L))
  n_total_orig       <- sum(block_lengths_orig)

  if (is.null(q)) q <- max(5L, floor(n_total_orig / 200L))
  if (is.null(h)) {
    q_tmp <- as.integer(q)
    h_raw <- 60L * q_tmp
    h <- min(max(h_raw, 10L * q_tmp), floor(n_total_orig / 4L))
    message(sprintf("q = %d, h = %d (auto)", q, h))
  }
  q <- as.integer(q)
  h <- as.integer(h)
  if (h <= 3L) stop("h must be > 3.")
  if (q < 1L)  stop("q must be >= 1.")

  min_len <- q + h - 1L

  too_short <- which(block_lengths_orig < min_len)
  if (length(too_short) > 0L) {
    block_names <- if (!is.null(names(blocks))) names(blocks)[too_short] else as.character(too_short)
    warning(sprintf(
      "Dropping %d block(s) with fewer than q + h - 1 = %d observations: %s",
      length(too_short), min_len, paste(block_names, collapse = ", ")
    ))
    blocks <- blocks[-too_short]
    if (length(blocks) == 0L)
      stop("No blocks remain after dropping blocks shorter than q + h - 1.")
  }

  K             <- length(blocks)
  block_lengths <- vapply(blocks, function(b) length(b$x1), integer(1L))
  n_total       <- sum(block_lengths)
  ma_kernel     <- rep(1 / q, q)

  block_fits <- lapply(blocks, function(b) {
    x1        <- as.numeric(b$x1)
    x2        <- as.numeric(b$x2)
    ma1       <- as.numeric(stats::filter(x1, ma_kernel, sides = 2))
    ma2       <- as.numeric(stats::filter(x2, ma_kernel, sides = 2))
    trend_hat <- (ma1 + ma2) / 2
    valid_t   <- which(!is.na(trend_hat))
    list(x1 = x1, x2 = x2, ma1 = ma1, ma2 = ma2,
         trend_hat = trend_hat, valid_t = valid_t,
         resid1 = x1[valid_t] - trend_hat[valid_t],
         resid2 = x2[valid_t] - trend_hat[valid_t])
  })

  sep <- rep(NA_real_, max_pq)
  concat_resids <- function(field) {
    vecs <- lapply(block_fits, `[[`, field)
    Reduce(function(a, b) c(a, sep, b), vecs)
  }

  fit1 <- .select_arma(concat_resids("resid1"), max_pq)
  fit2 <- .select_arma(concat_resids("resid2"), max_pq)

  all_trends <- unlist(lapply(block_fits, function(bf) bf$trend_hat[bf$valid_t]))
  all_resid1 <- unlist(lapply(block_fits, `[[`, "resid1"))
  all_resid2 <- unlist(lapply(block_fits, `[[`, "resid2"))
  var_trend  <- stats::var(all_trends, na.rm = TRUE)
  lambda1    <- var_trend / stats::var(all_resid1, na.rm = TRUE)
  lambda2    <- var_trend / stats::var(all_resid2, na.rm = TRUE)

  z0 <- atanh(rho0)

  block_corr <- lapply(block_fits, function(bf) {
    n_k <- length(bf$x1)
    R_k <- rep(NA_real_, n_k)
    p_k <- rep(NA_real_, n_k)
    for (t in seq_len(n_k)) {
      start <- t - h + 1L
      if (start < 1L) next
      w <- start:t
      y1w <- bf$ma1[w]; y2w <- bf$ma2[w]
      if (any(is.na(y1w)) || any(is.na(y2w))) next
      if (stats::sd(y1w) == 0 || stats::sd(y2w) == 0) next
      r <- suppressWarnings(stats::cor(y1w, y2w))
      if (!is.finite(r) || abs(r) >= 1) next
      R_k[t] <- r
      z_r    <- 0.5 * log((1 + r) / (1 - r))
      p_k[t] <- stats::pnorm((z_r - z0) * sqrt(h - 3L), lower.tail = TRUE)
    }
    list(R = R_k, p = p_k)
  })

  m_total <- sum(vapply(block_corr, function(bc) sum(!is.na(bc$R)), integer(1L)))
  if (m_total < 5L)
    stop("Too few valid rolling correlation windows across all blocks.")

  alpha_eff <- alpha / sum(1 / seq_len(m_total))
  x_eff     <- tanh(z0 + stats::qnorm(alpha_eff) / sqrt(h - 3L))

  block_state <- lapply(block_corr, function(bc) {
    n_k       <- length(bc$R)
    I_k       <- rep(NA_integer_, n_k)
    valid_k   <- which(!is.na(bc$R))
    I_k[valid_k] <- as.integer(bc$R[valid_k] < x_eff)
    list(I = I_k, valid_idx = valid_k)
  })

  I_seqs <- lapply(seq_len(K), function(k) {
    block_state[[k]]$I[block_state[[k]]$valid_idx]
  })

  all_I          <- unlist(I_seqs)
  frac_state_obs <- mean(all_I, na.rm = TRUE)

  total_steps01   <- 0L
  total_entries   <- 0L
  all_run_lengths <- integer(0L)

  for (k in seq_len(K)) {
    I_v <- I_seqs[[k]]
    m_k <- length(I_v)
    if (m_k < 2L) next
    total_steps01   <- total_steps01 + sum(I_v[-m_k] == 0L)
    total_entries   <- total_entries + sum(I_v[-1L] == 1L & I_v[-m_k] == 0L)
    runs_k          <- rle(I_v)
    all_run_lengths <- c(all_run_lengths, runs_k$lengths[runs_k$values == 1L])
  }

  entry_rate_obs      <- if (total_steps01 > 0L) total_entries / total_steps01 else NA_real_
  n_entries           <- length(all_run_lengths)
  mean_run_length_obs <- if (n_entries > 0L) mean(all_run_lengths) else NA_real_

  W_all <- unlist(lapply(seq_len(K), function(k) {
    bf  <- block_fits[[k]]
    idx <- block_state[[k]]$valid_idx
    (bf$ma1 * bf$ma2)[idx]
  }))
  W_all <- W_all[is.finite(W_all)]
  Nw    <- length(W_all)

  neff_approx <- NA_real_
  tau2_approx <- NA_real_

  if (Nw > 10L) {
    Wc      <- W_all - mean(W_all)
    gamma0  <- stats::var(Wc)
    K_lag   <- min(floor(Nw^(1/3)), 50L)
    if (K_lag >= 1L && is.finite(gamma0) && gamma0 > 0) {
      gamma_k <- vapply(seq_len(K_lag),
                        function(k) mean(Wc[seq_len(Nw - k)] * Wc[(k + 1L):Nw]),
                        numeric(1L))
      w_k        <- 1 - seq_len(K_lag) / (K_lag + 1L)
      sigma2_LR  <- gamma0 + 2 * sum(w_k * gamma_k)
      if (is.finite(sigma2_LR) && sigma2_LR > 0) {
        tau2_approx <- sigma2_LR / gamma0
        neff_approx <- h / tau2_approx
      }
    }
  }

  entry_rate_th      <- NA_real_
  mean_run_length_th <- NA_real_
  frac_state_th      <- NA_real_

  if (is.finite(neff_approx) && neff_approx > 0) {
    adj_pairs <- lapply(seq_len(K), function(k) {
      R_v <- block_corr[[k]]$R[block_state[[k]]$valid_idx]
      m_k <- length(R_v)
      if (m_k < 2L) return(NULL)
      list(r1 = R_v[-m_k], r2 = R_v[-1L])
    })
    adj_pairs <- Filter(Negate(is.null), adj_pairs)
    phi_R <- suppressWarnings(
      stats::cor(unlist(lapply(adj_pairs, `[[`, "r1")),
                 unlist(lapply(adj_pairs, `[[`, "r2")),
                 use = "complete.obs")
    )
    if (!is.finite(phi_R)) phi_R <- 0

    z_x_asym  <- (x_eff - rho0) * sqrt(neff_approx)
    p_x       <- stats::pnorm(z_x_asym)
    frac_state_th <- p_x

    if (p_x > 0 && p_x < 1 && requireNamespace("mvtnorm", quietly = TRUE)) {
      sigma_mat  <- matrix(c(1, phi_R, phi_R, 1), 2, 2)
      joint_prob <- as.numeric(mvtnorm::pmvnorm(
        upper = c(z_x_asym, z_x_asym), mean = c(0, 0), sigma = sigma_mat
      ))
      pi11_th            <- joint_prob / p_x
      pi01_th            <- (p_x - joint_prob) / (1 - p_x)
      entry_rate_th      <- pi01_th
      mean_run_length_th <- if (pi11_th < 1) 1 / (1 - pi11_th) else NA_real_
    }
  }

  message(sprintf(
    "K = %d blocks | ARMA(%d,%d)/ARMA(%d,%d) | n_eff ~ %s | x_eff = %s | n_entries = %d/%d",
    K, fit1$order[1], fit1$order[3], fit2$order[1], fit2$order[3],
    signif(neff_approx, 3), signif(x_eff, 3), n_entries, m_total
  ))

  offsets <- c(0L, cumsum(block_lengths[-K]))

  block_out <- lapply(seq_len(K), function(k) {
    list(trend_hat = block_fits[[k]]$trend_hat,
         ma1 = block_fits[[k]]$ma1, ma2 = block_fits[[k]]$ma2,
         R = block_corr[[k]]$R, I = block_state[[k]]$I,
         valid_idx = block_state[[k]]$valid_idx)
  })
  if (!is.null(names(blocks))) names(block_out) <- names(blocks)

  cat_field      <- function(f) unlist(lapply(block_out, `[[`, f))
  valid_idx_glob <- unlist(lapply(seq_len(K), function(k) {
    block_out[[k]]$valid_idx + offsets[k]
  }))

  list(
    method  = "state",
    null_model = list(
      series1 = list(order = fit1$order, ar = fit1$ar, ma = fit1$ma,
                     sigma2 = fit1$sigma2, lambda = lambda1),
      series2 = list(order = fit2$order, ar = fit2$ar, ma = fit2$ma,
                     sigma2 = fit2$sigma2, lambda = lambda2),
      q = q, h = h, alpha = alpha, rho0 = rho0,
      neff_approx = neff_approx, tau2_approx = tau2_approx
    ),
    observed = list(entry_rate = entry_rate_obs,
                    mean_run_length = mean_run_length_obs,
                    frac_state = frac_state_obs, n_entries = n_entries),
    expected_asymptotic = list(entry_rate = entry_rate_th,
                               mean_run_length = mean_run_length_th,
                               frac_state = frac_state_th),
    inputs = list(n = n_total, n_total = n_total,
                  block_lengths = block_lengths, max_pq = max_pq),
    trend_hat = cat_field("trend_hat"),
    ma1 = cat_field("ma1"), ma2 = cat_field("ma2"),
    R = cat_field("R"),
    thresholds = list(alpha_eff = alpha_eff, x_eff = x_eff),
    I = cat_field("I"),
    valid_idx = valid_idx_glob,
    blocks = block_out
  )
}
