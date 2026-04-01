#' Convert a data frame with a block-ID column into a list of series pairs
#'
#' @param df Data frame.
#' @param x1_col Character. Column name for the first series.
#' @param x2_col Character. Column name for the second series.
#' @param block_col Character. Column name for the block identifier.
#'
#' @return A named list of length K, each element a list with `x1` and `x2`
#'   numeric vectors.
#'
#' @export
blocks_from_df <- function(df, x1_col, x2_col, block_col) {
  ids <- unique(df[[block_col]])
  out <- lapply(ids, function(id) {
    sub <- df[df[[block_col]] == id, ]
    list(x1 = as.numeric(sub[[x1_col]]),
         x2 = as.numeric(sub[[x2_col]]))
  })
  setNames(out, as.character(ids))
}


# Internal: AIC-based ARMA order selection (mirrors lomad_fit's local helper).
# Accepts a residual vector that may contain NAs (handled by arima via Kalman).
.select_arma <- function(resid, max_pq) {
  best_aic   <- Inf
  best_fit   <- NULL
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
  coef <- best_fit$coef
  list(order  = best_order,
       ar     = as.numeric(coef[grepl("^ar", names(coef))]),
       ma     = as.numeric(coef[grepl("^ma", names(coef))]),
       sigma2 = best_fit$sigma2)
}


#' Fit the null model for local correlation decoupling across multiple blocks
#'
#' A multi-block variant of [lomad_fit()]. Treats each block as an independent
#' realization of the same generating process. Trend smoothing, rolling
#' correlations, and state classification are performed **within blocks**; no
#' computation crosses a block boundary. ARMA residual models are estimated
#' jointly by concatenating per-block detrended residuals with NA separators
#' (approximating independent block likelihoods via the Kalman filter in
#' [stats::arima()]). Decoupling statistics are aggregated across blocks,
#' excluding cross-block transitions.
#'
#' @param blocks Named or unnamed list of length K. Each element must be a list
#'   with numeric vectors `x1` and `x2` of equal length. Use
#'   [blocks_from_df()] to construct this from a tidy data frame.
#' @param q Integer or NULL. MA window for shared trend smoothing within each
#'   block. If `NULL`, chosen automatically as `max(5, floor(n_total / 200))`.
#' @param h Integer or NULL. Rolling correlation window length. If `NULL`,
#'   chosen automatically to target `n_eff ~ 30`.
#' @param alpha Numeric. FDR level for the BY multiple testing correction
#'   (default 0.05).
#' @param rho0 Numeric. Null correlation to test against (default 0).
#' @param max_pq Integer. Maximum AR/MA order for AIC-based model selection
#'   (default 2).
#'
#' @return A named list with the same structure as [lomad_fit()], compatible
#'   with [lomad_test_boot()], [lomad_test_mc()], and
#'   [lomad_test_analytic()]. Additional fields:
#'   \describe{
#'     \item{inputs$block_lengths}{Integer vector of length K (block lengths
#'       after dropping short blocks).}
#'     \item{inputs$n_total}{Total observations across all retained blocks.}
#'     \item{blocks}{List of K per-block lists, each containing `trend_hat`,
#'       `ma1`, `ma2`, `R`, `I`, and `valid_idx` (local indices within the
#'       block).}
#'   }
#'   The flat fields `trend_hat`, `ma1`, `ma2`, `R`, `I`, and `valid_idx`
#'   are block-concatenated for compatibility with [plot_lomad_fit()].
#'
#' @seealso [lomad_fit()], [blocks_from_df()]
#'
#' @export
lomad_fit_blocks <- function(blocks,
                             q      = NULL,
                             h      = NULL,
                             alpha  = 0.05,
                             rho0   = 0,
                             max_pq = 2) {

  # --- Input validation ---

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

  # --- Window selection (based on total n across all blocks) ---

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

  # --- Drop blocks shorter than the minimum viable length ---

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

  # --- Step 1: Per-block trend and smoothed series ---

  block_fits <- lapply(blocks, function(b) {
    x1        <- as.numeric(b$x1)
    x2        <- as.numeric(b$x2)
    ma1       <- as.numeric(stats::filter(x1, ma_kernel, sides = 2))
    ma2       <- as.numeric(stats::filter(x2, ma_kernel, sides = 2))
    trend_hat <- (ma1 + ma2) / 2
    valid_t   <- which(!is.na(trend_hat))
    list(x1        = x1,
         x2        = x2,
         ma1       = ma1,
         ma2       = ma2,
         trend_hat = trend_hat,
         valid_t   = valid_t,
         resid1    = x1[valid_t] - trend_hat[valid_t],
         resid2    = x2[valid_t] - trend_hat[valid_t])
  })

  # --- Step 2: Joint ARMA fitting via NA-separated residual concatenation ---
  # NA gaps cause the Kalman filter in arima() to restart state estimates at
  # each block boundary, approximating independent per-block likelihoods.

  sep <- rep(NA_real_, max_pq)

  concat_resids <- function(field) {
    vecs <- lapply(block_fits, `[[`, field)
    Reduce(function(a, b) c(a, sep, b), vecs)
  }

  fit1 <- .select_arma(concat_resids("resid1"), max_pq)
  fit2 <- .select_arma(concat_resids("resid2"), max_pq)

  # Pooled signal-to-noise
  all_trends <- unlist(lapply(block_fits, function(bf) bf$trend_hat[bf$valid_t]))
  all_resid1 <- unlist(lapply(block_fits, `[[`, "resid1"))
  all_resid2 <- unlist(lapply(block_fits, `[[`, "resid2"))
  var_trend  <- stats::var(all_trends, na.rm = TRUE)
  lambda1    <- var_trend / stats::var(all_resid1, na.rm = TRUE)
  lambda2    <- var_trend / stats::var(all_resid2, na.rm = TRUE)

  # --- Step 3: Per-block rolling correlations and Fisher z p-values ---

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

  # BY threshold from pooled valid count across all blocks
  m_total <- sum(vapply(block_corr, function(bc) sum(!is.na(bc$R)), integer(1L)))
  if (m_total < 5L)
    stop("Too few valid rolling correlation windows across all blocks. Check q and h.")

  alpha_eff <- alpha / sum(1 / seq_len(m_total))
  x_eff     <- tanh(z0 + stats::qnorm(alpha_eff) / sqrt(h - 3L))

  # --- Step 4: Per-block state classification ---

  block_state <- lapply(block_corr, function(bc) {
    n_k       <- length(bc$R)
    I_k       <- rep(NA_integer_, n_k)
    valid_k   <- which(!is.na(bc$R))
    I_k[valid_k] <- as.integer(bc$R[valid_k] < x_eff)
    list(I = I_k, valid_idx = valid_k)
  })

  # --- Step 5: Block-aware aggregation of observed statistics ---
  # Cross-block transitions are excluded throughout.

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
    # NOTE: runs reaching the block boundary are included at face value
    # (conservative). Censored-run correction is a planned future improvement.
  }

  entry_rate_obs      <- if (total_steps01 > 0L) total_entries / total_steps01 else NA_real_
  n_entries           <- length(all_run_lengths)
  mean_run_length_obs <- if (n_entries > 0L) mean(all_run_lengths) else NA_real_

  # --- Step 6: HAC n_eff (pooled within-block W = ma1 * ma2) ---

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

  # --- Asymptotic expectations under the null ---

  entry_rate_th      <- NA_real_
  mean_run_length_th <- NA_real_
  frac_state_th      <- NA_real_

  if (is.finite(neff_approx) && neff_approx > 0) {
    # Estimate lag-1 autocorrelation of R using within-block adjacent pairs only
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
    K,
    fit1$order[1], fit1$order[3],
    fit2$order[1], fit2$order[3],
    signif(neff_approx, 3),
    signif(x_eff, 3),
    n_entries, m_total
  ))

  # --- Assemble per-block output and concatenate for compatibility ---

  offsets <- c(0L, cumsum(block_lengths[-K]))

  block_out <- lapply(seq_len(K), function(k) {
    list(
      trend_hat = block_fits[[k]]$trend_hat,
      ma1       = block_fits[[k]]$ma1,
      ma2       = block_fits[[k]]$ma2,
      R         = block_corr[[k]]$R,
      I         = block_state[[k]]$I,
      valid_idx = block_state[[k]]$valid_idx   # local (within-block) indices
    )
  })
  if (!is.null(names(blocks))) names(block_out) <- names(blocks)

  cat_field      <- function(f) unlist(lapply(block_out, `[[`, f))
  valid_idx_glob <- unlist(lapply(seq_len(K), function(k) {
    block_out[[k]]$valid_idx + offsets[k]
  }))

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
    inputs = list(
      n             = n_total,        # for compatibility with test functions
      n_total       = n_total,
      block_lengths = block_lengths,
      max_pq        = max_pq
    ),
    trend_hat  = cat_field("trend_hat"),
    ma1        = cat_field("ma1"),
    ma2        = cat_field("ma2"),
    R          = cat_field("R"),
    thresholds = list(alpha_eff = alpha_eff, x_eff = x_eff),
    I          = cat_field("I"),
    valid_idx  = valid_idx_glob,
    blocks     = block_out
  )
}
