#' Parametric bootstrap test for local correlation decoupling
#'
#' Runs a parametric bootstrap under the fitted null model from
#' [lomad_fit()] to obtain null distributions and p-values for the
#' entry rate, mean run length, and fraction of time in the low-correlation
#' state.
#'
#' @param fit List returned by [lomad_fit()].
#' @param B Integer. Number of bootstrap replicates (default 500).
#' @param seed Integer or NULL. RNG seed for reproducibility.
#' @param ncores Integer. Number of cores for parallel bootstrap via
#'   [parallel::mclapply()] (default 1, i.e. sequential). Ignored on Windows,
#'   where [lapply()] is used regardless. When `ncores > 1`, the
#'   `"L'Ecuyer-CMRG"` RNG is used so that results are reproducible given
#'   `seed`.
#'
#' @return A named list:
#'   \describe{
#'     \item{null_model}{Null model from `fit`, augmented with
#'       bootstrap-estimated `neff` and `tau2`.}
#'     \item{observed}{Observed statistics from `fit` (including `n_entries`).}
#'     \item{expected}{Bootstrap means under the null.}
#'     \item{expected_asymptotic}{Asymptotic expectations from `fit`.}
#'     \item{p_values}{One-sided bootstrap p-values (P(stat >= observed)).}
#'     \item{inputs}{`n`, `max_pq`, `B`, `seed`.}
#'     \item{trend_hat, ma1, ma2, R, thresholds, I, valid_idx}{Passed through
#'       from `fit` for plotting and downstream use.}
#'   }
#'
#' @export
lomad_test_boot <- function(fit,
                                 B      = 500,
                                 seed   = NULL,
                                 ncores = 1L) {

  use_parallel <- ncores > 1L && .Platform$OS.type != "windows"

  if (use_parallel) {
    old_kind <- RNGkind()[1L]
    on.exit(RNGkind(old_kind), add = TRUE)
    RNGkind("L'Ecuyer-CMRG")
  }
  if (!is.null(seed)) set.seed(seed)

  null_model <- fit$null_model
  n          <- fit$inputs$n
  trend_hat  <- fit$trend_hat
  x_eff      <- fit$thresholds$x_eff
  z0         <- atanh(null_model$rho0)

  q  <- null_model$q
  h  <- null_model$h

  ma_kernel <- rep(1 / q, q)

  # --- Helper: simulate one series pair under the null ---

  sim_pair <- function() {
    e1 <- as.numeric(stats::arima.sim(
      model = list(ar = null_model$series1$ar, ma = null_model$series1$ma),
      n = n, sd = sqrt(null_model$series1$sigma2)
    ))
    e2 <- as.numeric(stats::arima.sim(
      model = list(ar = null_model$series2$ar, ma = null_model$series2$ma),
      n = n, sd = sqrt(null_model$series2$sigma2)
    ))
    list(x1 = trend_hat + e1, x2 = trend_hat + e2)
  }

  # --- Helper: state statistics from one simulated pair ---

  sim_stats <- function(x1, x2) {
    s1 <- as.numeric(stats::filter(x1, ma_kernel, sides = 2))
    s2 <- as.numeric(stats::filter(x2, ma_kernel, sides = 2))
    R  <- rep(NA_real_, n)

    for (t in seq_len(n)) {
      start <- t - h + 1L
      if (start < 1L) next
      w <- start:t
      y1w <- s1[w]; y2w <- s2[w]
      if (any(is.na(y1w)) || any(is.na(y2w))) next
      if (stats::sd(y1w) == 0 || stats::sd(y2w) == 0) next
      r <- suppressWarnings(stats::cor(y1w, y2w))
      if (is.finite(r)) R[t] <- r
    }

    v <- which(!is.na(R))
    if (length(v) < 5L)
      return(list(entry_rate = NA_real_, mean_run_length = NA_real_,
                  frac_state = NA_real_, var_R = NA_real_))

    I_v        <- as.integer(R[v] < x_eff)
    m          <- length(I_v)
    frac_state <- mean(I_v)

    steps01    <- sum(I_v[-m] == 0L)
    entry_rate <- if (steps01 > 0L)
      sum(I_v[-1L] == 1L & I_v[-m] == 0L) / steps01 else NA_real_

    runs           <- rle(I_v)
    run_lengths    <- runs$lengths[runs$values == 1L]
    n_entries_sim  <- length(run_lengths)
    mean_run_length <- if (n_entries_sim > 0L) mean(run_lengths) else NA_real_

    list(entry_rate      = entry_rate,
         mean_run_length = mean_run_length,
         frac_state      = frac_state,
         n_entries       = as.numeric(n_entries_sim),
         var_R           = stats::var(R[v], na.rm = TRUE))
  }

  # --- Bootstrap loop ---

  map_fn <- if (use_parallel) {
    function(X, FUN) parallel::mclapply(X, FUN, mc.cores = ncores)
  } else {
    lapply
  }

  if (!is.null(fit$inputs$block_lengths)) {
    # --- Block-aware bootstrap ---
    block_lengths <- fit$inputs$block_lengths
    Kb            <- length(block_lengths)
    block_trends  <- lapply(fit$blocks, `[[`, "trend_hat")
    ma_kernel_b   <- rep(1 / q, q)

    sim_pair_block <- function(k) {
      n_k <- block_lengths[k]
      e1  <- as.numeric(stats::arima.sim(
        model = list(ar = null_model$series1$ar, ma = null_model$series1$ma),
        n = n_k, sd = sqrt(null_model$series1$sigma2)
      ))
      e2  <- as.numeric(stats::arima.sim(
        model = list(ar = null_model$series2$ar, ma = null_model$series2$ma),
        n = n_k, sd = sqrt(null_model$series2$sigma2)
      ))
      list(x1 = block_trends[[k]] + e1, x2 = block_trends[[k]] + e2)
    }

    sim_stats_block <- function(x1_k, x2_k) {
      n_k  <- length(x1_k)
      s1   <- as.numeric(stats::filter(x1_k, ma_kernel_b, sides = 2))
      s2   <- as.numeric(stats::filter(x2_k, ma_kernel_b, sides = 2))
      R_k  <- rep(NA_real_, n_k)
      for (t in seq_len(n_k)) {
        start <- t - h + 1L
        if (start < 1L) next
        w   <- start:t
        y1w <- s1[w]; y2w <- s2[w]
        if (any(is.na(y1w)) || any(is.na(y2w))) next
        if (stats::sd(y1w) == 0 || stats::sd(y2w) == 0) next
        r <- suppressWarnings(stats::cor(y1w, y2w))
        if (is.finite(r)) R_k[t] <- r
      }
      v <- which(!is.na(R_k))
      list(I_v = as.integer(R_k[v] < x_eff), var_R = stats::var(R_k[v], na.rm = TRUE))
    }

    sim_stats_blocks <- function() {
      # Simulate all blocks once, then aggregate statistics
      block_results <- lapply(seq_len(Kb), function(k) {
        pair <- sim_pair_block(k)
        sim_stats_block(pair$x1, pair$x2)
      })

      steps01 <- 0L; entries <- 0L; run_lengths <- integer(0L)
      var_R_v <- numeric(Kb)
      all_I   <- list()

      for (k in seq_len(Kb)) {
        I_v          <- block_results[[k]]$I_v
        var_R_v[[k]] <- block_results[[k]]$var_R
        all_I[[k]]   <- I_v
        m_k          <- length(I_v)
        if (m_k < 2L) next
        steps01     <- steps01 + sum(I_v[-m_k] == 0L)
        entries     <- entries + sum(I_v[-1L] == 1L & I_v[-m_k] == 0L)
        runs_k      <- rle(I_v)
        run_lengths <- c(run_lengths, runs_k$lengths[runs_k$values == 1L])
      }

      frac_state      <- mean(unlist(all_I), na.rm = TRUE)
      entry_rate      <- if (steps01 > 0L) entries / steps01 else NA_real_
      n_entries_sim   <- length(run_lengths)
      mean_run_length <- if (n_entries_sim > 0L) mean(run_lengths) else NA_real_

      list(entry_rate      = entry_rate,
           mean_run_length = mean_run_length,
           frac_state      = frac_state,
           n_entries       = as.numeric(n_entries_sim),
           var_R           = mean(var_R_v, na.rm = TRUE))
    }

    boot <- map_fn(seq_len(B), function(b) sim_stats_blocks())

  } else {
    # --- Single-series bootstrap (original path) ---
    boot <- map_fn(seq_len(B), function(b) {
      pair <- sim_pair()
      sim_stats(pair$x1, pair$x2)
    })
  }

  entry_v  <- Filter(is.finite, vapply(boot, `[[`, numeric(1), "entry_rate"))
  run_v    <- Filter(is.finite, vapply(boot, `[[`, numeric(1), "mean_run_length"))
  frac_v   <- Filter(is.finite, vapply(boot, `[[`, numeric(1), "frac_state"))
  nent_v   <- Filter(is.finite, vapply(boot, `[[`, numeric(1), "n_entries"))
  var_R_v  <- Filter(is.finite, vapply(boot, `[[`, numeric(1), "var_R"))

  # --- Summaries ---

  obs <- fit$observed

  expected <- list(
    entry_rate      = if (length(entry_v)) mean(entry_v) else NA_real_,
    mean_run_length = if (length(run_v))   mean(run_v)   else NA_real_,
    frac_state      = if (length(frac_v))  mean(frac_v)  else NA_real_,
    n_entries       = if (length(nent_v))  mean(nent_v)  else NA_real_
  )

  p_values <- list(
    entry_rate      = if (!is.na(obs$entry_rate)      && length(entry_v))
      mean(entry_v >= obs$entry_rate)      else NA_real_,
    mean_run_length = if (!is.na(obs$mean_run_length) && length(run_v))
      mean(run_v   >= obs$mean_run_length) else NA_real_,
    frac_state      = if (!is.na(obs$frac_state)      && length(frac_v))
      mean(frac_v  >= obs$frac_state)      else NA_real_,
    n_entries       = if (!is.na(obs$n_entries)       && length(nent_v))
      mean(nent_v  >= obs$n_entries)       else NA_real_
  )

  mean_var_R  <- if (length(var_R_v)) mean(var_R_v) else NA_real_
  neff_boot   <- if (is.finite(mean_var_R) && mean_var_R > 0) 1 / mean_var_R else NA_real_
  tau2_boot   <- if (is.finite(neff_boot)) h / neff_boot else NA_real_

  message(sprintf(
    "B = %d | p: entry_rate = %s, run_length = %s, frac_state = %s, n_entries = %s",
    B,
    signif(p_values$entry_rate, 3),
    signif(p_values$mean_run_length, 3),
    signif(p_values$frac_state, 3),
    signif(p_values$n_entries, 3)
  ))

  null_model$neff <- neff_boot
  null_model$tau2 <- tau2_boot
  null_model$seed <- seed

  list(
    null_model          = null_model,
    observed            = obs,
    expected            = expected,
    expected_asymptotic = fit$expected_asymptotic,
    p_values            = p_values,
    inputs              = list(n = fit$inputs$n, max_pq = fit$inputs$max_pq,
                               B = B, seed = seed, ncores = ncores),
    trend_hat  = fit$trend_hat,
    ma1        = fit$ma1,
    ma2        = fit$ma2,
    R          = fit$R,
    thresholds = fit$thresholds,
    I          = fit$I,
    valid_idx  = fit$valid_idx
  )
}
