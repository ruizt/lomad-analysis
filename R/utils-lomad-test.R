# Internal test implementations for lomad_test() dispatch
#
# .lomad_test_clt()      — pointwise CLT test (paper method)
# .lomad_test_boot()     — [LEGACY] parametric bootstrap
# .lomad_test_mc()       — [LEGACY] Markov-chain bootstrap
# .lomad_test_analytic() — [LEGACY] analytic CLT for frac_state
# .boot_progress()       — progress bar for bootstrap loops


# Pointwise CLT test: the paper method.
.lomad_test_clt <- function(fit, alpha) {
  if (alpha <= 0 || alpha >= 1)
    stop("`alpha` must be in (0, 1).")

  n         <- fit$inputs$n
  s         <- fit$inputs$s
  valid_idx <- fit$valid_idx
  m         <- length(valid_idx)

  if (m == 0L) stop("No valid time points in fit.")

  R   <- fit$R
  rho <- fit$rho
  V   <- fit$V

  Z      <- rep(NA_real_, n)
  p_raw  <- rep(NA_real_, n)

  se                  <- sqrt(V[valid_idx] / s)
  Z[valid_idx]        <- (R[valid_idx] - rho[valid_idx]) / se
  p_raw[valid_idx]    <- stats::pnorm(Z[valid_idx])

  alpha_eff           <- alpha / sum(1 / seq_len(m))
  p_adj               <- rep(NA_real_, n)
  p_adj[valid_idx]    <- pmin(p_raw[valid_idx] * sum(1 / seq_len(m)), 1)

  rejected            <- rep(NA, n)
  rejected[valid_idx] <- p_raw[valid_idx] <= alpha_eff

  n_rejected <- sum(rejected, na.rm = TRUE)
  message(sprintf("Rejected %d / %d time points at BY-FDR = %.2f (alpha_eff = %.4f)",
                  n_rejected, m, alpha, alpha_eff))

  list(
    Z         = Z,
    p_values  = p_raw,
    p_adj     = p_adj,
    rejected  = rejected,
    alpha_eff = alpha_eff,
    inputs    = list(alpha = alpha, s = s)
  )
}


# [LEGACY] Parametric bootstrap test.
.lomad_test_boot <- function(fit, B, seed, ncores, verbose) {

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
  q          <- null_model$q
  h          <- null_model$h
  ma_kernel  <- rep(1 / q, q)

  sim_pair <- function() {
    e1 <- as.numeric(stats::arima.sim(
      model = list(ar = null_model$series1$ar, ma = null_model$series1$ma),
      n = n, sd = sqrt(null_model$series1$sigma2)))
    e2 <- as.numeric(stats::arima.sim(
      model = list(ar = null_model$series2$ar, ma = null_model$series2$ma),
      n = n, sd = sqrt(null_model$series2$sigma2)))
    list(x1 = trend_hat + e1, x2 = trend_hat + e2)
  }

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
    I_v <- as.integer(R[v] < x_eff)
    m <- length(I_v)
    frac_state <- mean(I_v)
    steps01 <- sum(I_v[-m] == 0L)
    entry_rate <- if (steps01 > 0L)
      sum(I_v[-1L] == 1L & I_v[-m] == 0L) / steps01 else NA_real_
    runs <- rle(I_v)
    run_lengths <- runs$lengths[runs$values == 1L]
    n_entries_sim <- length(run_lengths)
    mean_run_length <- if (n_entries_sim > 0L) mean(run_lengths) else NA_real_
    list(entry_rate = entry_rate, mean_run_length = mean_run_length,
         frac_state = frac_state, n_entries = as.numeric(n_entries_sim),
         var_R = stats::var(R[v], na.rm = TRUE))
  }

  if (!is.null(fit$inputs$block_lengths)) {
    block_lengths <- fit$inputs$block_lengths
    Kb <- length(block_lengths)
    block_trends <- lapply(fit$blocks, `[[`, "trend_hat")
    ma_kernel_b <- rep(1 / q, q)

    sim_pair_block <- function(k) {
      n_k <- block_lengths[k]
      e1 <- as.numeric(stats::arima.sim(
        model = list(ar = null_model$series1$ar, ma = null_model$series1$ma),
        n = n_k, sd = sqrt(null_model$series1$sigma2)))
      e2 <- as.numeric(stats::arima.sim(
        model = list(ar = null_model$series2$ar, ma = null_model$series2$ma),
        n = n_k, sd = sqrt(null_model$series2$sigma2)))
      list(x1 = block_trends[[k]] + e1, x2 = block_trends[[k]] + e2)
    }

    sim_stats_block <- function(x1_k, x2_k) {
      n_k <- length(x1_k)
      s1 <- as.numeric(stats::filter(x1_k, ma_kernel_b, sides = 2))
      s2 <- as.numeric(stats::filter(x2_k, ma_kernel_b, sides = 2))
      R_k <- rep(NA_real_, n_k)
      for (t in seq_len(n_k)) {
        start <- t - h + 1L
        if (start < 1L) next
        w <- start:t
        y1w <- s1[w]; y2w <- s2[w]
        if (any(is.na(y1w)) || any(is.na(y2w))) next
        if (stats::sd(y1w) == 0 || stats::sd(y2w) == 0) next
        r <- suppressWarnings(stats::cor(y1w, y2w))
        if (is.finite(r)) R_k[t] <- r
      }
      v <- which(!is.na(R_k))
      list(I_v = as.integer(R_k[v] < x_eff),
           var_R = stats::var(R_k[v], na.rm = TRUE))
    }

    iter_fn <- function(b) {
      block_results <- lapply(seq_len(Kb), function(k) {
        pair <- sim_pair_block(k)
        sim_stats_block(pair$x1, pair$x2)
      })
      steps01 <- 0L; entries <- 0L; run_lengths <- integer(0L)
      var_R_v <- numeric(Kb); all_I <- list()
      for (k in seq_len(Kb)) {
        I_v <- block_results[[k]]$I_v
        var_R_v[[k]] <- block_results[[k]]$var_R
        all_I[[k]] <- I_v
        m_k <- length(I_v)
        if (m_k < 2L) next
        steps01 <- steps01 + sum(I_v[-m_k] == 0L)
        entries <- entries + sum(I_v[-1L] == 1L & I_v[-m_k] == 0L)
        runs_k <- rle(I_v)
        run_lengths <- c(run_lengths, runs_k$lengths[runs_k$values == 1L])
      }
      frac_state <- mean(unlist(all_I), na.rm = TRUE)
      entry_rate <- if (steps01 > 0L) entries / steps01 else NA_real_
      n_entries_sim <- length(run_lengths)
      mean_run_length <- if (n_entries_sim > 0L) mean(run_lengths) else NA_real_
      list(entry_rate = entry_rate, mean_run_length = mean_run_length,
           frac_state = frac_state, n_entries = as.numeric(n_entries_sim),
           var_R = mean(var_R_v, na.rm = TRUE))
    }
  } else {
    iter_fn <- function(b) {
      pair <- sim_pair()
      sim_stats(pair$x1, pair$x2)
    }
  }

  boot <- if (use_parallel) {
    chunk_size <- max(1L, ncores)
    chunks <- split(seq_len(B), ceiling(seq_len(B) / chunk_size))
    results <- vector("list", B)
    if (verbose) .boot_progress(0L, B)
    for (ch in chunks) {
      results[ch] <- parallel::mclapply(ch, iter_fn, mc.cores = ncores)
      if (verbose) .boot_progress(max(ch), B)
    }
    results
  } else if (verbose) {
    .boot_progress(0L, B)
    lapply(seq_len(B), function(b) { res <- iter_fn(b); .boot_progress(b, B); res })
  } else {
    lapply(seq_len(B), iter_fn)
  }

  entry_v <- Filter(is.finite, vapply(boot, `[[`, numeric(1), "entry_rate"))
  run_v   <- Filter(is.finite, vapply(boot, `[[`, numeric(1), "mean_run_length"))
  frac_v  <- Filter(is.finite, vapply(boot, `[[`, numeric(1), "frac_state"))
  nent_v  <- Filter(is.finite, vapply(boot, `[[`, numeric(1), "n_entries"))
  var_R_v <- Filter(is.finite, vapply(boot, `[[`, numeric(1), "var_R"))

  obs <- fit$observed
  expected <- list(
    entry_rate      = if (length(entry_v)) mean(entry_v) else NA_real_,
    mean_run_length = if (length(run_v))   mean(run_v)   else NA_real_,
    frac_state      = if (length(frac_v))  mean(frac_v)  else NA_real_,
    n_entries       = if (length(nent_v))  mean(nent_v)  else NA_real_
  )
  p_values <- list(
    entry_rate      = if (!is.na(obs$entry_rate) && length(entry_v))
      mean(entry_v >= obs$entry_rate) else NA_real_,
    mean_run_length = if (!is.na(obs$mean_run_length) && length(run_v))
      mean(run_v >= obs$mean_run_length) else NA_real_,
    frac_state      = if (!is.na(obs$frac_state) && length(frac_v))
      mean(frac_v >= obs$frac_state) else NA_real_,
    n_entries       = if (!is.na(obs$n_entries) && length(nent_v))
      mean(nent_v >= obs$n_entries) else NA_real_
  )

  mean_var_R <- if (length(var_R_v)) mean(var_R_v) else NA_real_
  neff_boot  <- if (is.finite(mean_var_R) && mean_var_R > 0) 1 / mean_var_R else NA_real_
  tau2_boot  <- if (is.finite(neff_boot)) h / neff_boot else NA_real_

  message(sprintf("B = %d | p: entry_rate = %s, run_length = %s, frac_state = %s, n_entries = %s",
                  B, signif(p_values$entry_rate, 3), signif(p_values$mean_run_length, 3),
                  signif(p_values$frac_state, 3), signif(p_values$n_entries, 3)))

  null_model$neff <- neff_boot
  null_model$tau2 <- tau2_boot
  null_model$seed <- seed

  list(
    null_model = null_model, observed = obs, expected = expected,
    expected_asymptotic = fit$expected_asymptotic, p_values = p_values,
    inputs = list(n = fit$inputs$n, max_pq = fit$inputs$max_pq,
                  B = B, seed = seed, ncores = ncores),
    trend_hat = fit$trend_hat, ma1 = fit$ma1, ma2 = fit$ma2,
    R = fit$R, thresholds = fit$thresholds, I = fit$I, valid_idx = fit$valid_idx
  )
}


# [LEGACY] Markov-chain bootstrap test.
.lomad_test_mc <- function(fit, B, T_sim, seed, verbose) {

  if (!is.null(seed)) set.seed(seed)

  observed  <- fit$observed
  expectedA <- fit$expected_asymptotic

  pi1_th  <- expectedA$frac_state
  L_th    <- expectedA$mean_run_length
  pi11_th <- if (!is.na(L_th) && L_th > 1) 1 - 1 / L_th else NA_real_
  pi0_th  <- if (!is.na(pi1_th)) 1 - pi1_th else NA_real_
  pi01_th <- if (!is.na(pi1_th) && !is.na(pi11_th) && pi1_th > 0 && pi1_th < 1)
    pi1_th * (1 - pi11_th) / (1 - pi1_th) else NA_real_
  pi00_th <- if (!is.na(pi01_th)) 1 - pi01_th else NA_real_

  if (any(!is.finite(c(pi1_th, pi11_th, pi01_th, pi00_th))))
    stop("Markov transition parameters are not all finite.")
  if (pi1_th <= 0 || pi1_th >= 1 || pi11_th < 0 || pi11_th > 1 ||
      pi01_th < 0 || pi01_th > 1)
    stop("Markov transition parameters fall outside [0, 1].")

  if (is.null(T_sim)) {
    T_sim <- if (!is.null(fit$inputs$block_lengths)) {
      as.integer(stats::median(fit$inputs$block_lengths))
    } else {
      length(fit$valid_idx)
    }
  }

  sim_chain <- function() {
    I <- integer(T_sim)
    I[1] <- stats::rbinom(1, 1, pi1_th)
    for (t in 2:T_sim)
      I[t] <- stats::rbinom(1, 1, if (I[t - 1] == 1L) pi11_th else pi01_th)
    I
  }

  chain_stats <- function(I_v) {
    m <- length(I_v)
    frac_state <- mean(I_v)
    steps01 <- sum(I_v[-m] == 0L)
    entry_rate <- if (steps01 > 0L)
      sum(I_v[-1L] == 1L & I_v[-m] == 0L) / steps01 else NA_real_
    runs <- rle(I_v)
    run_lengths <- runs$lengths[runs$values == 1L]
    n_entries_sim <- length(run_lengths)
    mean_run_length <- if (n_entries_sim > 0L) mean(run_lengths) else NA_real_
    list(entry_rate = entry_rate, mean_run_length = mean_run_length,
         frac_state = frac_state, n_entries = as.numeric(n_entries_sim))
  }

  if (verbose) {
    .boot_progress(0L, B)
    boot <- lapply(seq_len(B), function(b) {
      res <- chain_stats(sim_chain()); .boot_progress(b, B); res
    })
  } else {
    boot <- lapply(seq_len(B), function(b) chain_stats(sim_chain()))
  }

  entry_v <- Filter(is.finite, vapply(boot, `[[`, numeric(1), "entry_rate"))
  run_v   <- Filter(is.finite, vapply(boot, `[[`, numeric(1), "mean_run_length"))
  frac_v  <- Filter(is.finite, vapply(boot, `[[`, numeric(1), "frac_state"))
  nent_v  <- Filter(is.finite, vapply(boot, `[[`, numeric(1), "n_entries"))

  expected <- list(
    entry_rate      = if (length(entry_v)) mean(entry_v) else NA_real_,
    mean_run_length = if (length(run_v))   mean(run_v)   else NA_real_,
    frac_state      = if (length(frac_v))  mean(frac_v)  else NA_real_,
    n_entries       = if (length(nent_v))  mean(nent_v)  else NA_real_
  )
  p_values <- list(
    entry_rate      = if (!is.na(observed$entry_rate) && length(entry_v))
      mean(entry_v >= observed$entry_rate) else NA_real_,
    mean_run_length = if (!is.na(observed$mean_run_length) && length(run_v))
      mean(run_v >= observed$mean_run_length) else NA_real_,
    frac_state      = if (!is.na(observed$frac_state) && length(frac_v))
      mean(frac_v >= observed$frac_state) else NA_real_,
    n_entries       = if (!is.na(observed$n_entries) && length(nent_v))
      mean(nent_v >= observed$n_entries) else NA_real_
  )

  message(sprintf("B = %d | p: entry_rate = %s, run_length = %s, frac_state = %s, n_entries = %s",
                  B, signif(p_values$entry_rate, 3), signif(p_values$mean_run_length, 3),
                  signif(p_values$frac_state, 3), signif(p_values$n_entries, 3)))

  null_model <- fit$null_model
  null_model$transitions <- list(pi1 = pi1_th, pi0 = pi0_th,
                                  pi11 = pi11_th, pi01 = pi01_th, pi00 = pi00_th)

  list(
    null_model = null_model, observed = observed, expected = expected,
    expected_asymptotic = expectedA, p_values = p_values,
    inputs = list(n = fit$inputs$n, max_pq = fit$inputs$max_pq,
                  B = B, T_sim = T_sim, seed = seed),
    trend_hat = fit$trend_hat, ma1 = fit$ma1, ma2 = fit$ma2,
    R = fit$R, thresholds = fit$thresholds, I = fit$I, valid_idx = fit$valid_idx
  )
}


# [LEGACY] Analytic CLT test for frac_state.
.lomad_test_analytic <- function(fit, T_eff) {

  observed  <- fit$observed
  expectedA <- fit$expected_asymptotic

  pi1_th  <- expectedA$frac_state
  L_th    <- expectedA$mean_run_length
  pi11_th <- if (!is.na(L_th) && L_th > 1) 1 - 1 / L_th else NA_real_
  pi0_th  <- if (!is.na(pi1_th)) 1 - pi1_th else NA_real_
  pi01_th <- if (!is.na(pi1_th) && !is.na(pi11_th) && pi1_th > 0 && pi1_th < 1)
    pi1_th * (1 - pi11_th) / (1 - pi1_th) else NA_real_
  pi00_th <- if (!is.na(pi01_th)) 1 - pi01_th else NA_real_

  if (is.null(T_eff)) {
    T_eff <- if (!is.null(fit$inputs$block_lengths)) {
      as.integer(stats::median(fit$inputs$block_lengths))
    } else {
      length(fit$valid_idx)
    }
  }

  p_frac <- NA_real_
  if (is.finite(pi1_th) && pi1_th > 0 && pi1_th < 1 &&
      is.finite(pi11_th) && is.finite(pi00_th) && T_eff > 0) {
    theta <- pi11_th + pi00_th - 1
    if (is.finite(theta) && abs(theta) < 1) {
      var_CLT <- pi1_th * (1 - pi1_th) * (1 + theta) / (1 - theta) / T_eff
      if (var_CLT > 0)
        p_frac <- 1 - stats::pnorm((observed$frac_state - pi1_th) / sqrt(var_CLT))
    }
  }

  message(sprintf("p (frac_state, analytic CLT) = %s", signif(p_frac, 3)))

  null_model <- fit$null_model
  null_model$transitions <- list(pi1 = pi1_th, pi0 = pi0_th,
                                  pi11 = pi11_th, pi01 = pi01_th, pi00 = pi00_th)
  null_model$T_eff <- T_eff

  list(
    null_model = null_model, observed = observed,
    expected = expectedA, expected_asymptotic = expectedA,
    p_values = list(entry_rate = NA_real_, mean_run_length = NA_real_,
                    frac_state = p_frac, n_entries = NA_real_),
    inputs = list(n = fit$inputs$n, max_pq = fit$inputs$max_pq, T_eff = T_eff),
    trend_hat = fit$trend_hat, ma1 = fit$ma1, ma2 = fit$ma2,
    R = fit$R, thresholds = fit$thresholds, I = fit$I, valid_idx = fit$valid_idx
  )
}


# Bootstrap progress bar.
.boot_progress <- function(b, B, width = 40L) {
  filled <- round(b / B * width)
  bar    <- paste0(strrep("=", filled), strrep(" ", width - filled))
  cat(sprintf("\r  [%s] %3d%%  (%d/%d)", bar, round(b / B * 100), b, B))
  if (b == B) cat("\n")
  flush(stdout())
}
