## Calibration study — single-replicate function + local driver
##
## Usage (local):   source("run.R")
## Usage (HPC):     Rscript tide/array-job.R <array_id> <out_dir>
##
## The run_rep() function is the only entry point called by the HPC scripts.

library(lomad)

# ---- Single-replicate function ----------------------------------------------

#' Run one calibration replicate
#'
#' @param params Named list with elements: T (series length), phi (AR1 coef),
#'   snr (target SNR), seed (integer RNG seed).
#' @return Named list of estimands (see design.md).
run_rep <- function(params) {
  set.seed(params$seed)

  T_   <- params$T
  phi  <- params$phi
  snr  <- params$snr

  # Generate shared trend (d = 0)
  trends <- make_trends_dist(n = T_, d = 0, seed = params$seed)

  # Auto-select windows (mirrors lomad_fit_clt defaults)
  h_win <- max(5L, floor(T_ / 200L))
  s_win <- min(60L * h_win, floor(T_ / 4L))

  # Add AR(1) noise at target SNR
  sim <- add_noise(
    trends        = trends,
    h             = h_win,
    lambda_target = snr,
    ar.coefs      = phi,
    seed          = params$seed + 1L
  )

  # Fit + test
  fit <- tryCatch(
    lomad_fit_clt(sim$y1, sim$y2, h = h_win, s = s_win, max_pq = 3L),
    error = function(e) NULL
  )
  if (is.null(fit)) return(list(error = TRUE, params = params))

  tst <- tryCatch(
    lomad_test_clt(fit, alpha = 0.05),
    error = function(e) NULL
  )
  if (is.null(tst)) return(list(error = TRUE, params = params))

  vi <- fit$valid_idx

  list(
    error          = FALSE,
    T              = T_,
    phi            = phi,
    snr            = snr,
    seed           = params$seed,
    h              = fit$inputs$h,
    s              = fit$inputs$s,
    n_valid        = length(vi),
    rejection_rate = if (length(vi) > 0) mean(tst$rejected[vi]) else NA_real_,
    arma_p1        = fit$noise$series1$order[1],
    arma_q1        = fit$noise$series1$order[3],
    arma_p2        = fit$noise$series2$order[1],
    arma_q2        = fit$noise$series2$order[3],
    sigma_hat1     = fit$acov_sums$L1,   # placeholder: update when exposed
    sigma_hat2     = fit$acov_sums$L2
  )
}

# ---- Parameter grid ---------------------------------------------------------

grid <- expand.grid(
  T   = c(250L, 500L, 1000L),
  phi = c(0.0, 0.3, 0.5, 0.8),
  snr = c(0.5, 1.0, 2.0),
  stringsAsFactors = FALSE
)

# ---- Local driver (small S for development) ---------------------------------
#
# Set RUN_LOCAL=TRUE to execute; leave FALSE so sourcing the file is safe.

RUN_LOCAL <- FALSE

if (RUN_LOCAL) {
  S_local <- 20   # small S for quick local checks
  set.seed(1234)
  seeds <- sample.int(1e6, nrow(grid) * S_local)

  results <- vector("list", nrow(grid) * S_local)
  k <- 1L
  for (i in seq_len(nrow(grid))) {
    for (s in seq_len(S_local)) {
      params        <- as.list(grid[i, ])
      params$seed   <- seeds[k]
      results[[k]]  <- run_rep(params)
      k             <- k + 1L
    }
  }

  results_df <- do.call(rbind, lapply(results, as.data.frame))
  print(
    aggregate(rejection_rate ~ T + phi + snr, data = results_df, FUN = mean)
  )
}
