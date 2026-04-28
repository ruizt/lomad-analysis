## Power study — single-replicate function + local driver
##
## Usage (local):   source("run.R")
## Usage (HPC):     Rscript tide/array-job.R <array_id> <out_dir>
##
## The run_rep() function is the only entry point called by the HPC scripts.

library(lomad)

# ---- Single-replicate function ----------------------------------------------

#' Run one power replicate
#'
#' @param params Named list with elements: T, d, structure, phi, snr,
#'   struct_param (structure-specific: rate r or bandwidth b), seed.
#' @return Named list of estimands (see design.md).
run_rep <- function(params) {
  set.seed(params$seed)

  T_         <- params$T
  d          <- params$d
  structure  <- params$structure
  phi        <- params$phi
  snr        <- params$snr
  sparam     <- params$struct_param

  # Generate trends for the given structure and separation d
  trends <- switch(structure,
    unstructured = make_trends_dist(n = T_, d = d, seed = params$seed),
    event_rate   = make_trends_rate(n = T_, d = d, rate = sparam,
                                    seed = params$seed),
    repulsion    = make_trends_smooth(n = T_, d = d, bw = sparam,
                                      seed = params$seed),
    crossing     = make_trends_cross(n = T_, d = d, bw = sparam,
                                     seed = params$seed),
    stop("Unknown structure: ", structure)
  )

  # Auto-select windows
  h_win <- max(5L, floor(T_ / 200L))
  s_win <- min(60L * h_win, floor(T_ / 4L))

  # Add AR(1) noise
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

  vi       <- fit$valid_idx
  rejected <- tst$rejected[vi]

  # Coupling weight w_t (NA for unstructured — no temporal structure)
  w <- if (!is.null(trends$w)) trends$w[vi] else rep(NA_real_, length(vi))

  # Estimands
  detected    <- any(rejected, na.rm = TRUE)
  sensitivity <- if (!is.null(trends$w) && any(w < 0.5))
                   mean(rejected[w < 0.5], na.rm = TRUE) else NA_real_
  fdr_val     <- if (!is.null(trends$w) && any(rejected, na.rm = TRUE))
                   mean(w[rejected] >= 0.5, na.rm = TRUE) else NA_real_

  list(
    error       = FALSE,
    T           = T_,
    d           = d,
    structure   = structure,
    phi         = phi,
    snr         = snr,
    struct_param = sparam,
    seed        = params$seed,
    n_valid     = length(vi),
    detected    = detected,
    sensitivity = sensitivity,
    fdr         = fdr_val,
    n_flagged   = sum(rejected, na.rm = TRUE)
  )
}

# ---- Parameter grid ---------------------------------------------------------

grid <- expand.grid(
  T            = c(250L, 500L, 1000L),
  d            = c(0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0),
  structure    = c("unstructured", "event_rate", "repulsion", "crossing"),
  phi          = c(0.3, 0.8),
  snr          = c(1.0, 2.0),   # update λ_lo / λ_hi after calibration
  struct_param = c(NA, NA),     # filled per-structure below
  stringsAsFactors = FALSE
)
# TODO: fill struct_param per structure once levels are finalised

# ---- Local driver -----------------------------------------------------------

RUN_LOCAL <- FALSE

if (RUN_LOCAL) {
  # Spot-check: unstructured only, two d values, small S
  grid_dev <- expand.grid(
    T = 500L, d = c(0, 1.5), structure = "unstructured",
    phi = 0.5, snr = 1.0, struct_param = NA,
    stringsAsFactors = FALSE
  )
  S_local <- 10
  set.seed(8888)
  seeds <- sample.int(1e6, nrow(grid_dev) * S_local)

  results <- vector("list", nrow(grid_dev) * S_local)
  k <- 1L
  for (i in seq_len(nrow(grid_dev))) {
    for (s in seq_len(S_local)) {
      params        <- as.list(grid_dev[i, ])
      params$seed   <- seeds[k]
      results[[k]]  <- run_rep(params)
      k             <- k + 1L
    }
  }

  results_df <- do.call(rbind, lapply(results, as.data.frame))
  print(aggregate(detected ~ d, data = results_df, FUN = mean))
}
