# Validate the perturbation bound (Appendix)
#
# Under near-common trends (|s1_t - s2_t| < eps), rho should be a perturbation
# of the common-trend expression, with remainder O(eps * tau). We verify this
# by varying eps from 0 (shared trend) to large, comparing the empirical mean
# R_t against both the common-trend rho and the perturbation expression.

library(lomad)
library(parallel)

set.seed(2917)

# ---- Parallel setup ---------------------------------------------------------
# The n_rep loop below is embarrassingly parallel (each replicate is
# independent); fork on Unix/macOS (fast, no data copying), fall back to a
# PSOCK cluster on Windows.
n_cores <- max(1L, detectCores() - 1L)
is_windows <- .Platform$OS.type == "windows"
if (is_windows) {
  cl <- makeCluster(n_cores)
  on.exit(stopCluster(cl), add = TRUE)
}
run_reps <- function(FUN, n_rep, ...) {
  if (is_windows) parLapply(cl, seq_len(n_rep), FUN, ...)
  else mclapply(seq_len(n_rep), FUN, ..., mc.cores = n_cores)
}

# ---- Parameters ------------------------------------------------------------

n     <- 2000
h     <- 20
s     <- 150
n_rep <- 500

# Base shared trend
K <- 10
basis <- sapply(1:K, function(k) sin(2 * pi * k * (1:n) / n))
coefs <- 2 / (1:K)^1.2
trend_base <- as.numeric(basis %*% coefs)

# Perturbation: a smooth function with unit sup-norm
perturb <- sin(2 * pi * 3 * (1:n) / n)
perturb <- perturb / max(abs(perturb))

# White noise for simplicity
sigma_sq <- 0.05
sigma_eta <- sigma_sq / h  # smoothed noise variance

eps_vals <- c(0, 0.01, 0.05, 0.1, 0.2, 0.5, 1.0)

# ---- Simulate and compare --------------------------------------------------

cat(sprintf("%-8s  %-12s  %-12s  %-12s\n",
            "eps", "RMSE(common)", "RMSE(actual)", "mean|R-rho0|"))

for (eps in eps_vals) {
  trend1 <- trend_base + (eps / 2) * perturb
  trend2 <- trend_base - (eps / 2) * perturb

  ma_t1 <- as.numeric(stats::filter(trend1, rep(1/h, h), sides = 1))
  ma_t2 <- as.numeric(stats::filter(trend2, rep(1/h, h), sides = 1))

  # Common-trend rho (using midpoint as shared signal)
  ma_mid <- (ma_t1 + ma_t2) / 2
  tau_sq_common <- compute_tau_sq(ma_mid, s)
  rho_common    <- compute_rho(tau_sq_common, sigma_eta, sigma_eta)

  # "Actual" rho from true window moments
  rho_actual <- rep(NA_real_, n)
  for (t in s:n) {
    idx <- (t - s + 1):t
    v1 <- ma_t1[idx]; v2 <- ma_t2[idx]
    if (any(is.na(v1)) || any(is.na(v2))) next
    cov12 <- mean((v1 - mean(v1)) * (v2 - mean(v2)))
    var1  <- mean((v1 - mean(v1))^2) + sigma_eta
    var2  <- mean((v2 - mean(v2))^2) + sigma_eta
    rho_actual[t] <- cov12 / sqrt(var1 * var2)
  }

  # Empirical mean R_t. Each replicate is independent, so run the n_rep loop
  # in parallel and reduce the per-replicate sums/counts afterward. Arguments
  # are passed explicitly (rather than captured by lexical scope) so this
  # works identically whether workers are forks (Unix/macOS) or a PSOCK
  # cluster (Windows).
  one_rep <- function(r, trend1, trend2, n, h, s, sigma_sq) {
    y1 <- trend1 + rnorm(n, sd = sqrt(sigma_sq))
    y2 <- trend2 + rnorm(n, sd = sqrt(sigma_sq))
    m1 <- as.numeric(stats::filter(y1, rep(1/h, h), sides = 1))
    m2 <- as.numeric(stats::filter(y2, rep(1/h, h), sides = 1))
    r_sum   <- rep(0, n)
    r_count <- rep(0L, n)
    for (t in s:n) {
      idx <- (t - s + 1):t
      a <- m1[idx]; b <- m2[idx]
      if (any(is.na(a)) || any(is.na(b))) next
      r_sum[t]   <- cor(a, b)
      r_count[t] <- 1L
    }
    list(r_sum = r_sum, r_count = r_count)
  }

  rep_results <- run_reps(one_rep, n_rep, trend1 = trend1, trend2 = trend2,
                           n = n, h = h, s = s, sigma_sq = sigma_sq)
  R_sum   <- Reduce(`+`, lapply(rep_results, `[[`, "r_sum"))
  R_count <- Reduce(`+`, lapply(rep_results, `[[`, "r_count"))

  R_mean <- ifelse(R_count > 0, R_sum / R_count, NA_real_)
  valid  <- which(!is.na(R_mean) & !is.na(rho_common) & !is.na(rho_actual))

  rmse_common <- sqrt(mean((R_mean[valid] - rho_common[valid])^2))
  rmse_actual <- sqrt(mean((R_mean[valid] - rho_actual[valid])^2))
  mean_gap    <- mean(abs(R_mean[valid] - rho_common[valid]))

  cat(sprintf("%-8.3f  %-12.5f  %-12.5f  %-12.5f\n",
              eps, rmse_common, rmse_actual, mean_gap))
}

cat("\nExpected: RMSE(common) grows with eps; RMSE(actual) stays small.\n")
