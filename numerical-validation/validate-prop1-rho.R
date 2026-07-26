# Validate Proposition 1: rho approximation
#
# Under the common-trend model, the theoretical rho_t from Proposition 1
# should match the mean empirical R_t across replications. We also check
# convergence as s increases.

library(lomad)

set.seed(5102)

# ---- Parameters ------------------------------------------------------------

n     <- 3000
h     <- 20
s_vals <- c(60, 120, 250)
n_rep <- 500

# Shared trend
K <- 12
basis <- sapply(1:K, function(k) sin(2 * pi * k * (1:n) / n))
coefs <- 3 / (1:K)^1.3
trend <- as.numeric(basis %*% coefs)

# ARMA noise
ar1 <- 0.5; ma1 <- 0.2; innov_sd1 <- 0.7
ar2 <- 0.3; ma2 <- -0.1; innov_sd2 <- 0.9

acov_raw1 <- arma_acov(ar1, ma1, innov_sd1^2, lag_max = 200)
acov_raw2 <- arma_acov(ar2, ma2, innov_sd2^2, lag_max = 200)
acov_eta1 <- .ma_filter_acov(acov_raw1, h, lag_max = 200)
acov_eta2 <- .ma_filter_acov(acov_raw2, h, lag_max = 200)
sigma1_sq <- acov_eta1[1]
sigma2_sq <- acov_eta2[1]
acov_eta1 <- acov_eta1[!is.na(acov_eta1)]
acov_eta2 <- acov_eta2[!is.na(acov_eta2)]
ml <- min(length(acov_eta1), length(acov_eta2))
acov_eta1 <- acov_eta1[1:ml]; acov_eta2 <- acov_eta2[1:ml]

ma_trend <- as.numeric(stats::filter(trend, rep(1/h, h), sides = 1))

# ---- Simulate and compare -------------------------------------------------

cat("Running", n_rep, "replications for each s...\n\n")

for (s in s_vals) {
  tau_sq <- compute_tau_sq(ma_trend, s)
  rho_th <- compute_rho(tau_sq, sigma1_sq, sigma2_sq)

  # Accumulate R_t across replications
  R_sum   <- rep(0, n)
  R_count <- rep(0L, n)

  for (r in 1:n_rep) {
    z1 <- arima.sim(model = list(ar = ar1, ma = ma1), n = n, sd = innov_sd1)
    z2 <- arima.sim(model = list(ar = ar2, ma = ma2), n = n, sd = innov_sd2)
    y1 <- trend + z1
    y2 <- trend + z2
    m1 <- as.numeric(stats::filter(y1, rep(1/h, h), sides = 1))
    m2 <- as.numeric(stats::filter(y2, rep(1/h, h), sides = 1))

    for (t in s:n) {
      idx <- (t - s + 1):t
      a <- m1[idx]; b <- m2[idx]
      if (any(is.na(a)) || any(is.na(b))) next
      R_sum[t]   <- R_sum[t] + cor(a, b)
      R_count[t] <- R_count[t] + 1L
    }
  }

  R_mean <- ifelse(R_count > 0, R_sum / R_count, NA_real_)
  valid  <- which(!is.na(R_mean) & !is.na(rho_th))

  corr <- cor(R_mean[valid], rho_th[valid])
  rmse <- sqrt(mean((R_mean[valid] - rho_th[valid])^2))
  bias <- mean(R_mean[valid] - rho_th[valid])

  cat(sprintf("s = %3d:  cor = %.4f  RMSE = %.5f  bias = %.5f\n",
              s, corr, rmse, bias))
}
