# Validate the CLT for R_t (Theorems 1 & 2)
#
# Under the common-trend model with known ARMA noise, the standardized
# statistic sqrt(s)*(R_t - rho_t)/sqrt(V_t) should be approximately N(0,1).
# We verify this at a few selected time points via repeated simulation.

devtools::load_all(".")
library(ggplot2)

set.seed(3847)

# ---- Parameters ------------------------------------------------------------

n    <- 2000
h    <- 20
s_vals <- c(80, 150, 300)  # window sizes to test
n_rep <- 2000

# Shared trend: deterministic Fourier series
K <- 10
basis <- sapply(1:K, function(k) sin(2 * pi * k * (1:n) / n))
coefs <- 2 / (1:K)^1.2
trend <- as.numeric(basis %*% coefs)

# ARMA noise parameters (independent for series 1 and 2)
ar1 <- 0.6; ma1 <- 0.3; innov_sd1 <- 0.8
ar2 <- 0.4; ma2 <- -0.2; innov_sd2 <- 1.0

# Precompute noise autocovariances and covariance sums
acov_raw1 <- arma_acov(ar1, ma1, innov_sd1^2, lag_max = 200)
acov_raw2 <- arma_acov(ar2, ma2, innov_sd2^2, lag_max = 200)
acov_eta1 <- .ma_filter_acov(acov_raw1, h, lag_max = 200)
acov_eta2 <- .ma_filter_acov(acov_raw2, h, lag_max = 200)
sigma1_sq <- acov_eta1[1]
sigma2_sq <- acov_eta2[1]
# Truncate trailing NAs from filter edge effects
acov_eta1 <- acov_eta1[!is.na(acov_eta1)]
acov_eta2 <- acov_eta2[!is.na(acov_eta2)]
ml <- min(length(acov_eta1), length(acov_eta2))
acov_eta1 <- acov_eta1[1:ml]; acov_eta2 <- acov_eta2[1:ml]
sums      <- acov_sums(acov_eta1, acov_eta2)

# MA-smoothed trend
ma_trend <- as.numeric(stats::filter(trend, rep(1/h, h), sides = 1))

# Select evaluation time points (avoiding edges)
eval_points <- c(500, 1000, 1500)

# ---- Simulate --------------------------------------------------------------

cat("Running", n_rep, "replications...\n")

results <- list()

for (s in s_vals) {
  # Theoretical quantities at eval points
  tau_sq <- compute_tau_sq(ma_trend, s)
  rho_th <- compute_rho(tau_sq, sigma1_sq, sigma2_sq)
  V_th   <- compute_V(tau_sq, sigma1_sq, sigma2_sq,
                       sums$L1, sums$L2, sums$Q1, sums$Q2, sums$Q12)

  # Collect R_t at eval points across replications
  R_mat <- matrix(NA_real_, n_rep, length(eval_points))

  for (r in 1:n_rep) {
    z1 <- arima.sim(model = list(ar = ar1, ma = ma1), n = n, sd = innov_sd1)
    z2 <- arima.sim(model = list(ar = ar2, ma = ma2), n = n, sd = innov_sd2)
    y1 <- trend + z1
    y2 <- trend + z2
    ma1_s <- as.numeric(stats::filter(y1, rep(1/h, h), sides = 1))
    ma2_s <- as.numeric(stats::filter(y2, rep(1/h, h), sides = 1))

    for (j in seq_along(eval_points)) {
      t0  <- eval_points[j]
      idx <- (t0 - s + 1):t0
      m1  <- ma1_s[idx]; m2 <- ma2_s[idx]
      if (any(is.na(m1)) || any(is.na(m2))) next
      R_mat[r, j] <- cor(m1, m2)
    }
  }

  for (j in seq_along(eval_points)) {
    t0   <- eval_points[j]
    R_j  <- R_mat[, j]
    R_j  <- R_j[!is.na(R_j)]
    rho_j <- rho_th[t0]
    V_j   <- V_th[t0]
    if (is.na(rho_j) || is.na(V_j) || V_j <= 0) next

    Z_j <- sqrt(s) * (R_j - rho_j) / sqrt(V_j)

    sw <- shapiro.test(sample(Z_j, min(length(Z_j), 5000)))
    cov90 <- mean(abs(Z_j) < qnorm(0.95))
    cov95 <- mean(abs(Z_j) < qnorm(0.975))

    results[[length(results) + 1]] <- data.frame(
      s = s, t = t0, n_obs = length(Z_j),
      mean_Z = mean(Z_j), sd_Z = sd(Z_j),
      shapiro_p = sw$p.value,
      cov90 = cov90, cov95 = cov95
    )
  }
}

res_df <- do.call(rbind, results)
cat("\n=== CLT Validation Summary ===\n")
cat("Target: mean_Z ≈ 0, sd_Z ≈ 1, cov90 ≈ 0.90, cov95 ≈ 0.95\n\n")
print(res_df, digits = 3, row.names = FALSE)

# ---- QQ plot at largest s ---------------------------------------------------

s_qq <- max(s_vals)
tau_sq_qq <- compute_tau_sq(ma_trend, s_qq)
rho_qq    <- compute_rho(tau_sq_qq, sigma1_sq, sigma2_sq)
V_qq      <- compute_V(tau_sq_qq, sigma1_sq, sigma2_sq,
                        sums$L1, sums$L2, sums$Q1, sums$Q2, sums$Q12)

# Reuse last simulation batch (s = s_qq)
t0 <- eval_points[2]
R_j <- R_mat[, 2]; R_j <- R_j[!is.na(R_j)]
Z_qq <- sqrt(s_qq) * (R_j - rho_qq[t0]) / sqrt(V_qq[t0])

qq_df <- data.frame(z = sort(Z_qq))
qq_df$theoretical <- qnorm(ppoints(nrow(qq_df)))

p_qq <- ggplot(qq_df, aes(x = theoretical, y = z)) +
  geom_point(alpha = 0.3, size = 0.8) +
  geom_abline(slope = 1, intercept = 0, color = "firebrick") +
  labs(x = "Theoretical N(0,1)", y = "Empirical",
       title = sprintf("QQ plot: standardized R_t (s=%d, t=%d)", s_qq, t0)) +
  theme_minimal()

ggsave("dev/numerical-validation/clt-qq.png", p_qq, width = 5, height = 5)
cat("\nQQ plot saved to dev/numerical-validation/clt-qq.png\n")
