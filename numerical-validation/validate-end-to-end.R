# End-to-end validation: full pipeline vs oracle
#
# Generates data under known parameters, runs the full lomad_fit() + lomad_test()
# pipeline (which estimates trends, noise, etc.), and compares:
#   1. Estimated rho_hat vs oracle rho (from known parameters)
#   2. Estimated V_hat vs oracle V
#   3. CLT coverage of the estimated test statistic
#   4. Type I error rate under H0 (common trend, d = 0)
#
# This checks whether estimation error degrades the CLT at practical sample sizes.

library(lomad)
library(ggplot2)

set.seed(7742)

# ---- Parameters ------------------------------------------------------------

n     <- 2000
h     <- 20
s     <- 150
n_rep <- 500
alpha <- 0.05

# AR(1) noise parameters.
# NOTE: The pipeline estimates AR(1) noise. If the true noise is ARMA(1,1),
# the AR(1) misspecification can severely overestimate the noise variance,
# clamping tau_sq to 0 and making rho_hat = 0 everywhere. Use AR(1) noise
# here for a fair test of the CLT machinery; ARMA misspecification is a
# separate issue for the estimation pipeline.
ar1 <- 0.5; ma1 <- numeric(0); innov1 <- 0.15
ar2 <- 0.3; ma2 <- numeric(0); innov2 <- 0.20

# Shared trend (H0: d = 0)
K <- 10
basis <- sapply(1:K, function(k) sin(2 * pi * k * (1:n) / n))
coefs <- 2 / (1:K)^1.2
trend <- as.numeric(basis %*% coefs)

# ---- Oracle quantities -----------------------------------------------------

acov_raw1  <- arma_acov(ar1, ma1, innov1^2, lag_max = h + 100)
acov_raw2  <- arma_acov(ar2, ma2, innov2^2, lag_max = h + 100)
# Truncate trailing NAs from filter edge effects
acov_filt1 <- .ma_filter_acov(acov_raw1, h, lag_max = 100)
acov_filt2 <- .ma_filter_acov(acov_raw2, h, lag_max = 100)
acov_filt1 <- acov_filt1[!is.na(acov_filt1)]
acov_filt2 <- acov_filt2[!is.na(acov_filt2)]
ml <- min(length(acov_filt1), length(acov_filt2))
acov_filt1 <- acov_filt1[1:ml]; acov_filt2 <- acov_filt2[1:ml]

sigma1_oracle <- acov_filt1[1]
sigma2_oracle <- acov_filt2[1]
sums_oracle   <- acov_sums(acov_filt1, acov_filt2)

ma_trend <- as.numeric(stats::filter(trend, rep(1/h, h), sides = 1))
# Subtract noise bias from tau_sq (as the pipeline does)
noise_bias    <- (sigma1_oracle + sigma2_oracle) / 4
tau_sq_oracle <- pmax(0, compute_tau_sq(ma_trend, s) - noise_bias)
rho_oracle    <- compute_rho(tau_sq_oracle, sigma1_oracle, sigma2_oracle)
V_oracle      <- compute_V(tau_sq_oracle, sigma1_oracle, sigma2_oracle,
                            sums_oracle$L1, sums_oracle$L2,
                            sums_oracle$Q1, sums_oracle$Q2, sums_oracle$Q12)

# ---- Evaluation points (avoid edges) --------------------------------------

eval_pts <- c(400, 700, 1000, 1300, 1600)

# ---- Run replications ------------------------------------------------------

cat(sprintf("Running %d replications (n=%d, h=%d, s=%d)...\n", n_rep, n, h, s))

# Store results at eval points
R_oracle_mat <- matrix(NA_real_, n_rep, length(eval_pts))   # R_t from data
rho_est_mat  <- matrix(NA_real_, n_rep, length(eval_pts))   # estimated rho
V_est_mat    <- matrix(NA_real_, n_rep, length(eval_pts))   # estimated V
Z_est_mat    <- matrix(NA_real_, n_rep, length(eval_pts))   # pipeline Z stat
rej_mat      <- matrix(NA,       n_rep, length(eval_pts))   # pipeline rejection

for (r in 1:n_rep) {
  z1 <- arima.sim(model = list(ar = ar1, ma = ma1), n = n, sd = innov1)
  z2 <- arima.sim(model = list(ar = ar2, ma = ma2), n = n, sd = innov2)
  y1 <- trend + z1
  y2 <- trend + z2

  fit <- suppressMessages(lomad_fit(y1, y2, method = "clt", h = h, s = s))
  tst <- suppressMessages(lomad_test(fit, alpha = alpha))

  for (j in seq_along(eval_pts)) {
    t0 <- eval_pts[j]
    R_oracle_mat[r, j] <- fit$R[t0]
    rho_est_mat[r, j]  <- fit$rho[t0]
    V_est_mat[r, j]    <- fit$V[t0]
    Z_est_mat[r, j]    <- tst$Z[t0]
    rej_mat[r, j]       <- tst$rejected[t0]
  }

  if (r %% 100 == 0) cat(sprintf("  %d / %d\n", r, n_rep))
}

# ---- Summary ---------------------------------------------------------------

cat("\n=== End-to-end validation ===\n\n")

cat("1. Estimated rho vs oracle rho (mean across reps):\n")
for (j in seq_along(eval_pts)) {
  t0 <- eval_pts[j]
  cat(sprintf("  t=%d: rho_oracle=%.4f  mean(rho_hat)=%.4f  sd(rho_hat)=%.4f\n",
              t0, rho_oracle[t0],
              mean(rho_est_mat[, j], na.rm = TRUE),
              sd(rho_est_mat[, j], na.rm = TRUE)))
}

cat("\n2. Estimated V vs oracle V (mean across reps):\n")
for (j in seq_along(eval_pts)) {
  t0 <- eval_pts[j]
  V_or <- V_oracle[t0]
  if (is.na(V_or)) next
  V_hat_mean <- mean(V_est_mat[, j], na.rm = TRUE)
  cat(sprintf("  t=%d: V_oracle=%.4f  mean(V_hat)=%.4f  ratio=%.3f\n",
              t0, V_or, V_hat_mean, V_hat_mean / V_or))
}

cat("\n3. CLT coverage (using pipeline Z, against N(0,1)):\n")
for (j in seq_along(eval_pts)) {
  t0 <- eval_pts[j]
  # Oracle standardization: sqrt(s)*(R - rho_oracle)/sqrt(V_oracle)
  R_j <- R_oracle_mat[, j]; V_or <- V_oracle[t0]; rho_or <- rho_oracle[t0]
  if (is.na(V_or) || V_or <= 0 || is.na(rho_or)) next
  Z_oracle <- sqrt(s) * (R_j - rho_or) / sqrt(V_or)
  Z_oracle <- Z_oracle[!is.na(Z_oracle)]

  # Pipeline Z (uses estimated rho, V)
  Z_pipe <- Z_est_mat[, j]; Z_pipe <- Z_pipe[!is.na(Z_pipe)]

  cov95_oracle <- mean(abs(Z_oracle) < qnorm(0.975))
  cov95_pipe   <- mean(abs(Z_pipe) < qnorm(0.975))

  cat(sprintf("  t=%d: cov95_oracle=%.3f  cov95_pipeline=%.3f  sd_Z_oracle=%.3f  sd_Z_pipe=%.3f\n",
              t0, cov95_oracle, cov95_pipe,
              sd(Z_oracle), sd(Z_pipe)))
}

cat("\n4. Type I error rate under H0 (should be ≈ 0 with BY-FDR correction):\n")
for (j in seq_along(eval_pts)) {
  t0 <- eval_pts[j]
  rej_j <- rej_mat[, j]; rej_j <- rej_j[!is.na(rej_j)]
  cat(sprintf("  t=%d: rejection rate = %.4f (%d / %d)\n",
              t0, mean(rej_j), sum(rej_j), length(rej_j)))
}

# ---- QQ plot: oracle vs pipeline standardization ---------------------------

j_mid <- which(eval_pts == 1000)
R_j <- R_oracle_mat[, j_mid]
V_or <- V_oracle[eval_pts[j_mid]]; rho_or <- rho_oracle[eval_pts[j_mid]]

Z_or <- sqrt(s) * (R_j - rho_or) / sqrt(V_or); Z_or <- Z_or[!is.na(Z_or)]
Z_pi <- Z_est_mat[, j_mid]; Z_pi <- Z_pi[!is.na(Z_pi)]

qq_df <- data.frame(
  z = c(sort(Z_or), sort(Z_pi)),
  theoretical = rep(qnorm(ppoints(length(Z_or))), 2),
  type = rep(c("Oracle", "Pipeline"), each = length(Z_or))
)

ggplot(qq_df, aes(x = theoretical, y = z, color = type)) +
  geom_point(alpha = 0.3, size = 0.8) +
  geom_abline(slope = 1, intercept = 0, color = "grey40") +
  labs(x = "Theoretical N(0,1)", y = "Empirical Z",
       title = sprintf("QQ: oracle vs pipeline (t=%d, s=%d)", eval_pts[j_mid], s),
       color = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("numerical-validation/end-to-end-qq.png", width = 6, height = 5)
cat("\nQQ plot saved to numerical-validation/end-to-end-qq.png\n")
