# Validate Proposition 1: V approximation
#
# The empirical variance of R_t across replications, scaled by s, should
# match the theoretical V_t from compute_V() at each time point.

devtools::load_all(".")
library(ggplot2)

set.seed(8314)

# ---- Parameters ------------------------------------------------------------

n     <- 2000
h     <- 20
s     <- 200
n_rep <- 2000

# Shared trend
K <- 10
basis <- sapply(1:K, function(k) sin(2 * pi * k * (1:n) / n))
coefs <- 2.5 / (1:K)^1.2
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
sums      <- acov_sums(acov_eta1, acov_eta2)

ma_trend <- as.numeric(stats::filter(trend, rep(1/h, h), sides = 1))

# Theoretical V
tau_sq <- compute_tau_sq(ma_trend, s)
V_th   <- compute_V(tau_sq, sigma1_sq, sigma2_sq,
                     sums$L1, sums$L2, sums$Q1, sums$Q2, sums$Q12)

# ---- Simulate --------------------------------------------------------------

cat("Running", n_rep, "replications (s =", s, ")...\n")

# Evaluate at a grid of time points
eval_pts <- seq(s + h, n, by = 20)

R_mat <- matrix(NA_real_, n_rep, length(eval_pts))

for (r in 1:n_rep) {
  z1 <- arima.sim(model = list(ar = ar1, ma = ma1), n = n, sd = innov_sd1)
  z2 <- arima.sim(model = list(ar = ar2, ma = ma2), n = n, sd = innov_sd2)
  y1 <- trend + z1
  y2 <- trend + z2
  m1 <- as.numeric(stats::filter(y1, rep(1/h, h), sides = 1))
  m2 <- as.numeric(stats::filter(y2, rep(1/h, h), sides = 1))

  for (j in seq_along(eval_pts)) {
    t0  <- eval_pts[j]
    idx <- (t0 - s + 1):t0
    a <- m1[idx]; b <- m2[idx]
    if (any(is.na(a)) || any(is.na(b))) next
    R_mat[r, j] <- cor(a, b)
  }
}

# Empirical s * Var(R_t) at each evaluation point
V_emp <- s * apply(R_mat, 2, var, na.rm = TRUE)
V_theory <- V_th[eval_pts]

valid <- which(!is.na(V_emp) & !is.na(V_theory) & V_theory > 0)

cat(sprintf("\nCor(V_emp, V_theory): %.4f\n",
            cor(V_emp[valid], V_theory[valid])))
cat(sprintf("Median ratio V_emp/V_theory: %.3f\n",
            median(V_emp[valid] / V_theory[valid])))

# ---- Plot ------------------------------------------------------------------

df <- data.frame(t = eval_pts[valid],
                 V_emp = V_emp[valid],
                 V_theory = V_theory[valid])

p <- ggplot(df, aes(x = V_theory, y = V_emp)) +
  geom_point(alpha = 0.4, size = 1.2) +
  geom_abline(slope = 1, intercept = 0, color = "firebrick") +
  labs(x = "Theoretical V", y = "Empirical s·Var(R_t)",
       title = sprintf("V validation (s=%d, %d reps)", s, n_rep)) +
  theme_minimal()

ggsave("dev/numerical-validation/prop1-V.png", p, width = 5, height = 5)
cat("Plot saved to dev/numerical-validation/prop1-V.png\n")
