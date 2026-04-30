# Validate the Sigma_{(3,4,5)} matrix entries (Prop 1 proof)
#
# Monte Carlo verification that the Sigma entries computed via Wick's theorem
# in centered coordinates match empirical covariances of the moment process.
# Tests with both white noise and ARMA noise.

set.seed(6291)

# ---- Helpers ---------------------------------------------------------------

# Compute empirical Sigma_{(3,4,5)} via MC
mc_sigma_345 <- function(s_fixed, sigma1, sigma2, n_rep,
                         ar1 = numeric(0), ma1 = numeric(0),
                         ar2 = numeric(0), ma2 = numeric(0)) {
  s_len <- length(s_fixed)
  H_bars <- matrix(NA, n_rep, 3)

  for (r in 1:n_rep) {
    if (length(ar1) > 0 || length(ma1) > 0) {
      eta1 <- arima.sim(model = list(ar = ar1, ma = ma1), n = s_len, sd = sigma1)
    } else {
      eta1 <- rnorm(s_len, sd = sigma1)
    }
    if (length(ar2) > 0 || length(ma2) > 0) {
      eta2 <- arima.sim(model = list(ar = ar2, ma = ma2), n = s_len, sd = sigma2)
    } else {
      eta2 <- rnorm(s_len, sd = sigma2)
    }
    Y1 <- s_fixed + eta1
    Y2 <- s_fixed + eta2
    H_bars[r, ] <- c(mean(Y1^2), mean(Y2^2), mean(Y1 * Y2))
  }

  s_len * cov(H_bars)
}

# Analytic Sigma_{(3,4,5)} in centered coordinates (from the proof)
analytic_sigma_345 <- function(tau2, L1, L2, Q1, Q2, Q12) {
  S <- matrix(0, 3, 3)
  S[1,1] <- 2 * Q1 + 4 * tau2 * L1
  S[2,2] <- 2 * Q2 + 4 * tau2 * L2
  S[3,3] <- Q12 + tau2 * (L1 + L2)
  S[1,3] <- S[3,1] <- 2 * tau2 * L1
  S[2,3] <- S[3,2] <- 2 * tau2 * L2
  S
}

# ---- Test 1: White noise ---------------------------------------------------

cat("=== Test 1: White noise, centered trend ===\n\n")

s_len  <- 500
tau2   <- 2
sigma1 <- 1; sigma2 <- sqrt(3)
n_rep  <- 50000

# Centered trend: zero mean, variance = tau2
s_fixed <- seq(-sqrt(3 * tau2), sqrt(3 * tau2), length.out = s_len)
s_fixed <- s_fixed - mean(s_fixed)  # ensure exactly centered
actual_tau2 <- mean(s_fixed^2)
cat(sprintf("Trend: mean = %.4f, tau2 = %.4f\n", mean(s_fixed), actual_tau2))

# White noise: L = sigma^2, Q = sigma^4
L1 <- sigma1^2; Q1 <- sigma1^4
L2 <- sigma2^2; Q2 <- sigma2^4
Q12 <- sigma1^2 * sigma2^2

S_analytic <- analytic_sigma_345(actual_tau2, L1, L2, Q1, Q2, Q12)
S_mc       <- mc_sigma_345(s_fixed, sigma1, sigma2, n_rep)

cat("\nAnalytic Sigma_{(3,4,5)}:\n")
print(round(S_analytic, 3))
cat("\nMC Sigma_{(3,4,5)}:\n")
print(round(S_mc, 3))
cat("\nRelative error (%):\n")
rel_err <- 100 * (S_mc - S_analytic) / S_analytic
rel_err[S_analytic == 0] <- NA
print(round(rel_err, 2))

# ---- Test 2: ARMA noise ---------------------------------------------------

cat("\n=== Test 2: AR(1) noise, centered trend ===\n\n")

devtools::load_all(".")

ar_coef <- 0.6
innov1 <- 0.8; innov2 <- 1.0
n_rep2 <- 30000

# Compute analytic covariance sums for AR(1) noise
acov1 <- arma_acov(ar_coef, numeric(0), innov1^2, lag_max = 200)
acov2 <- arma_acov(ar_coef, numeric(0), innov2^2, lag_max = 200)
csums  <- acov_sums(acov1, acov2)

S_analytic2 <- analytic_sigma_345(actual_tau2,
                                   csums$L1, csums$L2,
                                   csums$Q1, csums$Q2, csums$Q12)
S_mc2       <- mc_sigma_345(s_fixed, innov1, innov2, n_rep2,
                             ar1 = ar_coef, ar2 = ar_coef)

cat("Analytic Sigma_{(3,4,5)}:\n")
print(round(S_analytic2, 3))
cat("\nMC Sigma_{(3,4,5)}:\n")
print(round(S_mc2, 3))
cat("\nRelative error (%):\n")
rel_err2 <- 100 * (S_mc2 - S_analytic2) / S_analytic2
rel_err2[S_analytic2 == 0] <- NA
print(round(rel_err2, 2))
