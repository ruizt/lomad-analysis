# Can we do better than difference-based AR(1) now that leakage is gone?
#
# The test does not consume phi. It consumes the covariance sums of the
# MA-filtered noise: L_k = sum gamma_k(l), Q_k = sum gamma_k(l)^2, and Q_12.
# V is dominated by L, so accuracy is scored on L, not on phi.
#
# Three estimators, all applied directly to the observed series (no trend
# estimate), under correctly specified and misspecified noise:
#
#   ar1   difference-based AR(1)          -- what the package does now
#   npar  nonparametric ACVF from the variogram, gamma(l) = V(Lref) - V(l)
#   arp   AR(p) by Yule-Walker on that nonparametric ACVF, order by AIC
#
# Usage (from the repo root):
#   Rscript numerical-checks/assess-noise-estimators.R

suppressPackageStartupMessages({ library(lomad); library(dplyr) })

N <- 600L; H <- 5L; NREP <- 300L; LAG_MAX <- 60L; LREF <- 40L

# Nonparametric ACVF from the sample variogram of the raw series. V(l) =
# gamma(0) - gamma(l) for stationary noise; a trend inflates V(l) smoothly in
# l, so LREF is taken past the noise correlation length but short enough that
# the trend contribution is still small.
npar_acov <- function(y, lag_max = LAG_MAX, lref = LREF) {
  n <- length(y)
  V <- vapply(seq_len(lref), function(l) mean((y[(l + 1):n] - y[1:(n - l)])^2) / 2,
              numeric(1))
  g0 <- mean(V[max(1, lref - 4):lref])          # plateau = gamma(0)
  g  <- pmax(0, g0 - V[seq_len(min(lag_max, lref))])
  out <- numeric(lag_max + 1L)
  out[1] <- g0
  k <- min(length(g), lag_max)
  out[2:(k + 1L)] <- g[seq_len(k)]      # zero-padded past lref
  out
}

arp_from_acov <- function(g, p_max = 6L) {
  best <- NULL
  for (p in seq_len(p_max)) {
    fit <- tryCatch(solve(toeplitz(g[1:p]), g[2:(p + 1)]), error = function(e) NULL)
    if (is.null(fit)) next
    s2  <- g[1] - sum(fit * g[2:(p + 1)])
    if (!is.finite(s2) || s2 <= 0) next
    aic <- length(g) * log(s2) + 2 * p
    if (is.null(best) || aic < best$aic)
      best <- list(ar = fit, sigma2 = s2, p = p, aic = aic)
  }
  best
}

# L = sum over all integer lags = gamma(0) + 2 * sum_{l >= 1}
L_of <- function(acov_raw) {
  g <- lomad:::.ma_filter_acov(acov_raw, H, LAG_MAX)
  g[1] + 2 * sum(g[-1])
}

truths <- list(
  "AR(1) phi=0.5"      = list(ar = 0.5,          ma = numeric(0)),
  "AR(1) phi=0.8"      = list(ar = 0.8,          ma = numeric(0)),
  "AR(2) .6,-.3"       = list(ar = c(0.6, -0.3), ma = numeric(0)),
  "ARMA(1,1) .7,.4"    = list(ar = 0.7,          ma = 0.4),
  "MA(1) .6"           = list(ar = numeric(0),   ma = 0.6)
)

cat("\nscored on L = sum_l gamma_eta(l), the dominant term in V\n")
cat("relative error in L_hat / L_true, 300 reps, trend present (d = 0.75)\n\n")

out <- list()
for (nm in names(truths)) {
  sp <- truths[[nm]]
  est <- t(vapply(seq_len(NREP), function(r) {
    set.seed(4000L + r)
    tr  <- sim_trends(n = N, d = 0.75, method = "smooth", bw = 50,
                      coupling = 0.8, seed = 4000L + r)
    sim <- sim_noise_pair(tr, h = H, lambda_target = 1.0,
                          ar.coefs = sp$ar, ma.coefs = sp$ma, seed = 5000L + r)
    y <- sim$y1
    # sim_noise_pair() rescales the noise to hit lambda_target, so the truth
    # has to be evaluated at the realized innovation variance, not at 1.
    z  <- sim$y1 - sim$x1
    s2 <- if (length(sp$ar)) var(stats::filter(z, c(1, -sp$ar), sides = 1,
                                               method = "convolution"),
                                na.rm = TRUE) else var(z)
    if (length(sp$ma)) s2 <- s2 / (1 + sum(sp$ma^2))
    Lt <- L_of(arma_acov(sp$ar, sp$ma, s2, lag_max = LAG_MAX + H))
    # 1. difference-based AR(1)
    n1 <- suppressWarnings(lomad:::.variogram_ar1(y))
    L1 <- L_of(arma_acov(n1$ar, numeric(0), n1$sigma2, lag_max = LAG_MAX + H))
    # 2. nonparametric ACVF
    gn <- npar_acov(y)
    L2 <- L_of(gn)
    # 3. AR(p) on the nonparametric ACVF
    ap <- arp_from_acov(gn)
    L3 <- if (is.null(ap)) NA_real_
          else L_of(arma_acov(ap$ar, numeric(0), ap$sigma2, lag_max = LAG_MAX + H))
    c(L1 / Lt, L2 / Lt, L3 / Lt, if (is.null(ap)) NA_real_ else ap$p)
  }, numeric(4)))

  b <- 100 * (colMeans(est[, 1:3], na.rm = TRUE) - 1)
  s <- 100 * apply(est[, 1:3], 2, sd, na.rm = TRUE)
  cat(sprintf("%-18s | ar1 %+7.1f%% (sd %5.1f%%) | npar %+7.1f%% (sd %5.1f%%) | arp %+7.1f%% (p=%.1f)\n",
              nm, b[1], s[1], b[2], s[2], b[3], mean(est[, 4], na.rm = TRUE)))
  out[[nm]] <- list(bias = b, sd = s)
}

saveRDS(out, "numerical-checks/_noise-estimators.rds")
cat("\nwrote numerical-checks/_noise-estimators.rds\n")
