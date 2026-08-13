# Monte Carlo validation of Proposition 1
#
# Checks every quantity the proposition asserts, under a null of local affine
# similarity with two distinct signals:
#
#   1. Sigma_{(3,4,5)} -- in particular the (5,5) cross-pairing
#      tau_1^2 L_2 + tau_2^2 L_1, and the (3,4) zero.
#   2. rho = r tau_1 tau_2 / sqrt(BD), including r != 1.
#   3. V under H0, via the MC variance of the local sample correlation.
#   4. The attenuation identity rho = sqrt(1 - delta^2) rho^(0).
#   5. The symmetric case tau1 == tau2, r == 1, where the cross-pairing
#      collapses to tau^2 (L_1 + L_2) and rho to tau^2 / sqrt(BD).
#
# Usage (from the repo root):
#   Rscript numerical-checks/validate-prop1.R

suppressPackageStartupMessages(library(lomad))
set.seed(6291)

TOL_PCT <- 2.0   # MC tolerance, percent relative error

pass <- function(ok, label, detail = "") {
  cat(sprintf("  [%s] %s%s\n", if (ok) "PASS" else "FAIL", label,
              if (nzchar(detail)) paste0("  ", detail) else ""))
  if (!ok) assign("ANY_FAIL", TRUE, envir = globalenv())
}
ANY_FAIL <- FALSE

# Two centred, mean-square-orthonormal signal patterns, so that
#   mean(s1^2) = tau1^2, mean(s2^2) = tau2^2, mean(s1*s2) = r*tau1*tau2.
make_signals <- function(n, tau1, tau2, r) {
  t  <- seq_len(n)
  e1 <- sin(2 * pi * t / n); e2 <- cos(2 * pi * t / n)
  e1 <- e1 - mean(e1); e2 <- e2 - mean(e2)
  e1 <- e1 / sqrt(mean(e1^2)); e2 <- e2 / sqrt(mean(e2^2))
  # Gram-Schmidt guards against residual correlation from the centring
  e2 <- e2 - mean(e1 * e2) * e1
  e2 <- e2 / sqrt(mean(e2^2))
  list(s1 = tau1 * e1,
       s2 = tau2 * (r * e1 + sqrt(1 - r^2) * e2))
}

sim_noise_acov <- function(n, ar, ma, sigma2) {
  if (length(ar) == 0 && length(ma) == 0)
    return(stats::rnorm(n, sd = sqrt(sigma2)))
  as.numeric(stats::arima.sim(list(ar = ar, ma = ma), n = n,
                              sd = sqrt(sigma2)))
}


# ============================================================================
# 1. Sigma_{(3,4,5)} with tau1 != tau2 and r != 1
# ============================================================================

cat("\n=== 1. Sigma_(3,4,5), two distinct signals ===\n\n")

n_win <- 400L; n_rep <- 60000L
tau1 <- sqrt(2.0); tau2 <- sqrt(0.5); r_sig <- 0.6

# Deliberately asymmetric noise so L1 != L2 and the cross-pairing is visible
sp1 <- list(ar = 0.5,  ma = numeric(0), sigma2 = 1.0)
sp2 <- list(ar = -0.2, ma = numeric(0), sigma2 = 0.4)

g1 <- arma_acov(sp1$ar, sp1$ma, sp1$sigma2, lag_max = 200L)
g2 <- arma_acov(sp2$ar, sp2$ma, sp2$sigma2, lag_max = 200L)
sums <- acov_sums(g1, g2)
L1 <- sums$L1; L2 <- sums$L2; Q1 <- sums$Q1; Q2 <- sums$Q2; Q12 <- sums$Q12
cat(sprintf("L1 = %.4f, L2 = %.4f  (ratio %.2f)\n", L1, L2, L1 / L2))

sig <- make_signals(n_win, tau1, tau2, r_sig)
t1_sq <- mean(sig$s1^2); t2_sq <- mean(sig$s2^2)
A     <- mean(sig$s1 * sig$s2)
cat(sprintf("realized tau1^2 = %.4f, tau2^2 = %.4f, A = r*tau1*tau2 = %.4f\n",
            t1_sq, t2_sq, A))

Hb <- matrix(NA_real_, n_rep, 3)
for (i in seq_len(n_rep)) {
  Y1 <- sig$s1 + sim_noise_acov(n_win, sp1$ar, sp1$ma, sp1$sigma2)
  Y2 <- sig$s2 + sim_noise_acov(n_win, sp2$ar, sp2$ma, sp2$sigma2)
  Hb[i, ] <- c(mean(Y1^2), mean(Y2^2), mean(Y1 * Y2))
}
S_mc <- n_win * cov(Hb)

S_an <- matrix(0, 3, 3)
S_an[1, 1] <- 2 * Q1 + 4 * t1_sq * L1
S_an[2, 2] <- 2 * Q2 + 4 * t2_sq * L2
S_an[3, 3] <- Q12 + t1_sq * L2 + t2_sq * L1      # cross-paired
S_an[1, 3] <- S_an[3, 1] <- 2 * A * L1
S_an[2, 3] <- S_an[3, 2] <- 2 * A * L2

cat("\nanalytic:\n"); print(round(S_an, 4))
cat("MC:\n");        print(round(S_mc, 4))

rel <- 100 * abs(S_mc - S_an) / pmax(abs(S_an), 1e-8)
for (ij in list(c(1,1), c(2,2), c(3,3), c(1,3), c(2,3))) {
  i <- ij[1]; j <- ij[2]
  pass(rel[i, j] < TOL_PCT, sprintf("Sigma[%d,%d]", i + 2L, j + 2L),
       sprintf("analytic %.4f, MC %.4f, rel err %.2f%%",
               S_an[i, j], S_mc[i, j], rel[i, j]))
}
# (3,4) should be zero: it follows from eta_1 independent of eta_2, not from
# any shared-trend assumption, so it must hold at tau1 != tau2 as well.
se34 <- sd(Hb[, 1] * Hb[, 2]) / sqrt(n_rep)
pass(abs(S_mc[1, 2]) < 4 * n_win * se34 || abs(S_mc[1, 2]) < 0.02 * max(abs(S_an)),
     "Sigma[3,4] == 0", sprintf("MC %.5f", S_mc[1, 2]))

# The discriminating check: the un-crossed alternative must be rejected.
S55_wrong <- Q12 + t1_sq * L1 + t2_sq * L2
cat(sprintf("\n  cross-paired  (5,5) = %.4f   <- Proposition 1\n", S_an[3, 3]))
cat(sprintf("  un-crossed    (5,5) = %.4f   <- the transcription error\n",
            S55_wrong))
cat(sprintf("  MC            (5,5) = %.4f\n", S_mc[3, 3]))
pass(abs(S_mc[3, 3] - S_an[3, 3]) < abs(S_mc[3, 3] - S55_wrong),
     "MC favours the cross-paired (5,5)")


# ============================================================================
# 2 & 4. rho with r != 1, and the attenuation identity
# ============================================================================

cat("\n=== 2/4. rho and the attenuation identity ===\n\n")

n_rep2 <- 40000L
for (r_val in c(1, 0.8, 0.5, 0.0)) {
  sg <- make_signals(n_win, tau1, tau2, r_val)
  Rs <- numeric(n_rep2)
  for (i in seq_len(n_rep2)) {
    Y1 <- sg$s1 + sim_noise_acov(n_win, sp1$ar, sp1$ma, sp1$sigma2)
    Y2 <- sg$s2 + sim_noise_acov(n_win, sp2$ar, sp2$ma, sp2$sigma2)
    Rs[i] <- cor(Y1, Y2)
  }
  rho_an <- compute_rho(mean(sg$s1^2), mean(sg$s2^2), g1[1], g2[1], r = r_val)
  rho_0  <- compute_rho(mean(sg$s1^2), mean(sg$s2^2), g1[1], g2[1], r = 1)
  delta  <- sqrt(1 - r_val^2)
  cat(sprintf("r = %.2f: analytic %.5f, MC %.5f | rho0 = %.5f, ",
              r_val, rho_an, mean(Rs), rho_0))
  pass(abs(mean(Rs) - rho_an) < 0.004,
       sprintf("rho at r = %.2f", r_val))
  pass(abs(rho_an - sqrt(1 - delta^2) * rho_0) < 1e-12,
       sprintf("  attenuation identity at r = %.2f", r_val))
}


# ============================================================================
# 3. V under H0
# ============================================================================

cat("\n=== 3. V under H0 (r = 1) ===\n\n")

for (cfg in list(c(2.0, 0.5), c(1.0, 1.0), c(0.3, 3.0))) {
  tt1 <- sqrt(cfg[1]); tt2 <- sqrt(cfg[2])
  sg  <- make_signals(n_win, tt1, tt2, 1)
  Rs  <- numeric(n_rep2)
  for (i in seq_len(n_rep2)) {
    Y1 <- sg$s1 + sim_noise_acov(n_win, sp1$ar, sp1$ma, sp1$sigma2)
    Y2 <- sg$s2 + sim_noise_acov(n_win, sp2$ar, sp2$ma, sp2$sigma2)
    Rs[i] <- cor(Y1, Y2)
  }
  V_an <- compute_V(mean(sg$s1^2), mean(sg$s2^2), g1[1], g2[1],
                    L1, L2, Q1, Q2, Q12)
  V_mc <- n_win * var(Rs)
  rel  <- 100 * abs(V_mc - V_an) / V_an
  pass(rel < 6, sprintf("V at tau1^2 = %.2f, tau2^2 = %.2f", cfg[1], cfg[2]),
       sprintf("analytic %.5f, MC %.5f, rel err %.1f%%", V_an, V_mc, rel))
}


# ============================================================================
# 5. Back-compatibility with the single-tau formulae
# ============================================================================

cat("\n=== 5. tau1 == tau2 reproduces the pre-0.1.0 formulae ===\n\n")

old_rho <- function(tau_sq, s1, s2) tau_sq / sqrt((tau_sq + s1) * (tau_sq + s2))
old_V <- function(tau_sq, s1, s2, L1, L2, Q1, Q2, Q12) {
  B <- tau_sq + s1; D <- tau_sq + s2
  Q12/(B*D) + tau_sq*s1^2*L1/(B^3*D) + tau_sq*s2^2*L2/(B*D^3) +
    tau_sq^2*Q1/(2*B^3*D) + tau_sq^2*Q2/(2*B*D^3)
}
mx_r <- 0; mx_v <- 0
for (tsq in c(0, 0.01, 0.5, 2, 50)) {
  mx_r <- max(mx_r, abs(compute_rho(tsq, tsq, g1[1], g2[1]) -
                        old_rho(tsq, g1[1], g2[1])))
  mx_v <- max(mx_v, abs(compute_V(tsq, tsq, g1[1], g2[1], L1, L2, Q1, Q2, Q12) -
                        old_V(tsq, g1[1], g2[1], L1, L2, Q1, Q2, Q12)))
}
pass(mx_r < 1e-14, "rho matches old formula", sprintf("max abs diff %.2e", mx_r))
pass(mx_v < 1e-14, "V matches old formula",   sprintf("max abs diff %.2e", mx_v))

cat("\n", if (isTRUE(ANY_FAIL)) "SOME CHECKS FAILED\n" else "all checks passed\n",
    sep = "")
