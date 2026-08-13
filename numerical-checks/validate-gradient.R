# Validate the gradient of g(theta) and the centering argument
#
# 1. Verify analytic gradient nabla g(theta) against numerical finite
#    differences at multiple parameter configurations.
# 2. Verify that the full 5x5 quadratic form V = nabla_g^T Sigma nabla_g
#    equals the reduced 3x3 form from the Prop 1 proof (centering argument).
#
# Both series carry their own signal: means m1, m2, variances tau1^2, tau2^2,
# and cross term A = r tau1 tau2. The centering argument is a statement about
# the mean components of nabla g dropping out, and holds for any signal pair --
# testing it only at m1 = m2 and tau1 = tau2 would leave that untested.

# ---- Helpers ---------------------------------------------------------------

g <- function(th) {
  A <- th[5] - th[1] * th[2]
  B <- th[3] - th[1]^2
  D <- th[4] - th[2]^2
  A / sqrt(B * D)
}

grad_analytic <- function(th) {
  A <- th[5] - th[1] * th[2]
  B <- th[3] - th[1]^2
  D <- th[4] - th[2]^2
  c(
    -th[2] * (B*D)^(-1/2) + th[1] * A * D * (B*D)^(-3/2),
    -th[1] * (B*D)^(-1/2) + th[2] * A * B * (B*D)^(-3/2),
    -0.5 * A * D * (B*D)^(-3/2),
    -0.5 * A * B * (B*D)^(-3/2),
    (B*D)^(-1/2)
  )
}

grad_numeric <- function(f, x, eps = 1e-7) {
  sapply(seq_along(x), function(i) {
    xp <- xm <- x
    xp[i] <- x[i] + eps; xm[i] <- x[i] - eps
    (f(xp) - f(xm)) / (2 * eps)
  })
}

# Full 5x5 Sigma for white Gaussian noise, two distinct signals
build_Sigma_full <- function(m1, m2, tau1_sq, tau2_sq, A,
                             L1, L2, Q1, Q2, Q12) {
  p1 <- m1^2 + tau1_sq
  p2 <- m2^2 + tau2_sq
  S <- matrix(0, 5, 5)
  S[1,1] <- L1; S[2,2] <- L2
  S[1,3] <- S[3,1] <- 2 * m1 * L1
  S[2,4] <- S[4,2] <- 2 * m2 * L2
  S[1,5] <- S[5,1] <- m2 * L1
  S[2,5] <- S[5,2] <- m1 * L2
  S[3,3] <- 2 * Q1 + 4 * p1 * L1
  S[4,4] <- 2 * Q2 + 4 * p2 * L2
  # the cross-pairing: each series' signal power meets the *other* series'
  # long-run noise sum, which is the term the two-tau revision introduced
  S[5,5] <- Q12 + p1 * L2 + p2 * L1
  S[3,5] <- S[5,3] <- 2 * (m1 * m2 + A) * L1
  S[4,5] <- S[5,4] <- 2 * (m1 * m2 + A) * L2
  S
}

# The paper's reduced 3x3 Sigma: the same matrix in centered coordinates, so
# tau_k^2 replaces m_k^2 + tau_k^2 and A replaces m1 m2 + A
build_Sigma_paper <- function(tau1_sq, tau2_sq, A, L1, L2, Q1, Q2, Q12) {
  S <- matrix(0, 3, 3)
  S[1,1] <- 2 * Q1 + 4 * tau1_sq * L1
  S[2,2] <- 2 * Q2 + 4 * tau2_sq * L2
  S[3,3] <- Q12 + tau1_sq * L2 + tau2_sq * L1
  S[1,3] <- S[3,1] <- 2 * A * L1
  S[2,3] <- S[3,2] <- 2 * A * L2
  S
}

# ---- Test configurations ---------------------------------------------------

# r is the signal correlation; r = 1 with tau1 == tau2 is the old shared-trend
# case, kept as the first configuration so the reduction is exercised there too.
configs <- list(
  list(m1 = 0,  m2 = 0,  t1 = 2,   t2 = 2,   r = 1.0, s1sq = 1,   s2sq = 3),
  list(m1 = 5,  m2 = 5,  t1 = 2,   t2 = 0.5, r = 1.0, s1sq = 1,   s2sq = 3),
  list(m1 = 10, m2 = -4, t1 = 0.5, t2 = 3,   r = 0.7, s1sq = 2,   s2sq = 4),
  list(m1 = 3,  m2 = 8,  t1 = 5,   t2 = 1,   r = -0.4, s1sq = 0.5, s2sq = 1.5)
)

# Covariance sums (white noise for simplicity: L = sigma^2, Q = sigma^4)
cat("=== Gradient Verification ===\n\n")

all_pass <- TRUE

for (i in seq_along(configs)) {
  cfg <- configs[[i]]
  m1 <- cfg$m1; m2 <- cfg$m2; s1sq <- cfg$s1sq; s2sq <- cfg$s2sq
  tau1_sq <- cfg$t1; tau2_sq <- cfg$t2
  A <- cfg$r * sqrt(tau1_sq * tau2_sq)

  theta <- c(m1, m2,
             m1^2 + tau1_sq + s1sq,
             m2^2 + tau2_sq + s2sq,
             m1 * m2 + A)

  ga <- grad_analytic(theta)
  gn <- grad_numeric(g, theta)

  grad_ok <- all.equal(ga, gn, tolerance = 1e-5)

  cat(sprintf("Config %d (m=%.0f/%.0f, tau^2=%.1f/%.1f, r=%.1f, sigma^2=%.1f/%.1f):\n",
              i, m1, m2, tau1_sq, tau2_sq, cfg$r, s1sq, s2sq))
  cat(sprintf("  Gradient match: %s\n", if (isTRUE(grad_ok)) "OK" else grad_ok))

  # V comparison: full 5x5 vs paper 3x3
  L1 <- s1sq; L2 <- s2sq; Q1 <- s1sq^2; Q2 <- s2sq^2; Q12 <- s1sq * s2sq

  Sf <- build_Sigma_full(m1, m2, tau1_sq, tau2_sq, A, L1, L2, Q1, Q2, Q12)
  Sp <- build_Sigma_paper(tau1_sq, tau2_sq, A, L1, L2, Q1, Q2, Q12)

  V_full  <- as.numeric(t(ga) %*% Sf %*% ga)
  V_paper <- as.numeric(t(ga[3:5]) %*% Sp %*% ga[3:5])

  V_ok <- all.equal(V_full, V_paper, tolerance = 1e-12)
  cat(sprintf("  V_full = %.8f, V_paper = %.8f, match: %s\n\n",
              V_full, V_paper, if (isTRUE(V_ok)) "OK" else V_ok))

  if (!isTRUE(grad_ok) || !isTRUE(V_ok)) all_pass <- FALSE
}

cat(sprintf("Overall: %s\n", if (all_pass) "ALL PASSED" else "SOME FAILED"))
