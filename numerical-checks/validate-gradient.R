# Validate the gradient of g(theta) and the centering argument
#
# 1. Verify analytic gradient nabla g(theta) against numerical finite
#    differences at multiple parameter configurations.
# 2. Verify that the full 5x5 quadratic form V = nabla_g^T Sigma nabla_g
#    equals the reduced 3x3 form from the Prop 1 proof (centering argument).

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

# Build the full correct 5x5 Sigma for white Gaussian noise
# under the centered common-trend model
build_Sigma_full <- function(m, tau2, sigma1_sq, sigma2_sq,
                             L1, L2, Q1, Q2, Q12) {
  m2pt <- m^2 + tau2
  S <- matrix(0, 5, 5)
  S[1,1] <- L1; S[2,2] <- L2
  S[1,3] <- S[3,1] <- 2 * m * L1
  S[2,4] <- S[4,2] <- 2 * m * L2
  S[1,5] <- S[5,1] <- m * L1
  S[2,5] <- S[5,2] <- m * L2
  S[3,3] <- 2 * Q1 + 4 * m2pt * L1
  S[4,4] <- 2 * Q2 + 4 * m2pt * L2
  S[5,5] <- Q12 + m2pt * (L1 + L2)
  S[3,5] <- S[5,3] <- 2 * m2pt * L1
  S[4,5] <- S[5,4] <- 2 * m2pt * L2
  S
}

# Build the paper's reduced 3x3 Sigma (centered, tau2 not m2+tau2)
build_Sigma_paper <- function(tau2, L1, L2, Q1, Q2, Q12) {
  S <- matrix(0, 3, 3)
  S[1,1] <- 2 * Q1 + 4 * tau2 * L1
  S[2,2] <- 2 * Q2 + 4 * tau2 * L2
  S[3,3] <- Q12 + tau2 * (L1 + L2)
  S[1,3] <- S[3,1] <- 2 * tau2 * L1
  S[2,3] <- S[3,2] <- 2 * tau2 * L2
  S
}

# ---- Test configurations ---------------------------------------------------

configs <- list(
  list(m = 0,  tau2 = 2,   s1sq = 1,   s2sq = 3),
  list(m = 5,  tau2 = 2,   s1sq = 1,   s2sq = 3),
  list(m = 10, tau2 = 0.5, s1sq = 2,   s2sq = 4),
  list(m = 3,  tau2 = 5,   s1sq = 0.5, s2sq = 1.5)
)

# Covariance sums (white noise for simplicity: L = sigma^2, Q = sigma^4)
cat("=== Gradient Verification ===\n\n")

all_pass <- TRUE

for (i in seq_along(configs)) {
  cfg <- configs[[i]]
  m <- cfg$m; tau2 <- cfg$tau2; s1sq <- cfg$s1sq; s2sq <- cfg$s2sq

  theta <- c(m, m, m^2 + tau2 + s1sq, m^2 + tau2 + s2sq, m^2 + tau2)

  ga <- grad_analytic(theta)
  gn <- grad_numeric(g, theta)

  grad_ok <- all.equal(ga, gn, tolerance = 1e-5)

  cat(sprintf("Config %d (m=%.0f, tau2=%.1f, s1sq=%.1f, s2sq=%.1f):\n",
              i, m, tau2, s1sq, s2sq))
  cat(sprintf("  Gradient match: %s\n", if (isTRUE(grad_ok)) "OK" else grad_ok))

  # V comparison: full 5x5 vs paper 3x3
  L1 <- s1sq; L2 <- s2sq; Q1 <- s1sq^2; Q2 <- s2sq^2; Q12 <- s1sq * s2sq

  Sf <- build_Sigma_full(m, tau2, s1sq, s2sq, L1, L2, Q1, Q2, Q12)
  Sp <- build_Sigma_paper(tau2, L1, L2, Q1, Q2, Q12)

  V_full  <- as.numeric(t(ga) %*% Sf %*% ga)
  V_paper <- as.numeric(t(ga[3:5]) %*% Sp %*% ga[3:5])

  V_ok <- all.equal(V_full, V_paper, tolerance = 1e-12)
  cat(sprintf("  V_full = %.8f, V_paper = %.8f, match: %s\n\n",
              V_full, V_paper, if (isTRUE(V_ok)) "OK" else V_ok))

  if (!isTRUE(grad_ok) || !isTRUE(V_ok)) all_pass <- FALSE
}

cat(sprintf("Overall: %s\n", if (all_pass) "ALL PASSED" else "SOME FAILED"))
