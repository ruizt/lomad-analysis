# Tests for CLT building blocks: arma_acov, acov_sums, compute_rho,
# compute_V, compute_tau_sq, .ma_filter_acov

# ---- arma_acov ----

test_that("arma_acov: white noise returns [sigma2, 0, ...]", {
  g <- arma_acov(sigma2 = 2, lag_max = 5)
  expect_length(g, 6)
  expect_equal(g[1], 2)
  expect_equal(g[-1], rep(0, 5))
})

test_that("arma_acov: AR(1) matches exact formula", {
  phi <- 0.7
  sigma2 <- 1.5
  g <- unname(arma_acov(ar = phi, sigma2 = sigma2, lag_max = 10))

  gamma0_exact <- sigma2 / (1 - phi^2)
  # ARMAtoMA truncation introduces small errors at gamma(0)
  expect_equal(g[1], gamma0_exact, tolerance = 1e-2)

  for (l in 1:10) {
    expect_equal(g[l + 1], gamma0_exact * phi^l, tolerance = 1e-2,
                 label = sprintf("gamma(%d)", l))
  }
})

test_that("arma_acov: MA(1) has correct structure", {
  theta <- 0.6
  sigma2 <- 2
  g <- unname(arma_acov(ma = theta, sigma2 = sigma2, lag_max = 5))

  expect_equal(g[1], sigma2 * (1 + theta^2), tolerance = 1e-6)
  expect_equal(g[2], sigma2 * theta, tolerance = 1e-6)
  expect_true(all(abs(g[3:6]) < 1e-8))
})


# ---- acov_sums ----

test_that("acov_sums: symmetry in Q12", {
  a <- arma_acov(ar = 0.5, sigma2 = 1, lag_max = 50)
  b <- arma_acov(ar = 0.3, sigma2 = 2, lag_max = 50)
  s1 <- acov_sums(a, b)
  s2 <- acov_sums(b, a)
  expect_equal(s1$Q12, s2$Q12, tolerance = 1e-10)
})

test_that("acov_sums: single-series shortcut", {
  a <- arma_acov(ar = 0.5, sigma2 = 1, lag_max = 50)
  s <- acov_sums(a)
  expect_equal(s$L1, s$L2)
  expect_equal(s$Q1, s$Q2)
  expect_equal(s$Q1, s$Q12)
})


# ---- compute_rho ----

test_that("compute_rho: tau_sq = 0 returns 0", {
  expect_equal(compute_rho(0, 1, 1), 0)
})

test_that("compute_rho: exact formula", {
  tau <- 3; s1 <- 1; s2 <- 2
  expected <- tau / sqrt((tau + s1) * (tau + s2))
  expect_equal(compute_rho(tau, s1, s2), expected)
})

test_that("compute_rho: vectorised over tau_sq", {
  rho <- compute_rho(c(0, 1, 10), 1, 1)
  expect_length(rho, 3)
  expect_equal(rho[1], 0)
  expect_true(all(rho >= 0 & rho < 1))
  expect_true(rho[3] > rho[2])
})


# ---- compute_V ----

test_that("compute_V: returns positive for valid inputs", {
  V <- compute_V(tau_sq = 2, sigma1_sq = 1, sigma2_sq = 1,
                 L1 = 3, L2 = 3, Q1 = 2, Q2 = 2, Q12 = 2)
  expect_true(V > 0)
})

test_that("compute_V: vectorised over tau_sq", {
  V <- compute_V(tau_sq = c(1, 5, 10), sigma1_sq = 1, sigma2_sq = 1,
                 L1 = 3, L2 = 3, Q1 = 2, Q2 = 2, Q12 = 2)
  expect_length(V, 3)
  expect_true(all(V > 0))
})


# ---- compute_tau_sq ----

test_that("compute_tau_sq: constant trend returns 0", {
  tau <- compute_tau_sq(rep(5, 100), s = 10)
  expect_true(all(tau[10:100] == 0))
})

test_that("compute_tau_sq: NA for first s-1 positions", {
  tau <- compute_tau_sq(rnorm(50), s = 20)
  expect_true(all(is.na(tau[1:19])))
  expect_true(!any(is.na(tau[20:50])))
})


# ---- .ma_filter_acov ----

test_that(".ma_filter_acov: h=1 returns original acov", {
  acov_raw <- unname(arma_acov(ar = 0.5, sigma2 = 1, lag_max = 20))
  result   <- .ma_filter_acov(acov_raw, h = 1, lag_max = 10)
  expect_equal(result, acov_raw[1:11], tolerance = 1e-10)
})

test_that(".ma_filter_acov: white noise with h>1", {
  sigma2   <- 2
  acov_raw <- unname(arma_acov(sigma2 = sigma2, lag_max = 30))
  h        <- 5
  result   <- .ma_filter_acov(acov_raw, h = h, lag_max = 5)
  expect_equal(result[1], sigma2 / h, tolerance = 1e-10,
               label = "gamma_eta(0) = sigma2/h for white noise")
})
