# Tests for lomad_test_identity()

test_that("lomad_test_identity: returns expected structure", {
  set.seed(2291)
  n <- 200
  trend <- cumsum(rnorm(n, sd = 0.1))
  x1 <- trend + rnorm(n, sd = 0.5)
  x2 <- trend + rnorm(n, sd = 0.5)

  res <- suppressMessages(lomad_test_identity(x1, x2, q = 11))

  expect_named(res, c("I", "p_raw", "p_adj", "D", "se_D", "m_eff",
                      "global_p", "T_stat", "fit_diff"))
  expect_length(res$I, n)
  expect_true(res$se_D > 0)
  expect_true(res$m_eff > 0)
  expect_true(res$m_eff <= n)
  expect_true(res$global_p >= 0 && res$global_p <= 1)
})

test_that("lomad_test_identity: noise_override bypasses estimation", {
  set.seed(7142)
  n <- 200
  trend <- cumsum(rnorm(n, sd = 0.1))
  x1 <- trend + rnorm(n, sd = 0.5)
  x2 <- trend + rnorm(n, sd = 0.5)

  res <- suppressMessages(
    lomad_test_identity(x1, x2, q = 11,
                        noise_override = list(ar = 0.0, sigma2 = 0.25))
  )

  expect_named(res, c("I", "p_raw", "p_adj", "D", "se_D", "m_eff",
                      "global_p", "T_stat", "fit_diff"))
  expect_true(res$se_D > 0)
  expect_true(res$global_p >= 0 && res$global_p <= 1)
  expect_null(res$fit_diff)
})

test_that("lomad_test_identity: noise_override with AR(1) works", {
  tr <- sim_trends(n = 300, d = 0, seed = 5519)
  sp <- suppressMessages(
    sim_noise_pair(tr, h = 10, ar.coefs = 0.5, lambda_target = 2, seed = 5519)
  )
  true_sigma2 <- mean(c(sp$noise$series1$sigma, sp$noise$series2$sigma)^2)

  res <- suppressMessages(
    lomad_test_identity(sp$y1, sp$y2, q = 10, alpha = 0.05,
                        noise_override = list(ar = 0.5, sigma2 = true_sigma2))
  )

  expect_length(res$I, 300)
  expect_true(res$se_D > 0)
})

test_that("lomad_test_identity: few rejections under null (d = 0)", {
  tr <- sim_trends(n = 500, d = 0, seed = 3378)
  sp <- suppressMessages(
    sim_noise_pair(tr, h = 10, ar.coefs = 0.3, lambda_target = 2, seed = 3378)
  )

  res <- suppressMessages(lomad_test_identity(sp$y1, sp$y2, q = 11))
  frac_rejected <- mean(res$I, na.rm = TRUE)
  expect_true(frac_rejected < 0.2,
              label = sprintf("frac rejected = %.3f", frac_rejected))
})
