test_that("lomad_fit returns expected structure", {
  set.seed(42)
  n <- 200
  x1 <- cumsum(rnorm(n)) + rnorm(n, sd = 0.5)
  x2 <- cumsum(rnorm(n)) + rnorm(n, sd = 0.5)

  fit <- lomad_fit(x1, x2, q = 5, h = 30)

  expect_type(fit, "list")
  expect_named(fit, c(
    "null_model", "observed", "expected_asymptotic", "inputs",
    "trend_hat", "ma1", "ma2", "R", "p_series", "thresholds", "I", "valid_idx"
  ))
  expect_length(fit$trend_hat, n)
  expect_length(fit$R, n)
  expect_true(is.numeric(fit$observed$frac_state))
})

test_that("lomad_fit errors on invalid rho0", {
  x <- rnorm(100)
  expect_error(lomad_fit(x, x, rho0 = 1.5))
  expect_error(lomad_fit(x, x, rho0 = -1))
})

test_that("lomad_fit errors on mismatched series lengths", {
  expect_error(lomad_fit(rnorm(100), rnorm(50)))
})
