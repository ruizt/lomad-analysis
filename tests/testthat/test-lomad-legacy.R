# Tests for legacy state-process pipeline
# lomad_fit(method="state"), lomad_test(method="boot"/"mc"/"analytic")

test_that("lomad_fit(method='state'): returns expected structure", {
  set.seed(4411)
  n <- 200
  x1 <- cumsum(rnorm(n)) + rnorm(n, sd = 0.5)
  x2 <- cumsum(rnorm(n)) + rnorm(n, sd = 0.5)

  fit <- suppressMessages(lomad_fit(x1, x2, method = "state"))

  expect_equal(fit$method, "state")
  expect_true("null_model" %in% names(fit))
  expect_true("observed" %in% names(fit))
  expect_true("I" %in% names(fit))
  expect_length(fit$R, n)
  expect_true(is.numeric(fit$observed$frac_state))
})

test_that("lomad_fit(method='state'): errors on invalid rho0", {
  expect_error(
    lomad_fit(rnorm(100), rnorm(100), method = "state", rho0 = 1.5),
    "rho0"
  )
})

test_that("lomad_test(method='analytic'): returns p-value for frac_state", {
  skip_on_cran()

  set.seed(4411)
  n <- 300
  x1 <- cumsum(rnorm(n)) + rnorm(n, sd = 0.5)
  x2 <- cumsum(rnorm(n)) + rnorm(n, sd = 0.5)

  fit <- suppressMessages(lomad_fit(x1, x2, method = "state"))
  tst <- suppressMessages(lomad_test(fit, method = "analytic"))

  expect_true("p_values" %in% names(tst))
  expect_true(is.numeric(tst$p_values$frac_state))
})
