# Tests for lomad(), lomad_fit(), lomad_test()
# Primarily structural tests — statistical calibration is handled by
# simulation studies in dev/sims/.

# ---- Shared test data ----

# Stable null data: d = 0, high SNR, low phi, long series
make_null_data <- function(seed = 6247) {
  tr <- sim_trends(n = 500, d = 0, method = "dist", seed = seed)
  suppressMessages(
    sim_noise_pair(tr, h = 10, ar.coefs = 0.3, lambda_target = 2, seed = seed)
  )
}

# ---- lomad_fit (CLT) ----

test_that("lomad_fit: returns expected structure for CLT method", {
  sp <- make_null_data()
  fit <- suppressMessages(lomad_fit(sp$y1, sp$y2))

  expect_type(fit, "list")
  expect_equal(fit$method, "clt")

  expected_names <- c("method", "trend", "ma1", "ma2", "noise", "acov_sums",
                      "tau_sq", "rho", "V", "R", "valid_idx", "inputs")
  expect_true(all(expected_names %in% names(fit)))

  n <- length(sp$y1)
  expect_length(fit$trend, n)
  expect_length(fit$rho, n)
  expect_length(fit$V, n)
  expect_length(fit$R, n)
})

test_that("lomad_fit: rho in [0, 1] where defined", {
  sp <- make_null_data()
  fit <- suppressMessages(lomad_fit(sp$y1, sp$y2))

  rho_valid <- fit$rho[fit$valid_idx]
  expect_true(all(rho_valid >= 0 & rho_valid <= 1))
})

test_that("lomad_fit: V positive where defined", {
  sp <- make_null_data()
  fit <- suppressMessages(lomad_fit(sp$y1, sp$y2))

  V_valid <- fit$V[fit$valid_idx]
  expect_true(all(V_valid > 0))
})

test_that("lomad_fit: valid_idx is nonempty subset of 1:n", {
  sp <- make_null_data()
  fit <- suppressMessages(lomad_fit(sp$y1, sp$y2))

  expect_true(length(fit$valid_idx) > 0)
  expect_true(all(fit$valid_idx >= 1 & fit$valid_idx <= length(sp$y1)))
})

test_that("lomad_fit: noise contains expected fields", {
  sp <- make_null_data()
  fit <- suppressMessages(lomad_fit(sp$y1, sp$y2))

  for (series in c("series1", "series2")) {
    expect_true(all(c("ar", "ma", "sigma2") %in% names(fit$noise[[series]])))
    expect_true(fit$noise[[series]]$sigma2 > 0)
  }
})

test_that("lomad_fit: errors on mismatched lengths", {
  expect_error(lomad_fit(rnorm(100), rnorm(50)), "same length")
})

test_that("lomad_fit: errors on missing input", {
  expect_error(lomad_fit(), "`x1` and `x2` are required")
})


# ---- lomad_test ----

test_that("lomad_test: returns expected structure", {
  sp <- make_null_data()
  fit <- suppressMessages(lomad_fit(sp$y1, sp$y2))
  tst <- suppressMessages(lomad_test(fit))

  expect_named(tst, c("Z", "p_values", "p_adj", "rejected", "alpha_eff", "inputs"))

  n <- length(sp$y1)
  expect_length(tst$Z, n)
  expect_length(tst$p_values, n)
  expect_length(tst$rejected, n)
})

test_that("lomad_test: p-values in [0, 1]", {
  sp <- make_null_data()
  fit <- suppressMessages(lomad_fit(sp$y1, sp$y2))
  tst <- suppressMessages(lomad_test(fit))

  p_valid <- tst$p_values[fit$valid_idx]
  expect_true(all(p_valid >= 0 & p_valid <= 1))
})

test_that("lomad_test: p_adj >= p_values (BY is conservative)", {
  sp <- make_null_data()
  fit <- suppressMessages(lomad_fit(sp$y1, sp$y2))
  tst <- suppressMessages(lomad_test(fit))

  idx <- fit$valid_idx
  expect_true(all(tst$p_adj[idx] >= tst$p_values[idx] - 1e-10))
})

test_that("lomad_test: alpha_eff < alpha", {
  sp <- make_null_data()
  fit <- suppressMessages(lomad_fit(sp$y1, sp$y2))
  tst <- suppressMessages(lomad_test(fit, alpha = 0.05))
  expect_true(tst$alpha_eff < 0.05)
})

test_that("lomad_test: rejected is logical", {
  sp <- make_null_data()
  fit <- suppressMessages(lomad_fit(sp$y1, sp$y2))
  tst <- suppressMessages(lomad_test(fit))
  expect_true(is.logical(tst$rejected[fit$valid_idx]))
})

test_that("lomad_test: errors on invalid alpha", {
  sp <- make_null_data()
  fit <- suppressMessages(lomad_fit(sp$y1, sp$y2))
  expect_error(lomad_test(fit, alpha = 0), "in \\(0, 1\\)")
  expect_error(lomad_test(fit, alpha = 1), "in \\(0, 1\\)")
})


# ---- noise_override ----

test_that("lomad_fit: noise_override (single spec) bypasses estimation", {
  sp <- make_null_data()
  fit <- suppressMessages(
    lomad_fit(sp$y1, sp$y2,
              noise_override = list(ar = 0.3, sigma2 = 0.01))
  )

  expect_equal(fit$method, "clt")
  expect_equal(fit$noise$series1$ar, 0.3)
  expect_equal(fit$noise$series2$ar, 0.3)
  expect_equal(fit$noise$series1$sigma2, 0.01)
  expect_equal(fit$noise$series1$ma, numeric(0))
})

test_that("lomad_fit: noise_override (per-series) uses distinct specs", {
  sp <- make_null_data()
  fit <- suppressMessages(
    lomad_fit(sp$y1, sp$y2,
              noise_override = list(
                list(ar = 0.2, sigma2 = 0.005),
                list(ar = 0.4, ma = 0.1, sigma2 = 0.015)
              ))
  )

  expect_equal(fit$noise$series1$ar, 0.2)
  expect_equal(fit$noise$series2$ar, 0.4)
  expect_equal(fit$noise$series2$ma, 0.1)
  expect_equal(fit$noise$series1$sigma2, 0.005)
  expect_equal(fit$noise$series2$sigma2, 0.015)
})

test_that("lomad_fit: noise_override produces valid fit for lomad_test", {
  sp <- make_null_data()
  fit <- suppressMessages(
    lomad_fit(sp$y1, sp$y2,
              noise_override = list(ar = 0.3, sigma2 = 0.01))
  )
  tst <- lomad_test(fit)

  expect_true(length(fit$valid_idx) > 0)
  p_valid <- tst$p_values[fit$valid_idx]
  expect_true(all(p_valid >= 0 & p_valid <= 1))
})


# ---- lomad() integration ----

test_that("lomad: end-to-end structural check", {
  sp <- make_null_data()
  res <- suppressMessages(lomad(sp$y1, sp$y2))

  expect_named(res, c("fit", "test"))
  expect_equal(res$fit$method, "clt")
  expect_true(length(res$fit$valid_idx) > 0)
  expect_true(all(c("Z", "p_values", "rejected") %in% names(res$test)))
})

test_that("lomad: runs on short series (n = 100)", {
  tr <- sim_trends(n = 100, d = 0, seed = 8881)
  sp <- suppressMessages(
    sim_noise_pair(tr, h = 5, ar.coefs = 0.3, lambda_target = 2, seed = 8881)
  )
  res <- suppressMessages(lomad(sp$y1, sp$y2))
  expect_type(res, "list")
})


# ---- Null calibration (slow, skip on CRAN) ----

test_that("lomad: null rejection rate < 3*alpha", {
  skip_on_cran()

  S <- 50
  alpha <- 0.05
  rejections <- 0L

  for (i in seq_len(S)) {
    tr <- sim_trends(n = 2000, d = 0, method = "dist", seed = 1000 + i)
    sp <- suppressMessages(
      sim_noise_pair(tr, h = 10, ar.coefs = 0.3, lambda_target = 2,
                     seed = 2000 + i)
    )
    res <- suppressMessages(lomad(sp$y1, sp$y2, alpha = alpha))
    if (any(res$test$rejected, na.rm = TRUE)) {
      rejections <- rejections + 1L
    }
  }

  rejection_rate <- rejections / S
  expect_true(rejection_rate < 3 * alpha,
              label = sprintf("Rejection rate %.3f < 3*alpha = %.3f",
                              rejection_rate, 3 * alpha))
})
