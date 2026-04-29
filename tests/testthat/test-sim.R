# Tests for sim_trends(), sim_noise(), sim_noise_pair(), and key internals

# ---- sim_trends ----

test_that("sim_trends: returns correct structure and length", {
  tr <- sim_trends(n = 200, d = 1, seed = 4831)
  expect_type(tr, "list")
  expect_named(tr, c("x1", "x2", "x_mean", "w"))
  expect_length(tr$x1, 200)
  expect_length(tr$w, 200)
})

test_that("sim_trends: L2 distance equals d for all methods", {
  for (d in c(0.5, 2, 5)) {
    tr <- sim_trends(200, d = d, seed = 7713)
    dist <- sqrt(sum((tr$x1 - tr$x2)^2))
    expect_equal(dist, d, tolerance = 1e-8,
                 label = sprintf("dist method, d = %g", d))
  }

  tr_s <- sim_trends(200, d = 2, method = "smooth", bw = 30,
                     coupling = 0.8, seed = 7713)
  expect_equal(sqrt(sum((tr_s$x1 - tr_s$x2)^2)), 2, tolerance = 1e-8)

  tr_c <- sim_trends(200, d = 2, method = "cross", bw = 30,
                     coupling = 0.8, seed = 7713)
  expect_equal(sqrt(sum((tr_c$x1 - tr_c$x2)^2)), 2, tolerance = 1e-8)

  tr_r <- sim_trends(200, d = 2, method = "rate", rate = 0.02, seed = 7713)
  expect_equal(sqrt(sum((tr_r$x1 - tr_r$x2)^2)), 2, tolerance = 1e-8)
})

test_that("sim_trends: method='dist' gives w = 0", {
  tr <- sim_trends(100, d = 1, method = "dist", seed = 5122)
  expect_true(all(tr$w == 0))
})

test_that("sim_trends: method='smooth' gives w in (0, 1)", {
  tr <- sim_trends(500, d = 1, method = "smooth", bw = 50,
                   coupling = 0.8, seed = 5122)
  expect_true(all(tr$w > 0 & tr$w < 1))
})

test_that("sim_trends: method='rate' gives w in [0, 1]", {
  tr <- sim_trends(500, d = 1, method = "rate", rate = 0.01, seed = 5122)
  expect_true(all(tr$w >= 0 & tr$w <= 1))
})

test_that("sim_trends: custom w (vector and function)", {
  w_vec <- rep(c(0, 1), each = 100)
  tr_v <- sim_trends(200, d = 2, w = w_vec, seed = 3481)
  expect_equal(tr_v$w, w_vec)
  expect_equal(sqrt(sum((tr_v$x1 - tr_v$x2)^2)), 2, tolerance = 1e-8)

  tr_f <- sim_trends(200, d = 2, w = \(n) rep(0.5, n), seed = 3481)
  expect_equal(tr_f$w, rep(0.5, 200))
})

test_that("sim_trends: errors on bad inputs", {
  expect_error(sim_trends(100, nb = 24), "`nb` must be odd")
  expect_error(sim_trends(100, w = rep(0, 50)), "`w` must have length")
})

test_that("sim_trends: seed reproducibility", {
  tr1 <- sim_trends(100, d = 1, seed = 8192)
  tr2 <- sim_trends(100, d = 1, seed = 8192)
  expect_identical(tr1, tr2)
})


# ---- .generate_coef_pair ----

test_that(".generate_coef_pair: distance between coefs equals d", {
  for (d in c(0.5, 1, 3)) {
    cp <- .generate_coef_pair(nb = 25, d = d, seed = 2244)
    dist <- sqrt(sum((cp$coef2 - cp$coef1)^2))
    expect_equal(dist, d, tolerance = 1e-8,
                 label = sprintf("coef pair distance, d = %g", d))
  }
})


# ---- .make_w_smooth ----

test_that(".make_w_smooth: output in (0, 1)", {
  set.seed(6012)
  w <- .make_w_smooth(500, bw = 50, coupling = 0.8)
  expect_length(w, 500)
  expect_true(all(w > 0 & w < 1))
})


# ---- .make_w_rate ----

test_that(".make_w_rate: output in [0, 1]", {
  w <- .make_w_rate(500, rate = 0.01)
  expect_length(w, 500)
  expect_true(all(w >= 0 & w <= 1))
})


# ---- sim_noise / sim_noise_pair ----

test_that("sim_noise_pair: returns expected structure", {
  tr <- sim_trends(200, d = 1, seed = 9312)
  sp <- suppressMessages(
    sim_noise_pair(tr, h = 10, ar.coefs = 0.3, lambda_target = 2, seed = 9312)
  )
  expect_named(sp, c("y1", "y2", "x1", "x2", "noise"))
  expect_length(sp$y1, 200)
  expect_true(sp$noise$series1$sigma > 0)
})

test_that("sim_noise_pair: seed reproducibility", {
  tr <- sim_trends(200, d = 1, seed = 1142)
  sp1 <- suppressMessages(
    sim_noise_pair(tr, h = 10, ar.coefs = 0.3, lambda_target = 1, seed = 4401)
  )
  sp2 <- suppressMessages(
    sim_noise_pair(tr, h = 10, ar.coefs = 0.3, lambda_target = 1, seed = 4401)
  )
  expect_identical(sp1$y1, sp2$y1)
})
