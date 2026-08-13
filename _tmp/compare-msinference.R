## Does an existing multiscale trend-comparison method hold up under local
## affine similarity? One series, two framings.
suppressPackageStartupMessages({library(lomad); library(MSinference)})

n <- 1250L; s <- 50L; h <- 5L; nb <- 151L
SNR <- 0.5; PHI <- 0.5; ALPHA <- 0.05; SIM_RUNS <- 200L

set.seed(6001)
tr <- sim_trends(n = n, d = 0, method = "rs", bw = 50, nb = nb, seed = 6001,
                 affine_s = s, affine_cap = 0.015)
sim <- suppressMessages(sim_noise_pair(tr, h = h, lambda_target = SNR,
                                       ar.coefs = PHI, seed = 1001))
y1 <- sim$y1; y2 <- sim$y2

cat(sprintf("H_0: nu2 = a_t + b_t nu1, n = %d, s = %d, SNR = %.1f, phi = %.1f\n",
            n, s, SNR, PHI))
cat(sprintf("realized affine drift: b in [%.3f, %.3f], a/sd(x1) in [%.3f, %.3f]\n\n",
            min(tr$b), max(tr$b), min(tr$a)/sd(tr$x1), max(tr$a)/sd(tr$x1)))

## ---- lomad, for reference --------------------------------------------------
f  <- suppressMessages(suppressWarnings(lomad_fit(y1, y2, h = h, s = s)))
ts <- suppressMessages(lomad_test(f, alpha = ALPHA))
vi <- f$valid_idx
cat(sprintf("lomad          : %d / %d windows rejected (%.4f)\n\n",
            sum(ts$rejected[vi], na.rm = TRUE), length(vi),
            mean(ts$rejected[vi], na.rm = TRUE)))

## ---- MSinference -----------------------------------------------------------
grid <- MSinference::construct_grid(n)
cat(sprintf("MSinference grid: %d (location, bandwidth) pairs\n\n", nrow(grid$gset)))

run_ms <- function(a, b, label) {
  lrv <- c(MSinference::estimate_lrv(a, q = 25, r_bar = 10, p = 1)$lrv,
           MSinference::estimate_lrv(b, q = 25, r_bar = 10, p = 1)$lrv)
  t0 <- Sys.time()
  r <- MSinference::multiscale_test(data = cbind(a, b), sigma_vec = sqrt(lrv),
                                    n_ts = 2, grid = grid, alpha = ALPHA,
                                    sim_runs = SIM_RUNS)
  el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  gv <- r$gset_with_values[[1]]
  flagged <- sum(gv$test != 0)
  cat(sprintf("%-15s: stat = %.3f  quant = %.3f  reject = %-5s  %d / %d intervals flagged  [%.0f s]\n",
              label, r$stat, r$quant, r$stat > r$quant, flagged, nrow(gv), el))
  invisible(r)
}

## Easy: remove one global affine map, then compare trends directly.
zs <- function(v) (v - mean(v)) / sd(v)
run_ms(zs(y1), zs(y2), "global recentre")

## Hard: remove the map locally, the way the method under study does, then hand
## the realigned pair to the same simultaneous procedure. This asks whether the
## estimation error the realignment introduces is what breaks the inference.
roll_ab <- function(x, y, s) {
  n <- length(x); a <- numeric(n); b <- numeric(n)
  half <- s %/% 2L
  for (t in seq_len(n)) {
    lo <- max(1L, t - half); hi <- min(n, lo + s - 1L); lo <- max(1L, hi - s + 1L)
    w <- lo:hi
    vx <- stats::var(x[w])
    bt <- if (vx > 0) stats::cov(x[w], y[w]) / vx else 1
    if (!is.finite(bt) || abs(bt) < 1e-3) bt <- 1
    b[t] <- bt
    a[t] <- mean(y[w]) - bt * mean(x[w])
  }
  list(a = a, b = b)
}
ab <- roll_ab(y1, y2, s)
y2_re <- (y2 - ab$a) / ab$b
cat(sprintf("\nlocal realignment: b_hat in [%.3f, %.3f]\n", min(ab$b), max(ab$b)))
run_ms(zs(y1), zs(y2_re), "local realign")
