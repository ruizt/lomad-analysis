# Does trend leakage into the noise estimate survive the two-tau pipeline?
#
# The v0.0.1 pipeline formed noise residuals as y_k - (ma1 + ma2)/2. When the
# trends separate (d > 0) the discarded difference (nu_k - nu_bar) is real
# signal, so it leaked into the residual, inflated phi_hat, and made the test
# conservative -- power *fell* with d at phi = 0.8 (see power/design.md).
#
# The two-tau pipeline differences the observed series directly, so no trend
# estimate enters at all. This script compares the two on identical data.
#
# Metric of interest is not bias per se but its SLOPE IN d: a bias that is
# constant in d shifts the whole power curve, a bias that grows with d bends it.
#
# Usage (from the repo root):
#   Rscript numerical-checks/assess-trend-leakage.R

suppressPackageStartupMessages({ library(lomad); library(dplyr) })

N     <- 600L
NREP  <- 200L
H_VALS <- c(5L, 20L)   # 5 = the power sweep at n=600; 20 = the validation study

phi_of <- function(x) suppressWarnings(lomad:::.variogram_ar1(x)$ar)

one_rep <- function(d, phi, snr, seed, H_WIN) {
  set.seed(seed)
  tr  <- sim_trends(n = N, d = d, method = "smooth", bw = 50,
                    coupling = 0.8, seed = seed)
  sim <- sim_noise_pair(tr, h = H_WIN, lambda_target = snr,
                        ar.coefs = phi, seed = seed + 1L)

  ma1 <- as.numeric(stats::filter(sim$y1, rep(1 / H_WIN, H_WIN), sides = 1))
  ma2 <- as.numeric(stats::filter(sim$y2, rep(1 / H_WIN, H_WIN), sides = 1))
  shared <- (ma1 + ma2) / 2
  ok <- which(!is.na(shared))

  c(old1 = phi_of(sim$y1[ok] - shared[ok]),      # v0.0.1: shared-trend resid
    old2 = phi_of(sim$y2[ok] - shared[ok]),
    new1 = phi_of(sim$y1),                       # two-tau: direct differencing
    new2 = phi_of(sim$y2))
}

grid <- expand.grid(d = c(0, 0.5, 1.0, 1.5), phi = c(0.3, 0.5, 0.8),
                    snr = c(0.5, 1.5), h = H_VALS)

res <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  g <- grid[i, ]
  m <- vapply(seq_len(NREP),
              function(r) one_rep(g$d, g$phi, g$snr, 9000L + 37L * i + r, g$h),
              numeric(4))
  data.frame(h = g$h, d = g$d, phi = g$phi, snr = g$snr,
             old = mean(c(m["old1", ], m["old2", ])),
             new = mean(c(m["new1", ], m["new2", ])))
}))
res$bias_old <- res$old - res$phi
res$bias_new <- res$new - res$phi

cat("\n=== phi_hat bias, shared-trend residuals (v0.0.1) vs direct differencing\n\n")
cat(sprintf("%5s %5s %6s | %8s %8s | %8s %8s\n",
            "phi", "snr", "d", "old", "bias", "new", "bias"))
for (i in seq_len(nrow(res))) with(res[i, ],
  cat(sprintf("%5.1f %5.1f %6.2f | %8.3f %+8.3f | %8.3f %+8.3f\n",
              phi, snr, d, old, bias_old, new, bias_new)))

cat("\n=== slope of phi_hat bias in d (the mechanism that bent power curves)\n\n")
slopes <- res |>
  group_by(h, phi, snr) |>
  summarise(slope_old = coef(lm(bias_old ~ d))[2],
            slope_new = coef(lm(bias_new ~ d))[2], .groups = "drop")
cat(sprintf("%3s %5s %5s | %11s %11s   %s\n", "h", "phi", "snr",
            "old d-slope", "new d-slope", "reduction"))
for (i in seq_len(nrow(slopes))) with(slopes[i, ]
  , cat(sprintf("%3d %5.1f %5.1f | %+11.4f %+11.4f   %.1fx\n",
                h, phi, snr, slope_old, slope_new,
                abs(slope_old) / max(abs(slope_new), 1e-6))))

saveRDS(list(res = res, slopes = slopes),
        "numerical-checks/_trend-leakage.rds")
cat("\nwrote numerical-checks/_trend-leakage.rds\n")
