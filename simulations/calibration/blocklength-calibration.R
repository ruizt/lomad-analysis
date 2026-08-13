# blocklength-calibration.R -- does FDR control survive short blocks?
#
# The Morro Bay analysis admits a block only if its presmoothed length clears
# min_len = k*s + h. Since m = n - s - h + 2, that floor is a floor on m/s:
# min_len = k*s + h gives m/s = k - 1 + 2/s. The power sweep fixed n = 4s in
# every cell (s = min(60h, floor(n/4))), so everything ever simulated sits at
# m/s ~ 2.95. The analysis has been running at k = 3, i.e. m/s ~ 2.03, and we
# want to know whether k = 2.5 (m/s ~ 1.53) still controls FDR.
#
# Nothing about multiplicity is at issue: the correction is applied once,
# globally, over the pooled family, so shrinking a block does not change how
# many tests there are. What degrades is per-block estimation -- rho-hat and
# V-hat are both built from within-block windows -- and that shows up as raw
# p-values that are no longer uniform under the null, which no correction can
# repair.
#
# Two estimands, both under the global null (d = 0, so every window is null):
#
#   (1) per-block raw p-value calibration.  P(p <= 0.05) should be 0.05 at
#       every block length. This isolates the estimation layer.
#   (2) end-to-end global-null FDR.  Blocks drawn with the length distribution
#       the real threshold admits, p-values pooled, one BY step-up. With all
#       nulls true, FDR = P(at least one rejection), so this should be <= alpha.
#
# Usage:
#   Rscript simulations/calibration/blocklength-calibration.R
# Output:
#   simulations/calibration/results/blocklength-calibration.rds

suppressPackageStartupMessages({
  library(lomad); library(dplyr); library(tidyr); library(readr); library(purrr)
  library(parallel)
})

OUT_DIR <- "simulations/calibration/results"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- Operating point --------------------------------------------------------
# Matched to the Morro Bay fits: h and s are the analysis values, and the AR(1)
# coefficients are the medians over the 26 blocks that clear k = 2.5. lambda is
# set so the median local correlation lands near the 0.61 seen in the data.

S_WIN  <- 60L
H_WIN  <- 4L
ALPHA  <- 0.05
PHI1   <- 0.30
PHI2   <- 0.25
LAMBDA <- 1.5
NREP   <- 1000L
NCORE  <- max(1L, parallel::detectCores() - 2L)
pmap_ <- function(x, f) parallel::mclapply(x, f, mc.cores = NCORE)

# min_len for a given k, and the m/s it implies
min_len <- function(k) as.integer(k * S_WIN) + H_WIN
ms_of   <- function(k) k - 1 + 2 / S_WIN
K_GRID  <- c(4, 3, 2.5, 2)

# ---- One null replicate -----------------------------------------------------
# d = 0 gives two trends with zero L2 separation, so the local correlation sits
# at its noise-determined ceiling everywhere and every window is a true null.

null_block <- function(n, seed) {
  tr <- sim_trends(n, d = 0, method = "rs", bw = 50, seed = seed)
  sm <- suppressMessages(sim_noise_pair(
    tr, h = H_WIN, lambda_target = LAMBDA,
    ar.coefs = PHI1, s = S_WIN, seed = seed + 1L))
  sm2 <- suppressMessages(sim_noise_pair(
    tr, h = H_WIN, lambda_target = LAMBDA,
    ar.coefs = PHI2, s = S_WIN, seed = seed + 2L))
  fit <- suppressWarnings(suppressMessages(
    lomad_fit(sm$y1, sm2$y2, h = H_WIN, s = S_WIN)))
  tst <- suppressWarnings(suppressMessages(lomad_test(fit, alpha = ALPHA)))
  list(p = tst$p_values[fit$valid_idx], m = length(fit$valid_idx),
       rho = median(fit$rho, na.rm = TRUE))
}

# ---- (1) Per-block calibration across lengths -------------------------------

N_GRID <- c(min_len(2), min_len(2.5), min_len(3), min_len(4), 364L, 600L)

message("using ", NCORE, " cores")
message("(1) per-block raw p-value calibration, ", NREP, " reps x ",
        length(N_GRID), " lengths")

calib <- map_dfr(N_GRID, function(n) {
  reps <- pmap_(seq_len(NREP), \(i) null_block(n, seed = 10000L + 7L * i))
  p    <- unlist(map(reps, "p"))
  m    <- median(map_dbl(reps, "m"))
  tibble(n = n, m = m, ms = m / S_WIN,
         rho = median(map_dbl(reps, "rho"), na.rm = TRUE),
         n_p = length(p),
         p05 = mean(p <= 0.05), p01 = mean(p <= 0.01),
         p001 = mean(p <= 0.001),
         ks_p = suppressWarnings(ks.test(p, "punif")$p.value))
})
print(calib)

# ---- (2) End-to-end global-null FDR at each threshold -----------------------
# The ensemble is the real one: the presmoothed lengths of the blocks that the
# threshold actually admits, so the pooled family size M matches the analysis.

blk <- read_csv("_mb-data/ph_o2_blocks.csv", show_col_types = FALSE) |>
  count(location, block_id, name = "n_hr") |>
  mutate(n_6h = n_hr %/% 6L)

message("\n(2) global-null FDR by threshold, ", NREP, " reps")

fdr <- map_dfr(K_GRID, function(k) {
  lens <- blk$n_6h[blk$n_6h >= min_len(k)]
  res <- bind_rows(pmap_(seq_len(NREP), function(i) {
    out <- map(seq_along(lens),
               \(j) null_block(lens[j], seed = 500000L + 977L * i + 13L * j))
    p_all <- unlist(map(out, "p"))
    M     <- length(p_all)
    rej   <- p.adjust(p_all, method = "BY") <= ALPHA
    tibble(any_rej = any(rej), n_rej = sum(rej), M = M)
  }))
  tibble(k = k, min_len = min_len(k), ms = ms_of(k), n_blocks = length(lens),
         M = median(res$M),
         fdr = mean(res$any_rej),
         fdr_se = sqrt(mean(res$any_rej) * (1 - mean(res$any_rej)) / NREP),
         mean_rej = mean(res$n_rej))
})
print(fdr)

saveRDS(list(calibration = calib, fdr = fdr,
             settings = list(s = S_WIN, h = H_WIN, alpha = ALPHA,
                             phi = c(PHI1, PHI2), lambda = LAMBDA,
                             nrep = NREP)),
        file.path(OUT_DIR, "blocklength-calibration.rds"))

cat("\n==== verdict ====\n")
cat(sprintf("nominal alpha = %.3f\n", ALPHA))
for (i in seq_len(nrow(fdr)))
  cat(sprintf("  k = %.1f (m/s %.2f, %2d blocks): global-null FDR %.3f (SE %.3f) %s\n",
              fdr$k[i], fdr$ms[i], fdr$n_blocks[i], fdr$fdr[i], fdr$fdr_se[i],
              ifelse(fdr$fdr[i] <= ALPHA + 2 * fdr$fdr_se[i], "OK", "EXCEEDS")))
