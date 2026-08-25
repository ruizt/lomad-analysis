## sensitivity.R -- are the Morro Bay results stable near h = 4, s = 60?
##
## Inputs
##   _mb-data/ph_o2_blocks.csv          from mb-analysis/process_blocks.R
##
## Outputs -> mb-analysis/_tbl/
##   tbl-mb-sensitivity.csv             one row per (h, s) cell, h = 3 excluded
##
## The block set is held fixed at the length the reported analysis requires, so
## cells differ only in the fit and not in the data they see, and the reference
## cell reproduces Section 4. A cell needs only h + s - 1 observations to yield a
## window, so every cell is fittable on this set. Presmoothing does not depend on
## h or s and is done once.
##
## Usage (from the repo root):
##   Rscript mb-analysis/sensitivity.R

suppressPackageStartupMessages({
  library(tidyverse)
  library(lomad)
})
source("mb-analysis/utils.R")

tbl_out <- "mb-analysis/_tbl"; fs::dir_create(tbl_out)

ALPHA  <- 0.05
H_GRID <- c(3L, 4L, 5L)   # h = 3 is fitted but not tabulated
S_GRID <- c(50L, 60L, 70L)
H_REF  <- 4L; S_REF <- 60L
H_TBL  <- c(4L, 5L)       # rows kept for the supplement table
CV_K   <- 50L; CV_TEST <- 0.2
BLOCK_MIN <- as.integer(2.5 * S_REF) + H_REF   # matches MIN_LEN in analysis.R

# ---- Presmooth once ---------------------------------------------------------

block_data <- read_csv("_mb-data/ph_o2_blocks.csv", show_col_types = FALSE)

presm_all <- lapply(split(block_data, block_data$location), function(loc_dat) {
  out <- loc_dat |> group_split(block_id) |>
    lapply(presmooth_tidal, cols = c("o2", "ph"), step = "6 hours")
  setNames(out, vapply(out, \(b) as.character(b$block_id[[1]]), character(1)))
})
presm_all <- lapply(presm_all, function(bp)
  bp[vapply(bp, \(b) length(b$o2) >= BLOCK_MIN, logical(1))])

cat(sprintf("fixed block set (>= %d obs): %s\n", BLOCK_MIN,
            paste(sprintf("%s = %d", names(presm_all), lengths(presm_all)),
                  collapse = ", ")))

# ---- Fit one cell -----------------------------------------------------------
# One BY step-up over the pooled p-values, matching analysis.R.

run_cell <- function(h_win, s_win) {
  fits <- lapply(names(presm_all), function(loc) {
    out <- lapply(presm_all[[loc]], function(b) tryCatch({
      f <- lomad_fit(b$o2, b$ph, h = h_win, s = s_win)
      list(fit = f, tst = lomad_test(f, alpha = ALPHA),
           datetime = b$datetime, o2 = b$o2, ph = b$ph)
    }, error = function(e) NULL))
    Filter(Negate(is.null), out)
  })
  names(fits) <- names(presm_all)

  pooled <- bind_rows(lapply(names(fits), function(loc)
    bind_rows(lapply(names(fits[[loc]]), function(nm) {
      x <- fits[[loc]][[nm]]
      i <- x$fit$valid_idx
      tibble(location = loc, block_id = nm, datetime = x$datetime[i],
             o2 = x$o2[i], ph = x$ph[i], R = x$fit$R[i], rho = x$fit$rho[i],
             p = x$tst$p_values[i])
    }))))

  M <- nrow(pooled)
  k <- sum(p.adjust(pooled$p, method = "BY") <= ALPHA)
  p_star <- ALPHA * max(k, 1L) / (M * sum(1 / seq_len(M)))
  mutate(pooled, h = h_win, s = s_win, rejected = p <= p_star)
}

cells <- expand_grid(h = H_GRID, s = S_GRID)
key   <- sprintf("h%d_s%d", cells$h, cells$s)
res   <- setNames(pmap(cells, \(h, s) {
  cat(sprintf("  fitting h = %d, s = %d\n", h, s)); run_cell(h, s)
}), key)
ref <- res[[sprintf("h%d_s%d", H_REF, S_REF)]]

# ---- Summarise --------------------------------------------------------------

PH_SCALE <- read_csv("_mb-data/scale_constants.csv", show_col_types = FALSE) |>
  filter(variable == "ph") |> pull(scale)

# detection intensity by calendar month, divided by the cell's overall rate so
# cells with different rates are comparable
monthly_rel <- function(d) {
  m <- d |> mutate(m = month(datetime)) |>
    summarise(rate = mean(rejected), .by = m) |>
    complete(m = 1:12, fill = list(rate = NA_real_)) |> arrange(m) |> pull(rate)
  m / mean(m, na.rm = TRUE)
}
ref_rel <- monthly_rel(ref)

roll_sum <- function(v, n) {
  cs <- c(0, cumsum(v)); k <- length(v); o <- rep(NA_real_, k)
  if (k >= n) o[n:k] <- cs[(n:k) + 1] - cs[(n:k) - n + 1]
  o
}

# local RMSPE ratio, as in analysis.R: per-station fit on unflagged training
# observations, each window scored on its held-out observations only
rmspe_ratio <- function(d) {
  d <- arrange(d, location, block_id, datetime)
  w <- d$h[1] + d$s[1] - 1L
  key <- paste(d$location, d$block_id)
  set.seed(4471)
  acc <- matrix(NA_real_, nrow(d), CV_K)
  for (j in seq_len(CV_K)) {
    test <- runif(nrow(d)) < CV_TEST
    f <- d[!test, ] |> group_by(location) |>
      summarise(b0 = coef(lm(ph ~ o2))[1], b1 = coef(lm(ph ~ o2))[2], .groups = "drop")
    m  <- match(d$location, f$location)
    e2 <- ifelse(test, (d$ph - (f$b0[m] + f$b1[m] * d$o2))^2, 0)
    for (k in unique(key)) {
      i <- which(key == k); if (length(i) < w) next
      num <- roll_sum(e2[i], w); den <- roll_sum(as.numeric(test[i]), w)
      acc[i, j] <- ifelse(!is.na(den) & den > 0, sqrt(num / den), NA_real_)
    }
  }
  r <- PH_SCALE * rowMeans(acc, na.rm = TRUE)
  ok <- is.finite(r)
  mean(r[ok & d$rejected]) / mean(r[ok & !d$rejected])
}

summarise_cell <- function(d) {
  runs <- unlist(lapply(unique(paste(d$location, d$block_id)), function(k) {
    x <- d$rejected[paste(d$location, d$block_id) == k]
    r <- rle(x); r$lengths[r$values]
  }))
  rel <- monthly_rel(d)
  tibble(h = d$h[1], s = d$s[1],
         smooth_h = d$h[1] * 6, window_d = d$s[1] * 6 / 24,
         windows  = nrow(d),
         rate_BM  = mean(d$rejected[d$location == "BM1"]),
         rate_BH  = mean(d$rejected[d$location == "BS1"]),
         rate_ratio = rate_BM / rate_BH,
         n_runs   = length(runs),
         run_med_d = median(runs) * 6 / 24,
         monthly_pearson  = suppressWarnings(cor(rel, ref_rel, use = "complete.obs")),
         monthly_spearman = suppressWarnings(
           cor(rel, ref_rel, method = "spearman", use = "complete.obs")),
         rmspe_ratio = rmspe_ratio(d))
}

out <- map_dfr(res, summarise_cell) |> filter(h %in% H_TBL) |> arrange(h, s)
write_csv(out, file.path(tbl_out, "tbl-mb-sensitivity.csv"))
print(as.data.frame(out), digits = 3, row.names = FALSE)
cat(sprintf("\nWrote %s\n", file.path(tbl_out, "tbl-mb-sensitivity.csv")))
