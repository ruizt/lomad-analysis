# station-comparison.R -- is Bay South actually more stable than Bay Mouth?
#
# Usage (from the repo root):
#   Rscript mb-analysis/station-comparison.R
#
# Input
#   _mb-data/lomad_windows.rds    written by full_analysis.R
#
# Output
#   mb-analysis/_img/station-comparison.png
#   mb-analysis/_tbl/station-comparison.csv
#
# Three questions, none of which the side-by-side time courses answer well.
#
#   1. Rate.  BM1 flags a much larger share of its windows than BS1. Both the
#      naive difference and its uncertainty are wrong if computed pointwise:
#      windows overlap by construction (s = 60 at a 6-hour step, so successive
#      windows share 59/60 of their data) and rejections arrive in runs. The
#      block is the independent unit, so every interval here comes from a
#      bootstrap over blocks. The rate is also confounded with exposure -- the
#      two stations were up at different times -- so it is recomputed on the
#      concurrent support alone.
#
#   2. Alignment.  Do the two stations decouple at the same moments? Measured
#      as lift, P(both flagged) / (P(BM1) * P(BS1)), on the support where both
#      were under analysis, profiled across lags so a lead/lag relationship
#      would show up as an off-zero peak. The null rotates each block's flag
#      vector circularly within the block, which preserves run structure, the
#      marginal rate and the exposure pattern exactly, and destroys only the
#      cross-station alignment.
#
#   3. Season.  Flagged fraction by month against observed exposure, with the
#      same within-block rotation null and a first-harmonic amplitude statistic.

suppressPackageStartupMessages({
  library(tidyverse); library(patchwork)
})

set.seed(20260727)
img_out <- "mb-analysis/_img"; fs::dir_create(img_out)
tbl_out <- "mb-analysis/_tbl"; fs::dir_create(tbl_out)

ggthm <- theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "black", linewidth = 0.1))

STEP_H <- 6                       # presmoothed sampling interval, hours
B      <- 2000L                   # bootstrap / permutation replicates

w <- readRDS("_mb-data/lomad_windows.rds") |>
  # rejected was NA-filled to FALSE upstream; a window is only under test where
  # the fit produced both a correlation and a benchmark for it
  mutate(valid = !is.na(R) & !is.na(rho),
         blk   = paste(location, block_id))

stopifnot(all(w$rejected[!w$valid] == FALSE))

cat(sprintf("%d windows, %d under test, %d blocks\n",
            nrow(w), sum(w$valid), n_distinct(w$blk)))

# =============================================================================
# 1. Rate
# =============================================================================

# Bootstrap over blocks: resample blocks with replacement within a station and
# recompute the pooled fraction. Blocks differ enormously in length (156 to
# 2742 windows), so this also propagates the fact that BS1's exposure is
# dominated by a single very long block.

boot_rate <- function(d, B) {
  blocks <- split(d, d$blk)
  vapply(seq_len(B), function(i) {
    s <- blocks[sample.int(length(blocks), replace = TRUE)]
    sum(vapply(s, \(b) sum(b$rejected), numeric(1))) /
      sum(vapply(s, \(b) sum(b$valid),  numeric(1)))
  }, numeric(1))
}

wv <- filter(w, valid)

rate_tbl <- wv |>
  group_by(location) |>
  group_modify(~{
    bs <- boot_rate(.x, B)
    tibble(n_blocks = n_distinct(.x$blk), n_valid = nrow(.x),
           n_rej = sum(.x$rejected), rate = sum(.x$rejected) / nrow(.x),
           lo = quantile(bs, 0.025), hi = quantile(bs, 0.975))
  }) |> ungroup()

cat("\n==== 1a. overall rate (block bootstrap, B =", B, ") ====\n")
print(rate_tbl |> mutate(across(c(rate, lo, hi), \(x) round(100 * x, 1))))

# Rate ratio, resampling both stations jointly
bm <- filter(wv, location == "BM1"); bs <- filter(wv, location == "BS1")
rr_boot <- boot_rate(bm, B) / boot_rate(bs, B)
cat(sprintf("\nrate ratio BM1/BS1 = %.2f  (95%% CI %.2f to %.2f)\n",
            rate_tbl$rate[1] / rate_tbl$rate[2],
            quantile(rr_boot, 0.025), quantile(rr_boot, 0.975)))

# Block-level: how many blocks detect anything at all? Coarse, but it does not
# lean on the within-block dependence structure at all.
blk_tbl <- wv |> group_by(location, blk) |>
  summarise(n = n(), rej = sum(rejected), frac = rej / n, .groups = "drop")
det <- blk_tbl |> group_by(location) |>
  summarise(blocks = n(), detecting = sum(rej > 0), .groups = "drop")
ft <- fisher.test(matrix(c(det$detecting, det$blocks - det$detecting), nrow = 2))
cat(sprintf("\nblocks with >=1 detection: BM1 %d/%d, BS1 %d/%d (Fisher p = %.3f)\n",
            det$detecting[1], det$blocks[1], det$detecting[2], det$blocks[2],
            ft$p.value))

# ---- 1b. Concurrent support only --------------------------------------------
# The stations were not up at the same times. Restricting to timestamps where
# BOTH were under test removes exposure as an explanation entirely.

grid <- wv |> select(location, datetime, blk, rejected) |>
  pivot_wider(id_cols = datetime, names_from = location,
              values_from = c(rejected, blk),
              values_fn = list(rejected = any, blk = first)) |>
  drop_na(rejected_BM1, rejected_BS1) |>
  arrange(datetime)

cat(sprintf("\n==== 1b. concurrent support: %d timestamps (%.0f days), %s .. %s ====\n",
            nrow(grid), nrow(grid) * STEP_H / 24,
            as.Date(min(grid$datetime)), as.Date(max(grid$datetime))))

conc_boot <- function(g, col, blkcol, B) {
  blocks <- split(g, g[[blkcol]])
  vapply(seq_len(B), function(i) {
    s <- blocks[sample.int(length(blocks), replace = TRUE)]
    mean(unlist(lapply(s, \(b) b[[col]])))
  }, numeric(1))
}
cb_bm <- conc_boot(grid, "rejected_BM1", "blk_BM1", B)
cb_bs <- conc_boot(grid, "rejected_BS1", "blk_BS1", B)

conc_tbl <- tibble(
  location = c("BM1", "BS1"),
  rate = c(mean(grid$rejected_BM1), mean(grid$rejected_BS1)),
  lo   = c(quantile(cb_bm, 0.025), quantile(cb_bs, 0.025)),
  hi   = c(quantile(cb_bm, 0.975), quantile(cb_bs, 0.975)))
print(conc_tbl |> mutate(across(c(rate, lo, hi), \(x) round(100 * x, 1))))
cat(sprintf("concurrent rate ratio = %.2f  (95%% CI %.2f to %.2f)\n",
            conc_tbl$rate[1] / conc_tbl$rate[2],
            quantile(cb_bm / cb_bs, 0.025), quantile(cb_bm / cb_bs, 0.975)))

# =============================================================================
# 2. Alignment
# =============================================================================

lift_of <- function(a, b) {
  ok <- !is.na(a) & !is.na(b)
  if (sum(ok) < 100) return(NA_real_)
  pa <- mean(a[ok]); pb <- mean(b[ok])
  if (pa == 0 || pb == 0) return(NA_real_)
  mean(a[ok] & b[ok]) / (pa * pb)
}

# Full 6-hourly grid spanning the record, NA where a station is not under test
full_grid <- tibble(datetime = seq(min(wv$datetime), max(wv$datetime),
                                   by = paste(STEP_H, "hours")))
series <- lapply(c("BM1", "BS1"), function(loc) {
  full_grid |>
    left_join(wv |> filter(location == loc) |>
                group_by(datetime) |> summarise(r = any(rejected),
                                                blk = first(blk), .groups = "drop"),
              by = "datetime")
})
names(series) <- c("BM1", "BS1")

a0 <- series$BM1$r; b0 <- series$BS1$r
obs_lift <- lift_of(a0, b0)
n_both   <- sum(!is.na(a0) & !is.na(b0) & a0 & b0)

# Lag profile: shift BS1 forward/backward and recompute
lags <- seq(-480L, 480L, by = 4L)          # +/- 120 days in 1-day steps
shift <- function(x, k) if (k == 0) x else if (k > 0)
  c(rep(NA, k), head(x, -k)) else c(tail(x, k), rep(NA, -k))
prof <- tibble(lag_steps = lags,
               lag_days  = lags * STEP_H / 24,
               lift = vapply(lags, \(k) lift_of(a0, shift(b0, k)), numeric(1)))

# Null: rotate each block's flag vector circularly within that block. Preserves
# run lengths, the marginal rate and the exposure pattern; destroys only the
# alignment between stations.
rotate_within_blocks <- function(s) {
  idx <- which(!is.na(s$r))
  out <- s$r
  for (b in split(idx, s$blk[idx])) {
    if (length(b) > 1L)
      out[b] <- s$r[b][(seq_along(b) + sample.int(length(b), 1L) - 1L) %% length(b) + 1L]
  }
  out
}
null_lift <- vapply(seq_len(B), \(i) lift_of(a0, rotate_within_blocks(series$BS1)),
                    numeric(1))
p_align <- (1 + sum(null_lift >= obs_lift, na.rm = TRUE)) / (1 + sum(!is.na(null_lift)))

cat("\n==== 2. alignment ====\n")
# Whether this question is answerable at all: how many co-flagged timestamps
# would independence predict? If that number is small, observing zero says
# nothing, and the lift profile below is decoration.
exp_both <- nrow(grid) * mean(grid$rejected_BM1) * mean(grid$rejected_BS1)
cat(sprintf("concurrent timestamps %d; BM1 flags %d, BS1 flags %d\n",
            nrow(grid), sum(grid$rejected_BM1), sum(grid$rejected_BS1)))
cat(sprintf("co-flagged: observed %d, expected under independence %.1f\n",
            n_both, exp_both))
cat(sprintf("  -> P(0 | independence) = %.2f; the concurrent record cannot\n",
            dpois(0, exp_both)))
cat("     separate 'never together' from 'unrelated'.\n")
cat(sprintf("lift at lag 0 = %.2f  (null median %.2f, 95%% %.2f-%.2f), p = %.3f\n",
            obs_lift, median(null_lift, na.rm = TRUE),
            quantile(null_lift, 0.025, na.rm = TRUE),
            quantile(null_lift, 0.975, na.rm = TRUE), p_align))
best <- prof |> filter(!is.na(lift)) |> slice_max(lift, n = 1)
cat(sprintf("largest lift over the profile: %.2f at lag %+.0f days\n",
            best$lift, best$lag_days))

# =============================================================================
# 3. Season
# =============================================================================

seas <- wv |>
  mutate(month = month(datetime), doy = yday(datetime)) |>
  group_by(location, month) |>
  summarise(n = n(), rej = sum(rejected), frac = rej / n, .groups = "drop")

# First-harmonic amplitude of the flagged indicator on day-of-year, against a
# within-block rotation null (which holds exposure and run structure fixed).
harm_amp <- function(r, doy) {
  ok <- !is.na(r)
  th <- 2 * pi * doy[ok] / 365.25
  sqrt(mean(r[ok] * cos(th))^2 + mean(r[ok] * sin(th))^2) / mean(r[ok])
}

cat("\n==== 3. season ====\n")
seas_test <- map_dfr(c("BM1", "BS1"), function(loc) {
  s <- series[[loc]] |> mutate(doy = yday(datetime))
  obs <- harm_amp(s$r, s$doy)
  nul <- vapply(seq_len(B), \(i) harm_amp(rotate_within_blocks(s), s$doy), numeric(1))
  tibble(location = loc, amp = obs, null_med = median(nul, na.rm = TRUE),
         p = (1 + sum(nul >= obs, na.rm = TRUE)) / (1 + sum(!is.na(nul))))
})
print(seas_test |> mutate(across(c(amp, null_med), \(x) round(x, 3))))

peak <- seas |> group_by(location) |> slice_max(frac, n = 1) |> ungroup()
cat("\npeak month by station:\n")
print(peak |> transmute(location, month = month.abb[month],
                        frac = round(100 * frac, 1), n))

# =============================================================================
# Figure
# =============================================================================

pal <- c(BM1 = "#C44E52", BS1 = "#4C72B0")

p_rate <- bind_rows(
    rate_tbl |> transmute(location, rate, lo, hi, panel = "all analysed time"),
    conc_tbl |> mutate(panel = "concurrent only")) |>
  mutate(panel = factor(panel, c("all analysed time", "concurrent only"))) |>
  ggplot(aes(location, 100 * rate, colour = location)) +
  geom_point(data = blk_tbl, aes(y = 100 * frac, size = n),
             alpha = 0.25, position = position_nudge(x = 0.22),
             inherit.aes = TRUE, show.legend = FALSE) +
  geom_errorbar(aes(ymin = 100 * lo, ymax = 100 * hi), width = 0.08,
                linewidth = 0.6) +
  geom_point(size = 2.6) +
  facet_wrap(~panel) +
  scale_colour_manual(values = pal, guide = "none") +
  scale_size_area(max_size = 3) +
  ggthm + labs(x = NULL, y = "windows flagged (%)")

# Phenology raster: day-of-year across, year down, one lane per station. Shows
# exposure, cross-station alignment and season in one panel, without stacking
# time courses. Anything the lag profile could have said is visible here as
# whether the two lanes light up at the same x within a row.
ras <- wv |>
  mutate(year = year(datetime), doy = yday(datetime),
         lane = year + ifelse(location == "BM1", -0.19, 0.19))

p_ras <- ggplot(ras, aes(doy, lane)) +
  geom_tile(fill = "grey86", height = 0.34, width = 1) +
  geom_tile(data = filter(ras, rejected), aes(fill = location),
            height = 0.34, width = 1) +
  scale_fill_manual(values = pal) +
  scale_x_continuous(breaks = cumsum(c(1, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30)),
                     labels = month.abb, expand = c(0.01, 0)) +
  scale_y_reverse(breaks = 2020:2025) +
  ggthm + theme(legend.position = "bottom", legend.title = element_blank(),
                panel.grid.major.y = element_blank()) +
  labs(x = NULL, y = NULL,
       subtitle = "upper lane BM1, lower lane BS1; grey = under test")

# Monthly rate, with the exposure that produced it. BM1's May and June rates
# rest on 73 and 46 windows against 438 in April and 452 in September, so the
# point sizes are load-bearing: the tall thin months are the unreliable ones.
p_seas <- seas |>
  ggplot(aes(month, 100 * frac, colour = location)) +
  geom_hline(data = rate_tbl, aes(yintercept = 100 * rate, colour = location),
             linetype = "dashed", linewidth = 0.3) +
  geom_line(linewidth = 0.4, alpha = 0.6) +
  geom_point(aes(size = n)) +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  scale_colour_manual(values = pal, guide = "none") +
  scale_size_area(max_size = 5, name = "windows under test",
                  breaks = c(50, 200, 500)) +
  ggthm + theme(legend.position = "bottom") +
  labs(x = NULL, y = "windows flagged (%)")

plt <- (p_rate | p_seas) / p_ras + plot_layout(heights = c(1, 1.15)) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")
ggsave(file.path(img_out, "station-comparison.png"), plt,
       width = 12, height = 7, dpi = 200)

write_csv(
  bind_rows(
    rate_tbl |> transmute(quantity = "rate_all", location, n_valid, n_rej,
                          est = rate, lo, hi),
    conc_tbl |> transmute(quantity = "rate_concurrent", location,
                          n_valid = nrow(grid), n_rej = NA_integer_,
                          est = rate, lo, hi),
    seas_test |> transmute(quantity = "season_harmonic_amp", location,
                           n_valid = NA_integer_, n_rej = NA_integer_,
                           est = amp, lo = NA_real_, hi = p),
    tibble(quantity = "alignment_lift_lag0", location = "both",
           n_valid = sum(!is.na(a0) & !is.na(b0)), n_rej = n_both,
           est = obs_lift, lo = NA_real_, hi = p_align)),
  file.path(tbl_out, "station-comparison.csv"))

cat("\nWrote", file.path(img_out, "station-comparison.png"), "and the table.\n")
