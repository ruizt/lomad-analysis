library(tidyverse)
library(lomad)
source('mb-analysis/utils.R')   # presmooth_tidal()

# output directory for figures
img_out <- 'mb-analysis/_img'
fs::dir_create(img_out)

# ggplot theming
ggthm <- theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = 'black', linewidth = 0.1))

# --- Data -------------------------------------------------------------------

block_data <- read_csv('_mb-data/ph_o2_blocks.csv')

block_data |>
  count(location, block_id) |>
  mutate(days = n / 24) |>
  arrange(location, block_id)

# Split by location
location_data <- split(block_data, block_data$location)

# --- Presmooth ---------------------------------------------------------------
# Remove tidal periodicity (~25h cycle at hourly) and downsample to 6h

loc_results <- lapply(names(location_data), function(loc) {
  loc_dat <- location_data[[loc]]

  blocks_presm <- loc_dat |>
    group_split(block_id) |>
    lapply(presmooth_tidal, cols = c('o2', 'ph'), step = '6 hours')

  names(blocks_presm) <- vapply(
    blocks_presm, \(b) as.character(b$block_id[[1]]), character(1)
  )

  list(loc = loc, blocks_presm = blocks_presm)
})
names(loc_results) <- names(location_data)

# --- Fit per block, test globally --------------------------------------------
# Two stages. Stage 1 fits every block at both locations and computes raw
# pointwise p-values. Stage 2 pools the raw p-values across ALL blocks and
# applies one Benjamini-Yekutieli step-up correction to the pooled family, so
# FDR is controlled over every window tested in the study, not per block.
# Decisions map back to each block through the realized threshold p_star:
# for a step-up procedure the rejection set is exactly {p_raw <= p_star} with
# p_star = alpha * k / (M * c(M)), which keeps lomad_plot()'s tolerance band
# consistent with the flags.

h_win <- 4
s_win <- 60
alpha <- 0.05

# Minimum length: m >= 1.5s, i.e. n >= 2.5s + h. Since m = n - s - h + 2, a
# floor of k*s + h is a floor of m/s ~ k - 1. The bare fitting minimum
# (n >= 2h + s) admits blocks whose every window overlaps every other one --
# roughly one effective look -- so some floor is needed.
#
# Where to put it is settled by simulation rather than judgement, because the
# threshold and the result move together: the blocks in the 38-46 day band are
# almost all BM1, so the choice preferentially re-arms the station the analysis
# is about. simulations/calibration/blocklength-calibration.R runs the global
# null (d = 0) at this operating point and pools exactly the block ensemble each
# threshold admits. Global-null FDR against a nominal 0.05:
#
#     k = 4.0  (m/s 3.03, 11 blocks)   0.005
#     k = 3.0  (m/s 2.03, 15 blocks)   0.008
#     k = 2.5  (m/s 1.53, 26 blocks)   0.013   <- here
#     k = 2.0  (m/s 1.03, 29 blocks)   0.018
#
# Control holds throughout with room to spare; per-block p-values stay
# conservative at every length (P(p <= 0.05) = 0.020 at k = 2.5 against 0.017
# at k = 3). Shorter blocks do erode the margin monotonically, which is why
# this stops at 2.5 rather than 2: the extra half-step buys only 3 blocks and
# 2,594 hours, having already recovered 11 blocks and 11,292.
#
# Note the power sweep fixed n = 4s in every cell, so k = 4 is the only row
# above that was ever covered by the main simulations; the rest is why this
# study exists.
min_len <- as.integer(2.5 * s_win) + h_win

for (loc in names(loc_results)) {
  blocks_presm <- loc_results[[loc]]$blocks_presm
  blocks <- lapply(blocks_presm, \(b) list(x1 = b$o2, x2 = b$ph))

  too_short <- vapply(blocks, \(b) length(b$x1) < min_len, logical(1))
  if (any(too_short))
    cat(sprintf('  %s: dropping %d blocks shorter than %d observations (m < 1.5s)\n',
                loc, sum(too_short), min_len))
  blocks       <- blocks[!too_short]
  blocks_presm <- blocks_presm[!too_short]

  # Fit each block (variogram-based AR(1) noise); per-block lomad_test() is
  # only a container for the raw p-values here -- its decisions are
  # overwritten by the global stage below.
  block_fits <- lapply(names(blocks), \(nm) {
    b <- blocks[[nm]]
    tryCatch({
      fit <- lomad_fit(b$x1, b$x2, h = h_win, s = s_win)
      tst <- lomad_test(fit, alpha = alpha)
      list(fit = fit, tst = tst)
    }, error = function(e) { message(paste("  Block", nm, "failed:", e$message)); NULL })
  })
  names(block_fits) <- names(blocks)
  block_fits <- Filter(Negate(is.null), block_fits)

  cat(sprintf('=== %s === %d / %d blocks fit\n', loc, length(block_fits), length(blocks)))

  loc_results[[loc]]$blocks       <- blocks
  loc_results[[loc]]$blocks_presm <- blocks_presm
  loc_results[[loc]]$block_fits   <- block_fits
}

# --- Global BY across all blocks ---------------------------------------------

pooled <- bind_rows(lapply(names(loc_results), function(loc) {
  bf <- loc_results[[loc]]$block_fits
  if (is.null(bf) || length(bf) == 0) return(NULL)
  bind_rows(lapply(names(bf), function(nm) {
    fit <- bf[[nm]]$fit; tst <- bf[[nm]]$tst
    tibble(loc = loc, blk = nm, idx = fit$valid_idx,
           p = tst$p_values[fit$valid_idx])
  }))
}))

M      <- nrow(pooled)
p_adj  <- p.adjust(pooled$p, method = 'BY')
k      <- sum(p_adj <= alpha)
c_M    <- sum(1 / seq_len(M))
p_star <- alpha * max(k, 1L) / (M * c_M)

cat(sprintf('\nGlobal BY: M = %d pooled tests over %d blocks; %d rejected (p_star = %.3g)\n',
            M, n_distinct(paste(pooled$loc, pooled$blk)), k, p_star))

# Map global decisions back into each block's tst so every downstream
# consumer (plots, summaries) reflects the global correction.
pooled$p_adj <- p_adj
for (loc in names(loc_results)) {
  bf <- loc_results[[loc]]$block_fits
  for (nm in names(bf)) {
    fit <- bf[[nm]]$fit
    rows <- pooled$loc == loc & pooled$blk == nm
    stopifnot(sum(rows) == length(fit$valid_idx))
    tst <- bf[[nm]]$tst
    tst$p_adj[fit$valid_idx]    <- pooled$p_adj[rows]
    tst$rejected[fit$valid_idx] <- pooled$p[rows] <= p_star
    tst$alpha_eff               <- p_star
    loc_results[[loc]]$block_fits[[nm]]$tst <- tst
  }
}

# --- Rejection summary by location ------------------------------------------

rejection_summary <- bind_rows(lapply(names(loc_results), function(loc) {
  bf <- loc_results[[loc]]$block_fits
  if (is.null(bf)) return(NULL)
  bind_rows(lapply(names(bf), function(nm) {
    tst <- bf[[nm]]$tst
    n_valid    <- sum(!is.na(tst$rejected))
    n_rejected <- sum(tst$rejected, na.rm = TRUE)
    data.frame(location = loc, block_id = nm,
               n_valid = n_valid, n_rejected = n_rejected,
               frac_rejected = n_rejected / n_valid,
               detected = n_rejected > 0)
  }))
}))

cat('\n========== Rejection summary ==========\n')
print(rejection_summary, row.names = FALSE)

rejection_summary |>
  group_by(location) |>
  summarise(n_blocks = n(),
            blocks_detected = sum(detected),
            total_valid = sum(n_valid),
            total_rejected = sum(n_rejected),
            frac_rejected = total_rejected / total_valid,
            .groups = "drop") |>
  print()

# --- ggplot lomad-fit visualizations -----------------------------------------

library(patchwork)

make_lomad_plot_data <- function(loc_name, results) {
  loc <- results[[loc_name]]
  bf  <- loc$block_fits

  block_dfs    <- list()
  shade_list   <- list()
  trigger_list <- list()

  add_shade <- function(flag_vec, dates, bid, panel_label) {
    r  <- rle(flag_vec)
    en <- cumsum(r$lengths)
    st <- en - r$lengths + 1L
    ds <- st[r$values]; de <- en[r$values]
    if (length(ds) == 0) return(NULL)
    tibble(block_id = bid, xmin = dates[ds], xmax = dates[de], panel = panel_label)
  }

  for (nm in names(bf)) {
    presm <- loc$blocks_presm[[nm]]
    fit   <- bf[[nm]]$fit
    tst   <- bf[[nm]]$tst
    n     <- nrow(presm)
    dates <- presm$datetime
    bid   <- presm$block_id[[1]]
    h     <- fit$inputs$h

    rejected <- tst$rejected
    rejected[is.na(rejected)] <- FALSE

    block_dfs[[nm]] <- tibble(
      datetime = dates,
      block_id = bid,
      o2       = loc$blocks[[nm]]$x1,
      ph       = loc$blocks[[nm]]$x2,
      ma1      = fit$ma1,
      ma2      = fit$ma2,
      R        = fit$R,
      rho      = fit$rho,
      rejected = rejected
    )

    # Back-shifted shading for upper panel
    rej_shifted <- rep(FALSE, n)
    for (t in which(rejected))
      rej_shifted[max(1L, t - (h - 1L)):t] <- TRUE

    shade_list[[paste0(nm, '_up')]] <- add_shade(rej_shifted, dates, bid, 'upper')
    shade_list[[paste0(nm, '_lo')]] <- add_shade(rejected, dates, bid, 'lower')

    # Trigger points
    r_rej     <- rle(rejected)
    en_rej    <- cumsum(r_rej$lengths)
    st_rej    <- en_rej - r_rej$lengths + 1L
    entry_pos <- st_rej[r_rej$values]
    trig_t    <- pmax(1L, entry_pos - (h - 1L))
    if (length(trig_t) > 0)
      trigger_list[[nm]] <- tibble(block_id = bid, datetime = dates[trig_t])
  }

  list(
    main     = bind_rows(block_dfs),
    shade    = bind_rows(shade_list),
    triggers = bind_rows(trigger_list)
  )
}

make_lomad_ggplot <- function(pd) {
  shade_up  <- pd$shade |> filter(panel == 'upper')
  shade_lo  <- pd$shade |> filter(panel == 'lower')
  shade_col <- rgb(0.7, 0.85, 1, 0.4)

  p_up <- ggplot(pd$main, aes(x = datetime)) +
    geom_rect(data = shade_up, inherit.aes = FALSE,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              fill = shade_col) +
    geom_vline(data = pd$triggers, aes(xintercept = datetime),
               color = 'grey50', linetype = 'dashed', linewidth = 0.5) +
    geom_line(aes(y = o2), color = rgb(0, 0, 1, 0.2)) +
    geom_line(aes(y = ph), color = rgb(1, 0, 0, 0.2)) +
    geom_line(aes(y = ma1), color = 'blue', linewidth = 0.3) +
    geom_line(aes(y = ma2), color = 'red',  linewidth = 0.3) +
    facet_grid(~block_id, scales = 'free_x', space = 'free_x') +
    scale_x_datetime(breaks = function(x) mean(x), date_labels = '%b %Y') +
    ggthm +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          strip.text   = element_blank()) +
    labs(x = NULL, y = 'series')

  p_lo <- ggplot(pd$main, aes(x = datetime)) +
    geom_rect(data = shade_lo, inherit.aes = FALSE,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              fill = shade_col) +
    geom_hline(yintercept = 0, color = 'grey80', linewidth = 0.3) +
    geom_line(aes(y = R),   color = 'grey40') +
    geom_line(aes(y = rho), color = 'grey30', linetype = 'dashed') +
    facet_grid(~block_id, scales = 'free_x', space = 'free_x') +
    scale_x_datetime(breaks = function(x) mean(x), date_labels = '%b %Y') +
    ggthm +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          strip.text  = element_blank()) +
    labs(x = NULL, y = 'correlation')

  p_up / p_lo + plot_layout(heights = c(2, 1))
}

# Generate ggplots for each location
for (loc in names(loc_results)) {
  if (is.null(loc_results[[loc]]$block_fits)) next
  pd  <- make_lomad_plot_data(loc, loc_results)
  plt <- make_lomad_ggplot(pd)
  ggsave(paste0(img_out, '/', loc, '-lomad-fit.png'), plt, width = 12, height = 3)
  print(plt)
}


# --- Both stations on one calendar axis ---------------------------------------
# The per-location figures above lay blocks out side by side with free x scales,
# which compresses gaps and makes the two stations impossible to compare: block
# 7 at BM1 and block 7 at BS1 are different periods. Here x is real time,
# shared, so a rejection at one station can be read against what the other was
# doing at that moment. Blocks appear as segments with gaps between them, which
# is what the record actually looks like.

aligned <- bind_rows(lapply(names(loc_results), function(loc) {
  bf <- loc_results[[loc]]$block_fits
  if (is.null(bf) || length(bf) == 0) return(NULL)
  bind_rows(lapply(names(bf), function(nm) {
    presm <- loc_results[[loc]]$blocks_presm[[nm]]
    fit   <- bf[[nm]]$fit
    tst   <- bf[[nm]]$tst
    tibble(location = loc,
           block_id = presm$block_id[[1]],
           datetime = presm$datetime,
           o2       = loc_results[[loc]]$blocks[[nm]]$x1,
           ph       = loc_results[[loc]]$blocks[[nm]]$x2,
           ma1      = fit$ma1,
           ma2      = fit$ma2,
           R        = fit$R,
           rho      = fit$rho,
           rejected = replace_na(tst$rejected, FALSE))
  }))
}))

# Flagged runs as rectangles, so shading does not break at every observation
episodes <- aligned |>
  group_by(location, block_id) |>
  group_modify(~{
    r  <- rle(.x$rejected)
    en <- cumsum(r$lengths); st <- en - r$lengths + 1L
    if (!any(r$values)) return(tibble(xmin = as.POSIXct(character()),
                                      xmax = as.POSIXct(character())))
    tibble(xmin = .x$datetime[st[r$values]], xmax = .x$datetime[en[r$values]])
  }) |>
  ungroup()

# Facet on shared coverage windows rather than on blocks. Blocks cannot be
# faceted directly and still align: BM1 block 7 and BS1 block 7 are unrelated
# periods. Instead, merge the two stations' block extents into common intervals
# -- a window runs from where either station starts recording to where both
# have stopped for longer than `gap_days` -- so within a facet the two rows
# share an x range and are directly comparable. `space = "free_x"` makes panel
# width proportional to duration, so time is to scale within and across facets,
# while the dead stretches between windows are dropped.

gap_days <- 30

spans <- aligned |>
  group_by(location, block_id) |>
  summarise(start = min(datetime), end = max(datetime), .groups = "drop") |>
  arrange(start)

# merge overlapping-or-close spans across BOTH stations into shared windows
win <- spans |>
  mutate(new = start > lag(cummax(as.numeric(end)), default = -Inf) +
                 gap_days * 86400,
         window = cumsum(replace_na(new, TRUE))) |>
  group_by(window) |>
  summarise(wstart = min(start), wend = max(end), .groups = "drop") |>
  mutate(label = paste(format(wstart, "%b %Y"), format(wend, "%b %Y"), sep = " - "))

assign_window <- function(x) {
  i <- vapply(x, function(d) which(d >= win$wstart & d <= win$wend)[1], integer(1))
  factor(win$label[i], levels = win$label)
}
aligned  <- aligned  |> mutate(window = assign_window(datetime))
episodes <- episodes |> mutate(window = assign_window(xmin))

cat(sprintf("\n%d shared coverage windows (gap threshold %d days)\n",
            nrow(win), gap_days))

shade_col <- rgb(0.7, 0.85, 1, 0.5)

plt_aligned <- ggplot(aligned, aes(x = datetime, group = block_id)) +
  geom_rect(data = episodes, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = shade_col) +
  geom_line(aes(y = ma1), colour = "blue", linewidth = 0.3) +
  geom_line(aes(y = ma2), colour = "red",  linewidth = 0.3) +
  facet_grid(location ~ window, scales = "free_x", space = "free_x") +
  scale_x_datetime(date_breaks = "3 months", date_labels = "%b %Y",
                   expand = expansion(mult = 0.03)) +
  ggthm +
  theme(panel.spacing.x = unit(2, "pt"),
        strip.text.x = element_blank(),
        axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7)) +
  labs(x = NULL, y = "smoothed series")

ggsave(paste0(img_out, "/stations-aligned.png"), plt_aligned,
       width = 14, height = 4.2, dpi = 200)
print(plt_aligned)

# --- Persist the fitted window-level output -----------------------------------
# One row per analysed window, carrying the test decision under the global
# threshold. Kept for ad hoc work; everything below runs from `wv` directly.

saveRDS(aligned, "_mb-data/lomad_windows.rds")
message("Wrote _mb-data/lomad_windows.rds (", nrow(aligned), " windows)")


# =============================================================================
# Station comparison
# =============================================================================
# Two claims: the stations differ greatly in how much of their record gets
# flagged, and they nonetheless share a seasonal pattern.
#
# Every interval below is a bootstrap over BLOCKS, never over windows. At
# s = 60 on a 6-hourly step, successive windows share 59/60 of their data and
# rejections arrive in runs, so a pointwise standard error would be fiction.
# For the same reason the seasonality nulls rotate each block's flag vector
# circularly within that block: run lengths, the marginal rate and the exposure
# pattern are all held fixed, and only the calendar alignment is destroyed.

set.seed(20260727)
B <- 3000L

wv <- aligned |>
  # `rejected` was NA-filled to FALSE upstream; a window is only under test
  # where the fit produced both a correlation and a benchmark for it
  mutate(valid = !is.na(R) & !is.na(rho), blk = paste(location, block_id)) |>
  filter(valid)

tbl_out <- 'mb-analysis/_tbl'
fs::dir_create(tbl_out)

# --- Rate ---------------------------------------------------------------------

boot_rate <- function(d, B) {
  blocks <- split(d, d$blk)
  vapply(seq_len(B), function(i) {
    s <- blocks[sample.int(length(blocks), replace = TRUE)]
    sum(vapply(s, \(b) sum(b$rejected), numeric(1))) /
      sum(vapply(s, \(b) nrow(b), numeric(1)))
  }, numeric(1))
}

rate_tbl <- wv |>
  group_by(location) |>
  group_modify(~{
    bs <- boot_rate(.x, B)
    tibble(n_blocks = n_distinct(.x$blk), n_valid = nrow(.x),
           n_rej = sum(.x$rejected), rate = sum(.x$rejected) / nrow(.x),
           lo = quantile(bs, 0.025), hi = quantile(bs, 0.975))
  }) |> ungroup()

cat('\n========== Rate by station ==========\n')
print(rate_tbl |> mutate(across(c(rate, lo, hi), \(x) round(100 * x, 1))))

rr <- boot_rate(filter(wv, location == 'BM1'), B) /
      boot_rate(filter(wv, location == 'BS1'), B)
cat(sprintf('rate ratio BM1/BS1 = %.2f  (95%% CI %.2f to %.2f)\n',
            rate_tbl$rate[1] / rate_tbl$rate[2],
            quantile(rr, 0.025), quantile(rr, 0.975)))

# The two stations were up at different times, so the headline ratio could in
# principle be an exposure artifact. Restricted to timestamps where both were
# under test, it is not.
conc <- wv |>
  select(location, datetime, rejected) |>
  pivot_wider(id_cols = datetime, names_from = location,
              values_from = rejected, values_fn = any) |>
  drop_na(BM1, BS1)
cat(sprintf('concurrent support: %d timestamps (%.0f days) -- BM1 %.1f%%, BS1 %.1f%%\n',
            nrow(conc), nrow(conc) * 6 / 24,
            100 * mean(conc$BM1), 100 * mean(conc$BS1)))

# Leave-one-block-out. The cleanest robustness statement available, because it
# makes no distributional assumption at all: if the two ranges do not overlap,
# no single block is carrying the difference.
blk_tbl <- wv |> group_by(location, blk) |>
  summarise(n = n(), rej = sum(rejected), frac = rej / n, .groups = 'drop')
for (loc in c('BM1', 'BS1')) {
  d  <- filter(blk_tbl, location == loc)
  lo <- vapply(seq_len(nrow(d)), \(i) sum(d$rej[-i]) / sum(d$n[-i]), numeric(1))
  cat(sprintf('%s leave-one-block-out: %.1f%% to %.1f%%\n',
              loc, 100 * min(lo), 100 * max(lo)))
}

# The conservative counterweight: a block-level test that ignores within-block
# structure entirely. It does not reach significance, and the write-up should
# say so -- BS1's whole signal is two blocks.
det <- blk_tbl |> group_by(location) |>
  summarise(blocks = n(), detecting = sum(rej > 0), .groups = 'drop')
ft <- fisher.test(matrix(c(det$detecting, det$blocks - det$detecting), nrow = 2))
cat(sprintf('blocks detecting: BM1 %d/%d, BS1 %d/%d (Fisher p = %.3f)\n',
            det$detecting[1], det$blocks[1],
            det$detecting[2], det$blocks[2], ft$p.value))

# --- Seasonality --------------------------------------------------------------

rotate_within_blocks <- function(r, blk) {
  out <- r
  for (b in split(seq_along(r), blk))
    if (length(b) > 1L)
      out[b] <- r[b][(seq_along(b) + sample.int(length(b), 1L) - 1L) %% length(b) + 1L]
  out
}

# Spring contrast rather than a harmonic. A single sinusoid is the wrong shape
# here: BM1 has two peaks about five months apart, so the fitted annual maximum
# lands in July, between them, and the test loses most of its power (BM1
# p = 0.046 on harmonic 1 but its estimated peak is meaningless; pooled, the
# annual harmonic gives p = 0.27). The semiannual harmonic does better at BM1
# (p = 0.026) but has no interpretation at BS1, which has only one peak.
#
# What the two stations actually share is a window, not a waveform: both are
# silent in midwinter and midsummer and both concentrate their flags in March
# to June. A one-degree-of-freedom contrast on that window is the powerful test.
spring_contrast <- function(d, B) {
  s    <- month(d$datetime) %in% 3:6
  stat <- function(r) mean(r[s]) - mean(r[!s])
  obs  <- stat(d$rejected)
  nul  <- vapply(seq_len(B),
                 \(i) stat(rotate_within_blocks(d$rejected, d$blk)), numeric(1))
  tibble(spring = 100 * mean(d$rejected[s]), rest = 100 * mean(d$rejected[!s]),
         diff = 100 * obs, p = (1 + sum(nul >= obs)) / (1 + B))
}

cat('\n========== Seasonality: March-June vs rest of year ==========\n')
seas_test <- bind_rows(
  wv |> group_by(location) |> group_modify(~spring_contrast(.x, B)) |> ungroup(),
  spring_contrast(wv, B) |> mutate(location = 'pooled')) |>
  relocate(location)
print(seas_test |> mutate(across(c(spring, rest, diff), \(x) round(x, 1)),
                          p = round(p, 4)))

# Pooling is what makes this work, and the reason is worth recording: BM1 alone
# has a wide rotation null (its 13.7% flagged fraction and second, autumn peak
# mean random rotations land in spring almost as often as the data do), while
# BS1 alone is too sparse. Together the shared window clears.
cat('\nshare of each station\'s flags falling in March-June:\n')
print(wv |> filter(rejected) |> mutate(sp = month(datetime) %in% 3:6) |>
      group_by(location) |>
      summarise(spring = sum(sp), other = sum(!sp),
                pct_spring = round(100 * mean(sp)), .groups = 'drop'))

seas <- wv |> mutate(month = month(datetime)) |>
  group_by(location, month) |>
  summarise(n = n(), rej = sum(rejected), frac = rej / n, .groups = 'drop')

cat('\nmonths with no detection at either station: ')
cat(paste(month.abb[setdiff(1:12, unique(seas$month[seas$rej > 0]))],
          collapse = ', '), '\n')

write_csv(bind_rows(
  rate_tbl |> transmute(quantity = 'rate', location, n = n_valid, rej = n_rej,
                        est = rate, lo, hi),
  seas_test |> transmute(quantity = 'spring_contrast', location, n = NA_integer_,
                         rej = NA_integer_, est = diff, lo = NA_real_, hi = p)),
  file.path(tbl_out, 'station-comparison.csv'))

# --- Figure -------------------------------------------------------------------

pal <- c(BM1 = '#C44E52', BS1 = '#4C72B0')

# Phenology raster: day of year across, year down, one lane per station, grey
# where under test. This is the evidence that the spring concentration recurs
# rather than resting on one episode, and it shows the exposure that produced
# every rate below it.
ras <- wv |>
  mutate(year = year(datetime), doy = yday(datetime),
         lane = year + ifelse(location == 'BM1', -0.19, 0.19))

p_ras <- ggplot(ras, aes(doy, lane)) +
  geom_tile(fill = 'grey86', height = 0.34, width = 1) +
  geom_tile(data = filter(ras, rejected), aes(fill = location),
            height = 0.34, width = 1) +
  scale_fill_manual(values = pal) +
  scale_x_continuous(
    breaks = cumsum(c(1, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30)),
    labels = month.abb, expand = c(0.01, 0)) +
  scale_y_reverse(breaks = 2020:2025) +
  ggthm + theme(legend.position = 'bottom', legend.title = element_blank(),
                panel.grid.major.y = element_blank()) +
  labs(x = NULL, y = NULL,
       subtitle = 'upper lane BM1, lower lane BS1; grey = under test')

# Point size is load-bearing: BM1's May and June rates rest on 73 and 46
# windows against 438 in April and 452 in September.
p_seas <- ggplot(seas, aes(month, 100 * frac, colour = location)) +
  geom_hline(data = rate_tbl, aes(yintercept = 100 * rate, colour = location),
             linetype = 'dashed', linewidth = 0.3) +
  geom_line(linewidth = 0.4, alpha = 0.6) +
  geom_point(aes(size = n)) +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  scale_colour_manual(values = pal, guide = 'none') +
  scale_size_area(max_size = 5, name = 'windows under test',
                  breaks = c(50, 200, 500)) +
  ggthm + theme(legend.position = 'bottom') +
  labs(x = NULL, y = 'windows flagged (%)')

plt_seas <- p_ras / p_seas + plot_layout(heights = c(1.25, 1)) +
  plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_suffix = ')')
ggsave(paste0(img_out, '/seasonality.png'), plt_seas,
       width = 10, height = 7, dpi = 200)
print(plt_seas)
