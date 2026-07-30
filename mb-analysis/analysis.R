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
# Fit every block, then apply one Benjamini-Yekutieli step-up correction to the
# pooled p-values so FDR is controlled study-wide rather than per block.
# Decisions map back through p_star = alpha * k / (M * c(M)), the realized
# threshold: for a step-up procedure the rejection set is {p_raw <= p_star}.

h_win <- 4
s_win <- 60
alpha <- 0.05

# Minimum block length, m >= 1.5s. This is below the regime the power sweep
# covered (n = 4s throughout), so it was checked against the global null in
# simulations/calibration/blocklength-calibration.R: FDR 0.013 here, 0.008 at
# the previous m >= 2s, against a nominal 0.05.
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

  # Per-block lomad_test() is only a container for the raw p-values; its
  # decisions are overwritten by the global stage below.
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

# Returns the two panels rather than a composed plot so both stations can stack.
STATION_NAME <- c(BM1 = 'Bay Mouth (BM)', BS1 = 'Bay Head (BH)')

make_lomad_ggplot <- function(pd, loc) {
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
          strip.text   = element_blank(),
          axis.title.y = element_text(size = 8),
          axis.text.y  = element_text(size = 7),
          plot.title   = element_text(size = 9, face = 'plain')) +
    labs(x = NULL, y = 'moving averages', title = STATION_NAME[[loc]])

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
    theme(axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1,
                                      size = 6),
          axis.title.y = element_text(size = 8),
          axis.text.y  = element_text(size = 7),
          strip.text   = element_blank()) +
    labs(x = NULL, y = 'correlation')

  list(up = p_up, lo = p_lo)
}

# Both stations stacked, series over correlation. Independent x scales: the
# stations have different blocks, so this is a per-block view.
panels <- unlist(lapply(names(loc_results), function(loc) {
  if (is.null(loc_results[[loc]]$block_fits)) return(NULL)
  make_lomad_ggplot(make_lomad_plot_data(loc, loc_results), loc)
}), recursive = FALSE)

plt_fits <- wrap_plots(panels, ncol = 1) +
  plot_layout(heights = rep(c(2, 1), length(panels) / 2))
ggsave(paste0(img_out, '/sfig-mb-detections.png'), plt_fits,
       width = 10, height = 5, dpi = 200)
print(plt_fits)


# --- Window-level output ------------------------------------------------------
# One row per analysed window; input to the station comparison below.

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


# --- Persist the fitted window-level output -----------------------------------
# Kept for ad hoc work; everything below runs from `aligned` directly.

saveRDS(aligned, "_mb-data/lomad_windows.rds")
message("Wrote _mb-data/lomad_windows.rds (", nrow(aligned), " windows)")


# =============================================================================
# Station comparison
# =============================================================================
# Descriptive only. A seasonal hypothesis would have to be specified after
# seeing these data, and a rate interval would be inference on inference. What
# follows the point estimates is sensitivity, not uncertainty.

wv <- aligned |>
  # `rejected` was NA-filled to FALSE upstream; a window is only under test
  # where the fit produced both a correlation and a benchmark for it
  mutate(valid = !is.na(R) & !is.na(rho), blk = paste(location, block_id)) |>
  filter(valid)

tbl_out <- 'mb-analysis/_tbl'
fs::dir_create(tbl_out)

# --- Rate ---------------------------------------------------------------------

rate_tbl <- wv |>
  group_by(location) |>
  summarise(n_blocks = n_distinct(blk), n_valid = n(),
            n_rej = sum(rejected), rate = n_rej / n_valid, .groups = 'drop')

cat('\n========== Rate by station ==========\n')
print(rate_tbl |> mutate(rate = round(100 * rate, 1)))
cat(sprintf('rate ratio BM1/BS1 = %.2f\n',
            rate_tbl$rate[rate_tbl$location == 'BM1'] /
            rate_tbl$rate[rate_tbl$location == 'BS1']))

# --- Sensitivity of the ratio -------------------------------------------------
# Does the ratio depend on any one block, or on the stations having been up at
# different times?

ratio_of <- function(d) {
  r <- tapply(d$rejected, d$location, mean)
  unname(r['BM1'] / r['BS1'])
}

# (i) drop each block in turn
blk_tbl <- wv |> group_by(location, blk) |>
  summarise(n = n(), rej = sum(rejected), frac = rej / n, .groups = 'drop')

loo <- vapply(unique(wv$blk), \(b) ratio_of(filter(wv, blk != b)), numeric(1))
cat(sprintf('\nleave-one-block-out ratio: %.2f to %.2f over %d blocks\n',
            min(loo), max(loo), length(loo)))
for (loc in c('BM1', 'BS1')) {
  d  <- filter(blk_tbl, location == loc)
  lo <- vapply(seq_len(nrow(d)), \(i) sum(d$rej[-i]) / sum(d$n[-i]), numeric(1))
  cat(sprintf('  %s rate spans %.1f%% to %.1f%%\n',
              loc, 100 * min(lo), 100 * max(lo)))
}

# (ii) restrict to timestamps where both stations were under test, which removes
# differing exposure as an explanation
conc <- wv |>
  select(location, datetime, rejected) |>
  pivot_wider(id_cols = datetime, names_from = location,
              values_from = rejected, values_fn = any) |>
  drop_na(BM1, BS1)
cat(sprintf('concurrent support (%d timestamps, %.0f days): BM1 %.1f%%, BS1 %.1f%%, ratio %.2f\n',
            nrow(conc), nrow(conc) * 6 / 24,
            100 * mean(conc$BM1), 100 * mean(conc$BS1),
            mean(conc$BM1) / mean(conc$BS1)))

# --- Season -------------------------------------------------------------------

seas <- wv |> mutate(month = month(datetime)) |>
  group_by(month) |>
  summarise(n = n(), rej = sum(rejected), frac = rej / n, .groups = 'drop')

seas_loc <- wv |> mutate(month = month(datetime)) |>
  group_by(location, month) |>
  summarise(n = n(), rej = sum(rejected), .groups = 'drop')

cat('\n========== Detections by month, both stations pooled ==========\n')
print(seas |> mutate(pct = round(100 * frac, 1)) |>
      transmute(month = month.abb[month], n, rej, pct), n = 12)

cat('\nmonths with no detection at either station: ')
cat(paste(month.abb[setdiff(1:12, unique(seas_loc$month[seas_loc$rej > 0]))],
          collapse = ', '), '\n')
cat('share of flags in Mar-Jun: ')
cat(paste(wv |> filter(rejected) |> group_by(location) |>
          summarise(p = sprintf('%s %.0f%%', location[1], 100 * mean(month(datetime) %in% 3:6)),
                    .groups = 'drop') |> pull(p), collapse = ', '), '\n')

write_csv(
  rate_tbl |> transmute(quantity = 'rate', location, n = n_valid, rej = n_rej,
                        est = rate),
  file.path(tbl_out, 'station-rates.csv'))

# --- Figure -------------------------------------------------------------------

# Phenology raster: day of year across, year down, upper lane BM, lower lane BS.
# Grey is every window the test reached a decision on, which excludes the
# leading s + h - 2 points of each block.
pal <- c(BM = '#C44E52', BH = '#4C72B0')

# Shared x geometry. Tiles are centred on integer doy with width 1, so a
# month's visual centre sits half a day right of its arithmetic midpoint.
MONTH_END   <- cumsum(c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31))
MONTH_START <- c(0, head(MONTH_END, -1))
MONTH_MID   <- (MONTH_START + MONTH_END + 1) / 2
X_EXPAND    <- expansion(mult = c(0.045, 0.01))   # gutter holds the lane labels
X_LIM       <- c(-2, 366)   # shared, not just breaks: expand() is relative to
                            # each panel's own range, which differs

ras <- wv |>
  mutate(year    = year(datetime), doy = yday(datetime),
         lane    = year + ifelse(location == 'BM1', -0.19, 0.19),
         station = c(BM1 = 'BM', BS1 = 'BH')[location])

lanes <- distinct(ras, lane, station)

p_ras <- ggplot(ras, aes(doy, lane)) +
  geom_tile(fill = 'grey86', height = 0.34, width = 1) +
  geom_tile(data = filter(ras, rejected), aes(fill = station),
            height = 0.34, width = 1) +
  geom_text(data = lanes, aes(x = -2, y = lane, label = station),
            inherit.aes = FALSE, hjust = 1, size = 2.5, colour = 'grey35') +
  scale_fill_manual(values = pal,
                    labels = paste(names(pal), 'detection')) +
  scale_x_continuous(breaks = MONTH_MID, labels = month.abb,
                     limits = X_LIM, oob = scales::oob_keep,
                     expand = X_EXPAND) +
  scale_y_reverse(breaks = 2020:2025) +
  ggthm + theme(legend.position = 'bottom', legend.title = element_blank(),
                panel.grid.major.y = element_blank(),
                axis.ticks.length = unit(0, 'in'),
                panel.border = element_blank(),
                axis.text.y = element_text(angle = 90, hjust = 0.5)) +
  labs(x = NULL, y = "Year")

# Pooled monthly rate as a marginal strip above the raster, on the same axis.
p_seas <- seas |>
  mutate(mid = MONTH_MID[month]) |>
  ggplot(aes(mid, 100 * frac)) +
  geom_line(linewidth = 0.4, colour = 'grey25') +
  geom_point(size = 1.3, colour = 'grey15') +
  scale_x_continuous(breaks = MONTH_MID, labels = month.abb,
                     limits = X_LIM, oob = scales::oob_keep,
                     expand = X_EXPAND) +
  scale_y_continuous(breaks = c(0, 10, 20)) +
  ggthm + theme(axis.text.x  = element_blank(),
                axis.ticks = element_blank(),
                panel.border = element_blank(),
                panel.grid.major.y = element_line(linewidth = 0.1, color = 'darkgrey')) +
  labs(x = NULL, y = 'Detections \n(%)')

plt_seas <- p_seas / p_ras + plot_layout(heights = c(1.5, 5))
ggsave(paste0(img_out, '/fig-mb-seasonality.png'), plt_seas,
       width = 6, height = 4, dpi = 200)
print(plt_seas)
