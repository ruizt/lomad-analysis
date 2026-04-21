library(tidyverse)
devtools::load_all()
source('dev/mb-analysis/utils.R')   # presmooth_tidal()

# output directory for figures
img_out <- '_img'
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

# Split by location; block_id is not globally unique
location_data <- split(block_data, block_data$location)

# --- Presmooth --------------------------------------------------------------
# Remove tidal periodicity (~25 h cycle at hourly resolution) and downsample
# to 6-hour resolution. Applied independently within each block.

loc_results <- lapply(names(location_data), function(loc) {
  loc_dat <- location_data[[loc]]

  blocks_presm <- loc_dat |>
    group_split(block_id) |>
    lapply(presmooth_tidal, cols = c('o2', 'ph'), cycle_len = 25, step = '6 hours')

  names(blocks_presm) <- vapply(
    blocks_presm, \(b) as.character(b$block_id[[1]]), character(1)
  )

  list(loc = loc, blocks_presm = blocks_presm)
})
names(loc_results) <- names(location_data)

# --- Fit null model, separately per location --------------------------------

block_fit_view <- function(fit, k) {
  blk <- fit$blocks[[k]]
  fit$trend_hat <- blk$trend_hat
  fit$ma1       <- blk$ma1
  fit$ma2       <- blk$ma2
  fit$R         <- blk$R
  fit$I         <- blk$I
  fit$valid_idx <- blk$valid_idx   # local indices, correct for this block's length
  fit
}

for (loc in names(loc_results)) {
  blocks_presm <- loc_results[[loc]]$blocks_presm
  blocks       <- lapply(blocks_presm, \(b) list(x1 = b$o2, x2 = b$ph))

  fit <- lomad_fit_blocks(blocks,
                          rho0   = 0,
                          q      = 5 * 4,
                          h      = 10 * 4,
                          max_pq = 2)

  cat('\n===', loc, '===\n')
  print(fit$null_model)
  print(fit$observed)
  print(fit$expected_asymptotic)

  # --- Visualize individual blocks ------------------------------------------

  pdf(paste0(img_out, '/', loc, '-blocks-fit.pdf'), width = 5, height = 4)
  for (nm in names(fit$blocks)) {
    plot_lomad_fit(blocks[[nm]]$x1, blocks[[nm]]$x2,
                   fit   = block_fit_view(fit, nm),
                   alpha = 0.4)
    title(main = paste(loc, 'block', nm), line = 0.5)
  }
  dev.off()

  # --- Inference ------------------------------------------------------------

  # Asymptotic test (rough approximation)
  test_analytic <- lomad_test_analytic(fit)
  cat('\n--- Analytic p-values (', loc, ') ---\n')
  print(test_analytic$p_values)

  # Markov-chain bootstrap (fast)
  test_mc <- lomad_test_mc(fit, B = 1000, seed = 4721,
                           verbose = TRUE)
  cat('\n--- MC bootstrap p-values (', loc, ') ---\n')
  print(test_mc$p_values)

  # Full parametric bootstrap (use parallel cores to speed up)
  test_boot <- lomad_test_boot(fit, B = 1000, seed = 4721,
                               ncores = parallel::detectCores() - 1,
                               verbose = TRUE)
  cat('\n--- Parametric bootstrap p-values (', loc, ') ---\n')
  print(test_boot$p_values)
  print(test_boot$expected)
  print(test_boot$observed)

  # Store results back
  loc_results[[loc]]$blocks <- blocks
  loc_results[[loc]]$fit    <- fit
  loc_results[[loc]]$test_analytic <- test_analytic
  loc_results[[loc]]$test_mc       <- test_mc
  loc_results[[loc]]$test_boot     <- test_boot
}

# --- Summary table: p-values by location ------------------------------------

pval_summary <- function(results, test_slot) {
  stats <- c('entry_rate', 'mean_run_length', 'frac_state', 'n_entries')
  rows <- lapply(names(results), function(loc) {
    pv <- results[[loc]][[test_slot]]$p_values
    data.frame(
      location = loc,
      statistic = stats,
      p_value   = signif(unlist(pv[stats]), 3),
      row.names = NULL
    )
  })
  bind_rows(rows) |>
    pivot_wider(names_from = location, values_from = p_value)
}

cat('\n========== P-value summary ==========\n')
cat('\n-- Analytic --\n');   print(pval_summary(loc_results, 'test_analytic'))
cat('\n-- MC bootstrap --\n'); print(pval_summary(loc_results, 'test_mc'))
cat('\n-- Parametric bootstrap --\n'); print(pval_summary(loc_results, 'test_boot'))

# --- ggplot lomad-fit visualizations (BM1 and BS1) ---------------------------

library(patchwork)

# Assemble per-location data frames from loc_results for ggplot rendering.
# Returns main time-series data, shading intervals, trigger points, and the
# BY-corrected correlation threshold.
make_lomad_plot_data <- function(loc_name, results) {
  loc   <- results[[loc_name]]
  fit   <- loc$fit
  h     <- fit$null_model$h
  x_eff <- fit$thresholds$x_eff

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

  for (nm in names(fit$blocks)) {
    presm <- loc$blocks_presm[[nm]]
    blk   <- fit$blocks[[nm]]
    n     <- nrow(presm)
    dates <- presm$datetime
    I     <- blk$I
    valid <- blk$valid_idx
    bid   <- presm$block_id[[1]]

    block_dfs[[nm]] <- tibble(
      datetime  = dates,
      block_id  = bid,
      o2        = loc$blocks[[nm]]$x1,
      ph        = loc$blocks[[nm]]$x2,
      ma1       = blk$ma1,
      ma2       = blk$ma2,
      trend_hat = blk$trend_hat,
      R         = blk$R,
      I         = I
    )

    # Back-shifted shading for upper panel (expand each detected point t back
    # to [t-(h-1), t] to align with the data window that drove the detection)
    I_shifted <- rep(FALSE, n)
    for (t in which(!is.na(I) & I == 1L))
      I_shifted[max(1L, t - (h - 1L)):t] <- TRUE

    shade_list[[paste0(nm, '_up')]] <- add_shade(I_shifted, dates, bid, 'upper')
    # Lower panel shading: no back-shift, just where R actually crossed threshold
    shade_list[[paste0(nm, '_lo')]] <- add_shade(!is.na(I) & I == 1L, dates, bid, 'lower')

    # Trigger points: left edge of the back-shifted window for each new episode
    I_v       <- I[valid]
    m2        <- length(I_v)
    entry_pos <- which(I_v == 1L & c(0L, I_v[-m2]) != 1L)
    trig_t    <- pmax(1L, valid[entry_pos] - (h - 1L))
    if (length(trig_t) > 0)
      trigger_list[[nm]] <- tibble(block_id = bid, datetime = dates[trig_t])
  }

  list(
    main     = bind_rows(block_dfs),
    shade    = bind_rows(shade_list),
    triggers = bind_rows(trigger_list),
    x_eff    = x_eff
  )
}

make_lomad_ggplot <- function(pd, title_str) {
  shade_up  <- pd$shade |> filter(panel == 'upper')
  shade_lo  <- pd$shade |> filter(panel == 'lower')
  shade_col <- rgb(0.7, 0.85, 1, 0.4)

  # Upper panel: raw series (transparent), smoothed MA lines, shared trend,
  # back-shifted decoupling shading, and dashed trigger-point lines
  p_up <- ggplot(pd$main, aes(x = datetime)) +
    geom_rect(data = shade_up, inherit.aes = FALSE,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              fill = shade_col) +
    geom_vline(data = pd$triggers, aes(xintercept = datetime),
               color = 'grey50', linetype = 'dashed', linewidth = 0.5) +
    geom_line(aes(y = o2),        color = rgb(0, 0, 1, 0.2)) +
    geom_line(aes(y = ph),        color = rgb(1, 0, 0, 0.2)) +
    # geom_line(aes(y = trend_hat), color = rgb(0.4, 0.4, 0.4, 0.8), linewidth = 0.2) +
    geom_line(aes(y = ma1),       color = 'blue', linewidth = 0.3) +
    geom_line(aes(y = ma2),       color = 'red',  linewidth = 0.3) +
    facet_grid(~block_id, scales = 'free_x', space = 'free_x') +
    scale_x_datetime(breaks = function(x) mean(x), date_labels = '%b %Y') +
    ggthm +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          strip.text   = element_text(size = 7)) +
    labs(x = NULL, y = 'series', title = title_str)

  # Lower panel: rolling correlation with BY-corrected threshold
  p_lo <- ggplot(pd$main, aes(x = datetime, y = R)) +
    geom_rect(data = shade_lo, inherit.aes = FALSE,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              fill = shade_col) +
    geom_hline(yintercept = 0,        color = 'grey80', linewidth = 0.3) +
    geom_hline(yintercept = pd$x_eff, color = 'grey30', linetype = 'dashed') +
    geom_line(color = 'grey40') +
    facet_grid(~block_id, scales = 'free_x', space = 'free_x') +
    scale_x_datetime(breaks = function(x) mean(x), date_labels = '%b %Y') +
    ggthm +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          strip.text  = element_blank()) +
    labs(x = NULL, y = 'correlation')

  # p_up / p_lo + plot_layout(heights = c(2, 1))
  p_up
}

pd_bm1  <- make_lomad_plot_data('BM1', loc_results)
pd_bs1  <- make_lomad_plot_data('BS1', loc_results)

plt_bm1 <- make_lomad_ggplot(pd_bm1, 'BM1 \u2013 lomad fit')
plt_bs1 <- make_lomad_ggplot(pd_bs1, 'BS1 \u2013 lomad fit')

wrap_elements(plt_bm1) / wrap_elements(plt_bs1)

ggsave(paste0(img_out, '/BM1-lomad-fit.pdf'), plt_bm1, width = 16, height = 5)
ggsave(paste0(img_out, '/BS1-lomad-fit.pdf'), plt_bs1, width = 16, height = 5)

print(plt_bm1)
print(plt_bs1)
