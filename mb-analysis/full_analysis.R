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

# Minimum length: m >= 2s, i.e. n >= 3s + h. The bare fitting minimum
# (n >= 2h + s) admits blocks whose every window overlaps every other one --
# roughly one effective look -- and per-block BY then hands the weakest
# evidence the laxest threshold. Requiring two window-lengths of valid test
# points drops those.
min_len <- 3L * s_win + h_win

for (loc in names(loc_results)) {
  blocks_presm <- loc_results[[loc]]$blocks_presm
  blocks <- lapply(blocks_presm, \(b) list(x1 = b$o2, x2 = b$ph))

  too_short <- vapply(blocks, \(b) length(b$x1) < min_len, logical(1))
  if (any(too_short))
    cat(sprintf('  %s: dropping %d blocks shorter than %d observations (m < 2s)\n',
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

make_lomad_ggplot <- function(pd, title_str) {
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
    labs(x = NULL, y = 'series', title = title_str)

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
  plt <- make_lomad_ggplot(pd, paste0(loc, ' - lomad fit (AR(1) noise, global BY)'))
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

xlim <- range(aligned$datetime)
shade_col <- rgb(0.7, 0.85, 1, 0.5)

p_series <- ggplot(aligned, aes(x = datetime, group = block_id)) +
  geom_rect(data = episodes, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = shade_col) +
  geom_line(aes(y = ma1), colour = 'blue', linewidth = 0.3) +
  geom_line(aes(y = ma2), colour = 'red',  linewidth = 0.3) +
  facet_grid(location ~ .) +
  scale_x_datetime(limits = xlim, date_labels = '%b %Y', expand = c(0.01, 0)) +
  ggthm +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  labs(x = NULL, y = 'smoothed series',
       title = 'Morro Bay: both stations on a common time axis',
       subtitle = 'blue = dissolved oxygen, red = pH; shading = flagged windows')

p_corr <- ggplot(aligned, aes(x = datetime, group = block_id)) +
  geom_rect(data = episodes, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = shade_col) +
  geom_hline(yintercept = 0, colour = 'grey80', linewidth = 0.3) +
  geom_line(aes(y = rho), colour = 'grey30', linetype = 'dashed') +
  geom_line(aes(y = R),   colour = 'grey20') +
  facet_grid(location ~ .) +
  scale_x_datetime(limits = xlim, date_labels = '%b %Y', expand = c(0.01, 0)) +
  ggthm +
  labs(x = NULL, y = 'correlation')

plt_aligned <- p_series / p_corr + plot_layout(heights = c(1, 1))
ggsave(paste0(img_out, '/stations-aligned.png'), plt_aligned,
       width = 12, height = 7, dpi = 200)
print(plt_aligned)
