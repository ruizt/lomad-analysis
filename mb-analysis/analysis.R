## analysis.R -- Morro Bay: fit, test, and every figure the paper uses
##
## Inputs
##   _mb-data/ph_o2_blocks.csv          from mb-analysis/process_blocks.R
##   lomad::morro_bay                   the installed package's example block
##   _map/tl_2023_06079_areawater.shp   TIGER/Line, downloaded on first run
##
## Outputs -> mb-analysis/_img/
##   fig-mb-sites-coupling.png   site map, both stations, coupled -> decoupled
##   fig-affine-similarity.png   local level and amplitude differences
##   fig-mb-example.png          the vignette's worked example
##   fig-mb-seasonality.png      monthly detection rate over a phenology raster
##   sfig-mb-detections.png      every fitted block, both stations
##
## Outputs -> _mb-data/
##   lomad_windows.rds           one row per analysed window, for ad hoc work
##
## Rates and summaries are printed, not written.
##
## Usage (from the repo root), either of:
##   Rscript mb-analysis/analysis.R
##   source("mb-analysis/analysis.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(sf)
  library(lomad)
})
source("mb-analysis/utils.R")   # presmooth_tidal(), VAR_PAL, LW_*, figure-theme

img_out <- "mb-analysis/_img"; fs::dir_create(img_out)
map_dir <- "_map";             fs::dir_create(map_dir)

H_WIN <- 4L
S_WIN <- 60L
ALPHA <- 0.05

# m >= 1.5s. Below the regime the power sweep covered (n = 4s throughout), so
# it was checked against the global null in
# simulations/calibration/blocklength-calibration.R: FDR 0.013, against 0.05.
MIN_LEN <- as.integer(2.5 * S_WIN) + H_WIN

STATION_NAME <- c(BM1 = "Bay Mouth (BM)", BS1 = "Bay Head (BH)")
STATION_ABBR <- c(BM1 = "BM", BS1 = "BH")

# theme_minimal like every other figure, with the panel border kept explicitly:
# the raster reads better framed.
ggthm <- theme_minimal(base_size = PT$title) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "black", linewidth = 0.1),
        panel.border = element_rect(fill = NA, colour = "grey40", linewidth = 0.3)) +
  fig_sizes()


# =============================================================================
# Presmooth
# Notch out tidal periodicity (~25 h cycle at hourly) and downsample to 6 h.
# =============================================================================

block_data <- read_csv("_mb-data/ph_o2_blocks.csv", show_col_types = FALSE)

presm_all <- lapply(split(block_data, block_data$location), function(loc_dat) {
  out <- loc_dat |>
    group_split(block_id) |>
    lapply(presmooth_tidal, cols = c("o2", "ph"), step = "6 hours")
  setNames(out, vapply(out, \(b) as.character(b$block_id[[1]]), character(1)))
})


# =============================================================================
# Fit per block, test globally
# One Benjamini-Yekutieli step-up correction over the pooled p-values, so FDR
# is controlled study-wide rather than per block.
# =============================================================================

loc_results <- lapply(names(presm_all), function(loc) {
  blocks_presm <- presm_all[[loc]]
  blocks <- lapply(blocks_presm, \(b) list(x1 = b$o2, x2 = b$ph))

  too_short <- vapply(blocks, \(b) length(b$x1) < MIN_LEN, logical(1))
  if (any(too_short))
    cat(sprintf("  %s: dropping %d blocks shorter than %d observations\n",
                loc, sum(too_short), MIN_LEN))
  blocks       <- blocks[!too_short]
  blocks_presm <- blocks_presm[!too_short]

  # lomad_test() here is only a container for the raw p-values; its decisions
  # are overwritten by the global stage below.
  block_fits <- lapply(names(blocks), \(nm) {
    b <- blocks[[nm]]
    tryCatch({
      fit <- lomad_fit(b$x1, b$x2, h = H_WIN, s = S_WIN)
      list(fit = fit, tst = lomad_test(fit, alpha = ALPHA))
    }, error = function(e) { message("  Block ", nm, " failed: ", e$message); NULL })
  })
  names(block_fits) <- names(blocks)
  block_fits <- Filter(Negate(is.null), block_fits)

  cat(sprintf("=== %s === %d / %d blocks fit\n",
              loc, length(block_fits), length(blocks)))

  list(loc = loc, blocks = blocks, blocks_presm = blocks_presm,
       block_fits = block_fits)
})
names(loc_results) <- names(presm_all)

pooled <- bind_rows(lapply(names(loc_results), function(loc) {
  bf <- loc_results[[loc]]$block_fits
  if (length(bf) == 0) return(NULL)
  bind_rows(lapply(names(bf), function(nm) {
    tibble(loc = loc, blk = nm, idx = bf[[nm]]$fit$valid_idx,
           p = bf[[nm]]$tst$p_values[bf[[nm]]$fit$valid_idx])
  }))
}))

# For a step-up procedure the rejection set is {p_raw <= p_star}, where p_star
# is the realized threshold.
M      <- nrow(pooled)
p_adj  <- p.adjust(pooled$p, method = "BY")
k      <- sum(p_adj <= ALPHA)
c_M    <- sum(1 / seq_len(M))
p_star <- ALPHA * max(k, 1L) / (M * c_M)

cat(sprintf("\nGlobal BY: M = %d pooled tests over %d blocks; %d rejected (p_star = %.3g)\n",
            M, n_distinct(paste(pooled$loc, pooled$blk)), k, p_star))

# Map the global decisions back so every downstream consumer sees them.
pooled$p_adj <- p_adj
for (loc in names(loc_results)) {
  bf <- loc_results[[loc]]$block_fits
  for (nm in names(bf)) {
    fit  <- bf[[nm]]$fit
    rows <- pooled$loc == loc & pooled$blk == nm
    stopifnot(sum(rows) == length(fit$valid_idx))
    tst <- bf[[nm]]$tst
    tst$p_adj[fit$valid_idx]    <- pooled$p_adj[rows]
    tst$rejected[fit$valid_idx] <- pooled$p[rows] <= p_star
    tst$alpha_eff               <- p_star
    loc_results[[loc]]$block_fits[[nm]]$tst <- tst
  }
}

rejection_summary <- bind_rows(lapply(names(loc_results), function(loc) {
  bf <- loc_results[[loc]]$block_fits
  bind_rows(lapply(names(bf), function(nm) {
    tst <- bf[[nm]]$tst
    data.frame(location = loc, block_id = nm,
               n_valid = sum(!is.na(tst$rejected)),
               n_rejected = sum(tst$rejected, na.rm = TRUE),
               frac_rejected = sum(tst$rejected, na.rm = TRUE) / sum(!is.na(tst$rejected)),
               detected = sum(tst$rejected, na.rm = TRUE) > 0)
  }))
}))

cat("\n========== Rejection summary ==========\n")
print(rejection_summary, row.names = FALSE)

print(rejection_summary |>
  group_by(location) |>
  summarise(n_blocks = n(), blocks_detected = sum(detected),
            total_valid = sum(n_valid), total_rejected = sum(n_rejected),
            frac_rejected = total_rejected / total_valid, .groups = "drop"))

# One row per analysed window; everything below runs from this.
aligned <- bind_rows(lapply(names(loc_results), function(loc) {
  bf <- loc_results[[loc]]$block_fits
  if (length(bf) == 0) return(NULL)
  bind_rows(lapply(names(bf), function(nm) {
    presm <- loc_results[[loc]]$blocks_presm[[nm]]
    fit   <- bf[[nm]]$fit
    tibble(location = loc,
           block_id = presm$block_id[[1]],
           datetime = presm$datetime,
           o2       = loc_results[[loc]]$blocks[[nm]]$x1,
           ph       = loc_results[[loc]]$blocks[[nm]]$x2,
           ma1      = fit$ma1,
           ma2      = fit$ma2,
           R        = fit$R,
           rho      = fit$rho,
           rejected = replace_na(bf[[nm]]$tst$rejected, FALSE))
  }))
}))

saveRDS(aligned, "_mb-data/lomad_windows.rds")
message("Wrote _mb-data/lomad_windows.rds (", nrow(aligned), " windows)")


# =============================================================================
# Station comparison
# Descriptive only. A seasonal hypothesis would have to be specified after
# seeing these data, and a rate interval would be inference on inference. What
# follows the point estimates is sensitivity, not uncertainty.
# =============================================================================

# `rejected` was NA-filled to FALSE upstream; a window is under test only where
# the fit produced both a correlation and a benchmark for it.
wv <- aligned |>
  mutate(valid = !is.na(R) & !is.na(rho), blk = paste(location, block_id)) |>
  filter(valid)

rate_tbl <- wv |>
  group_by(location) |>
  summarise(n_blocks = n_distinct(blk), n_valid = n(),
            n_rej = sum(rejected), rate = n_rej / n_valid, .groups = "drop")

cat("\n========== Rate by station ==========\n")
print(rate_tbl |> mutate(rate_pct = round(100 * rate, 1)))
cat(sprintf("rate ratio BM1/BS1 = %.2f\n",
            rate_tbl$rate[rate_tbl$location == "BM1"] /
            rate_tbl$rate[rate_tbl$location == "BS1"]))

# Does the ratio depend on any one block, or on the stations having been up at
# different times?
ratio_of <- function(d) {
  r <- tapply(d$rejected, d$location, mean)
  unname(r["BM1"] / r["BS1"])
}

blk_tbl <- wv |> group_by(location, blk) |>
  summarise(n = n(), rej = sum(rejected), frac = rej / n, .groups = "drop")

loo <- vapply(unique(wv$blk), \(b) ratio_of(filter(wv, blk != b)), numeric(1))
cat(sprintf("\nleave-one-block-out ratio: %.2f to %.2f over %d blocks\n",
            min(loo), max(loo), length(loo)))
for (loc in c("BM1", "BS1")) {
  d  <- filter(blk_tbl, location == loc)
  lo <- vapply(seq_len(nrow(d)), \(i) sum(d$rej[-i]) / sum(d$n[-i]), numeric(1))
  cat(sprintf("  %s rate spans %.1f%% to %.1f%%\n",
              loc, 100 * min(lo), 100 * max(lo)))
}

# Restricting to timestamps where both stations were under test removes
# differing exposure as an explanation.
conc <- wv |>
  select(location, datetime, rejected) |>
  pivot_wider(id_cols = datetime, names_from = location,
              values_from = rejected, values_fn = any) |>
  drop_na(BM1, BS1)
cat(sprintf("concurrent support (%d timestamps, %.0f days): BM1 %.1f%%, BS1 %.1f%%, ratio %.2f\n",
            nrow(conc), nrow(conc) * 6 / 24,
            100 * mean(conc$BM1), 100 * mean(conc$BS1),
            mean(conc$BM1) / mean(conc$BS1)))

seas <- wv |> mutate(month = month(datetime)) |>
  group_by(month) |>
  summarise(n = n(), rej = sum(rejected), frac = rej / n, .groups = "drop")

seas_loc <- wv |> mutate(month = month(datetime)) |>
  group_by(location, month) |>
  summarise(n = n(), rej = sum(rejected), .groups = "drop")

cat("\n========== Detections by month, both stations pooled ==========\n")
print(seas |> mutate(pct = round(100 * frac, 1)) |>
      transmute(month = month.abb[month], n, rej, pct), n = 12)

cat("\nmonths with no detection at either station: ")
cat(paste(month.abb[setdiff(1:12, unique(seas_loc$month[seas_loc$rej > 0]))],
          collapse = ", "), "\n")
cat("share of flags in Mar-Jun: ")
cat(paste(wv |> filter(rejected) |> group_by(location) |>
          summarise(p = sprintf("%s %.0f%%", location[1],
                                100 * mean(month(datetime) %in% 3:6)),
                    .groups = "drop") |> pull(p), collapse = ", "), "\n")


# =============================================================================
# fig-mb-sites-coupling.png
# Site map and both stations over the example block, above Bay Mouth's moving
# averages for the same window.
# =============================================================================

STN <- tibble(
  station = c("Bay Mouth (BM)", "Bay Head (BH)"),
  lon = c(-(120 + 51/60 + 32.04/3600), -(120 + 50/60 + 50.28/3600)),
  lat = c(  35 + 22/60 + 15.2394/3600,   35 + 20/60 + 1.6794/3600))

BOX <- list(x = c(-120.885, -120.788), y = c(35.305, 35.395))
PANEL_FILL <- "grey92"

compass_star <- function(cx, cy, r) {
  k <- 1 / cos(cy * pi / 180)          # degrees of longitude per degree of latitude
  th <- seq(90, 405, by = 45) * pi / 180
  rad <- rep(c(r, r * 0.36), length.out = length(th))
  data.frame(x = cx + k * rad * cos(th), y = cy + rad * sin(th))
}
STAR <- compass_star(-120.8035, 35.3195, 0.0075)

# `grid` takes a colour or FALSE. panel.grid.major has to be set by name:
# fig_theme() sets it explicitly, and ggplot does not let a parent element
# override an explicitly-set child.
fig_panel <- function(fill = NA, grid = FALSE, border = TRUE) fig_theme() +
  theme(panel.grid.major = if (isFALSE(grid)) element_blank()
                           else element_line(linewidth = 0.1, colour = grid),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = fill, colour = NA),
        panel.border = if (border) element_rect(fill = NA, colour = "grey40",
                                                linewidth = 0.3)
                       else element_blank())

shp <- file.path(map_dir, "tl_2023_06079_areawater.shp")
if (!file.exists(shp)) {
  zip <- file.path(map_dir, "slo_water.zip")
  download.file(paste0("https://www2.census.gov/geo/tiger/TIGER2023/",
                       "AREAWATER/tl_2023_06079_areawater.zip"), zip, quiet = TRUE)
  unzip(zip, exdir = map_dir)
}
water <- st_read(shp, quiet = TRUE) |> filter(MTFCC %in% c("H2051", "H2053"))
stn   <- st_as_sf(STN, coords = c("lon", "lat"), crs = 4269)
ca    <- st_as_sf(maps::map("state", "california", plot = FALSE, fill = TRUE))

inset <- ggplot() +
  geom_sf(data = ca, fill = "grey80", colour = "grey40", linewidth = 0.2) +
  annotate("rect", xmin = BOX$x[1] - 0.5, xmax = BOX$x[2] + 0.5,
           ymin = BOX$y[1] - 0.5, ymax = BOX$y[2] + 0.5,
           fill = NA, colour = "#B22222", linewidth = 0.6) +
  coord_sf(expand = FALSE) +
  theme_void() +
  theme(panel.background = element_rect(fill = "white", colour = "grey40",
                                        linewidth = 0.3),
        plot.margin = margin(1, 1, 1, 1))

p_map <- ggplot() +
  geom_sf(data = water, fill = "white", colour = "grey60", linewidth = 0.2) +
  geom_sf(data = stn, size = 2.4, colour = "#B22222") +
  # labels below their points, leaving the top-right corner for the inset
  geom_sf_text(data = stn, aes(label = station), nudge_y = -0.0045,
               size = ANNOT, colour = "#B22222", fontface = "bold") +
  geom_polygon(data = STAR, aes(x, y), fill = "grey25", colour = "grey25",
               linewidth = 0.2) +
  annotate("text", x = -120.8035, y = 35.3315, label = "N", size = ANNOT,
           fontface = "bold", colour = "grey25") +
  coord_sf(xlim = BOX$x, ylim = BOX$y, expand = FALSE) +
  scale_x_continuous(breaks = c(-120.87, -120.83)) +
  scale_y_continuous(breaks = seq(35.32, 35.38, by = 0.02)) +
  fig_theme() +
  theme(panel.grid = element_line(linewidth = 0.1, colour = "grey85"),
        panel.background = element_rect(fill = "grey88", colour = NA),
        panel.border = element_rect(fill = NA, colour = "grey40", linewidth = 0.3),
        axis.title = element_blank()) +
  inset_element(inset, left = 0.60, bottom = 0.50, right = 0.99, top = 0.99)

# Every block overlapping the example window, not just the analysed ones.
W <- range(morro_bay$datetime)
ser_dat <- imap(presm_all, \(blks, loc) {
    imap(blks, \(b, nm) mutate(b, location = loc, blk = nm)) |> bind_rows()
  }) |>
  bind_rows() |>
  filter(datetime >= W[1], datetime <= W[2]) |>
  mutate(station = factor(STATION_NAME[location], levels = STATION_NAME)) |>
  pivot_longer(c(o2, ph), names_to = "var", values_to = "z") |>
  mutate(var = c(o2 = "DO", ph = "pH")[var])

p_ser <- ser_dat |>
  # grouped on the block too, so lines break at gaps instead of interpolating
  ggplot(aes(datetime, z, colour = var, group = interaction(var, station, blk))) +
  geom_line(linewidth = LW_OBS, alpha = 0.35) +
  facet_grid(station ~ .) +
  scale_colour_manual(values = VAR_PAL, name = NULL) +
  scale_x_datetime(date_breaks = "2 weeks", date_labels = "%d %b") +
  guides(colour = "none") +
  # values are standardized anomalies, so a faint grid but no tick labels
  fig_panel(grid = "grey88", border = FALSE) +
  theme(strip.text.y = element_text(angle = -90),
        strip.background = element_blank(),
        axis.text = element_blank()) +
  labs(x = NULL, y = NULL)

cpl_fit <- lomad_fit(morro_bay$o2, morro_bay$ph, h = H_WIN, s = S_WIN)
brk     <- as.POSIXct("2022-09-25", tz = "UTC")   # the flagged run ends 26 Sep

cpl_d <- tibble(datetime = morro_bay$datetime, DO = morro_bay$o2,
                pH = morro_bay$ph, ma_DO = cpl_fit$ma1, ma_pH = cpl_fit$ma2)
cpl_raw <- cpl_d |> select(datetime, DO, pH) |>
  pivot_longer(-datetime, names_to = "var", values_to = "z")
cpl_ma  <- cpl_d |> select(datetime, ma_DO, ma_pH) |>
  pivot_longer(-datetime, names_to = "var", values_to = "z") |>
  mutate(var = sub("^ma_", "", var))

ay  <- max(cpl_d$DO, cpl_d$pH, na.rm = TRUE) + 0.5   # annotation row
gap <- as.difftime(1, units = "days"); len <- as.difftime(5.5, units = "days")
ann <- tibble(x0 = c(brk - gap, brk + gap), x1 = c(brk - gap - len, brk + gap + len),
              label = c("decoupled", "coupled"))
ann$mid <- ann$x0 + (ann$x1 - ann$x0) / 2

p_cpl <- ggplot() +
  geom_vline(xintercept = brk, linetype = "dashed", colour = "grey35", linewidth = 0.45) +
  geom_line(data = cpl_raw, aes(datetime, z, colour = var), linewidth = LW_OBS, alpha = 0.3) +
  geom_line(data = cpl_ma,  aes(datetime, z, colour = var), linewidth = LW_MA) +
  geom_segment(data = ann, aes(x = x0, xend = x1, y = ay, yend = ay),
               arrow = arrow(length = unit(0.055, "in"), type = "closed"),
               colour = "grey35", linewidth = 0.35) +
  geom_text(data = ann, aes(mid, ay + 0.4, label = label), size = ANNOT, colour = "grey20") +
  scale_colour_manual(values = VAR_PAL, name = NULL) +
  scale_x_datetime(date_breaks = "2 weeks", date_labels = "%d %b") +
  # a duplicated axis with no breaks or labels is just a right-hand title,
  # matching the facet strips in the row above
  scale_y_continuous(sec.axis = dup_axis(name = "Bay Mouth (BM)",
                                         breaks = NULL, labels = NULL)) +
  coord_cartesian(ylim = c(min(cpl_d$DO, cpl_d$pH, na.rm = TRUE), ay + 0.85)) +
  fig_panel(PANEL_FILL, grid = "grey65") +
  theme(legend.position = c(0.995, 0.98), legend.justification = c(1, 1),
        legend.direction = "vertical",
        legend.background = element_rect(fill = alpha("white", 0.75), colour = NA),
        legend.key.width = unit(0.22, "in"),
        axis.text.y = element_blank(),
        # vjust nudges a rotated title horizontally; the gutter's width is set
        # by the map's latitude labels, which patchwork matches across rows
        axis.title.y       = element_text(margin = margin(r = 1), vjust = 0),
        axis.title.y.right = element_text(margin = margin(l = 1))) +
  labs(x = NULL, y = "Moving average")

# The lower panel carries its own legend, so the site row drops one rather than
# collecting it, which would take width from the map.
ggsave(file.path(img_out, "fig-mb-sites-coupling.png"),
       (p_map + p_ser + plot_layout(widths = c(1, 1.7))) / p_cpl +
         plot_layout(heights = c(1, 1)),
       width = 6, height = 5, dpi = 450)


# =============================================================================
# fig-affine-similarity.png
# Global standardization does not remove local differences in level and
# amplitude, and the affine map that would remove them is local, not global.
# Motivates the affine-invariant null against a pointwise-equality null.
#
# Two adjacent windows of equal width in the longest block. W2 is 140 points to
# match W1 rather than the 98 that maximise the contrast: at 98 its correlation
# is 0.834 and the map removes 58% of the RMS separation, against 0.774 and 40%
# here. The stronger 140-point windows nearby all overlap W1.
# =============================================================================

AFF_LOC <- "BS1"; AFF_BLK <- 23
AFF_W   <- list(W1 = c(1235L, 1374L), W2 = c(1375L, 1514L))
AFF_CTX <- c(1235L, 1612L)
AFF_FILL <- c(W1 = "#D9EAD3", W2 = "#EAD9F0")   # clear of the DO/pH colours

aff <- aligned |>
  filter(location == AFF_LOC, block_id == AFF_BLK) |>
  arrange(datetime)

# a_t, b_t from OLS of the pH moving average on the DO moving average over the
# window. The method never estimates b; this shows what an affine map absorbs.
aff_fit <- function(ix) {
  d  <- aff[ix[1]:ix[2], ]
  ok <- is.finite(d$ma1) & is.finite(d$ma2)
  x  <- d$ma1[ok]; y <- d$ma2[ok]
  m  <- stats::lm(y ~ x)
  a  <- unname(coef(m)[1]); b <- unname(coef(m)[2])
  list(a = a, b = b, corr = cor(x, y), kappa = sd(y) / sd(x),
       t1 = d$datetime[1], t2 = d$datetime[nrow(d)],
       rms0 = sqrt(mean((y - x)^2)),
       rms1 = sqrt(mean(((y - a) / b - x)^2)))
}
aff_long <- function(d, ph = d$ma2) {
  bind_rows(tibble(datetime = d$datetime, value = d$ma1, Series = "DO"),
            tibble(datetime = d$datetime, value = ph,    Series = "pH")) |>
    filter(is.finite(value))
}

aff_f <- lapply(AFF_W, aff_fit)

cat("\n========== Affine similarity windows ==========\n")
for (k in names(aff_f)) {
  f <- aff_f[[k]]
  cat(sprintf("%s  %s to %s  corr %.3f | kappa %.2f | a_t %+.2f | b_t %.2f | RMS %.3f -> %.3f (%.0f%%)\n",
              k, as.Date(f$t1), as.Date(f$t2), f$corr, f$kappa, f$a, f$b,
              f$rms0, f$rms1, 100 * (1 - f$rms1 / f$rms0)))
}

aff_band <- tibble(Region = names(aff_f),
                   xmin = as.POSIXct(sapply(aff_f, \(f) f$t1), origin = "1970-01-01"),
                   xmax = as.POSIXct(sapply(aff_f, \(f) f$t2), origin = "1970-01-01"))

p_aff_ctx <- ggplot() +
  geom_rect(data = aff_band,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = Region),
            alpha = 0.55) +
  geom_line(data = aff_long(aff[AFF_CTX[1]:AFF_CTX[2], ]),
            aes(datetime, value, colour = Series), linewidth = LW_MA) +
  geom_text(data = aff_band,
            aes(x = xmin + (xmax - xmin) / 2, y = Inf, label = Region),
            vjust = 1.4, size = ANNOT, colour = "grey25") +
  scale_fill_manual(values = AFF_FILL, guide = "none") +
  scale_colour_manual(values = VAR_PAL) +
  scale_x_datetime(date_labels = "%b %d") +
  labs(x = NULL, y = "Standardized units", colour = NULL,
       title = "Globally standardized 24h moving averages") +
  fig_theme() +
  theme(legend.position = c(0.99, 0.02), legend.justification = c(1, 0),
        legend.direction = "horizontal", legend.background = element_blank(),
        legend.key.size = unit(0.35, "cm"))

aff_row <- function(k) {
  f <- aff_f[[k]]; z <- aff[AFF_W[[k]][1]:AFF_W[[k]][2], ]
  raw <- ggplot(aff_long(z), aes(datetime, value, colour = Series)) +
    geom_line(linewidth = LW_MA) +
    scale_colour_manual(values = VAR_PAL, guide = "none") +
    scale_x_datetime(date_labels = "%b %d", breaks = scales::breaks_pretty(3)) +
    labs(x = NULL, y = "Standardized units",
         title = sprintf("Window %s (original)", k), subtitle = " ") +
    fig_theme() + theme(plot.subtitle = element_text(size = PT$annot))
  # Coefficients are quoted inside the expression so plotmath keeps a trailing
  # zero; unquoted, 1.70 is evaluated and renders as 1.7.
  adj <- ggplot(aff_long(z, (z$ma2 - f$a) / f$b), aes(datetime, value, colour = Series)) +
    geom_line(linewidth = LW_MA) +
    scale_colour_manual(values = VAR_PAL, guide = "none") +
    scale_x_datetime(date_labels = "%b %d", breaks = scales::breaks_pretty(3)) +
    labs(x = NULL, y = NULL, title = sprintf("Window %s (realigned)", k),
         subtitle = bquote(hat(a)[t] == .(sprintf("%.2f", f$a)) * "," ~
                           hat(b)[t] == .(sprintf("%.2f", f$b)))) +
    fig_theme() +
    theme(plot.subtitle = element_text(size = PT$annot, colour = "grey25"))
  list(raw = raw, adj = adj)
}
aff_r1 <- aff_row("W1"); aff_r2 <- aff_row("W2")

# Tags mark rows, not panels: patchwork's tag_levels would letter all five.
ggsave(file.path(img_out, "fig-affine-similarity.png"),
       (p_aff_ctx + labs(tag = "A")) /
         ((aff_r1$raw + labs(tag = "B")) | aff_r1$adj) /
         ((aff_r2$raw + labs(tag = "C")) | aff_r2$adj) +
         plot_layout(heights = c(1, 1, 1)) &
         theme(plot.tag = element_text(size = PT$ltitle)),
       width = 5, height = 5, dpi = 450)


# =============================================================================
# fig-mb-example.png
# The vignette's worked example. Drawn in ggplot rather than with lomad_plot()
# so it matches the other figures; the two-panel layout and the shading
# convention are lomad_plot()'s.
#
# Reads the INSTALLED package's data, so the fit is checked against what the
# vignette reports: a stale install would silently draw a different block.
# =============================================================================

ex_fit <- lomad_fit(morro_bay$o2, morro_bay$ph, h = H_WIN, s = S_WIN)
ex_tst <- lomad_test(ex_fit, alpha = ALPHA)

stopifnot(
  nrow(morro_bay) == 208L,
  as.Date(min(morro_bay$datetime)) == as.Date("2022-08-24"),
  sum(ex_tst$rejected, na.rm = TRUE) == 15L
)

ex_rej <- replace(ex_tst$rejected, is.na(ex_tst$rejected), FALSE)
t_idx  <- morro_bay$datetime

runs <- function(flag) {
  r <- rle(flag); en <- cumsum(r$lengths); st <- en - r$lengths + 1L
  tibble(xmin = t_idx[st[r$values]], xmax = t_idx[en[r$values]])
}
# A rejection at t concerns the window ending at t, so the upper panel shades
# back to t - s + 1. The span comes from the package, so this figure and the
# vignette cannot disagree about which window a flag refers to.
shade_up <- runs(lomad:::.rejected_window_span(ex_rej, S_WIN))
shade_lo <- runs(ex_rej)

ex_d <- tibble(datetime = t_idx, DO = ex_fit$ma1, pH = ex_fit$ma2,
               R = ex_fit$R, rho = ex_fit$rho,
               crit = ex_fit$rho + qnorm(ex_tst$alpha_eff) * sqrt(ex_fit$V / S_WIN))
ex_ma <- ex_d |> select(datetime, DO, pH) |>
  pivot_longer(-datetime, names_to = "var", values_to = "z")

p_ex_up <- ggplot() +
  geom_rect(data = shade_up, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = SHADE, alpha = 0.45) +
  geom_line(data = ex_ma, aes(datetime, z, colour = var), linewidth = LW_MA) +
  scale_colour_manual(values = VAR_PAL, name = NULL) +
  scale_x_datetime(date_breaks = "2 weeks", date_labels = "%d %b") +
  fig_theme() +
  theme(legend.position = c(0.995, 0.98), legend.justification = c(1, 1),
        legend.direction = "vertical",
        legend.background = element_rect(fill = alpha("white", 0.75), colour = NA),
        legend.key.width = unit(0.22, "in"),
        axis.text.x = element_blank(), axis.text.y = element_blank()) +
  labs(x = NULL, y = "Moving average")

# The band runs between rho and the critical value below which R is flagged.
# Both move with t, which is why the deepest dip in R need not be the flagged
# one. Correlation keeps its tick labels: the value is directly interpretable.
p_ex_lo <- ggplot() +
  geom_rect(data = shade_lo, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = SHADE, alpha = 0.45) +
  geom_ribbon(data = filter(ex_d, !is.na(crit), !is.na(rho)),
              aes(datetime, ymin = crit, ymax = rho), fill = "grey55", alpha = 0.25) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.2) +
  geom_line(data = ex_d, aes(datetime, rho), colour = "grey30", linetype = "dashed",
            linewidth = LW_MA * 0.8) +
  geom_line(data = ex_d, aes(datetime, R), colour = "grey15", linewidth = LW_MA) +
  scale_x_datetime(date_breaks = "2 weeks", date_labels = "%d %b") +
  fig_theme() +
  labs(x = NULL, y = "Correlation series")

ggsave(file.path(img_out, "fig-mb-example.png"),
       p_ex_up / p_ex_lo + plot_layout(heights = c(1.35, 1)),
       width = 5, height = 3, dpi = 450)


# =============================================================================
# sfig-mb-detections.png
# Every fitted block, both stations, series over correlation. Independent x
# scales: the stations have different blocks, so this is a per-block view.
# =============================================================================

make_lomad_plot_data <- function(loc_name, results) {
  loc <- results[[loc_name]]
  bf  <- loc$block_fits

  block_dfs <- shade_list <- list()

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
    dates <- presm$datetime
    bid   <- presm$block_id[[1]]

    rejected <- replace_na(bf[[nm]]$tst$rejected, FALSE)

    block_dfs[[nm]] <- tibble(
      datetime = dates, block_id = bid,
      o2 = loc$blocks[[nm]]$x1, ph = loc$blocks[[nm]]$x2,
      ma1 = fit$ma1, ma2 = fit$ma2,
      R = fit$R, rho = fit$rho, rejected = rejected
    )

    # A rejection at t concerns the window {t - s + 1, ..., t}, so shade back
    # to its start. Same helper fig-mb-example uses.
    rej_shifted <- lomad:::.rejected_window_span(rejected, fit$inputs$s)

    shade_list[[paste0(nm, "_up")]] <- add_shade(rej_shifted, dates, bid, "upper")
    shade_list[[paste0(nm, "_lo")]] <- add_shade(rejected, dates, bid, "lower")
  }

  list(main = bind_rows(block_dfs), shade = bind_rows(shade_list))
}

# Returns the two panels rather than a composed plot, so both stations stack.
make_lomad_ggplot <- function(pd, loc) {
  shade_col <- rgb(0.7, 0.85, 1, 0.4)

  p_up <- ggplot(pd$main, aes(x = datetime)) +
    geom_rect(data = filter(pd$shade, panel == "upper"), inherit.aes = FALSE,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              fill = shade_col) +
    geom_line(aes(y = o2), color = rgb(0, 0, 1, 0.2)) +
    geom_line(aes(y = ph), color = rgb(1, 0, 0, 0.2)) +
    geom_line(aes(y = ma1), color = "blue", linewidth = 0.3) +
    geom_line(aes(y = ma2), color = "red",  linewidth = 0.3) +
    facet_grid(~block_id, scales = "free_x", space = "free_x") +
    scale_x_datetime(breaks = function(x) mean(x), date_labels = "%b %Y") +
    ggthm +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          strip.text   = element_blank(),
          plot.title = element_text(face = "plain")) +
    labs(x = NULL, y = "moving average", title = STATION_NAME[[loc]])

  p_lo <- ggplot(pd$main, aes(x = datetime)) +
    geom_rect(data = filter(pd$shade, panel == "lower"), inherit.aes = FALSE,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              fill = shade_col) +
    geom_hline(yintercept = 0, color = "grey80", linewidth = 0.3) +
    geom_line(aes(y = R),   color = "grey40") +
    geom_line(aes(y = rho), color = "grey30", linetype = "dashed") +
    facet_grid(~block_id, scales = "free_x", space = "free_x") +
    scale_x_datetime(breaks = function(x) mean(x), date_labels = "%b %Y") +
    ggthm +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          strip.text  = element_blank()) +
    labs(x = NULL, y = "correlation")

  list(up = p_up, lo = p_lo)
}

panels <- unlist(lapply(names(loc_results), function(loc) {
  if (length(loc_results[[loc]]$block_fits) == 0) return(NULL)
  make_lomad_ggplot(make_lomad_plot_data(loc, loc_results), loc)
}), recursive = FALSE)

plt_fits <- wrap_plots(panels, ncol = 1) +
  plot_layout(heights = rep(c(2, 1), length(panels) / 2))
ggsave(file.path(img_out, "sfig-mb-detections.png"), plt_fits,
       width = 10, height = 5, dpi = 450)


# =============================================================================
# fig-mb-seasonality.png
# Pooled monthly rate over a phenology raster: day of year across, year down,
# upper lane BM, lower lane BH. Grey is every window the test reached a
# decision on, which excludes the leading s + h - 2 points of each block.
# =============================================================================

pal <- c(BM = "#C44E52", BH = "#4C72B0")

# Tiles are centred on integer doy with width 1, so a month's visual centre
# sits half a day right of its arithmetic midpoint.
MONTH_END   <- cumsum(c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31))
MONTH_START <- c(0, head(MONTH_END, -1))
MONTH_MID   <- (MONTH_START + MONTH_END + 1) / 2
X_EXPAND    <- expansion(mult = c(0.01, 0.045))   # right gutter holds the lane labels
X_LIM       <- c(0, 367)    # shared, not just breaks: expand() is relative to
                            # each panel's own range, which differs
LANE_X      <- 367          # just past the last tile's edge

ras <- wv |>
  mutate(year    = year(datetime), doy = yday(datetime),
         lane    = year + ifelse(location == "BM1", -0.19, 0.19),
         station = STATION_ABBR[location])

lanes <- distinct(ras, lane, station)

p_ras <- ggplot(ras, aes(doy, lane)) +
  geom_tile(fill = "grey86", height = 0.34, width = 1) +
  geom_tile(data = filter(ras, rejected), aes(fill = station),
            height = 0.34, width = 1) +
  geom_text(data = lanes, aes(x = LANE_X, y = lane, label = station),
            inherit.aes = FALSE, hjust = 0, size = ANNOT, colour = "grey35") +
  # Labelled by a function, not a vector: values= is matched by name but
  # labels= is matched by position, so a literal vector would follow the
  # scale's own (alphabetical) break order and swap the two stations.
  scale_fill_manual(values = pal, labels = \(x) paste(x, "detection")) +
  scale_x_continuous(breaks = MONTH_MID, labels = month.abb,
                     limits = X_LIM, oob = scales::oob_keep,
                     expand = X_EXPAND) +
  scale_y_reverse(breaks = 2020:2025) +
  ggthm + theme(legend.position = "bottom", legend.title = element_blank(),
                panel.grid.major.y = element_blank(),
                axis.ticks.length = unit(0, "in"),
                panel.border = element_blank(),
                axis.text.y = element_text(angle = 90, hjust = 0.5)) +
  labs(x = NULL, y = "Year")

# Marginal strip above the raster, on the same x axis. Its height is set by the
# rotated y title, which is bounded by panel height.
p_seas <- seas |>
  mutate(mid = MONTH_MID[month]) |>
  ggplot(aes(mid, 100 * frac)) +
  geom_line(linewidth = 0.4, colour = "grey25") +
  geom_point(size = 1.3, colour = "grey15") +
  scale_x_continuous(breaks = MONTH_MID, labels = month.abb,
                     limits = X_LIM, oob = scales::oob_keep,
                     expand = X_EXPAND) +
  scale_y_continuous(breaks = c(0, 10, 20)) +
  ggthm + theme(axis.text.x  = element_blank(),
                axis.ticks = element_blank(),
                panel.border = element_blank(),
                panel.grid.major.y = element_line(linewidth = 0.1, color = "darkgrey")) +
  labs(x = NULL, y = "Detections (%)")

plt_seas <- p_seas / p_ras + plot_layout(heights = c(1.9, 5))
ggsave(file.path(img_out, "fig-mb-seasonality.png"), plt_seas,
       width = 5, height = 3.5, dpi = 450)

# Only when someone is watching: printing under Rscript opens a device and
# leaves an Rplots.pdf behind.
if (interactive()) {
  print(plt_fits)
  print(plt_seas)
}

cat("\nFigures written to ", img_out, "\n", sep = "")
