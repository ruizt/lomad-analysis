## analysis.R -- Morro Bay: fit, test, and every figure the paper uses
##
## Inputs
##   _mb-data/ph_o2_blocks.csv          from mb-analysis/process_blocks.R
##   _mb-data/scale_constants.csv       ditto; puts pH results in measured units
##   lomad::morro_bay                   the installed package's example block
##   _map/tl_2023_06079_areawater.shp   TIGER/Line, downloaded on first run
##
## Outputs -> mb-analysis/_img/
##   fig-mb-sites-coupling.png   site map, both stations, coupled -> decoupled
##   fig-affine-similarity.png   local level and amplitude differences
##   fig-mb-example.png          the vignette's worked example
##   fig-mb-seasonality.png      monthly detection rate over a phenology raster
##   fig-mb-reconstruction.png   pooled vs local DO->pH relationship, by state
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

# m >= 1.5s, below the power sweep's regime; checked separately against the
# global null.
MIN_LEN <- as.integer(2.5 * S_WIN) + H_WIN

STATION_NAME <- c(BM1 = "Bay Mouth (BM)", BS1 = "Bay Head (BH)")
STATION_ABBR <- c(BM1 = "BM", BS1 = "BH")

# panel border kept explicitly; the raster reads better framed
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
# Fit per block; one BY step-up over the pooled p-values, so FDR is
# controlled study-wide.
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

  # a container for the raw p-values; decisions are overwritten below
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
# Descriptive only: point estimates and sensitivity, not inference.
# =============================================================================

# `rejected` was NA-filled to FALSE upstream; keep only tested windows.
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

# Leave-one-block-out on the rate ratio.
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

# Timestamps where both stations were under test.
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
# Site map and both stations over the example block, above BM's moving averages.
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

# panel.grid.major must be set by name: fig_theme() sets it explicitly and a
# parent element cannot override an explicitly-set child.
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
  labs(x = NULL, y = NULL)

# The lower panel carries its own legend, so the site row drops one rather than
# collecting it, which would take width from the map.
ggsave(file.path(img_out, "fig-mb-sites-coupling.png"),
       (p_map + p_ser + plot_layout(widths = c(1, 1.7))) / p_cpl +
         plot_layout(heights = c(1, 1)),
       width = 5, height = 4.5, dpi = 450)

# =============================================================================
# fig-affine-similarity.png
# Two windows a season apart, 84 points each: 28 Jan - 18 Feb and September.
# =============================================================================

AFF_LOC <- "BS1"; AFF_BLK <- 23
AFF_W   <- list(W1 = c(1393L, 1476L), W2 = c(2292L, 2375L))   # 84 points each
AFF_CTX <- c(1284L, 2619L)   # 2025-01-01 to 2025-12-01, spanning both windows
AFF_FILL <- c(W1 = "#D9EAD3", W2 = "#EAD9F0")   # clear of the DO/pH colours

aff <- aligned |>
  filter(location == AFF_LOC, block_id == AFF_BLK) |>
  arrange(datetime)

# a_t, b_t from OLS of the pH moving average on the DO moving average.
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
# Each window's map applied to the other.
for (k in names(aff_f)) {
  o <- setdiff(names(aff_f), k)
  d <- aligned |> filter(location == AFF_LOC, block_id == AFF_BLK) |>
    arrange(datetime) |> slice(AFF_W[[k]][1]:AFF_W[[k]][2])
  ok <- is.finite(d$ma1) & is.finite(d$ma2); x <- d$ma1[ok]; y <- d$ma2[ok]
  g <- aff_f[[o]]
  cat(sprintf("  %s under %s's map: RMS %.3f -> %.3f (%+.0f%%)\n", k, o,
              sqrt(mean((y - x)^2)), sqrt(mean(((y - g$a) / g$b - x)^2)),
              100 * (1 - sqrt(mean(((y - g$a) / g$b - x)^2)) / sqrt(mean((y - x)^2)))))
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
# The vignette's worked example, drawn in ggplot to match the other figures.
# Reads the INSTALLED package's data.
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
# A rejection at t concerns the window ending at t; shade back to t - s + 1.
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
  theme(legend.position = c(1, 1), legend.justification = c(1, 1),
        legend.direction = "horizontal",
        legend.background = element_rect(fill = NA, colour = NA),
        legend.key.width = unit(0.22, "in"),
        axis.text.x = element_blank(), axis.text.y = element_blank()) +
  labs(x = NULL, y = expression(Y[t]))

# Band between rho and the critical value below which R is flagged; both move
# with t.
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
  labs(x = NULL, y = expression(R[t]))

ggsave(file.path(img_out, "fig-mb-example.png"),
       p_ex_up / p_ex_lo + plot_layout(heights = c(1.35, 1)),
       width = 3.5, height = 2, dpi = 450)

# =============================================================================
# sfig-mb-detections.png
# Every fitted block, both stations, series over correlation; independent x.
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
# Pooled monthly rate over a phenology raster: day of year across, year down.
# Grey is every window the test decided on.
# =============================================================================

pal <- c(BM = "#C44E52", BH = "#4C72B0")

# Tiles centred on integer doy with width 1.
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
  scale_fill_manual(values = pal) +
  scale_x_continuous(breaks = MONTH_MID, labels = month.abb,
                     limits = X_LIM, oob = scales::oob_keep,
                     expand = X_EXPAND) +
  scale_y_reverse(breaks = 2020:2025) +
  ggthm + theme(legend.position = "none",
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

# Pooled rates as a diverging strip, standing in for the legend. Bar lengths
# are comparable to each other, not to an absolute reference.
SEAS_W   <- 5      # export width, inches
PANEL_IN <- 4.2    # bar row once patchwork aligns it to the raster
GAP      <- 0.006  # hairline between the two bars
TITLE    <- "Pooled detection rates:"

rate_bm <- rate_tbl$rate[rate_tbl$location == "BM1"]
rate_bh <- rate_tbl$rate[rate_tbl$location == "BS1"]

grDevices::pdf(NULL)
title_in <- grid::convertWidth(grid::grobWidth(grid::textGrob(
  TITLE, gp = grid::gpar(fontsize = PT$annot))), "in", valueOnly = TRUE)
grDevices::dev.off()

t_x  <- 2 * title_in / PANEL_IN            # title width, panel spans 2
ceil <- (rate_bm + rate_bh) / (2 - t_x - 0.10)
f_bm <- rate_bm / ceil; f_bh <- rate_bh / ceil
cx   <- -1 + t_x + 0.10 + f_bm             # centre, so BM starts after the title

p_rate <- ggplot() +
  geom_text(aes(-1, 0, label = TITLE), hjust = 0, size = ANNOT, colour = "grey25") +
  geom_rect(data = tibble(xmin = c(cx - GAP - f_bm, cx + GAP),
                          xmax = c(cx - GAP, cx + GAP + f_bh),
                          fill = c(pal[["BM"]], pal[["BH"]])),
            aes(xmin = xmin, xmax = xmax, ymin = -0.5, ymax = 0.5, fill = fill)) +
  geom_text(data = tibble(x = c(cx - GAP - f_bm, cx + GAP + f_bh), h = c(-0.35, 1.35),
                          lab = sprintf("%.1f%%", 100 * c(rate_bm, rate_bh))),
            aes(x, 0, label = lab, hjust = h), colour = "white",
            size = ANNOT, fontface = "bold") +
  geom_text(data = tibble(x = c(cx - GAP - f_bm / 2, cx + GAP + f_bh / 2),
                          lab = c("BM", "BH"), col = c(pal[["BM"]], pal[["BH"]])),
            aes(x, -0.5, label = lab, colour = col), size = ANNOT, vjust = 1.5) +
  scale_fill_identity() + scale_colour_identity() +
  coord_cartesian(xlim = c(-1, 1), ylim = c(-1.6, 0.7), expand = FALSE) +
  labs(x = NULL, y = NULL) + theme_void(base_size = PT$title)

plt_seas <- p_seas / p_ras / p_rate + plot_layout(heights = c(1.9, 5, 0.85))
ggsave(file.path(img_out, "fig-mb-seasonality.png"), plt_seas,
       width = SEAS_W, height = 3.9, dpi = 450)

# =============================================================================
# fig-mb-reconstruction.png
# A pooled DO->pH fit hides the local breakdown that raises prediction error.
# =============================================================================

CV_K      <- 50L     # random splits
CV_TEST   <- 0.2     # held out per split
BW_ADJ    <- 1.5     # kde bandwidth; windows overlap, so the default undersmooths
# R_t is a correlation of moving averages, so its support runs back a further
# h - 1 observations than the s the correlation is taken over.
W_RAW     <- S_WIN + H_WIN - 1L
PH_SCALE  <- read_csv("_mb-data/scale_constants.csv", show_col_types = FALSE) |>
  filter(variable == "ph") |> pull(scale)
COUPLED   <- "grey45"
DECOUP    <- setNames(c("#C44E52", "#4C72B0"), STATION_NAME)

rec <- aligned |> filter(!is.na(o2), !is.na(ph)) |>
  arrange(location, block_id, datetime)
rec_key <- paste(rec$location, rec$block_id)

# rolling sum over the W_RAW rows ending at each position
roll_sum <- function(v, s) {
  cs <- c(0, cumsum(v)); n <- length(v); out <- rep(NA_real_, n)
  if (n >= s) out[s:n] <- cs[(s:n) + 1] - cs[(s:n) - s + 1]
  out
}

# within-window correlation of the raw pair
rec$w_r <- NA_real_
for (k in unique(rec_key)) {
  i <- which(rec_key == k); if (length(i) < W_RAW) next
  x <- rec$o2[i]; y <- rec$ph[i]; r <- rep(NA_real_, length(i))
  for (t in W_RAW:length(i)) {
    w <- (t - W_RAW + 1L):t
    if (sd(x[w]) > 0 && sd(y[w]) > 0) r[t] <- cor(x[w], y[w])
  }
  rec$w_r[i] <- r
}

# Per-station line fit on training observations only; each window's error uses
# only the held-out observations inside it, so nothing informs its own window.
set.seed(4471)
cv_mat <- matrix(NA_real_, nrow(rec), CV_K)
for (j in seq_len(CV_K)) {
  test <- runif(nrow(rec)) < CV_TEST
  fit  <- rec[!test, ] |> group_by(location) |>
    summarise(b0 = coef(lm(ph ~ o2))[1], b1 = coef(lm(ph ~ o2))[2], .groups = "drop")
  m  <- match(rec$location, fit$location)
  e2 <- ifelse(test, (rec$ph - (fit$b0[m] + fit$b1[m] * rec$o2))^2, 0)
  for (k in unique(rec_key)) {
    i <- which(rec_key == k); if (length(i) < W_RAW) next
    num <- roll_sum(e2[i], W_RAW); den <- roll_sum(as.numeric(test[i]), W_RAW)
    cv_mat[i, j] <- ifelse(!is.na(den) & den > 0, sqrt(num / den), NA_real_)
  }
}
rec$cv_rmspe <- PH_SCALE * rowMeans(cv_mat, na.rm = TRUE)

rec <- rec |> filter(is.finite(cv_rmspe), !is.na(w_r)) |>
  mutate(station = factor(STATION_NAME[location], levels = STATION_NAME),
         state   = factor(ifelse(rejected, "decoupled", "coupled"),
                          levels = c("coupled", "decoupled")),
         col     = ifelse(state == "coupled", COUPLED, DECOUP[as.character(station)]),
         lcol    = ifelse(state == "coupled", "grey65", DECOUP[as.character(station)]))

rec_pooled <- rec |> group_by(station) |> summarise(r = cor(o2, ph), .groups = "drop")

rec_lab <- rec |> group_by(station, state, col) |>
  summarise(r = cor(o2, ph), n = n(), .groups = "drop") |>
  mutate(lab = sprintf("%s:  r = %.2f  (n = %s)", state, r, format(n, big.mark = ",")),
         vj  = ifelse(state == "coupled", -8.2, -6.7))

p_rec_scatter <- ggplot(rec, aes(o2, ph, colour = col)) +
  geom_point(data = ~filter(.x, state == "coupled"),   size = 0.55, stroke = 0, alpha = 0.18) +
  geom_point(data = ~filter(.x, state == "decoupled"), size = 0.70, stroke = 0, alpha = 0.55) +
  geom_smooth(aes(colour = lcol), method = "lm", formula = y ~ x,
              se = FALSE, linewidth = 0.6) +
  geom_text(data = rec_lab, aes(x = Inf, y = -Inf, label = lab, vjust = vj),
            hjust = 1.04, size = ANNOT, show.legend = FALSE) +
  facet_wrap(~station) +
  scale_colour_identity() +
  labs(x = NULL, y = "pH ~ DO") +
  fig_theme() +
  theme(axis.text = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.text = element_text(size = PT$title))

p_rec_dens <- ggplot(rec, aes(w_r, fill = col, colour = col)) +
  geom_density(alpha = 0.35, linewidth = 0.35, adjust = BW_ADJ) +
  geom_vline(data = rec_pooled, aes(xintercept = r), inherit.aes = FALSE,
             linetype = "dashed", colour = "grey25", linewidth = 0.4) +
  geom_text(data = rec_pooled, aes(x = r, y = Inf, label = "pooled"),
            inherit.aes = FALSE, hjust = 1.12, vjust = 1.6,
            size = ANNOT, colour = "grey25") +
  facet_wrap(~station) +
  scale_colour_identity() + scale_fill_identity() +
  # baseline drawn in data space so it stops at the limits, not the panel edge
  annotate("segment", x = -1, xend = 1, y = 0, yend = 0,
           linewidth = 0.25, colour = "grey40") +
  scale_x_continuous(limits = c(-1, 1), breaks = c(-1, 0, 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "local correlation") +
  fig_theme() +
  theme(strip.text = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(linewidth = 0.25, colour = "grey40"),
        axis.ticks.length.x = unit(1.5, "pt"))

rec_ratio <- rec |> group_by(station, state) |>
  summarise(m = mean(cv_rmspe), .groups = "drop") |>
  pivot_wider(names_from = state, values_from = m) |>
  mutate(lab = sprintf("%.2f\u00d7 increase", decoupled / coupled),
         col = DECOUP[as.character(station)]) |>
  # rests on the top edge of the decoupled box
  left_join(rec |> filter(state == "decoupled") |> group_by(station) |>
              summarise(y = as.numeric(quantile(cv_rmspe, 0.75)), .groups = "drop"),
            by = "station")

p_rec_err <- ggplot(rec, aes(state, cv_rmspe, fill = col)) +
  geom_boxplot(outlier.size = 0.2, outlier.alpha = 0.2, linewidth = 0.3,
               width = 0.5, colour = "grey25") +
  geom_text(data = rec_ratio, aes(x = 1, y = y, label = lab, colour = col),
            inherit.aes = FALSE, hjust = -0.06, vjust = -0.6, size = ANNOT) +
  facet_wrap(~station) +
  scale_fill_identity() + scale_colour_identity() +
  scale_x_discrete(limits = c("decoupled", "coupled")) +
  scale_y_continuous(breaks = c(0.05, 0.10, 0.15)) +
  labs(x = NULL, y = "local RMSPE (pH)") +
  fig_theme() +
  theme(strip.text = element_blank(), panel.grid.major.x = element_blank())

plt_rec <- (p_rec_scatter / p_rec_dens / p_rec_err) +
  plot_layout(heights = c(1.1, 0.8, 0.8))
ggsave(file.path(img_out, "fig-mb-reconstruction.png"), plt_rec,
       width = 5.5, height = 6, units = "in", dpi = 450)

cat("\nLocal RMSPE (pH units), decoupled / coupled:\n")
for (k in seq_len(nrow(rec_ratio)))
  cat(sprintf("  %-16s %.4f / %.4f = %.2f\n", rec_ratio$station[k],
              rec_ratio$decoupled[k], rec_ratio$coupled[k],
              rec_ratio$decoupled[k] / rec_ratio$coupled[k]))
cat(sprintf("  %-16s %.4f / %.4f = %.2f\n", "pooled",
            mean(rec$cv_rmspe[rec$state == "decoupled"]),
            mean(rec$cv_rmspe[rec$state == "coupled"]),
            mean(rec$cv_rmspe[rec$state == "decoupled"]) /
              mean(rec$cv_rmspe[rec$state == "coupled"])))

# Only when someone is watching: printing under Rscript opens a device and
# leaves an Rplots.pdf behind.
if (interactive()) {
  print(plt_fits)
  print(plt_seas)
  print(plt_rec)
}

cat("\nFigures written to ", img_out, "\n", sep = "")
