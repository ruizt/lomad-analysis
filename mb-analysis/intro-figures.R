# intro-figures.R -- site map and the coupling/decoupling illustration
#
# Usage (from the repo root):
#   Rscript mb-analysis/intro-figures.R
#
# Output
#   mb-analysis/_img/fig-mb-sites.png      map + both stations over one block
#   mb-analysis/_img/fig-mb-coupling.png   the same block, coupled -> decoupled
#
# Both cover exactly the block that ships with the package as morro_bay, so the
# introduction, this pair and fig-mb-example.png all show the same series.
#
# Shapefile: TIGER/Line area water for San Luis Obispo County. Downloaded on
# first run into _map/ (gitignored) since it is 441 KB of third-party data.

suppressPackageStartupMessages({
  library(tidyverse); library(sf); library(patchwork); library(lomad)
})
source("mb-analysis/utils.R")

img_out <- "mb-analysis/_img"; fs::dir_create(img_out)
map_dir <- "_map";             fs::dir_create(map_dir)

# DO blue, pH red, the same pure hues example-figure.R and analysis.R use for
# these two variables. Deliberately not the muted pair the seasonality raster
# uses for the two stations -- same family, but not the same colours, so a
# reader cannot carry "red = pH" across into a figure where red means Bay Mouth.
# Palette, line widths and theme come from utils.R, shared with example-figure.R
STN <- tibble(
  station = c("Bay Mouth (BM)", "Bay Head (BH)"),
  lon = c(-(120 + 51/60 + 32.04/3600), -(120 + 50/60 + 50.28/3600)),
  lat = c(  35 + 22/60 + 15.2394/3600,   35 + 20/60 + 1.6794/3600))

BOX <- list(x = c(-120.885, -120.788), y = c(35.305, 35.395))

compass_star <- function(cx, cy, r) {
  k <- 1 / cos(cy * pi / 180)          # degrees of longitude per degree of latitude
  th <- seq(90, 405, by = 45) * pi / 180
  rad <- rep(c(r, r * 0.36), length.out = length(th))
  data.frame(x = cx + k * rad * cos(th), y = cy + rad * sin(th))
}
STAR <- compass_star(-120.8035, 35.3195, 0.0075)

ggthm <- fig_theme()

# The sub-window the lower row of the composite shows: a week trimmed off each
# end of the block, which keeps the episode centred with context either side.
# Snapped to observation times, not to the nominal date. With no expansion the
# lower panel spans its first and last observation, so a nominal window would
# leave the box a few hours wider than the panel it marks.
TRIM      <- local({
  nom <- range(morro_bay$datetime) + c(7, -7) * 86400
  range(morro_bay$datetime[morro_bay$datetime >= nom[1] &
                           morro_bay$datetime <= nom[2]])
})
TRIM_FILL <- "grey92"

# No grid, with a border, for the composite. A grid across a shaded panel reads
# as two competing backgrounds.
# panel.grid.major has to be blanked by name: fig_theme() sets it explicitly,
# and ggplot does not let a parent element override an explicitly-set child.
fig_theme_framed <- function(fill = NA, grid = FALSE, border = TRUE) fig_theme() +
  theme(panel.grid.major = if (grid) element_line(linewidth = 0.1, colour = "grey65")
                           else element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = fill, colour = NA),
        panel.border = if (border) element_rect(fill = NA, colour = "grey40",
                                                linewidth = 0.3)
                       else element_blank())

# ---- shapefile --------------------------------------------------------------

shp <- file.path(map_dir, "tl_2023_06079_areawater.shp")
if (!file.exists(shp)) {
  zip <- file.path(map_dir, "slo_water.zip")
  download.file(paste0("https://www2.census.gov/geo/tiger/TIGER2023/",
                       "AREAWATER/tl_2023_06079_areawater.zip"), zip, quiet = TRUE)
  unzip(zip, exdir = map_dir)
}
water <- st_read(shp, quiet = TRUE) |> filter(MTFCC %in% c("H2051", "H2053"))
stn   <- st_as_sf(STN, coords = c("lon", "lat"), crs = 4269)

# ---- map, with a California locator inset ----------------------------------

ca <- st_as_sf(maps::map("state", "california", plot = FALSE, fill = TRUE))

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
  # labels below their points, so the top-right corner stays clear for the inset
  geom_sf_text(data = stn, aes(label = station), nudge_y = -0.0045,
               size = 2.8, colour = "#B22222", fontface = "bold") +
  geom_polygon(data = STAR, aes(x, y), fill = "grey25", colour = "grey25",
               linewidth = 0.2) +
  annotate("text", x = -120.8035, y = 35.3315, label = "N", size = 2.6,
           fontface = "bold", colour = "grey25") +
  coord_sf(xlim = BOX$x, ylim = BOX$y, expand = FALSE) +
  scale_x_continuous(breaks = c(-120.87, -120.83)) +
  scale_y_continuous(breaks = seq(35.32, 35.38, by = 0.02)) +
  theme_minimal(base_size = 9) +
  theme(panel.grid = element_line(linewidth = 0.1, colour = "grey85"),
        panel.background = element_rect(fill = "grey88", colour = NA),
        panel.border = element_rect(fill = NA, colour = "grey40", linewidth = 0.3),
        axis.title = element_blank(), axis.text = element_text(size = 6.5)) +
  inset_element(inset, left = 0.60, bottom = 0.50, right = 0.99, top = 0.99)

# ---- both stations over the example block ----------------------------------

W <- range(morro_bay$datetime)
pres <- read_csv("_mb-data/ph_o2_blocks.csv", show_col_types = FALSE) |>
  group_split(location, block_id) |>
  lapply(function(d) {
    if (max(d$datetime) < W[1] || min(d$datetime) > W[2]) return(NULL)
    suppressWarnings(presmooth_tidal(d, cols = c("o2", "ph"), step = "6 hours")) |>
      mutate(location = d$location[1], blk = d$block_id[1])
  }) |> bind_rows() |> filter(datetime >= W[1], datetime <= W[2])

ser_dat <- pres |>
  mutate(station = factor(c(BM1 = "Bay Mouth (BM)", BS1 = "Bay Head (BH)")[location],
                          levels = c("Bay Mouth (BM)", "Bay Head (BH)"))) |>
  pivot_longer(c(o2, ph), names_to = "var", values_to = "z") |>
  mutate(var = c(o2 = "DO", ph = "pH")[var])

make_series <- function(alpha = 1) ser_dat |>
  # grouped on the block too, so the line breaks at gaps instead of
  # interpolating across them
  ggplot(aes(datetime, z, colour = var, group = interaction(var, station, blk))) +
  geom_line(linewidth = LW_OBS, alpha = alpha) +
  facet_grid(station ~ .) +
  scale_colour_manual(values = VAR_PAL, name = NULL) +
  scale_x_datetime(date_breaks = "2 weeks", date_labels = "%d %b") +
  ggthm + theme(legend.position = "top",
                strip.text.y = element_text(angle = -90),
                strip.background = element_blank(),
                axis.text.y = element_blank()) +
  labs(x = NULL, y = NULL)

p_ser <- make_series()

ggsave(file.path(img_out, "fig-mb-sites.png"),
       p_map + p_ser + plot_layout(widths = c(1, 1.7)),
       width = FIG_W, height = 3.1, dpi = 450)

# ---- coupling / decoupling -------------------------------------------------

fit <- lomad_fit(morro_bay$o2, morro_bay$ph, h = 4, s = 60)
brk <- as.POSIXct("2022-09-25", tz = "UTC")   # the flagged run ends 26 Sep

d <- tibble(datetime = morro_bay$datetime, DO = morro_bay$o2, pH = morro_bay$ph,
            ma_DO = fit$ma1, ma_pH = fit$ma2)
raw <- d |> select(datetime, DO, pH) |> pivot_longer(-datetime, names_to = "var", values_to = "z")
ma  <- d |> select(datetime, ma_DO, ma_pH) |>
  pivot_longer(-datetime, names_to = "var", values_to = "z") |>
  mutate(var = sub("^ma_", "", var))

ytop <- max(d$DO, d$pH, na.rm = TRUE); ay <- ytop + 0.5
gap <- as.difftime(1, units = "days"); len <- as.difftime(5.5, units = "days")
ann <- tibble(x0 = c(brk - gap, brk + gap), x1 = c(brk - gap - len, brk + gap + len),
              label = c("decoupled", "coupled"))
ann$mid <- ann$x0 + (ann$x1 - ann$x0) / 2

# `win` filters rather than clips, and the scale takes no expansion, so the panel
# spans exactly the window the box above it marks. coord_cartesian alone left
# default padding on both sides, which is why the two did not line up.
make_coupling <- function(win = NULL, framed = FALSE) {
  keep <- function(x) if (is.null(win)) x else
    dplyr::filter(x, datetime >= win[1], datetime <= win[2])
  ggplot() +
  geom_vline(xintercept = brk, linetype = "dashed", colour = "grey35", linewidth = 0.45) +
  geom_line(data = keep(raw), aes(datetime, z, colour = var), linewidth = LW_OBS, alpha = 0.3) +
  geom_line(data = keep(ma),  aes(datetime, z, colour = var), linewidth = LW_MA) +
  geom_segment(data = ann, aes(x = x0, xend = x1, y = ay, yend = ay),
               arrow = arrow(length = unit(0.055, "in"), type = "closed"),
               colour = "grey35", linewidth = 0.35) +
  geom_text(data = ann, aes(mid, ay + 0.4, label = label), size = 3.2, colour = "grey20") +
  scale_colour_manual(values = VAR_PAL, name = NULL) +
  scale_x_datetime(date_breaks = "2 weeks", date_labels = "%d %b",
                   expand = if (is.null(win)) waiver() else expansion(0)) +
  coord_cartesian(xlim = win,
                  ylim = c(min(d$DO, d$pH, na.rm = TRUE), ay + 0.85)) +
  (if (framed) fig_theme_framed(TRIM_FILL, grid = TRUE) else ggthm) +
  theme(legend.position = c(0.995, 0.98), legend.justification = c(1, 1),
        legend.direction = "vertical",
        legend.background = element_rect(fill = alpha("white", 0.75), colour = NA),
        legend.key.width = unit(0.22, "in"),
        axis.text.y = element_blank()) +
  labs(x = NULL, y = "Moving averages (BM)")
}

ggsave(file.path(img_out, "fig-mb-coupling.png"), make_coupling(),
       width = FIG_W, height = 3.2, dpi = 450)

# ---- the two stacked ---------------------------------------------------------
# The site series go faint here, matching the background of the coupling panel
# below: the top row is context for the block, the bottom row is the point.
# Its legend is dropped rather than collected, since the coupling panel already
# carries one inside its own bounds.

# ---- the two stacked, in two variants --------------------------------------
# v1: both rows over the whole block, gridded, unframed. Reads as two views of
#     the same window.
# v2: the upper row marks a sub-window in the fill the lower panel is drawn on
#     and the lower panel is that window, so it reads as a detail figure. The
#     box sits only in the Bay Mouth facet, since Bay Mouth is what the panel
#     below shows, and is outlined so it reads as a marked region rather than a
#     change of background.
#
# Both drop the site legend rather than collecting it: the coupling panel
# carries one inside its own bounds, and collecting would take width from the
# map.

faint <- function() make_series(alpha = 0.35) + guides(colour = "none")

ggsave(file.path(img_out, "fig-mb-sites-coupling-v1.png"),
       (p_map + faint() + plot_layout(widths = c(1, 1.7))) /
         make_coupling() + plot_layout(heights = c(2, 1)),
       width = FIG_W, height = 6.6, dpi = 450)

box <- tibble(x1 = TRIM[1], x2 = TRIM[2],
              station = factor("Bay Mouth (BM)", levels = levels(ser_dat$station)))

p_ser_box <- faint() +
  geom_rect(data = box, inherit.aes = FALSE,
            aes(xmin = x1, xmax = x2, ymin = -Inf, ymax = Inf),
            fill = TRIM_FILL, colour = "grey40", linewidth = 0.3) +
  geom_line(linewidth = LW_OBS, alpha = 0.35) +   # redrawn over the box
  fig_theme_framed(border = FALSE) +
  theme(legend.position = "none", strip.text.y = element_text(angle = -90),
        strip.background = element_blank(), axis.text = element_blank()) +
  labs(x = NULL, y = NULL)

ggsave(file.path(img_out, "fig-mb-sites-coupling-v2.png"),
       (p_map + p_ser_box + plot_layout(widths = c(1, 1.7))) /
         make_coupling(win = TRIM, framed = TRUE) + plot_layout(heights = c(2, 1)),
       width = FIG_W, height = 6.6, dpi = 450)

cat("Wrote fig-mb-sites.png, fig-mb-coupling.png,",
    "fig-mb-sites-coupling-v1.png, fig-mb-sites-coupling-v2.png\n")
