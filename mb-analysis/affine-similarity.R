## affine-similarity.R -- global standardisation does not remove local
## differences in level and amplitude
##
## Motivates the affine-invariant null against a pointwise-equality null: much
## of the apparent separation between the pH and DO trends is affine in origin,
## and the affine map that removes it is local, not global.
##
## Inputs
##   _mb-data/ph_o2_blocks.csv    from mb-analysis/process_blocks.R
##
## Outputs -> mb-analysis/_img/
##   fig-affine-similarity.png
##
## Region indices come from a Python prototype and are half-open there, so
## `[i1:i2]` is `(i1 + 1):i2` here. Region A matches the prototype only on that
## reading -- closed indexing shifts corr from 0.615 to 0.627.
##
## Usage (from the repo root):
##   Rscript mb-analysis/affine-similarity.R

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})
source("mb-analysis/utils.R")   # presmooth_tidal(), VAR_PAL, LW_*, figure theme

img_out <- "mb-analysis/_img"; fs::dir_create(img_out)

H_WIN <- 4L      # 24-h moving average on the 6-hourly series
S_WIN <- 60L     # 15-day correlation window
LOC   <- "BS1"
BLOCK <- 23      # longest block

A_IDX <- c(1235L, 1374L)   # location offset dominant
B_IDX <- c(1375L, 1514L)   # scale difference dominant; widened to match W1
CTX   <- c(1235L, 1612L)   # display range for the context panel

FILL_A <- "#D9EAD3"        # light green / light purple: distinct from the
FILL_B <- "#EAD9F0"        # blue and red the series themselves use

# ---- Data -------------------------------------------------------------------

blk <- readr::read_csv("_mb-data/ph_o2_blocks.csv", show_col_types = FALSE) |>
  filter(location == LOC, block_id == BLOCK) |>
  arrange(datetime)

# Notch the tidal band, then keep one observation per 6-hour bin.
ps <- presmooth_tidal(blk, cols = c("o2", "ph")) |> arrange(datetime)

kern <- rep(1 / H_WIN, H_WIN)
ps <- ps |> mutate(
  ma_o2 = as.numeric(stats::filter(o2, kern, sides = 1)),
  ma_ph = as.numeric(stats::filter(ph, kern, sides = 1)),
  idx   = row_number()
)

message(sprintf("%s block %d: %d raw hourly rows -> %d at 6-hourly, %s to %s",
                LOC, BLOCK, nrow(blk), nrow(ps),
                as.Date(min(ps$datetime)), as.Date(max(ps$datetime))))

# ---- Region fits ------------------------------------------------------------
# b_hat, a_hat from OLS of the pH moving average on the DO moving average over
# the region; kappa is the ratio of their standard deviations. The method never
# estimates b -- this is here to show what an affine map would absorb.

fit_region <- function(i1, i2, label) {
  d  <- ps[i1:i2, ]
  ok <- is.finite(d$ma_o2) & is.finite(d$ma_ph)
  x  <- d$ma_o2[ok]; y <- d$ma_ph[ok]
  m  <- stats::lm(y ~ x)
  a  <- unname(coef(m)[1]); b <- unname(coef(m)[2])
  list(label = label, i1 = i1, i2 = i2,
       t1 = d$datetime[1], t2 = d$datetime[nrow(d)],
       n = length(x), corr = cor(x, y), kappa = sd(y) / sd(x), a = a, b = b,
       rms0 = sqrt(mean((y - x)^2)),
       rms1 = sqrt(mean(((y - a) / b - x)^2)))
}

regA <- fit_region(A_IDX[1], A_IDX[2], "W1")
regB <- fit_region(B_IDX[1], B_IDX[2], "W2")

# Rolling correlation of the two moving averages, for the range checks quoted
# in the caption. Reported, not plotted.
R_t <- rep(NA_real_, nrow(ps))
for (t in S_WIN:nrow(ps)) {
  w <- (t - S_WIN + 1L):t
  if (all(is.finite(ps$ma_o2[w])) && all(is.finite(ps$ma_ph[w])))
    R_t[t] <- cor(ps$ma_o2[w], ps$ma_ph[w])
}

for (r in list(regA, regB)) {
  rr <- R_t[r$i1:r$i2]; rr <- rr[is.finite(rr)]
  cat(sprintf(
    "\nwindow %s  idx %d-%d  %s to %s  (n = %d)\n  corr %.3f | kappa %.3f | b_hat %.3f | a_hat %.3f\n  RMS %.3f -> %.3f (%.0f%% reduction)\n  rolling R_t median %.3f, range [%.2f, %.2f]\n",
    r$label, r$i1, r$i2, as.Date(r$t1), as.Date(r$t2), r$n,
    r$corr, r$kappa, r$b, r$a, r$rms0, r$rms1, 100 * (1 - r$rms1 / r$rms0),
    median(rr), min(rr), max(rr)))
}

# ---- Panels -----------------------------------------------------------------

long_ma <- function(d, o2 = "ma_o2", ph = "ma_ph") {
  bind_rows(
    tibble(datetime = d$datetime, value = d[[o2]], Series = "DO"),
    tibble(datetime = d$datetime, value = d[[ph]], Series = "pH")
  ) |> filter(is.finite(value))
}

ctx <- ps[CTX[1]:CTX[2], ]
band <- tibble(
  xmin = c(regA$t1, regB$t1), xmax = c(regA$t2, regB$t2),
  Region = c("W1", "W2")
)

p_a <- ggplot() +
  geom_rect(data = band, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf,
                             fill = Region), alpha = 0.55) +
  geom_line(data = long_ma(ctx), aes(datetime, value, colour = Series),
            linewidth = LW_MA) +
  geom_text(data = band, aes(x = xmin + (xmax - xmin) / 2, y = Inf, label = Region),
            vjust = 1.4, size = ANNOT, colour = "grey25") +
  scale_fill_manual(values = c(W1 = FILL_A, W2 = FILL_B), guide = "none") +
  scale_colour_manual(values = VAR_PAL) +
  scale_x_datetime(date_labels = "%b %d") +
  labs(x = NULL, y = "Standardized units", colour = NULL,
       title = "Globally standardized 24h moving averages") +
  fig_theme() +
  theme(legend.position = c(0.99, 0.02), legend.justification = c(1, 0),
        legend.direction = "horizontal", legend.background = element_blank(),
        legend.key.size = unit(0.35, "cm"))

# Each region gets an as-observed panel and a realigned one. b_hat is fitted
# over the whole region, not per window; the method never estimates b.
region_row <- function(r, show_y) {
  z <- ps[r$i1:r$i2, ]

  p_raw <- ggplot(long_ma(z), aes(datetime, value, colour = Series)) +
    geom_line(linewidth = LW_MA) +
    scale_colour_manual(values = VAR_PAL, guide = "none") +
    scale_x_datetime(date_labels = "%b %d") +
    labs(x = NULL, y = if (show_y) "Standardized units" else NULL,
         title = sprintf("Window %s (original)", r$label),
         subtitle = " ") +
    fig_theme() +
    theme(plot.subtitle = element_text(size = PT$annot))

  z_adj <- z |> mutate(ma_ph = (ma_ph - r$a) / r$b)
  p_adj <- ggplot(long_ma(z_adj), aes(datetime, value, colour = Series)) +
    geom_line(linewidth = LW_MA) +
    scale_colour_manual(values = VAR_PAL, guide = "none") +
    scale_x_datetime(date_labels = "%b %d") +
    # In the subtitle, not annotated inside: region B's trace reaches the top
    # left corner. Numbers are quoted so plotmath keeps the trailing zero --
    # unquoted, `1.70` is evaluated and renders as 1.7.
    labs(x = NULL, y = NULL,
         title = sprintf("Window %s (realigned)", r$label),
         subtitle = bquote(hat(a)[t] == .(sprintf("%.2f", r$a)) * "," ~
                           hat(b)[t] == .(sprintf("%.2f", r$b)))) +
    fig_theme() +
    theme(plot.subtitle = element_text(size = PT$annot, colour = "grey25"))

  list(raw = p_raw, adj = p_adj)
}

rowA <- region_row(regA, TRUE)
rowB <- region_row(regB, TRUE)

# Tags mark rows, not panels, and are set by hand as in
# simulations/simulation-results.R -- patchwork's tag_levels would letter all
# five panels.
fig <- (p_a + labs(tag = "A")) /
  ((rowA$raw + labs(tag = "B")) | rowA$adj) /
  ((rowB$raw + labs(tag = "C")) | rowB$adj) +
  plot_layout(heights = c(1, 1, 1)) &
  theme(plot.tag = element_text(size = PT$ltitle))

ggsave(file.path(img_out, "fig-affine-similarity.png"), fig,
       width = 6.5, height = 6.6, dpi = 450)
cat(sprintf("\nWrote %s\n", file.path(img_out, "fig-affine-similarity.png")))
