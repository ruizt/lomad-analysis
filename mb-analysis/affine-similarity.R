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
## Both examples are exactly S_WIN = 60 points, one 15-day analysis window, so
## the figure shows the scale the test operates on. They were selected by
## sweeping every 60-point window in the record rather than chosen by eye, and
## each isolates one half of the affine map: the level example has an amplitude
## ratio near 1, the amplitude example an offset near 0. A window with both at
## once reads as a level shift and buries the amplitude difference.
##
## The two come from different blocks, so each carries its own context panel.
##
## Usage (from the repo root):
##   Rscript mb-analysis/affine-similarity.R

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})
source("mb-analysis/utils.R")   # presmooth_tidal(), VAR_PAL, LW_*, figure theme

img_out <- "mb-analysis/_img"; fs::dir_create(img_out)

H_WIN   <- 4L      # 24-h moving average on the 6-hourly series
S_WIN   <- 60L     # 15-day correlation window
CTX_PAD <- 120L    # context extends this far each side of the window

FILL <- "#D9EAD3"

EXAMPLES <- list(
  list(loc = "BM1", blk = 33, i1 =  78L, tag = "level offset"),
  list(loc = "BS1", blk = 22, i1 = 217L, tag = "amplitude difference")
)

# ---- Preparation ------------------------------------------------------------

raw <- readr::read_csv("_mb-data/ph_o2_blocks.csv", show_col_types = FALSE)

prep <- function(loc, blk) {
  b  <- raw |> filter(location == loc, block_id == blk) |> arrange(datetime)
  ps <- presmooth_tidal(b, cols = c("o2", "ph")) |> arrange(datetime)
  kern <- rep(1 / H_WIN, H_WIN)
  ps |> mutate(ma_o2 = as.numeric(stats::filter(o2, kern, sides = 1)),
               ma_ph = as.numeric(stats::filter(ph, kern, sides = 1)))
}

# b_hat, a_hat from OLS of the pH moving average on the DO moving average over
# the window; kappa is the ratio of their standard deviations. The method never
# estimates b -- this is here to show what an affine map would absorb.
fit_window <- function(ps, i1, i2) {
  d  <- ps[i1:i2, ]
  ok <- is.finite(d$ma_o2) & is.finite(d$ma_ph)
  x  <- d$ma_o2[ok]; y <- d$ma_ph[ok]
  m  <- stats::lm(y ~ x)
  a  <- unname(coef(m)[1]); b <- unname(coef(m)[2])
  list(a = a, b = b, corr = cor(x, y), kappa = sd(y) / sd(x),
       offset = mean(y) - mean(x),
       rms0 = sqrt(mean((y - x)^2)),
       rms1 = sqrt(mean(((y - a) / b - x)^2)),
       t1 = d$datetime[1], t2 = d$datetime[nrow(d)])
}

long_ma <- function(d) {
  bind_rows(
    tibble(datetime = d$datetime, value = d$ma_o2, Series = "DO"),
    tibble(datetime = d$datetime, value = d$ma_ph, Series = "pH")
  ) |> filter(is.finite(value))
}

# ---- Panels -----------------------------------------------------------------

build <- function(ex, show_legend) {
  ps <- prep(ex$loc, ex$blk)
  i1 <- ex$i1; i2 <- i1 + S_WIN - 1L
  f  <- fit_window(ps, i1, i2)

  cat(sprintf(
    "\n%s -- %s block %d, idx %d-%d, %s to %s\n  corr %.3f | kappa %.2f | offset %+.2f | b_hat %.2f | a_hat %+.2f\n  RMS %.3f -> %.3f (%.0f%% reduction)\n",
    ex$tag, ex$loc, ex$blk, i1, i2, as.Date(f$t1), as.Date(f$t2),
    f$corr, f$kappa, f$offset, f$b, f$a, f$rms0, f$rms1,
    100 * (1 - f$rms1 / f$rms0)))

  ctx <- ps[max(1L, i1 - CTX_PAD):min(nrow(ps), i2 + CTX_PAD), ]
  p_ctx <- ggplot() +
    geom_rect(data = tibble(x1 = f$t1, x2 = f$t2),
              aes(xmin = x1, xmax = x2, ymin = -Inf, ymax = Inf),
              fill = FILL, alpha = 0.7) +
    geom_line(data = long_ma(ctx), aes(datetime, value, colour = Series),
              linewidth = LW_MA) +
    scale_colour_manual(values = VAR_PAL) +
    scale_x_datetime(date_labels = "%b %d") +
    labs(x = NULL, y = "Standardized units", colour = NULL,
         title = sprintf("%s, block %d -- shaded: the %s window below",
                         ex$loc, ex$blk, ex$tag)) +
    fig_theme() +
    theme(legend.position = if (show_legend) c(0.99, 0.02) else "none",
          legend.justification = c(1, 0), legend.direction = "horizontal",
          legend.background = element_blank(),
          legend.key.size = unit(0.35, "cm"))

  z <- ps[i1:i2, ]
  p_raw <- ggplot(long_ma(z), aes(datetime, value, colour = Series)) +
    geom_line(linewidth = LW_MA) +
    scale_colour_manual(values = VAR_PAL, guide = "none") +
    scale_x_datetime(date_labels = "%b %d") +
    labs(x = NULL, y = "Standardized units", title = "As observed",
         subtitle = bquote(kappa == .(sprintf("%.2f", f$kappa)) * "," ~
                           offset == .(sprintf("%+.2f", f$offset)))) +
    fig_theme() +
    theme(plot.subtitle = element_text(size = PT$annot, colour = "grey25"))

  z_adj <- z |> mutate(ma_ph = (ma_ph - f$a) / f$b)
  p_adj <- ggplot(long_ma(z_adj), aes(datetime, value, colour = Series)) +
    geom_line(linewidth = LW_MA) +
    scale_colour_manual(values = VAR_PAL, guide = "none") +
    scale_x_datetime(date_labels = "%b %d") +
    # Numbers are quoted so plotmath keeps the trailing zero: unquoted, 1.70 is
    # evaluated and renders as 1.7.
    labs(x = NULL, y = NULL, title = "Realigned",
         subtitle = bquote(hat(a) == .(sprintf("%+.2f", f$a)) * "," ~
                           hat(b) == .(sprintf("%.2f", f$b)) * "," ~
                           RMS ~ .(sprintf("-%.0f%%", 100*(1 - f$rms1/f$rms0))))) +
    fig_theme() +
    theme(plot.subtitle = element_text(size = PT$annot, colour = "grey25"))

  list(ctx = p_ctx, raw = p_raw, adj = p_adj)
}

e1 <- build(EXAMPLES[[1]], TRUE)
e2 <- build(EXAMPLES[[2]], FALSE)

fig <- e1$ctx / (e1$raw | e1$adj) / e2$ctx / (e2$raw | e2$adj) +
  plot_layout(heights = c(0.8, 1, 0.8, 1)) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")",
                  theme = theme(plot.tag = element_text(size = PT$ltitle)))

ggsave(file.path(img_out, "fig-affine-similarity.png"), fig,
       width = 6.5, height = 7.6, dpi = 450)
cat(sprintf("\nWrote %s\n", file.path(img_out, "fig-affine-similarity.png")))
