# example-figure.R -- the vignette's worked example, as a paper figure
#
# Usage (from the repo root):
#   Rscript mb-analysis/example-figure.R
#
# Output
#   mb-analysis/_img/fig-mb-example.png
#
# Same data, h, s and alpha as vignette("lomad"), so the paper and the package
# report the same fit. Drawn in ggplot rather than with lomad_plot() so it
# matches the other Morro Bay figures: palette, line widths and theme come from
# utils.R. The two-panel layout and the shading convention are lomad_plot()'s.
#
# Note this reads the INSTALLED package's data, not anything in this repo. A
# stale install silently draws a different block, so the fit is checked against
# what the vignette reports before anything is written.

suppressPackageStartupMessages({
  library(tidyverse); library(patchwork); library(lomad)
})
source("mb-analysis/utils.R")

img_out <- "mb-analysis/_img"; fs::dir_create(img_out)
H <- 4L; S <- 60L; ALPHA <- 0.05

fit <- lomad_fit(morro_bay$o2, morro_bay$ph, h = H, s = S)
tst <- lomad_test(fit, alpha = ALPHA)

stopifnot(
  nrow(morro_bay) == 208L,
  as.Date(min(morro_bay$datetime)) == as.Date("2022-08-24"),
  sum(tst$rejected, na.rm = TRUE) == 52L
)

rejected <- replace(tst$rejected, is.na(tst$rejected), FALSE)
t_idx    <- morro_bay$datetime

# A rejection at t concerns the window ending at t, so the upper panel shades
# back to t - s + 1 while the lower panel shades t itself. The span comes from
# the package, so this figure and the vignette cannot disagree about which
# window a flag refers to.
runs <- function(flag) {
  r <- rle(flag); en <- cumsum(r$lengths); st <- en - r$lengths + 1L
  tibble(xmin = t_idx[st[r$values]], xmax = t_idx[en[r$values]])
}
shade_up <- runs(lomad:::.rejected_window_span(rejected, S))
shade_lo <- runs(rejected)

d <- tibble(datetime = t_idx, DO = fit$ma1, pH = fit$ma2, trend = fit$trend,
            R = fit$R, rho = fit$rho,
            crit = fit$rho + qnorm(tst$alpha_eff) * sqrt(fit$V / S))
ma <- d |> select(datetime, DO, pH) |>
  pivot_longer(-datetime, names_to = "var", values_to = "z")

p_up <- ggplot() +
  geom_rect(data = shade_up, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = SHADE, alpha = 0.45) +
  geom_line(data = d, aes(datetime, trend), colour = "grey45",
            linewidth = LW_MA * 0.8) +
  geom_line(data = ma, aes(datetime, z, colour = var), linewidth = LW_MA) +
  scale_colour_manual(values = VAR_PAL, name = NULL) +
  scale_x_datetime(date_breaks = "2 weeks", date_labels = "%d %b") +
  fig_theme() +
  theme(legend.position = c(0.995, 0.98), legend.justification = c(1, 1),
        legend.direction = "vertical",
        legend.background = element_rect(fill = alpha("white", 0.75), colour = NA),
        legend.key.width = unit(0.22, "in"),
        axis.text.x = element_blank(), axis.text.y = element_blank()) +
  labs(x = NULL, y = "Moving averages")

# The band is the region between rho and the critical value below which R is
# flagged. Both move with t, which is why the deepest dip in R need not be the
# flagged one. Correlation keeps its tick labels: unlike a standardized
# anomaly, the value is directly interpretable.
p_lo <- ggplot() +
  geom_rect(data = shade_lo, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = SHADE, alpha = 0.45) +
  geom_ribbon(data = filter(d, !is.na(crit), !is.na(rho)),
              aes(datetime, ymin = crit, ymax = rho), fill = "grey55", alpha = 0.25) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.2) +
  geom_line(data = d, aes(datetime, rho), colour = "grey30", linetype = "dashed",
            linewidth = LW_MA * 0.8) +
  geom_line(data = d, aes(datetime, R), colour = "grey15", linewidth = LW_MA) +
  scale_x_datetime(date_breaks = "2 weeks", date_labels = "%d %b") +
  fig_theme() +
  labs(x = NULL, y = "Correlation series")

ggsave(file.path(img_out, "fig-mb-example.png"),
       p_up / p_lo + plot_layout(heights = c(1.35, 1)),
       width = 5, height = 3, dpi = 450)

cat(sprintf("Wrote fig-mb-example.png -- %d rows, %s to %s, %d of %d flagged\n",
            nrow(morro_bay), as.Date(min(t_idx)), as.Date(max(t_idx)),
            sum(rejected), length(fit$valid_idx)))
