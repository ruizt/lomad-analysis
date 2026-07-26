library(lomad)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# ---- Parameters --------------------------------------------------------

n        <- 500
d        <- 2       # L2 separation, held constant across all four panels
bw       <- 50      # bandwidth b for stochastic methods
coupling <- 0.8     # coupling fraction c (proportion of time w_t > 1/2)
rate     <- 0.01    # event rate r (events per unit time)

seed_coef <- 2847   # seed for Fourier base (shared across all panels)

# ---- Generate trend pairs via sim_trends() ------------------------------

tr_dist   <- sim_trends(n, d = d, method = "dist",   seed = seed_coef)
tr_rate   <- sim_trends(n, d = d, method = "rate",   seed = seed_coef, 
                        rate = rate)
tr_smooth <- sim_trends(n, d = d, method = "smooth", seed = seed_coef,
                        bw = bw, coupling = coupling)
tr_cross  <- sim_trends(n, d = d, method = "cross",  seed = seed_coef,
                        bw = bw, coupling = coupling)

# ---- Common y-axis range -----------------------------------------------

y_lim <- range(c(
  tr_dist$x1,   tr_dist$x2,
  tr_rate$x1,   tr_rate$x2,   tr_rate$x_mean,
  tr_smooth$x1, tr_smooth$x2, tr_smooth$x_mean,
  tr_cross$x1,  tr_cross$x2,  tr_cross$x_mean
))
y_pad <- diff(y_lim) * 0.05
y_lim <- y_lim + c(-y_pad, y_pad)

# ---- Shared themes -----------------------------------------------------

theme_top <- theme_minimal(base_size = 9) +
  theme(
    legend.position  = "none",
    axis.title.x     = element_blank(),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    plot.title       = element_text(size = 9, face = "bold"),
    panel.grid.minor = element_blank()
  )

theme_bot <- theme_minimal(base_size = 9) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank()
  )

# ---- Panel builders ----------------------------------------------------

# Unstructured panel: w = 0, so x1/x2 are the raw Fourier trends.
panel_top_dist <- function(tr, title) {
  x_mean_loc <- (tr$x1 + tr$x2) / 2
  df <- data.frame(t = seq_len(n), x1 = tr$x1, x2 = tr$x2,
                   x_mean = x_mean_loc) |>
    pivot_longer(c(x1, x2, x_mean), names_to = "series", values_to = "value") |>
    mutate(series = factor(series, levels = c("x_mean", "x1", "x2")))
  ggplot(df, aes(t, value, color = series, linewidth = series,
                 linetype = series, alpha = series)) +
    geom_line() +
    scale_color_manual(values = c(x_mean = "gray50", x1 = "#0072B2", x2 = "#D55E00")) +
    scale_linewidth_manual(values = c(x_mean = 0.35, x1 = 0.45, x2 = 0.45)) +
    scale_linetype_manual(values = c(x_mean = "dashed", x1 = "solid", x2 = "solid")) +
    scale_alpha_manual(values = c(x_mean = 1, x1 = 0.9, x2 = 0.9)) +
    scale_y_continuous(limits = y_lim) +
    labs(title = title, y = expression(nu[it])) +
    theme_top
}

# Structured panel: show x1, x2, and x_mean.
panel_top_struct <- function(tr, title) {
  df <- data.frame(t = seq_len(n), x1 = tr$x1, x2 = tr$x2,
                   x_mean = tr$x_mean) |>
    pivot_longer(c(x1, x2, x_mean), names_to = "series", values_to = "value") |>
    mutate(series = factor(series, levels = c("x_mean", "x1", "x2")))
  ggplot(df, aes(t, value, color = series, linewidth = series,
                 linetype = series, alpha = series)) +
    geom_line() +
    scale_color_manual(values = c(x_mean = "gray50", x1 = "#0072B2", x2 = "#D55E00")) +
    scale_linewidth_manual(values = c(x_mean = 0.35, x1 = 0.45, x2 = 0.45)) +
    scale_linetype_manual(values = c(x_mean = "dashed", x1 = "solid", x2 = "solid")) +
    scale_alpha_manual(values = c(x_mean = 1, x1 = 0.9, x2 = 0.9)) +
    scale_y_continuous(limits = y_lim) +
    labs(title = title, y = expression(nu[it])) +
    theme_top
}

# Coupling weight sub-panel.
panel_wt <- function(tr, ref_lines = c(0, 1), ylim = NULL) {
  if (is.null(ylim)) {
    rng  <- range(tr$w)
    pad  <- diff(rng) * 0.1
    ylim <- c(rng[1] - pad, rng[2] + pad)
  }
  data.frame(t = seq_len(n), w = tr$w) |>
    ggplot(aes(t, w)) +
    geom_hline(yintercept = ref_lines, linetype = "dashed",
               color = "gray65", linewidth = 0.3) +
    geom_line(linewidth = 0.4, color = "gray20") +
    scale_y_continuous(limits = ylim, breaks = ref_lines) +
    labs(y = expression(w[t]), x = "t") +
    theme_bot
}

# ---- Composite panels (series / w_t) -----------------------------------

comp_dist <- panel_top_dist(tr_dist, "(a) Unstructured") /
  panel_wt(tr_dist, ref_lines = c(0, 1), ylim = c(-0.1, 1.1)) +
  plot_layout(heights = c(3, 1))

comp_rate <- panel_top_struct(tr_rate, "(b) Event rate") /
  panel_wt(tr_rate) +
  plot_layout(heights = c(3, 1))

comp_smooth <- panel_top_struct(tr_smooth, "(c) Stochastic repulsion") /
  panel_wt(tr_smooth) +
  plot_layout(heights = c(3, 1))

comp_cross <- panel_top_struct(tr_cross, "(d) Stochastic crossing") /
  panel_wt(tr_cross, ref_lines = c(0, 1)) +
  plot_layout(heights = c(3, 1))

# ---- 2×2 figure --------------------------------------------------------

fig_2x2 <- (comp_dist | comp_rate) / (comp_smooth | comp_cross)
print(fig_2x2)

# ---- 1×4 figure --------------------------------------------------------

fig_1x4 <- comp_dist | comp_rate | comp_smooth | comp_cross
print(fig_1x4)

ggsave("_img/fig_trend_construction.png",
       fig_1x4, width = 8, height = 2.5, units = "in", dpi = 400)
