devtools::load_all()
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# ---- shared helpers --------------------------------------------------------

# Two-panel diagnostic for methods that return lambda
plot_trends <- function(sim, title = NULL) {
  n  <- length(sim$y1)
  df <- data.frame(t = seq_len(n), y1 = sim$y1, y2 = sim$y2,
                   x_mean = sim$x_mean, lambda = sim$lambda)

  p_series <- df |>
    pivot_longer(c(y1, y2, x_mean), names_to = "series", values_to = "value") |>
    ggplot(aes(x = t, y = value, color = series)) +
    geom_line(linewidth = 0.45) +
    scale_color_manual(values = c(y1 = "steelblue", y2 = "tomato",
                                  x_mean = "gray60")) +
    labs(title = title, y = NULL, color = NULL) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "none", axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          plot.title = element_text(size = 9))

  p_lambda <- ggplot(df, aes(x = t, y = lambda)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.5, ymax = 1.5,
             fill = "gray90", alpha = 0.6) +
    geom_hline(yintercept = c(0, 1), linetype = "dashed",
               color = "gray70", linewidth = 0.3) +
    geom_line(linewidth = 0.4) +
    labs(y = "λ", x = "t") +
    theme_minimal(base_size = 10)

  p_series / p_lambda + plot_layout(heights = c(2, 1))
}

# Single-panel plot for make_trends_dist (no lambda)
plot_dist <- function(sim, title = NULL) {
  n <- length(sim$x1)
  data.frame(t = seq_len(n), y1 = sim$x1, y2 = sim$x2) |>
    pivot_longer(c(y1, y2), names_to = "series", values_to = "value") |>
    ggplot(aes(x = t, y = value, color = series)) +
    geom_line(linewidth = 0.45) +
    scale_color_manual(values = c(y1 = "steelblue", y2 = "tomato")) +
    labs(title = title, y = NULL, x = "t", color = NULL) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "none", plot.title = element_text(size = 9))
}

ds <- c(0.5, 1, 2, 4)


# ---- make_trends_dist ------------------------------------------------------

plots_dist <- lapply(ds, \(d)
  plot_dist(make_trends_dist(n = 500, d = d, seed = 7),
            title = sprintf("d = %.1f", d))
)

wrap_plots(plots_dist, ncol = 2) +
  plot_annotation(title = "make_trends_dist(): varying d")


# ---- make_trends_rate ------------------------------------------------------
# rate = 0.01 (5 events), vary d

plots_rate <- lapply(ds, \(d)
  plot_trends(make_trends_rate(n = 500, d = d, rate = 0.01, seed = 7),
              title = sprintf("d = %.1f", d))
)

wrap_plots(plots_rate, ncol = 2) +
  plot_annotation(title = "make_trends_rate(): varying d  (rate = 0.01)")


# ---- make_trends_smooth ----------------------------------------------------
# bw = 50, coupling = 0.8, vary d

plots_smooth <- lapply(ds, \(d)
  plot_trends(make_trends_smooth(n = 500, d = d, bw = 50,
                                 coupling = 0.8, seed = 7),
              title = sprintf("d = %.1f", d))
)

wrap_plots(plots_smooth, ncol = 2) +
  plot_annotation(title = "make_trends_smooth(): varying d  (bw = 50, coupling = 0.8)")


# ---- make_trends_cross -----------------------------------------------------
# bw = 50, coupling = 0.8, vary d

plots_cross <- lapply(ds, \(d)
  plot_trends(make_trends_cross(n = 500, d = d, bw = 50,
                                coupling = 0.8, seed = 7),
              title = sprintf("d = %.1f", d))
)

wrap_plots(plots_cross, ncol = 2) +
  plot_annotation(title = "make_trends_cross(): varying d  (bw = 50, coupling = 0.8)")
