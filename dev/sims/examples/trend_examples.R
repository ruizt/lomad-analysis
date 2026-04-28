devtools::load_all()
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

seed_coef   <- 2847   # seed for Fourier base (shared across structured panels)
seed_smooth <- 2222   # seed for stochastic repulsion G_t
seed_cross  <- 3333   # seed for stochastic crossing G_t

# ---- Base trends mu1, mu2 and midpoint mu_bar --------------------------
#
# Generate mu1, mu2 once at unit distance (d=1) in coefficient space.
# mu_bar = (mu1 + mu2) / 2 is the trend midpoint shown in structured panels.
# The structured methods mix mu1, mu2 via the coupling weight w_t; see
# eqn:w-mixing in the paper.

base   <- make_trends_dist(n = n, d = 1, seed = seed_coef)
mu1    <- base$x1
mu2    <- base$x2
mu_bar <- (mu1 + mu2) / 2

# ---- Coupling-weight generators ----------------------------------------

# Event-rate coupling weight (paper eqn:w-rate):
#   w_t = [1 - sum_j f(t - tau_j)]_+
# where f is a gamma density (shape=2) normalised to unit peak, and
# tau_1,...,tau_m are m = floor(r*T) evenly spaced event times.
make_w_rate <- function(n, rate) {
  m               <- round(rate * n)
  gap             <- n / m
  delta           <- 1 / m
  decouple_length <- round(gap * delta)
  decouple_str    <- m
  events <- round(seq(from = gap / 2, by = gap, length.out = m))
  t      <- seq_len(n)
  w      <- rep(1, n)
  for (i in events) {
    bump <- stats::dgamma(t - i, shape = 2, scale = decouple_length * 0.5)
    bump <- bump / max(bump)
    w    <- w - delta * decouple_str * bump
  }
  pmax(0, pmin(1, w))
}

# Stochastic repulsion coupling weight (paper eqn:w-smooth):
#   w_t = Phi(Phi^{-1}(c) - G_t)
# where G_t is a standardised kernel-smoothed Gaussian process with bandwidth b
# and c controls the fraction of time w_t > 1/2 (strong coupling).
make_w_smooth <- function(n, bw, coupling, seed) {
  set.seed(seed)
  G_raw    <- stats::rnorm(n)
  G_smooth <- stats::ksmooth(seq_len(n), G_raw, kernel = "normal",
                             bandwidth = bw, x.points = seq_len(n))$y
  G        <- (G_smooth - mean(G_smooth)) / stats::sd(G_smooth)
  stats::pnorm(stats::qnorm(coupling) - G)
}

# Stochastic crossing coupling weight (paper eqn:w-crossing):
#   w_t = 1 - sigma_w * G_t,  sigma_w = 0.5 / Phi^{-1}((1+c)/2)
# Allows trends to cross mu_bar; P(|w_t - 1| < 1/2) = c.
make_w_cross <- function(n, bw, coupling, seed) {
  set.seed(seed)
  G_raw    <- stats::rnorm(n)
  G_smooth <- stats::ksmooth(seq_len(n), G_raw, kernel = "normal",
                             bandwidth = bw, x.points = seq_len(n))$y
  G        <- (G_smooth - mean(G_smooth)) / stats::sd(G_smooth)
  sigma_w  <- 0.5 / stats::qnorm((1 + coupling) / 2)
  1 - sigma_w * G
}

# ---- Apply w_t and rescale to target d ---------------------------------
#
# Constructs nu1, nu2 via eqn:w-mixing:
#   nu_{it} = mu_bar_t + (1 - w_t)(mu_{it} - mu_bar_t)
# then rescales so that ||nu1 - nu2||_2 = d.
apply_w <- function(w, d) {
  nu1_unit <- mu_bar + (1 - w) * (mu1 - mu_bar)
  nu2_unit <- mu_bar + (1 - w) * (mu2 - mu_bar)
  sc       <- d / sqrt(sum((nu1_unit - nu2_unit)^2))
  list(
    nu1    = mu_bar + sc * (nu1_unit - mu_bar),
    nu2    = mu_bar + sc * (nu2_unit - mu_bar),
    mu_bar = mu_bar,
    w      = w
  )
}

# ---- Assemble trend lists ----------------------------------------------

tr_dist   <- make_trends_dist(n = n, d = d, seed = seed_coef)
tr_rate   <- apply_w(make_w_rate(n, rate),                        d)
tr_smooth <- apply_w(make_w_smooth(n, bw, coupling, seed_smooth), d)
tr_cross  <- apply_w(make_w_cross(n, bw, coupling, seed_cross),   d)

# ---- Common y-axis range -----------------------------------------------

y_lim <- range(c(
  tr_dist$x1,    tr_dist$x2,
  tr_rate$nu1,   tr_rate$nu2,   tr_rate$mu_bar,
  tr_smooth$nu1, tr_smooth$nu2, tr_smooth$mu_bar,
  tr_cross$nu1,  tr_cross$nu2,  tr_cross$mu_bar
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

# Unstructured panel: show mu1, mu2, and their midpoint mu_bar.
# (w_t = 0 so nu_it = mu_it; no rescaling needed beyond d.)
panel_top_dist <- function(tr, title) {
  mu_bar_loc <- (tr$x1 + tr$x2) / 2
  df <- data.frame(t = seq_len(n), mu1 = tr$x1, mu2 = tr$x2,
                   mu_bar = mu_bar_loc) |>
    pivot_longer(c(mu1, mu2, mu_bar), names_to = "series", values_to = "value") |>
    mutate(series = factor(series, levels = c("mu_bar", "mu1", "mu2")))
  ggplot(df, aes(t, value, color = series, linewidth = series,
                 linetype = series, alpha = series)) +
    geom_line() +
    scale_color_manual(values = c(mu_bar = "gray50", mu1 = "#0072B2", mu2 = "#D55E00")) +
    scale_linewidth_manual(values = c(mu_bar = 0.35, mu1 = 0.45, mu2 = 0.45)) +
    scale_linetype_manual(values = c(mu_bar = "dashed", mu1 = "solid", mu2 = "solid")) +
    scale_alpha_manual(values = c(mu_bar = 1, mu1 = 0.9, mu2 = 0.9)) +
    scale_y_continuous(limits = y_lim) +
    labs(title = title, y = expression(nu[it])) +
    theme_top
}

# Structured panel: show nu1, nu2, and mu_bar.
panel_top_struct <- function(tr, title) {
  df <- data.frame(t = seq_len(n), nu1 = tr$nu1, nu2 = tr$nu2,
                   mu_bar = tr$mu_bar) |>
    pivot_longer(c(nu1, nu2, mu_bar), names_to = "series", values_to = "value") |>
    mutate(series = factor(series, levels = c("mu_bar", "nu1", "nu2")))
  ggplot(df, aes(t, value, color = series, linewidth = series,
                 linetype = series, alpha = series)) +
    geom_line() +
    scale_color_manual(values = c(mu_bar = "gray50", nu1 = "#0072B2", nu2 = "#D55E00")) +
    scale_linewidth_manual(values = c(mu_bar = 0.35, nu1 = 0.45, nu2 = 0.45)) +
    scale_linetype_manual(values = c(mu_bar = "dashed", nu1 = "solid", nu2 = "solid")) +
    scale_alpha_manual(values = c(mu_bar = 1, nu1 = 0.9, nu2 = 0.9)) +
    scale_y_continuous(limits = y_lim) +
    labs(title = title, y = expression(nu[it])) +
    theme_top
}

# Coupling weight sub-panel: plot w_t.
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
  panel_wt(list(w = rep(0, n)), ref_lines = c(0, 1), ylim = c(-0.1, 1.1)) +
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

ggsave("dev/sims/examples/_img/fig_trend_construction.png",
       fig_1x4, width = 8, height = 2.5, units = "in", dpi = 400)
ggsave("../lomad-paper/img/fig_trend_construction.png",
       fig_1x4, width = 8, height = 2.5, units = "in", dpi = 400)
