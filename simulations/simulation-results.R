## simulation-results.R -- figures and tables for the simulation studies
##
## Reads compiled summaries only, never raw per-job files, so it runs in
## seconds. The expensive stages are simulations/<study>/collect-results.R and
## simulations/power/collect-results.R.
##
## Inputs
##   simulations/power/results/simulations-power-curves.rds
##   simulations/power/results/simulations-power-roc.rds
##   simulations/power/results/simulations-power-auc.rds
##   simulations/validation/results/simulations-validation-results.rds
##   (fig-trends.png simulates its own data)
##
## Outputs -> simulations/_img/
##   fig-trends.png             methods of simulating trend separation
##   fig-power.png              local power and classification accuracy
##   fig-validation.png         finite-sample accuracy of the CLT
##
## Outputs -> simulations/_tbl/
##   tbl-localization-auc.csv         AUC and concordance per design cell
##
## Usage (from the repo root):
##   Rscript simulations/simulation-results.R

suppressPackageStartupMessages({
  library(lomad)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(ggplot2)
  library(ggh4x)
  library(patchwork)
})
source(here::here("figure-theme.R"))   # PT, ANNOT, fig_sizes()

IMG_DIR <- "simulations/_img"
TBL_DIR <- "simulations/_tbl"
alpha   <- 0.05
dir.create(IMG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TBL_DIR, showWarnings = FALSE, recursive = TRUE)

# Keyed by internal code, not display label: a named vector spliced into
# another named vector via c() gets its names silently mangled.
STRUCT_FULL <- c(fr = "Fixed Rate", rs = "Random Separation",
                  rm = "Random Mixing")
STRUCT_ABBR <- c(fr = "FR", rs = "RS", rm = "RM")
STRUCT_HEX  <- c(fr = "#009E73", rs = "#0072B2", rm = "#D55E00")
STRUCT_PAL  <- setNames(STRUCT_HEX, STRUCT_ABBR[names(STRUCT_HEX)])

# =============================================================================
# fig-trends.png -- ways of distributing separation across a series pair
# =============================================================================

local({

n        <- 500
d        <- 2       # base amplitude of the distinct component
bw       <- 50      # bandwidth b
coupling <- 0.8     # coupling fraction c
rate     <- 0.01    # event rate r

seed_coef <- 2847   # Fourier base, shared across panels

# The affine layer is illustrative, not the study's calibration: the cap is
# 5% per window against the sweep's 1.5%, and map_seed is picked so a_t rises
# and b_t falls monotonically. At the sweep's settings the drift is too small
# to see at this scale.
cap      <- 0.05
a_mult   <- 0.9     # a_t amplitude, as a multiple of a trend's sd
map_seed <- 40

# One d for every structure, unlike the power study, which scales d per
# structure to put them on a common delta_t. Here the base pair is drawn behind
# each panel, and a per-structure d would leave the mixed trends short of it
# even where w_t = 0 -- FR reaches only 45% of the backdrop at d = 0.45 d. The
# columns are then not on a common delta_t, which an illustration of *where*
# separation falls does not need.
struct_args <- list(
  dist = list(),
  fr   = list(rate = rate, bump = "gaussian"),
  rs   = list(bw = bw, coupling = coupling),
  rm   = list(bw = bw, coupling = coupling))

amp <- sd(sim_trends(n, d = d, method = "dist", seed = seed_coef)$x1)

set.seed(map_seed)
a_t <- a_mult * amp * lomad:::.make_affine_walk(n, 20L, bw = 0.5, cap = cap)
b_t <- 1 + lomad:::.make_affine_walk(n, 20L, bw = 0.5, cap = cap)

# Each structure twice: once without the affine layer, once with. The map runs
# after the mixing, so the middle row is the pre-map state and only x2 differs
# between the two. `dist` has w == 0 throughout, giving the uncoupled pair the
# other columns are drawn against.
trends <- lapply(setNames(nm = names(struct_args)), function(code) {
  base <- list(n, d = d, method = code, seed = seed_coef)
  list(pre  = do.call(sim_trends, c(base, struct_args[[code]])),
       post = do.call(sim_trends, c(base, list(affine_a = a_t, affine_b = b_t),
                                    struct_args[[code]])))
})

y_lim <- range(unlist(lapply(trends, function(x)
  c(x$pre$x1, x$pre$x2, x$post$x2))))
y_pad <- diff(y_lim) * 0.06
y_lim <- y_lim + c(-y_pad, y_pad)

theme_panel <- theme_minimal(base_size = PT$title) +
  theme(
    legend.position  = "none",
    axis.title.x     = element_blank(),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    axis.title.y     = element_text(margin = margin(r = 1)),
    plot.title       = element_text(face = "plain"),
    panel.grid.minor = element_blank()
  ) +
  fig_sizes()

# Only the first column is labelled; every row shares its scale across columns.
bare_y <- theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
                axis.ticks.y = element_blank())

# The uncoupled pair, drawn behind every mixed panel. Mixing pulls the trends in
# toward their mean, so the coloured pair meets the backdrop exactly where
# w_t = 0 and collapses onto the midline where w_t = 1.
backdrop <- data.frame(t = seq_len(n),
                       x1 = trends[["dist"]]$pre$x1,
                       x2 = trends[["dist"]]$pre$x2) |>
  pivot_longer(c(x1, x2), names_to = "series", values_to = "value")

# x1/x2 share the structure's colour: the panel is about the structure, not
# about which series is which.
panel_trends <- function(x1, x2, colour, ylab, bg = FALSE, key = FALSE) {
  p <- data.frame(t = seq_len(n), x1 = x1, x2 = x2) |>
    pivot_longer(c(x1, x2), names_to = "series", values_to = "value") |>
    ggplot(aes(t, value, group = series))
  if (bg)
    p <- p + geom_line(data = backdrop, colour = "grey65",
                       linewidth = 0.35, alpha = 0.55)
  p <- p + geom_line(colour = colour, linewidth = 0.45, alpha = 0.85) +
    scale_y_continuous(limits = y_lim) +
    labs(y = ylab) +
    theme_panel
  # coloured to match the backdrop, so it needs no swatch
  if (key)
    p <- p + annotate("text", x = 0.98 * n, y = y_lim[2],
                      label = "'grey = base trends'~mu[it]", parse = TRUE,
                      hjust = 1, vjust = 1, size = ANNOT, colour = "grey55")
  p
}

panel_wt <- function(tr, colour, title) {
  rng <- range(tr$w)
  pad <- diff(rng) * 0.1
  ggplot(data.frame(t = seq_len(n), w = tr$w), aes(t, w)) +
    geom_hline(yintercept = c(0, 1), linetype = "dashed",
               colour = "gray65", linewidth = 0.3) +
    geom_line(linewidth = 0.4, colour = colour, alpha = 0.85) +
    scale_y_continuous(limits = c(rng[1] - pad, rng[2] + pad), breaks = c(0, 1)) +
    labs(title = title, y = expression(w[t])) +
    theme_panel
}

# Ticks at the limits only: the strips are too short for an interior break to
# read, and the pair of endpoints is all the scale a reader needs here.
panel_coef <- function(y, ylab, limits) {
  ggplot(data.frame(t = seq_len(n), y = y), aes(t, y)) +
    geom_line(linewidth = 0.4, colour = "grey25", alpha = 0.9) +
    scale_y_continuous(breaks = limits, limits = limits, expand = expansion(0.08)) +
    labs(y = ylab) +
    theme_panel
}

struct_title <- function(code) paste0(STRUCT_FULL[[code]], " (", STRUCT_ABBR[[code]], ")")

# The a_t / b_t strips are the same series in all three columns. Repeating them
# above each rescaled panel rather than drawing them once costs nothing and
# keeps every column readable on its own. They lead the bottom half so each
# column reads as map then effect, the order the generator applies them in.
#
# Titles ride on the w_t row so they head the whole column.
#
# nu* is mixed but not yet rescaled; the bottom row's nu is the final trend.
column <- function(code) {
  tr     <- trends[[code]]
  colour <- STRUCT_HEX[[code]]

  p <- panel_wt(tr$pre, colour, struct_title(code)) /
    panel_trends(tr$pre$x1, tr$pre$x2, colour, expression(tilde(nu)[it]),
                 bg = TRUE, key = code == "fr") /
    panel_coef(b_t, expression(b[t]), c(0.2, 1.0)) /
    panel_coef(a_t, expression(a[t]), c(0, 0.08)) /
    panel_trends(tr$post$x1, tr$post$x2, colour, expression(nu[it])) +
    plot_layout(heights = c(1, 3, 0.6, 0.6, 3))

  if (code == "fr") p else p & bare_y
}

# No `widths`: patchwork allocates the null (panel) space and leaves the axis
# furniture outside it, so equal weights already give equal panels. Weighting
# the labelled column up skews them instead.
fig_trends <- column("fr") | column("rs") | column("rm")

ggsave(file.path(IMG_DIR, "fig-trends.png"),
       fig_trends, width = 6.5, height = 4, units = "in", dpi = 450)

})


# =============================================================================
# fig-power.png, panel A -- local power against realized separation
# =============================================================================

Z_LEVEL <- qnorm(0.995)          # pointwise 99%

p_power <- local({
  cv <- readRDS("simulations/power/results/simulations-power-curves.rds") |>
    mutate(Structure = factor(STRUCT_ABBR[struct], levels = STRUCT_ABBR),
           lo = pmax(0, rejection - Z_LEVEL * se),
           hi = pmin(1, rejection + Z_LEVEL * se))

  # delta_t = 1 is an atom, not the end of the mesh: under b > 0 every window
  # whose trends reverse maps there exactly. Drawn as a separate point so the
  # curve is not read as turning up at the boundary.
  cont <- filter(cv, !atom)
  at   <- filter(cv, atom) |> mutate(delta_mean = 1)

  ggplot(cont, aes(delta_mean, rejection, colour = Structure)) +
    geom_hline(yintercept = alpha, linetype = "dashed",
               colour = "grey60", linewidth = 0.4) +
    geom_ribbon(aes(ymin = lo, ymax = hi, fill = Structure),
                alpha = 0.25, colour = NA, show.legend = FALSE) +
    geom_line(linewidth = 0.45, alpha = 0.85) +
    geom_linerange(data = at, aes(ymin = lo, ymax = hi),
                   linewidth = 0.4, alpha = 0.85, show.legend = FALSE) +
    geom_point(data = at, aes(fill = Structure), size = 1.4, shape = 21,
               stroke = 0.35, colour = "white", alpha = 0.85,
               show.legend = FALSE) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    scale_colour_manual(values = STRUCT_PAL) +
    scale_fill_manual(values = STRUCT_PAL) +
    guides(colour = guide_legend(override.aes = list(
             linewidth = 0.5, shape = 16, size = 1.4, linetype = "solid",
             alpha = 1))) +
    facet_nested(phi ~ snr + s, labeller = labeller(
      phi = \(x) paste0("phi == ", x),
      snr = \(x) paste0("SNR == ", x),
      s   = \(x) paste0("s[T] == ", x),
      .default = label_parsed)) +
    labs(x = expression("Effect size " * delta[t]), y = "Rejection rate",
         colour = "Structure") +
    theme_minimal(base_size = PT$title) +
    theme(legend.position = "right",
          legend.key = element_blank(),
          panel.spacing = unit(8, "pt"),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(linewidth = 0.1, colour = "darkgray"))
})


# =============================================================================
# fig-power.png, panel B -- classification accuracy against the decoupling
# threshold, and tbl-localization-auc.csv
# =============================================================================

# Read against panel A, not on its own. Precision and NPV both condition on the
# decision, so they answer "given a flag, was it right" and are silent on what
# was missed. Where power is low the surviving rejections concentrate on the
# largest separations, which holds precision up: the s = 50, phi = 0.7 cells
# score around 0.9 while flagging under 3% of windows.
p_local <- local({
  roc <- readRDS("simulations/power/results/simulations-power-roc.rds") |>
    mutate(Structure = factor(STRUCT_ABBR[struct], levels = STRUCT_ABBR),
           one_npv   = 1 - npv)
  au  <- readRDS("simulations/power/results/simulations-power-auc.rds")

  write.csv(au[, c("struct", "s", "n", "snr", "phi", "windows", "rejection",
                   "auc")],
            file.path(TBL_DIR, "tbl-localization-auc.csv"), row.names = FALSE)
  cat(sprintf("Wrote %s\n", file.path(TBL_DIR, "tbl-localization-auc.csv")))

  # one label per panel: the mean AUC over the three structures
  ann <- au |>
    group_by(phi, snr, s) |>
    summarise(x = 0.97, y = 0.10,
              label = sprintf("AUC = %.3f", mean(auc)), .groups = "drop")

  ggplot(roc, aes(one_npv, prec, colour = Structure)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                colour = "grey70", linewidth = 0.3) +
    geom_path(linewidth = 0.45, alpha = 0.85) +
    geom_text(data = ann, aes(x = x, y = y, label = label), inherit.aes = FALSE,
              hjust = 1, vjust = 0, size = ANNOT, colour = "grey30") +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    scale_colour_manual(values = STRUCT_PAL) +
    facet_nested(phi ~ snr + s, labeller = labeller(
      phi = \(x) paste0("phi == ", x),
      snr = \(x) paste0("SNR == ", x),
      s   = \(x) paste0("s[T] == ", x),
      .default = label_parsed)) +
    labs(x = "1 - NPV", y = "Precision", colour = "Structure") +
    theme_minimal(base_size = PT$title) +
    theme(legend.position = "right",
          legend.key = element_blank(),
          panel.spacing = unit(8, "pt"),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(linewidth = 0.1, colour = "darkgray"))
})


# =============================================================================
# fig-power.png -- the two panel sets stacked
# =============================================================================

local({
  # guides = "collect" merges only identical guides and panel A also maps
  # fill, so suppress on B and C instead -- via guides(), which the trailing
  # `&` cannot override.
  drop_guide <- guides(colour = "none", fill = "none")

  composite <-
    (p_power + guides(fill = "none")) /
    (p_local + drop_guide) +
    plot_layout(guides = "collect", heights = c(1, 1)) +
    plot_annotation(tag_levels = "A") &
    theme(legend.position = "bottom")

  out <- file.path(IMG_DIR, "fig-power.png")
  ggsave(out, composite, width = 6, height = 7, dpi = 450)
  cat(sprintf("\nWrote %s\n", out))
})


# =============================================================================
# fig-validation.png -- finite-sample accuracy of the CLT and of the plug-in
# =============================================================================

local({
  results <- readRDS("simulations/validation/results/simulations-validation-results.rds")

base_theme <- theme_minimal(base_size = PT$title) +
  theme(panel.grid.minor = element_blank()) +
  fig_sizes()

col_th  <- "firebrick"
col_emp <- "grey60"
col_oracle   <- "#0072B2"
col_pipeline <- "#D55E00"

# ---- CLT QQ plots, oracle, faceted by s -------------------------------------

clt_ids <- grep("^clt-", names(results), value = TRUE)

qq_list <- lapply(clt_ids, function(id) {
  r     <- results[[id]]
  R_vec <- r$R_vec[!is.na(r$R_vec)]
  Z     <- sqrt(r$s) * (R_vec - r$rho_oracle) / sqrt(r$V_oracle)
  nn    <- length(Z)

  data.frame(
    theoretical = qnorm(ppoints(nn)),
    empirical   = sort(Z),
    s_label     = sprintf("s = %d", r$s),
    s_num       = r$s
  )
})
qq_df <- do.call(rbind, qq_list)
qq_df$s_label <- reorder(qq_df$s_label, qq_df$s_num)

p_clt <- ggplot(qq_df, aes(theoretical, empirical)) +
  geom_abline(slope = 1, intercept = 0, colour = col_th, linewidth = 0.5) +
  geom_point(colour = "black", alpha = 0.15, size = 0.9) +
  facet_wrap(~ s_label, nrow = 1) +
  labs(x = "N(0, 1) quantile", y = expression(G[t] ~ "quantile")) +
  base_theme

# ---- The two trends ---------------------------------------------------------
# Regenerated, not stored: sim_trends() is deterministic given the seed. The
# scale factor comes from the results rather than being hard-coded, so the
# panel cannot drift from the study that produced the numbers beside it.
# nu_2 = b * nu_1 with b > 0 satisfies H_0, so every window is a true null
# despite the two series differing in amplitude.

b_scale <- results[["e2e-s150"]]$b_scale
if (is.null(b_scale)) stop("results predate b_scale; re-run the validation study")

trend_v <- sim_trends(n = 2000, d = 0, nb = 25, sd0 = 50, p = 1.5, seed = 5381)$x1

# Shared with rho, which is undefined over the first window: without a common
# limit a given t lands at a different x in each panel.
TLIM <- c(1, length(trend_v))

trend_df <- rbind(
  data.frame(t = seq_along(trend_v), nu = trend_v,          k = "1"),
  data.frame(t = seq_along(trend_v), nu = b_scale * trend_v, k = "2")
)

# Line type, not colour. Red is theory and grey is empirical throughout the
# rest of this figure, and blue/orange are oracle and end-to-end; a fifth and
# sixth colour here would collide with all of that. The two trends are the
# same curve up to scale, so a linetype contrast separates them without
# adding to the palette.
# Both solid, same colour, no legend: which line is which follows from the
# title and the amplitudes, so a key would only add ink.
p_trend <- ggplot(trend_df, aes(t, nu, group = k)) +
  geom_line(linewidth = 0.3, colour = col_th) +
  coord_cartesian(xlim = TLIM) +
  base_theme +
  theme(axis.text.x = element_blank(), plot.margin = margin(5.5, 5.5, 6, 5.5),
        legend.position = "none") +
  labs(x = NULL, y = expression(nu[it]),
       title = NULL)

# ---- Proposition 1 moment accuracy, oracle, s = 150 -------------------------

r_rho <- results[["rho-s150"]]
valid <- which(!is.na(r_rho$R_mean) & !is.na(r_rho$rho_th))

rho_df <- data.frame(
  t   = rep(valid, 2),
  rho = c(r_rho$rho_th[valid], r_rho$R_mean[valid]),
  type = rep(c("theoretical", "empirical"), each = length(valid))
)

p_rho <- ggplot() +
  geom_line(data = rho_df[rho_df$type == "empirical", ],
            aes(t, rho), colour = "grey30", linewidth = 0.5) +
  geom_line(data = rho_df[rho_df$type == "theoretical", ],
            aes(t, rho), colour = col_th, linewidth = 0.6) +
  # empty layer, drawn only to build the legend
  geom_line(data = data.frame(
              x = c(NA, NA), y = c(NA, NA),
              label = factor(c("Empirical", "Theoretical"),
                             levels = c("Empirical", "Theoretical"))),
            aes(x, y, colour = label)) +
  scale_colour_manual(values = c("Empirical" = "grey30",
                                  "Theoretical" = col_th)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(NA, 1)) +
  coord_cartesian(xlim = TLIM) +
  labs(x = "Time", y = expression(rho[t]), title = NULL) +
  base_theme +
  theme(
    legend.position = c(1,1),
    legend.justification = c(1, 1),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.title = element_blank(),
    legend.direction = "horizontal",
    legend.key.size = unit(0.4, "cm"),
    plot.margin = margin(3, 5.5, 5.5, 5.5)
  )

r_var   <- results[["var-s150"]]
V_emp   <- r_var$s * apply(r_var$R_mat, 2, var, na.rm = TRUE)
V_theory <- r_var$V_theory
valid_v  <- which(!is.na(V_emp) & !is.na(V_theory) & V_theory > 0)

v_df <- data.frame(
  V_theory = V_theory[valid_v],
  V_emp    = V_emp[valid_v]
)
rng <- range(c(v_df$V_theory, v_df$V_emp))

p_v <- ggplot(v_df, aes(V_theory, V_emp)) +
  geom_abline(slope = 1, intercept = 0, colour = col_th, linewidth = 0.6) +
  geom_point(colour = "black", size = 1.2, alpha = 0.4) +
  coord_equal(xlim = rng, ylim = rng) +
  labs(x = expression("Theoretical" ~ V[t]),
       y = expression(s %.% Var(R[t])),
       title = NULL) +
  base_theme

# ---- End to end, s = 150 ----------------------------------------------------

r_e2e    <- results[["e2e-s150"]]
s_val    <- r_e2e$s
eval_pts <- r_e2e$eval_pts
n_pts    <- length(eval_pts)

# Every evaluation point, not one of them. Showing a single window hid three
# and made the choice of which look arbitrary.
qq_e <- do.call(rbind, lapply(seq_along(eval_pts), function(j) {
  Zo <- sqrt(s_val) * (r_e2e$R_mat[, j] - r_e2e$rho_oracle[j]) / sqrt(r_e2e$V_oracle[j])
  Zo <- Zo[is.finite(Zo)]
  Zp <- r_e2e$Z_est_mat[, j]; Zp <- Zp[is.finite(Zp)]
  nn <- min(length(Zo), length(Zp))
  data.frame(
    theoretical = rep(qnorm(ppoints(nn)), 2),
    empirical   = c(sort(Zo[seq_len(nn)]), sort(Zp[seq_len(nn)])),
    type        = rep(c("Oracle", "End-to-end"), each = nn),
    t           = eval_pts[j]
  )
}))
qq_e$t <- factor(qq_e$t, levels = sort(eval_pts),
                 labels = paste("t =", sort(eval_pts)))

QQ_LIM   <- 3.5
QQ_CLIP_N <- sum(abs(qq_e$empirical) > QQ_LIM)
cat(sprintf("QQ: %d of %d points (%.3f%%) clipped at |Z| > %.1f\n",
            QQ_CLIP_N, nrow(qq_e), 100 * QQ_CLIP_N / nrow(qq_e), QQ_LIM))

p_qq <- ggplot(qq_e, aes(theoretical, empirical, colour = type)) +
  geom_abline(slope = 1, intercept = 0, colour = col_th) +
  geom_point(size = 0.35, alpha = 0.5) +
  facet_wrap(~ t, nrow = 2) +
  scale_colour_manual(values = c(Oracle = col_oracle, `End-to-end` = col_pipeline)) +
  # Square panels on a common symmetric range. A handful of lower-tail points
  # fall outside and are clipped rather than allowed to set the scale for all
  # four facets; QQ_CLIP_N below reports how many, for the caption.
  coord_cartesian(xlim = c(-QQ_LIM, QQ_LIM), ylim = c(-QQ_LIM, QQ_LIM)) +
  labs(x = "N(0, 1) quantile", y = expression(G[t] ~ "quantile")) +
  base_theme +
  theme(aspect.ratio = 1, legend.position = "none")

# lomad_test forms p = pnorm(Z) and rejects when p <= alpha, so only the lower
# tail fires. Two-sided coverage would let an inflated lower tail cancel a
# deflated upper one; at s = 80 in panel A that cancellation is total.
cov_rows <- vector("list", 2 * n_pts)
for (j in seq_len(n_pts)) {
  rho_j <- r_e2e$rho_oracle[j]
  V_j   <- r_e2e$V_oracle[j]

  Zo  <- sqrt(s_val) * (r_e2e$R_mat[, j] - rho_j) / sqrt(V_j)
  Zo  <- Zo[!is.na(Zo)]
  c_o <- mean(Zo < qnorm(alpha))
  n_o <- length(Zo)

  Zp  <- r_e2e$Z_est_mat[, j]
  Zp  <- Zp[!is.na(Zp)]
  c_p <- mean(Zp < qnorm(alpha))
  n_p <- length(Zp)

  cov_rows[[j]]         <- data.frame(
    t = eval_pts[j], cov = c_o,
    lo = c_o - 1.96 * sqrt(c_o * (1 - c_o) / n_o),
    hi = c_o + 1.96 * sqrt(c_o * (1 - c_o) / n_o),
    type = "Oracle"
  )
  cov_rows[[j + n_pts]] <- data.frame(
    t = eval_pts[j], cov = c_p,
    lo = c_p - 1.96 * sqrt(c_p * (1 - c_p) / n_p),
    hi = c_p + 1.96 * sqrt(c_p * (1 - c_p) / n_p),
    type = "End-to-end"
  )
}
cov_df <- do.call(rbind, cov_rows)

dodge <- position_dodge(width = 50)

# The dense grid is computed by the validation runner and travels with the
# results, so the band is reproducible from the repo.
band_df <- r_e2e$band
band_df$se05 <- sqrt(band_df$t05 * (1 - band_df$t05) / band_df$n)

band_layers <- if (!is.null(band_df)) list(
  geom_ribbon(data = band_df, inherit.aes = FALSE,
              aes(x = t, ymin = t05 - 1.96 * se05, ymax = t05 + 1.96 * se05,
                  fill = type), alpha = 0.10),
  geom_line(data = band_df, inherit.aes = FALSE,
            aes(x = t, y = t05, colour = type), alpha = 0.30, linewidth = 0.3),
  scale_fill_manual(values = c(Oracle = col_oracle, `End-to-end` = col_pipeline),
                    guide = "none")
) else list()

# Local SNR over t, above the error panel on a shared axis, as the trend sits
# above rho in row B. The rejection rate tracks SNR, so the strip is what makes
# that panel readable rather than a scatter of four points. Both series
# are drawn: lambda_2 is a constant multiple of lambda_1 by construction, and
# showing only one would imply a single SNR governs the test when the
# experiment cannot separate them.
lam_df <- rbind(
  data.frame(t = seq_along(r_e2e$lambda1), lam = r_e2e$lambda1, k = "1"),
  data.frame(t = seq_along(r_e2e$lambda2), lam = r_e2e$lambda2, k = "2")
)
lam_df <- lam_df[is.finite(lam_df$lam) & lam_df$lam > 0, ]
lam_df$lam <- log10(lam_df$lam)
lam_ratio <- round(median(r_e2e$lambda2 / r_e2e$lambda1, na.rm = TRUE), 2)

p_lam <- ggplot(lam_df, aes(t, lam, group = k)) +
  geom_line(linewidth = 0.3, colour = col_th) +
  # Three breaks, two decades apart, spanning the realized range.
  scale_y_continuous(breaks = c(-3, -1, 1)) +
  coord_cartesian(xlim = TLIM) +
  base_theme +
  theme(axis.text.x = element_blank(), plot.margin = margin(5.5, 5.5, 6, 5.5),
        legend.position = "none") +
  labs(x = NULL, y = expression(log[10] ~ lambda[kt]),
       title = NULL)

p_cov <- ggplot(cov_df, aes(x = t, y = cov, colour = type)) +
  band_layers +
  geom_hline(yintercept = alpha, linetype = "dashed", colour = col_th) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 55,
                position = dodge, linewidth = 0.5) +
  geom_point(size = 1.5, position = dodge) +
  scale_colour_manual(values = c(Oracle = col_oracle, `End-to-end` = col_pipeline)) +
  coord_cartesian(xlim = TLIM) +
  labs(x = "Time", y = "Type I error") +
  base_theme +
  theme(
    # The two panels of row C share a colour scale, so one key serves both. It
    # goes here rather than on the QQ, whose square panels have no free corner.
    # Nudged in from the right edge: the t = 1900 intervals reach 0.085 and
    # the key sat on them.
    legend.position = c(0.86, 0.99),
    legend.justification = c(1, 1),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.title = element_blank(),
    legend.key.size = unit(0.3, "cm"),
    legend.text = element_text(size = 7)
  ) +
  guides(colour = guide_legend(override.aes = list(size = 1.2)))

# ---- Composite --------------------------------------------------------------
# Three rows of equal height, trend flush above rho on a shared time axis.

full_fig <- (p_clt + labs(tag = "A")) + (p_trend + labs(tag = "B")) + p_rho +
  p_v + (p_lam + labs(tag = "C")) + p_cov + p_qq +
  plot_layout(design = c(
    area(1,  1,  6, 6),   # A  CLT QQ facets
    area(7,  1,  9, 4),   # B  trend
    area(10, 1, 12, 4),   #    rho
    area(7,  5, 12, 6),   #    V
    # Row C runs taller than A and B: the QQ facets are square, so their size
    # is set by whichever of width or height binds first, and at an equal span
    # it was height -- the panels came out a third of their width.
    area(13,  1, 15, 4),  # C  local SNR, flush above coverage
    area(16,  1, 21, 4),  #    coverage, on panel B's time axis
    area(13,  5, 21, 6)   #    end-to-end QQ, faceted 2x2 by evaluation point
  ))

ggsave(file.path(IMG_DIR, "fig-validation.png"),
         full_fig, width = 6, height = 6.9, dpi = 450)
})

cat("All figures written to ", IMG_DIR, "\n", sep = "")
