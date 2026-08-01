## simulation-results.R -- figures and tables for the simulation studies
##
## The cheap stage of the pipeline: reads compiled summaries and precomputed
## intermediates only, never the raw per-job files, so it runs in seconds and
## can be re-run freely while drafting. The expensive stages are
## simulations/<study>/collect-results.R and
## simulations/power/localization-sweep.R.
##
## Inputs
##   simulations/power/results/simulations-power-summary.rds
##   simulations/power/results/simulations-power-localization.rds
##   simulations/validation/results/simulations-validation-results.rds
##   (fig-trends.png simulates its own data)
##
## Outputs -> simulations/_img/
##   fig-trends.png             methods of simulating trend separation
##   fig-power-composite.png    power, concordance and profile stacked
##   fig-validation.png         finite-sample accuracy of the CLT
##
## Outputs -> simulations/_tbl/
##   tbl-localization-auc.csv         concordance AUC per design cell
##   tbl-localization-resolution.csv  separation at which rejection hits .50/.95
##
## Both output directories are _-prefixed and untracked: everything here
## regenerates in seconds from the compiled results, which are tracked.
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
STRUCT_FULL <- c(rate = "Fixed Rate", smooth = "Random Separation",
                  cross = "Random Mixing")
STRUCT_ABBR <- c(rate = "FR", smooth = "RS", cross = "RM")
STRUCT_HEX  <- c(rate = "#009E73", smooth = "#0072B2", cross = "#D55E00")
STRUCT_PAL  <- setNames(STRUCT_HEX, STRUCT_ABBR[names(STRUCT_HEX)])

# Mirrors localization-sweep.R. The power summary carries T but not s_T.
h_win_fcn <- function(n) max(5L, floor(n / 200L))
s_win_fcn <- function(n) min(60L * h_win_fcn(n), floor(n / 4L))

ESTIMATION_LAB <- c(estimated = "Lomad", oracle = "Oracle")
ESTIMATION_LTY <- c(estimated = "solid", oracle = "dashed")


# =============================================================================
# fig-trends.png -- ways of distributing separation across a series pair
# =============================================================================

local({

n        <- 500
d        <- 2       # L2 separation, held constant across all four panels
bw       <- 50      # bandwidth b
coupling <- 0.8     # coupling fraction c
rate     <- 0.01    # event rate r

seed_coef <- 2847   # Fourier base, shared across panels

tr_dist   <- sim_trends(n, d = d, method = "dist",   seed = seed_coef)
tr_rate   <- sim_trends(n, d = d, method = "rate",   seed = seed_coef,
                        rate = rate, bump = "gaussian")
tr_smooth <- sim_trends(n, d = d, method = "smooth", seed = seed_coef,
                        bw = bw, coupling = coupling)
tr_cross  <- sim_trends(n, d = d, method = "cross",  seed = seed_coef,
                        bw = bw, coupling = coupling)

y_lim <- range(c(
  tr_dist$x1,   tr_dist$x2,
  tr_rate$x1,   tr_rate$x2,   tr_rate$x_mean,
  tr_smooth$x1, tr_smooth$x2, tr_smooth$x_mean,
  tr_cross$x1,  tr_cross$x2,  tr_cross$x_mean
))
y_pad <- diff(y_lim) * 0.05
y_lim <- y_lim + c(-y_pad, y_pad)

theme_top <- theme_minimal(base_size = PT$title) +
  theme(
    legend.position  = "none",
    axis.title.x     = element_blank(),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    plot.title       = element_text(face = "plain"),
    panel.grid.minor = element_blank()
  ) +
  fig_sizes()

theme_bot <- theme_minimal(base_size = PT$title) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank()
  ) +
  fig_sizes()

# w = 0 here, so x1/x2 are the raw Fourier trends and share one colour.
panel_top_dist <- function(tr, title) {
  x_mean_loc <- (tr$x1 + tr$x2) / 2
  df <- data.frame(t = seq_len(n), x1 = tr$x1, x2 = tr$x2,
                   x_mean = x_mean_loc) |>
    pivot_longer(c(x1, x2, x_mean), names_to = "series", values_to = "value") |>
    mutate(series = factor(series, levels = c("x_mean", "x1", "x2")))
  ggplot(df, aes(t, value, color = series, linewidth = series,
                 linetype = series, alpha = series)) +
    geom_line() +
    scale_color_manual(values = c(x_mean = "gray50", x1 = "black", x2 = "black")) +
    scale_linewidth_manual(values = c(x_mean = 0.35, x1 = 0.45, x2 = 0.45)) +
    scale_linetype_manual(values = c(x_mean = "dashed", x1 = "solid", x2 = "solid")) +
    scale_alpha_manual(values = c(x_mean = 1, x1 = 0.9, x2 = 0.9)) +
    scale_y_continuous(limits = y_lim) +
    labs(title = title, y = expression(mu[it])) +
    theme_top
}

# x1/x2 share the structure's colour: the panel is about the structure, not
# about which series is which.
panel_top_struct <- function(tr, title, colour) {
  df <- data.frame(t = seq_len(n), x1 = tr$x1, x2 = tr$x2,
                   x_mean = tr$x_mean) |>
    pivot_longer(c(x1, x2, x_mean), names_to = "series", values_to = "value") |>
    mutate(series = factor(series, levels = c("x_mean", "x1", "x2")))
  ggplot(df, aes(t, value, color = series, linewidth = series,
                 linetype = series, alpha = series)) +
    geom_line() +
    scale_color_manual(values = c(x_mean = "gray50", x1 = colour, x2 = colour)) +
    scale_linewidth_manual(values = c(x_mean = 0.35, x1 = 0.45, x2 = 0.45)) +
    scale_linetype_manual(values = c(x_mean = "dashed", x1 = "solid", x2 = "solid")) +
    scale_alpha_manual(values = c(x_mean = 1, x1 = 0.9, x2 = 0.9)) +
    scale_y_continuous(limits = y_lim) +
    labs(title = title, y = expression(nu[it])) +
    theme_top
}

panel_wt <- function(tr, ref_lines = c(0, 1), ylim = NULL, colour = "gray20") {
  if (is.null(ylim)) {
    rng  <- range(tr$w)
    pad  <- diff(rng) * 0.1
    ylim <- c(rng[1] - pad, rng[2] + pad)
  }
  data.frame(t = seq_len(n), w = tr$w) |>
    ggplot(aes(t, w)) +
    geom_hline(yintercept = ref_lines, linetype = "dashed",
               color = "gray65", linewidth = 0.3) +
    geom_line(linewidth = 0.4, color = colour) +
    scale_y_continuous(limits = ylim, breaks = ref_lines) +
    labs(y = expression(w[t]), x = "t") +
    theme_bot
}

# Abbreviation on its own line: "Random Separation (RS)" on one line is wider
# than the panel and ggplot truncates it silently. The unstructured panel takes
# a blank second line so all four titles are the same height.
struct_title <- function(code) paste0(STRUCT_FULL[[code]], "\n(", STRUCT_ABBR[[code]], ")")

comp_dist <- panel_top_dist(tr_dist, "Base trends\n") /
  panel_wt(tr_dist, ref_lines = c(0, 1), ylim = c(-0.1, 1.1)) +
  plot_layout(heights = c(3, 1))

comp_rate <- panel_top_struct(tr_rate, struct_title("rate"),
                               colour = STRUCT_HEX[["rate"]]) /
  panel_wt(tr_rate, colour = STRUCT_HEX[["rate"]]) +
  plot_layout(heights = c(3, 1))

comp_smooth <- panel_top_struct(tr_smooth, struct_title("smooth"),
                                 colour = STRUCT_HEX[["smooth"]]) /
  panel_wt(tr_smooth, colour = STRUCT_HEX[["smooth"]]) +
  plot_layout(heights = c(3, 1))

comp_cross <- panel_top_struct(tr_cross, struct_title("cross"),
                                colour = STRUCT_HEX[["cross"]]) /
  panel_wt(tr_cross, ref_lines = c(0, 1), colour = STRUCT_HEX[["cross"]]) +
  plot_layout(heights = c(3, 1))

fig_1x4 <- comp_dist | comp_rate | comp_smooth | comp_cross
ggsave(file.path(IMG_DIR, "fig-trends.png"),
         fig_1x4, width = 6.5, height = 2, units = "in", dpi = 400)

})


# =============================================================================
# fig-power-composite.png, panel A -- power against L2 separation d
# =============================================================================

p_power <- local({
  results_summary <- readRDS("simulations/power/results/simulations-power-summary.rds")

# Oracle runs exist only at phi = 0.8, so the solid/dashed contrast gets an
# inline note there rather than a legend covering panels it does not apply to.
oracle_note <- data.frame(phi = 0.8, snr = 0.5, n = 200,
                           d = 0.05, detection = 0.97,
                           label = "dashed =\noracle")   # one line overruns

p <- results_summary |>
  mutate(struct = factor(STRUCT_ABBR[struct], levels = STRUCT_ABBR)) |>
  ggplot(aes(d, detection, colour = struct,
           linetype = method,
           group = interaction(struct, method))) +
  geom_hline(yintercept = alpha, linetype = "dashed",
             colour = "grey60", linewidth = 0.5) +
  geom_line(linewidth = 0.5, alpha = 0.8) +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi, fill = struct),
              alpha = 0.2, colour = NA) +
  geom_text(data = oracle_note, aes(x = d, y = detection, label = label),
            inherit.aes = FALSE, hjust = 0, vjust = 1, lineheight = 0.95,
            size = ANNOT, colour = "grey30") +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  scale_x_continuous(breaks = 0:2) +
  scale_colour_manual(values = STRUCT_PAL) +
  scale_fill_manual(values = STRUCT_PAL) +
  scale_linetype_manual(values = ESTIMATION_LTY) +
  guides(linetype = guide_none()) +
  facet_nested(phi ~ snr + n, labeller = labeller(
    phi = \(x) paste0("phi == ", x),
    snr = \(x) paste0("SNR == ", x),
    n   = \(x) paste0("s[T] == ", vapply(as.integer(x), s_win_fcn, numeric(1))),
    .default = label_parsed
  )) +
  labs(x = "Separation (d)", y = "Power",
       colour = "Structure", fill = "Structure") +
  theme_minimal(base_size = PT$title) +
  theme(legend.position = "right",
        panel.spacing = unit(8, "pt"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.1, color = "darkgray"))

  p
})


# =============================================================================
# fig-power-composite.png, panel B -- concordance between rejections and
# true local separation
# =============================================================================

p_local <- local({
  L      <- readRDS("simulations/power/results/simulations-power-localization.rds")
  sweep  <- L$sweep
  ORIENT <- L$orient
  C_MAX   <- 0.30
  MIN_N   <- 10000L

  if (ORIENT == "conventional") {
    sweep$xx <- 1 - sweep$spec; sweep$yy <- sweep$sens
    xlab <- "1 - Specificity"
    ylab <- "Sensitivity"
  } else {
    sweep$xx <- 1 - sweep$npv;  sweep$yy <- sweep$prec
    xlab <- "1 - NPV"
    ylab <- "Precision"
  }
sw <- sweep |>
  filter(c <= C_MAX, n_above >= MIN_N, n_below >= MIN_N) |>
  mutate(Structure = factor(STRUCT_ABBR[struct], levels = STRUCT_ABBR),
         method = factor(method, levels = c("estimated", "oracle")))

# Integrated over the full c grid, not the trimmed display range: a partial
# area would not read as a concordance.
auc <- sweep |>
  group_by(struct, phi, snr, s_win, method) |>
  arrange(c, .by_group = TRUE) |>
  summarise(auc = {
    x <- 1 - npv; y <- prec; o <- order(x)
    sum(diff(x[o]) * (y[o][-1] + head(y[o], -1)) / 2, na.rm = TRUE)
  }, .groups = "drop")

write.csv(auc, file.path(TBL_DIR, "tbl-localization-auc.csv"), row.names = FALSE)
cat(sprintf("Wrote %s\n", file.path(TBL_DIR, "tbl-localization-auc.csv")))

# Mean over exactly the solid curves drawn in each panel: the semi-join keeps
# the annotation tied to what is visible if the filters above change.
auc_panel <- auc |>
  filter(method == "estimated") |>
  semi_join(distinct(sw, struct, phi, snr, s_win, method),
            by = c("struct", "phi", "snr", "s_win", "method")) |>
  group_by(phi, snr, s_win) |>
  summarise(xx = 0.95, yy = 0.2, k = n(),
            label = sprintf("AUC = %.3f", mean(auc)), .groups = "drop")

p <- ggplot(sw, aes(xx, yy, colour = Structure, linetype = method,
                    group = interaction(struct, method))) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted",
              colour = "grey70", linewidth = 0.3) +
  geom_path(linewidth = 0.6, alpha = 0.85) +
  geom_text(data = auc_panel, aes(x = xx, y = yy, label = label),
            inherit.aes = FALSE, hjust = 1, size = ANNOT, colour = "grey20") +
  facet_nested(phi ~ snr + s_win, labeller = labeller(
    phi   = function(x) paste0("phi == ", x),
    snr   = function(x) paste0("SNR == ", x),
    s_win = function(x) paste0("s[T] == ", x),
    .default = label_parsed)) +
  scale_colour_manual(values = STRUCT_PAL) +
  scale_linetype_manual(values = ESTIMATION_LTY) +
  guides(linetype = guide_none()) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  labs(x = xlab, y = ylab, colour = "Structure") +
  theme_minimal(base_size = PT$title) +
  theme(legend.position = "right",
        panel.spacing = unit(8, "pt"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.1, color = "darkgray"))

p
})


# =============================================================================
# fig-power-composite.png, panel C -- rejection probability against true
# windowed separation, and tbl-localization-resolution.csv
# =============================================================================

p_profile <- local({
  sweep <- readRDS("simulations/power/results/simulations-power-localization.rds")$sweep

  BW    <- 0.01               # separation bin width
  C_MAX <- 0.30               # matches panel B's plotted range
  EDGES <- seq(BW, C_MAX, by = BW)

  # sens(c) * n_above(c) counts rejected windows above the cut, so differencing
  # adjacent cuts gives exact within-bin counts.
  prof <- sweep |>
    mutate(A = sens * n_above) |>
    filter(c %in% round(EDGES, 3)) |>
    group_by(struct, phi, snr, s_win, method) |>
    arrange(c, .by_group = TRUE) |>
    reframe(lo    = head(c, -1),
            n_bin = head(n_above, -1) - tail(n_above, -1),
            n_rej = head(A, -1)       - tail(A, -1)) |>
    mutate(mid = lo + BW / 2, rate = n_rej / n_bin) |>
    mutate(Structure = factor(STRUCT_ABBR[struct], levels = STRUCT_ABBR),
           method    = factor(method, levels = c("estimated", "oracle")))

  # No error bars: neighbouring windows overlap almost completely, so effective
  # sample size is far below the bin counts and a binomial interval would be
  # badly overconfident.
  p <- ggplot(prof, aes(mid, rate, colour = Structure, linetype = method,
                        group = interaction(struct, method))) +
    geom_hline(yintercept = 0.5, linetype = "dotted", colour = "grey70",
               linewidth = 0.3) +
    geom_line(linewidth = 0.6, alpha = 0.9) +
    facet_nested(phi ~ snr + s_win, labeller = labeller(
      phi   = function(x) paste0("phi == ", x),
      snr   = function(x) paste0("SNR == ", x),
      s_win = function(x) paste0("s[T] == ", x),
      .default = label_parsed)) +
    scale_colour_manual(values = STRUCT_PAL) +
    scale_linetype_manual(values = ESTIMATION_LTY) +
    guides(linetype = guide_none()) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    labs(x = "Windowed maximum separation", y = "Rejection probability",
         colour = "Structure") +
    theme_minimal(base_size = PT$title) +
    theme(legend.position = "right",
            panel.grid.minor = element_blank(),
          panel.grid.major = element_line(linewidth = 0.1, color = "darkgray"))

  res <- prof |>
    group_by(struct, phi, snr, s_win, method) |>
    summarise(sep_50 = if (any(rate >= 0.50)) mid[which(rate >= 0.50)[1]] else NA_real_,
              sep_95 = if (any(rate >= 0.95)) mid[which(rate >= 0.95)[1]] else NA_real_,
              .groups = "drop")

  write.csv(res, file.path(TBL_DIR, "tbl-localization-resolution.csv"),
            row.names = FALSE)
  cat(sprintf("Wrote %s\n", file.path(TBL_DIR, "tbl-localization-resolution.csv")))

  cat("\nSeparation at which rejection probability reaches 0.50 (NA: not within",
      C_MAX, "):\n")
  print(as.data.frame(res |> filter(method == "estimated") |>
    mutate(sep_50 = ifelse(is.na(sep_50), "--", sprintf("%.3f", sep_50))) |>
    tidyr::pivot_wider(names_from = s_win, values_from = sep_50,
                       names_prefix = "s=", id_cols = c(struct, snr, phi)) |>
    arrange(struct, snr, phi)), row.names = FALSE)

  p
})


# =============================================================================
# fig-power-composite.png -- the three panel sets stacked
# =============================================================================

local({
  # guides = "collect" merges only identical guides, and panel A maps fill as
  # well as colour, so its guide never matches. Suppress on B and C instead,
  # via guides() so the trailing `&` cannot override it.
  drop_guide <- guides(colour = "none", fill = "none")

  composite <-
    (p_power   + guides(fill = "none")) /
    (p_local   + drop_guide) /
    (p_profile + drop_guide) +
    plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "A") &
    theme(legend.position = "bottom")

  out <- file.path(IMG_DIR, "fig-power-composite.png")
  ggsave(out, composite, width = 6.5, height = 8, dpi = 450)
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
  labs(x = "Theoretical N(0,1)", y = expression("Empirical" ~ Z[t])) +
  base_theme

# ---- The common trend -------------------------------------------------------
# Regenerated rather than stored: sim_trends() is deterministic given the seed.
# d = 0, so both series share this trend and every window is a true null.

trend_v <- sim_trends(n = 2000, d = 0, nb = 25, sd0 = 50, p = 1.5, seed = 5381)$x1

# Shared with rho below. rho is undefined over the first window, so without a
# common limit a given t would land at a different x in each panel.
TLIM <- c(1, length(trend_v))

p_trend <- ggplot(data.frame(t = seq_along(trend_v), nu = trend_v),
                  aes(t, nu)) +
  geom_line(colour = col_th, linewidth = 0.3) +
  coord_cartesian(xlim = TLIM) +
  scale_y_continuous(breaks = c(-2, 0, 2)) +
  base_theme +
  theme(axis.text.x = element_blank(), plot.margin = margin(5.5, 5.5, 0, 5.5)) +
  labs(x = NULL, y = expression(nu[t]), title = "Shared trend")

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
  labs(x = "Time", y = expression(rho[t]), title = "Local correlation (s = 150)") +
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
       y = expression(s %.% Var(R[t]))) +
  base_theme

# ---- End to end, s = 150 ----------------------------------------------------

r_e2e    <- results[["e2e-s150"]]
s_val    <- r_e2e$s
eval_pts <- r_e2e$eval_pts
n_pts    <- length(eval_pts)

j_mid    <- which(eval_pts == 1000)
R_j      <- r_e2e$R_mat[, j_mid]
rho_or   <- r_e2e$rho_oracle[j_mid]
V_or     <- r_e2e$V_oracle[j_mid]
Z_oracle <- sqrt(s_val) * (R_j - rho_or) / sqrt(V_or)
Z_oracle <- Z_oracle[!is.na(Z_oracle)]
Z_pipe   <- r_e2e$Z_est_mat[, j_mid]
Z_pipe   <- Z_pipe[!is.na(Z_pipe)]

nn_e <- min(length(Z_oracle), length(Z_pipe))
qq_e <- data.frame(
  theoretical = rep(qnorm(ppoints(nn_e)), 2),
  empirical   = c(sort(Z_oracle[seq_len(nn_e)]), sort(Z_pipe[seq_len(nn_e)])),
  type        = rep(c("Oracle", "End to end"), each = nn_e)
)

p_qq <- ggplot(qq_e, aes(theoretical, empirical, colour = type)) +
  geom_abline(slope = 1, intercept = 0, colour = col_th) +
  geom_point(size = 0.6, alpha = 0.6) +
  scale_colour_manual(values = c(Oracle = col_oracle, `End to end` = col_pipeline)) +
  labs(x = "Theoretical N(0,1)", y = expression("Empirical" ~ Z[t])) +
  base_theme +
  theme(
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.title = element_blank(),
    legend.key.size = unit(0.35, "cm")
  ) +
  guides(colour = guide_legend(override.aes = list(size = 1.5, alpha = 1)))

cov_rows <- vector("list", 2 * n_pts)
for (j in seq_len(n_pts)) {
  rho_j <- r_e2e$rho_oracle[j]
  V_j   <- r_e2e$V_oracle[j]

  Zo  <- sqrt(s_val) * (r_e2e$R_mat[, j] - rho_j) / sqrt(V_j)
  Zo  <- Zo[!is.na(Zo)]
  c_o <- mean(abs(Zo) < qnorm(0.975))
  n_o <- length(Zo)

  Zp  <- r_e2e$Z_est_mat[, j]
  Zp  <- Zp[!is.na(Zp)]
  c_p <- mean(abs(Zp) < qnorm(0.975))
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
    type = "End to end"
  )
}
cov_df <- do.call(rbind, cov_rows)

dodge <- position_dodge(width = 50)

p_cov <- ggplot(cov_df, aes(x = t, y = cov, colour = type)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", colour = col_th) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 55,
                position = dodge, linewidth = 0.5) +
  geom_point(size = 1.5, position = dodge) +
  scale_colour_manual(values = c(Oracle = col_oracle, `End to end` = col_pipeline)) +
  labs(x = "Time", y = "95% coverage") +
  base_theme +
  theme(legend.position = "none")

# ---- Composite --------------------------------------------------------------
# Three rows of equal height. The trend sits flush on top of rho, sharing its
# time axis, so rho can be read against the trend that generates it.

full_fig <- (p_clt + labs(tag = "A")) + (p_trend + labs(tag = "B")) + p_rho +
  p_v + (p_qq + labs(tag = "C")) + p_cov +
  plot_layout(design = c(
    area(1,  1,  6, 6),   # A  CLT QQ facets
    area(7,  1,  9, 4),   # B  trend
    area(10, 1, 12, 4),   #    rho
    area(7,  5, 12, 6),   #    V
    area(13, 1, 18, 3),   # C  end-to-end QQ
    area(13, 4, 18, 6)    #    coverage
  ))

ggsave(file.path(IMG_DIR, "fig-validation.png"),
         full_fig, width = 6, height = 5.25, dpi = 300)
})

cat("All figures written to ", IMG_DIR, "\n", sep = "")
