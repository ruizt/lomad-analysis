## simulation-figures.R — every figure in the paper, from compiled results
##
## This is the cheap stage of the pipeline. It reads only compiled summaries
## and precomputed intermediates, never the raw per-job files, so it runs in
## seconds and can be re-run freely while drafting. The expensive stages are
## simulations/<study>/collect-results.R (assemble + summarise) and
## simulations/power/localization-sweep.R (the localization sweep).
##
## Inputs
##   simulations/power/results/simulations-power-summary.rds
##   simulations/validation/results/simulations-validation-results.rds
##   simulations/power/results/simulations-power-localization.rds
##   (the trend-construction figure needs no inputs; it simulates its own)
##
## Outputs -> simulations/_img/
##   fig-trends.png                methods of simulating trend separation
##   fig-power.png                 detection rate vs separation d
##   fig-localization.png          localization threshold sweep
##   fig-validation.png            finite-sample accuracy of the CLT
##
## Usage (from the repo root):
##   Rscript simulations/simulation-figures.R

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

IMG_DIR <- "simulations/_img"
alpha   <- 0.05
dir.create(IMG_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- Shared structure naming/colour convention (all figures) ---------------
# Internal codes (as stored in the data) map to the same display name,
# abbreviation, and colour everywhere: trend construction, power curves, and
# localization. Keep the palette keyed by the internal code (not by a display
# label) — a display label used as a vector name gets silently mangled if it
# is itself named when spliced into another named vector via c().
STRUCT_FULL <- c(rate = "Fixed Rate", smooth = "Stochastic Modulation",
                  cross = "Stochastic Blending")
STRUCT_ABBR <- c(rate = "FR", smooth = "SM", cross = "SB")
STRUCT_HEX  <- c(rate = "#009E73", smooth = "#0072B2", cross = "#D55E00")

# Full display label used only for the trend-construction panel titles.
STRUCT_LABELS <- STRUCT_FULL
# Abbreviation-keyed palette, for the Structure legend in the other figures.
STRUCT_PAL <- setNames(STRUCT_HEX, STRUCT_ABBR[names(STRUCT_HEX)])

# Rolling-window length s_T as a function of series length T. Mirrors
# localization-sweep.R exactly; fig-power.png doesn't carry a window-size
# column of its own, so T's s_T is derived here for the facet labels.
h_win_fcn <- function(n) max(5L, floor(n / 200L))
s_win_fcn <- function(n) min(60L * h_win_fcn(n), floor(n / 4L))

# Shared legend for the solid/dashed (estimated vs. oracle noise) contrast,
# used identically in the power-curve and localization figures.
ESTIMATION_LAB <- c(estimated = "Lomad", oracle = "Oracle")
ESTIMATION_LTY <- c(estimated = "solid", oracle = "dashed")


# =============================================================================
# Trend construction
# Methods of distributing separation across the series (paper Fig. 1).
# =============================================================================

local({

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

# Unstructured panel: w = 0, so x1/x2 are the raw Fourier trends. Both lines
# black — there is no structure/colour to distinguish here.
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
    labs(title = title, y = expression(nu[it])) +
    theme_top
}

# Structured panel: show x1, x2, and x_mean. x1/x2 share a single colour
# (the structure's colour, from STRUCT_PAL) rather than being distinguished
# from each other — the panel is about the structure, not which series is
# which.
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

# Coupling weight sub-panel. Colour-matched to the structure above it (gray
# for the unstructured panel, where there is no structure colour).
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

# ---- Composite panels (series / w_t) -----------------------------------

struct_title <- function(code) paste0(STRUCT_FULL[[code]], " (", STRUCT_ABBR[[code]], ")")

comp_dist <- panel_top_dist(tr_dist, "Unstructured") /
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

# ---- 2×2 figure --------------------------------------------------------

fig_1x4 <- comp_dist | comp_rate | comp_smooth | comp_cross
ggsave(file.path(IMG_DIR, "fig-trends.png"),
         fig_1x4, width = 8.6, height = 2.5, units = "in", dpi = 400)

})


# =============================================================================
# Power curves
# Detection rate as a function of L2 separation d, by structure/phi/SNR/T.
# =============================================================================

local({
  results_summary <- readRDS("simulations/power/results/simulations-power-summary.rds")

library(ggh4x)

# Oracle only exists at phi = 0.8 (dashed vs. solid is meaningless elsewhere
# in this grid), so it gets a one-off inline note there instead of a legend
# that would otherwise apply, misleadingly, to every panel.
oracle_note <- data.frame(phi = 0.8, snr = 0.5, n = 200,
                           d = 0.05, detection = 0.97,
                           label = "dashed = oracle")

results_summary |>
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
            inherit.aes = FALSE, hjust = 0, size = 2.8, colour = "grey30") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_colour_manual(values = STRUCT_PAL) +
  scale_fill_manual(values = STRUCT_PAL) +
  scale_linetype_manual(values = ESTIMATION_LTY) +
  guides(linetype = guide_none()) +
  facet_nested(phi ~ snr + n, labeller = labeller(
    phi = \(x) paste0("phi == ", x),
    snr = \(x) paste0("SNR == ", x),
    n   = \(x) paste0("T==", x, "*', '~s[T]==", vapply(as.integer(x), s_win_fcn, numeric(1))),
    .default = label_parsed
  )) +
  labs(x = "Separation (d)", y = "Power",
       colour = "Structure", fill = "Structure") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right",
        axis.text = element_text(size = 8),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.1, color = "darkgray"))


ggsave(file.path(IMG_DIR, "fig-power.png"),
         width = 9, height = 4, dpi = 450)

})


# =============================================================================
# Localization
# Threshold sweep relating rejections to true local separation.
# =============================================================================

local({
  L      <- readRDS("simulations/power/results/simulations-power-localization.rds")
  sweep  <- L$sweep
  ORIENT <- L$orient
  C_MAX   <- 0.30
  MIN_N   <- 10000L
  # The rolling-max estimator's localization collapses to near-chance only for
  # "FR" in the two hardest cells (phi = 0.8, SNR = 1.5, s_T = 100 and 150) —
  # dropped there and only there; FR is kept everywhere else.
  DROP_STRUCT <- "rate"; DROP_PHI <- 0.8; DROP_SNR <- 1.5
  DROP_S_WIN  <- c(100L, 150L)

  # Axis pair depends on the orientation recorded by localization-sweep.R.
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
  filter(!(struct == DROP_STRUCT & phi == DROP_PHI & snr == DROP_SNR &
             s_win %in% DROP_S_WIN)) |>
  mutate(Structure = factor(STRUCT_ABBR[struct], levels = STRUCT_ABBR),
         method = factor(method, levels = c("estimated", "oracle")))

# Oracle only exists at phi = 0.8 here too (same as the power-curve figure),
# so it gets the same one-off inline note instead of a legend that would
# otherwise apply, misleadingly, to every panel.
oracle_note_loc <- data.frame(phi = 0.8, snr = 0.5, s_win = 50,
                               xx = 0.95, yy = 0.05, label = "dashed = oracle")

p <- ggplot(sw, aes(xx, yy, colour = Structure, linetype = method,
                    group = interaction(struct, method))) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted",
              colour = "grey70", linewidth = 0.3) +
  geom_path(linewidth = 0.6, alpha = 0.85) +
  geom_text(data = oracle_note_loc, aes(x = xx, y = yy, label = label),
            inherit.aes = FALSE, hjust = 1, size = 2.8, colour = "grey30") +
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
  theme_minimal(base_size = 12) +
  theme(legend.position = "right",
        axis.text = element_text(size = 8),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.1, color = "darkgray"))

ggsave(file.path(IMG_DIR, "fig-localization.png"), p,
       width = 9, height = 4, dpi = 450)
cat(sprintf("\nWrote %s\n", file.path(IMG_DIR, "fig-localization.png")))

# ---- Concordance AUC ---------------------------------------------------------
# Area under the curve above, which in the predictive orientation is a genuine
# ROC: true separation m_t is the score, rejection status the class label. So
# the area reads as P(a rejected window has larger m_t than an unrejected one).
#
# Integrated over the FULL c grid, not the plotted range: C_MAX and MIN_N are
# display choices that trim the tails, and a partial area would not carry the
# concordance reading. The curve is anchored at both corners without them --
# every window is above the cut as c -> 0, none as c -> max.
auc <- sweep |>
  group_by(struct, phi, snr, s_win, method) |>
  arrange(c, .by_group = TRUE) |>
  summarise(auc = {
    x <- 1 - npv; y <- prec; o <- order(x)
    sum(diff(x[o]) * (y[o][-1] + head(y[o], -1)) / 2, na.rm = TRUE)
  }, .groups = "drop")

cat("\nLocalization concordance AUC:\n")
print(as.data.frame(auc |>
  mutate(auc = sprintf("%.3f", auc)) |>
  tidyr::pivot_wider(names_from = s_win, values_from = auc,
                     names_prefix = "s=") |>
  arrange(method, struct, snr, phi)), row.names = FALSE)

})


# =============================================================================
# Validation composite
# Finite-sample accuracy of the CLT approximation and the plug-in pipeline.
# =============================================================================

local({
  results <- readRDS("simulations/validation/results/simulations-validation-results.rds")
# ---- Shared theme ------------------------------------------------------------

base_theme <- theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey92")
  )

col_th  <- "firebrick"
col_emp <- "grey60"
col_oracle   <- "#0072B2"
col_pipeline <- "#D55E00"

# ==============================================================================
# Row (a): CLT QQ plots (oracle, faceted by s)
# ==============================================================================

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

row_a <- ggplot(qq_df, aes(theoretical, empirical)) +
  geom_abline(slope = 1, intercept = 0, colour = col_th, linewidth = 0.5) +
  geom_point(colour = "black", alpha = 0.15, size = 0.9) +
  facet_wrap(~ s_label, nrow = 1) +
  labs(x = "Theoretical N(0,1)", y = expression("Empirical" ~ Z[t])) +
  base_theme

# ==============================================================================
# Row (b): Proposition 1 moment accuracy (oracle, s = 150)
# ==============================================================================

# ---- Left panel: rho --------------------------------------------------------

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
  # Dummy layer for legend
  geom_line(data = data.frame(
              x = c(NA, NA), y = c(NA, NA),
              label = factor(c("Empirical", "Theoretical"),
                             levels = c("Empirical", "Theoretical"))),
            aes(x, y, colour = label)) +
  scale_colour_manual(values = c("Empirical" = "grey30",
                                  "Theoretical" = col_th)) +
  labs(x = "Time", y = expression(rho[t])) +
  base_theme +
  theme(
    legend.position = c(1,1),
    legend.justification = c(1, 1),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.4, "cm")
  )

# ---- Right panel: V ----------------------------------------------------------

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

row_b <- p_rho + p_v

# ==============================================================================
# Row (c): End-to-end pipeline validation (s = 150)
# ==============================================================================

r_e2e    <- results[["e2e-s150"]]
s_val    <- r_e2e$s
eval_pts <- r_e2e$eval_pts
n_pts    <- length(eval_pts)

# ---- Left panel: QQ at t = 1000 ---------------------------------------------

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
  type        = rep(c("Oracle", "Pipeline"), each = nn_e)
)

p_qq <- ggplot(qq_e, aes(theoretical, empirical, colour = type)) +
  geom_abline(slope = 1, intercept = 0, colour = col_th, linetype = "dashed") +
  geom_point(size = 0.6, alpha = 0.6) +
  scale_colour_manual(values = c(Oracle = col_oracle, Pipeline = col_pipeline)) +
  labs(x = "Theoretical N(0,1)", y = expression("Empirical" ~ Z[t])) +
  base_theme +
  theme(
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.35, "cm")
  ) +
  guides(colour = guide_legend(override.aes = list(size = 1.5, alpha = 1)))

# ---- Right panel: coverage ---------------------------------------------------

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
    type = "Pipeline"
  )
}
cov_df <- do.call(rbind, cov_rows)

dodge <- position_dodge(width = 50)

p_cov <- ggplot(cov_df, aes(x = t, y = cov, colour = type)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", colour = col_th) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 55,
                position = dodge, linewidth = 0.5) +
  geom_point(size = 1.5, position = dodge) +
  scale_colour_manual(values = c(Oracle = col_oracle, Pipeline = col_pipeline)) +
  labs(x = "Time", y = "95% coverage") +
  base_theme +
  theme(legend.position = "none")

row_c <- p_qq + p_cov

# ==============================================================================
# Composite figure
# ==============================================================================

# Use design layout so patchwork can align axes across rows.
# Row tags are added via labs(tag) on the first panel of each row.
row_a <- row_a + labs(tag = "a")
p_rho <- p_rho + labs(tag = "b")
p_qq  <- p_qq  + labs(tag = "c")

design <- "
AAAAAA
AAAAAA
BBBBCC
BBBBCC
DDDEEE
DDDEEE
DDDEEE
"

full_fig <- row_a + p_rho + p_v + p_qq + p_cov +
  plot_layout(design = design) +
  plot_annotation(
    theme = theme(plot.tag = element_text(size = 12, face = "bold"))
  )

ggsave(file.path(IMG_DIR, "fig-validation.png"),
         full_fig, width = 6, height = 6, dpi = 300)
})

cat("All figures written to ", IMG_DIR, "\n", sep = "")
