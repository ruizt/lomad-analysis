## simulation-results.R -- figures and tables for the simulation studies
##
## Reads compiled summaries only; the collect-results.R scripts do the work.
##
## Inputs
##   simulations/power/results/simulations-power-curves.rds
##   simulations/power/results/simulations-power-roc.rds
##   simulations/power/results/simulations-power-auc.rds
##   simulations/power/results/simulations-power-fdr.rds
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
##   tbl-realized-fdr.csv             realized FDR under the BY adjustment
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

# Keyed by internal code: c() mangles names when splicing named vectors.
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

# Illustrative settings, not the sweep's: the cap is 5% per window, and
# map_seed gives a monotone a_t and b_t.
cap      <- 0.05
a_mult   <- 0.9     # a_t amplitude, as a multiple of a trend's sd
map_seed <- 40

# One d for every structure, so each column meets the backdrop at w_t = 0.
struct_args <- list(
  dist = list(),
  fr   = list(rate = rate, bump = "gaussian"),
  rs   = list(bw = bw, coupling = coupling),
  rm   = list(bw = bw, coupling = coupling))

amp <- sd(sim_trends(n, d = d, method = "dist", seed = seed_coef)$x1)

set.seed(map_seed)
a_t <- a_mult * amp * lomad:::.make_affine_walk(n, 20L, bw = 0.5, cap = cap)
b_t <- 1 + lomad:::.make_affine_walk(n, 20L, bw = 0.5, cap = cap)

# Each structure twice, without and with the affine layer. `dist` has w == 0
# throughout and supplies the backdrop.
trends <- lapply(setNames(nm = names(struct_args)), function(code) {
  base <- list(n, d = d, method = code, seed = seed_coef)
  list(pre  = do.call(sim_trends, c(base, struct_args[[code]])),
       post = do.call(sim_trends, c(base, list(affine_a = a_t, affine_b = b_t),
                                    struct_args[[code]])))
})

y_lim <- range(unlist(lapply(trends, function(x)
  c(x$pre$x1, x$pre$x2, x$post$x2))))
y_pad <- diff(y_lim) * 0.2
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

# The uncoupled pair, drawn behind every mixed panel.
backdrop <- data.frame(t = seq_len(n),
                       x1 = trends[["dist"]]$pre$x1,
                       x2 = trends[["dist"]]$pre$x2) |>
  pivot_longer(c(x1, x2), names_to = "series", values_to = "value")

# x1/x2 share the structure's colour.
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
  # coloured to match the backdrop
  if (key)
    p <- p + annotate("text", x = 0.0 * n, y = y_lim[1],
                      label = "'grey = base trends'~mu[it]", parse = TRUE,
                      hjust = 0, vjust = 0, size = ANNOT, colour = "grey55")
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

# Ticks at the limits only; the strips are too short for interior breaks.
panel_coef <- function(y, ylab, limits) {
  ggplot(data.frame(t = seq_len(n), y = y), aes(t, y)) +
    geom_line(linewidth = 0.4, colour = "grey25", alpha = 0.9) +
    scale_y_continuous(breaks = limits, limits = limits, expand = expansion(0.08)) +
    labs(y = ylab) +
    theme_panel
}

struct_title <- function(code) paste0(STRUCT_FULL[[code]], " (", STRUCT_ABBR[[code]], ")")

# a_t / b_t repeat in every column. Titles ride on the w_t row so they head
# the column. nu* is mixed but not rescaled; nu is the final trend.
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

# No `widths`: patchwork sizes the null (panel) space, so equal weights
# already give equal panels.
fig_trends <- column("fr") | column("rs") | column("rm")

ggsave(file.path(IMG_DIR, "fig-trends.png"),
       fig_trends, width = 5.5, height = 3, units = "in", dpi = 450)

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

  # delta_t = 1 is an atom; drawn as a point, not the end of the curve.
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
  # guides = "collect" merges identical guides only, and panel A also maps
  # fill; suppress on B and C via guides(), which a trailing `&` cannot undo.
  drop_guide <- guides(colour = "none", fill = "none")

  composite <-
    (p_power + guides(fill = "none")) /
    (p_local + drop_guide) +
    plot_layout(guides = "collect", heights = c(1, 1)) +
    plot_annotation(tag_levels = "A") &
    theme(legend.position = "bottom")

  out <- file.path(IMG_DIR, "fig-power.png")
  ggsave(out, composite, width = 5.5, height = 6, dpi = 450)
  cat(sprintf("\nWrote %s\n", out))
})

# =============================================================================
# fig-validation.png -- finite-sample accuracy of the CLT and of the plug-in
# =============================================================================

local({
  results <- readRDS("simulations/validation/results/simulations-validation-results.rds")

# Experiments are keyed <type>-s<window>; look up by type.
exp_of <- function(type) {
  id <- grep(paste0("^", type, "-"), names(results), value = TRUE)
  if (length(id) != 1L)
    stop("expected exactly one '", type, "' experiment, found ", length(id))
  results[[id]]
}

base_theme <- theme_minimal(base_size = PT$title) +
  theme(panel.grid.minor = element_blank()) +
  fig_sizes()

col_th  <- "firebrick"
col_emp <- "grey60"
col_oracle   <- "#0072B2"
col_pipeline <- "#D55E00"

# ---- The two trends ---------------------------------------------------------
# Regenerated from the seed; b_scale is read from the results, not hard-coded.

b_scale <- exp_of("e2e")$b_scale
if (is.null(b_scale)) stop("results predate b_scale; re-run the validation study")

trend_v <- sim_trends(n = 2000, d = 0, nb = 25, sd0 = 50, p = 1.5, seed = 5381)$x1

# Shared with rho, which is undefined over the first window.
TLIM <- c(1, length(trend_v))

trend_df <- rbind(
  data.frame(t = seq_along(trend_v), nu = trend_v,          k = "1"),
  data.frame(t = seq_along(trend_v), nu = b_scale * trend_v, k = "2")
)

p_trend <- ggplot(trend_df, aes(t, nu, group = k)) +
  geom_line(linewidth = 0.3, colour = col_th) +
  coord_cartesian(xlim = TLIM) +
  base_theme +
  theme(axis.text.x = element_blank(), plot.margin = margin(5.5, 5.5, 6, 5.5),
        legend.position = "none") +
  labs(x = NULL, y = expression(nu[it]),
       title = NULL)

# ---- Proposition 1 moment accuracy ------------------------------------------

# Both moments come from the end-to-end experiment's dense grid.
mom <- exp_of("e2e")$band
mom <- mom[mom$type == "Oracle" & is.finite(mom$R_mean) & is.finite(mom$R_var), ]

rho_df <- data.frame(
  t    = rep(mom$t, 2),
  rho  = c(mom$rho, mom$R_mean),
  type = rep(c("theoretical", "empirical"), each = nrow(mom))
)

p_rho <- ggplot() +
  geom_line(data = rho_df[rho_df$type == "empirical", ],
            aes(t, rho), colour = "grey30", linewidth = 0.4, alpha = 0.8) +
  geom_line(data = rho_df[rho_df$type == "theoretical", ],
            aes(t, rho), colour = col_th, linewidth = 0.2, alpha = 0.8) +
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

v_df <- data.frame(V_theory = mom$V_th, V_emp = exp_of("e2e")$s * mom$R_var)
v_df <- v_df[is.finite(v_df$V_theory) & is.finite(v_df$V_emp) & v_df$V_theory > 0, ]
rng <- range(c(v_df$V_theory, v_df$V_emp))

p_v <- ggplot(v_df, aes(V_theory, V_emp)) +
  geom_abline(slope = 1, intercept = 0, colour = col_th, linewidth = 0.6) +
  geom_point(colour = "black", size = 1.2, alpha = 0.4) +
  coord_equal(xlim = rng, ylim = rng) +
  labs(x = expression("Theoretical" ~ V[t]),
       y = expression(s %.% Var(R[t])),
       title = NULL) +
  base_theme

# ---- End to end -------------------------------------------------------------

r_e2e    <- exp_of("e2e")
s_val    <- r_e2e$s
eval_pts <- r_e2e$eval_pts
n_pts    <- length(eval_pts)

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
  # Common symmetric range; QQ_CLIP_N reports the points this clips.
  coord_cartesian(xlim = c(-QQ_LIM, QQ_LIM), ylim = c(-QQ_LIM, QQ_LIM)) +
  labs(x = "N(0, 1) quantile", y = expression(G[t] ~ "quantile")) +
  base_theme +
  theme(aspect.ratio = 1, legend.position = "none")

# Lower tail only: lomad_test rejects on p = pnorm(Z) <= alpha.
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

# Local SNR over t, on the coverage panel's axis. Both series are drawn.
lam_df <- rbind(
  data.frame(t = seq_along(r_e2e$lambda1), lam = r_e2e$lambda1, k = "1"),
  data.frame(t = seq_along(r_e2e$lambda2), lam = r_e2e$lambda2, k = "2")
)
lam_df <- lam_df[is.finite(lam_df$lam) & lam_df$lam > 0, ]
lam_df$lam <- log10(lam_df$lam)
lam_ratio <- round(median(r_e2e$lambda2 / r_e2e$lambda1, na.rm = TRUE), 2)

p_lam <- ggplot(lam_df, aes(t, lam, group = k)) +
  geom_line(linewidth = 0.3, colour = col_th) +

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
  scale_y_continuous(breaks = c(0.05, 0.10, 0.15)) +
  coord_cartesian(xlim = TLIM, ylim = c(0.025, 0.15)) +
  labs(x = "Time", y = "Type I error") +
  base_theme +
  theme(
    # One key for both panels of the row; inset to clear the t = 1900 bars.
    legend.position = c(0.99, 0.99),
    legend.justification = c(1, 1),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.title = element_blank(),
    legend.key.size = unit(0.4, "cm")
  ) +
  guides(colour = guide_legend(override.aes = list(size = 1.2)))

# ---- Composite --------------------------------------------------------------
# Two rows; trend flush above rho on a shared time axis.

full_fig <- (p_trend + labs(tag = "A")) + p_rho + p_v +
  (p_lam + labs(tag = "B")) + p_cov + p_qq +
  plot_layout(design = c(
    area(1,  1,  3, 4),   # A  trend
    area(4,  1,  6, 4),   #    rho
    area(1,  5,  6, 6),   #    V
    # Row B is taller: the QQ facets are square and height binds first.
    area(7,  1,  9, 4),   # B  local SNR, flush above coverage
    area(10, 1, 15, 4),   #    coverage, on panel A's time axis
    area(7,  5, 15, 6)    #    end-to-end QQ, faceted 2x2 by evaluation point
  ))

ggsave(file.path(IMG_DIR, "fig-validation.png"),
         full_fig, width = 5.5, height = 4, dpi = 450)
})

# =============================================================================
# tbl-validation-params.csv, tbl-oracle-quantities.csv -- the paper's Tables 1
# and 2, from the compiled results rather than re-derived
# =============================================================================

local({
  results <- readRDS("simulations/validation/results/simulations-validation-results.rds")
  e <- results[[grep("^e2e-", names(results), value = TRUE)]]

  params <- data.frame(
    series     = 1:2,
    phi        = c(0.5, 0.3),
    sigma2     = c(0.64, 0.64),
    sigma2_eta = round(e$sigma2_eta, 3))
  write.csv(params, file.path(TBL_DIR, "tbl-validation-params.csv"), row.names = FALSE)

  p <- e$eval_pts
  oracle <- data.frame(
    s      = e$s,
    t      = p,
    tau2_1 = round(e$tau2_1[p], 4),
    tau2_2 = round(e$tau2_2[p], 4),
    lambda1 = round(e$lambda1[p], 3),
    lambda2 = round(e$lambda2[p], 3),
    rho    = round(e$rho_oracle, 3),
    V      = round(e$V_oracle, 3))
  write.csv(oracle, file.path(TBL_DIR, "tbl-oracle-quantities.csv"), row.names = FALSE)

  cat(sprintf("Wrote %s and %s\n",
              file.path(TBL_DIR, "tbl-validation-params.csv"),
              file.path(TBL_DIR, "tbl-oracle-quantities.csv")))
})

cat("All figures written to ", IMG_DIR, "\n", sep = "")

# =============================================================================
# tbl-realized-fdr.csv -- realized false discovery rate under the BY adjustment
# =============================================================================

local({
  fdr <- readRDS("simulations/power/results/simulations-power-fdr.rds")

  out <- fdr |>
    transmute(epsilon    = eps,
              null_share = round(null_share, 3),
              fdr        = round(fdr, 4),
              se         = round(se, 4))

  write.csv(out, file.path(TBL_DIR, "tbl-realized-fdr.csv"), row.names = FALSE)

  cat(sprintf("\nRealized FDR over %s datasets (BY at the simulation's alpha):\n",
              format(fdr$datasets[1], big.mark = ",")))
  print(as.data.frame(out), row.names = FALSE)
})
