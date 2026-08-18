## realignment-results.R -- figure and tables for the realignment study
##
## Inputs
##   simulations/realignment/results/realignment-results.rds
##
## Outputs
##   simulations/_img/sfig-realignment.png    simulated series and estimated maps
##   simulations/_tbl/stbl-realignment*.csv   rejections by scenario x method x framing
##
## `stat`/`quant` stay in the compiled .rds and out of the table: NA for lomad,
## and on the multiscale-corrected scale for MSinference. `flagged`/`total`
## count windows for lomad and grid points for MSinference, so only `reject`
## compares across methods.
##
## Usage (from the repo root):
##   Rscript simulations/realignment/realignment-results.R

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})
source(here::here("figure-theme.R"))   # PT, ANNOT, fig_sizes()

IMG_DIR <- "simulations/_img"
TBL_DIR <- "simulations/_tbl"
dir.create(IMG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TBL_DIR, showWarnings = FALSE, recursive = TRUE)

res <- readRDS("simulations/realignment/results/realignment-results.rds")
n   <- res$params$n

SERIES_PAL <- c(`1` = "#0072B2", `2` = "#D55E00")
# Disjoint from SERIES_PAL: these are estimates of the map, not series.
REALIGN_PAL <- c(raw = "#CC79A7", smoothed = "#009E73")

# Panels are stacked and aligned on time, as in fig-trends: only the bottom
# one carries the x axis.
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

# ---- the observed series, on a shared scale ---------------------------------
y_lim <- range(unlist(lapply(res$scenarios, function(sc) c(sc$y1, sc$y2))))

panel_series <- function(nm, title) {
  sc <- res$scenarios[[nm]]
  bind_rows(
    tibble(t = seq_len(n), value = sc$y1, series = "1"),
    tibble(t = seq_len(n), value = sc$y2, series = "2")) |>
    ggplot(aes(t, value, colour = series)) +
    geom_line(linewidth = 0.2) +
    scale_colour_manual(values = SERIES_PAL) +
    scale_y_continuous(limits = y_lim) +
    labs(title = title, y = expression(X[it])) +
    theme_panel
}

# ---- the map, as strips between the aligned and misaligned pair -------------
# Ticks at the limits only; the strips are too short for interior breaks.
panel_coef <- function(y, ylab, digits, title = NULL) {
  lim <- range(y)
  ggplot(data.frame(t = seq_len(n), y = y), aes(t, y)) +
    geom_line(linewidth = 0.4, colour = "grey25", alpha = 0.9) +
    scale_y_continuous(breaks = lim, limits = lim,
                       labels = function(v) sprintf(paste0("%.", digits, "f"), v),
                       expand = expansion(0.08)) +
    labs(y = ylab, title = title) +
    theme_panel
}

# ---- what the correction recovers -------------------------------------------
# Error rather than level, so the true map is the zero line. The band spans
# the full range of the true b_t.
b_drift <- diff(range(res$scenarios$af$b))
p_err <- bind_rows(
  tibble(t = seq_len(n), value = res$realignment$raw$b - res$scenarios$af$b,
         source = "raw"),
  tibble(t = seq_len(n), value = res$realignment$smoothed$b - res$scenarios$af$b,
         source = "smoothed")) |>
  mutate(source = factor(source, levels = names(REALIGN_PAL))) |>
  ggplot(aes(t, value, colour = source)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -b_drift, ymax = b_drift,
           fill = "grey50", alpha = 0.18) +
  geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey80") +
  geom_line(linewidth = 0.25) +
  scale_colour_manual(values = REALIGN_PAL, name = "Realignment scale") +
  labs(x = NULL, y = expression(hat(b)[t] - b[t]),
       title = "Local realignment error") +
  theme_panel +
  theme(legend.position       = "inside",
        legend.position.inside = c(0.01, 0.02),
        legend.justification  = c(0, 0),
        legend.background     = element_rect(fill = NA, colour = NA),
        legend.key.size       = unit(0.5, "lines"),
        legend.margin         = margin(1, 3, 1, 3),
        axis.title.x          = element_text(),
        axis.text.x           = element_text(),
        axis.ticks.x          = element_line(linewidth = 0.2))

fig <- panel_series("eq", expression("Identity"~(nu[1] == nu[2]))) /
  panel_coef(res$scenarios$af$b, expression(b[t]), 2,
             title = "Affine transformation") /
  panel_coef(res$scenarios$af$a, expression(a[t]), 3) /
  panel_series("af", expression("Local affine similarity"~
                                (nu[2] == a[t] + b[t] * nu[1]))) /
  p_err +
  plot_layout(heights = c(3, 1.0, 0.85, 3, 3))

ggsave(file.path(IMG_DIR, "sfig-realignment.png"), fig,
       width = 5, height = 6, units = "in", dpi = 450)

# ---- tables -----------------------------------------------------------------
SCEN_TBL <- c(eq = "Equality", af = "Affine similarity")
FRAM_TBL <- c(global         = "Global standardization",
              `realign-raw`  = "Local realignment (raw)",
              `realign-smooth` = "Local realignment (smoothed)")

format_table <- function(tab) {
  tab |>
    mutate(
      scenario = unname(SCEN_TBL[scenario]),
      framing  = ifelse(framing %in% names(FRAM_TBL),
                        unname(FRAM_TBL[framing]), framing),
      reject   = ifelse(reject, "Yes", "No")) |>
    select(scenario, method, framing, reject, flagged, total)
}

write.csv(format_table(res$table), file.path(TBL_DIR, "stbl-realignment.csv"),
          row.names = FALSE)

# Table only for the second draw; the figure would be a duplicate layout.
ALT_SEED <- "-seed7307.2411"
alt_path <- sprintf("simulations/realignment/results/realignment-results%s.rds",
                    ALT_SEED)
if (file.exists(alt_path)) {
  write.csv(format_table(readRDS(alt_path)$table),
            file.path(TBL_DIR, sprintf("stbl-realignment%s.csv", ALT_SEED)),
            row.names = FALSE)
} else {
  message("no alternate draw at ", alt_path, "; skipping its table")
}

message("wrote figure and tables to ", IMG_DIR, " and ", TBL_DIR)
