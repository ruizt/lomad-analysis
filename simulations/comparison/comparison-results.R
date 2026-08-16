## comparison-results.R -- figure and table for the comparison study
##
## Inputs
##   simulations/comparison/results/comparison-results.rds
##
## Outputs
##   simulations/_img/fig-comparison.png    the three simulated series
##   simulations/_tbl/tbl-comparison.csv    rejections by scenario x method x framing
##
## The table drops `stat`/`quant`: they are NA for lomad, and for MSinference
## they sit on the multiscale-corrected scale (see design.md), which invites
## cross-method comparison that is not meaningful. Both remain in the compiled
## .rds. `flagged`/`total` count windows for lomad and (location, bandwidth)
## grid points for MSinference; the denominators are not comparable across
## methods, and only `reject` is.
##
## Usage (from the repo root):
##   Rscript simulations/comparison/comparison-results.R

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
})
source(here::here("figure-theme.R"))   # PT, ANNOT, fig_sizes()

IMG_DIR <- "simulations/_img"
TBL_DIR <- "simulations/_tbl"
dir.create(IMG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TBL_DIR, showWarnings = FALSE, recursive = TRUE)

res <- readRDS("simulations/comparison/results/comparison-results.rds")
n   <- res$params$n

SERIES_PAL <- c(`1` = "#0072B2", `2` = "#D55E00")
SCEN_LAB <- c(
  eq = "Equality:~nu[2]==nu[1]",
  af = "Affine~similarity:~nu[2]==a[t]+b[t]*nu[1]")

# ---- fig-comparison ---------------------------------------------------------
series_df <- function(field) {
  bind_rows(lapply(names(res$scenarios), function(nm) {
    sc <- res$scenarios[[nm]]
    bind_rows(
      tibble(t = seq_len(n), value = sc[[paste0(field, "1")]], series = "1"),
      tibble(t = seq_len(n), value = sc[[paste0(field, "2")]], series = "2")) |>
      mutate(scenario = factor(nm, levels = names(SCEN_LAB)))
  }))
}

fig <- ggplot(series_df("y"), aes(t, value, color = series)) +
  geom_line(linewidth = 0.2, alpha = 0.3) +
  geom_line(data = series_df("x"), linewidth = 0.5) +
  facet_wrap(~scenario, ncol = 1, scales = "free_y",
             labeller = as_labeller(SCEN_LAB, label_parsed)) +
  scale_color_manual(values = SERIES_PAL, name = "Series") +
  labs(x = "Time", y = "Value") +
  fig_theme()

ggsave(file.path(IMG_DIR, "fig-comparison.png"), fig,
       width = 6.5, height = 5, units = "in", dpi = 300)

# ---- tbl-comparison ---------------------------------------------------------
res$table |>
  select(scenario, method, framing, reject, flagged, total) |>
  write.csv(file.path(TBL_DIR, "tbl-comparison.csv"), row.names = FALSE)

message("wrote ", file.path(IMG_DIR, "fig-comparison.png"),
        " and ", file.path(TBL_DIR, "tbl-comparison.csv"))
