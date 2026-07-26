## collect.R — assemble per-job .rds files into combined results
##
## Run this after fetch.sh has copied the per-job .rds files locally.
## By default reads from results/raw/ (where fetch.sh writes).
##
## Usage (from the repo root):
##   Rscript sims/tide-example/tide/collect.R
## or interactively in RStudio — the defaults should work as-is.
##
## Override the source directory:
##   RAW_DIR=/some/other/path Rscript sims/tide-example/tide/collect.R

library(dplyr)
library(ggplot2)

RAW_DIR <- Sys.getenv("RAW_DIR", "sims/tide-example/results/raw")
OUT_DIR <- "sims/tide-example/results"
alpha   <- 0.05

# ---- Collect -----------------------------------------------------------------

files <- list.files(RAW_DIR, pattern = "\\.rds$", full.names = TRUE)

if (length(files) == 0) {
  stop("No .rds files found in: ", RAW_DIR,
       "\nCheck that the PVC is mounted and all jobs have completed.")
}

results <- lapply(files, function(f) {
  obj <- readRDS(f)
  obj$results
}) |> bind_rows()

# ---- Summary -----------------------------------------------------------------

results_summary <- results |>
  group_by(n) |>
  summarise(
    S        = n(),
    coverage = mean(covered),
    se       = sqrt(coverage * (1 - coverage) / S),
    .groups  = "drop"
  )

# ---- Save --------------------------------------------------------------------

saveRDS(results, file.path(OUT_DIR, "results.rds"))
saveRDS(results_summary, file.path(OUT_DIR, "results_summary.rds"))

# ---- Plot --------------------------------------------------------------------

fig <- ggplot(results_summary, aes(x = factor(n), y = coverage)) +
  geom_hline(yintercept = 1 - alpha, linetype = "dashed", colour = "grey60") +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = coverage - 1.96 * se,
                    ymax = coverage + 1.96 * se),
                width = 0.15) +
  scale_y_continuous(limits = c(0.85, 1),
                     labels = scales::percent_format(accuracy = 1)) +
  labs(x        = "Sample size (n)",
       y        = "Coverage",
       title    = "Hotelling T² confidence region — empirical coverage",
       subtitle = sprintf("nominal = %.0f%%  |  dashed: target", (1 - alpha) * 100)) +
  theme_minimal(base_size = 12)

ggsave(file.path(OUT_DIR, "coverage_plot.png"),
       fig, width = 6, height = 4, dpi = 150)
