## collect.R — assemble per-job .rds files into combined results
##
## Run this after fetch.sh has copied the per-job .rds files locally.
## By default reads from results/raw/ (where fetch.sh writes).
##
## Usage (from the repo root):
##   Rscript dev/sims/power/tide/collect.R

library(dplyr)
library(ggplot2)

RAW_DIR <- Sys.getenv("RAW_DIR", "dev/sims/power/results/raw")
OUT_DIR <- "dev/sims/power/results"
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
  group_by(structure, d) |>
  summarise(
    S           = n(),
    detection   = mean(detected, na.rm = TRUE),
    sensitivity = mean(sensitivity, na.rm = TRUE),
    fdr         = mean(fdr, na.rm = TRUE),
    .groups     = "drop"
  )

# ---- Save --------------------------------------------------------------------

saveRDS(results, file.path(OUT_DIR, "results.rds"))
saveRDS(results_summary, file.path(OUT_DIR, "results_summary.rds"))

# ---- Plot: power curves by structure -----------------------------------------

fig <- ggplot(results_summary, aes(d, detection)) +
  geom_hline(yintercept = alpha, linetype = "dashed", colour = "grey60") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  scale_y_continuous(limits = c(0, 1),
                     labels = scales::percent_format(accuracy = 1)) +
  facet_wrap(~structure) +
  labs(x      = expression(paste(italic(d), "  (L"^2, " separation)")),
       y      = "Detection rate",
       title  = "Power — CLT test by trend structure",
       subtitle = sprintf("alpha = %.2f  |  dashed: nominal level", alpha)) +
  theme_minimal(base_size = 12)

ggsave(file.path(OUT_DIR, "power_curves.png"),
       fig, width = 8, height = 6, dpi = 150)
