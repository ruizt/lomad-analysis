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

# Summary files only; series files are kept on disk for downstream analysis
files <- list.files(RAW_DIR, pattern = "\\.rds$", full.names = TRUE)
files <- files[!grepl("-series\\.rds$", files)]

if (length(files) == 0) {
  stop("No .rds files found in: ", RAW_DIR,
       "\nCheck that the PVC is mounted and all jobs have completed.")
}

results <- lapply(files, function(f) {
  is_oracle <- grepl("-oracle\\.rds$", f)
  readRDS(f)$results |> mutate(method = if (is_oracle) "oracle" else "estimated")
}) |> bind_rows()

# ---- Summary -----------------------------------------------------------------

results_summary <- results |>
  group_by(struct, d, n, snr, phi, method) |>
  summarise(
    S         = n(),
    detection = mean(detected, na.rm = TRUE),
    .groups   = "drop"
  )

# ---- Save --------------------------------------------------------------------

saveRDS(results, file.path(OUT_DIR, "results.rds"))
saveRDS(results_summary, file.path(OUT_DIR, "results_summary.rds"))

# ---- Plot: power curves faceted by snr x phi, one PNG per n -----------------

for (n_val in sort(unique(results_summary$n))) {
  df <- filter(results_summary, n == n_val)

  fig <- ggplot(df, aes(d, detection, colour = struct,
                         linetype = method, group = interaction(struct, method))) +
    geom_hline(yintercept = alpha, linetype = "dashed", colour = "grey60") +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.5) +
    scale_y_continuous(limits = c(0, 1),
                       labels = scales::percent_format(accuracy = 1)) +
    scale_colour_manual(values = c("smooth" = "#0072B2",
                                   "cross"  = "#D55E00",
                                   "rate"   = "#009E73")) +
    scale_linetype_manual(values = c("estimated" = "solid",
                                     "oracle"    = "dashed")) +
    facet_grid(snr ~ phi, labeller = label_both) +
    labs(x        = expression(paste(italic(d), "  (L"^2, " separation)")),
         y        = "Detection rate",
         colour   = "Structure",
         linetype = "Method",
         title    = sprintf("Power — CLT test  |  n = %d", n_val),
         subtitle = sprintf("alpha = %.2f  |  dashed lines: oracle noise params", alpha)) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "top")

  ggsave(file.path(OUT_DIR, sprintf("power_curves_n%d.png", n_val)),
         fig, width = 8, height = 6, dpi = 150)
}
