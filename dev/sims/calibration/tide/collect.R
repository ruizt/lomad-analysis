## collect.R — assemble per-job .rds files into combined results
##
## Run this after fetch.sh has copied the per-job .rds files locally.
## By default reads from results/raw/ (where fetch.sh writes).
##
## Usage (from the repo root):
##   Rscript dev/sims/calibration/tide/collect.R
## or interactively in RStudio — the defaults should work as-is.
##
## Override the source directory:
##   RAW_DIR=/some/other/path Rscript dev/sims/calibration/tide/collect.R

library(dplyr)
library(ggplot2)

RAW_DIR <- Sys.getenv("RAW_DIR", "dev/sims/calibration/results/raw")
OUT_DIR <- "dev/sims/calibration/results"
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
  group_by(d) |>
  summarise(
    S             = n(),
    clt_rate      = mean(clt_rejected,      na.rm = TRUE),
    oracle_rate   = mean(oracle_rejected,   na.rm = TRUE),
    identity_rate = mean(identity_rejected, na.rm = TRUE),
    .groups = "drop"
  )

# ---- Save --------------------------------------------------------------------

saveRDS(results, file.path(OUT_DIR, "results.rds"))
saveRDS(results_summary, file.path(OUT_DIR, "results_summary.rds"))

# ---- Plot --------------------------------------------------------------------

fig <- results_summary |>
  tidyr::pivot_longer(c(clt_rate, oracle_rate, identity_rate),
                      names_to = "method", values_to = "rate") |>
  mutate(method = factor(method,
                         levels = c("clt_rate", "oracle_rate", "identity_rate"),
                         labels = c("CLT (estimated)", "CLT (oracle)",
                                    "Identity (oracle)"))) |>
  ggplot(aes(d, rate, colour = method, group = method)) +
  geom_hline(yintercept = alpha, linetype = "dashed", colour = "grey60") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  scale_y_continuous(limits = c(0, 1),
                     labels = scales::percent_format(accuracy = 1)) +
  scale_colour_manual(values = c("CLT (estimated)"  = "#0072B2",
                                 "CLT (oracle)"      = "#009E73",
                                 "Identity (oracle)" = "#D55E00")) +
  labs(x      = expression(paste(italic(d), "  (L"^2, " separation)")),
       y      = "Rejection rate",
       colour = NULL,
       title  = "Calibration — CLT test vs identity test",
       subtitle = sprintf("alpha = %.2f  |  dashed: nominal level", alpha)) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")

ggsave(file.path(OUT_DIR, "power_curve.png"),
       fig, width = 6, height = 4, dpi = 150)
