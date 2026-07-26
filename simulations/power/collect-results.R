## collect-results.R — assemble per-job .rds files into compiled results + summary
##
## Run this after fetch.sh has copied the per-job .rds files locally. Produces
## the two compiled artifacts the figure stage consumes, so the batch pipeline
## goes raw per-job files -> compiled results -> results summary in one step.
## Deliberately does no plotting: figures are built by
## simulations/simulation-figures.R.
##
## Usage (from the repo root):
##   Rscript simulations/power/collect-results.R

library(dplyr)

RAW_DIR <- Sys.getenv("RAW_DIR", "simulations/power/results/_raw")
OUT_DIR <- "simulations/power/results"

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

# ---- Summarise ---------------------------------------------------------------

alpha <- 0.05

results_summary <- results |>
  group_by(struct, d, n, snr, phi, method) |>
  summarise(
    S         = n(),
    detection = mean(detected, na.rm = TRUE),
    ci_lo     = qbeta(0.025, sum(detected), S - sum(detected) + 1),
    ci_hi     = qbeta(0.975, sum(detected) + 1, S - sum(detected)),
    .groups   = "drop"
  )

# ---- Save --------------------------------------------------------------------

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
saveRDS(results,         file.path(OUT_DIR, "simulations-power-results.rds"))
saveRDS(results_summary, file.path(OUT_DIR, "simulations-power-summary.rds"))

cat(sprintf("Wrote %d compiled rows and %d summary rows to %s\n",
            nrow(results), nrow(results_summary), OUT_DIR))
