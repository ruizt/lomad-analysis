## collect.R — assemble per-job .rds files into a single results data frame
##
## Run this after fetch.sh has copied the per-job .rds files locally.
## By default reads from results/raw/ (where fetch.sh writes).
##
## Usage (from the repo root):
##   Rscript sims/power/tide/collect.R

library(dplyr)

RAW_DIR <- Sys.getenv("RAW_DIR", "sims/power/results/_raw")
OUT_DIR <- "sims/power/results"

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

# ---- Save --------------------------------------------------------------------

saveRDS(results, file.path(OUT_DIR, "results.rds"))
