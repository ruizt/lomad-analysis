## collect.R — assemble per-job .rds files into a compiled results object
##
## Run after fetch.sh has copied per-job .rds files locally. Deliberately does
## no plotting: the composite figure is built by
## simulations/simulation-figures.R from the object written here.
##
## Usage (from the repo root):
##   Rscript simulations/validation/collect.R
##
## Override the source directory:
##   RAW_DIR=/some/other/path Rscript simulations/validation/collect.R


RAW_DIR <- Sys.getenv("RAW_DIR", "simulations/validation/results/raw")
OUT_DIR <- "simulations/validation/results"

# ---- Load all results --------------------------------------------------------

files <- list.files(RAW_DIR, pattern = "\\.rds$", full.names = TRUE)

if (length(files) == 0) {
  stop("No .rds files found in: ", RAW_DIR,
       "\nCheck that the PVC is mounted and all jobs have completed.")
}

results <- lapply(files, readRDS)
names(results) <- vapply(results, function(r) r$experiment, character(1))

cat(sprintf("Loaded %d experiment files: %s\n",
            length(results), paste(names(results), collapse = ", ")))

# ---- Save --------------------------------------------------------------------

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
saveRDS(results,
        file.path(OUT_DIR, "simulations-validation-results.rds"))

cat(sprintf("Wrote compiled results for %d experiments to %s\n",
            length(results), OUT_DIR))
