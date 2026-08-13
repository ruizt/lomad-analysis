## collect-results.R — assemble per-cell .rds files into the compiled artifacts
##
## Run after fetch.sh has copied the per-cell files locally. Produces everything
## the figure stage consumes, so the pipeline is raw per-cell files -> compiled
## artifacts in one step. Deliberately does no plotting: figures are built by
## simulations/simulation-results.R.
##
## Three artifacts:
##   simulations-power-results.rds — one row per replicate
##   simulations-power-curves.rds  — rejection rate by delta_t bin (local power)
##   simulations-power-auc.rds     — concordance between rejection and delta_t
##
## The study reports rejection against *realized* local separation delta_t, so
## d is pooled over rather than plotted: it is a generator knob scaled per
## structure, and the same d means different separations for different
## structures. Pooling is what makes the curves comparable across them.
##
## This replaces localization-sweep.R, which had to reread every full-length
## series file to recover (delta_t, rejected) pairs. tide/sim.R now emits those
## directly, so the expensive stage is gone.
##
## Usage (from the repo root):
##   Rscript simulations/power/collect-results.R

library(dplyr)

RAW_DIR <- Sys.getenv("RAW_DIR", "simulations/power/results/_raw")
OUT_DIR <- "simulations/power/results"

# Bin edges for the local power curve. Finer at the bottom, where the curve
# turns, and where most windows sit.
DELTA_BREAKS <- c(0, 0.02, 0.05, 0.10, 0.20, 0.30, 0.45, 0.60, 0.80, 1.00)

# ---- Collect ----------------------------------------------------------------

all_files <- list.files(RAW_DIR, pattern = "\\.rds$", full.names = TRUE)
sum_files <- all_files[!grepl("-windows\\.rds$", all_files)]
win_files <- all_files[grepl("-windows\\.rds$", all_files)]

if (length(sum_files) == 0) {
  stop("No .rds files found in: ", RAW_DIR,
       "\nCheck that the PVC is mounted and all jobs have completed.")
}

tag <- function(f) if (grepl("-oracle", f)) "oracle" else "estimated"

results <- lapply(sum_files, function(f) {
  readRDS(f)$results |> mutate(method = tag(f))
}) |> bind_rows()

# Window files carry no cell metadata of their own -- they are keyed by seed,
# which the summary rows already identify -- so the cell is joined back on.
cells <- results |> distinct(seed, struct, d, s, n, snr, phi, method)

windows <- lapply(win_files, function(f) {
  readRDS(f) |> mutate(method = tag(f))
}) |> bind_rows() |> inner_join(cells, by = c("seed", "method"))

# ---- Local power curves -----------------------------------------------------

# Pooled over d and replicate. A window is the unit; d only decides which part
# of the delta_t range gets populated.
curves <- windows |>
  filter(is.finite(delta), is.finite(rejected)) |>
  mutate(bin = cut(delta, DELTA_BREAKS, include.lowest = TRUE)) |>
  group_by(struct, s, n, snr, phi, method, bin) |>
  summarise(
    windows    = n(),
    rejection  = mean(rejected),
    delta_mid  = median(delta),
    lambda_med = median(lambda, na.rm = TRUE),
    # binomial interval on the window count, which overstates precision because
    # windows overlap by s - 1; it is a floor on the uncertainty, not a CI
    se         = sqrt(rejection * (1 - rejection) / n()),
    .groups    = "drop"
  )

# ---- Concordance between rejection and true separation ----------------------

# P(a randomly chosen rejected window is more separated than a randomly chosen
# non-rejected one), by Mann-Whitney. Pooled over d for the same reason.
auc_of <- function(score, label) {
  n1 <- sum(label); n0 <- sum(!label)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  (mean(rank(score)[label]) - (n1 + 1) / 2) / n0
}

auc <- windows |>
  filter(is.finite(delta), is.finite(rejected)) |>
  group_by(struct, s, n, snr, phi, method) |>
  summarise(
    windows   = n(),
    rejection = mean(rejected),
    auc       = auc_of(delta, rejected),
    .groups   = "drop"
  )

# ---- Save -------------------------------------------------------------------

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
saveRDS(results, file.path(OUT_DIR, "simulations-power-results.rds"))
saveRDS(curves,  file.path(OUT_DIR, "simulations-power-curves.rds"))
saveRDS(auc,     file.path(OUT_DIR, "simulations-power-auc.rds"))

cat(sprintf("%d replicates, %d windows -> %d curve rows, %d cells\n",
            nrow(results), nrow(windows), nrow(curves), nrow(auc)))
