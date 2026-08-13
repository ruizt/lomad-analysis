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

results <- lapply(sum_files, function(f) readRDS(f)$results) |> bind_rows()

# The sweep holds ~194 million window rows, which cannot be bound into one
# frame: reading them all at once exhausts memory well past 16 GB. Nothing
# downstream needs the rows themselves, only counts per delta_t bin, so each
# file is reduced as it is read and only the tallies are kept.
#
# NBIN is the resolution the concordance is computed at; the curve uses the
# coarser DELTA_BREAKS. 1000 bins over [0, 1] makes the tie correction below
# negligible.
NBIN <- 1000L

# Window files carry no cell metadata of their own, and cannot be joined on
# `seed`: sim.R derives its seeds from `d` alone, so every (s, snr, phi) cell
# sharing a `d` draws the *same* seeds. Each window file is therefore paired
# with the summary file of the same name, and the cell taken from there.
tally <- lapply(win_files, function(f) {
  meta <- readRDS(sub("-windows\\.rds$", ".rds", f))
  w    <- readRDS(f)
  w    <- w[is.finite(w$delta) & is.finite(w$rejected), ]

  coarse <- cut(w$delta, DELTA_BREAKS, include.lowest = TRUE)
  fine   <- pmin(NBIN, pmax(1L, ceiling(w$delta * NBIN)))

  list(
    cell = data.frame(struct = meta$struct, s = meta$s, n = meta$n,
                      snr = meta$snr, phi = meta$phi),
    coarse = data.frame(
      bin        = levels(coarse),
      windows    = as.integer(table(coarse)),
      rejected   = as.integer(tapply(w$rejected, coarse, sum, default = 0)),
      delta_sum  = as.numeric(tapply(w$delta, coarse, sum, default = 0)),
      lambda_sum = as.numeric(tapply(w$lambda, coarse,
                                     function(v) sum(v, na.rm = TRUE), default = 0))
    ),
    rej_fine = as.integer(table(factor(fine[w$rejected], levels = 1:NBIN))),
    non_fine = as.integer(table(factor(fine[!w$rejected], levels = 1:NBIN)))
  )
})

cell_key <- vapply(tally, function(x) paste(x$cell, collapse = "|"), character(1))

# ---- Local power curves -----------------------------------------------------

# Pooled over d and replicate. A window is the unit; d only decides which part
# of the delta_t range gets populated.
curves <- lapply(split(tally, cell_key), function(g) {
  cs <- lapply(g, `[[`, "coarse") |> bind_rows() |>
    group_by(bin) |>
    summarise(windows = sum(windows), rejected = sum(rejected),
              delta_sum = sum(delta_sum), lambda_sum = sum(lambda_sum),
              .groups = "drop")
  cbind(g[[1]]$cell, cs)
}) |> bind_rows() |>
  filter(windows > 0) |>
  mutate(
    rejection  = rejected / windows,
    delta_mean = delta_sum / windows,
    lambda_mean = lambda_sum / windows,
    # binomial interval on the window count, which overstates precision because
    # windows overlap by s - 1; it is a floor on the uncertainty, not a CI
    se         = sqrt(rejection * (1 - rejection) / windows)
  ) |>
  select(-delta_sum, -lambda_sum)

# ---- Concordance between rejection and true separation ----------------------

# P(a randomly chosen rejected window is more separated than a randomly chosen
# non-rejected one). Computed from the binned counts rather than by ranking
# 194 million values: with r_k rejected and m_k non-rejected in bin k,
#   U = sum_k r_k (sum_{j<k} m_j + m_k / 2),
# which is Mann-Whitney with ties inside a bin taking the usual half credit.
auc <- lapply(split(tally, cell_key), function(g) {
  r <- Reduce(`+`, lapply(g, `[[`, "rej_fine"))
  m <- Reduce(`+`, lapply(g, `[[`, "non_fine"))
  # as.numeric throughout: R * M reaches ~1e11 here, past integer range
  R <- as.numeric(sum(r)); M <- as.numeric(sum(m))
  below <- cumsum(as.numeric(m)) - m
  cbind(g[[1]]$cell,
        data.frame(windows   = R + M,
                   rejection = R / (R + M),
                   auc       = if (R == 0 || M == 0) NA_real_
                               else sum(r * (below + m / 2)) / (R * M)))
}) |> bind_rows()

# ---- Save -------------------------------------------------------------------

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
saveRDS(results, file.path(OUT_DIR, "simulations-power-results.rds"))
saveRDS(curves,  file.path(OUT_DIR, "simulations-power-curves.rds"))
saveRDS(auc,     file.path(OUT_DIR, "simulations-power-auc.rds"))

cat(sprintf("%d replicates, %s windows -> %d curve rows, %d cells\n",
            nrow(results), format(sum(auc$windows), big.mark = ","), nrow(curves), nrow(auc)))
