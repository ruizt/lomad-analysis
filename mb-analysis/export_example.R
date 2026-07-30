# export_example.R — export the vignette example subset for the R package
#
# Usage (from the repo root):
#   Rscript mb-analysis/export_example.R
#
# Input
#   _mb-data/ph_o2_blocks.csv     (produced by mb-analysis/process_blocks.R)
#
# Output
#   mb-analysis/_export/morro_bay.rda
#
# Copy that file into lomad-package/data/ by hand. Neither repository reads
# from the other: this script writes locally, the package ships the .rda it was
# given. Re-run only when the example blocks should change.

suppressPackageStartupMessages({
  library(tidyverse)
  library(lomad)
})
source("mb-analysis/utils.R")   # presmooth_tidal()

# The block is identified by the period it covers, not by block_id: ids are
# assigned by position within a station and shift whenever the record is
# extended or the QA changes. The previous choice (Bay Mouth, spring 2023) is
# gone for exactly that reason -- correcting morro_param_map() applied 23 pH
# flag windows that had been silently skipped, one of which splits that block
# into three pieces too short to fit.
#
# One Bay Mouth block, late summer 2022, containing a single sustained
# decoupling episode of about two weeks. Chosen so that both series fit without
# hitting a variogram boundary: several candidates detect just as clearly but
# clamp phi_hat at its lower bound, and example data should not greet a user
# with a warning that the noise model may be misspecified.
#
# Bay Mouth rather than Bay Head because Bay Head records detections in only
# two of its analysable blocks, both very long; a short Bay Head example would
# show the method finding nothing.
KEEP_LOCATION <- "BM1"
PERIODS   <- list(c("2022-08-20", "2022-10-20"))
MIN_HOURS <- 1000L   # excludes short fragments sharing a window

blocks_csv <- "_mb-data/ph_o2_blocks.csv"
if (!file.exists(blocks_csv))
  stop(blocks_csv, " not found; run mb-analysis/process_blocks.R first.")

all_blocks <- read_csv(blocks_csv, show_col_types = FALSE) |>
  filter(location == KEEP_LOCATION)

pick <- vapply(PERIODS, function(w) {
  hit <- all_blocks |>
    group_by(block_id) |>
    summarise(n = n(), start = min(datetime), end = max(datetime),
              .groups = "drop") |>
    filter(start >= as.POSIXct(w[1], tz = "UTC"),
           end   <= as.POSIXct(w[2], tz = "UTC"),
           n >= MIN_HOURS)
  if (nrow(hit) != 1L)
    stop("expected exactly one ", KEEP_LOCATION, " block within ",
         w[1], " .. ", w[2], "; found ", nrow(hit))
  hit$block_id
}, numeric(1))

morro_bay <- all_blocks |>
  filter(block_id %in% pick) |>
  group_split(block_id) |>
  lapply(presmooth_tidal, cols = c("o2", "ph"), step = "6 hours") |>
  bind_rows() |>
  transmute(datetime = datetime, o2 = o2, ph = ph) |>
  arrange(datetime) |>
  as.data.frame()
rownames(morro_bay) <- NULL

dir.create("mb-analysis/_export", showWarnings = FALSE, recursive = TRUE)
out <- "mb-analysis/_export/morro_bay.rda"
save(morro_bay, file = out, compress = "xz")

cat(sprintf("Wrote %s: %d rows, %s .. %s (%.0f KB)\n", out, nrow(morro_bay),
            as.Date(min(morro_bay$datetime)), as.Date(max(morro_bay$datetime)),
            file.size(out) / 1024))

# The exported data must be something the noise model can actually fit -- the
# previous example shipped with a lag-2/lag-1 variogram ratio above 2, which no
# stationary AR(1) can produce, so lomad_fit() clamped phi_hat on every fit.
vg <- function(x, l) mean((x[(l + 1):length(x)] - x[1:(length(x) - l)])^2) / 2
f <- suppressWarnings(suppressMessages(
  lomad::lomad_fit(morro_bay$o2, morro_bay$ph, h = 4, s = 60)))
vt <- which(!is.na(f$trend)); r <- morro_bay$o2[vt] - f$trend[vt]
tst <- suppressMessages(lomad::lomad_test(f, alpha = 0.05))
cat(sprintf("  V2/V1 = %.2f, phi_hat = (%.3f, %.3f), rejections = %d\n",
            vg(r, 2) / vg(r, 1), f$noise$series1$ar, f$noise$series2$ar,
            sum(tst$rejected, na.rm = TRUE)))
cat("  (V2/V1 must be below 2 for a stationary AR(1) to fit)\n")
cat("\nNow copy it across:\n")
cat("  cp", out, "../lomad-package/data/morro_bay.rda\n")
