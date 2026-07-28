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

# Blocks are identified by the period they cover, not by block_id: ids are
# assigned by position within a station and shift whenever the record is
# extended or the QA changes.
#
# Two consecutive Bay Mouth blocks: winter 2021-22, where the test flags
# nothing, and the spring 2022 block running directly on from it, which holds
# one sustained episode. Bay Mouth rather than Bay South because under the
# current pipeline Bay South records no detections in any of its analysable
# blocks, so a Bay South example would show the method finding nothing.
KEEP_LOCATION <- "BM1"
PERIODS   <- list(c("2021-10-01", "2022-01-31"),
                  c("2022-02-01", "2022-04-30"))
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
  transmute(block    = as.integer(block_id),
            datetime = datetime,
            o2       = o2,
            ph       = ph) |>
  arrange(block, datetime) |>
  as.data.frame()
rownames(morro_bay) <- NULL

dir.create("mb-analysis/_export", showWarnings = FALSE, recursive = TRUE)
out <- "mb-analysis/_export/morro_bay.rda"
save(morro_bay, file = out, compress = "xz")

cat(sprintf("Wrote %s: %d rows, blocks %s (%.0f KB)\n", out, nrow(morro_bay),
            paste(unique(morro_bay$block), collapse = ", "),
            file.size(out) / 1024))

# The exported data must be something the noise model can actually fit -- the
# previous example shipped with a lag-2/lag-1 variogram ratio above 2, which no
# stationary AR(1) can produce, so lomad_fit() clamped phi_hat on every fit.
vg <- function(x, l) mean((x[(l + 1):length(x)] - x[1:(length(x) - l)])^2) / 2
for (b in unique(morro_bay$block)) {
  d <- morro_bay[morro_bay$block == b, ]
  f <- suppressWarnings(suppressMessages(
    lomad::lomad_fit(d$o2, d$ph, h = 4, s = 60)))
  vt <- which(!is.na(f$trend)); r <- d$o2[vt] - f$trend[vt]
  cat(sprintf("  block %2d: n = %3d, %s .. %s, V2/V1 = %.2f, phi_hat = %.3f\n",
              b, nrow(d), as.Date(min(d$datetime)), as.Date(max(d$datetime)),
              vg(r, 2) / vg(r, 1), f$noise$series1$ar))
}
cat("  (V2/V1 must be below 2 for a stationary AR(1) to fit)\n")
cat("\nNow copy it across:\n")
cat("  cp", out, "../lomad-package/data/morro_bay.rda\n")
