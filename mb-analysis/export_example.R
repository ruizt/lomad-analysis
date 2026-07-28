# export_example.R — export the vignette example subset for the R package
#
# Usage (from the repo root):
#   Rscript mb-analysis/export_example.R
#
# Input
#   _mb-data/ph_o2_blocks.csv     (produced by mb-analysis/process_blocks.R)
#
# Output
#   ../lomad-package/data-raw/mb-example.rds
#
# This is a one-way hand-off, not a build dependency. The .rds is committed to
# the package repository, and the package rebuilds `morro_bay` from its own
# local copy -- it never reads anything from this repo. Re-run this script only
# when the example blocks should change.

suppressPackageStartupMessages({
  library(tidyverse)
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

example <- all_blocks |>
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
rownames(example) <- NULL

out <- "../lomad-package/data-raw/mb-example.rds"
if (!dir.exists(dirname(out)))
  stop("lomad-package/data-raw not found at ", dirname(out),
       "\nExpected lomad-package as a sibling of this repo.")

saveRDS(example, out, compress = "xz")
cat(sprintf("Wrote %s: %d rows, blocks %s (%.0f KB)\n", out, nrow(example),
            paste(unique(example$block), collapse = ", "),
            file.size(out) / 1024))

for (b in unique(example$block)) {
  d <- example[example$block == b, ]
  cat(sprintf("  block %2d: n = %3d, %s .. %s\n", b, nrow(d),
              as.Date(min(d$datetime)), as.Date(max(d$datetime))))
}
