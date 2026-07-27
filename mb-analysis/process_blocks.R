# process_blocks.R — build the analysis blocks from the QA'd sensor record
#
# Usage (from the repo root):
#   Rscript mb-analysis/process_blocks.R
#
# Input
#   ../mb-qartod/_data/clean/wp_data.parquet  (copied into _mb-data/ on the way)
#
# Output
#   _mb-data/wp_data.parquet      the QA'd record, verbatim
#   _mb-data/ph_o2_blocks.csv     hourly, blocked, standardized; the analysis input
#
# This step used to live in mb-qartod as inst/scripts/run_block_pipeline.R. It
# belongs here: every choice in it — the bin width, the gap threshold, the
# minimum block length, which variables are standardized — is an analysis
# decision, not a quality-control one. Leaving it upstream also meant the
# pipeline's `select(time_utc, o2, ph, location)` silently discarded the tide
# and pressure channels, which we later needed to evaluate tidal filtering.
#
# The block-construction logic reproduces mb-qartod's morro_make_blocks()
# exactly, with one addition: tide and pressure are carried through. They are
# deliberately NOT part of `block_vars`, so their missingness (tide is ~50%
# complete) cannot influence where blocks begin and end.

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
})
# arrow masks lubridate::duration() with an Arrow type constructor, so it is
# used qualified rather than attached.
requireNamespace("arrow", quietly = TRUE) ||
  stop("package 'arrow' is required to read wp_data.parquet")

# ---- Configuration ----------------------------------------------------------

block_vars    <- c("o2", "ph")   # variables that jointly determine blocks
carry_vars    <- c("tide", "press")
max_lag_hours <- 24              # gap threshold
min_days      <- 5               # minimum block length
bin_size      <- "1 hour"

mb_clean <- fs::path(here::here(), "../mb-qartod/_data/clean")
out_dir  <- here::here("_mb-data")
fs::dir_create(out_dir)

# ---- Import -----------------------------------------------------------------

src <- fs::path(mb_clean, "wp_data.parquet")
if (!fs::file_exists(src))
  stop("wp_data.parquet not found at: ", src,
       "\nCheck that the mb-qartod repo is a sibling of this one and that ",
       "inst/scripts/run_qartod_pipeline.R has been run.")

dst <- fs::path(out_dir, "wp_data.parquet")
fs::file_copy(src, dst, overwrite = TRUE)
message("Copied wp_data.parquet to ", out_dir)

wp <- arrow::read_parquet(dst, col_select = c(time_utc, file_name, oxygen_mg_l,
                                       p_h_internal, tide_m_mllw,
                                       pressure_dbar)) |>
  mutate(location = substr(file_name, 9, 11),
         o2       = oxygen_mg_l,
         ph       = p_h_internal,
         tide     = tide_m_mllw,
         press    = pressure_dbar) |>
  select(time_utc, location, all_of(c(block_vars, carry_vars)))

# ---- Bin to hourly means ----------------------------------------------------

hourly <- wp |>
  mutate(bin = floor_date(time_utc, unit = bin_size)) |>
  group_by(location, bin) |>
  summarise(across(all_of(c(block_vars, carry_vars)),
                   \(x) mean(x, na.rm = TRUE)), .groups = "drop") |>
  rename(datetime = bin) |>
  # mean(NA, na.rm = TRUE) is NaN; restore it to NA so drop_na() sees it
  mutate(across(all_of(c(block_vars, carry_vars)),
                \(x) ifelse(is.nan(x), NA_real_, x))) |>
  drop_na(all_of(block_vars))

# ---- Standardize ------------------------------------------------------------
# Global, matching the upstream pipeline: one location/scale per variable over
# the whole record rather than per block or per station. Blocks are compared
# with each other downstream, so a per-block standardization would remove
# exactly the between-block differences the analysis looks at.

hourly <- hourly |>
  mutate(across(all_of(block_vars), \(x) as.numeric(scale(x))))

# ---- Identify and filter blocks ---------------------------------------------
# A new block starts at the first observation of a location or after any gap
# longer than max_lag_hours. Missing rows were already dropped above, so
# missingness enters only through the gaps it creates.

blocks <- hourly |>
  arrange(location, datetime) |>
  group_by(location) |>
  mutate(
    lag_hours = as.numeric(difftime(datetime, lag(datetime), units = "hours")),
    new_block = is.na(lag_hours) | lag_hours > max_lag_hours,
    block_id  = cumsum(new_block)
  ) |>
  ungroup() |>
  select(-lag_hours, -new_block)

min_rows <- min_days * 24 / as.numeric(lubridate::duration(bin_size), "hours")
keep <- blocks |>
  count(location, block_id) |>
  filter(n >= min_rows) |>
  select(location, block_id)

blocks <- semi_join(blocks, keep, by = c("location", "block_id")) |>
  select(datetime, location, all_of(c(block_vars, carry_vars)), block_id)

# ---- Report and write -------------------------------------------------------

summary_tbl <- blocks |>
  group_by(location, block_id) |>
  summarise(start = min(as.Date(datetime)), end = max(as.Date(datetime)),
            n_obs = n(), tide_pct = round(100 * mean(!is.na(tide))),
            .groups = "drop") |>
  mutate(days = as.numeric(end - start))

message(sprintf("%d blocks over %d rows (%d dropped below %d days)",
                nrow(summary_tbl), nrow(blocks),
                n_distinct(paste(blocks$location, blocks$block_id)) -
                  nrow(summary_tbl) + 0L, min_days))
print(summary_tbl, n = Inf)

out_path <- fs::path(out_dir, "ph_o2_blocks.csv")
write_csv(blocks, out_path)
message("Wrote ", out_path)
