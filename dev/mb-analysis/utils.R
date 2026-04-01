library(data.table)
library(lubridate)
library(dplyr)

#' Presmooth a time series to remove tidal (or other periodic) fluctuations
#'
#' Applies a centered rolling mean of width `cycle_len` to each column in
#' `cols`, drops rows with any resulting NAs, then optionally downsamples by
#' keeping the earliest observation within each floored time bin.
#'
#' @param df Data frame containing at least a datetime column and the series
#'   columns to smooth.
#' @param cols Character vector of column names to smooth.
#' @param datetime_col Name of the datetime column (default `"datetime"`).
#' @param cycle_len Integer. Window width in samples equal to one period of the
#'   fluctuation to remove (e.g. 25 for hourly data with a ~25 h tidal cycle).
#' @param step Character string passed to [lubridate::floor_date()] for
#'   downsampling (e.g. `"6 hours"`). Pass `NULL` to skip downsampling.
#'
#' @return The input data frame with smoothed values in `cols`, fewer rows (NAs
#'   dropped and downsampled), and otherwise the same columns.
presmooth_tidal <- function(df,
                            cols,
                            datetime_col = "datetime",
                            cycle_len,
                            step = "6 hours") {
  df <- df |>
    mutate(across(
      .cols = all_of(cols),
      .fns  = ~data.table::frollmean(.x, n = cycle_len, fill = NA, align = "center")
    )) |>
    drop_na(all_of(cols))

  if (!is.null(step)) {
    df <- df |>
      group_by(time_rounded = floor_date(.data[[datetime_col]], unit = step)) |>
      slice_min(.data[[datetime_col]], n = 1, with_ties = FALSE) |>
      ungroup() |>
      select(-time_rounded)
  }

  df
}
