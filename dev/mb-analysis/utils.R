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
      dplyr::select(-time_rounded)
  }

  df
}


#' Convert a data frame with a block-ID column into a list of series pairs
#'
#' @param df Data frame.
#' @param x1_col Character. Column name for the first series.
#' @param x2_col Character. Column name for the second series.
#' @param block_col Character. Column name for the block identifier.
#'
#' @return A named list of length K, each element a list with `x1` and `x2`
#'   numeric vectors.
blocks_from_df <- function(df, x1_col, x2_col, block_col) {
  ids <- unique(df[[block_col]])
  out <- lapply(ids, function(id) {
    sub <- df[df[[block_col]] == id, ]
    list(x1 = as.numeric(sub[[x1_col]]),
         x2 = as.numeric(sub[[x2_col]]))
  })
  stats::setNames(out, as.character(ids))
}


#' Fit the legacy state-process model across multiple independent blocks
#'
#' This is the multi-block fitting pipeline formerly in the lomad package
#' as `.lomad_fit_blocks()`. It pools ARMA noise estimates across blocks
#' and computes rolling correlations + Fisher-z thresholds per block.
#'
#' @param blocks Named list of `list(x1, x2)` pairs (one per block).
#' @param q Integer. MA smoothing window.
#' @param h Integer. Rolling correlation window.
#' @param alpha Numeric. FDR level (default 0.05).
#' @param rho0 Numeric. Null correlation (default 0).
#' @param max_pq Integer. Maximum ARMA order (default 2).
#'
#' @return A list compatible with `lomad_test()` (method = "state").
fit_blocks_state <- function(blocks, q = NULL, h = NULL,
                             alpha = 0.05, rho0 = 0, max_pq = 2) {
  # Delegate to the internal package function (still exists, marked legacy)
  lomad:::.lomad_fit_blocks(blocks = blocks, q = q, h = h,
                            alpha = alpha, rho0 = rho0, max_pq = max_pq)
}
