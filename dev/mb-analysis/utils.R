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


#' Estimate pooled ARMA noise across multiple blocks
#'
#' Takes a list of initial `lomad_fit()` results (one per block) and the
#' corresponding block data. For each block, calls `estimate_arma_noise()`
#' using the trend from the initial fit. The per-block ARMA estimates are
#' then pooled via weighted averaging of AR coefficients and innovation
#' variances (weighted by block length).
#'
#' @param fits List of `lomad_fit()` results (method = "clt").
#' @param blocks List of `list(x1, x2)` pairs, same order as `fits`.
#' @param p_max Maximum AR order for BIC selection (default 5).
#'
#' @return A list of two specs (one per series) suitable for `noise_override`.
pool_arma_noise <- function(fits, blocks, p_max = 5L) {
  n_blocks <- length(fits)
  n_vec <- vapply(fits, \(f) f$inputs$n, numeric(1))

  # Estimate ARMA noise per block
  block_noise <- mapply(function(f, b) {
    estimate_arma_noise(b$x1, b$x2, f$trend, p_max = p_max)
  }, fits, blocks, SIMPLIFY = FALSE)

  pool_one <- function(slot) {
    # Determine max AR order across blocks
    orders  <- vapply(block_noise, \(bn) bn[[slot]]$order[1], integer(1))
    max_p   <- max(orders)

    # Pad shorter AR vectors with zeros
    ar_mat <- t(sapply(block_noise, function(bn) {
      ar <- bn[[slot]]$ar
      c(ar, rep(0, max_p - length(ar)))
    }))

    sigma2s <- vapply(block_noise, \(bn) bn[[slot]]$sigma2, numeric(1))

    list(ar     = as.numeric(apply(ar_mat, 2, weighted.mean, w = n_vec)),
         sigma2 = weighted.mean(sigma2s, n_vec))
  }

  list(pool_one("series1"), pool_one("series2"))
}
