library(data.table)
library(lubridate)
library(dplyr)

#' Remove tidal periodicity from a time series by spectral notching
#'
#' Zeros the Fourier components in narrow bands around the named tidal
#' constituents, then optionally downsamples by keeping the earliest
#' observation within each floored time bin.
#'
#' This replaces an earlier implementation that used a centered rolling mean of
#' one tidal period. That is a low-pass filter: it removes the tide, but also
#' everything faster, and a boxcar's attenuation extends well below its nominal
#' cutoff. The result at 6-hourly sampling was a series band-limited far below
#' its own Nyquist frequency -- smooth by construction, with lag-2 variogram
#' more than twice the lag-1 value, which no stationary AR(1) can produce. The
#' noise model could not be fitted at all (`phi_hat` clamped on 36 of 37
#' blocks). Notching removes only the tidal bands and leaves the rest of the
#' spectrum intact, which both fits the noise model (`phi_hat` ~ 0.28, no
#' clamping) and removes the tide more completely (residual spectral line 4.4x
#' background against 8.3x for the rolling mean).
#'
#' @param df Data frame with a datetime column and the series to filter.
#' @param cols Character vector of column names to filter.
#' @param datetime_col Name of the datetime column (default `"datetime"`).
#' @param periods Numeric vector of tidal periods in hours. Defaults to the
#'   five dominant constituents: M2 (principal lunar semidiurnal), S2
#'   (principal solar semidiurnal), N2 (larger lunar elliptic semidiurnal), K1
#'   (lunisolar declinational diurnal) and O1 (principal lunar diurnal).
#' @param half_width Notch half-width in cycles per hour. The default 1/300
#'   spans the M2/S2 beat at ~355 h (the spring-neap cycle), so amplitude
#'   modulation is removed with the constituent rather than left behind.
#' @param step Passed to [lubridate::floor_date()] for downsampling. `NULL`
#'   skips downsampling.
#'
#' @return The input data frame with filtered values in `cols`, fewer rows if
#'   downsampled, otherwise the same columns.
TIDAL_PERIODS <- c(M2 = 12.4206, S2 = 12.0000, N2 = 12.6583,
                   K1 = 23.9345, O1 = 25.8193)

presmooth_tidal <- function(df,
                            cols,
                            datetime_col = "datetime",
                            periods      = TIDAL_PERIODS,
                            half_width   = 1/300,
                            step         = "6 hours") {
  notch <- function(x) {
    keep <- !is.na(x)
    if (sum(keep) < 4L) return(x)
    y  <- x[keep]; n <- length(y); mu <- mean(y)
    # frequency in cycles per sample, folded to [0, 1/2]; the series is hourly,
    # so cycles per sample == cycles per hour
    f  <- (seq_len(n) - 1) / n
    f  <- pmin(f, 1 - f)
    X  <- stats::fft(y - mu)
    X[Reduce(`|`, lapply(periods, \(p) abs(f - 1/p) <= half_width))] <- 0
    x[keep] <- Re(stats::fft(X, inverse = TRUE)) / n + mu
    x
  }

  df <- df |>
    mutate(across(.cols = all_of(cols), .fns = notch)) |>
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
