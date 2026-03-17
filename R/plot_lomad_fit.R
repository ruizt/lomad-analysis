#' Plot detected decoupling periods from a lomad fit
#'
#' Two-panel base R plot. The upper panel overlays the raw and smoothed series
#' with shaded regions marking detected decoupling periods (back-shifted by
#' `h - 1` to align with the data window) and dashed vertical lines at each
#' trigger point (the moment the rolling correlation first crossed the
#' threshold). The lower panel shows the rolling correlation series with the
#' BY-corrected threshold and the same shaded regions.
#'
#' @param x1 Numeric vector. First raw time series.
#' @param x2 Numeric vector. Second raw time series, same length as `x1`.
#' @param fit List returned by [lomad_fit()] or [lomad()].
#' @param dates Optional vector of dates or axis labels, same length as `x1`.
#'   If `NULL` (default), plots against a 1-based integer index.
#' @param alpha Numeric. Fill transparency for shaded regions (default 0.25).
#'
#' @return Invisibly returns `NULL`. Called for its side effect (base R plot).
#'
#' @export
plot_lomad_fit <- function(x1, x2, fit, dates = NULL, alpha = 0.25) {

  n     <- length(x1)
  h     <- fit$null_model$h
  I     <- fit$I
  R     <- fit$R
  x_eff <- fit$thresholds$x_eff

  t_idx <- if (is.null(dates)) seq_len(n) else dates

  # --- Build back-shifted, merged decoupling intervals ---
  # Expand each detected point t back to [t-(h-1), t], then merge overlaps.
  I_shifted <- rep(FALSE, n)
  for (t in which(!is.na(I) & I == 1L)) {
    s <- max(1L, t - (h - 1L))
    I_shifted[s:t] <- TRUE
  }
  r_bs   <- rle(I_shifted)
  ends   <- cumsum(r_bs$lengths)
  starts <- ends - r_bs$lengths + 1L
  decouple_starts <- starts[r_bs$values]
  decouple_ends   <- ends[r_bs$values]

  # --- Trigger points: start of the data window that triggered each episode ---
  # R[t] is right-aligned over [t-h+1, t], so the first data point involved
  # is at t_start - (h-1), i.e. the left edge of the back-shifted region.
  I_v       <- I[fit$valid_idx]
  m         <- length(I_v)
  entry_pos <- which(I_v == 1L & c(0L, I_v[-m]) != 1L)
  trigger_t <- pmax(1L, fit$valid_idx[entry_pos] - (h - 1L))

  # --- Colors ---
  col1_raw  <- grDevices::rgb(0, 0, 1, 0.2)
  col2_raw  <- grDevices::rgb(1, 0, 0, 0.2)
  col_trend <- grDevices::rgb(0.4, 0.4, 0.4, 0.8)
  shade_col <- grDevices::rgb(0.7, 0.85, 1, alpha)
  trig_col  <- grDevices::rgb(0.3, 0.3, 0.3, 0.5)

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))

  graphics::par(mfrow = c(2, 1), oma = c(3, 0, 0, 0))

  # ---- Upper panel: raw + smoothed series ----
  graphics::par(mar = c(0, 4, 2, 1))

  ylim_top <- c(-1.1, 1.1) * max(abs(c(x1, x2)), na.rm = TRUE)

  graphics::plot(t_idx, x1,
                 type = "l", col = col1_raw, ylim = ylim_top,
                 xlab = "", ylab = "series", xaxt = "n")
  graphics::lines(t_idx, x2, col = col2_raw)

  .shade_intervals(t_idx, decouple_starts, decouple_ends, shade_col)

  for (tt in trigger_t)
    graphics::abline(v = t_idx[tt], col = trig_col, lty = 2, lwd = 0.8)

  # Redraw series on top of shading
  graphics::lines(t_idx, x1,        col = col1_raw)
  graphics::lines(t_idx, x2,        col = col2_raw)
  graphics::lines(t_idx, fit$ma1,   col = "blue", lwd = 1.5)
  graphics::lines(t_idx, fit$ma2,   col = "red",  lwd = 1.5)
  graphics::lines(t_idx, fit$trend_hat, col = col_trend, lwd = 1.2)

  # ---- Lower panel: rolling correlation ----
  graphics::par(mar = c(0, 4, 0, 1))

  R_range  <- range(R, na.rm = TRUE)
  ylim_bot <- c(min(R_range[1], x_eff, -0.1) - 0.05,
                max(R_range[2],  0.1)          + 0.05)

  graphics::plot(t_idx, R,
                 type = "l", col = "grey40", ylim = ylim_bot,
                 xlab = "", ylab = "correlation", xaxt = "n")

  # Shade periods where R is actually below the threshold (no back-shift)
  below_thresh <- !is.na(I) & I == 1L
  r_bt   <- rle(below_thresh)
  ends_bt   <- cumsum(r_bt$lengths)
  starts_bt <- ends_bt - r_bt$lengths + 1L
  .shade_intervals(t_idx, starts_bt[r_bt$values], ends_bt[r_bt$values], shade_col)

  graphics::lines(t_idx, R, col = "grey40")
  graphics::abline(h = x_eff, col = "grey30", lty = 2, lwd = 1)
  graphics::abline(h = 0,     col = "grey80", lty = 1, lwd = 0.5)

  graphics::axis(1)

  invisible(NULL)
}

# Internal helper: shade rectangular regions between start/end index pairs.
.shade_intervals <- function(t_idx, starts, ends, col) {
  if (length(starts) == 0L) return(invisible(NULL))
  usr <- graphics::par("usr")
  for (k in seq_along(starts)) {
    graphics::rect(
      xleft   = t_idx[starts[k]],
      ybottom = usr[3],
      xright  = t_idx[ends[k]],
      ytop    = usr[4],
      col     = col,
      border  = NA
    )
  }
}
