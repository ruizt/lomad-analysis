#' Plot detected decoupling periods from a lomad fit
#'
#' Two-panel base R plot. The upper panel shows the smoothed series with shaded
#' regions marking detected decoupling periods. The lower panel shows the
#' rolling correlation series.
#'
#' For CLT fits, the upper panel plots the MA-smoothed series (`ma1`, `ma2`)
#' and shared trend, with shading where `tst$rejected` is `TRUE`. The lower
#' panel plots `R_t` with theoretical `rho_t` as a dashed reference.
#'
#' For legacy state fits, the upper panel plots raw series `x1`, `x2` with
#' the MA lines overlaid, and the lower panel shows `R_t` with the
#' BY-corrected Fisher-z threshold.
#'
#' @param fit List returned by [lomad_fit()].
#' @param tst List returned by [lomad_test()] (required for CLT fits, ignored
#'   for state fits).
#' @param x1 Numeric vector. First raw series (state fits only).
#' @param x2 Numeric vector. Second raw series (state fits only).
#' @param dates Optional vector of dates or axis labels.
#' @param alpha Numeric. Fill transparency for shaded regions (default 0.25).
#'
#' @return Invisibly returns `NULL`. Called for its side effect (base R plot).
#'
#' @export
lomad_plot <- function(fit, tst = NULL, x1 = NULL, x2 = NULL,
                       dates = NULL, alpha = 0.25) {

  method <- fit$method %||% "state"

  if (method == "clt") {
    .lomad_plot_clt(fit, tst, dates, alpha)
  } else {
    .lomad_plot_state(x1, x2, fit, dates, alpha)
  }
}


# CLT plot implementation
.lomad_plot_clt <- function(fit, tst, dates, alpha) {

  if (is.null(tst))
    stop("lomad_plot() requires a `tst` argument (lomad_test result) for CLT fits.")

  n     <- fit$inputs$n
  h     <- fit$inputs$h
  R     <- fit$R
  rho   <- fit$rho
  ma1   <- fit$ma1
  ma2   <- fit$ma2
  trend <- fit$trend

  rejected <- tst$rejected
  rejected[is.na(rejected)] <- FALSE

  t_idx <- if (is.null(dates)) seq_len(n) else dates

  # Back-shift rejected regions for upper panel
  rej_shifted <- rep(FALSE, n)
  for (t in which(rejected)) {
    s <- max(1L, t - (h - 1L))
    rej_shifted[s:t] <- TRUE
  }
  r_bs   <- rle(rej_shifted)
  ends   <- cumsum(r_bs$lengths)
  starts <- ends - r_bs$lengths + 1L
  shade_starts <- starts[r_bs$values]
  shade_ends   <- ends[r_bs$values]

  # Lower panel shading: no back-shift
  r_lo   <- rle(rejected)
  ends_lo   <- cumsum(r_lo$lengths)
  starts_lo <- ends_lo - r_lo$lengths + 1L
  shade_lo_starts <- starts_lo[r_lo$values]
  shade_lo_ends   <- ends_lo[r_lo$values]

  shade_col <- grDevices::rgb(0.7, 0.85, 1, alpha)
  col_trend <- grDevices::rgb(0.4, 0.4, 0.4, 0.8)

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))

  graphics::par(mfrow = c(2, 1), oma = c(3, 0, 0, 0))

  # Upper panel: smoothed series + trend
  graphics::par(mar = c(0, 4, 2, 1))
  ylim_top <- range(c(ma1, ma2), na.rm = TRUE)
  ylim_top <- ylim_top + c(-1, 1) * diff(ylim_top) * 0.05
  graphics::plot(t_idx, ma1, type = "n", ylim = ylim_top,
                 xlab = "", ylab = "smoothed series", xaxt = "n")
  .shade_intervals(t_idx, shade_starts, shade_ends, shade_col)
  graphics::lines(t_idx, ma1, col = "blue", lwd = 1.5)
  graphics::lines(t_idx, ma2, col = "red", lwd = 1.5)
  graphics::lines(t_idx, trend, col = col_trend, lwd = 1.2)

  # Lower panel: R_t with rho_t reference
  graphics::par(mar = c(0, 4, 0, 1))
  R_range  <- range(c(R, rho), na.rm = TRUE)
  ylim_bot <- c(min(R_range[1], -0.1) - 0.05, max(R_range[2], 0.1) + 0.05)
  graphics::plot(t_idx, R, type = "n", ylim = ylim_bot,
                 xlab = "", ylab = "correlation", xaxt = "n")
  .shade_intervals(t_idx, shade_lo_starts, shade_lo_ends, shade_col)
  graphics::lines(t_idx, R, col = "grey40")
  graphics::lines(t_idx, rho, col = "grey30", lty = 2, lwd = 1)
  graphics::abline(h = 0, col = "grey80", lty = 1, lwd = 0.5)
  graphics::axis(1)

  invisible(NULL)
}


# [LEGACY] State-process plot implementation
.lomad_plot_state <- function(x1, x2, fit, dates, alpha) {

  if (is.null(x1) || is.null(x2))
    stop("lomad_plot() requires `x1` and `x2` arguments for state fits.")

  n     <- length(x1)
  h     <- fit$null_model$h
  I     <- fit$I
  R     <- fit$R
  x_eff <- fit$thresholds$x_eff

  t_idx <- if (is.null(dates)) seq_len(n) else dates

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

  I_v       <- I[fit$valid_idx]
  m         <- length(I_v)
  entry_pos <- which(I_v == 1L & c(0L, I_v[-m]) != 1L)
  trigger_t <- pmax(1L, fit$valid_idx[entry_pos] - (h - 1L))

  col1_raw  <- grDevices::rgb(0, 0, 1, 0.2)
  col2_raw  <- grDevices::rgb(1, 0, 0, 0.2)
  col_trend <- grDevices::rgb(0.4, 0.4, 0.4, 0.8)
  shade_col <- grDevices::rgb(0.7, 0.85, 1, alpha)
  trig_col  <- grDevices::rgb(0.3, 0.3, 0.3, 0.5)

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))

  graphics::par(mfrow = c(2, 1), oma = c(3, 0, 0, 0))

  graphics::par(mar = c(0, 4, 2, 1))
  ylim_top <- c(-1.1, 1.1) * max(abs(c(x1, x2)), na.rm = TRUE)
  graphics::plot(t_idx, x1, type = "l", col = col1_raw, ylim = ylim_top,
                 xlab = "", ylab = "series", xaxt = "n")
  graphics::lines(t_idx, x2, col = col2_raw)
  .shade_intervals(t_idx, decouple_starts, decouple_ends, shade_col)
  for (tt in trigger_t)
    graphics::abline(v = t_idx[tt], col = trig_col, lty = 2, lwd = 0.8)
  graphics::lines(t_idx, x1, col = col1_raw)
  graphics::lines(t_idx, x2, col = col2_raw)
  graphics::lines(t_idx, fit$ma1, col = "blue", lwd = 1.5)
  graphics::lines(t_idx, fit$ma2, col = "red", lwd = 1.5)
  graphics::lines(t_idx, fit$trend_hat, col = col_trend, lwd = 1.2)

  graphics::par(mar = c(0, 4, 0, 1))
  R_range  <- range(R, na.rm = TRUE)
  ylim_bot <- c(min(R_range[1], x_eff, -0.1) - 0.05,
                max(R_range[2], 0.1) + 0.05)
  graphics::plot(t_idx, R, type = "l", col = "grey40", ylim = ylim_bot,
                 xlab = "", ylab = "correlation", xaxt = "n")
  below_thresh <- !is.na(I) & I == 1L
  r_bt   <- rle(below_thresh)
  ends_bt   <- cumsum(r_bt$lengths)
  starts_bt <- ends_bt - r_bt$lengths + 1L
  .shade_intervals(t_idx, starts_bt[r_bt$values], ends_bt[r_bt$values], shade_col)
  graphics::lines(t_idx, R, col = "grey40")
  graphics::abline(h = x_eff, col = "grey30", lty = 2, lwd = 1)
  graphics::abline(h = 0, col = "grey80", lty = 1, lwd = 0.5)
  graphics::axis(1)

  invisible(NULL)
}


# Shade rectangular regions between start/end index pairs.
.shade_intervals <- function(t_idx, starts, ends, col) {
  if (length(starts) == 0L) return(invisible(NULL))
  usr <- graphics::par("usr")
  for (k in seq_along(starts)) {
    graphics::rect(xleft = t_idx[starts[k]], ybottom = usr[3],
                   xright = t_idx[ends[k]], ytop = usr[4],
                   col = col, border = NA)
  }
}
