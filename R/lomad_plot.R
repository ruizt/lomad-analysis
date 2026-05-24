#' Plot detected decoupling periods from a lomad fit
#'
#' Two-panel base R plot. The upper panel shows the MA-smoothed series
#' (`ma1`, `ma2`) and shared trend, with shading where `tst$rejected` is
#' `TRUE`. The lower panel plots the rolling correlation `R_t` with
#' theoretical `rho_t` as a dashed reference.
#'
#' @param fit List returned by [lomad_fit()].
#' @param tst List returned by [lomad_test()].
#' @param dates Optional vector of dates or axis labels.
#' @param alpha Numeric. Fill transparency for shaded regions (default 0.25).
#'
#' @return Invisibly returns `NULL`. Called for its side effect (base R plot).
#'
#' @export
lomad_plot <- function(fit, tst, dates = NULL, alpha = 0.25) {

  if (is.null(tst))
    stop("lomad_plot() requires a `tst` argument (lomad_test result).")

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
