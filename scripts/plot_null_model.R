plot_null_model <- function(x1, x2, fit, dates,
                            alpha = 0.3,
                            labels = list(x1 = "x1 (smoothed)",
                                                 x2 = "x2 (smoothed)"),
                            legend_pos = "topleft") {
  
  stopifnot(length(x1) == length(x2),
            length(x1) == length(dates))
  
  I     <- fit$I
  h     <- fit$null_model$h_corr
  dates <- as.POSIXct(dates)
  shade_col = rgb(0.8, 0.8, 0.8, alpha)
  
  # fallback defaults if user omits one label
  if (is.null(labels$x1)) labels$x1 <- "x1 (smoothed)"
  if (is.null(labels$x2)) labels$x2 <- "x2 (smoothed)"
  
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  par(mfrow = c(1, 1), mar = c(4, 4, 2, 1))
  
  ylim_top <- c(-1.1, 1.1) * max(abs(c(x1, x2)), na.rm = TRUE)
  
  ## ---------------------------
  ## Main plot
  ## ---------------------------
  plot(dates, x1,
       type = "l",
       col  = rgb(0, 0, 1, 0.2),
       ylim = ylim_top,
       xlab = "",
       ylab = "")
  
  ## --- Shade regions where I_t = 1, back-shifted by h ---
  r      <- rle(I == 1)
  ends   <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1
  
  for (k in which(r$values)) {
    s_idx_raw <- starts[k] - (h - 1L)
    s_idx     <- max(s_idx_raw, 1L)
    e_idx     <- min(ends[k], length(dates))
    
    rect(xleft   = dates[s_idx],
         ybottom = par("usr")[3],
         xright  = dates[e_idx],
         ytop    = par("usr")[4],
         col     = shade_col,
         border  = NA)
  }
  
  ## --- Redraw raw + smoothed series ---
  lines(dates, x1, type = "l", col = rgb(0, 0, 1, 0.2))
  lines(dates, x2, type = "l", col = rgb(1, 0, 0, 0.2))
  
  ## Smoothed versions
  lines(dates, fit$ma1,       type = "l", col = "red")
  lines(dates, fit$ma2,       type = "l", col = "blue")
  
  ## Shared trend
  lines(dates, fit$trend_hat, type = "l", col = "black")
  
  ## ---------------------------
  ## Legend (red + blue only)
  ## ---------------------------
  legend(
    legend_pos,
    legend = c(labels$x1, labels$x2),
    col    = c("red", "blue"),
    lty    = 1,
    lwd    = 2,
    bg     = rgb(1, 1, 1, 0.6)
  )
}