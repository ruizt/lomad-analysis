# example-figure.R -- the vignette's worked example, as a paper figure
#
# Usage (from the repo root):
#   Rscript mb-analysis/example-figure.R
#
# Output
#   mb-analysis/_img/fig-mb-example.png
#
# Same data, h, s and alpha as vignette("lomad"), so the paper and the package
# show the reader the same fit. The two panels are drawn here rather than by
# lomad_plot() only because the paper wants a colour legend and different axis
# labels, neither of which that function exposes; the layout, colours and
# shading convention are otherwise its.
#
# Note this reads the INSTALLED package's data, not anything in this repo. A
# stale install silently draws a different block, so the fit is checked against
# what the vignette reports before anything is written.

suppressPackageStartupMessages(library(lomad))

img_out <- "mb-analysis/_img"
fs::dir_create(img_out)

H <- 4L; S <- 60L; ALPHA <- 0.05

fit <- lomad_fit(morro_bay$o2, morro_bay$ph, h = H, s = S)
tst <- lomad_test(fit, alpha = ALPHA)

stopifnot(
  nrow(morro_bay) == 208L,
  as.Date(min(morro_bay$datetime)) == as.Date("2022-08-24"),
  sum(tst$rejected, na.rm = TRUE) == 52L
)

rejected <- replace(tst$rejected, is.na(tst$rejected), FALSE)
t_idx    <- morro_bay$datetime

# A rejection at t concerns the window ending at t, so the upper panel shades
# back to t - s + 1 while the lower panel shades t itself. Taken from the
# package so the two figures cannot disagree about which window a flag means.
rej_upper <- lomad:::.rejected_window_span(rejected, S)

crit <- fit$rho + qnorm(tst$alpha_eff) * sqrt(fit$V / S)

shade_runs <- function(flag) {
  r  <- rle(flag); en <- cumsum(r$lengths); st <- en - r$lengths + 1L
  cbind(st[r$values], en[r$values])
}
# Shade to the panel's own limits. An arbitrary large y range (+/-1e6) silently
# draws nothing where the panel's user span is small -- the correlation panel
# spans about 1.1 units against the upper panel's 4.3, and the same 1e6 maps to
# a device coordinate several times larger there and is dropped.
draw_shade <- function(flag, col) {
  yr <- par("usr")[3:4]
  for (i in seq_len(nrow(sp <- shade_runs(flag))))
    rect(t_idx[sp[i, 1]], yr[1], t_idx[sp[i, 2]], yr[2], col = col, border = NA)
}

# Fill transparency, not the FDR level -- lomad_plot() takes these as separate
# arguments and reusing ALPHA here shaded at 0.05 instead of 0.25.
shade_col <- rgb(0.7, 0.85, 1, 0.25)
# Match intro-figures.R. lwd is in 1/96 inch (0.2646 mm), so a width given in
# mm has to be converted; ggplot's linewidth is mm already.
MM   <- 1 / 0.2646
LWD_MA <- 0.45 * MM
col_do    <- "blue"
col_ph    <- "red"
col_trend <- rgb(0.4, 0.4, 0.4, 0.8)

png(file.path(img_out, "fig-mb-example.png"),
    width = 7, height = 4.5, units = "in", res = 450)
op <- par(no.readonly = TRUE)
par(mfrow = c(2, 1), oma = c(3, 0, 0, 0))

# ---- upper: moving averages -------------------------------------------------
par(mar = c(0, 4, 2, 1))
yl <- range(c(fit$ma1, fit$ma2), na.rm = TRUE)
plot(t_idx, fit$ma1, type = "n", ylim = yl + c(-1, 1) * diff(yl) * 0.05,
     xlab = "", ylab = "Moving averages", xaxt = "n")
draw_shade(rej_upper, shade_col)
lines(t_idx, fit$ma1,  col = col_do,    lwd = LWD_MA)
lines(t_idx, fit$ma2,  col = col_ph,    lwd = LWD_MA)
lines(t_idx, fit$trend, col = col_trend, lwd = LWD_MA * 0.8)
legend("bottomleft", legend = c("DO", "pH"), col = c(col_do, col_ph),
       lwd = 1.5, horiz = TRUE, bty = "n", cex = 0.9)

# ---- lower: correlation -----------------------------------------------------
par(mar = c(0, 4, 0, 1))
rng <- range(c(fit$R, fit$rho, crit), na.rm = TRUE)
plot(t_idx, fit$R, type = "n",
     ylim = c(min(rng[1], -0.1) - 0.05, max(rng[2], 0.1) + 0.05),
     xlab = "", ylab = "Correlation series", xaxt = "n")
draw_shade(rejected, shade_col)
ok <- which(!is.na(crit) & !is.na(fit$rho))
polygon(c(t_idx[ok], rev(t_idx[ok])), c(crit[ok], rev(fit$rho[ok])),
        col = rgb(0.55, 0.55, 0.55, 0.22), border = NA)
lines(t_idx, fit$R,   col = "grey40")
lines(t_idx, fit$rho, col = "grey30", lty = 2)
abline(h = 0, col = "grey80", lwd = 0.5)
axis.POSIXct(1, x = t_idx)

par(op); invisible(dev.off())

cat(sprintf("Wrote %s/fig-mb-example.png -- %d rows, %s to %s, %d of %d flagged\n",
            img_out, nrow(morro_bay),
            as.Date(min(t_idx)), as.Date(max(t_idx)),
            sum(rejected), length(fit$valid_idx)))
