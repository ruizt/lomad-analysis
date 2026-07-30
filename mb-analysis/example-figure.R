# example-figure.R -- the vignette's worked example, as a paper figure
#
# Usage (from the repo root):
#   Rscript mb-analysis/example-figure.R
#
# Output
#   mb-analysis/_img/fig-mb-example.png
#
# Reproduces vignette("lomad") exactly -- same data, h, s, alpha and plot call --
# so the paper and the package show the reader the same thing.
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

out <- file.path(img_out, "fig-mb-example.png")
png(out, width = 7, height = 4.5, units = "in", res = 450)
lomad_plot(fit, tst, dates = morro_bay$datetime)
invisible(dev.off())

cat(sprintf("Wrote %s -- %d rows, %s to %s, %d of %d windows flagged\n",
            out, nrow(morro_bay),
            as.Date(min(morro_bay$datetime)), as.Date(max(morro_bay$datetime)),
            sum(tst$rejected, na.rm = TRUE), length(fit$valid_idx)))
