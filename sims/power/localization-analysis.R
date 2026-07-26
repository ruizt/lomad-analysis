## localization-analysis.R — localization figure for the paper
##
## Question: how well do the test's rejections align with genuine local trend
## separation, and at what magnitude of separation does that alignment hold?
##
## DESIGN (threshold sweep). For each window W_t = {t - s + 1, ..., t} the
## true local separation is summarised by the windowed maximum
##       m_t = max_{u in W_t} |nu_1u - nu_2u| ,
## computed from the true simulated trends. A window counts as decoupled when
## m_t > c. Sweeping c traces a curve: at small c almost every window counts
## as decoupled, at large c almost none, and the two error rates trade off
## against each other in between. Each point on the curve is one value of c.
##
##   sens(c) = P(rejected  | m_t >  c)     "of the windows separated by at
##                                          least c, how many did we flag?"
##   spec(c) = P(!rejected | m_t <= c)     "of the windows separated by less
##                                          than c, how many did we pass?"
##
## ORIENTATION. Set ORIENT below:
##   "predictive" (default) — precision P(m_t > c | rejected) against 1 - NPV,
##                            conditioning on the test's decision. Because the
##                            sweep asks how well the rejections line up with
##                            separation of at least magnitude c, conditioning
##                            on the flags is the natural choice. This is the
##                            orientation used in postprocess_roc.qmd (where
##                            the two quantities were labelled sens/spec).
##   "conventional"         — sens/spec as defined above, conditioning on the
##                            truth instead.
##
## In the predictive orientation both axes are P(m_t > c | .) conditioned on
## rejection and on non-rejection, so sweeping c traces a genuine ROC: the
## true separation m_t is the score and rejection status is the class label.
## Its area therefore carries the standard reading — the probability that a
## randomly chosen rejected window has larger true separation than a randomly
## chosen non-rejected one — a concordance measure of localization. The area
## is not plotted here but follows directly from the saved sweep.
##
## Curves pool over separation d; they are separated by trend structure,
## AR(1) coefficient phi, SNR, and rolling-window size s (equivalently T).
##
## Usage (from the repo root):
##   Rscript sims/power/localization-analysis.R
##
## Outputs:
##   sims/power/results/_img/fig_localization.png
##   sims/power/results/localization.rds  (both orientations + profile)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggh4x)
})

RAW_DIR <- "sims/power/results/_raw"
OUT_DIR <- "sims/power/results"
IMG_DIR <- file.path(OUT_DIR, "_img")

ORIENT   <- "predictive"            # or "conventional"
MMAX     <- 1.0                     # top of the separation grid
NBIN     <- 1000L                   # separation-grid resolution
C_MAX    <- 0.30                    # plot the sweep over c in [0, C_MAX]
C_MARKS  <- c(0.01, 0.02, 0.05, 0.10, 0.20)  # c values annotated on the curve
MIN_N    <- 10000L    # drop sweep points with fewer than this many windows on
                      # either side of c: at large c the "separated" class
                      # empties out and the rates become pure noise

dir.create(IMG_DIR, showWarnings = FALSE, recursive = TRUE)

h_win_fcn <- function(n) max(5L, floor(n / 200L))
s_win_fcn <- function(n) min(60L * h_win_fcn(n), floor(n / 4L))

## Trailing rolling maximum, van Herk / Gil-Werman (O(n), verified against
## a brute-force implementation).
runmax_trailing <- function(x, s) {
  n <- length(x); if (s <= 1L) return(x)
  m  <- ceiling(n / s) * s
  xp <- c(x, rep(-Inf, m - n))
  M  <- matrix(xp, nrow = s)
  Fw <- as.vector(matrixStats::colCummaxs(M))
  Bw <- as.vector(matrixStats::colCummaxs(M[s:1, , drop = FALSE])[s:1, , drop = FALSE])
  out <- rep(NA_real_, n); idx <- s:n
  out[idx] <- pmax(Bw[idx - s + 1L], Fw[idx])
  out
}

parse_meta <- function(f) {
  bn <- sub("-series\\.rds$", "", basename(f))
  oracle <- grepl("-oracle", bn); bn <- sub("-oracle", "", bn)
  p <- regmatches(bn, regexec(
    "^([a-z]+)_d([0-9-]+)_n([0-9]+)_snr([0-9-]+)_phi([0-9-]+)$", bn))[[1]]
  if (length(p) < 6L) return(NULL)
  data.frame(file = f, struct = p[2],
             d = as.numeric(sub("-", ".", p[3])), n = as.integer(p[4]),
             snr = as.numeric(sub("-", ".", p[5])),
             phi = as.numeric(sub("-", ".", p[6])),
             oracle = oracle, stringsAsFactors = FALSE)
}

series_files <- list.files(RAW_DIR, pattern = "-series\\.rds$", full.names = TRUE)
meta <- do.call(rbind, lapply(series_files, parse_meta))
## Exclude d = 0. Those series are null everywhere, so they contribute a large
## mass of zero-separation windows that are almost never rejected; including
## them drags P(m > c | not rejected) toward 0 and flatters the curves.
meta <- meta[!meta$oracle & meta$d > 0, ]
message(sprintf("Series files (estimated noise, d > 0): %d", nrow(meta)))

## For each cell, tally windows by binned m_t, split by the test's decision.
## Every sweep quantity follows from these two vectors by cumulative sums.
cells <- list(); corrupt <- character(0)
edges <- seq(0, MMAX, length.out = NBIN + 1L)

for (i in seq_len(nrow(meta))) {
  m <- meta[i, ]
  x <- tryCatch(readRDS(m$file), error = function(e) NULL)
  if (is.null(x)) { corrupt <- c(corrupt, basename(m$file)); next }

  s_win <- s_win_fcn(m$n)
  key <- paste(m$struct, m$phi, m$snr, s_win, sep = "|")
  if (is.null(cells[[key]]))
    cells[[key]] <- list(rej = numeric(NBIN), non = numeric(NBIN))
  cl <- cells[[key]]

  for (rep in x) {
    if (is.null(rep)) next
    vi <- rep$vi; if (!length(vi)) next
    mx <- runmax_trailing(rep$sep, s_win)[vi]
    rj <- rep$rejected[vi]
    ok <- is.finite(mx) & !is.na(rj)
    if (!any(ok)) next
    mx <- mx[ok]; rj <- rj[ok]
    b  <- pmin(pmax(findInterval(mx, edges, rightmost.closed = TRUE), 1L), NBIN)
    if (any(rj))  cl$rej <- cl$rej + tabulate(b[rj],  nbins = NBIN)
    if (any(!rj)) cl$non <- cl$non + tabulate(b[!rj], nbins = NBIN)
  }
  cells[[key]] <- cl
  if (i %% 25L == 0L) message(sprintf("  ... %d / %d files", i, nrow(meta)))
}
if (length(corrupt))
  warning("Skipped unreadable file(s): ", paste(corrupt, collapse = ", "))

## ---- sweep ---------------------------------------------------------------
## At threshold c = edges[k+1]: bins 1..k are "m_t <= c", bins k+1..NBIN are
## "m_t > c".  A = rejected & above, B = rejected & below,
##             C = passed   & above, D = passed   & below.
sweep <- bind_rows(lapply(names(cells), function(key) {
  cl <- cells[[key]]; k <- strsplit(key, "|", fixed = TRUE)[[1]]
  cum_rej <- cumsum(cl$rej); cum_non <- cumsum(cl$non)
  tot_rej <- sum(cl$rej);    tot_non <- sum(cl$non)
  B <- cum_rej; D <- cum_non
  A <- tot_rej - B; C <- tot_non - D
  data.frame(
    struct = k[1], phi = as.numeric(k[2]), snr = as.numeric(k[3]),
    s_win = as.integer(k[4]),
    c = edges[-1],
    sens = A / pmax(1, A + C),          # P(rejected  | m > c)
    spec = D / pmax(1, D + B),          # P(!rejected | m <= c)
    prec = A / pmax(1, A + B),          # P(m > c | rejected)
    npv  = D / pmax(1, D + C),          # P(m <= c | !rejected)
    n_above = A + C, n_below = B + D,
    stringsAsFactors = FALSE
  )
}))

saveRDS(list(sweep = sweep, orient = ORIENT, corrupt = corrupt),
        file.path(OUT_DIR, "localization.rds"))

if (ORIENT == "conventional") {
  sweep$xx <- 1 - sweep$spec; sweep$yy <- sweep$sens
  xlab <- "1 - specificity   P(rejected | separation <= c)"
  ylab <- "Sensitivity   P(rejected | separation > c)"
} else {
  sweep$xx <- 1 - sweep$npv;  sweep$yy <- sweep$prec
  xlab <- "1 - NPV   P(separation > c | not rejected)"
  ylab <- "Precision   P(separation > c | rejected)"
}

cat(sprintf("\n== sweep at selected c (orientation: %s) ==\n", ORIENT))
print(as.data.frame(
  sweep |> filter(c %in% sapply(C_MARKS, function(z) edges[-1][which.min(abs(edges[-1] - z))])) |>
    mutate(across(c(sens, spec, prec, npv), \(z) round(z, 3))) |>
    select(struct, phi, snr, s_win, c, sens, spec, prec, npv) |>
    arrange(struct, phi, snr, s_win, c)), row.names = FALSE)

## ---- figure ---------------------------------------------------------------
lab_struct <- c(smooth = "Smooth", cross = "Cross", rate = "Rate")
pal <- c(Smooth = "#0072B2", Cross = "#D55E00", Rate = "#009E73")

sw <- sweep |>
  filter(c <= C_MAX, n_above >= MIN_N, n_below >= MIN_N) |>
  mutate(Structure = factor(lab_struct[struct], levels = names(pal)))

marks <- bind_rows(lapply(C_MARKS, function(z) {
  sw |> group_by(struct, phi, snr, s_win) |>
    slice_min(abs(c - z), n = 1, with_ties = FALSE) |>
    ungroup() |> mutate(c_lab = z)
}))

p <- ggplot(sw, aes(xx, yy, colour = Structure)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey70", linewidth = 0.3) +
  geom_path(linewidth = 0.6) +
  geom_point(data = marks, size = 1.4) +
  facet_nested(phi ~ snr + s_win, labeller = labeller(
    phi   = function(x) paste0("phi == ", x),
    snr   = function(x) paste0("SNR == ", x),
    s_win = function(x) paste0("s == ", x),
    .default = label_parsed)) +
  scale_colour_manual(values = pal) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  labs(x = xlab, y = ylab, colour = NULL,
       caption = paste0(
         "Each curve sweeps the separation threshold c defining a decoupled ",
         "window (windowed maximum |nu_1 - nu_2|).\nPoints mark c = ",
         paste(C_MARKS, collapse = ", "),
         ". Curves pool over separation d.")) +
  theme_bw(base_size = 11) +
  theme(legend.position  = "bottom",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.15, colour = "grey85"),
        strip.background = element_rect(fill = "grey95", colour = NA),
        axis.title       = element_text(size = 9),
        plot.caption     = element_text(size = 7, hjust = 0, colour = "grey30"))

ggsave(file.path(IMG_DIR, "fig_localization.png"), p,
       width = 11, height = 6.5, dpi = 400)
cat(sprintf("\nWrote %s\n", file.path(IMG_DIR, "fig_localization.png")))
