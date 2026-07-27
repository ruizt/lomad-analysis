## localization-sweep.R — compute stage of the localization analysis
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
## Both the estimated-noise pipeline and the oracle pipeline are swept (the
## oracle series only exist at phi = 0.8, mirroring the power-curve figure),
## distinguished by a `method` column ("estimated" / "oracle").
##
## This is the expensive stage: it reads every -series.rds file in the power
## study (~1.4 GB across 360 files) and takes several minutes. It is separated
## from figure generation so that redrawing a figure never re-runs it.
##
## Usage (from the repo root):
##   Rscript simulations/power/localization-sweep.R
##
## Outputs:
##   simulations/power/results/_img/fig-localization.png
##   simulations/power/results/localization.rds  (both orientations + profile)

suppressPackageStartupMessages({
  library(dplyr)
})

RAW_DIR <- "simulations/power/results/_raw"
OUT_DIR <- "simulations/power/results"

ORIENT   <- "predictive"            # or "conventional"
MMAX     <- 1.0                     # top of the separation grid
NBIN     <- 1000L                   # separation-grid resolution
C_MAX    <- 0.30                    # plot the sweep over c in [0, C_MAX]
C_MARKS  <- c(0.01, 0.02, 0.05, 0.10, 0.20)  # c values annotated on the curve
MIN_N    <- 10000L    # drop sweep points with fewer than this many windows on
                      # either side of c: at large c the "separated" class
                      # empties out and the rates become pure noise

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

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
meta <- meta[meta$d > 0, ]
message(sprintf("Series files (estimated + oracle, d > 0): %d", nrow(meta)))

## For each cell, tally windows by binned m_t, split by the test's decision.
## Every sweep quantity follows from these two vectors by cumulative sums.
## Cells are keyed separately by method (estimated / oracle) so the two
## pipelines never mix.
cells <- list(); corrupt <- character(0)
edges <- seq(0, MMAX, length.out = NBIN + 1L)

for (i in seq_len(nrow(meta))) {
  m <- meta[i, ]
  x <- tryCatch(readRDS(m$file), error = function(e) NULL)
  if (is.null(x)) { corrupt <- c(corrupt, basename(m$file)); next }

  s_win  <- s_win_fcn(m$n)
  method <- if (m$oracle) "oracle" else "estimated"
  key <- paste(m$struct, m$phi, m$snr, s_win, method, sep = "|")
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
    s_win = as.integer(k[4]), method = k[5],
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
        file.path(OUT_DIR, "simulations-power-localization.rds"))
cat(sprintf("Wrote %s\n",
            file.path(OUT_DIR, "simulations-power-localization.rds")))
