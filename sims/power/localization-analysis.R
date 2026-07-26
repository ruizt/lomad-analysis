## localization-analysis.R — localization figure for the paper
##
## Question: how much *local* trend separation does the pointwise test need
## before it flags a window, and does it stay calibrated where the trends are
## genuinely identical?
##
## TRUTH.  At each time t the local RMS separation over the test's own window
## W_t = {t - s + 1, ..., t} is
##       delta_t = sqrt( mean_{u in W_t} |nu_1u - nu_2u|^2 ),
## the paper's local RMS separation. This is computed from the true simulated
## trends, so it is exact.
##
## PRIMARY OUTPUT: the detection profile — empirical rejection rate as a
## function of delta_t. It is cutoff-free and reads directly:
##   * the delta_t = 0 bin (trends exactly identical) gives the pointwise
##     false positive rate, which should sit at or below alpha;
##   * the rise across delta_t > 0 gives the resolution of the procedure, i.e.
##     how small a local separation it reliably detects.
##
## WHY NOT AN ROC.  An ROC needs a binary truth, which forces a cutoff c on
## delta_t. Because separation is continuous and mostly small, any moderate c
## labels genuinely separated windows as negatives and understates the method
## badly (at c = 0.05 the apparent AUC is ~0.74, while the profile below shows
## near-perfect detection by delta_t = 0.03). The ROC is still computed and
## saved for reference, but the profile is the honest summary.
##
## Profiles pool over separation d and series length T; they are separated by
## trend structure, AR(1) coefficient phi, and SNR.
##
## Memory note: ~40M pointwise records, so results are accumulated as binned
## counts rather than materialized as rows.
##
## Usage (from the repo root):
##   Rscript sims/power/localization-analysis.R
##
## Outputs:
##   sims/power/results/_img/fig_localization.png
##   sims/power/results/localization.rds   (profile, ROC, AUC, operating points)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

RAW_DIR <- "sims/power/results/_raw"
OUT_DIR <- "sims/power/results"
IMG_DIR <- file.path(OUT_DIR, "_img")

ALPHA  <- 0.05            # nominal FDR level used in the simulations
C_ROC  <- 0.05            # truth cutoff, reference ROC only
NBINS  <- 500L            # p-value grid for the reference ROC

## delta_t bins: an exact-zero class, then increasing separation
DBRK <- c(-1e-9, 1e-9, 0.005, 0.01, 0.02, 0.03, 0.05, 0.08, 0.12, Inf)
DLAB <- c("0", "(0,.005]", "(.005,.01]", "(.01,.02]", "(.02,.03]",
          "(.03,.05]", "(.05,.08]", "(.08,.12]", ">.12")
DMID <- c(0, 0.0025, 0.0075, 0.015, 0.025, 0.04, 0.065, 0.10, 0.15)
nb   <- length(DLAB)

dir.create(IMG_DIR, showWarnings = FALSE, recursive = TRUE)

h_win_fcn <- function(n) max(5L, floor(n / 200L))
s_win_fcn <- function(n) min(60L * h_win_fcn(n), floor(n / 4L))

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
meta <- meta[!meta$oracle, ]
message(sprintf("Series files (estimated noise): %d", nrow(meta)))

new_cell <- function() list(
  n_tot = numeric(nb), n_rej = numeric(nb),          # detection profile
  dec = numeric(NBINS), sim = numeric(NBINS)         # reference ROC
)
cells <- list(); corrupt <- character(0)

for (i in seq_len(nrow(meta))) {
  m <- meta[i, ]
  x <- tryCatch(readRDS(m$file), error = function(e) NULL)
  if (is.null(x)) { corrupt <- c(corrupt, basename(m$file)); next }

  s_win <- s_win_fcn(m$n); kern <- rep(1 / s_win, s_win)
  key <- paste(m$struct, m$phi, m$snr, sep = "|")
  if (is.null(cells[[key]])) cells[[key]] <- new_cell()
  cl <- cells[[key]]

  for (rep in x) {
    if (is.null(rep)) next
    vi <- rep$vi; if (!length(vi)) next
    rms <- sqrt(pmax(0, as.numeric(
      stats::filter(rep$sep^2, kern, sides = 1))))[vi]
    rj <- rep$rejected[vi]; p <- rep$p_raw[vi]
    ok <- is.finite(rms) & !is.na(rj)
    if (!any(ok)) next
    rms <- rms[ok]; rj <- rj[ok]; p <- p[ok]

    db <- as.integer(cut(rms, DBRK, labels = FALSE))
    cl$n_tot <- cl$n_tot + tabulate(db, nbins = nb)
    if (any(rj)) cl$n_rej <- cl$n_rej + tabulate(db[rj], nbins = nb)

    fin <- is.finite(p)
    if (any(fin)) {
      pb <- pmin(floor(p[fin] * NBINS) + 1L, NBINS)
      isd <- rms[fin] > C_ROC
      if (any(isd))  cl$dec <- cl$dec + tabulate(pb[isd],  nbins = NBINS)
      if (any(!isd)) cl$sim <- cl$sim + tabulate(pb[!isd], nbins = NBINS)
    }
  }
  cells[[key]] <- cl
  if (i %% 25L == 0L) message(sprintf("  ... %d / %d files", i, nrow(meta)))
}
if (length(corrupt))
  warning("Skipped unreadable file(s): ", paste(corrupt, collapse = ", "))

kk <- function(key) strsplit(key, "|", fixed = TRUE)[[1]]

profile <- bind_rows(lapply(names(cells), function(key) {
  cl <- cells[[key]]; k <- kk(key)
  data.frame(struct = k[1], phi = as.numeric(k[2]), snr = as.numeric(k[3]),
             bin = factor(DLAB, levels = DLAB), delta = DMID,
             n = cl$n_tot, n_rej = cl$n_rej,
             rate = ifelse(cl$n_tot > 0, cl$n_rej / cl$n_tot, NA_real_),
             stringsAsFactors = FALSE)
}))

trapz_auc <- function(fpr, tpr) {
  o <- order(fpr, tpr); x <- c(0, fpr[o], 1); y <- c(0, tpr[o], 1)
  sum(diff(x) * (utils::head(y, -1) + utils::tail(y, -1)) / 2)
}
roc <- bind_rows(lapply(names(cells), function(key) {
  cl <- cells[[key]]; k <- kk(key)
  nd <- sum(cl$dec); ns <- sum(cl$sim); if (nd == 0 || ns == 0) return(NULL)
  data.frame(struct = k[1], phi = as.numeric(k[2]), snr = as.numeric(k[3]),
             alpha = seq_len(NBINS) / NBINS,
             tpr = cumsum(cl$dec) / nd, fpr = cumsum(cl$sim) / ns,
             stringsAsFactors = FALSE)
}))
auc_tab <- roc |> group_by(struct, phi, snr) |>
  summarise(auc = trapz_auc(fpr, tpr), .groups = "drop")

saveRDS(list(profile = profile, roc = roc, auc = auc_tab,
             alpha = ALPHA, c_roc = C_ROC, corrupt = corrupt),
        file.path(OUT_DIR, "localization.rds"))

cat("\n== False positive rate where trends are exactly identical (delta_t = 0) ==\n")
print(as.data.frame(profile |> filter(bin == "0") |>
        mutate(fpr = round(rate, 4)) |> select(struct, phi, snr, fpr, n) |>
        arrange(struct, phi, snr)), row.names = FALSE)

cat("\n== Detection profile pooled over cells ==\n")
print(as.data.frame(profile |> group_by(bin) |>
        summarise(rate = round(sum(n_rej) / sum(n), 3), n = sum(n),
                  .groups = "drop")), row.names = FALSE)

cat(sprintf("\n== Reference ROC AUC (truth cutoff %.2f; see header caveat) ==\n", C_ROC))
print(as.data.frame(auc_tab |> mutate(auc = round(auc, 3)) |>
        arrange(struct, phi, snr)), row.names = FALSE)

## ---- figure -------------------------------------------------------------------
lab_struct <- c(smooth = "Smooth", cross = "Cross", rate = "Rate")
pal <- c(Smooth = "#0072B2", Cross = "#D55E00", Rate = "#009E73")

pdat <- profile |> filter(!is.na(rate), n >= 500) |>
  mutate(Structure = factor(lab_struct[struct], levels = names(pal)))

p <- ggplot(pdat, aes(delta, rate, colour = Structure)) +
  geom_hline(yintercept = ALPHA, linetype = "dashed",
             colour = "grey60", linewidth = 0.3) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.1) +
  facet_grid(phi ~ snr, labeller = labeller(
    phi = function(x) paste0("phi == ", x),
    snr = function(x) paste0("SNR == ", x),
    .default = label_parsed)) +
  scale_colour_manual(values = pal) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, ALPHA, 0.5, 1),
                     labels = c("0", ".05", ".5", "1")) +
  labs(x = expression("Local RMS trend separation " * delta[t]),
       y = "Rejection rate", colour = NULL,
       caption = paste0(
         "Leftmost point (delta_t = 0) is the false positive rate where the ",
         "trends are exactly identical; dashed line marks alpha = 0.05.\n",
         "Profiles pool over separation d and series length T.")) +
  theme_bw(base_size = 11) +
  theme(legend.position  = "bottom",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.15, colour = "grey85"),
        strip.background = element_rect(fill = "grey95", colour = NA),
        plot.caption     = element_text(size = 7, hjust = 0, colour = "grey30"))

ggsave(file.path(IMG_DIR, "fig_localization.png"), p,
       width = 7.5, height = 5.5, dpi = 400)
cat(sprintf("\nWrote %s\n", file.path(IMG_DIR, "fig_localization.png")))
