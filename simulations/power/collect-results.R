## collect-results.R — assemble per-cell .rds files into the compiled artifacts
##
## Run after fetch.sh has copied the per-cell files locally. Does no plotting;
## figures are built by simulations/simulation-results.R.
##
## Four artifacts, all keyed by the cell (struct, s, n, snr, phi), pooled over d:
##   simulations-power-results.rds — one row per replicate
##   simulations-power-curves.rds  — rejection rate by delta_t (local power)
##   simulations-power-roc.rds     — precision and NPV vs the threshold c
##   simulations-power-auc.rds     — area under the precision / 1 - NPV curve
##   simulations-power-fdr.rds     — realized FDR under the BY adjustment
##
## Usage (from the repo root):
##   Rscript simulations/power/collect-results.R

library(dplyr)

RAW_DIR <- Sys.getenv("RAW_DIR", "simulations/power/results/_raw")
OUT_DIR <- "simulations/power/results"

MESH     <- 100L     # mesh cells for the power curve
NBIN     <- 1000L    # mesh cells for the threshold sweep
MIN_CELL <- 10000L   # minimum per class for a cut to be kept

# delta_t is never exactly zero under the trend construction, so a window
# counts as a true null when delta_t <= EPS. Reported across a range because
# the share of windows so classified rises with it.
EPS <- c(0.001, 0.005, 0.010, 0.020, 0.050)

# ---- Collect ----------------------------------------------------------------

all_files <- list.files(RAW_DIR, pattern = "\\.rds$", full.names = TRUE)
sum_files <- all_files[!grepl("-windows\\.rds$", all_files)]
win_files <- all_files[grepl("-windows\\.rds$", all_files)]

if (length(sum_files) == 0) {
  stop("No .rds files found in: ", RAW_DIR,
       "\nCheck that the PVC is mounted and all jobs have completed.")
}

results <- lapply(sum_files, function(f) readRDS(f)$results) |> bind_rows()

# Mesh index of delta: 1..B over [0, 1), and B + 1 for delta == 1, which is an
# atom and so gets a cell of its own.
mesh_of <- function(x, B)
  ifelse(x >= 1, B + 1L, pmin(B, pmax(1L, as.integer(ceiling(x * B)))))

# Reduce each window file to counts as it is read: the sweep holds ~194M rows,
# too many to bind into one frame.
#
# Each window file is paired with the summary file of the same name to recover
# its cell. Do not join on `seed` instead -- sim.R derives seeds from `d` alone,
# so every cell sharing a `d` draws the same ones.
tally <- lapply(win_files, function(f) {
  meta <- readRDS(sub("-windows\\.rds$", ".rds", f))
  w    <- readRDS(f)
  w    <- w[is.finite(w$delta) & is.finite(w$rejected), ]

  coarse <- mesh_of(w$delta, MESH)
  fine   <- mesh_of(w$delta, NBIN)

  # counts per (replicate, mesh cell), kept split for the variance below
  rep_id <- match(w$seed, unique(w$seed))
  k      <- length(unique(w$seed))
  NCELL  <- MESH + 1L                              # mesh cells plus the atom
  idx    <- (rep_id - 1L) * NCELL + coarse         # replicate-major flat index

  list(
    cell = data.frame(struct = meta$struct, s = meta$s, n = meta$n,
                      snr = meta$snr, phi = meta$phi),
    k         = k,
    n_ij      = matrix(tabulate(idx, nbins = k * NCELL), nrow = k, byrow = TRUE),
    r_ij      = matrix(tabulate(idx[w$rejected], nbins = k * NCELL),
                       nrow = k, byrow = TRUE),
    delta_sum = as.numeric(tapply(w$delta, factor(coarse, levels = seq_len(NCELL)),
                                  sum, default = 0)),
    # NBIN + 1 so the delta == 1 atom lands in its own top rank rather than
    # being discarded by tabulate()
    rej_fine  = tabulate(fine[w$rejected],  nbins = NBIN + 1L),
    non_fine  = tabulate(fine[!w$rejected], nbins = NBIN + 1L),

    # Per-replicate false discovery proportion, V / max(R, 1), summed here so
    # the pooled mean and its standard error reconstruct without holding every
    # replicate. Uses the stored decisions, which are the BY step-up applied
    # once per replicate over that replicate's full family of windows.
    fdr = do.call(rbind, lapply(EPS, function(e) {
      nul <- w$delta <= e
      R   <- tapply(w$rejected,       rep_id, sum)
      V   <- tapply(w$rejected & nul, rep_id, sum)
      fdp <- as.numeric(V) / pmax(as.numeric(R), 1)
      data.frame(eps = e, replicates = k,
                 sum_fdp = sum(fdp), sum_fdp2 = sum(fdp^2),
                 sum_null = sum(as.numeric(tapply(nul, rep_id, mean))),
                 windows = nrow(w), rejected = sum(w$rejected))
    }))
  )
})

cell_key <- vapply(tally, function(x) paste(x$cell, collapse = "|"), character(1))

# ---- Local power curves -----------------------------------------------------

# Rejection rate per mesh cell as a ratio estimator over replicates, with the
# cluster variance for unequal cluster sizes:
#
#   p_hat = sum_i r_i / sum_i n_i,
#   Var   = k/(k-1) * sum_i (r_i - p_hat n_i)^2 / (sum_i n_i)^2,
#
# where r_i and n_i are replicate i's rejected and total counts in the cell.
curves <- lapply(split(tally, cell_key), function(g) {
  n_ij <- do.call(rbind, lapply(g, `[[`, "n_ij"))     # replicates x mesh
  r_ij <- do.call(rbind, lapply(g, `[[`, "r_ij"))
  k    <- nrow(n_ij)

  n_j <- colSums(n_ij)
  r_j <- colSums(r_ij)
  p_j <- ifelse(n_j > 0, r_j / n_j, NA_real_)

  resid <- r_ij - rep(p_j, each = k) * n_ij
  se_j  <- ifelse(n_j > 0,
                  sqrt(k / (k - 1) * colSums(resid^2)) / n_j,
                  NA_real_)

  cbind(g[[1]]$cell,
        data.frame(
          mesh       = seq_len(MESH + 1L),
          atom       = seq_len(MESH + 1L) == MESH + 1L,
          windows    = n_j,
          replicates = k,
          rejection  = p_j,
          delta_mean = Reduce(`+`, lapply(g, `[[`, "delta_sum")) / pmax(n_j, 1),
          se         = se_j))
}) |> bind_rows() |> filter(windows > 0)

# ---- Classification accuracy against the decoupling threshold ---------------

# Sweep the threshold c. At each cut, with truth = {delta_t > c} and prediction
# = rejection,
#
#   prec = P(delta > c | reject),   npv = P(delta <= c | no reject).
#
# Cuts leaving fewer than MIN_CELL windows in either class are dropped.
roc <- lapply(split(tally, cell_key), function(g) {
  r <- Reduce(`+`, lapply(g, `[[`, "rej_fine"))
  m <- Reduce(`+`, lapply(g, `[[`, "non_fine"))
  R <- as.numeric(sum(r)); M <- as.numeric(sum(m))

  # counts strictly above each cut, taken over the cut points between cells
  r_above <- rev(cumsum(rev(as.numeric(r))))[-1]
  m_above <- rev(cumsum(rev(as.numeric(m))))[-1]
  n_above <- r_above + m_above
  n_below <- (R + M) - n_above

  keep <- n_above >= MIN_CELL & n_below >= MIN_CELL
  cbind(g[[1]]$cell,
        data.frame(c    = (seq_along(r_above) / NBIN)[keep],
                   prec = (r_above / R)[keep],
                   npv  = ((M - m_above) / M)[keep]))
}) |> bind_rows()

# Add the c = 0 and c = 1 endpoints, which the sweep cannot evaluate because
# one class is empty there. Both are exact.
roc <- roc |>
  group_by(struct, s, n, snr, phi) |>
  group_modify(~ bind_rows(data.frame(c = 0, prec = 1, npv = 0),
                           arrange(.x, c),
                           data.frame(c = 1, prec = 0, npv = 1))) |>
  ungroup()

# ---- Concordance ------------------------------------------------------------

# Total windows and overall rejection rate per cell.
totals <- lapply(split(tally, cell_key), function(g) {
  R <- as.numeric(sum(Reduce(`+`, lapply(g, `[[`, "rej_fine"))))
  M <- as.numeric(sum(Reduce(`+`, lapply(g, `[[`, "non_fine"))))
  cbind(g[[1]]$cell,
        data.frame(windows = R + M, rejection = R / (R + M)))
}) |> bind_rows()

# Trapezoid area under the precision / 1 - NPV curve above.
auc <- roc |>
  group_by(struct, s, n, snr, phi) |>
  summarise(auc = { x <- 1 - npv; y <- prec; o <- order(x)
                    sum(diff(x[o]) * (y[o][-1] + head(y[o], -1)) / 2) },
            .groups = "drop") |>
  right_join(totals, by = c("struct", "s", "n", "snr", "phi"))

# ---- Realized FDR under the BY adjustment ------------------------------------

# Pooled over every replicate in the design, not averaged over cells: the
# dataset is the unit, so cell averages would weight unequal replicate counts.
fdr <- lapply(tally, function(x) cbind(x$cell, x$fdr)) |>
  bind_rows() |>
  group_by(eps) |>
  summarise(datasets   = sum(replicates),
            null_share = sum(sum_null) / sum(replicates),
            fdr        = sum(sum_fdp) / sum(replicates),
            sd         = sqrt((sum(sum_fdp2) - sum(sum_fdp)^2 / sum(replicates)) /
                              (sum(replicates) - 1)),
            .groups    = "drop") |>
  mutate(se = sd / sqrt(datasets))

# ---- Save -------------------------------------------------------------------

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
saveRDS(results, file.path(OUT_DIR, "simulations-power-results.rds"))
saveRDS(curves,  file.path(OUT_DIR, "simulations-power-curves.rds"))
saveRDS(roc,     file.path(OUT_DIR, "simulations-power-roc.rds"))
saveRDS(auc,     file.path(OUT_DIR, "simulations-power-auc.rds"))
saveRDS(fdr,     file.path(OUT_DIR, "simulations-power-fdr.rds"))

cat(sprintf("%d replicates, %s windows -> %d curve rows, %d cells\n",
            nrow(results), format(sum(auc$windows), big.mark = ","), nrow(curves), nrow(auc)))
