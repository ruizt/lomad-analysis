## collect-results.R — assemble per-cell .rds files into the compiled artifacts
##
## Run after fetch.sh has copied the per-cell files locally. Produces everything
## the figure stage consumes, so the pipeline is raw per-cell files -> compiled
## artifacts in one step. Deliberately does no plotting: figures are built by
## simulations/simulation-results.R.
##
## Three artifacts:
##   simulations-power-results.rds — one row per replicate
##   simulations-power-curves.rds  — rejection rate by delta_t (local power)
##   simulations-power-auc.rds     — concordance between rejection and delta_t
##
## The study reports rejection against *realized* local separation delta_t, so
## d is pooled over rather than plotted: it is a generator knob scaled per
## structure, and the same d means different separations for different
## structures. Pooling is what makes the curves comparable across them.
##
## This replaces localization-sweep.R, which had to reread every full-length
## series file to recover (delta_t, rejected) pairs. tide/sim.R now emits those
## directly, so the expensive stage is gone.
##
## Usage (from the repo root):
##   Rscript simulations/power/collect-results.R

library(dplyr)

RAW_DIR <- Sys.getenv("RAW_DIR", "simulations/power/results/_raw")
OUT_DIR <- "simulations/power/results"

# Resolution of the local power curve: a uniform mesh over [0, 1] rather than a
# handful of bins. At this sweep's size every mesh cell still holds tens of
# thousands of windows -- 20k at the thinnest -- so the curve is effectively
# continuous without any smoothing choice to defend.
MESH <- 100L

# ---- Collect ----------------------------------------------------------------

all_files <- list.files(RAW_DIR, pattern = "\\.rds$", full.names = TRUE)
sum_files <- all_files[!grepl("-windows\\.rds$", all_files)]
win_files <- all_files[grepl("-windows\\.rds$", all_files)]

if (length(sum_files) == 0) {
  stop("No .rds files found in: ", RAW_DIR,
       "\nCheck that the PVC is mounted and all jobs have completed.")
}

results <- lapply(sum_files, function(f) readRDS(f)$results) |> bind_rows()

# The sweep holds ~194 million window rows, which cannot be bound into one
# frame: reading them all at once exhausts memory well past 16 GB. Nothing
# downstream needs the rows themselves, only counts per delta_t bin, so each
# file is reduced as it is read and only the tallies are kept.
#
# NBIN is the resolution the concordance is computed at, finer than the curve's:
# 1000 bins over [0, 1] makes the tie correction below negligible.
NBIN <- 1000L

# Mesh index of a delta value: 1..B over [0, 1), and B + 1 for delta == 1.
#
# delta = 1 is an atom, not the end of a continuum. Under the constraint b > 0
# the minimiser sits at the boundary whenever the window correlation is
# non-positive, so every window whose trends reverse maps to exactly 1 -- about
# one in eight. Folding those into the top mesh cell would mix a point mass
# with a sliver of genuine spread and read as the curve turning up at the edge.
mesh_of <- function(x, B)
  ifelse(x >= 1, B + 1L, pmin(B, pmax(1L, as.integer(ceiling(x * B)))))

# Window files carry no cell metadata of their own, and cannot be joined on
# `seed`: sim.R derives its seeds from `d` alone, so every (s, snr, phi) cell
# sharing a `d` draws the *same* seeds. Each window file is therefore paired
# with the summary file of the same name, and the cell taken from there.
tally <- lapply(win_files, function(f) {
  meta <- readRDS(sub("-windows\\.rds$", ".rds", f))
  w    <- readRDS(f)
  w    <- w[is.finite(w$delta) & is.finite(w$rejected), ]

  coarse <- mesh_of(w$delta, MESH)
  fine   <- mesh_of(w$delta, NBIN)

  # Per (replicate, mesh cell) counts, not pooled. Windows overlap by s - 1, so
  # they are not independent and a binomial interval on their raw count
  # understates the uncertainty severalfold. The replicate is the independent
  # unit -- its own trends, its own noise -- so the counts are kept split by
  # seed and the variance is formed across replicates further down.
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
    lam_sum   = as.numeric(tapply(w$lambda, factor(coarse, levels = seq_len(NCELL)),
                                  function(v) sum(v, na.rm = TRUE), default = 0)),
    # NBIN + 1 so the delta == 1 atom lands in its own top rank rather than
    # being discarded by tabulate()
    rej_fine  = tabulate(fine[w$rejected],  nbins = NBIN + 1L),
    non_fine  = tabulate(fine[!w$rejected], nbins = NBIN + 1L)
  )
})

cell_key <- vapply(tally, function(x) paste(x$cell, collapse = "|"), character(1))

# ---- Local power curves -----------------------------------------------------

# Pooled over d, which only decides which part of the delta_t range gets
# populated, but *not* over replicate: the replicate is the sampling unit.
#
# Windows overlap by s - 1, so a binomial interval on their raw count treats
# near-duplicates as independent observations and understates the standard
# error by a factor of two to four. The estimate is instead a ratio estimator
# over replicates, with the usual cluster variance for unequal cluster sizes,
#
#   p_hat = sum_i r_i / sum_i n_i,
#   Var   = k/(k-1) * sum_i (r_i - p_hat n_i)^2 / (sum_i n_i)^2,
#
# where r_i and n_i are replicate i's rejected and total window counts in the
# mesh cell. This makes no assumption about the dependence within a replicate,
# and reduces to the binomial form when there is none.
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
          mesh        = seq_len(MESH + 1L),
          atom        = seq_len(MESH + 1L) == MESH + 1L,
          windows     = n_j,
          replicates  = k,
          rejection   = p_j,
          delta_mean  = Reduce(`+`, lapply(g, `[[`, "delta_sum")) / pmax(n_j, 1),
          lambda_mean = Reduce(`+`, lapply(g, `[[`, "lam_sum")) / pmax(n_j, 1),
          se          = se_j))
}) |> bind_rows() |> filter(windows > 0)

# ---- Concordance between rejection and true separation ----------------------

# P(a randomly chosen rejected window is more separated than a randomly chosen
# non-rejected one). Computed from the binned counts rather than by ranking
# 194 million values: with r_k rejected and m_k non-rejected in bin k,
#   U = sum_k r_k (sum_{j<k} m_j + m_k / 2),
# which is Mann-Whitney with ties inside a bin taking the usual half credit.
auc <- lapply(split(tally, cell_key), function(g) {
  r <- Reduce(`+`, lapply(g, `[[`, "rej_fine"))
  m <- Reduce(`+`, lapply(g, `[[`, "non_fine"))
  # as.numeric throughout: R * M reaches ~1e11 here, past integer range
  R <- as.numeric(sum(r)); M <- as.numeric(sum(m))
  below <- cumsum(as.numeric(m)) - m
  cbind(g[[1]]$cell,
        data.frame(windows   = R + M,
                   rejection = R / (R + M),
                   auc       = if (R == 0 || M == 0) NA_real_
                               else sum(r * (below + m / 2)) / (R * M)))
}) |> bind_rows()

# ---- Save -------------------------------------------------------------------

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
saveRDS(results, file.path(OUT_DIR, "simulations-power-results.rds"))
saveRDS(curves,  file.path(OUT_DIR, "simulations-power-curves.rds"))
saveRDS(auc,     file.path(OUT_DIR, "simulations-power-auc.rds"))

cat(sprintf("%d replicates, %s windows -> %d curve rows, %d cells\n",
            nrow(results), format(sum(auc$windows), big.mark = ","), nrow(curves), nrow(auc)))
