## run-realignment.R -- does removing the affine map first make trend
## comparison work under local affine similarity? See design.md for the cells.
##
## Outputs
##   results/_raw/ms-{scenario}-{framing}.rds   per-cell MSinference cache
##   results/realignment-results.rds             compiled study object
##
## Usage (from the repo root):
##   REALIGNMENT_DRY=1 Rscript simulations/realignment/run-realignment.R  # fast
##   Rscript simulations/realignment/run-realignment.R                   # full

suppressPackageStartupMessages({library(lomad); library(MSinference)})

N <- 1250L; H <- 5L; NB <- 151L
AFFINE_S <- 50L    # window the affine drift cap applies over (data generation)
REALIGN_S <- 50L   # window for the rolling realignment given to MSinference
SNR <- 0.5; PHI <- 0.5; ALPHA <- 0.05

## Caches and compiled output are namespaced by seed, so draws never collide.
## Second draw: REALIGNMENT_TREND_SEED=7307 REALIGNMENT_NOISE_SEED=2411
TREND_SEED <- as.integer(Sys.getenv("REALIGNMENT_TREND_SEED", "6001"))
NOISE_SEED <- as.integer(Sys.getenv("REALIGNMENT_NOISE_SEED", "1001"))

SIM_RUNS <- as.integer(Sys.getenv("SIM_RUNS", "1000"))
DRY      <- nzchar(Sys.getenv("REALIGNMENT_DRY"))

OUT_DIR <- "simulations/realignment/results"
RAW_DIR <- file.path(OUT_DIR, "_raw")
dir.create(RAW_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- scenarios --------------------------------------------------------------
# `eq` reuses the `af` base trend with the affine layer removed, so the two
# differ only in the map. Both satisfy the null.
#
# The integer literals are load-bearing: sim_trends() consumes RNG differently
# for `affine_s = 50` than for `affine_s = 50L`.
tr_af <- sim_trends(n = N, d = 0, method = "rs", bw = 50, nb = NB,
                    seed = TREND_SEED, affine_s = AFFINE_S, affine_cap = 0.015)
tr_eq <- tr_af
tr_eq$x2 <- tr_eq$x1
tr_eq$a  <- rep(0, N)
tr_eq$b  <- rep(1, N)

sim_pair <- function(tr) {
  sim <- suppressMessages(sim_noise_pair(tr, h = H, lambda_target = SNR,
                                         ar.coefs = PHI, seed = NOISE_SEED))
  list(y1 = sim$y1, y2 = sim$y2, tr = tr)
}
scen <- list(eq = sim_pair(tr_eq), af = sim_pair(tr_af))

# ---- lomad ------------------------------------------------------------------
# Run across window sizes; s = 50 is the hardest cell of the power study.
LOMAD_CELLS <- list(
  list(scenario = "eq", s = 50L),
  list(scenario = "af", s = 50L),
  list(scenario = "af", s = 80L),
  list(scenario = "af", s = 100L)
)

run_lomad <- function(cell) {
  sc <- scen[[cell$scenario]]
  f  <- suppressMessages(suppressWarnings(lomad_fit(sc$y1, sc$y2, h = H, s = cell$s)))
  ts <- suppressMessages(lomad_test(f, alpha = ALPHA))
  vi <- f$valid_idx
  nrej <- sum(ts$rejected[vi], na.rm = TRUE)
  list(scenario = cell$scenario, method = "lomad",
       framing = sprintf("s = %d", cell$s), pending = FALSE,
       reject = nrej > 0, flagged = nrej, total = length(vi))
}

# ---- realignment ------------------------------------------------------------
# Centered rolling regression of y on x over windows of width REALIGN_S;
# unconstrained apart from a finiteness guard.
roll_ab <- function(x, y, s) {
  n <- length(x); a <- numeric(n); b <- numeric(n)
  half <- s %/% 2L
  for (t in seq_len(n)) {
    lo <- max(1L, t - half); hi <- min(n, lo + s - 1L); lo <- max(1L, hi - s + 1L)
    w  <- (lo:hi)[stats::complete.cases(x[lo:hi], y[lo:hi])]
    vx <- stats::var(x[w])
    bt <- if (isTRUE(vx > 0)) stats::cov(x[w], y[w]) / vx else NA_real_
    if (!is.finite(bt) || bt == 0) bt <- 1
    b[t] <- bt
    a[t] <- mean(y[w]) - bt * mean(x[w])
  }
  list(a = a, b = b)
}

ma <- function(x, h) as.numeric(stats::filter(x, rep(1 / h, h), sides = 1))

realign <- function(sc, from = c("raw", "smoothed")) {
  from <- match.arg(from)
  map <- if (from == "raw") roll_ab(sc$y1, sc$y2, REALIGN_S)
         else roll_ab(ma(sc$y1, H), ma(sc$y2, H), REALIGN_S)
  list(y2_re = (sc$y2 - map$a) / map$b,
       map = map,
       diag = list(b_mean = mean(map$b), b_range = range(map$b),
                   b_sd = stats::sd(map$b), b_negative = sum(map$b < 0),
                   true_b_range = range(sc$tr$b)))
}

# ---- MSinference cells ------------------------------------------------------
zs   <- function(v) (v - mean(v)) / stats::sd(v)
grid <- MSinference::construct_grid(N)

run_ms <- function(y1, y2) {
  a <- zs(y1); b <- zs(y2)
  lrv <- c(MSinference::estimate_lrv(a, q = 25, r_bar = 10, p = 1)$lrv,
           MSinference::estimate_lrv(b, q = 25, r_bar = 10, p = 1)$lrv)
  t0 <- Sys.time()
  r <- MSinference::multiscale_test(data = cbind(a, b), sigma_vec = sqrt(lrv),
                                    n_ts = 2, grid = grid, alpha = ALPHA,
                                    sim_runs = SIM_RUNS)
  gv <- r$gset_with_values[[1]]
  # r$stat is max() over a pairwise matrix padded with structural zeros, which
  # floors it at 0; report the true (1,2) pairwise statistic instead. The test
  # decision is unaffected since the critical value is positive.
  sp <- r$stat_pairwise
  stat <- if (is.matrix(sp) && all(dim(sp) == 2)) sp[1, 2] else r$stat
  list(stat = stat, quant = r$quant, reject = r$stat > r$quant,
       flagged = sum(gv$test != 0), total = nrow(gv),
       elapsed_s = as.numeric(difftime(Sys.time(), t0, units = "secs")))
}

MS_CELLS <- list(
  list(scenario = "eq", framing = "global",         prep = function(sc) list(y1 = sc$y1, y2 = sc$y2)),
  list(scenario = "af", framing = "global",         prep = function(sc) list(y1 = sc$y1, y2 = sc$y2)),
  list(scenario = "af", framing = "realign-raw",    prep = function(sc) { re <- realign(sc, "raw");      list(y1 = sc$y1, y2 = re$y2_re, diag = re$diag) }),
  list(scenario = "af", framing = "realign-smooth", prep = function(sc) { re <- realign(sc, "smoothed"); list(y1 = sc$y1, y2 = re$y2_re, diag = re$diag) })
)

# Default-seed caches keep their original names; other seeds get a suffix.
seed_tag <- {
  if (TREND_SEED == 6001L && NOISE_SEED == 1001L) ""
  else sprintf("-seed%d.%d", TREND_SEED, NOISE_SEED)
}

run_cell <- function(cell, idx) {
  cache <- file.path(RAW_DIR, sprintf("ms-%s-%s%s.rds", cell$scenario,
                                      cell$framing, seed_tag))
  if (file.exists(cache)) {
    out <- readRDS(cache)
    if (out$sim_runs < SIM_RUNS)
      warning(sprintf("%s was cached at sim_runs = %d < %d; delete %s to rerun",
                      basename(cache), out$sim_runs, SIM_RUNS, cache), call. = FALSE)
    return(out)
  }
  prep <- cell$prep(scen[[cell$scenario]])
  if (DRY)
    return(list(scenario = cell$scenario, framing = cell$framing,
                method = "MSinference", pending = TRUE, diag = prep$diag))
  set.seed(TREND_SEED + idx)
  message(sprintf("running MSinference: %s / %s (sim_runs = %d) ...",
                  cell$scenario, cell$framing, SIM_RUNS))
  out <- c(list(scenario = cell$scenario, framing = cell$framing,
                method = "MSinference", pending = FALSE,
                sim_runs = SIM_RUNS, diag = prep$diag,
                trend_seed = TREND_SEED, noise_seed = NOISE_SEED,
                msinference_version = as.character(utils::packageVersion("MSinference")),
                run_date = Sys.Date()),
           run_ms(prep$y1, prep$y2))
  saveRDS(out, cache)
  message(sprintf("  done in %.0f s -> %s", out$elapsed_s, cache))
  out
}

# ---- run and compile --------------------------------------------------------
lomad_rows <- lapply(LOMAD_CELLS, run_lomad)
ms_rows    <- Map(run_cell, MS_CELLS, seq_along(MS_CELLS))

row_df <- function(r) data.frame(
  scenario = r$scenario,
  # cache files written before `method` was recorded hold MSinference cells only
  method   = if (is.null(r$method)) "MSinference" else r$method,
  framing  = r$framing,
  reject   = if (isTRUE(r$pending)) NA else r$reject,
  flagged  = if (isTRUE(r$pending)) NA_integer_ else r$flagged,
  total    = if (isTRUE(r$pending)) NA_integer_ else r$total,
  stat     = if (is.null(r$stat))  NA_real_ else r$stat,
  quant    = if (is.null(r$quant)) NA_real_ else r$quant)
tab <- do.call(rbind, lapply(c(lomad_rows, ms_rows), row_df))

results <- list(
  params = list(n = N, h = H, nb = NB, affine_s = AFFINE_S,
                realign_s = REALIGN_S, snr = SNR, phi = PHI,
                alpha = ALPHA, sim_runs = SIM_RUNS,
                trend_seed = TREND_SEED, noise_seed = NOISE_SEED),
  table  = tab,
  cells  = c(lomad_rows, ms_rows),
  # Recomputed rather than read from the cells, which predate this field.
  realignment = lapply(c(raw = "raw", smoothed = "smoothed"), function(k) {
    re <- realign(scen$af, k)
    list(a = re$map$a, b = re$map$b, diag = re$diag)
  }),
  scenarios = lapply(scen, function(sc)
    list(y1 = sc$y1, y2 = sc$y2, x1 = sc$tr$x1, x2 = sc$tr$x2,
         a = sc$tr$a, b = sc$tr$b)))
saveRDS(results, file.path(OUT_DIR, sprintf("realignment-results%s.rds", seed_tag)))

print(tab, row.names = FALSE, digits = 3)
if (any(is.na(tab$reject))) message("\npending cells remain; rerun without REALIGNMENT_DRY=1")
