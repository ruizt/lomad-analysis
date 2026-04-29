library(lomad)

# Parameters (override via environment variables)
# SIM_H=0 and SIM_S=0 mean "auto-select" (recommended)
n      <- as.integer(Sys.getenv("SIM_N",      "500"))
d      <- as.numeric(Sys.getenv("SIM_D",      "5"))
seed   <- as.integer(Sys.getenv("SIM_SEED",   "32026"))
h_env  <- as.integer(Sys.getenv("SIM_H",      "0"))   # 0 = auto
s_env  <- as.integer(Sys.getenv("SIM_S",      "0"))   # 0 = auto
max_pq <- as.integer(Sys.getenv("SIM_MAX_PQ", "3"))

# Pass NULL to lomad_fit_clt() to trigger auto window selection
h_arg <- if (h_env == 0L) NULL else h_env
s_arg <- if (s_env == 0L) NULL else s_env

cat(sprintf("=== lomad CLT simulation ===\n"))
cat(sprintf("n=%d  d=%.1f  seed=%d  max_pq=%d\n", n, d, seed, max_pq))
cat(sprintf("h=%s  s=%s\n",
            if (is.null(h_arg)) "auto" else as.character(h_arg),
            if (is.null(s_arg)) "auto" else as.character(s_arg)))

# Step 1: generate trends
trends <- make_trends_dist(n = n, d = d, seed = seed)

# Step 2: add ARMA noise.
# h for add_noise (MA window for noise injection) uses the same auto rule as
# lomad_fit_clt so the noise scale is consistent with what the CLT estimator
# sees. scale=5 gives a moderate SNR; lambda_target=4 sets the noise level.
h_noise <- if (!is.null(h_arg)) h_arg else max(5L, floor(n / 200L))
sim <- add_noise(trends,
                 h             = h_noise,
                 lambda_target = 4,
                 scale         = 5,
                 order         = c(2, 1),
                 s             = 100,
                 n_start       = 30,
                 seed          = seed)

cat("Noise summary:\n")
print(sim$noise)

# Step 3: fit the CLT model.
# lomad_fit_clt() estimates the local population correlation rho_t and its
# asymptotic variance V_t at each time point using a CLT approximation —
# much faster than bootstrap and no random variation across runs.
fit <- lomad_fit_clt(sim$y1, sim$y2, h = h_arg, s = s_arg, max_pq = max_pq)

cat(sprintf("Fitted h=%d  s=%d\n", fit$inputs$h, fit$inputs$s))
cat(sprintf("ARMA orders: series1 ARMA(%d,%d)  series2 ARMA(%d,%d)\n",
            fit$noise$series1$order[1], fit$noise$series1$order[3],
            fit$noise$series2$order[1], fit$noise$series2$order[3]))

# Step 4: pointwise CLT test with Benjamini-Yekutieli FDR control.
# Tests H0: R_t >= rho_t (series are locally similar) at each valid time
# point. Rejection means the observed local correlation has dropped
# significantly below the expected baseline — evidence of decoupling.
tst <- lomad_test_clt(fit, alpha = 0.05)

n_valid    <- length(fit$valid_idx)
n_rejected <- sum(tst$rejected, na.rm = TRUE)
cat(sprintf("Rejected %d / %d valid time points (alpha_eff = %.5f)\n",
            n_rejected, n_valid, tst$alpha_eff))

# Save full output (fit + test together)
out_dir <- Sys.getenv("SIM_OUT_DIR", ".")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_file <- file.path(out_dir, sprintf("out_n%d_d%.1f_seed%d.rds", n, d, seed))
saveRDS(list(fit = fit, tst = tst), out_file)
cat(sprintf("Saved -> %s\n", out_file))

# Append one summary row to the shared CSV log.
# Columns match dev/sims/calibration/run.R so results can be compared directly.
log_row <- data.frame(
  timestamp   = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  seed        = seed,
  n           = n,
  d           = d,
  h           = fit$inputs$h,
  s           = fit$inputs$s,
  max_pq      = max_pq,
  arma_p1     = fit$noise$series1$order[1],
  arma_q1     = fit$noise$series1$order[3],
  arma_p2     = fit$noise$series2$order[1],
  arma_q2     = fit$noise$series2$order[3],
  n_valid     = n_valid,
  n_rejected  = n_rejected,
  rej_rate    = round(n_rejected / max(n_valid, 1L), 6),
  alpha_eff   = round(tst$alpha_eff, 6),
  out_file    = out_file
)

log_file <- file.path(out_dir, "lomad_results_log.csv")
write.table(
  log_row,
  file      = log_file,
  sep       = ",",
  row.names = FALSE,
  col.names = !file.exists(log_file),
  append    = TRUE
)
cat(sprintf("Logged  -> %s\n", log_file))
