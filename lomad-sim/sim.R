library(lomad)

# Parameters (override via environment variables)
n      <- as.integer(Sys.getenv("SIM_N",      "500"))
d      <- as.numeric(Sys.getenv("SIM_D",      "5"))
seed   <- as.integer(Sys.getenv("SIM_SEED",   "32026"))
h      <- as.integer(Sys.getenv("SIM_H",      "30"))
q      <- as.integer(Sys.getenv("SIM_Q",      "30"))
B      <- as.integer(Sys.getenv("SIM_B",      "1000"))
method <- Sys.getenv("SIM_METHOD", "boot")
ncores <- as.integer(Sys.getenv("SIM_NCORES", "1"))

cat(sprintf("=== lomad simulation ===\n"))
cat(sprintf("n=%d  d=%.1f  h=%d  q=%d  B=%d  method=%s  seed=%d\n",
            n, d, h, q, B, method, seed))

# Step 1: generate trends (test-script.R L3)
trends <- make_trends_dist(n = n, d = d, seed = seed)

# Step 2: add ARMA noise (test-script.R L9-17)
sim <- add_noise(trends,
                 h             = h,
                 lambda_target = 4,
                 scale         = 5,
                 order         = c(2, 1),
                 s             = 100,
                 n_start       = 30,
                 seed          = seed)

cat("Noise summary:\n")
print(sim$noise)

# Step 3: fit and test (test-script.R L30)
out <- lomad(sim$y1, sim$y2,
             q      = q,
             h      = 50,
             B      = B,
             seed   = seed,
             method = method,
             ncores = ncores)

cat("p_values:\n");  print(out$p_values)
cat("observed:\n");  print(out$observed)
cat("expected:\n");  print(out$expected)

# Save full output
out_dir <- Sys.getenv("SIM_OUT_DIR", ".")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_file <- file.path(out_dir, sprintf("out_n%d_d%.1f_seed%d.rds", n, d, seed))
saveRDS(out, out_file)
cat(sprintf("Saved -> %s\n", out_file))

# Append summary row to CSV log
log_row <- data.frame(
  timestamp       = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  seed            = seed,
  n               = n,
  d               = d,
  h               = h,
  q               = q,
  B               = B,
  method          = method,
  p_entry_rate    = round(out$p_values$entry_rate,      6),
  p_run_length    = round(out$p_values$mean_run_length, 6),
  p_frac_state    = round(out$p_values$frac_state,      6),
  p_n_entries     = round(out$p_values$n_entries,       6),
  obs_entry_rate  = round(out$observed$entry_rate,      6),
  obs_run_length  = round(out$observed$mean_run_length, 6),
  obs_frac_state  = round(out$observed$frac_state,      6),
  out_file        = out_file
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
