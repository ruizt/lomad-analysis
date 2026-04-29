## sim.R — container entrypoint for MVN coverage simulation
##
## Reads parameters from environment variables, runs S replicates for one
## value of n, and saves the results as an .rds file.
##
## Environment variables:
##   SIM_N       — sample size (default: 30)
##   SIM_S       — number of replicates (default: 200)
##   SIM_SEED    — base seed for reproducibility (default: 7291)
##   SIM_OUT_DIR — output directory (default: /jobs/output)
##
## Test locally:
##   SIM_N=30 SIM_S=5 Rscript dev/sims/tide-example/tide/sim.R

library(mvtnorm)
library(dplyr)

# ---- Read environment variables ----------------------------------------------

n       <- as.integer(Sys.getenv("SIM_N",       "30"))
S       <- as.integer(Sys.getenv("SIM_S",       "200"))
seed0   <- as.integer(Sys.getenv("SIM_SEED",    "7291"))
out_dir <- Sys.getenv("SIM_OUT_DIR", "/jobs/output")

# ---- Fixed parameters (must match template.R) --------------------------------

p     <- 3
mu    <- c(1, -0.5, 2)
Sigma <- matrix(c(
  1.0, 0.5, 0.3,
  0.5, 2.0, 0.4,
  0.3, 0.4, 1.5
), nrow = p, byrow = TRUE)
alpha <- 0.05

# ---- Single-replicate function (copied from template.R) ----------------------

run_rep <- function(n, seed) {
  set.seed(seed)

  X     <- rmvnorm(n, mean = mu, sigma = Sigma)
  x_bar <- colMeans(X)
  S_hat <- cov(X)

  diff  <- x_bar - mu
  T2    <- n * drop(diff %*% solve(S_hat) %*% diff)

  F_stat   <- T2 * (n - p) / (p * (n - 1))
  crit     <- qf(1 - alpha, df1 = p, df2 = n - p)
  covered  <- F_stat <= crit

  data.frame(n = n, seed = seed, T2 = T2, F_stat = F_stat, covered = covered)
}

# ---- Simulation loop ---------------------------------------------------------

set.seed(seed0 + n)
seeds <- sample.int(1e6, S)

results <- bind_rows(lapply(seq_len(S), function(i) {
  run_rep(n, seeds[i])
}))

# ---- Summary -----------------------------------------------------------------

results_summary <- results |>
  group_by(n) |>
  summarise(
    S        = n(),
    coverage = mean(covered),
    se       = sqrt(coverage * (1 - coverage) / S),
    .groups  = "drop"
  )

# ---- Save results ------------------------------------------------------------

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
filename <- sprintf("n%d.rds", n)
saveRDS(list(n = n, S = S, seed0 = seed0,
             results = results, results_summary = results_summary),
        file.path(out_dir, filename))
