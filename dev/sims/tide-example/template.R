## template.R — local proof-of-concept for MVN coverage simulation
##
## Draws S replicates of n i.i.d. MVN(mu, Sigma) samples for each sample size
## in n_vals. For each replicate, builds a 95% confidence ellipsoid for the
## mean and checks whether the true mu is inside. Reports coverage rates.
##
## Usage (from repo root):
##   source("dev/sims/tide-example/template.R")

library(mvtnorm)
library(dplyr)

# ---- Parameters --------------------------------------------------------------

p     <- 3                         # dimension
mu    <- c(1, -0.5, 2)            # true mean
Sigma <- matrix(c(                 # true covariance
  1.0, 0.5, 0.3,
  0.5, 2.0, 0.4,
  0.3, 0.4, 1.5
), nrow = p, byrow = TRUE)

alpha  <- 0.05                     # nominal miscoverage
n_vals <- c(10, 30, 100, 500)      # sample sizes to sweep
S      <- 200                      # replicates per n

# ---- Single-replicate function -----------------------------------------------

run_rep <- function(n, seed) {
  set.seed(seed)

  X     <- rmvnorm(n, mean = mu, sigma = Sigma)
  x_bar <- colMeans(X)
  S_hat <- cov(X)

  # Hotelling's T^2 statistic for H0: mean = mu
  diff  <- x_bar - mu
  T2    <- n * drop(diff %*% solve(S_hat) %*% diff)

  # Under H0: T^2 ~ p(n-1)/(n-p) * F(p, n-p)
  F_stat   <- T2 * (n - p) / (p * (n - 1))
  crit     <- qf(1 - alpha, df1 = p, df2 = n - p)
  covered  <- F_stat <= crit

  data.frame(n = n, seed = seed, T2 = T2, F_stat = F_stat, covered = covered)
}

# testing
run_rep(n = 50, seed = 123)

# ---- Main loop ---------------------------------------------------------------

set.seed(7291)
all_seeds <- sample.int(1e6, max(n_vals) * S)

results <- lapply(n_vals, function(n) {
  seeds <- all_seeds[seq_len(S)]
  bind_rows(lapply(seeds, function(s) run_rep(n, s)))
}) |> bind_rows()

# ---- Summary -----------------------------------------------------------------

results_summary <- results |>
  group_by(n) |>
  summarise(
    S        = n(),
    coverage = mean(covered),
    se       = sqrt(coverage * (1 - coverage) / S),
    .groups  = "drop"
  )

# ---- Save --------------------------------------------------------------------

saveRDS(results, "dev/sims/tide-example/results/template-results.rds")
saveRDS(results_summary, "dev/sims/tide-example/results/template-results-summary.rds")
