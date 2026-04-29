## template.R — local proof-of-concept for power study
##
## Runs S replicates for each (structure, d) combination and reports
## detection rates (and sensitivity/FDR for structured methods).
##
## Usage (from repo root):
##   source("dev/sims/power/template.R")

library(lomad)
library(dplyr)

# ---- Parameters --------------------------------------------------------------

n     <- 500L
phi   <- 0.5
snr   <- 1.5
alpha <- 0.05
S     <- 20

d_vals     <- c(0, 0.5, 1.0, 1.5, 2.0)
structures <- c("dist", "smooth", "cross", "rate")

h_win <- max(5L, floor(n / 200L))
s_win <- min(60L * h_win, floor(n / 4L))

# Structure-specific parameters
struct_params <- list(
  dist   = list(),
  smooth = list(bw = 50),
  cross  = list(bw = 50),
  rate   = list(rate = 0.01)
)

# ---- Single-replicate function -----------------------------------------------

run_rep <- function(d, structure, seed) {
  set.seed(seed)

  # Generate trends
  trends <- do.call(sim_trends,
    c(list(n = n, d = d, method = structure, seed = seed),
      struct_params[[structure]]))

  # Add noise
  sim <- sim_noise_pair(trends, h = h_win, lambda_target = snr,
                        ar.coefs = phi, seed = seed + 1L)

  # Fit + test
  fit <- tryCatch(
    lomad_fit(sim$y1, sim$y2, h = h_win, s = s_win),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(data.frame(d = d, structure = structure, seed = seed,
                      detected = NA, sensitivity = NA_real_, fdr = NA_real_,
                      n_flagged = NA_integer_))
  }

  tst <- lomad_test(fit, alpha = alpha)
  vi  <- fit$valid_idx
  rejected <- tst$rejected[vi]

  # Coupling weight (NULL for dist/unstructured)
  w <- if (!is.null(trends$w)) trends$w[vi] else NULL

  detected    <- any(rejected, na.rm = TRUE)
  sensitivity <- if (!is.null(w) && any(w < 0.5))
                   mean(rejected[w < 0.5], na.rm = TRUE) else NA_real_
  fdr_val     <- if (!is.null(w) && any(rejected, na.rm = TRUE))
                   mean(w[rejected] >= 0.5, na.rm = TRUE) else NA_real_

  data.frame(d = d, structure = structure, seed = seed,
             detected = detected, sensitivity = sensitivity,
             fdr = fdr_val, n_flagged = sum(rejected, na.rm = TRUE))
}

# ---- Main loop ---------------------------------------------------------------

set.seed(2847)
all_seeds <- sample.int(1e6, S)

results <- lapply(structures, function(struct) {
  lapply(d_vals, function(d) {
    bind_rows(lapply(all_seeds, function(s) run_rep(d, struct, s)))
  }) |> bind_rows()
}) |> bind_rows()

# ---- Summary -----------------------------------------------------------------

results_summary <- results |>
  group_by(structure, d) |>
  summarise(
    S           = n(),
    detection   = mean(detected, na.rm = TRUE),
    sensitivity = mean(sensitivity, na.rm = TRUE),
    fdr         = mean(fdr, na.rm = TRUE),
    .groups     = "drop"
  )

# ---- Save --------------------------------------------------------------------

saveRDS(results, "dev/sims/power/results/results.rds")
saveRDS(results_summary, "dev/sims/power/results/results_summary.rds")
