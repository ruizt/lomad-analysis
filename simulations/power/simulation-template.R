## simulation-template.R — local illustration of the power study
##
## Runs S replicates for each (structure, d, n, phi, snr) combination so the
## simulation can be inspected and reasoned about locally, at a scale that runs
## in seconds rather than on the cluster.
##
## NOT the source of truth. tide/sim.R is what actually runs on Tide and
## produces the archived results; run_rep() below deliberately mirrors it so
## the illustration is faithful. If you change the simulation logic, change
## tide/sim.R first, then mirror it here. The two differ only in how they take
## parameters -- arguments here, environment variables there -- and in the loop
## that drives them.
##
## Writes nothing: only tide/sim.R and collect-results.R touch results/.
##
## Usage (from repo root):
##   source("simulations/power/simulation-template.R")

library(lomad)
library(dplyr)
library(tidyverse)
library(patchwork)

# ---- Parameters --------------------------------------------------------------

n     <- 500L
phi   <- 0.5
snr   <- 1.5
alpha <- 0.05
S     <- 20

d_vals     <- c(0, 0.25, 0.5, 0.75, 1, 1.5)
structs <- c("dist", "smooth", "cross", "rate")
n_vals   <- c(200L, 400L, 600L)
phi_vals <- c(0.3, 0.5, 0.8)
snr_vals <- c(0.5, 1.5)

h_win <- max(5L, floor(n / 200L))
s_win <- min(60L * h_win, floor(n / 4L))

# Structure-specific parameters
struct_params <- list(
  dist   = list(),
  smooth = list(bw = 50),
  cross  = list(bw = 50),
  # bump = "gaussian": the shape-2 gamma default puts a corner at each event
  # onset, which difference-based noise estimation cannot cancel (Hall and Van
  # Keilegom 2003 require a bounded derivative). The leftover biases the
  # residual autocovariance upward, and near the unit root that bias is
  # amplified by ~2/(1-phi)^2, so at phi = 0.8 the fixed-rate structure loses
  # almost all detection. The gaussian pulse is matched on width and smooth at
  # onset; it still trails the stochastic structures but is no longer anomalous.
  rate   = list(rate = 0.01, bump = "gaussian")
)

# ---- Single-replicate function -----------------------------------------------

run_rep <- function(d, struct, n, phi, snr, seed, oracle = FALSE) {
  set.seed(seed)

  # Generate trends
  trends <- do.call(sim_trends,
                    c(list(n = n, d = d, method = struct, seed = seed),
                      struct_params[[struct]]))

  # Coupling weight (NULL for dist/unstructured)
  w <- if (!is.null(trends$w)) trends$w else NULL

  # Add noise
  sim <- sim_noise_pair(trends, h = h_win, lambda_target = snr,
                        ar.coefs = phi, seed = seed + 1L)

  # Oracle: bypass noise estimation with the true AR params
  noise_ov <- NULL
  if (oracle) {
    z1 <- sim$y1 - sim$x1
    innov1 <- z1[-1] - phi * z1[-length(z1)]
    noise_ov <- list(ar = phi, sigma2 = var(innov1))
  }

  # Fit
  fit <- tryCatch(
    lomad_fit(sim$y1, sim$y2, h = h_win, s = s_win, noise_override = noise_ov),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(list(
      summary = data.frame(d = d, struct = struct, n = n,
                           phi = phi, snr = snr, seed = seed,
                           detected = NA),
      series = NULL
    ))
  }

  # Test
  tst <- lomad_test(fit, alpha = alpha)

  # Outputs. Schema matches tide/sim.R exactly: collect-results.R and
  # localization-sweep.R both assume it. `sep` is the pointwise true trend
  # separation and is what localization-sweep.R measures rejections against.
  list(
    summary = data.frame(d = d, n = n, phi = phi, snr = snr,
                         struct = struct, seed = seed,
                         detected = any(tst$rejected, na.rm = TRUE)),
    series = list(w = w,
                  sep = abs(trends$x1 - trends$x2),
                  vi = fit$valid_idx,
                  p_raw = tst$p_values,
                  p_adj = tst$p_adj,
                  rejected = tst$rejected)
  )
}

# summary output
run_rep(d=0.5, struct='rate', n=500, phi=0.5, snr=1.5, seed=123)$summary

# series output
run_rep(d=0.5, struct='rate', n=500, phi=0.5, snr=1.5, seed=123)$series |>
  str()

# ---- Main loop ---------------------------------------------------------------

set.seed(2847)
all_seeds <- sample.int(1e6, S)


results <- lapply(structs, function(struct) {
  lapply(d_vals, function(d) {
    lapply(n_vals, function(n) {
      lapply(phi_vals, function(phi) {
        lapply(snr_vals, function(snr) {
          reps <- lapply(all_seeds, function(s) run_rep(d, struct, n, phi, snr, s))
          
          bind_rows(lapply(reps, `[[`, "summary"))
        }) |> bind_rows()
      }) |> bind_rows()
    }) |> bind_rows()
  }) |> bind_rows()
}) |> bind_rows()


# ---- Summary -----------------------------------------------------------------

results_summary <- results |>
  group_by(struct, d, n, phi, snr) |>
  summarise(
    S           = n(),
    detection   = mean(detected, na.rm = TRUE),
    .groups     = "drop"
  )

# Nothing is written here by design: simulation-template.R is a local proof-of-concept for
# the simulation logic. Only tide/sim.R (on the cluster) and collect-results.R
# write into results/.

#------ Plot -------------------------------------------------------------------

# example
n_val <- 400
phi_val <- 0.5
snr_val <- 1.5

# Graph of Detection Rate
results_summary |>
  mutate(
    var = (detection * (1 - detection)) / S,
    se = sqrt(var)
  ) |>
  filter(
    n == n_val,
    phi == phi_val,
    snr == snr_val
  ) |>
  ggplot(aes(x = d, y = detection)) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(
      ymin = detection - se,
      ymax = detection + se
    ),
    width = 0.2
  ) +
  geom_smooth(se= FALSE) +
  facet_wrap(~ struct, ncol=2) +   
  theme_minimal(base_size = 18) +
  labs(
    title = sprintf("Detection Rate by Structure 
(n = %d, phi = %.1f, snr = %.1f)", 
                    n_val, phi_val, snr_val),
    x = "Distance",
    y = "Detection Rate"
  ) +
  scale_x_continuous(breaks = scales::breaks_width(0.5)) +
  scale_y_continuous(limits = c(0, 1)) +
  theme(
    plot.title = element_text(size = 22),
    strip.text = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    panel.grid.minor = element_blank()
  )
