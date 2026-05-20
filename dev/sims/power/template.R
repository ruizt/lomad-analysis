## template.R — local proof-of-concept for power study
##
## Runs S replicates for each (structure, d, n, phi, snr) combination and reports
## detection rates, sensitivity, and FDR for structured methods.
##
## Usage (from repo root):
##   source("dev/sims/power/template.R")

library(lomad)
library(dplyr)
library(tidyverse)
devtools::load_all()
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
  rate   = list(rate = 0.01)
)

# ---- Single-replicate function -----------------------------------------------

run_rep <- function(d, struct, n, phi, snr, seed) {
  set.seed(seed)
  
  # Generate trends
  trends <- do.call(sim_trends,
                    c(list(n = n, d = d, method = struct, seed = seed),
                      struct_params[[struct]]))
  
  # Add noise
  sim <- sim_noise_pair(trends, h = h_win, lambda_target = snr,
                        ar.coefs = phi, seed = seed + 1L)
  
  # Fit + test
  fit <- tryCatch(
    lomad_fit(sim$y1, sim$y2, h = h_win, s = s_win),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(data.frame(d = d, 
                      struct = struct, 
                      n = n,
                      phi = phi,
                      snr = snr,
                      seed = seed,
                      detected = NA, 
                      sensitivity = NA_real_, 
                      fdr = NA_real_,
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
  
  data.frame(d = d, 
             n = n,
             phi = phi,
             snr = snr,
             struct = struct, 
             seed = seed,
             detected = detected, 
             sensitivity = sensitivity,
             fdr = fdr_val, 
             n_flagged = sum(rejected, na.rm = TRUE))
}

run_rep(d=0.5, struct='rate', n=500, phi=0.5, snr=1.5, seed=123)

# ---- Main loop ---------------------------------------------------------------

set.seed(2847)
all_seeds <- sample.int(1e6, S)

results <- lapply(structs, function(struct) {
  lapply(d_vals, function(d) {
    lapply(n_vals, function(n) {
      lapply(phi_vals, function(phi) {
        lapply(snr_vals, function(snr) {
          bind_rows(lapply(all_seeds, function(s) run_rep(d, struct, n, phi, snr, s)))
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
    sensitivity = mean(sensitivity, na.rm = TRUE),
    fdr         = mean(fdr, na.rm = TRUE),
    n_flagged   = mean(n_flagged, na.rm = TRUE),
    .groups     = "drop"
  )

# ---- Save --------------------------------------------------------------------

saveRDS(results, "dev/sims/power/results/results.rds")
saveRDS(results_summary, "dev/sims/power/results/results_summary.rds")

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


# Graph of Sensitivity 
results_summary |>
  mutate(
    var = (sensitivity * (1 - sensitivity)) / S,
    se  = sqrt(var)
  ) |>
  filter(
    n == n_val,
    phi == phi_val,
    snr == snr_val
  ) |>
  ggplot(aes(x = d, y = sensitivity)) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(
      ymin = sensitivity - se,
      ymax = sensitivity + se,
    ),
    width = 0.2
  ) +
  geom_smooth(se= FALSE) +
  facet_wrap(~ struct, ncol = 2) +   
  theme_minimal(base_size = 18) +
  labs(
    title = sprintf("Sensitivity by Structure 
(n = %d, phi = %.1f, snr = %.1f)", 
                    n_val, phi_val, snr_val),
    x = "Distance",
    y = "Sensitivity"
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


# Graph of FDR 
results_summary |>
  mutate(
    var = (fdr * (1 - fdr)) / S,
    se  = sqrt(var)
  ) |>
  filter(
    n == n_val,
    phi == phi_val,
    snr == snr_val
  ) |>
  ggplot(aes(x = d, y = fdr)) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(
      ymin = fdr - se,
      ymax = fdr + se,
    ),
    width = 0.2
  ) +
  geom_smooth(se=FALSE) +
  facet_wrap(~ struct, ncol = 2) +   
  theme_minimal(base_size = 18) +
  labs(
    title = sprintf("FDR by Structure 
(n = %d, phi = %.1f, snr = %.1f)", 
    n_val, phi_val, snr_val),
    x = "Distance",
    y = "FDR"
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

# Graph of n flagged 
results_summary |>
  mutate(
    var = (n_flagged * (1 - n_flagged)) / S,
    se  = sqrt(var)
  ) |>
  filter(
    n == n_val,
    phi == phi_val,
    snr == snr_val
  ) |>
  ggplot(aes(x = d, y = n_flagged)) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(
      ymin = n_flagged - se,
      ymax = n_flagged + se,
    ),
    width = 0.2
  ) +
  geom_smooth(se=FALSE) +
  facet_wrap(~ struct, ncol = 2) +   
  theme_minimal(base_size = 18) +
  labs(
    title = sprintf("Number of Values Flagged by Structure
(n = %d, phi = %.1f, snr = %.1f)", 
                    n_val, phi_val, snr_val),
    x = "Distance",
    y = "Values Flagged"
  ) +
  scale_x_continuous(breaks = scales::breaks_width(0.5)) +
  theme(
    plot.title = element_text(size = 22),
    strip.text = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    panel.grid.minor = element_blank()
  )

