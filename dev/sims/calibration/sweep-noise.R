## Noise-type sweep: general AR(p) estimator diagnostics
##
## Verifies that estimate_arma_noise() (multi-lag variogram + Yule-Walker + BIC)
## recovers noise parameters and maintains calibration across four noise types:
##   AR(1), AR(2), ARMA(1,1), MA(1).
##
## For each noise type we check:
##   - Estimated vs true innovation variance (sigma2)
##   - Selected AR order
##   - Calibration (rejection rate) at d = 0
##   - Power at d = 1.5

devtools::load_all()
library(dplyr)
library(ggplot2)

# ---- Parameters ------------------------------------------------------------

n      <- 1000
snr    <- 1
alpha  <- 0.05
h_win  <- 10
s_win  <- 50

noise_types <- list(
  "AR(1)"     = list(ar = 0.5,           ma = numeric(0)),
  "AR(2)"     = list(ar = c(0.6, -0.3),  ma = numeric(0)),
  "ARMA(1,1)" = list(ar = 0.5,           ma = 0.3),
  "MA(1)"     = list(ar = numeric(0),     ma = 0.5)
)

d_sweep    <- c(0, 0.5, 1.5)
S_sweep    <- 50L
seed_sweep <- 6241

# ---- Single-replicate function ---------------------------------------------

run_rep_arma <- function(d, noise_spec, seed) {
  trends <- sim_trends(n = n, d = d, method = "dist", seed = seed)
  sim    <- suppressMessages(
    sim_noise_pair(trends, h = h_win, s = s_win, lambda_target = snr,
                   ar.coefs = noise_spec$ar, ma.coefs = noise_spec$ma,
                   seed = seed + 1L)
  )
  true_sigma2 <- mean(c(sim$noise$series1$sigma, sim$noise$series2$sigma)^2)

  # Estimate noise with the general estimator
  tr    <- estimate_trends(sim$y1, sim$y2, h_win)
  noise <- estimate_arma_noise(sim$y1, sim$y2, tr$trend)

  est_sigma2 <- mean(c(noise$series1$sigma2, noise$series2$sigma2))
  est_order  <- mean(c(noise$series1$order[1], noise$series2$order[1]))

  # Build the CLT fit, then re-derive test quantities with the general estimator
  fit <- suppressMessages(
    lomad_fit(sim$y1, sim$y2, h = h_win, s = s_win)
  )
  raw_lag_max <- fit$inputs$lag_max + h_win - 1L
  acov1_raw <- arma_acov(ar = noise$series1$ar, ma = noise$series1$ma,
                         sigma2 = noise$series1$sigma2,
                         lag_max = raw_lag_max)
  acov2_raw <- arma_acov(ar = noise$series2$ar, ma = noise$series2$ma,
                         sigma2 = noise$series2$sigma2,
                         lag_max = raw_lag_max)
  acov1 <- lomad:::.ma_filter_acov(acov1_raw, h_win, fit$inputs$lag_max)
  acov2 <- lomad:::.ma_filter_acov(acov2_raw, h_win, fit$inputs$lag_max)
  sums  <- acov_sums(acov1, acov2)

  sigma1_sq  <- acov1[1L]
  sigma2_sq  <- acov2[1L]
  noise_bias <- (sigma1_sq + sigma2_sq) / 4
  tau_sq     <- pmax(0, compute_tau_sq(tr$trend, s_win) - noise_bias)
  rho        <- compute_rho(tau_sq, sigma1_sq, sigma2_sq)
  V          <- compute_V(tau_sq, sigma1_sq, sigma2_sq,
                          sums$L1, sums$L2, sums$Q1, sums$Q2, sums$Q12)

  valid_idx <- which(is.finite(fit$R) & is.finite(rho) & is.finite(V) & V > 0)
  if (length(valid_idx) == 0L) {
    rejected <- NA
  } else {
    se <- sqrt(V[valid_idx] / s_win)
    Z  <- (fit$R[valid_idx] - rho[valid_idx]) / se
    p_raw <- pnorm(Z)
    m <- length(valid_idx)
    alpha_eff <- alpha / sum(1 / seq_len(m))
    rejected  <- any(p_raw <= alpha_eff, na.rm = TRUE)
  }

  data.frame(d = d, seed = seed,
             true_sigma2 = true_sigma2,
             est_sigma2  = est_sigma2,
             est_order   = est_order,
             rejected    = rejected)
}

# ---- Simulation ------------------------------------------------------------

set.seed(seed_sweep)
sweep_seeds <- sample.int(1e6, S_sweep)

sweep_results <- list()
for (noise_name in names(noise_types)) {
  cat(sprintf("\n=== %s ===\n", noise_name))
  noise_spec <- noise_types[[noise_name]]
  for (dv in d_sweep) {
    cat(sprintf("  d = %.1f: ", dv))
    reps <- lapply(sweep_seeds, function(s) run_rep_arma(dv, noise_spec, s))
    df   <- bind_rows(reps)
    df$noise_type <- noise_name
    sweep_results <- c(sweep_results, list(df))
    cat(sprintf("rej = %.2f, sigma2_ratio = %.2f, mean_order = %.1f\n",
                mean(df$rejected, na.rm = TRUE),
                mean(df$est_sigma2 / df$true_sigma2, na.rm = TRUE),
                mean(df$est_order, na.rm = TRUE)))
  }
}

sweep_results <- bind_rows(sweep_results)

# ---- Summary and plot ------------------------------------------------------

sweep_summary <- sweep_results |>
  group_by(noise_type, d) |>
  summarise(
    rejection_rate = mean(rejected, na.rm = TRUE),
    sigma2_ratio   = mean(est_sigma2 / true_sigma2, na.rm = TRUE),
    mean_order     = mean(est_order, na.rm = TRUE),
    .groups = "drop"
  )

print(sweep_summary, n = 50)

sweep_summary |>
  ggplot(aes(d, rejection_rate, colour = noise_type, group = noise_type)) +
  geom_hline(yintercept = alpha, linetype = "dashed", colour = "grey60") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  scale_y_continuous(limits = c(0, 1),
                     labels = scales::percent_format(accuracy = 1)) +
  labs(x = expression(paste(italic(d), "  (L"^2, " separation)")),
       y = "Rejection rate", colour = "Noise type",
       title = "Calibration — general AR(p) estimator across noise types",
       subtitle = sprintf("n = %d, SNR = %.1f, S = %d per d",
                          n, snr, S_sweep)) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")
