## Calibration simulation — small-scale proof of concept
##
## Three methods are compared across a range of L2 separations d:
##
##   CLT test (estimated): lomad_fit() + lomad_test() with BY-FDR.
##   CLT test (oracle): lomad_fit(noise_override = ...) + lomad_test().
##   Identity test (oracle): L2 norm of MA-smoothed difference.
##
## See templates.R for a visual walkthrough of the simulation setup.

devtools::load_all()
library(dplyr)
library(ggplot2)

# ---- Parameters ------------------------------------------------------------

n      <- 1000
phi    <- 0.5
snr    <- 1
S      <- 50          # replicates per d value
d_vals <- c(0, 0.2, 0.5, 1, 1.5)
alpha  <- 0.05
seed0  <- 4853

h_win <- 10
s_win <- 50

# ---- Single-replicate function ---------------------------------------------

run_rep <- function(d, seed) {
  trends <- sim_trends(n = n, d = d, method = "dist", seed = seed)
  sim    <- suppressMessages(
    sim_noise_pair(trends, h = h_win, s = s_win, lambda_target = snr,
                   ar.coefs = phi, seed = seed + 1L)
  )
  true_sigma2 <- sim$noise$series1$sigma^2

  # 1. CLT test: estimated pipeline
  fit <- suppressMessages(
    lomad_fit(sim$y1, sim$y2, h = h_win, s = s_win)
  )
  tst <- suppressMessages(lomad_test(fit, alpha = alpha))
  clt_frac <- if (length(fit$valid_idx) > 0) {
    mean(tst$rejected[fit$valid_idx], na.rm = TRUE)
  } else NA_real_

  # 2. CLT test: oracle (true noise parameters)
  fit_orc <- suppressMessages(
    lomad_fit(sim$y1, sim$y2, h = h_win, s = s_win,
              noise_override = list(ar = phi, sigma2 = true_sigma2))
  )
  tst_orc <- lomad_test(fit_orc, alpha = alpha)
  oracle_frac <- if (length(fit_orc$valid_idx) > 0) {
    mean(tst_orc$rejected[fit_orc$valid_idx], na.rm = TRUE)
  } else NA_real_

  # 3. Identity test (oracle contrast)
  sigma2_innov <- mean(c(sim$noise$series1$sigma,
                         sim$noise$series2$sigma)^2)
  ident <- suppressMessages(
    lomad_test_identity(sim$y1, sim$y2, q = h_win, alpha = alpha,
                        noise_override = list(ar = phi, sigma2 = sigma2_innov))
  )
  identity_frac <- mean(ident$I, na.rm = TRUE)

  data.frame(d = d, seed = seed,
             clt_frac      = clt_frac,
             oracle_frac   = oracle_frac,
             identity_frac = identity_frac)
}

# ---- Simulation ------------------------------------------------------------

set.seed(seed0)
seeds <- sample.int(1e6, S)

results <- vector("list", length(d_vals))

for (i in seq_along(d_vals)) {
  dv  <- d_vals[i]
  cat(sprintf("d = %.1f  ", dv))
  reps         <- lapply(seeds, function(seed) run_rep(dv, seed))
  results[[i]] <- bind_rows(reps)
  cat(sprintf("clt: %.3f  oracle: %.3f  identity: %.3f\n",
              mean(results[[i]]$clt_frac,      na.rm = TRUE),
              mean(results[[i]]$oracle_frac,   na.rm = TRUE),
              mean(results[[i]]$identity_frac, na.rm = TRUE)))
}

results <- bind_rows(results)

# ---- Summary and plot ------------------------------------------------------

summary_tbl <- results |>
  group_by(d) |>
  summarise(clt_rate      = mean(clt_frac,      na.rm = TRUE),
            oracle_rate   = mean(oracle_frac,   na.rm = TRUE),
            identity_rate = mean(identity_frac, na.rm = TRUE),
            .groups = "drop")

print(summary_tbl)

summary_tbl |>
  tidyr::pivot_longer(c(clt_rate, oracle_rate, identity_rate),
                      names_to = "method", values_to = "rate") |>
  mutate(method = factor(method,
                         levels = c("clt_rate", "oracle_rate", "identity_rate"),
                         labels = c("CLT (estimated)", "CLT (oracle)",
                                    "Identity (oracle)"))) |>
  ggplot(aes(d, rate, colour = method, group = method)) +
  geom_hline(yintercept = alpha, linetype = "dashed", colour = "grey60") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  scale_y_continuous(limits = c(0, 1),
                     labels = scales::percent_format(accuracy = 1)) +
  scale_colour_manual(values = c("CLT (estimated)"  = "#0072B2",
                                 "CLT (oracle)"      = "#009E73",
                                 "Identity (oracle)" = "#D55E00")) +
  labs(x = expression(paste(italic(d), "  (L"^2, " separation)")),
       y = "Fraction rejected", colour = NULL,
       title  = "Calibration — fraction of time points rejected",
       subtitle = sprintf("n = %d, phi = %.1f, SNR = %.1f, S = %d per d",
                          n, phi, snr, S)) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")
