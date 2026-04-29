## Calibration simulation — small-scale proof of concept
##
## Three methods are compared across a range of L2 separations d:
##
##   CLT test (estimated): lomad_fit() + lomad_test() with BY-FDR.
##   CLT test (oracle): lomad_fit(noise_override = ...) + lomad_test().
##   Identity test (oracle): L2 norm of MA-smoothed difference.
##
## See settings.R for a visual walkthrough of the simulation setup.

devtools::load_all()
library(dplyr)
library(ggplot2)

# ---- Parameters ------------------------------------------------------------

n      <- 1000
phi    <- 0.5
snr    <- 1
S      <- 50          # replicates per d value
d_vals <- c(0, 0.2, 0.5, 1)
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
  # 1. CLT test: estimated pipeline
  fit <- suppressMessages(
    lomad_fit(sim$y1, sim$y2, h = h_win, s = s_win)
  )
  tst <- suppressMessages(lomad_test(fit, alpha = alpha))
  clt_rejected <- if (length(fit$valid_idx) > 0) {
    any(tst$rejected[fit$valid_idx], na.rm = TRUE)
  } else NA

  # 2. CLT test: oracle (true per-series noise parameters)
  fit_orc <- suppressMessages(
    lomad_fit(sim$y1, sim$y2, h = h_win, s = s_win,
              noise_override = list(
                list(ar = phi, sigma2 = sim$noise$series1$sigma^2),
                list(ar = phi, sigma2 = sim$noise$series2$sigma^2)
              ))
  )
  tst_orc <- lomad_test(fit_orc, alpha = alpha)
  oracle_rejected <- if (length(fit_orc$valid_idx) > 0) {
    any(tst_orc$rejected[fit_orc$valid_idx], na.rm = TRUE)
  } else NA

  # 3. Identity test (oracle, global Bonferroni p-value)
  sigma2_innov <- mean(c(sim$noise$series1$sigma,
                         sim$noise$series2$sigma)^2)
  ident <- suppressMessages(
    lomad_test_identity(sim$y1, sim$y2, q = h_win, alpha = alpha,
                        noise_override = list(ar = phi, sigma2 = sigma2_innov))
  )
  identity_rejected <- ident$global_p < alpha

  data.frame(d = d, seed = seed,
             clt_rejected      = clt_rejected,
             oracle_rejected   = oracle_rejected,
             identity_rejected = identity_rejected)
}

# testing
run_rep(0, 1)

# ---- Simulation ------------------------------------------------------------

set.seed(seed0)
seeds <- sample.int(1e6, S)

results <- vector("list", length(d_vals))

for (i in seq_along(d_vals)) {
  dv  <- d_vals[i]
  cat(sprintf("d = %.1f  ", dv))
  reps         <- lapply(seeds, function(seed) run_rep(dv, seed))
  results[[i]] <- bind_rows(reps)
  cat(sprintf("clt: %.2f  oracle: %.2f  identity: %.2f\n",
              mean(results[[i]]$clt_rejected,      na.rm = TRUE),
              mean(results[[i]]$oracle_rejected,    na.rm = TRUE),
              mean(results[[i]]$identity_rejected, na.rm = TRUE)))
}

# collate
results <- bind_rows(results)

# export
write_rds(results, file = 'dev/sims/calibration/results/template-result.rds')

# ---- Summary and plot ------------------------------------------------------

summary_tbl <- results |>
  group_by(d) |>
  summarise(clt_rate      = mean(clt_rejected,      na.rm = TRUE),
            oracle_rate   = mean(oracle_rejected,    na.rm = TRUE),
            identity_rate = mean(identity_rejected, na.rm = TRUE),
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
       y = "Rejection rate", colour = NULL,
       title  = "Calibration — CLT test vs identity test",
       subtitle = sprintf("n = %d, phi = %.1f, SNR = %.1f, S = %d per d",
                          n, phi, snr, S)) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")
