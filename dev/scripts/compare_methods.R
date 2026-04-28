devtools::load_all()
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)

# ---- Design notes -------------------------------------------------------
#
# SIGNAL STRUCTURE: make_trends_dist with default Fourier parameters (p=2.5,
# nb=25) concentrates power at k=1 (period = n = 1000).  Within any h=50
# window the trends are nearly flat, so the local (rolling) correlation
# between the two MA_q-smoothed series is noise-dominated at ~0.13-0.15
# regardless of d.  This is the mechanism behind the "too rigid" contrast:
#
#   lomad (rho0 = 0): flags only windows where correlation drops *below zero*.
#   With rolling correlations firmly positive at every d, lomad correctly
#   does not flag decoupling — the series remain locally co-moving even though
#   they are globally offset by d.  Type I error ≈ 0% (conservative but
#   controlled).
#
#   Full-oracle identity test: tests H0: d = 0 via the L2 norm of the
#   MA-smoothed difference D_t = MA_q(y1 - y2).  Uses the true AR(1)
#   coefficient (phi) and true innovation variance (sigma2_innov) to compute
#   Var(D_t) exactly, giving a well-calibrated ~5% level at d = 0.  Power
#   is monotone in d.  The test rejects any location/scale difference without
#   distinguishing whether the series remain locally correlated.
#
# NOISE: fixed AR(1) coefficient (phi = 0.5) so the noise spectrum is
# consistent across replicates and the oracle variance estimator is exact.

n      <- 1000
snr    <- 2.0
phi    <- 0.5
rho0   <- 0.0
q      <- 10
h      <- 50
S      <- 200
seed0  <- 7421

# d values for the power curve — stop before full-oracle saturates (~d = 4)
d_vals <- c(0, 0.5, 1, 1.5, 2, 2.5)

# ---- Full-oracle identity test ------------------------------------------
#
# Under H0 (d = 0): D_t = MA_q(e1 - e2) where e_i ~ AR(1, phi).
# Var(D_t) is computed from the exact AR(1) autocovariance and the MA_q
# filter weights.  T_stat = (mean(D^2) - Var(D)) / se, tested as N(0,1).
#
# sigma2_innov: if supplied, uses the true innovation variance (full oracle).

test_identity_oracle2 <- function(y1, y2, q, phi, sigma2_innov = NULL) {
  q  <- as.integer(q)
  u  <- y1 - y2
  D  <- as.numeric(stats::filter(u, rep(1/q, q), sides = 2))
  ok <- !is.na(D)
  Dv <- D[ok]
  nv <- length(Dv)

  ls       <- -(q - 1L):(q - 1L)
  wts_ls   <- q - abs(ls)
  acf_unit <- function(k) (1/q^2) * sum(wts_ls * phi^abs(ls + k))

  var_D_unit <- acf_unit(0)

  if (is.null(sigma2_innov)) {
    sigma2_u <- var(Dv) / var_D_unit
  } else {
    sigma2_u <- (2 * sigma2_innov) / (1 - phi^2)
  }

  var_D      <- sigma2_u * var_D_unit
  max_lag    <- min(nv - 1L, 10L * q)
  rho_sq_sum <- sum(vapply(seq_len(max_lag),
                           function(k) (acf_unit(k) / var_D_unit)^2,
                           numeric(1)))
  se_D2   <- sqrt(2 * var_D^2 / nv * (1 + 2 * rho_sq_sum))
  T_stat  <- (mean(Dv^2) - var_D) / se_D2
  p_val   <- pnorm(T_stat, lower.tail = FALSE)

  list(p_val = p_val, T_stat = T_stat, var_D = var_D, D = D)
}

# ---- Helper: run both methods on one replicate -------------------------

run_one <- function(dv, seed) {
  trends <- make_trends_dist(n = n, d = dv, seed = seed)
  sim    <- suppressMessages(
    add_noise(trends, h = h, lambda_target = snr,
              ar.coefs = c(phi), s = 2 * h, seed = seed + 1)
  )
  out_l  <- suppressMessages(
    lomad(sim$y1, sim$y2, q = q, h = h, rho0 = rho0, method = "analytic")
  )
  true_sigma2 <- mean(c(sim$noise$series1$sigma, sim$noise$series2$sigma)^2)
  out_fo <- test_identity_oracle2(sim$y1, sim$y2, q = q, phi = phi,
                                  sigma2_innov = true_sigma2)
  list(trends = trends, sim = sim, out_l = out_l, out_fo = out_fo)
}

# ---- Example plots: one column per d value ----------------------------

make_example_plot <- function(dv, seed) {
  r   <- run_one(dv, seed)
  sim <- r$sim

  df_ex <- data.frame(
    t  = seq_len(n),
    y1 = sim$y1, y2 = sim$y2,
    x1 = r$trends$x1, x2 = r$trends$x2,
    R  = r$out_l$R,
    D  = r$out_fo$D
  )

  lomad_p <- r$out_l$p_values$frac_state
  fo_p    <- r$out_fo$p_val

  p_ser <- ggplot(df_ex, aes(t)) +
    geom_line(aes(y = y1), color = "steelblue", linewidth = 0.25, alpha = 0.4) +
    geom_line(aes(y = y2), color = "tomato",    linewidth = 0.25, alpha = 0.4) +
    geom_line(aes(y = x1), color = "steelblue", linewidth = 1.0) +
    geom_line(aes(y = x2), color = "tomato",    linewidth = 1.0) +
    labs(
      title = sprintf("d = %.1f  |  lomad p = %.3f  |  full-oracle p = %.3f",
                      dv, lomad_p, fo_p),
      y = "Value"
    ) +
    theme_minimal(base_size = 10) +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
          plot.title = element_text(size = 8))

  p_r <- ggplot(df_ex, aes(t, R)) +
    geom_hline(yintercept = rho0, linetype = "dashed",
               color = "grey50", linewidth = 0.4) +
    geom_line(linewidth = 0.3, color = "grey30", na.rm = TRUE) +
    labs(y = expression(R[t])) +
    theme_minimal(base_size = 10) +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank())

  p_d <- ggplot(df_ex, aes(t, D)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_line(linewidth = 0.3, color = "grey30", na.rm = TRUE) +
    labs(y = expression(D[t]), x = "t") +
    theme_minimal(base_size = 10)

  p_ser / p_r / p_d + plot_layout(heights = c(2, 1, 1))
}

# Example plots at d = 0, 1, 2, 4
patchwork::wrap_plots(
  mapply(make_example_plot,
         c(0, 1, 2, 4),
         seed0 + c(0, 100, 200, 300),
         SIMPLIFY = FALSE),
  ncol = 4
)

# ---- Calibration check (type I error at d = 0) -------------------------

d_cal  <- 0       # H0: no separation
S_cal  <- 200
set.seed(4853)
seeds_cal <- sample.int(1e5, S_cal)

cat("\n--- Calibration check (d =", d_cal, ", S =", S_cal, "reps) ---\n")
cal_results <- vapply(seeds_cal, function(seed) {
  r <- run_one(d_cal, seed)
  c(lomad       = as.integer(r$out_l$p_values$frac_state < 0.05),
    full_oracle = as.integer(r$out_fo$p_val              < 0.05))
}, numeric(2))

cat(sprintf("  lomad        rejection rate: %.3f  (expected: near 0)\n",
            mean(cal_results["lomad", ])))
cat(sprintf("  full-oracle  rejection rate: %.3f  (expected: ~0.05)\n\n",
            mean(cal_results["full_oracle", ])))

# ---- Power curve: S reps per d value -----------------------------------

set.seed(6214)
seeds_sim <- sample.int(1e5, S)

cat("Running power curve:", S, "reps x", length(d_vals),
    "d-values (d =", paste(d_vals, collapse = ", "), ") ...\n")

sim_results <- mapply(
  function(dv, seed) {
    r <- run_one(dv, seed)
    data.frame(
      d           = dv,
      lomad       = as.integer(r$out_l$p_values$frac_state < 0.05),
      full_oracle = as.integer(r$out_fo$p_val              < 0.05)
    )
  },
  rep(d_vals, each = S),
  rep(seeds_sim, times = length(d_vals)),
  SIMPLIFY = FALSE
) |> bind_rows()

# ---- Detection rates ---------------------------------------------------

detection_rates <- sim_results |>
  group_by(d) |>
  summarise(
    lomad       = mean(lomad),
    full_oracle = mean(full_oracle),
    .groups = "drop"
  )

cat("\nRejection rates (S =", S, "reps, alpha = 0.05)\n")
cat(sprintf("  rho0 = %.1f  |  phi = %.1f  |  snr = %.1f  |  n = %d\n\n",
            rho0, phi, snr, n))
print(detection_rates)

# ---- Rejection rate plot -----------------------------------------------

detection_rates |>
  pivot_longer(c(lomad, full_oracle), names_to = "method", values_to = "rate") |>
  mutate(method = factor(
    method,
    levels = c("lomad", "full_oracle"),
    labels = c("lomad", "full-oracle identity")
  )) |>
  ggplot(aes(d, rate, color = method, group = method)) +
  geom_hline(yintercept = 0.05, linetype = "dashed",
             color = "grey60", linewidth = 0.4) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_y_continuous(limits = c(0, 1),
                     labels = scales::percent_format(accuracy = 1)) +
  scale_color_manual(values = c("lomad"                 = "steelblue",
                                "full-oracle identity"  = "tomato")) +
  labs(
    x = expression(paste(italic(d), "  (L"^2, " separation)")),
    y = "Rejection rate",
    color = NULL,
    title = sprintf(
      "lomad vs. full-oracle identity test  (n = %d, snr = %.1f, S = %d)",
      n, snr, S
    ),
    subtitle = sprintf(
      "lomad: rho0 = %.1f, h = %d, q = %d  |  noise: AR(1), phi = %.1f  |  dashed: alpha = 0.05",
      rho0, h, q, phi
    )
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")
