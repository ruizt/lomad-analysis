## Calibration study — visual walkthrough
##
## Shows what a single replicate looks like at each d value.
## Each section generates trends, adds noise, fits the CLT pipeline,
## and runs all three tests (estimated CLT, oracle CLT, oracle identity).

devtools::load_all()
library(ggplot2)
library(patchwork)

# ---- Shared parameters ------------------------------------------------------

n     <- 1000
phi   <- 0.5
snr   <- 1
alpha <- 0.05
seed  <- 4853

h_win <- 10
s_win <- 50

theme_strip <- theme_minimal(base_size = 10) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())

# ---- Helper: panel for one d value ------------------------------------------

make_panel <- function(d, seed) {
  trends <- sim_trends(n = n, d = d, method = "dist", seed = seed)
  sim    <- sim_noise_pair(trends, h = h_win, s = s_win,
                           lambda_target = snr, ar.coefs = phi,
                           seed = seed + 1L)

  # Estimated CLT
  fit <- lomad_fit(sim$y1, sim$y2, h = h_win, s = s_win)
  tst <- lomad_test(fit, alpha = alpha)
  vi  <- fit$valid_idx

  # Oracle CLT (true noise parameters)
  true_sigma2 <- sim$noise$series1$sigma^2
  fit_orc <- lomad_fit(sim$y1, sim$y2, h = h_win, s = s_win,
                       noise_override = list(ar = phi, sigma2 = true_sigma2))
  tst_orc <- lomad_test(fit_orc, alpha = alpha)

  # Oracle identity test (true noise parameters)
  sigma2_innov <- mean(c(sim$noise$series1$sigma, sim$noise$series2$sigma)^2)
  ident <- lomad_test_identity(sim$y1, sim$y2, q = h_win, alpha = alpha,
                               noise_override = list(ar = phi, sigma2 = sigma2_innov))

  n_rej     <- sum(tst$rejected[vi], na.rm = TRUE)
  n_rej_orc <- sum(tst_orc$rejected[fit_orc$valid_idx], na.rm = TRUE)
  n_rej_idt <- sum(ident$I, na.rm = TRUE)

  t_seq <- seq_len(n)
  df <- data.frame(t = t_seq, y1 = sim$y1, y2 = sim$y2,
                   x1 = trends$x1, x2 = trends$x2,
                   R = fit$R, rho = fit$rho,
                   Z = tst$Z, rej = tst$rejected,
                   D = ident$D)

  p_ser <- ggplot(df, aes(t)) +
    geom_line(aes(y = y1), colour = "steelblue", linewidth = 0.2, alpha = 0.35) +
    geom_line(aes(y = y2), colour = "tomato",    linewidth = 0.2, alpha = 0.35) +
    geom_line(aes(y = x1), colour = "steelblue", linewidth = 0.9) +
    geom_line(aes(y = x2), colour = "tomato",    linewidth = 0.9) +
    labs(title = sprintf("d = %.1f  |  est: %d  orc: %d  ident: %d",
                         d, n_rej, n_rej_orc, n_rej_idt),
         y = "Value") +
    theme_strip + theme(plot.title = element_text(size = 8))

  p_cor <- ggplot(df, aes(t)) +
    geom_line(aes(y = R),   colour = "grey40",  linewidth = 0.3, na.rm = TRUE) +
    geom_line(aes(y = rho), colour = "#0072B2", linewidth = 0.7, na.rm = TRUE) +
    labs(y = expression(R[t] ~ "/" ~ hat(rho)[t])) +
    theme_strip

  p_z <- ggplot(df, aes(t, Z)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_line(linewidth = 0.3, colour = "grey30", na.rm = TRUE) +
    geom_point(data = \(x) x[!is.na(x$rej) & x$rej, ],
               colour = "#D55E00", size = 1.2) +
    labs(y = expression(Z[t])) +
    theme_strip

  p_d <- ggplot(df, aes(t, D)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_line(linewidth = 0.3, colour = "grey30", na.rm = TRUE) +
    geom_point(data = \(x) x[!is.na(ident$I) & ident$I == 1L, ],
               aes(t, D), colour = "#D55E00", size = 1.2) +
    labs(y = expression(D[t]), x = "t") +
    theme_minimal(base_size = 10)

  p_ser / p_cor / p_z / p_d + plot_layout(heights = c(2, 1, 1, 1))
}

# ---- d = 0 -----------------------------------------------------------------

make_panel(d = 0, seed = seed)

# ---- d = 0.5 ---------------------------------------------------------------

make_panel(d = 0.5, seed = seed + 100L)

# ---- d = 1 -----------------------------------------------------------------

make_panel(d = 1, seed = seed + 200L)

# ---- d = 1.5 ---------------------------------------------------------------

make_panel(d = 2, seed = seed + 300L)

# ---- All together -----------------------------------------------------------

patchwork::wrap_plots(
  make_panel(0,   seed),
  make_panel(0.5, seed + 100L),
  make_panel(1,   seed + 200L),
  make_panel(1.5, seed + 300L),
  ncol = 4
)
