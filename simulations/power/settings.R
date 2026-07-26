## settings.R — visual walkthrough of power study replicates
##
## One section per trend structure, each with d adjustable at the top.
## For each structure, three panels are produced:
##   1. Example dataset (observed series + true trends)
##   2. Observed vs expected rolling correlation
##   3. Z statistic with rejected points highlighted
##
## Usage (from repo root):
##   source("simulations/power/settings.R")
##   — or step through one structure at a time interactively.

library(lomad)
library(ggplot2)
library(patchwork)

# ---- Shared parameters -------------------------------------------------------

n     <- 1000L # 250, 500, 1000
phi   <- 0.5 # 0.3, 0.8
snr   <- 0.5 # hi/low (TBD)
alpha <- 0.05 # stays fixed
seed  <- 6183 # stays fixed

h_win <- 5 # seemed reasonable earlier
s_win <- 50 # seemed reasonable earlier

theme_strip <- theme_minimal(base_size = 10) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())

# ---- Helper: panel for one (structure, d) ------------------------------------

make_panel <- function(structure, d, extra_params = list()) {
  set.seed(seed)

  trends <- do.call(sim_trends,
    c(list(n = n, d = d, method = structure, seed = seed), extra_params))

  sim <- sim_noise_pair(trends, h = h_win, lambda_target = snr,
                        ar.coefs = phi, seed = seed + 1L)

  fit <- lomad_fit(sim$y1, sim$y2, h = h_win, s = s_win)
  tst <- lomad_test(fit, alpha = alpha)
  vi  <- fit$valid_idx

  n_rej <- sum(tst$rejected, na.rm = TRUE)

  df <- data.frame(
    t   = seq_len(n),
    y1  = sim$y1, y2 = sim$y2,
    f1  = sim$x1, f2 = sim$x2,
    w   = trends$w,
    R   = fit$R,  rho = fit$rho,
    Z   = rep(NA_real_, n),
    rej = rep(NA, n)
  )
  df$Z[vi]   <- tst$Z
  df$rej[vi] <- tst$rejected

  title_str <- sprintf("%s, d = %.1f  |  %d / %d rejected",
                       structure, d, n_rej, length(vi))

  p_ser <- ggplot(df, aes(x = t)) +
    geom_line(aes(y = y1), colour = "steelblue", linewidth = 0.2, alpha = 0.35) +
    geom_line(aes(y = y2), colour = "tomato",    linewidth = 0.2, alpha = 0.35) +
    geom_line(aes(y = f1), colour = "steelblue", linewidth = 0.9) +
    geom_line(aes(y = f2), colour = "tomato",    linewidth = 0.9) +
    labs(title = title_str, y = "value") +
    theme_strip + theme(plot.title = element_text(size = 8))

  p_cor <- ggplot(df, aes(x = t)) +
    geom_line(aes(y = R),   colour = "grey40",  linewidth = 0.3, na.rm = TRUE) +
    geom_line(aes(y = rho), colour = "#0072B2", linewidth = 0.7, na.rm = TRUE) +
    labs(y = expression(R[t] ~ "/" ~ hat(rho)[t])) +
    theme_strip

  p_z <- ggplot(df, aes(x = t, y = Z)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_line(linewidth = 0.3, colour = "grey30", na.rm = TRUE) +
    geom_point(data = \(x) x[!is.na(x$rej) & x$rej, ],
               colour = "#D55E00", size = 1.2) +
    geom_hline(yintercept = qnorm(1 - tst$alpha_eff),
               linetype = "dashed", colour = "steelblue") +
    labs(y = expression(Z[t])) +
    theme_strip

  p_w <- ggplot(df, aes(x = t, y = w)) +
    geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey60") +
    geom_line(linewidth = 0.5, colour = "grey30") +
    labs(y = expression(w[t]), x = "t") +
    theme_minimal(base_size = 10)

  p_ser / p_cor / p_z / p_w + plot_layout(heights = c(2, 1, 1, 1))
}

# ---- smooth ------------------------------------------------------------------

d <- 1.5
make_panel("smooth", d = d, extra_params = list(bw = 50))

# ---- cross -------------------------------------------------------------------

d <- 1.5
make_panel("cross", d = d, extra_params = list(bw = 50))

# ---- rate --------------------------------------------------------------------

d <- 1.5
make_panel("rate", d = d, extra_params = list(rate = 0.005))
