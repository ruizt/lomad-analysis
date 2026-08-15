## Validation study — local proof-of-concept
##
## Mirrors tide/sim.R at small S and draws draft versions of the two panels of
## fig-validation.png. tide/sim.R is the source of truth; change it first, then
## mirror the change here.
##
## See design.md for the full study specification.

library(lomad)
library(ggplot2)
library(patchwork)

S      <- 50L      # tide/sim.R runs 1000
s_win  <- 100L
seed0  <- 7291

# ---- Common DGP --------------------------------------------------------------

n_obs <- 2000
h_win <- 5

# nu_2 = B_SCALE * nu_1, so H_0 holds everywhere while tau_2^2 = B_SCALE^2 tau_1^2
B_SCALE <- 2
tr      <- sim_trends(n = n_obs, d = 0, nb = 25, sd0 = 50, p = 1.5, seed = 5381)
trend   <- tr$x1
trend2  <- B_SCALE * trend

ma_trend  <- as.numeric(stats::filter(trend,  rep(1 / h_win, h_win), sides = 1))
ma_trend2 <- as.numeric(stats::filter(trend2, rep(1 / h_win, h_win), sides = 1))

# AR(1) noise
e2e_ar1 <- 0.5;  e2e_sd1 <- 0.8
e2e_ar2 <- 0.3;  e2e_sd2 <- 0.8

.ma_filter_acov <- lomad:::.ma_filter_acov

# The window the oracle quantities are built at; must match the s the replicates
# run at, or the oracle arm is standardized at the wrong width.
e2e_s <- s_win

e2e_acov_raw1  <- arma_acov(e2e_ar1, numeric(0), e2e_sd1^2, lag_max = h_win + 100)
e2e_acov_raw2  <- arma_acov(e2e_ar2, numeric(0), e2e_sd2^2, lag_max = h_win + 100)
e2e_acov_filt1 <- .ma_filter_acov(e2e_acov_raw1, h_win, lag_max = 100)
e2e_acov_filt2 <- .ma_filter_acov(e2e_acov_raw2, h_win, lag_max = 100)
e2e_acov_filt1 <- e2e_acov_filt1[!is.na(e2e_acov_filt1)]
e2e_acov_filt2 <- e2e_acov_filt2[!is.na(e2e_acov_filt2)]
e2e_ml         <- min(length(e2e_acov_filt1), length(e2e_acov_filt2))
e2e_acov_filt1 <- e2e_acov_filt1[1:e2e_ml]
e2e_acov_filt2 <- e2e_acov_filt2[1:e2e_ml]

e2e_sigma1     <- e2e_acov_filt1[1]
e2e_sigma2     <- e2e_acov_filt2[1]
e2e_sums       <- acov_sums(e2e_acov_filt1, e2e_acov_filt2)

e2e_tau_sq     <- compute_tau_sq(ma_trend,  e2e_s)
e2e_tau2_sq    <- compute_tau_sq(ma_trend2, e2e_s)
e2e_rho_oracle <- compute_rho(e2e_tau_sq, e2e_tau2_sq, e2e_sigma1, e2e_sigma2)
e2e_V_oracle   <- compute_V(e2e_tau_sq, e2e_tau2_sq, e2e_sigma1, e2e_sigma2,
                            e2e_sums$L1, e2e_sums$L2,
                            e2e_sums$Q1, e2e_sums$Q2, e2e_sums$Q12)

e2e_eval_pts <- c(400, 900, 1400, 1900)
e2e_band_pts <- seq(h_win + e2e_s, n_obs, by = 5L)

# ---- run_rep -----------------------------------------------------------------

run_rep_e2e <- function(s, seed, eval_pts, alpha = 0.05) {
  set.seed(seed)
  z1 <- arima.sim(model = list(ar = e2e_ar1), n = n_obs, sd = e2e_sd1)
  z2 <- arima.sim(model = list(ar = e2e_ar2), n = n_obs, sd = e2e_sd2)
  y1 <- trend + z1
  y2 <- trend2 + z2

  fit <- suppressMessages(lomad_fit(y1, y2, h = h_win, s = s))
  tst <- suppressMessages(lomad_test(fit, alpha = alpha))

  list(
    eval = data.frame(
      t        = eval_pts,
      R        = fit$R[eval_pts],
      rho_hat  = fit$rho[eval_pts],
      V_hat    = fit$V[eval_pts],
      Z_pipe   = tst$Z[eval_pts],
      rejected = tst$rejected[eval_pts]
    ),
    band = data.frame(
      R      = fit$R[e2e_band_pts],
      Z_pipe = tst$Z[e2e_band_pts]
    )
  )
}

# ---- Run ---------------------------------------------------------------------

cat(sprintf("=== e2e-s%d, S = %d ===\n", s_win, S))
set.seed(seed0 + s_win)
seeds <- sample.int(1e6, S)

reps <- lapply(seeds, function(sd) run_rep_e2e(s_win, sd, e2e_eval_pts))

bR <- do.call(rbind, lapply(reps, function(x) x$band$R))
bZ <- do.call(rbind, lapply(reps, function(x) x$band$Z_pipe))

band <- do.call(rbind, lapply(seq_along(e2e_band_pts), function(j) {
  tt <- e2e_band_pts[j]
  Zo <- sqrt(s_win) * (bR[, j] - e2e_rho_oracle[tt]) / sqrt(e2e_V_oracle[tt])
  Zo <- Zo[is.finite(Zo)]
  Zp <- bZ[, j]; Zp <- Zp[is.finite(Zp)]
  f <- function(Z, ty) data.frame(
    t = tt, rho = e2e_rho_oracle[tt],
    R_mean = mean(bR[, j], na.rm = TRUE),
    R_var  = stats::var(bR[, j], na.rm = TRUE),
    V_th   = e2e_V_oracle[tt],
    t05 = mean(Z < stats::qnorm(0.05)), n = length(Z), type = ty)
  rbind(f(Zo, "Oracle"), f(Zp, "End-to-end"))
}))

# ---- Draft panels ------------------------------------------------------------

mom <- band[band$type == "Oracle", ]

p_rho <- ggplot(mom, aes(t)) +
  geom_line(aes(y = R_mean), colour = "grey30") +
  geom_line(aes(y = rho), colour = "firebrick") +
  labs(x = "Time", y = expression(rho[t]), title = "Empirical (grey) vs theoretical") +
  theme_minimal()

p_v <- ggplot(mom, aes(V_th, s_win * R_var)) +
  geom_abline(slope = 1, intercept = 0, colour = "firebrick") +
  geom_point(alpha = 0.4) +
  labs(x = expression("Theoretical" ~ V[t]), y = expression(s %.% Var(R[t]))) +
  theme_minimal()

p_cov <- ggplot(band, aes(t, t05, colour = type)) +
  geom_hline(yintercept = 0.05, linetype = "dashed", colour = "firebrick") +
  geom_line() +
  labs(x = "Time", y = "Type I error", colour = NULL) +
  theme_minimal() + theme(legend.position = "bottom")

print((p_rho | p_v) / p_cov)

cat(sprintf("\nrho: max |emp - th| = %.4f\n", max(abs(mom$R_mean - mom$rho))))
cat(sprintf("V:   median emp/th   = %.3f\n", median(s_win * mom$R_var / mom$V_th)))
for (ty in unique(band$type))
  cat(sprintf("%-11s Type I error: median %.3f\n", ty, median(band$t05[band$type == ty])))
