## Validation study — local proof-of-concept
##
## Runs all seven experiments at small scale and produces draft versions of the
## three paper figures. The run_rep_*() functions defined here are the same ones
## used in tide/sim.R.
##
## See design.md for the full study specification.

library(lomad)
library(ggplot2)
library(patchwork)

# ---- Common DGP --------------------------------------------------------------

n_obs <- 2000
h_win <- 20

# Shared trend via sim_trends (common trend, d = 0)
tr     <- sim_trends(n = n_obs, d = 0, nb = 25, sd0 = 50, p = 1.5, seed = 5381)
trend  <- tr$x1

# Noiseless MA-smoothed trend — used for oracle tau_sq everywhere.
# Because d = 0 (shared trend), this equals filter(trend, ...) directly;
# no bias correction is needed.
ma_trend <- as.numeric(stats::filter(trend, rep(1 / h_win, h_win), sides = 1))

# ARMA(1,1) noise parameters (oracle experiments, Figures 1–2)
oracle_ar1 <- 0.6;  oracle_ma1 <- 0.3;   oracle_sd1 <- 0.8
oracle_ar2 <- 0.4;  oracle_ma2 <- -0.2;  oracle_sd2 <- 1.0

# AR(1) noise parameters (end-to-end, Figure 3)
e2e_ar1 <- 0.5;  e2e_sd1 <- 0.8
e2e_ar2 <- 0.3;  e2e_sd2 <- 0.8
e2e_s   <- 150

# Precompute oracle filtered autocovariances
.ma_filter_acov <- lomad:::.ma_filter_acov

acov_raw1  <- arma_acov(oracle_ar1, oracle_ma1, oracle_sd1^2, lag_max = 200)
acov_raw2  <- arma_acov(oracle_ar2, oracle_ma2, oracle_sd2^2, lag_max = 200)
acov_eta1  <- .ma_filter_acov(acov_raw1, h_win, lag_max = 200)
acov_eta2  <- .ma_filter_acov(acov_raw2, h_win, lag_max = 200)
sigma1_sq  <- acov_eta1[1]
sigma2_sq  <- acov_eta2[1]
acov_eta1  <- acov_eta1[!is.na(acov_eta1)]
acov_eta2  <- acov_eta2[!is.na(acov_eta2)]
ml         <- min(length(acov_eta1), length(acov_eta2))
acov_eta1  <- acov_eta1[1:ml]
acov_eta2  <- acov_eta2[1:ml]
cov_sums   <- acov_sums(acov_eta1, acov_eta2)

# E2E oracle quantities
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
# Use noiseless ma_trend directly — no bias correction needed in oracle setting
e2e_tau_sq     <- compute_tau_sq(ma_trend, e2e_s)
e2e_rho_oracle <- compute_rho(e2e_tau_sq, e2e_sigma1, e2e_sigma2)
e2e_V_oracle   <- compute_V(e2e_tau_sq, e2e_sigma1, e2e_sigma2,
                              e2e_sums$L1, e2e_sums$L2,
                              e2e_sums$Q1, e2e_sums$Q2, e2e_sums$Q12)
e2e_eval_pts   <- c(500, 850, 1000, 1400, 1600)

# ==============================================================================
# run_rep functions
# ==============================================================================

# ---- CLT experiment (Figure 1) -----------------------------------------------

run_rep_clt <- function(s, seed, eval_t = 600) {
  set.seed(seed)
  z1 <- arima.sim(model = list(ar = oracle_ar1, ma = oracle_ma1),
                  n = n_obs, sd = oracle_sd1)
  z2 <- arima.sim(model = list(ar = oracle_ar2, ma = oracle_ma2),
                  n = n_obs, sd = oracle_sd2)
  y1 <- trend + z1
  y2 <- trend + z2
  m1 <- as.numeric(stats::filter(y1, rep(1 / h_win, h_win), sides = 1))
  m2 <- as.numeric(stats::filter(y2, rep(1 / h_win, h_win), sides = 1))
  idx <- (eval_t - s + 1):eval_t
  if (any(is.na(m1[idx])) || any(is.na(m2[idx]))) return(NA_real_)
  cor(m1[idx], m2[idx])
}

# ---- Rho experiment (Figure 2A) ----------------------------------------------

run_rep_rho <- function(s, seed) {
  set.seed(seed)
  z1 <- arima.sim(model = list(ar = oracle_ar1, ma = oracle_ma1),
                  n = n_obs, sd = oracle_sd1)
  z2 <- arima.sim(model = list(ar = oracle_ar2, ma = oracle_ma2),
                  n = n_obs, sd = oracle_sd2)
  y1 <- trend + z1
  y2 <- trend + z2
  m1 <- as.numeric(stats::filter(y1, rep(1 / h_win, h_win), sides = 1))
  m2 <- as.numeric(stats::filter(y2, rep(1 / h_win, h_win), sides = 1))

  R <- rep(NA_real_, n_obs)
  for (t in s:n_obs) {
    idx <- (t - s + 1):t
    a <- m1[idx]; b <- m2[idx]
    if (any(is.na(a)) || any(is.na(b))) next
    R[t] <- cor(a, b)
  }
  R
}

# ---- Var experiment (Figure 2B) ----------------------------------------------

run_rep_var <- function(s, seed, eval_pts) {
  set.seed(seed)
  z1 <- arima.sim(model = list(ar = oracle_ar1, ma = oracle_ma1),
                  n = n_obs, sd = oracle_sd1)
  z2 <- arima.sim(model = list(ar = oracle_ar2, ma = oracle_ma2),
                  n = n_obs, sd = oracle_sd2)
  y1 <- trend + z1
  y2 <- trend + z2
  m1 <- as.numeric(stats::filter(y1, rep(1 / h_win, h_win), sides = 1))
  m2 <- as.numeric(stats::filter(y2, rep(1 / h_win, h_win), sides = 1))

  R <- numeric(length(eval_pts))
  for (j in seq_along(eval_pts)) {
    t0  <- eval_pts[j]
    idx <- (t0 - s + 1):t0
    a <- m1[idx]; b <- m2[idx]
    R[j] <- if (any(is.na(a)) || any(is.na(b))) NA_real_ else cor(a, b)
  }
  R
}

# ---- End-to-end experiment (Figure 3) ----------------------------------------

run_rep_e2e <- function(s, seed, eval_pts, alpha = 0.05) {
  set.seed(seed)
  z1 <- arima.sim(model = list(ar = e2e_ar1), n = n_obs, sd = e2e_sd1)
  z2 <- arima.sim(model = list(ar = e2e_ar2), n = n_obs, sd = e2e_sd2)
  y1 <- trend + z1
  y2 <- trend + z2

  fit <- suppressMessages(lomad_fit(y1, y2, method = "clt", h = h_win, s = s))
  tst <- suppressMessages(lomad_test(fit, alpha = alpha))

  data.frame(
    t        = eval_pts,
    R        = fit$R[eval_pts],
    rho_hat  = fit$rho[eval_pts],
    V_hat    = fit$V[eval_pts],
    Z_pipe   = tst$Z[eval_pts],
    rejected = tst$rejected[eval_pts]
  )
}

# ==============================================================================
# Run experiments at small scale
# ==============================================================================

seed0 <- 7291

# ---- Figure 1: CLT QQ -------------------------------------------------------

cat("=== Figure 1: CLT QQ ===\n")

S_clt   <- 200
s_vals  <- c(80, 150, 300)
eval_t  <- 600

fig1_data <- list()

for (s in s_vals) {
  cat(sprintf("  s = %d ...\n", s))
  set.seed(seed0 + s)
  seeds <- sample.int(1e6, S_clt)

  R_vec <- vapply(seeds, function(sd) run_rep_clt(s, sd, eval_t), numeric(1))

  tau_sq <- compute_tau_sq(ma_trend, s)
  rho_t  <- compute_rho(tau_sq, sigma1_sq, sigma2_sq)
  V_t    <- compute_V(tau_sq, sigma1_sq, sigma2_sq,
                       cov_sums$L1, cov_sums$L2,
                       cov_sums$Q1, cov_sums$Q2, cov_sums$Q12)

  rho_0 <- rho_t[eval_t]
  V_0   <- V_t[eval_t]
  Z     <- sqrt(s) * (R_vec - rho_0) / sqrt(V_0)
  Z     <- Z[!is.na(Z)]

  cov95 <- mean(abs(Z) < qnorm(0.975))

  fig1_data[[as.character(s)]] <- list(Z = Z, s = s, cov95 = cov95)
  cat(sprintf("    rho=%.3f  mean(Z)=%.3f  sd(Z)=%.3f  cov95=%.3f\n",
              rho_0, mean(Z), sd(Z), cov95))
}

# Draft Figure 1
qq_list <- lapply(fig1_data, function(d) {
  nn <- length(d$Z)
  data.frame(
    theoretical = qnorm(ppoints(nn)),
    empirical   = sort(d$Z),
    s_label     = sprintf("s = %d  (95%% cov = %.2f)", d$s, d$cov95),
    s_num       = d$s
  )
})
qq_df <- do.call(rbind, qq_list)
qq_df$s_label <- reorder(qq_df$s_label, qq_df$s_num)

ggplot(qq_df, aes(theoretical, empirical)) +
  geom_abline(slope = 1, intercept = 0, colour = "firebrick", linewidth = 0.5) +
  geom_point(alpha = 0.3, size = 0.7) +
  facet_wrap(~s_label, nrow = 1) +
  labs(x = "Theoretical N(0,1)", y = "Empirical Z") +
  theme_bw(base_size = 10)

# ---- Figure 2A: rho accuracy ------------------------------------------------

cat("=== Figure 2A: rho ===\n")

S_rho  <- 50
s_rho  <- c(80, 250)

fig2a_data <- list()

for (s in s_rho) {
  cat(sprintf("  s = %d ...\n", s))
  set.seed(seed0 + s + 1000)
  seeds <- sample.int(1e6, S_rho)

  R_accum <- rep(0, n_obs)
  R_count <- rep(0L, n_obs)

  for (sd in seeds) {
    R_rep <- run_rep_rho(s, sd)
    ok <- !is.na(R_rep)
    R_accum[ok] <- R_accum[ok] + R_rep[ok]
    R_count[ok] <- R_count[ok] + 1L
  }

  R_mean <- ifelse(R_count > 0, R_accum / R_count, NA_real_)
  tau_sq <- compute_tau_sq(ma_trend, s)
  rho_th <- compute_rho(tau_sq, sigma1_sq, sigma2_sq)

  fig2a_data[[as.character(s)]] <- list(
    R_mean = R_mean, rho_th = rho_th, s = s
  )
}

# Draft Figure 2A
rho_list <- lapply(fig2a_data, function(d) {
  valid <- which(!is.na(d$R_mean) & !is.na(d$rho_th))
  data.frame(
    t      = valid,
    R_mean = d$R_mean[valid],
    rho_th = d$rho_th[valid],
    s      = sprintf("s = %d", d$s)
  )
})
rho_df <- do.call(rbind, rho_list)

p2a <- ggplot(rho_df, aes(t)) +
  geom_line(aes(y = rho_th, colour = "Theoretical"), linewidth = 0.6) +
  geom_line(aes(y = R_mean, colour = "Empirical"), linewidth = 0.4, alpha = 0.7) +
  facet_wrap(~s, nrow = 1) +
  scale_colour_manual(values = c("Theoretical" = "firebrick",
                                  "Empirical" = "grey30")) +
  labs(y = expression(rho[t]), colour = NULL, x = "Time") +
  theme_bw(base_size = 10) +
  theme(legend.position = "top")

# ---- Figure 2B: V accuracy --------------------------------------------------

cat("=== Figure 2B: V ===\n")

S_var     <- 200
s_var     <- 200
var_grid  <- seq(s_var + h_win, n_obs, by = 20)

set.seed(seed0 + 2000)
seeds_var <- sample.int(1e6, S_var)

R_mat_var <- matrix(NA_real_, S_var, length(var_grid))
for (i in seq_along(seeds_var)) {
  R_mat_var[i, ] <- run_rep_var(s_var, seeds_var[i], var_grid)
}

V_emp    <- s_var * apply(R_mat_var, 2, var, na.rm = TRUE)
tau_sq_v <- compute_tau_sq(ma_trend, s_var)
V_th     <- compute_V(tau_sq_v, sigma1_sq, sigma2_sq,
                       cov_sums$L1, cov_sums$L2,
                       cov_sums$Q1, cov_sums$Q2, cov_sums$Q12)
V_theory <- V_th[var_grid]

valid_v  <- which(!is.na(V_emp) & !is.na(V_theory) & V_theory > 0)
cor_V    <- cor(V_emp[valid_v], V_theory[valid_v])

p2b <- ggplot(
  data.frame(V_theory = V_theory[valid_v], V_emp = V_emp[valid_v]),
  aes(V_theory, V_emp)
) +
  geom_abline(slope = 1, intercept = 0, colour = "firebrick", linewidth = 0.5) +
  geom_point(alpha = 0.4, size = 1) +
  annotate("text", x = min(V_theory[valid_v]), y = max(V_emp[valid_v]),
           label = sprintf("r = %.3f", cor_V), hjust = 0, vjust = 1, size = 3) +
  labs(x = expression("Theoretical " * V[t]),
       y = expression("Empirical " * s %.% Var(R[t]))) +
  theme_bw(base_size = 10)

# Combine Figure 2
p2a / p2b + plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")

# ---- Figure 3: End-to-end pipeline ------------------------------------------

cat("=== Figure 3: end-to-end ===\n")

S_e2e  <- 50
alpha  <- 0.05

set.seed(seed0 + 3000)
seeds_e2e <- sample.int(1e6, S_e2e)

e2e_results <- lapply(seeds_e2e, function(sd) {
  run_rep_e2e(e2e_s, sd, e2e_eval_pts, alpha)
})

# Stack into matrices
n_pts     <- length(e2e_eval_pts)
R_mat_e   <- matrix(NA_real_, S_e2e, n_pts)
Z_mat_e   <- matrix(NA_real_, S_e2e, n_pts)

for (i in seq_along(e2e_results)) {
  R_mat_e[i, ]   <- e2e_results[[i]]$R
  Z_mat_e[i, ]   <- e2e_results[[i]]$Z_pipe
}

# Panel A: QQ at t = 1000
j_mid   <- which(e2e_eval_pts == 1000)
R_j     <- R_mat_e[, j_mid]
rho_or  <- e2e_rho_oracle[e2e_eval_pts[j_mid]]
V_or    <- e2e_V_oracle[e2e_eval_pts[j_mid]]
Z_oracle <- sqrt(e2e_s) * (R_j - rho_or) / sqrt(V_or)
Z_oracle <- Z_oracle[!is.na(Z_oracle)]
Z_pipe   <- Z_mat_e[, j_mid]
Z_pipe   <- Z_pipe[!is.na(Z_pipe)]

nn_e  <- min(length(Z_oracle), length(Z_pipe))
qq_e  <- data.frame(
  theoretical = rep(qnorm(ppoints(nn_e)), 2),
  z           = c(sort(Z_oracle[1:nn_e]), sort(Z_pipe[1:nn_e])),
  type        = rep(c("Oracle", "Pipeline"), each = nn_e)
)

p3a <- ggplot(qq_e, aes(theoretical, z, colour = type)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey40", linewidth = 0.5) +
  geom_point(alpha = 0.4, size = 0.8) +
  scale_colour_manual(values = c("Oracle" = "#0072B2", "Pipeline" = "#D55E00")) +
  labs(x = "Theoretical N(0,1)", y = "Empirical Z", colour = NULL) +
  theme_bw(base_size = 10) +
  theme(legend.position = "top")

# Panel B: coverage at each eval point
cov_rows <- vector("list", 2 * n_pts)
for (j in seq_len(n_pts)) {
  rho_j <- e2e_rho_oracle[e2e_eval_pts[j]]
  V_j   <- e2e_V_oracle[e2e_eval_pts[j]]

  Zo <- sqrt(e2e_s) * (R_mat_e[, j] - rho_j) / sqrt(V_j)
  Zo <- Zo[!is.na(Zo)]
  c_o <- mean(abs(Zo) < qnorm(0.975))
  n_o <- length(Zo)

  Zp <- Z_mat_e[, j]
  Zp <- Zp[!is.na(Zp)]
  c_p <- mean(abs(Zp) < qnorm(0.975))
  n_p <- length(Zp)

  cov_rows[[j]]         <- data.frame(
    t = e2e_eval_pts[j], cov = c_o,
    se = sqrt(c_o * (1 - c_o) / n_o), type = "Oracle"
  )
  cov_rows[[j + n_pts]] <- data.frame(
    t = e2e_eval_pts[j], cov = c_p,
    se = sqrt(c_p * (1 - c_p) / n_p), type = "Pipeline"
  )
}
cov_df <- do.call(rbind, cov_rows)

p3b <- ggplot(cov_df, aes(factor(t), cov, colour = type)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", colour = "grey60") +
  geom_point(position = position_dodge(width = 0.4), size = 2) +
  geom_errorbar(aes(ymin = cov - 2 * se, ymax = cov + 2 * se),
                position = position_dodge(width = 0.4), width = 0.2) +
  scale_colour_manual(values = c("Oracle" = "#0072B2", "Pipeline" = "#D55E00")) +
  labs(x = "Time point", y = "95% coverage", colour = NULL) +
  theme_bw(base_size = 10) +
  theme(legend.position = "top")

# Combine Figure 3
p3a + p3b +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")

cat("\nDone. All draft figures should be visible in the Plots pane.\n")
