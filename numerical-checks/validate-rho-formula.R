# Numerical verification of the general rho_t formula (Eq. rho-general)
# against empirical local sample correlation R_t from simulated data.
#
# Tests two coupling methods (smooth, rate) and diagnoses the contribution
# of each simplifying assumption (drop cross-covariance, w-constant).

library(ggplot2)
library(lomad)

set.seed(4817)

# ---- Simulation parameters -------------------------------------------------

n        <- 4000
h        <- 20
s        <- 150
d        <- 15
sigma2   <- 0.01
sigma2_sm <- sigma2 / h  # smoothed noise variance

fourier_args <- list(sd0 = 10, p = 1.5, nb = 51)

# ---- Helper: compute R_t and rho approximations ----------------------------

compute_all_rho <- function(tr, y1, y2, h, s, sigma2_sm) {
  n <- length(y1)

  ma1    <- as.numeric(stats::filter(y1, rep(1/h, h), sides = 1))
  ma2    <- as.numeric(stats::filter(y2, rep(1/h, h), sides = 1))
  ma_nu1 <- as.numeric(stats::filter(tr$x1, rep(1/h, h), sides = 1))
  ma_nu2 <- as.numeric(stats::filter(tr$x2, rep(1/h, h), sides = 1))

  mu_bar_sm <- as.numeric(stats::filter((tr$x1 + tr$x2) / 2,
                                         rep(1/h, h), sides = 1))
  ma_hd     <- as.numeric(stats::filter((tr$x1 - tr$x2) / 2,
                                         rep(1/h, h), sides = 1))

  R_emp        <- rep(NA_real_, n)
  rho_direct   <- rep(NA_real_, n)
  rho_general  <- rep(NA_real_, n)
  rho_no_cross <- rep(NA_real_, n)

  for (t in s:n) {
    idx <- (t - s + 1):t
    m1 <- ma1[idx]; m2 <- ma2[idx]
    v1 <- ma_nu1[idx]; v2 <- ma_nu2[idx]
    mb <- mu_bar_sm[idx]; hd <- ma_hd[idx]
    if (any(is.na(m1)) || any(is.na(m2)) ||
        any(is.na(v1)) || any(is.na(v2)) ||
        any(is.na(mb)) || any(is.na(hd))) next

    # Empirical R_t
    R_emp[t] <- cor(m1, m2)

    # Direct oracle rho (exact window moments + known noise)
    cov12 <- mean((v1 - mean(v1)) * (v2 - mean(v2)))
    var1  <- mean((v1 - mean(v1))^2) + sigma2_sm
    var2  <- mean((v2 - mean(v2))^2) + sigma2_sm
    rho_direct[t] <- cov12 / sqrt(var1 * var2)

    var_mb  <- mean((mb - mean(mb))^2)
    var_hd  <- mean((hd - mean(hd))^2)
    cov_mh  <- mean((mb - mean(mb)) * (hd - mean(hd)))

    # General formula: drop cross-covariance, use exact variances
    rho_general[t] <- (var_mb - var_hd) / (var_mb + var_hd + sigma2_sm)

    # Diagnostic: drop cross-cov only (exact variances, equal noise)
    rho_no_cross[t] <- (var_mb - var_hd) /
                        (var_mb + var_hd + sigma2_sm)
  }

  list(R_emp = R_emp, rho_direct = rho_direct,
       rho_general = rho_general, rho_no_cross = rho_no_cross)
}

# sim_trends() no longer normalises to ||x1 - x2|| = d -- d scales the distinct
# component, and separation is linear in it. This check compares rho formulas at
# a fixed separation, so the multiplier that reaches `d` is solved for here
# rather than assumed. It differs by structure, which is exactly what the old
# normalisation was hiding.
at_sep <- function(args) {
  probe <- do.call(sim_trends, c(list(n = n, d = 1), args, fourier_args))
  scale <- d / sqrt(sum((probe$x1 - probe$x2)^2))
  do.call(sim_trends, c(list(n = n, d = scale), args, fourier_args))
}

report <- function(label, x, y) {
  valid <- which(!is.na(x) & !is.na(y))
  cat(sprintf("  %-45s cor = %.3f  RMSE = %.4f\n", label,
              cor(x[valid], y[valid]),
              sqrt(mean((x[valid] - y[valid])^2))))
}

# ---- Smooth coupling -------------------------------------------------------

cat("=== Random separation (bw = 150, coupling = 0.6) ===\n")
tr_sm <- at_sep(list(method = "rs", bw = 150, coupling = 0.6))

y1_sm <- tr_sm$x1 + rnorm(n, sd = sqrt(sigma2))
y2_sm <- tr_sm$x2 + rnorm(n, sd = sqrt(sigma2))

res_sm <- compute_all_rho(tr_sm, y1_sm, y2_sm, h, s, sigma2_sm)

report("Direct (oracle) vs R_emp",     res_sm$rho_direct,  res_sm$R_emp)
report("General formula vs R_emp",     res_sm$rho_general, res_sm$R_emp)
report("General formula vs Direct",    res_sm$rho_general, res_sm$rho_direct)

# ---- Rate coupling ----------------------------------------------------------

cat("\n=== Fixed rate (rate = 0.003) ===\n")
set.seed(7213)
tr_rt <- at_sep(list(method = "fr", rate = 0.003))

y1_rt <- tr_rt$x1 + rnorm(n, sd = sqrt(sigma2))
y2_rt <- tr_rt$x2 + rnorm(n, sd = sqrt(sigma2))

res_rt <- compute_all_rho(tr_rt, y1_rt, y2_rt, h, s, sigma2_sm)

report("Direct (oracle) vs R_emp",     res_rt$rho_direct,  res_rt$R_emp)
report("General formula vs R_emp",     res_rt$rho_general, res_rt$R_emp)
report("General formula vs Direct",    res_rt$rho_general, res_rt$rho_direct)

# ---- Diagnostic: w-constant approximation -----------------------------------

cat("\n=== Diagnostic: w-constant approximation (smooth coupling) ===\n")

# Recover base Delta (before coupling modulation)
w_safe <- pmax(tr_sm$w, 0.001)
delta_base <- (tr_sm$x1 - tr_sm$x2) / (1 - w_safe)
ma_db <- as.numeric(stats::filter(delta_base, rep(1/h, h), sides = 1))
mu_bar_sm <- as.numeric(stats::filter((tr_sm$x1 + tr_sm$x2) / 2,
                                       rep(1/h, h), sides = 1))

rho_w_param <- rep(NA_real_, n)
for (t in s:n) {
  idx <- (t - s + 1):t
  mb <- mu_bar_sm[idx]; db <- ma_db[idx]; ww <- tr_sm$w[idx]
  if (any(is.na(mb)) || any(is.na(db))) next
  tau_mb2 <- mean((mb - mean(mb))^2)
  tau_db2 <- mean((db - mean(db))^2)
  w_avg   <- mean(ww)
  lam <- tau_mb2 / sigma2_sm
  kap <- tau_db2 / (4 * sigma2_sm)
  rho_w_param[t] <- (lam - (1 - w_avg)^2 * kap) /
                     (lam + (1 - w_avg)^2 * kap + 1)
}

report("w-parametrized formula vs R_emp",  rho_w_param, res_sm$R_emp)
report("w-parametrized formula vs Direct", rho_w_param, res_sm$rho_direct)

# ---- Plots ------------------------------------------------------------------

make_plot <- function(res, title) {
  valid <- which(!is.na(res$R_emp) & !is.na(res$rho_general))
  df <- data.frame(t = valid, R_emp = res$R_emp[valid],
                   rho = res$rho_general[valid])
  ggplot(df, aes(x = t)) +
    geom_line(aes(y = R_emp), alpha = 0.4, linewidth = 0.3) +
    geom_line(aes(y = rho), color = "firebrick", linewidth = 0.6) +
    labs(x = "Time", y = expression(rho[t]), title = title) +
    theme_minimal()
}

p1 <- make_plot(res_sm, "Smooth coupling: R (grey) vs formula (red)")
p2 <- make_plot(res_rt, "Rate coupling: R (grey) vs formula (red)")
