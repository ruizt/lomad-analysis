# Approximation Quality Study Design

## Objective

Assess the empirical accuracy of the Proposition 1 approximations for the
local population correlation ρ and asymptotic variance V under the
signal-plus-noise model. The goal is to determine when the approximations
support the methodology and where they degrade, not to prove new theory.

## Background

Proposition 1 gives closed-form expressions for ρ and V in terms of the
signal variance τ² and the smoothed noise ACVF γ(l). These expressions are
derived under the assumption that the shared trend is nearly constant over
each local window. The study asks how well this approximation holds in
practice across a range of conditions, and whether any material degradation
occurs in the parameter regimes used by the paper.

## Quantities to compare

For each simulation condition, compute:

- **Empirical mean of Rₜ** — average of R_t across replicates at each *t*
- **Empirical variance of Rₜ** — variance of R_t across replicates at each *t*  
- **Theoretical ρₜ** — from Proposition 1 using true τ², σ₁², σ₂²
- **Theoretical V/s** — from Proposition 1 using true ACVF and τ²

Then report:
- Absolute error in the mean: |E[Rₜ] − ρₜ|
- Relative error in the mean: |E[Rₜ] − ρₜ| / ρₜ
- Absolute error in the variance: |Var(Rₜ) − Vₜ/s|
- Relative error in the variance: |Var(Rₜ) − Vₜ/s| / (Vₜ/s)

Summarise errors pooled over valid time points and across replicates.

## Simulation design

Use a focused factorial over the following factors. This is not an exhaustive
grid — choose representative combinations that span the range of conditions
relevant to the paper.

### Core factors

| Factor | Levels |
|--------|--------|
| Series length *T* | 500, 1000 |
| Smoothing window *h* | small (5), medium (10), large (25) |
| Correlation window *s* | small (50), medium (150), large (250) |
| Signal-to-noise ratio *λ* | 0.5, 1.5, 3.0 |
| Trend mismatch *δ* | 0 (shared), 0.25, 0.5 (via `d` in `make_trends_dist`) |
| Trend smoothness | slow (default Fourier basis), fast (higher-frequency basis) |

AR(1) noise with φ = 0.5 throughout. Use *R* = 500 replicates per condition
so that Monte Carlo error in the empirical mean and variance is clearly
smaller than any approximation errors of interest.

### Interpretation of *δ*

The approximation for ρ is derived under a shared-trend assumption. Setting
δ > 0 (via *d* > 0 in `make_trends_dist`) introduces mismatch between the
two trends, testing how far the approximation extends beyond the exact
shared-trend case. Small δ (d ≤ 0.25) is the regime where the CLT test is
used in practice; moderate δ (d ≈ 0.5) is near the edge of the null region.

---

## Available functions

All functions below are available after `devtools::load_all()`.

### Data generation

```r
# Generate trend pair with L² separation d
trends <- make_trends_dist(n, d, seed)
# Returns: list(x1, x2)  — trend vectors of length n

# Add calibrated AR(1) noise
sim <- add_noise(trends, h, lambda_target, ar.coefs, seed)
# Returns: list(y1, y2, x1, x2, noise)
#   noise$series1$sigma  — innovation SD (true parameter)
#   noise$series1$ar     — AR coefficients used
```

### True ACVF and covariance sums

```r
# Theoretical ACVF of AR(1) noise process
acov_raw <- arma_acov(ar = phi, ma = numeric(0), sigma2 = sigma2_innov,
                      lag_max = 100L + h - 1L)
# Returns: numeric vector gamma(0), gamma(1), ..., gamma(lag_max)

# Apply MA(h) filter to get smoothed noise ACVF
# (copy ma_filter_acov() from dev/sims/calibration/run.R — not exported)
acov_eta <- ma_filter_acov(acov_raw, h, lag_max = 100L)

# Compute covariance sums L1, L2, Q1, Q2, Q12
sums <- acov_sums(acov_eta, acov_eta)
# Returns: list(L1, L2, Q1, Q2, Q12)
```

### Theoretical approximations (Proposition 1)

```r
# Signal variance over each rolling window
tau_sq <- compute_tau_sq(trend, s)
# trend: shared trend (true or estimated), e.g. (x1 + x2) / 2
# Returns: numeric vector of length n, NA for first s-1 positions

# Local population correlation
rho <- compute_rho(tau_sq, sigma1_sq = acov_eta[1], sigma2_sq = acov_eta[1])
# Returns: numeric vector of same length as tau_sq

# Asymptotic variance V (SE² of R_t is V/s)
V <- compute_V(tau_sq, sigma1_sq = acov_eta[1], sigma2_sq = acov_eta[1],
               L1 = sums$L1, L2 = sums$L2,
               Q1 = sums$Q1, Q2 = sums$Q2, Q12 = sums$Q12)
```

### Computing Rₜ from data

The rolling correlation R_t must be computed manually; there is no exported
wrapper. Use the same calculation as in `dev/sims/calibration/run.R`:

```r
h_win <- h  # smoothing window
s_win <- s  # correlation window

ma_y1 <- as.numeric(stats::filter(sim$y1, rep(1/h_win, h_win), sides = 1))
ma_y2 <- as.numeric(stats::filter(sim$y2, rep(1/h_win, h_win), sides = 1))

R <- rep(NA_real_, n)
for (t in s_win:n) {
  w <- (t - s_win + 1L):t
  if (any(is.na(ma_y1[w])) || any(is.na(ma_y2[w]))) next
  if (sd(ma_y1[w]) == 0 || sd(ma_y2[w]) == 0) next
  r_val <- cor(ma_y1[w], ma_y2[w])
  if (is.finite(r_val) && abs(r_val) < 1) R[t] <- r_val
}
```

---

## Recommended protocol

### Single condition

```r
# 1. Fix condition parameters
n <- 1000; h <- 10; s <- 150; snr <- 1.5; phi <- 0.5; d <- 0

# 2. Run R replicates, collect R_t vectors
R_reps <- matrix(NA_real_, nrow = n, ncol = R_reps_count)
for (r in seq_len(R_reps_count)) {
  trends <- make_trends_dist(n = n, d = d, seed = base_seed + r)
  sim    <- add_noise(trends, h = h, lambda_target = snr,
                      ar.coefs = phi, seed = base_seed + r + 1L)
  # ... compute R_t as above, store in R_reps[, r]
}

# 3. Compute empirical mean and variance at each t
emp_mean <- rowMeans(R_reps, na.rm = TRUE)
emp_var  <- apply(R_reps, 1, var, na.rm = TRUE)

# 4. Compute theoretical rho and V/s using true parameters
sigma2_innov <- mean(c(sim$noise$series1$sigma, sim$noise$series2$sigma)^2)
acov_raw     <- arma_acov(ar = phi, ma = numeric(0), sigma2 = sigma2_innov,
                           lag_max = 100L + h - 1L)
acov_eta     <- ma_filter_acov(acov_raw, h, 100L)
sums         <- acov_sums(acov_eta, acov_eta)
true_trend   <- (trends$x1 + trends$x2) / 2
tau_sq       <- compute_tau_sq(
  as.numeric(stats::filter(true_trend, rep(1/h, h), sides = 1)), s)
rho_theory   <- compute_rho(tau_sq, acov_eta[1L], acov_eta[1L])
V_theory     <- compute_V(tau_sq, acov_eta[1L], acov_eta[1L],
                           sums$L1, sums$L2, sums$Q1, sums$Q2, sums$Q12)

# 5. Compute errors at valid time points
valid <- is.finite(rho_theory) & is.finite(V_theory) &
         is.finite(emp_mean)   & is.finite(emp_var)
err_rho <- abs(emp_mean[valid] - rho_theory[valid])
err_V   <- abs(emp_var[valid]  - V_theory[valid] / s)
```

### Across conditions

Loop over the factorial and store, for each condition, a summary row:
median absolute error in ρ, median relative error in ρ, median absolute
error in V/s, median relative error in V/s. Collect into a data frame for
table/figure production.

---

## Required outputs

### Tables
- Mean-approximation error (median |E[Rₜ] − ρₜ|) by factor level
- Variance-approximation error (median |Var(Rₜ) − Vₜ/s|) by factor level

### Figures
- Scatter: empirical E[Rₜ] vs theoretical ρₜ (one point per valid *t* × condition)
- Scatter: empirical Var(Rₜ) vs theoretical Vₜ/s
- Error vs *h*, error vs *s*, error vs SNR, error vs trend mismatch *δ*

---

## Student instructions

### Getting started

Source the package and verify the key functions are available:

```r
devtools::load_all()
# Check:
?arma_acov
?compute_rho
?compute_V
?compute_tau_sq
?acov_sums
```

`ma_filter_acov()` is not exported — copy it from
`dev/sims/calibration/run.R` into your script.

### Suggested development order

1. Implement the single-condition protocol above for *one* condition (d=0,
   h=10, s=150, λ=1.5, T=1000) with *R* = 50 replicates. Produce a scatter
   plot of E[Rₜ] vs ρₜ. It should fall near the diagonal.

2. Once working, wrap the single-condition code into a function
   `run_condition(n, h, s, snr, phi, d, R_reps, seed)` that returns a
   summary row (median errors).

3. Build the factorial loop and collect results.

4. Produce the required tables and figures.

### What to look for

- **ρ approximation**: errors should be small (< 0.05 in absolute terms)
  for the operating regime of the paper (h ≥ 5, s ≥ 50, λ ≥ 1.5, d ≤ 0.5).
  If errors grow materially with trend mismatch δ or with small h, note this.

- **V approximation**: this is the harder quantity. Expect larger relative
  errors, especially for small *s* or large *h*. If the V approximation is
  only qualitative (order-of-magnitude correct but not precise), that is a
  finding worth documenting.

- **Primary determinants**: identify which factor (h, s, SNR, δ) most
  strongly drives approximation error. This is the key interpretive result.

### Deliverable

Return a concise written summary (1–2 pages) covering:
1. In what regimes is the ρ approximation reliable?
2. Is ρ consistently better approximated than V?
3. When does V break down materially?
4. What are the primary determinants of error?
5. A recommendation: should the paper address approximation quality in the
   main text, appendix, or not at all?
