# Calibration Study Design

## Objective

Validate the CLT-based approximation (Proposition 1) by checking that the
pointwise test is well-calibrated under H₀ — i.e. the BY-adjusted rejection
rate is at or below the nominal level α = 0.05 across a range of conditions.

## Data-generating process

Both series share an identical trend (*d* = 0) with independent AR(1) noise:

```
X_{kt} = ν_t + Z_{kt},   Z_{kt} ~ AR(1, φ, σ²),   k = 1, 2
```

Trends are generated via `make_trends_dist()` with *d* = 0 (shared Fourier
basis, no separation). Seeds are drawn fresh for each replicate.

## Parameter grid

| Factor | Levels |
|--------|--------|
| Series length *T* | 250, 500, 1000 |
| AR(1) coefficient *φ* | 0.0, 0.3, 0.5, 0.8 |
| SNR *λ* | 0.5, 1.0, 2.0 |
| Smoothing window *h* | auto (`max(5, floor(T/200))`) |
| Correlation window *s* | auto (`min(60h, floor(T/4))`) |

*S* = 500 replicates per cell. Total cells: 3 × 4 × 3 = 36.
Total replicates: 18,000.

## Estimands

For each replicate, record:

- **rejection_rate**: proportion of valid time points flagged at α = 0.05
  (BY-FDR adjusted) — primary calibration diagnostic
- **arma_order**: fitted ARMA order (p, q) for each series — to track model
  selection behaviour
- **sigma_hat**: estimated smoothed noise variance γ̂η(0) — to assess noise
  estimation accuracy
- **n_valid**: number of valid time points used in inference

## Expected outputs

`results/summary.rds` — data frame with one row per replicate, columns for
all grid parameters plus the estimands above.

`figures/calibration_rejection_rate.png` — faceted plot of mean rejection rate
(± 2 SE) by *T*, *φ*, and *λ*, with α = 0.05 reference line.

`figures/calibration_arma_order.png` — heatmap of selected ARMA order
frequency by *φ*.

## Acceptance criterion

Mean rejection rate within [0, 0.08] across all cells, with no cell exceeding
0.15 at *S* = 500 replicates.
