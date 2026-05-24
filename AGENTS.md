# lomad Package — Source Conventions

This document describes the conventions used to organize the `R/` source
directory, write tests, and handle related code. It is intended for both human
developers and AI agents working on the package.

## File naming

### Exported functions

Exported functions live in `{theme}_{function}.R`. The theme prefix groups
related functions:

| Prefix | Theme | Files |
|--------|-------|-------|
| `sim_` | Simulation / data generation | `sim_trends.R`, `sim_noise.R` |
| `estimate_` | Parameter estimation | `estimate_trends.R`, `estimate_noise.R` |
| `lomad_` | Method implementations | `lomad.R`, `lomad_fit.R`, `lomad_test.R`, `lomad_test_identity.R`, `lomad_plot.R` |
| *(none)* | CLT building blocks | `arma_acov.R` |

### Internal utilities

Internal helpers live in `utils-{theme}.R`:

| File | Theme | Contents |
|------|-------|----------|
| `utils-sim.R` | Simulation | `.generate_fourier_coef`, `.generate_coef_pair`, `.make_basis_trends`, `.make_w_smooth`, `.make_w_cross`, `.make_w_rate`, `.apply_w`, `.pacf_to_arma_coefs` |
| `utils-estimate.R` | Estimation | `.variogram_ar1`, `.variogram_ar`, `.yule_walker`, `.fit_arma`, `.select_arma`, `.smooth_noise_var` |
| `utils-clt.R` | CLT theory | `.ma_filter_acov`, `.ma_filtered_var`, `compute_tau_sq`, `compute_rho`, `compute_V` |
| `utils-lomad-fit.R` | Fit implementation | `.lomad_fit_clt` |
| `utils-lomad-test.R` | Test implementation | `.lomad_test_clt` |

### Convention summary

- **All internals are dot-prefixed**: `.foo()`, not `foo()`.
- **No duplicated code**: shared helpers live in exactly one utils file.

## Function interfaces

### `lomad_fit()`

```r
lomad_fit(x1, x2, h = NULL, s = NULL, noise_override = NULL, lag_max = 100L)
```

Runs the CLT pipeline: MA trend estimation, noise fitting, and rolling
computation of R_t, rho_t, V_t.

- `noise_method = c("ar1", "arma")`: selects the noise model. `"ar1"`
  (default) uses variogram-based AR(1); `"arma"` uses BIC-selected AR(p)
  via `estimate_arma_noise()`. Use `"arma"` for data with periodic or
  complex noise structure (e.g. tidal signals).
- Optional `noise_override`: a list with `ar`, `ma` (optional), `sigma2`
  to bypass noise estimation entirely. Accepts a single spec (shared for
  both series) or a list of two specs (one per series). Used for oracle
  experiments and for externally estimated noise parameters.

### `lomad_test()`

```r
lomad_test(fit, alpha = 0.05)
```

BY-corrected pointwise Z-test on the fit.

### `lomad_test_identity()`

Standalone pointwise test of exact trend identity. Optional `noise_override`
(list with `ar`, `ma` (optional), `sigma2` for one noise series) bypasses
ARMA estimation for oracle experiments. The difference variance is computed
internally as twice the single-series variance.

### `lomad()`

Convenience wrapper: `lomad_fit()` → `lomad_test()`. Returns `list(fit, test)`.

### `lomad_plot()`

```r
lomad_plot(fit, tst, dates = NULL, alpha = 0.25)
```

Two-panel base R plot. Upper panel: smoothed series (ma1, ma2, trend) with
`tst$rejected` shading. Lower panel: rolling correlation R_t vs rho_t.

## Testing conventions

Tests live in `tests/testthat/` with one file per theme:

| File | Covers |
|------|--------|
| `test-clt.R` | `arma_acov`, `acov_sums`, `compute_rho`, `compute_V`, `compute_tau_sq`, `.ma_filter_acov` |
| `test-sim.R` | `sim_trends`, `sim_noise`, `sim_noise_pair`, key internals |
| `test-estimate.R` | `estimate_trends`, `estimate_ar1_noise`, `estimate_arma_noise`, `.variogram_*`, `.select_arma` |
| `test-lomad.R` | `lomad`, `lomad_fit`, `lomad_test` — structural tests |
| `test-lomad-identity.R` | `lomad_test_identity` |

### Testing internals

In `testthat`, test files run inside the package namespace, so dot-prefixed
internals are callable directly in tests without `:::`. Test internals directly
when they have nontrivial logic; test trivial wrappers indirectly via exports.

### CRAN rules

- Tests must complete in < 5 seconds total on CRAN (use small `n`).
- Wrap slow or stochastic tests in `skip_on_cran()`.
- Use `set.seed()` for reproducibility.
- Use generous tolerances for statistical checks.

## Data generation for simulations

The pipeline for generating test data is:

```r
sim_trends()  -->  sim_noise_pair()  -->  observed series (y1, y2)
```

- `sim_trends()` generates trend pairs with controlled L2 distance `d` and
  coupling weight `w` (methods: `"dist"`, `"smooth"`, `"cross"`, `"rate"`,
  or custom `w` vector/function).
- `sim_noise_pair()` adds calibrated ARMA noise to both series at a target SNR.

## R/ directory layout

```
R/
├── sim_trends.R              # sim_trends()
├── sim_noise.R               # sim_noise(), sim_noise_pair()
├── utils-sim.R               # .generate_*, .make_w_*, .apply_w, .pacf_to_arma_coefs
│
├── estimate_trends.R         # estimate_trends()
├── estimate_noise.R          # estimate_ar1_noise(), estimate_arma_noise()
├── utils-estimate.R          # .select_arma, .variogram_*, .yule_walker, .fit_arma, .smooth_noise_var
│
├── lomad.R                   # lomad()
├── lomad_fit.R               # lomad_fit()
├── lomad_test.R              # lomad_test()
├── lomad_test_identity.R     # lomad_test_identity()
├── lomad_plot.R              # lomad_plot(), .shade_intervals
├── utils-lomad-fit.R         # .lomad_fit_clt
├── utils-lomad-test.R        # .lomad_test_clt
│
├── arma_acov.R               # arma_acov(), acov_sums()
└── utils-clt.R               # .ma_filter_acov, .ma_filtered_var, compute_tau_sq, compute_rho, compute_V
```
