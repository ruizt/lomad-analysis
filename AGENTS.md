# lomad Package — Source Conventions

This document describes the conventions used to organize the `R/` source
directory, write tests, and handle legacy code. It is intended for both human
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
| `utils-lomad-fit.R` | Fit implementations | `.lomad_fit_clt`, `.lomad_fit_state` [LEGACY], `.lomad_fit_blocks` [LEGACY] |
| `utils-lomad-test.R` | Test implementations | `.lomad_test_clt`, `.lomad_test_boot` [LEGACY], `.lomad_test_mc` [LEGACY], `.lomad_test_analytic` [LEGACY], `.boot_progress` [LEGACY] |

### Convention summary

- **All internals are dot-prefixed**: `.foo()`, not `foo()`.
- **No duplicated code**: shared helpers live in exactly one utils file.
- **Legacy code is marked `[LEGACY]`** in file headers and function comments.
  See `dev/legacy-sunset.md` for removal instructions.

## Function dispatch

### `lomad_fit()`

Takes `method = c("clt", "state")`:
- `"clt"` (default): dispatches to `.lomad_fit_clt()` — the paper method.
- `"state"`: dispatches to `.lomad_fit_state()` — legacy state-process pipeline.
- If `blocks` argument is supplied, forces `"state"` and calls `.lomad_fit_blocks()`.
- Optional `noise_override` (CLT only): a list with `ar`, `ma` (optional), `sigma2`
  to bypass noise estimation. Accepts a single spec (shared for both series) or a
  list of two specs (one per series). Used for oracle experiments.

All fit objects include `$method` (`"clt"` or `"state"`) for automatic dispatch
by `lomad_test()`.

### `lomad_test()`

Takes `fit` and optional `method`:
- For CLT fits: auto-dispatches to `.lomad_test_clt()`. No `method` needed.
- For state fits: `method = c("boot", "mc", "analytic")` selects the test.

### `lomad_test_identity()`

Standalone pointwise test of exact trend identity. Optional `noise_override`
(list with `ar`, `ma` (optional), `sigma2` for one noise series) bypasses
ARMA estimation for oracle experiments. The difference variance is computed
internally as twice the single-series variance.

### `lomad()`

Convenience wrapper: `lomad_fit()` → `lomad_test()`. Returns `list(fit, test)`.

## Testing conventions

Tests live in `tests/testthat/` with one file per theme:

| File | Covers |
|------|--------|
| `test-clt.R` | `arma_acov`, `acov_sums`, `compute_rho`, `compute_V`, `compute_tau_sq`, `.ma_filter_acov` |
| `test-sim.R` | `sim_trends`, `sim_noise`, `sim_noise_pair`, key internals |
| `test-estimate.R` | `estimate_trends`, `estimate_ar1_noise`, `estimate_arma_noise`, `.variogram_*`, `.select_arma` |
| `test-lomad.R` | `lomad`, `lomad_fit(method="clt")`, `lomad_test` — structural tests |
| `test-lomad-identity.R` | `lomad_test_identity` |
| `test-lomad-legacy.R` | Legacy state pipeline |

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
├── lomad_fit.R               # lomad_fit(), blocks_from_df()
├── lomad_test.R              # lomad_test()
├── lomad_test_identity.R     # lomad_test_identity()
├── lomad_plot.R              # lomad_plot(), .shade_intervals
├── utils-lomad-fit.R         # .lomad_fit_clt, .lomad_fit_state [LEGACY], .lomad_fit_blocks [LEGACY]
├── utils-lomad-test.R        # .lomad_test_clt, .lomad_test_boot [LEGACY], .lomad_test_mc [LEGACY],
│                             #   .lomad_test_analytic [LEGACY], .boot_progress [LEGACY]
│
├── arma_acov.R               # arma_acov(), acov_sums()
└── utils-clt.R               # .ma_filter_acov, .ma_filtered_var, compute_tau_sq, compute_rho, compute_V
```
