# Legacy Code Sunset — COMPLETED

The legacy state-process pipeline was removed on 2026-05-24.

## What was removed

- `.lomad_fit_state()` — state-process fitting with Fisher-z thresholding
- `.lomad_fit_blocks()` — multi-block variant of the above
- `.lomad_test_boot()` — parametric bootstrap on the full correlation process
- `.lomad_test_mc()` — Markov-chain bootstrap on the binary state process
- `.lomad_test_analytic()` — CLT p-value for `frac_state`
- `.boot_progress()` — progress bar for bootstrap loops
- `.lomad_plot_state()` — base R plot for state fits
- `test-lomad-legacy.R` — test file for legacy pipeline
- `dev/scripts/test-script.R` — dev script using legacy API
- `dev/scripts/compare_methods.R` — dev script comparing legacy and CLT methods

## What was simplified

- `lomad_fit()` — removed `method`, `q`, `alpha`, `rho0`, `max_pq` parameters
- `lomad_test()` — removed `method`, `B`, `seed`, `ncores`, `verbose`, `T_sim`, `T_eff` parameters
- `lomad()` — removed `method`, `test_method` parameters
- `lomad_plot()` — removed state dispatch, `x1`/`x2` arguments

## What was retained

- `.select_arma()` and `.smooth_noise_var()` in `utils-estimate.R` — still
  used by `lomad_test_identity()`.

## What was added

- `noise_method = c("ar1", "arma")` parameter on `lomad_fit()` (default
  `"ar1"` for backward compatibility). Use `"arma"` for data with periodic
  or complex noise structure.

## Reduction

~900 lines of implementation code removed. File count unchanged (utils files
retained with only the CLT implementations).
