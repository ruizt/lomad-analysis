# Legacy Code Sunset Instructions

This document describes how to remove the legacy state-process pipeline
from the `lomad` package. The legacy code is marked `[LEGACY]` in source
comments and in `AGENTS.md`.

## What gets removed

The legacy pipeline consists of:
- `.lomad_fit_state()` — state-process fitting with Fisher-z thresholding
- `.lomad_fit_blocks()` — multi-block variant of the above
- `.lomad_test_boot()` — parametric bootstrap on the full correlation process
- `.lomad_test_mc()` — Markov-chain bootstrap on the binary state process
- `.lomad_test_analytic()` — CLT p-value for `frac_state`
- `.boot_progress()` — progress bar for bootstrap loops
- `.select_arma()` — AIC-based ARMA order selection (only if
  `lomad_test_identity` is also dropped; otherwise retain it)

Estimated reduction: ~650 lines of implementation code.

## Steps

### 1. Remove legacy implementations

In `R/utils-lomad-fit.R`:
- Delete `.lomad_fit_state()` and `.lomad_fit_blocks()`.
- Only `.lomad_fit_clt()` should remain. The file can be renamed to
  drop the `utils-` prefix if desired, or the function can be inlined
  into `lomad_fit.R`.

In `R/utils-lomad-test.R`:
- Delete `.lomad_test_boot()`, `.lomad_test_mc()`, `.lomad_test_analytic()`,
  and `.boot_progress()`.
- Only `.lomad_test_clt()` should remain.

### 2. Simplify `lomad_fit.R`

Remove the `method` argument and dispatch logic. Inline `.lomad_fit_clt()`
directly:

```r
lomad_fit <- function(x1, x2, h = NULL, s = NULL, lag_max = 100L) {
  # ... (body of .lomad_fit_clt) ...
}
```

Decide whether to keep `blocks_from_df()`. If block support is not planned
for the CLT pipeline, remove it.

### 3. Simplify `lomad_test.R`

Remove the `method` argument and state-pipeline dispatch. Inline
`.lomad_test_clt()`:

```r
lomad_test <- function(fit, alpha = 0.05) {
  # ... (body of .lomad_test_clt) ...
}
```

### 4. Simplify `lomad.R`

Remove `method` and `test_method` arguments. The wrapper becomes:

```r
lomad <- function(x1, x2, alpha = 0.05, ...) {
  fit <- lomad_fit(x1, x2, ...)
  tst <- lomad_test(fit, alpha = alpha)
  list(fit = fit, test = tst)
}
```

### 5. Clean up `utils-estimate.R`

If `lomad_test_identity()` is also being dropped:
- Remove `.select_arma()` and `.smooth_noise_var()` (no remaining consumers).

If keeping `lomad_test_identity()`:
- `.select_arma()` and `.smooth_noise_var()` must remain.

### 6. Delete empty files

After inlining, `utils-lomad-fit.R` and `utils-lomad-test.R` can be deleted
(or kept as thin files if you prefer not to inline).

### 7. Update documentation

- Remove `method` parameter docs from `lomad_fit`, `lomad_test`, `lomad`.
- Remove cross-references to legacy functions.
- Update `AGENTS.md`: remove `[LEGACY]` markers and legacy file descriptions.

### 8. Update tests

- Delete `tests/testthat/test-lomad-legacy.R`.
- Simplify `test-lomad.R` to remove `method = "state"` tests.

### 9. Verify

```r
devtools::document()
devtools::check()
```

No errors, no warnings about missing exports.

## After sunset

File count drops from 15 to 11-13. The three lomad files (`lomad.R`,
`lomad_fit.R`, `lomad_test.R`) become simple and self-contained with no
dispatch logic.
