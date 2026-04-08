# lomad

**Lo**cal **M**oving **A**verage **D**ecoupling — an R package for detecting
and characterizing periods of decoupling between two correlated time series.

## Overview

`lomad` provides a statistical framework for detecting time windows in which
two previously correlated series exhibit locally low correlation. The approach
uses rolling window correlations on moving-average smoothed series, applies
Fisher z-tests with Benjamini-Yekutieli multiple testing correction, and
models the resulting binary state process as a 2-state Markov chain.

Inference is built around a two-step workflow:

| Function | Role |
|---|---|
| `lomad_fit()` | Fit the null model to a single contiguous series pair |
| `lomad_fit_blocks()` | Fit the null model jointly across multiple independent data blocks |
| `lomad_test_boot()` | Test via full parametric bootstrap (p-values for all statistics) |
| `lomad_test_mc()` | Test via Markov-chain bootstrap on state process only (faster) |
| `lomad_test_analytic()` | Test via closed-form CLT (`frac_state` only) |
| `lomad()` | Convenience wrapper: runs both steps via `method = "boot"`, `"mc"`, or `"analytic"` |

## Installation

### Development version (from source)

Clone the repository and install with `devtools`:

```bash
git clone <repo-url>
cd lomad
```

```r
# Install dependencies, then install the package
devtools::install_deps()
devtools::install()
```

Or build a source tarball and install it manually:

```bash
R CMD build .
R CMD INSTALL lomad_*.tar.gz
```

To load the package in-place without installing (useful during development):

```r
devtools::load_all()
```

## Quickstart

### Single series

```r
library(lomad)

# Simulate paired trend series and add calibrated ARMA noise
trends <- make_trends_dist(n = 500, d = 5, seed = 1)
sim    <- add_noise(trends, h = 30, lambda_target = 4, scale = 5,
                    order = c(2, 1), seed = 2)

# Fit null model and test (method = "boot", "mc", or "analytic")
out <- lomad(sim$y1, sim$y2, q = 30, h = 50, B = 500, seed = 3, method = "boot")

out$observed   # observed decoupling statistics
out$p_values   # p-values

# Visualise
plot_lomad_fit(sim$y1, sim$y2, out)
```

### Multiple blocks

When data arrive as independent contiguous chunks (e.g. instrument deployments,
field seasons), use `lomad_fit_blocks()` to pool all blocks into a single null
model estimate. Block boundaries are fully respected: no smoothing, filter
state, or transition counting crosses a boundary.

```r
# blocks is a named list of list(x1, x2) pairs — one per block.
# Use blocks_from_df() to build it from a tidy data frame:
blocks <- blocks_from_df(df, x1_col = "o2", x2_col = "ph", block_col = "block_id")

fit <- lomad_fit_blocks(blocks, q = 56, h = 112, max_pq = 2)

fit$observed            # decoupling statistics aggregated across blocks
fit$expected_asymptotic

# All three test functions accept the output of lomad_fit_blocks() directly.
# lomad_test_boot() simulates one replicate per block at each block's observed
# length, rather than one long aggregate series.
test <- lomad_test_boot(fit, B = 1000, seed = 42,
                        ncores = parallel::detectCores() - 1)
test$p_values
```

## For contributors

Clone the repo and open `lomad.Rproj` in RStudio. Development notebooks and
scripts live in `dev/`. Load all package functions with:

```r
devtools::load_all()
```

See `dev/README.md` for the full contributor workflow.

## Project structure

```
lomad/
├── R/                     # Package functions
├── dev/
│   ├── mb-analysis/       # MB field data analysis scripts
│   ├── notebooks/         # Quarto simulation studies
│   └── scripts/           # Exploratory R scripts
├── tests/testthat/        # Unit tests
├── man/                   # Auto-generated documentation
├── DESCRIPTION
└── NAMESPACE
```
