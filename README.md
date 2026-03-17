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
| `lomad_fit()` | Fit the null model: smooth trends, fit ARMA residuals, compute rolling correlations, classify decoupling periods |
| `lomad_test()` | Test for decoupling: parametric bootstrap p-values for entry rate, mean run length, fraction of time decoupled, and number of episodes |
| `lomad()` | Convenience wrapper: runs both steps in one call |

Additional inference methods are in development and will eventually be
selectable via a `method` argument to `lomad()`:

| Method | Description |
|---|---|
| `"bootstrap"` | Full parametric bootstrap (current default) |
| `"mc_bootstrap"` | Bootstrap on the Markov state process only (faster) |
| `"analytic"` | Closed-form CLT via Markov chain theory |

## Installation

```r
# Install from source (development)
devtools::install()
```

## Quickstart

```r
library(lomad)

# Simulate paired trend series and add calibrated ARMA noise
trends <- make_trends_dist(n = 500, d = 5, seed = 1)
sim    <- add_noise(trends, h = 30, lambda_target = 4, scale = 5,
                    order = c(2, 1), seed = 2)

# Fit null model and run parametric bootstrap
out <- lomad(sim$y1, sim$y2, q = 30, h = 50, B = 500, seed = 3)

out$observed   # observed decoupling statistics
out$p_values   # bootstrap p-values

# Visualise
plot_lomad_fit(sim$y1, sim$y2, out)
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
├── R/                  # Package functions
├── dev/
│   ├── notebooks/      # Quarto simulation studies
│   └── scripts/        # Exploratory R scripts
├── tests/testthat/     # Unit tests
├── man/                # Auto-generated documentation
├── DESCRIPTION
└── NAMESPACE
```
