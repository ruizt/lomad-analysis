# lomad

**Lo**cal **M**oving **A**verage **D**ecoupling — an R package for detecting
and characterizing periods of decoupling between two correlated time series.

## Overview

`lomad` provides a statistical framework for detecting time windows in which
two previously correlated series exhibit locally low correlation. Given paired
series following a signal-plus-noise model, the method:

1. Smooths each series with a moving average filter to isolate trends.
2. Computes rolling window correlations on the smoothed series.
3. Tests pointwise whether observed correlations fall below their expected
   values under a shared-trend null, using a CLT-based test statistic with
   Benjamini–Yekutieli FDR correction.

Noise parameters (AR/ARMA) are estimated via a variogram-based approach that
is robust to trend contamination. The asymptotic variance of the local
correlation accounts for autocorrelation in the smoothed noise.

### Core workflow

| Function | Role |
|---|---|
| `lomad_fit()` | Fit the null model (estimate noise, compute expected correlations and asymptotic variance) |
| `lomad_test()` | Pointwise test for local decoupling with FDR correction |
| `lomad()` | Convenience wrapper: `lomad_fit()` → `lomad_test()` |

### Supporting functions

| Function | Role |
|---|---|
| `lomad_test_identity()` | Global test of exact trend identity (H₀: d = 0) |
| `lomad_plot()` | Visualise fit and test results |
| `estimate_trends()` | Extract trends via moving average |
| `estimate_ar1_noise()` | Estimate AR(1) noise parameters via variogram |
| `estimate_arma_noise()` | Estimate ARMA noise parameters via long-AR approximation |
| `sim_trends()` | Generate synthetic trend pairs at controlled L² separation |
| `sim_noise_pair()` | Add calibrated ARMA noise to trend pairs |

## Installation

### Development version (from source)

Clone the repository and install with `devtools`:

```bash
git clone <repo-url>
cd lomad
```

```r
devtools::install_deps()
devtools::install()
```

Or build and install from the command line:

```bash
R CMD build .
R CMD INSTALL lomad_*.tar.gz
```

To load the package in-place during development:

```r
devtools::load_all()
```

## Quickstart

```r
library(lomad)

# Simulate paired trends with controlled L² separation
trends <- sim_trends(n = 500, d = 2, method = "smooth", bw = 50, seed = 1)

# Add calibrated AR(1) noise at target SNR
sim <- sim_noise_pair(trends, h = 10, lambda_target = 1.5,
                      ar.coefs = 0.5, seed = 2)

# Fit null model and test for decoupling
fit <- lomad_fit(sim$y1, sim$y2, h = 10, s = 50)
tst <- lomad_test(fit, alpha = 0.05)

# Which time points show significant decoupling?
which(tst$rejected)

# Visualise
lomad_plot(fit, tst)
```

### Convenience wrapper

```r
out <- lomad(sim$y1, sim$y2, h = 10, s = 50, alpha = 0.05)
out$fit   # lomad_fit output
out$test  # lomad_test output
```

## For contributors

Clone the repo and open `lomad.Rproj` in RStudio. Load all package functions
with `devtools::load_all()`.

Development files live in `dev/`:

| Directory | Contents |
|-----------|----------|
| `dev/sims/` | Simulation studies (calibration, power, Tide HPC example) |
| `dev/mb-analysis/` | Morro Bay field data analysis |
| `dev/scripts/` | Exploratory R scripts |

See `AGENTS.md` for source conventions and `dev/sims/README.md` for the
simulation infrastructure.

## Project structure

```
lomad/
├── R/                     # Package source (themed: sim_*, estimate_*, lomad_*, utils-*)
├── tests/testthat/        # Unit tests (150 tests, themed by module)
├── man/                   # Auto-generated documentation
├── dev/
│   ├── sims/              # Simulation studies + Tide HPC scaffolding
│   │   ├── calibration/   # CLT test calibration study
│   │   ├── power/         # Power curves across trend structures
│   │   └── tide-example/  # Minimal MVN example to learn the Tide workflow
│   ├── mb-analysis/       # Morro Bay field data analysis
│   └── scripts/           # Exploratory scripts
├── DESCRIPTION
├── NAMESPACE
└── AGENTS.md              # Source conventions for developers and AI agents
```
