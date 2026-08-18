# lomad-analysis

Reproduces results in:

> Ruiz, T. D., Seifert, A. J., Hamilton, E., Mispagel, C. M., Hunt, O. P.,
> Garcia, J., and Bockmon, E. E. (2026). Inference for local trend similarity in
> nonstationary time series via rolling correlation, with application to
> assessing stability in an estuarine system. *Manuscript in preparation.*

The method itself is implemented in an R package in the companion repository
[`lomad-package`](https://github.com/ruizt/lomad-package). For the method, the
package vignette (`vignette("lomad")`) is the better starting point.

## Layout

```
lomad-analysis/
├── figure-theme.R              # shared plot sizing and theme helpers
├── simulations/
│   ├── simulation-results.R    # all paper figures + tables, from compiled results
│   ├── Dockerfile              # shared image for the cluster jobs
│   ├── build-image.sh          # builds and pushes that image
│   ├── power/                  # local power vs realized separation
│   ├── realignment/            # local realignment as a preprocessing strategy
│   └── validation/             # finite-sample accuracy of the CLT
├── mb-analysis/                # Morro Bay field data analysis
│   ├── process_blocks.R        # sensor record -> analysis blocks
│   ├── analysis.R              # presmoothing, fitting, pooled inference, figures
│   ├── sensitivity.R           # detection stability across h and s
│   ├── export_example.R        # regenerates the block shipped with the package
│   └── utils.R
└── numerical-checks/           # numerical verification of the theory
```

`_`-prefixed entries are gitignored throughout.

## Setup

```r
pak::pak("ruizt/lomad-package")   # or remotes::install_github(...)
library(lomad)
```

Scripts assume they are run **from the repository root**, e.g.
`source("simulations/power/simulation-template.R")`.

## Simulations

Two studies, each with its own directory:

- **`simulations/validation/` — does the asymptotic theory hold at finite $n$?**
- **`simulations/power/` — can the test find a separation, and where?**

In addition, a brief supplemental experiment:

- **`simulations/realignment/` — is local realignment a workable alternative?**

The first two run on the cluster and each hold:

- `design.md` describes the study design and implementation details
- `simulation-template.R` gives a local illustration at small scale (writes nothing)
- `tide/` containing files to execute study at full scale on a Kubernetes cluster
- `collect-results.R` assembles fetched results from cluster jobs
- `results/` containing raw and assembled results 

The cluster jobs use a shared image stored at `ghcr.io/ruizt/lomad-simulations:latest`.
The results of cluster jobs are archived on Zenodo, but everything needed to 
reproduce the paper's figures is contained in the assembled results stored in the 
repository. The following will generate results as they appear in the paper:

```bash
Rscript simulations/simulation-results.R    # rebuilds all paper figures and tables
```

## Numerical checks

Verification of certain derivations are stored in `numerical-checks` and check algebra
and analysis that appears in the theory. Nothing in the paper's figures depends on
them.

| Script | Paper result | What it checks |
|---|---|---|
| `validate-gradient.R` | Theorem 2 ($V = \nabla g^T \Sigma \nabla g$) | Analytic $\nabla g(\theta)$ against finite differences, and that the full $5\times5$ quadratic form equals the reduced $3\times3$ one in centered coordinates |
| `validate-prop1.R` | Proposition 1 | $\Sigma_{(3,4,5)}$ against Monte Carlo with two distinct signals, $\rho = r\tau_1\tau_2/\sqrt{BD}$, $V$ under $H_0$, and the attenuation identity $\rho = (1-\delta^2)^{1/2}\rho^{(0)}$ |
| `validate-rho-formula.R` | Appendix ($\rho_t$ under the simulation model) | The general $\rho_t$ formula against empirical $R_t$ under the `rs` and `fr` coupling structures |

## Morro Bay analysis

Analysis of CeNCOOS data from Morro Bay during the five-year period 2020-2025.

| Script | What it does |
|---|---|
| `process_blocks.R` | Builds the analysis blocks from `wp_data.parquet`: hourly binning, a 24-hour gap threshold, a 5-day minimum, and global standardization of pH and dissolved oxygen |
| `analysis.R` | Presmooths, fits and tests each block, pools p-values across blocks for one Benjamini--Yekutieli correction, and writes every figure the paper uses |
| `sensitivity.R` | Refits the analysis across a grid of smoothing bandwidths and window widths |
| `utils.R` | Spectral notch filter that removes tidal periodicity, plus the downsampling helper |
| `export_example.R` | Regenerates the one Bay Mouth block (late summer 2022) shipped with the package as `morro_bay` |

The scripts read and write `_mb-data/` at the repository root, which is
gitignored and not distributed. They assume `wp_data.parquet` — the
quality-controlled sensor record — has been placed there; nothing in this
repository reaches outside it to fetch the data. That record is produced by
[`mb-qartod`](https://github.com/ruizt/mb-qartod), which handles the
quality control.
