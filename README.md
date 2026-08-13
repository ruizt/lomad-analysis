# lomad-analysis

> Ruiz, T. D., Seifert, A. J., Hamilton, E., Mispagel, C. M., Hunt, O. P.,
> Garcia, J., and Bockmon, E. E. (2026). Inference for local trend similarity in
> nonstationary time series via rolling correlation, with application to
> assessing stability in an estuarine system. *Manuscript in preparation.*

Everything behind the **lomad** paper — the simulation studies, the Morro Bay
field analysis, and numerical checks of the theory not reported in the paper. 
The method itself is implemented in an R package in the companion repository
[`lomad-package`](https://github.com/ruizt/lomad-package).

For the method itself, the package vignette (`vignette("lomad")`) is the
better starting point.

## Layout

```
lomad-analysis/
├── figure-theme.R              # shared plot sizing and theme helpers
├── simulations/
│   ├── simulation-results.R    # all paper figures + tables, from compiled results
│   ├── Dockerfile              # shared image for the cluster jobs
│   ├── build-image.sh          # builds and pushes that image
│   ├── power/                  # local power vs realized separation
│   └── validation/             # finite-sample accuracy of the CLT
├── mb-analysis/                # Morro Bay field data analysis
│   ├── process_blocks.R        # sensor record -> analysis blocks
│   ├── analysis.R              # presmoothing, fitting, pooled inference, figures
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

Each study directory holds:

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

`mb-analysis/process_blocks.R` copies `wp_data.parquet` — the
quality-controlled sensor record — from the sibling
[`mb-qartod`](https://github.com/ruizt/mb-qartod) repository into `_mb-data/`
(gitignored) and builds the analysis blocks from it, so keep `mb-qartod`
checked out alongside this repo. Block construction lives here rather than
upstream because every choice in it is an analysis decision: hourly binning, a
24-hour gap threshold, a 5-day minimum, and global standardization of pH and
dissolved oxygen. Tide and pressure are carried through but excluded from the
block definition, so their missingness cannot move block boundaries.

`mb-analysis/analysis.R` then presmooths (spectral notch at the tidal bands,
downsample to 6-hourly), fits each block at $h = 4$ and $s = 60$, pools the raw
p-values across all blocks for a single Benjamini--Yekutieli correction, and
writes figures to `mb-analysis/_img/`. Presmoothing stays in the analysis script
rather than the processing script so that both the raw and presmoothed series
are available for figures.

Blocks shorter than $2.5s + h$ are dropped before fitting. This is a second,
stricter length restriction than the 5-day minimum applied upstream, and it is
an inference decision rather than a data one: a block barely longer than the
window contributes few windows, all of them heavily overlapping, and the FDR
over such a block is not controlled at the nominal level. Pooling p-values
across blocks means one short block degrades the correction for every other, so
the restriction is applied before pooling rather than after.

One Bay Mouth block (late summer 2022) ships with the package as `morro_bay`;
`mb-analysis/export_example.R` regenerates it.
