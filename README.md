# lomad-analysis

Analysis, simulation studies, and HPC scaffolding for the **lomad** method
(**Lo**cal **M**oving **A**verage **D**ecoupling) — detecting periods of
decoupling between two correlated time series.

The R package itself lives in the companion repository
[`lomad-package`](https://github.com/ruizt/lomad-package). This repository
contains everything built *on top of* that package: the paper's simulation
studies, the Morro Bay field-data analysis, numerical validation of the theory,
and exploratory notebooks and scripts.

## Setup

Install the package from GitHub, then load it at the top of any script or
notebook:

```r
remotes::install_github("ruizt/lomad-package")
library(lomad)
```

Re-run `remotes::install_github("ruizt/lomad-package")` to pick up upstream
changes to the package. (To iterate on the package and this analysis together,
clone `lomad-package` as a sibling directory and `devtools::load_all()` it
instead.)

Scripts and docs assume they are run **from the repository root** (e.g.
`source("simulations/power/template.R")`, `Rscript simulations/power/tide/collect.R`).

## Structure

```
lomad-analysis/
├── simulations/            # Simulation studies + Tide HPC scaffolding
│   ├── power/              # Power curves and localization across trend structures
│   ├── validation/         # Finite-sample accuracy of the CLT approximation
│   └── trend_examples.R    # Trend-construction illustrations
├── mb-analysis/            # Morro Bay field data analysis (CLT + ARMA noise)
└── numerical-validation/   # End-to-end numerical validation scripts
```

See `simulations/README.md` for the simulation infrastructure and Tide workflow.

## Data

`mb-analysis/import_mb_data.R` copies cleaned Morro Bay data from the sibling
[`mb-qartod`](https://github.com/ruizt/mb-qartod) repository into `_mb-data/`
(gitignored). Keep `mb-qartod` checked out as a sibling of this repo.

## Keeping files out of version control

Two options, both covered by `.gitignore`:

- **Prefix with `_`** — for individual files alongside tracked ones (e.g. `_my_scratch.R`)
- **`scratch/` folder** — drop anything in `scratch/` and it won't be tracked
