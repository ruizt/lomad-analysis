# Simulation Studies

This directory contains all simulation code for the lomad paper. Each study
has its own subdirectory with a self-contained design document, R code, and
HPC submission materials.

## Studies

| Directory | Purpose |
|-----------|---------|
| `calibration/` | Validate CLT approximation: rejection rates under H₀ across noise/window conditions |
| `power/` | Power curves as a function of separation *d* for each trend structure |
| `examples/` | Illustrative figures for the paper (not batch studies) |

## Conventions

- **`design.md`** — prose description of the simulation design, parameter grid,
  estimands, and intended outputs. This is the source of truth; code should
  match it, not the other way around.
- **`run.R`** — defines a single-replicate function `run_rep(params)` that
  accepts a named list of parameters and returns a named list of outputs.
  A local driver at the bottom runs a small grid for development.
- **`tide/`** — HPC submission materials for the Tide cluster (SLURM). Contains
  an array job script and a submission shell script. These call `run.R` via
  `Rscript`; no other files should be needed on the cluster.
- **`results/`** — gitignored. Populated by HPC runs. Intermediate `.rds` files
  (one per replicate) are combined by a `collect.R` script within each study.
- **`figures/`** — tracked. Final figures produced by `figures.R` within each
  study, for inclusion in the paper.

## Running on Tide

From the study directory (e.g. `calibration/`):

```bash
cd tide/
sbatch submit.sh
```

Results land in `results/` as individual `.rds` files named `rep_<array_id>.rds`.
Collect and summarise with:

```r
source("collect.R")   # produces results/summary.rds
source("figures.R")   # produces figures/
```
