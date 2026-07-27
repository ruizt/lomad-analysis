# lomad-analysis

Everything behind the **lomad** paper — the simulation studies, the Morro Bay
field analysis, and the numerical checks of the theory. The method itself is an
R package in the companion repository
[`lomad-package`](https://github.com/ruizt/lomad-package); this repository is
what was built on top of it.

## Setup

```r
pak::pak("ruizt/lomad-package")   # or remotes::install_github(...)
library(lomad)
```

Scripts assume they are run **from the repository root**, e.g.
`source("simulations/power/simulation-template.R")`.

## What to run

Everything you need to reproduce the paper's figures is in the repository — no
downloads required:

```bash
Rscript simulations/simulation-results.R    # rebuilds all paper figures and tables
```

To see a study at work locally, at a scale that runs in seconds rather than on
a cluster, source a study's `simulation-template.R`. For the method itself
rather than the studies, the package vignette (`vignette("lomad")`) is the
better starting point.

The `tide/` directories reproduce the full simulations on a Kubernetes cluster.
That is how the archived results were produced, and it is here for transparency
and for us — it needs cluster access and hours of compute, and you do not need
it to rebuild any figure.

## The simulation studies

Two studies, each with its own `design.md` giving the parameter grid,
estimands, and acceptance criteria. In brief:

**`simulations/power/` — can the test find a separation, and where?**
Trends are simulated at a controlled $L^2$ separation $d$ and the test is run
on the resulting series. The question is how detection depends on $d$, and on
the conditions that make it harder: how the separation is distributed over
time (three *trend structures* — evenly spaced events, irregular episodes, or
episodes where the trends cross), the series length, the signal-to-noise
ratio, and the autocorrelation of the noise. An **oracle arm** re-runs the
hardest cells with the true noise parameters supplied, which separates the
method's limits from the noise estimator's. Two figures come out: power curves
against $d$, and a localization sweep asking whether rejections land where the
trends actually separate.

**`simulations/validation/` — does the asymptotic theory hold at finite $n$?**
Everything here runs under $H_0$ with a single shared trend, so any rejection
is an error. Three things are checked: that the standardized statistic $Z_t$
is approximately standard normal, at three window sizes; that the
Proposition 1 expressions for $\rho_t$ and $V_t$ match their empirical
counterparts; and that the approximation survives replacing those quantities
with plug-in estimates. The composite figure is the paper's validation figure.

Both are run on a Kubernetes cluster because the grids are large — the power
study alone is 360 parameter combinations at 500 replicates. The compiled
results are in this repository, so none of that has to be re-run to use them.

## Layout

```
lomad-analysis/
├── simulations/
│   ├── simulation-results.R    # all paper figures + tables, from compiled results
│   ├── Dockerfile              # shared image for the cluster jobs
│   ├── power/                  # power curves and localization vs separation d
│   └── validation/             # finite-sample accuracy of the CLT
├── mb-analysis/                # Morro Bay field data analysis
└── numerical-checks/           # numerical verification of the theory
```

Each study directory holds `design.md` (the design, parameter grid, and
expected outputs), `simulation-template.R` (a local illustration; writes
nothing), `collect-results.R` (assembles fetched results), `results/`, and
`tide/`. The power study adds `localization-sweep.R`.

## How the pipeline fits together

The stages are split by cost, so redrawing a figure never re-runs a
computation:

| Stage | Script | Cost |
|---|---|---|
| Simulate | `<study>/tide/sim.R` (cluster) | hours |
| Assemble | `<study>/collect-results.R` | seconds |
| Analyse | `power/localization-sweep.R` | minutes (reads ~1.4 GB) |
| Draw | `simulations/simulation-results.R` | seconds |

`tide/` holds only what talks to the cluster. `collect-results.R` and
`localization-sweep.R` read local files, so they sit outside it.

`tide/sim.R` is the source of truth for the simulation logic — it is what runs
and what produced the archived results. `simulation-template.R` mirrors its
`run_rep()` so the local illustration stays faithful. **Change `sim.R` first,
then mirror.**

## Results, and what lives on Zenodo

Compiled results are tracked here (about 3.6 MB), so a fresh clone rebuilds
every figure with no downloads. The bulk — the per-job raw output and its
archive, roughly 2.6 GB — is on Zenodo and is only needed to re-derive the
compiled files from scratch.

The naming carries the rule, and it is the only rule `.gitignore` needs:

- **`_`-prefixed** — bulky, gitignored, archived on Zenodo (`_raw/`, the zips)
- **unprefixed** — small compiled artifacts, tracked here *and* archived
  (`simulations-power-{results,summary,localization}.rds`,
  `simulations-validation-results.rds`)

The same convention covers scratch work: prefix any file with `_`, or drop it
in `scratch/`, and it stays out of version control.

## Reproducing the simulations on Tide

Per study, from the repository root:

```bash
kubectl apply -n cal-poly-ruiz -f simulations/<study>/tide/pvc.yaml  # once
bash    simulations/<study>/tide/submit_sweep.sh                     # submit
kubectl get jobs -n cal-poly-ruiz -l app=lomad-<study>               # monitor
bash    simulations/<study>/tide/fetch.sh                            # retrieve
Rscript simulations/<study>/collect-results.R                        # assemble
kubectl delete jobs -n cal-poly-ruiz -l app=lomad-<study>            # clean up
```

`submit_sweep.sh` mounts `sim.R` into the container through a ConfigMap, so
changing the simulation needs no image rebuild. The Job spec lives once in
`tide/job.yaml` and is filled with `envsubst`; `power/tide/test_one_job.sh`
submits a single small job from the same template as a smoke test.

The image (`ghcr.io/ruizt/lomad-simulations`) is built with, after
`docker login ghcr.io`:

```bash
bash simulations/build-image.sh              # lomad from GitHub, push
bash simulations/build-image.sh --local      # lomad from ../lomad-package
bash simulations/build-image.sh --no-push    # build only
```

By default `lomad` is installed from GitHub at `LOMAD_REF` (default `main`),
which **requires `ruizt/lomad-package` to be public** — `pak` inside the
container has no credentials, so while the repo is private this fails and
`--local` is the working route.

`--local` builds `lomad` from the sibling checkout and installs that source
tree. Worth using even after the repo is public whenever the image should
match your working tree rather than a branch — as when re-running simulations
against a change that is not yet released. Either way the image is labelled
with the package version and commit it came from.

The GitHub Actions build is manual-only for the same reason; re-add its push
trigger once the package repo is public.

The ghcr package must be **public** or the cluster cannot pull it without an
image pull secret.

## Numerical checks

Verification of the *derivations*, as distinct from `simulations/validation/`,
which measures the finite-sample behavior of the method. These check algebra
and analysis — that an analytic gradient matches finite differences, that
Wick's-theorem covariance entries match Monte Carlo, that a stated bound holds
— and they run locally in seconds. Nothing in the paper's figures depends on
them.

The CLT approximation, the Proposition 1 moments, and the end-to-end pipeline
are *not* checked here — `simulations/validation/` measures all three at a
scale these scripts could not match, and duplicating them locally only
produced copies that drifted out of step.

Run from the repository root; scripts that plot save `.png` files alongside
themselves.

| Script | Paper result | What it checks |
|---|---|---|
| `validate-perturbation.R` | Perturbation bound | Under near-common trends, $\rho$ is a perturbation of the common-trend expression with remainder $O(\varepsilon\tau)$ |
| `validate-gradient.R` | Gradient remark (Thm 2) | Analytic $\nabla g(\theta)$ against finite differences; also that the full $5\times5$ $V$ equals the reduced $3\times3$ |
| `validate-sigma-matrix.R` | $\Sigma_{(3,4,5)}$ | Wick's-theorem entries against empirical covariances, for white and ARMA noise |
| `validate-rho-formula.R` | Appendix ($\rho_t$) | The general $\rho_t$ formula against empirical $R_t$ under smooth and rate coupling |

## Morro Bay data

`mb-analysis/process_blocks.R` copies `wp_data.parquet` — the
quality-controlled sensor record — from the sibling
[`mb-qartod`](https://github.com/ruizt/mb-qartod) repository into `_mb-data/`
(gitignored) and builds the analysis blocks from it, so keep `mb-qartod`
checked out alongside this repo. Block construction lives here rather than
upstream because every choice in it is an analysis decision: hourly binning, a
24-hour gap threshold, a 5-day minimum, and global standardization of pH and
dissolved oxygen. Tide and pressure are carried through but excluded from the
block definition, so their missingness cannot move block boundaries.

`mb-analysis/full_analysis.R` then presmooths (spectral notch at the tidal
bands, downsample to 6-hourly), fits each block, pools the raw p-values across
all blocks for a single Benjamini--Yekutieli correction, and writes figures to
`mb-analysis/_img/`. Presmoothing stays in the analysis script rather than the
processing script so that both the raw and presmoothed series are available for
figures.

The two blocks used in the paper also ship with the package as `morro_bay`.
