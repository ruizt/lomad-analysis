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
Rscript simulations/simulation-figures.R    # rebuilds all four paper figures
```

To see a study at work locally, at a scale that runs in seconds rather than on
a cluster, source a study's `simulation-template.R`. For the method itself
rather than the studies, the package vignette (`vignette("lomad")`) is the
better starting point.

The `tide/` directories reproduce the full simulations on a Kubernetes cluster.
That is how the archived results were produced, and it is here for transparency
and for us — it needs cluster access and hours of compute, and you do not need
it to rebuild any figure.

## Layout

```
lomad-analysis/
├── simulations/
│   ├── simulation-figures.R    # all paper figures, from compiled results
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
| Draw | `simulations/simulation-figures.R` | seconds |

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

The image (`ghcr.io/ruizt/lomad-sims`) installs `lomad` from GitHub at build
time. Rebuild it via the *Build and push lomad-sims image* workflow, or by hand
after `docker login ghcr.io`:

```bash
docker buildx build --platform linux/amd64 \
  -f simulations/Dockerfile \
  -t ghcr.io/ruizt/lomad-sims:latest --push .
```

`--platform linux/amd64` is required for the cluster even on Apple Silicon, and
`--build-arg LOMAD_REF=<ref>` pins the package version. The ghcr package must
be **public** or the cluster cannot pull it without an image pull secret.

## Numerical checks

Monte Carlo and numerical verification of the paper's theoretical results, one
script per result. Run from the repository root; scripts that plot save `.png`
files alongside themselves.

| Script | Paper result | What it checks |
|---|---|---|
| `validate-clt.R` | Theorems 1 & 2 | Sampling distribution of $R_t$ against the CLT approximation — QQ, Shapiro–Wilk, and interval coverage at several window sizes |
| `validate-prop1-rho.R` | Proposition 1 ($\rho$) | Theoretical $\rho_t$ against mean empirical $R_t$, and convergence as $s$ grows |
| `validate-prop1-V.R` | Proposition 1 ($V$) | Theoretical $V_t$ against $s \cdot \mathrm{Var}(R_t)$ |
| `validate-perturbation.R` | Perturbation bound | Under near-common trends, $\rho$ is a perturbation of the common-trend expression with remainder $O(\varepsilon\tau)$ |
| `validate-gradient.R` | Gradient remark (Thm 2) | Analytic $\nabla g(\theta)$ against finite differences; also that the full $5\times5$ $V$ equals the reduced $3\times3$ |
| `validate-sigma-matrix.R` | $\Sigma_{(3,4,5)}$ | Wick's-theorem entries against empirical covariances, for white and ARMA noise |
| `validate-rho-formula.R` | Appendix ($\rho_t$) | The general $\rho_t$ formula against empirical $R_t$ under smooth and rate coupling |
| `validate-end-to-end.R` | Full pipeline | `lomad_fit()` + `lomad_test()` on known parameters: $\hat\rho$, $\hat V$, and coverage against oracle values |

## Morro Bay data

`mb-analysis/import_mb_data.R` copies cleaned data from the sibling
[`mb-qartod`](https://github.com/ruizt/mb-qartod) repository into `_mb-data/`
(gitignored), so keep `mb-qartod` checked out alongside this repo. The two
blocks used in the paper also ship with the package as `morro_bay`.
