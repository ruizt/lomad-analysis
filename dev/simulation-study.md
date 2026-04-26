# Simulation Study Design

## Overview

The simulation study has three parts:

1.  **Part 1 — Calibration under similarity.** Show the method behaves correctly when trends are similar but not identical. This is not a power study — it validates that the method does not falsely flag near-identical trends as decoupled.

2.  **Part 2 — Power study: structured decoupling.** A full factorial experiment using `make_trends_rate`. The periodic structure of decoupling events is easy to understand and establishes the core power results.

3.  **Part 3 — Power study: naturalistic decoupling.** Same factorial using `make_trends_smooth` and/or `make_trends_cross`. More realistic alternative structure; directly comparable to Part 2.

------------------------------------------------------------------------

## Data generation

All simulated datasets consist of a **pair of trend series** with **calibrated ARMA noise** added to each. The full pipeline is:

```         
make_trends_*()  -->  add_noise()  -->  observed series (y1, y2)
```

### Trend generation

Four functions generate paired trend series. All are built on a shared Fourier basis and return the true coupling state `lambda` as ground truth (except `make_trends_dist`).

| Function | Structure | lambda range | Primary parameters |
|------------------|------------------|------------------|------------------|
| `make_trends_dist` | Static — constant separation | none | `d` |
| `make_trends_rate` | Periodic gamma-shaped decoupling events | [0, 1] | `d`, `rate` |
| `make_trends_smooth` | Stochastic, repel only | (0, 1) | `d`, `bw`, `coupling` |
| `make_trends_cross` | Stochastic, with crossings | unbounded | `d`, `bw`, `coupling` |

In all cases, **`d` is the Euclidean distance between the two output series**: `||y1 - y2||`. Lambda is the ground-truth coupling weight at each time point: lambda = 1 means fully coupled (both series track their shared mean); lambda = 0 means fully decoupled (maximum separation); lambda \> 1 means the series have crossed past the mean to opposite sides (`make_trends_cross` only).

### Noise

Noise is added via `add_noise()`, which wraps `make_noise()`. The noise process for each series is ARMA with innovation SD calibrated so that the mean local SNR matches `lambda_target`:

$$\lambda_{\text{target}} = \frac{\tau^2}{\sigma^2}$$

where $\tau^2$ is the smoothed signal variance (rolling variance of the MA-smoothed trend) and $\sigma^2$ is the smoothed noise variance. This calibration is done separately for each series. `add_noise()` also returns diagnostics including `mean_snr` and `noise_var_ratio` (empirical / theoretical smooth noise variance; should be near 1).

------------------------------------------------------------------------

## Part 1: Calibration under similarity

**Generator:** `make_trends_dist` with small `d`

**Goal:** Confirm the method does not flag near-but-not-identical trends as decoupled. Contrast with a test of exact trend identity, which would reject for any `d > 0`.

**Design:** Vary `d` over a small range, e.g. `d ∈ {0, 0.25, 0.5, 1}`, at two SNR levels (high and low). Run `S` replicates per cell. Record the fraction of time points flagged as decoupled — this should remain near zero across all `d` values.

------------------------------------------------------------------------

## Parts 2 and 3: Power studies

**Generators:** `make_trends_rate` (Part 2), `make_trends_smooth` / `make_trends_cross` (Part 3)

**Primary axis:** `d`. Use a **common absolute grid** across all conditions — e.g. `d ∈ {0, 0.5, 1, 2, 3, 4, 5}`. Power curves are allowed to saturate at different `d` values; the saturation point is itself informative (minimum detectable effect size under each condition).

### Factorial conditions

| Factor | Levels | Notes |
|------------------------|------------------------|------------------------|
| `d` | 0, …, max | Includes 0 = null |
| SNR (`lambda_target`) | high (\~4), low (\~0.5) | Shared across both series |
| Autocorrelation (AR1 `phi`) | weak (\~0.3), strong (\~0.8) |  |
| `rate` | sparse (\~2 events), moderate (\~5 events) | Part 2 only |
| `bw` / `coupling` | narrow/wide bandwidth; moderate/high coupling | Part 3 only; fixed at canonical values |
| `n` | 250, 500, 1000 |  |

`S` replicates per cell (target: MC standard error \< 0.05 on power estimates, roughly `S = 100`).

### Outputs per cell

| Output       | Description                                              |
|--------------|----------------------------------------------------------|
| Global power | P(any decoupling detected)                               |
| Sensitivity  | P(flagged \| truly decoupled), using ground-truth lambda |
| FDR          | P(truly coupled \| flagged), using ground-truth lambda   |

A time point is "truly decoupled" if `lambda < 0.5` (lambda below the midpoint of the coupled band).

------------------------------------------------------------------------

## Repository structure for the simulation study

All simulation work lives under `dev/simulation/`. This folder is not part of the package build.

```
dev/simulation/
├── grid.R          # Task 1: factorial grid definition
├── run_one.R       # Task 2: single-replicate pipeline function
├── pilot.R         # Task 0: pilot run script (coarse d grid, hardest/easiest conditions)
├── run_sim.R       # Task 3: single-row simulation script (HPC entrypoint)
├── results/        # saved .rds output files, one per grid row (gitignored)
└── analysis.qmd    # Task 4: aggregation and power curve plots
```

`run_one.R` defines a function; it should be `source()`d at the top of `pilot.R` and `run_sim.R`. `grid.R` defines the grid object and is `source()`d where needed. `pilot.R` uses `devtools::load_all()` and is run interactively in RStudio — it does not need to run on the cluster. `run_sim.R` should use `library(lomad)` since the cluster environment will have the package installed; the batch job design (see `lomad-sim/`) will handle the rest of the deployment details.

The key structural principle for cluster execution is **one job per grid row**: `run_sim.R` runs all `S` replicates for a single row and saves one result file. The row is identified by a command-line argument, which makes the script easy to test locally and straightforward for the batch job layer to parameterize:

```r
args   <- commandArgs(trailingOnly = TRUE)
row_id <- as.integer(args[1])
S      <- as.integer(args[2])
```

To run locally: `Rscript run_sim.R 5 100` (row 5, 100 replicates). The batch job wrapper passes the same arguments at scale.

Result files named by row index (e.g. `row_001.rds`, `row_042.rds`) make it easy to check which rows are complete and avoid long filenames encoding all parameter values.

## Student tasks: power study scaffolding

The primary student contribution is building the infrastructure that makes the power study runnable. This involves five components, intended to be tackled in order: calibrating the `d` grid from a pilot run, setting up the factorial design, writing a single-replicate pipeline, looping over replicates and aggregating, and visualizing results. A validation step (type I error check) should be run before the full study.

### Task 0: Pilot study and d grid calibration

Before running the full factorial, run a small pilot to determine the appropriate `d` grid. The goal is to find the range of `d` where the power curve transitions from near 0 to near 1, so the full study grid covers the interesting region without wasting replicates on values that are trivially powerful or trivially powerless.

**Design:** Run the pipeline at a coarse `d` grid (e.g. `d ∈ {0, 1, 2, 4, 6, 8, 10}`) with few replicates (`S = 30`) and only the **hardest condition**: lowest SNR (`lambda_target = 0.5`), strongest autocorrelation (`phi = 0.8`), smallest `n` (250). This is the condition where power will be lowest — if the curve saturates within this grid the other conditions certainly will too. Also run the **easiest condition** (high SNR, weak autocorrelation, large `n`) to confirm it saturates within the same grid.

**Output:** Rough power curves for the hardest and easiest conditions. From these, identify: - The approximate `d` where power first exceeds \~0.1 (lower bound of interesting range) - The approximate `d` where power reaches \~0.9 (upper bound) - A refined grid of 6–8 points covering that range

Report the chosen grid with a brief justification referencing the pilot curves.

### Task 1: Factorial design setup

Write a function or script that generates the full factorial grid as a data frame, one row per condition cell. For Part 2 (structured decoupling) this looks like:

``` r
expand.grid(
  d             = c(0, 0.5, 1, 2, 3, 4, 5),
  lambda_target = c(0.5, 4),
  phi           = c(0.3, 0.8),
  rate          = c(0.004, 0.01),
  n             = c(250, 500, 1000)
)
```

Each row defines one condition cell. The grid for Part 3 replaces `rate` with `bw` and `coupling`. Confirm the grid is sensible — total number of cells, and expected runtime per cell at `S` replicates — before running anything.

### Task 2: Single-replicate pipeline

Write a function `run_one(d, lambda_target, phi, rate, n, seed)` that runs the full pipeline for one replicate:

1.  Generate a trend pair: `make_trends_rate(n = n, d = d, rate = rate, seed = seed)`
2.  Add calibrated noise: `add_noise(trends, h = ..., lambda_target = lambda_target, order = c(1, 0), ar.coefs = phi)`
3.  Fit the model: `lomad_fit(y1, y2, ...)`
4.  Compute and return the three outputs as a one-row data frame:
    -   **Global power indicator**: 1 if any time point is flagged, 0 otherwise
    -   **Sensitivity**: fraction of truly decoupled time points (lambda \< 0.5) that are flagged
    -   **FDR**: fraction of flagged time points where lambda ≥ 0.5

Returning a one-row data frame makes it easy to collect results with `bind_rows()`.

### Task 3: Single-row simulation script (`run_sim.R`)

Write `run_sim.R` as the HPC entrypoint — it runs all `S` replicates for **one grid row** and saves one result file. It should not loop over the full grid; the cluster job sweep handles that by submitting one job per row.

```r
library(lomad)
source("grid.R")     # loads `grid` data frame
source("run_one.R")  # loads `run_one()` function

args    <- commandArgs(trailingOnly = TRUE)
row_id  <- as.integer(args[1])
S       <- as.integer(args[2])
out_dir <- "results"

row <- grid[row_id, ]

results <- lapply(seq_len(S), \(s)
  run_one(d             = row$d,
          lambda_target = row$lambda_target,
          phi           = row$phi,
          rate          = row$rate,
          n             = row$n,
          seed          = s)
) |> bind_rows()

out_file <- file.path(out_dir, sprintf("row_%03d.rds", row_id))
saveRDS(list(params = row, results = results), out_file)
```

For local aggregation after all jobs complete, load all result files and bind:

```r
files <- list.files("results", pattern = "row_.*\\.rds", full.names = TRUE)
all_results <- lapply(files, \(f) {
  x <- readRDS(f)
  cbind(x$params, x$results)
}) |> bind_rows()
```

Then summarise with `group_by()` / `summarise()` as in Task 1.

### Task 4: Power curve visualization

Write a plotting function that takes the aggregated results and produces power curves against `d`, with `±2 SE` ribbons. Facet or color by condition (SNR, autocorrelation, `n`, or rate) to show how operating conditions affect the curves.

### Validation: type I error check

Before running the full study, run Tasks 2–3 with `d = 0` only (`S = 200` replicates) across all SNR × autocorrelation × `n` combinations. The global power at `d = 0` is the empirical type I error rate — it should be near the nominal level (5%). Produce a table of these rates. If any cell is substantially inflated, investigate before proceeding to the full grid.

------------------------------------------------------------------------

## Open questions

-   Whether Part 3 uses `make_trends_smooth`, `make_trends_cross`, or both.
-   Number of replicates `S` per cell: target MC SE \< 0.05 suggests `S ≈ 100` for power estimates near 0.5; more may be needed at the tails.
-   Whether `bw` and `coupling` in Part 3 should be chosen to approximately match the expected fraction of time decoupled in Part 2 (for cleaner comparisons).
