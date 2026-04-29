# Calibration Study Design

## Objective

Verify that the CLT-based test controls type I error across the null region — i.e., for small *d* where trends remain locally similar — while a full-oracle identity test that targets global separation rejects quickly as *d* grows. Type I error control at *d* = 0 is necessary but not sufficient; the test should stay near the nominal level for all *d* small enough that local windows still reflect co-moving trends.

## Methods compared

-   **CLT test (estimated)**: the full estimation pipeline — `lomad_fit()` estimates the noise ARMA process and smoothed-noise ACVF, then `lomad_test()` applies Benjamini–Yekutieli FDR correction pointwise. The per-dataset rejection indicator is `any(rejected[valid_idx])`.

-   **Identity test (oracle)**: tests H₀: *d* = 0 via the L² norm of the MA(*h*)-smoothed difference, using the *true* AR(1) coefficient and innovation variance. Well-calibrated at *d* = 0 and monotone in *d*, but detects any global separation regardless of local co-movement. Included as a contrast to show the CLT test is selective.

## Data-generating process

Trends from `sim_trends(method = "dist")` (Fourier basis, target L² separation *d*), with AR(1) noise added via `sim_noise_pair()` at a fixed SNR.

## Parameters

| Parameter              | Value                         |
|------------------------|-------------------------------|
| Series length *T*      | 1000                          |
| AR(1) coefficient *φ*  | 0.5                           |
| Target SNR *λ*         | 1.5                           |
| Smoothing window *h*   | auto (`max(5, floor(T/200))`) |
| Correlation window *s* | auto (`min(60h, floor(T/4))`) |
| ARMA order bound       | `max_pq = 3`                  |
| *d* grid               | 0, 0.1, 0.2, 0.3, 0.4, 0.5    |
| Replicates *S*         | 10 (local) → 100–500 (Tide)   |
| Significance level *α* | 0.05                          |

The *d* grid focuses on the null region (d ≤ 0.5) where the CLT approximation is expected to hold. The identity test serves as a reference showing that non-zero *d* is detectable in principle.

## Estimand

For each (d, replicate) pair: does the test produce at least one rejection?

-   **CLT test**: `any(tst$rejected[fit$valid_idx], na.rm = TRUE)`
-   **Identity test**: single p-value \< α

The primary output is **rejection rate** (mean of 0/1 across *S* replicates) as a function of *d* for both methods.

## Acceptance criterion

-   CLT rejection rate ≤ 0.10 at *d* = 0 (type I error controlled near nominal).
-   CLT rejection rate ≤ 0.15 for all *d* ≤ 0.5 (null region).
-   Identity test rejection rate ≥ 0.50 at *d* = 0.5 (confirms the contrast is informative even in the null region).

------------------------------------------------------------------------

## Workflow

### Local development

Source `run-small.R` directly in RStudio. It loops over all *d* values, prints a summary table, and produces a rejection-rate plot. The parameter `S` at the top of the script controls the number of replicates; the default (10) is for fast iteration. Increase to 30–50 before concluding anything.

To inspect a single replicate interactively, source the file and then call:

``` r
run_rep(d = 0, seed = 12345)
```

------------------------------------------------------------------------

## Student assignment: scaling up on Tide

Tide runs Kubernetes. The goal is to parallelise the *d* sweep: one container per *d* value, each running *S* = 100–500 replicates, all jobs submitted simultaneously. The files in `tide/` provide the scaffolding.

### Files provided

| File | Status | Purpose |
|----|----|----|
| `tide/Dockerfile` | complete | Builds the container image from the lomad source |
| `tide/job.yaml` | complete | Single-job template (reference only) |
| `tide/submit_sweep.sh` | complete | Submits one job per *d* value |
| `tide/collect.R` | complete | Assembles per-job `.rds` files after completion |
| `tide/sim.R` | **your task** | Container entrypoint — one *d* per run |

### Your task: write `tide/sim.R`

The entrypoint must do four things:

1.  **Read parameters from environment variables.** The skeleton at the top of `tide/sim.R` already does this — do not change it.

2.  **Copy the `run_rep` function from `run-small.R`** into `tide/sim.R`. It depends only on the fixed parameters and the `lomad` package, so it works unchanged.

3.  **Run the simulation loop.** Draw *S* seeds and call `run_rep(d, seed)` for each. Collect results into a single data frame. Use `set.seed(seed0 + as.integer(d * 100))` before drawing seeds so each *d* value gets a distinct random stream.

4.  **Save results.** Write a list `(d, S, seed0, results)` to `file.path(out_dir, <filename>.rds)`. See the hint in `tide/sim.R` for how to derive the filename from *d*.

### Testing locally before building the image

Once `tide/sim.R` is complete, test it with a small *S* before building the container:

``` bash
SIM_D=0 SIM_S=5 SIM_SEED=4853 Rscript dev/sims/calibration/tide/sim.R
```

This runs five replicates at *d* = 0 and writes output to `/jobs/output/` (or `SIM_OUT_DIR` if you override it). Check that the `.rds` file appears and the rejection rates printed to the log look sensible.

### Building and pushing the container image

From the repo root (substitute your GitHub org):

``` bash
docker buildx build --platform linux/amd64 \
  -f dev/sims/calibration/tide/Dockerfile \
  -t ghcr.io/<org>/lomad-calib:latest --push .
```

Update the `IMAGE` variable in `submit_sweep.sh` to match.

### Submitting the sweep

Edit `SIM_S` in `submit_sweep.sh` before submitting. Start with 20 to confirm
everything works end-to-end, then delete the jobs and resubmit at 200–500.

```bash
bash dev/sims/calibration/tide/submit_sweep.sh
```

This submits six jobs in parallel (one per *d* value). Monitor progress:

```bash
kubectl get jobs -l app=lomad-calib          # completion status
kubectl logs job/lomad-calib-d0-3            # stdout for d = 0.3
```

Wait until all jobs show `Complete` before fetching results.

### Fetching results

Run `fetch.sh` from the repo root. It spins up a temporary accessor pod,
copies all `.rds` files from the PVC to `results/raw/` locally, then tears
the pod down automatically:

```bash
bash dev/sims/calibration/tide/fetch.sh
```

You should see one `.rds` file per *d* value in `dev/sims/calibration/results/raw/`.

### Assembling the summary

```bash
Rscript dev/sims/calibration/tide/collect.R
```

This reads `results/raw/*.rds`, prints rejection rates, and writes
`results/summary.rds` and `results/power_curve.png`.

### Full sequence at a glance

```bash
# 1. Build and push the image (once per code change)
docker buildx build --platform linux/amd64 \
  -f dev/sims/calibration/tide/Dockerfile \
  -t ghcr.io/<org>/lomad-calib:latest --push .

# 2. Submit all d values in parallel
bash dev/sims/calibration/tide/submit_sweep.sh

# 3. Monitor until all jobs show Complete
kubectl get jobs -l app=lomad-calib

# 4. Copy results from PVC to local
bash dev/sims/calibration/tide/fetch.sh

# 5. Assemble summary and plot
Rscript dev/sims/calibration/tide/collect.R
```
