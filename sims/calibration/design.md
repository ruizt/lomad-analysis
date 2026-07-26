# Calibration Study Design

## Objective

Verify that the CLT-based test controls type I error across the null region — i.e., for small *d* where trends remain locally similar — while a full-oracle identity test that targets global separation rejects quickly as *d* grows. Type I error control at *d* = 0 is necessary but not sufficient; the test should stay near the nominal level for all *d* small enough that local windows still reflect co-moving trends.

## Methods compared

-   **CLT test (estimated)**: the full estimation pipeline — `lomad_fit()` estimates the noise ARMA process and smoothed-noise ACVF, then `lomad_test()` applies Benjamini–Yekutieli FDR correction pointwise. The per-dataset rejection indicator is `any(rejected[valid_idx])`.

-   **CLT test (oracle)**: same as above but with true per-series AR(1) coefficient and innovation variance supplied via `noise_override`, bypassing noise estimation. Isolates whether the CLT framework itself is correctly calibrated from any estimation error.

-   **Identity test (oracle)**: tests H₀: *d* = 0 via a global excess-variance statistic on the MA(*h*)-smoothed difference, using the *true* AR(1) coefficient and innovation variance. The test statistic is *T* = (mean(*D*²) − Var(*D*)) / se(*T*), where se(*T*) accounts for autocorrelation in *D*² via the smoothed ACVF. Well-calibrated at *d* = 0 and saturates quickly as *d* grows. Included as a contrast to show the CLT test is selective.

## Data-generating process

Trends from `sim_trends(method = "dist")` (Fourier basis, target L² separation *d*), with AR(1) noise added via `sim_noise_pair()` at a fixed SNR.

## Parameters

| Parameter              | Value                   |
|------------------------|-------------------------|
| Series length *T*      | 1000                    |
| AR(1) coefficient *φ*  | 0.5                     |
| Target SNR *λ*         | 1                       |
| Smoothing window *h*   | 10                      |
| Correlation window *s* | 50                      |
| *d* grid               | 0, 0.2, 0.5, 1          |
| Replicates *S*         | 50 (local) → 500 (Tide) |
| Significance level *α* | 0.05                    |

The *d* grid spans the null (*d* = 0) through moderate separation (*d* = 1). The identity test serves as a contrast: it saturates quickly as *d* grows, while the CLT test remains well-calibrated.

## Estimand

For each (d, replicate) pair: does the test produce at least one rejection?

-   **CLT tests (estimated and oracle)**: `any(tst$rejected[fit$valid_idx], na.rm = TRUE)`
-   **Identity test**: `global_p < α`

The primary output is **rejection rate** (mean of 0/1 across *S* replicates) as a function of *d* for all three methods.

## Acceptance criterion

-   CLT rejection rate ≤ 0.10 at *d* = 0 (type I error controlled near nominal).
-   CLT rejection rate ≤ 0.15 for all *d* in the grid (null region for similarity).
-   Identity test rejection rate ≥ 0.50 at *d* = 0.5 (confirms the contrast is informative).

------------------------------------------------------------------------

## File layout

```         
calibration/
├── design.md           ← you are here
├── template.R          ← local proof-of-concept (defines run_rep())
├── settings.R          ← visual walkthrough of single replicates
├── results/            ← output (gitignored)
│   ├── results.rds
│   ├── results_summary.rds
│   ├── power_curve.png
│   └── raw/            ← per-job .rds files fetched from Tide
└── tide/               ← Kubernetes scaffolding
    ├── sim.R           ← simulation script (mounted into container)
    ├── submit.sh       ← runs the full pipeline (PVC → jobs → wait → fetch)
    ├── submit_sweep.sh ← submits one Job per d (called by submit.sh)
    ├── pvc.yaml        ← shared storage (create once)
    ├── job.yaml        ← Job template (reference only)
    ├── accessor.yaml   ← lightweight pod for file retrieval
    ├── fetch.sh        ← copies results from PVC to local machine
    └── collect.R       ← assembles per-job files into summary
```

------------------------------------------------------------------------

## Workflow

### Local development

Source `template.R` directly in RStudio. It loops over all *d* values, prints a summary table, and exports results to `results/`. The parameter `S` at the top of the script controls the number of replicates (default 50). See `settings.R` for a visual walkthrough of a single replicate at each *d* value.

To inspect a single replicate interactively:

``` r
run_rep(d = 0, seed = 12345)
```

### Writing `tide/sim.R`

The entrypoint reads parameters from environment variables and must do three things:

1.  **Copy `run_rep()` from `template.R`** — it depends only on the fixed parameters and the lomad package, so it works unchanged.

2.  **Run the simulation loop** — draw *S* seeds and call `run_rep(d, seed)` for each. Collect results into a single data frame. Use `set.seed(seed0 + as.integer(d * 100))` before drawing seeds so each *d* value gets a distinct random stream.

3.  **Save results** — write `list(d, S, seed0, results)` to `file.path(out_dir, <filename>.rds)`. See the hint in `tide/sim.R` for how to derive the filename from *d*.

### Testing locally before submitting

Once `tide/sim.R` is complete, test it with a small *S*:

``` bash
SIM_D=0 SIM_S=5 SIM_SEED=4853 SIM_OUT_DIR=sims/calibration/results/raw \
  Rscript sims/calibration/tide/sim.R
```

### Shared container image

All studies share a single Docker image that is pre-built and publicly available:

```         
ghcr.io/ruizt/lomad-sims:latest
```

The image contains R and all packages (including lomad) but no simulation scripts — `sim.R` is mounted into the container at `/scripts/sim.R` at runtime via a Kubernetes ConfigMap. See `../tide-example/building-containers.md` for details on building and updating the image.

### Running on Tide

#### One-shot pipeline

``` bash
bash sims/calibration/tide/submit.sh
```

This creates the PVC, submits one job per *d* value (with the ConfigMap), polls until all jobs complete, and fetches results to `results/raw/`.

#### Step-by-step

``` bash
# 1. Create PVC (once)
kubectl apply -n cal-poly-ruiz -f sims/calibration/tide/pvc.yaml

# 2. Submit all d values (creates ConfigMap + Jobs)
bash sims/calibration/tide/submit_sweep.sh

# 3. Monitor until all jobs show Complete
kubectl get jobs -n cal-poly-ruiz -l app=lomad-calib

# 4. Fetch results from PVC to local
bash sims/calibration/tide/fetch.sh

# 5. Assemble summary and plot
Rscript sims/calibration/tide/collect.R

# 6. Clean up
kubectl delete jobs -n cal-poly-ruiz -l app=lomad-calib
kubectl delete -n cal-poly-ruiz -f sims/calibration/tide/pvc.yaml
```

#### Monitoring and debugging

``` bash
# Job status
kubectl get jobs -n cal-poly-ruiz -l app=lomad-calib

# Container logs for a specific d
kubectl logs -n cal-poly-ruiz job/lomad-calib-d0-5

# Detailed events and error info
kubectl describe job -n cal-poly-ruiz lomad-calib-d0-5
```

#### Tips

-   Start with `SIM_S=20` in `submit_sweep.sh` to confirm everything works end-to-end, then delete the jobs and resubmit at 200–500.
-   If you update `sim.R`, re-running `submit_sweep.sh` will update the ConfigMap automatically — no image rebuild needed.
-   Results are saved to `results/results.rds` (all replicates) and `results/results_summary.rds` (rejection rates by *d*).
