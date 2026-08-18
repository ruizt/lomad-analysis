# Validation Study Design

## Experiments

One experiment, reported in two parts:

1. **Moment accuracy** checks Proposition 1's expressions for the limiting mean
    $\rho_t$ and variance $V_t$ against their empirical counterparts

2. **Calibration** checks that the full `lomad_fit()` + `lomad_test()` pipeline
    preserves approximation accuracy despite estimation error, against an
    oracle standardization using the true $\rho_t$ and $V_t$

### Data generation

A single random Fourier trend is drawn, and the second series is a fixed
affine map of it: $\nu_2 = 2\nu_1$. Local affine similarity therefore holds
everywhere by construction, while trend variances are related by 
$\tau_2^2 = 4\tau_1^2$.

Fixed parameters throughout both experiments are:

| Parameter | Value |
|-----------|-------|
| Series length *n* | 2000 |
| MA smoothing window *h* | 5 |
| Correlation window *s* | 100 |
| Trend | `sim_trends(n = 2000, d = 0, nb = 25, sd0 = 50, p = 1.5, seed = 5381)` |
| Affine map | $\nu_2 = 2\nu_1$ (`affine_a = 0`, `affine_b = 2`) |

Noise is independent AR(1). Here $\phi_k, \sigma^2_k$ are the AR coefficient
and innovation variance for series $k$, and $\sigma^2_{\eta_k}$ is the smoothed
noise variance that enters the CLT.

| Series $k$ | $\phi_k$ | $\sigma^2_k$ | $\sigma^2_{\eta_k}$ |
|-----------|----------|--------------|----------------------|
| 1 | 0.5 | 0.64 | 0.380 |
| 2 | 0.3 | 0.64 | 0.227 |

Approximations are checked at t ∈ {400, 900, 1400, 1900}, plus a dense grid at
spacing 5 over the full testable range. The dense grid carries the coverage
band and the moments of $R_t$ used for the $\rho_t$ and $V_t$ comparisons.
$\bar R_t$ is formed from the first 100 replicates only (`E2E_RHO_REPS`): over
all 1000 its Monte Carlo error is thinner than the plotted line.

### Job table

A single job:

| Experiment ID | *s* | Reps | Output |
|---------------|-----|------|--------|
| `e2e-s100` | 100 | 1000 | Oracle + pipeline quantities at 4 eval pts, and moments of R_t over the dense grid |

---

### Expected outputs

Per job (Tide) or per experiment (local):

- `{experiment}.rds` — one file per experiment in `results/_raw/`, named by the
  `SIM_EXPERIMENT` label (`e2e-s100`)

Aggregated (after `collect-results.R`):

- `results/simulations-validation-results.rds` — named list of all experiment
  objects, keyed by experiment label

------------------------------------------------------------------------

## File layout

```
validation/
├── design.md           ← you are here
├── simulation-template.R          ← local proof-of-concept, mirrors tide/sim.R
├── collect-results.R           ← assembles fetched per-job files into a compiled results object
├── results/
│   ├── simulations-validation-results.rds  ← tracked
│   └── _raw/           ← per-job .rds files fetched from Tide, Zenodo only
└── tide/               ← Kubernetes scaffolding
    ├── sim.R           ← container entrypoint (dispatches on SIM_EXPERIMENT)
    ├── submit.sh       ← full pipeline (PVC → jobs → wait → fetch)
    ├── submit_sweep.sh ← submits one Job per experiment
    ├── pvc.yaml        ← shared storage (create once)
    ├── accessor.yaml   ← lightweight pod for file retrieval
    └── fetch.sh        ← copies results from PVC to local machine
```

`_`-prefixed entries are gitignored and archived on Zenodo; the compiled `.rds`
files are small enough to track, so a fresh clone can rebuild every figure
without downloading anything.

`tide/` holds everything that talks to the cluster. `collect-results.R` sits
outside it because it only reads local files.

---

## Workflow

### Local testing

Source `simulation-template.R`; it runs the experiment at *S* = 50 and draws
draft versions of both panels. `tide/sim.R` is the source of truth: change it
first, then mirror the change into the template.

Before submitting at scale, test the container entrypoint locally with a small
number of replicates:

```bash
SIM_EXPERIMENT=e2e-s100 SIM_S=5 SIM_SEED=7291 \
  SIM_OUT_DIR=simulations/validation/results/_raw \
  Rscript simulations/validation/tide/sim.R
```

### Deployment

All studies share a single Docker image `ghcr.io/ruizt/lomad-simulations:latest`.
The image contains R and all packages (including lomad) but no simulation
scripts — `sim.R` is mounted at runtime via a Kubernetes ConfigMap.

A wrapper around the simulation steps provides a one-shot submission: 

```bash
bash simulations/validation/tide/submit.sh
```

This creates the PVC, submits the job, polls until all
complete, and fetches results to `results/_raw/`.

This can also be executed step by step:

```bash
# 1. Create PVC (once)
kubectl apply -n cal-poly-ruiz -f simulations/validation/tide/pvc.yaml

# 2. Submit all experiments (creates ConfigMap + Jobs)
bash simulations/validation/tide/submit_sweep.sh

# 3. Monitor until all jobs show Complete
kubectl get jobs -n cal-poly-ruiz -l app=lomad-valid

# 4. Fetch results from PVC to local
bash simulations/validation/tide/fetch.sh

# 5. Assemble results
Rscript simulations/validation/collect-results.R   # compiled results

# 6. Clean up
kubectl delete jobs -n cal-poly-ruiz -l app=lomad-valid
kubectl delete -n cal-poly-ruiz -f simulations/validation/tide/pvc.yaml
```
