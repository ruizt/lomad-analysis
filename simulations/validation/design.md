# Validation Study Design

## Experiments

Two validation experiments:

1. **Oracle experiment** checks that the asymptotic normal approximation for the 
    standardized rolling correlation is accurate at realistic sample sizes

2. **End-to-end experiment** checks that the full `lomad_fit()` + `lomad_test()` 
    pipeline preserves approximation accuracy despite estimation error

### Data generation

A single random Fourier trend is drawn, and the second series is a fixed
affine map of it: $\nu_2 = 2\nu_1$. Local affine similarity therefore holds
everywhere by construction, while trend variances are related by 
$\tau_2^2 = 4\tau_1^2$.

Fixed parameters throughout both experiments are:

| Parameter | Value |
|-----------|-------|
| Series length *n* | 2000 |
| MA smoothing window *h* | 20 |
| Trend | `sim_trends(n = 2000, d = 0, nb = 25, sd0 = 50, p = 1.5, seed = 5381)` |
| Affine map | $\nu_2 = 2\nu_1$ (`affine_a = 0`, `affine_b = 2`) |

Noise is independent ARMA, fixed within each experiment and differing between
them: ARMA(1,1) for the oracle experiment, AR(1) for the end-to-end. Here
$\phi_k, \theta_k, \sigma^2_k$ are the AR coefficient, MA coefficient and
innovation variance for series $k$, and $\sigma^2_{\eta_k}$ is the smoothed
noise variance that enters the CLT.

| Experiment | Series $k$ | $\phi_k$ | $\theta_k$ | $\sigma^2_k$ | $\sigma^2_{\eta_k}$ |
|------------|-----------|----------|------------|--------------|----------------------|
| Oracle | 1 | 0.6 | 0.3 | 0.64 | 0.305 |
| Oracle | 2 | 0.4 | −0.2 | 1.00 | 0.086 |
| End-to-end | 1 | 0.5 | — | 0.64 | 0.119 |
| End-to-end | 2 | 0.3 | — | 0.64 | 0.063 |

Each experiment also checks approximations at specific evaluation points:

- CLT experiments: t = 600 (where ρ ≈ 0.24–0.35 across s values)
- E2E experiment: t ∈ {400, 900, 1400, 1900}, plus a dense grid at
  spacing 5 over the full testable range for the coverage band

### Job table

Experiments are distributed across six jobs:

| Experiment ID | *s* | Reps | Output |
|---------------|-----|------|--------|
| `clt-s80` | 80 | 2000 | R_t at t = 600 per rep |
| `clt-s150` | 150 | 2000 | R_t at t = 600 per rep |
| `clt-s300` | 300 | 2000 | R_t at t = 600 per rep |
| `rho-s150` | 150 | 500 | Running mean R̄_t over time |
| `var-s150` | 150 | 2000 | R_t at every 20th eval point |
| `e2e-s150` | 150 | 1000 | Oracle + pipeline quantities at 4 eval pts |

---

### Expected outputs

Per job (Tide) or per experiment (local):

- `{experiment}.rds` — one file per experiment in `results/_raw/`, named by the
  `SIM_EXPERIMENT` label (`clt-s80`, `clt-s150`, `clt-s300`, `rho-s150`,
  `var-s150`, `e2e-s150`)

Aggregated (after `collect-results.R`):

- `results/simulations-validation-results.rds` — named list of all experiment
  objects, keyed by experiment label

------------------------------------------------------------------------

## File layout

```
validation/
├── design.md           ← you are here
├── simulation-template.R          ← local proof-of-concept (defines run_rep_*())
├── collect-results.R           ← assembles fetched per-job files into a compiled results object
├── results/
│   ├── simulations-validation-results.rds  ← tracked
│   ├── _simulations-validation.zip         ← archive of _raw/, Zenodo only
│   └── _raw/           ← per-job .rds files fetched from Tide, Zenodo only
└── tide/               ← Kubernetes scaffolding
    ├── sim.R           ← container entrypoint (dispatches on SIM_EXPERIMENT)
    ├── submit.sh       ← full pipeline (PVC → jobs → wait → fetch)
    ├── submit_sweep.sh ← submits one Job per experiment
    ├── pvc.yaml        ← shared storage (create once)
    ├── accessor.yaml   ← lightweight pod for file retrieval
    └── fetch.sh        ← copies results from PVC to local machine
```

---

## Workflow

### Local testing

Source `simulation-template.R`; it runs all experiments at small scale (reduced
reps). The `run_rep_*()` functions it defines are the same ones used in
`tide/sim.R`.

Before submitting at scale, test the container entrypoint locally with a small
number of replicates:

```bash
SIM_EXPERIMENT=clt-s80 SIM_S=5 SIM_SEED=7291 \
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

This creates the PVC, submits one job per experiment (6 jobs), polls until all
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
