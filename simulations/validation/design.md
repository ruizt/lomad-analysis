# Validation Study Design

## Objective

Empirically verify two things:

1. **Oracle CLT**: The asymptotic normal approximation for the standardized
   rolling correlation $\sqrt{s}(R_t - \rho_t)/\sqrt{V_t} \approx N(0,1)$ is
   accurate at practical sample sizes, with both the mean ($\rho_t$) and
   variance ($V_t$) formulas from Proposition 1 matching their empirical
   counterparts.

2. **Estimation layer**: The full `lomad_fit()` + `lomad_test()` pipeline —
   which estimates trends, noise, $\hat\rho$, and $\hat V$ — preserves CLT
   coverage at realistic sample sizes despite estimation error.

## Experiments

Six independent jobs, producing one composite three-row figure. Each job is
identified by a `SIM_EXPERIMENT` environment variable and produces one `.rds`
file.

### Data generation

A single random Fourier trend is drawn, and the second series is a fixed
affine map of it: $\nu_2 = 2\nu_1$. Local affine similarity therefore holds
everywhere by construction, while trend variances are related by 
$\tau_2^2 = 4\tau_1^2$.

| Parameter | Value |
|-----------|-------|
| Series length *n* | 2000 |
| MA smoothing window *h* | 20 |
| Trend | `sim_trends(n = 2000, d = 0, nb = 25, sd0 = 50, p = 1.5, seed = 5381)` [^1] |
| Affine map | $\nu_2 = 2\nu_1$ (`affine_a = 0`, `affine_b = 2`) |

The data-generating processes for noise are as follows
| Noise (series 1), oracle | ARMA(1,1): φ = 0.6, θ = 0.3, σ_ε = 0.8 |
| Noise (series 2), oracle | ARMA(1,1): φ = 0.4, θ = −0.2, σ_ε = 1.0 |
| Noise (series 1), e2e | AR(1): φ = 0.5, σ_ε = 0.8 |
| Noise (series 2), e2e | AR(1): φ = 0.3, σ_ε = 0.8 |

The oracle ARMA(1,1) noise is used for rows (a)–(b) (where all parameters are
known). The AR(1) noise is used for row (c), where the pipeline estimates
AR(1) — using a matched DGP tests the CLT machinery without confounding it with
model misspecification.

Eval points:
- CLT experiments (row a): t = 600 (where ρ ≈ 0.24–0.35 across s values)
- E2E experiment (row c): t ∈ {400, 900, 1400, 1900}, plus a dense grid at
  spacing 5 over the full testable range for the coverage band

### Job table

| Experiment ID | Figure | *s* | Reps | Output |
|---------------|--------|-----|------|--------|
| `clt-s80` | 1 | 80 | 2000 | R_t at t = 600 per rep |
| `clt-s150` | 1 | 150 | 2000 | R_t at t = 600 per rep |
| `clt-s300` | 1 | 300 | 2000 | R_t at t = 600 per rep |
| `rho-s150` | 2 (left) | 150 | 500 | Running mean R̄_t over time |
| `var-s150` | 2 (right) | 150 | 2000 | R_t at grid of eval points |
| `e2e-s150` | 3 | 150 | 1000 | Oracle + pipeline quantities at 5 eval pts |

---

## Expected outputs

### Per job (Tide) or per experiment (local)

- `{experiment}.rds` — one file per experiment in `results/_raw/`, named by the
  `SIM_EXPERIMENT` label (`clt-s80`, `clt-s150`, `clt-s300`, `rho-s80`,
  `rho-s150`, `rho-s250`, `var-s150`, `var-s200`, `e2e-s150`)

### Aggregated (after `collect-results.R`)

- `results/simulations-validation-results.rds` — named list of all experiment
  objects, keyed by experiment label

### Figures

Built by `simulations/simulation-results.R`, not by this study, and written to
`simulations/_img/fig-validation.png`.

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

Source `simulation-template.R` in RStudio. It runs all experiments at small scale (reduced
reps) and produces draft versions of all three figures. The `run_rep_*()`
functions it defines are the same ones used in `tide/sim.R`.

Before submitting to Tide, test the container entrypoint locally with a small
number of replicates:

```bash
SIM_EXPERIMENT=clt-s80 SIM_S=5 SIM_SEED=7291 \
  SIM_OUT_DIR=simulations/validation/results/_raw \
  Rscript simulations/validation/tide/sim.R
```

### Shared container image

All studies share a single Docker image:

```
ghcr.io/ruizt/lomad-simulations:latest
```

The image contains R and all packages (including lomad) but no simulation
scripts — `sim.R` is mounted at runtime via a Kubernetes ConfigMap. See
`../README.md` for details on building/updating.

### Running on Tide

#### One-shot pipeline

```bash
bash simulations/validation/tide/submit.sh
```

Creates the PVC, submits one job per experiment (6 jobs), polls until all
complete, and fetches results to `results/_raw/`.

#### Step-by-step

```bash
# 1. Create PVC (once)
kubectl apply -n cal-poly-ruiz -f simulations/validation/tide/pvc.yaml

# 2. Submit all experiments (creates ConfigMap + Jobs)
bash simulations/validation/tide/submit_sweep.sh

# 3. Monitor until all jobs show Complete
kubectl get jobs -n cal-poly-ruiz -l app=lomad-valid

# 4. Fetch results from PVC to local
bash simulations/validation/tide/fetch.sh

# 5. Assemble results, then build the figure
Rscript simulations/validation/collect-results.R   # compiled results
Rscript simulations/simulation-results.R       # composite figure

# 6. Clean up
kubectl delete jobs -n cal-poly-ruiz -l app=lomad-valid
kubectl delete -n cal-poly-ruiz -f simulations/validation/tide/pvc.yaml
```

#### Monitoring and debugging

```bash
# Job status
kubectl get jobs -n cal-poly-ruiz -l app=lomad-valid

# Container logs for a specific experiment
kubectl logs -n cal-poly-ruiz job/lomad-valid-clt-s80

# Detailed events and error info
kubectl describe job -n cal-poly-ruiz lomad-valid-clt-s80
```

#### Tips

- Start with small `SIM_S` (e.g., 20) in `submit_sweep.sh` to confirm
  everything works end-to-end, then delete jobs and resubmit at full scale.
- If you update `sim.R`, re-running `submit_sweep.sh` will update the
  ConfigMap automatically — no image rebuild needed.
