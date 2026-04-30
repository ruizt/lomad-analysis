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

### Common DGP (shared across all experiments)

Trends are generated using the package function `sim_trends(n = 2000, d = 0,
nb = 25, sd0 = 50, p = 1.5, seed = 5381)`, producing a shared random Fourier
trend with heterogeneous local variance.

| Parameter | Value |
|-----------|-------|
| Series length *n* | 2000 |
| MA smoothing window *h* | 20 |
| Trend | `sim_trends(d = 0, nb = 25, sd0 = 50, p = 1.5, seed = 5381)` [^1] |

[^1]: The spectral decay `p = 1.5` concentrates energy in low frequencies,
producing a heterogeneous trend with long near-flat stretches where τ² is close
to the noise floor. This creates a mix of easy windows (moderate ρ) and hard
windows (ρ ≈ 0 or near ceiling) — a realistic but challenging regime for the
CLT. Smaller `p` (e.g., 1.0) raises τ² but pushes ρ uniformly high, where the
bounded-correlation skew dominates. Testing at `p = 1.0` and `p = 1.25` showed
no meaningful improvement in the Z-score diagnostics; the dominant factor is the
oracle τ² computation, not the spectral decay.
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
- E2E experiment (row c): t ∈ {500, 850, 1000, 1400, 1600} (ρ spans 0.02–0.66)

### Job table

| Experiment ID | Figure | *s* | Reps | Output |
|---------------|--------|-----|------|--------|
| `clt-s80` | 1 | 80 | 2000 | R_t at t = 600 per rep |
| `clt-s150` | 1 | 150 | 2000 | R_t at t = 600 per rep |
| `clt-s300` | 1 | 300 | 2000 | R_t at t = 600 per rep |
| `rho-s150` | 2 (left) | 150 | 500 | Running mean R̄_t over time |
| `var-s150` | 2 (right) | 150 | 2000 | R_t at grid of eval points |
| `e2e-s150` | 3 | 150 | 1000 | Oracle + pipeline quantities at 5 eval pts |

## Composite figure (`validation-composite.png`)

A single three-row figure saved as both PNG (300 dpi) and PDF.

### Row (a) — CLT QQ plots (oracle, faceted by *s*)

Three-panel QQ plot of $\sqrt{s}(R_t - \rho_t)/\sqrt{V_t}$ vs $N(0,1)$, one
panel per $s \in \{80, 150, 300\}$, evaluated at $t = 600$. Black
semi-transparent points, firebrick reference line.

**Story**: The CLT approximation is already adequate at $s = 80$ and tightens
with increasing window size.

### Row (b) — Proposition 1 moment accuracy (oracle, $s = 150$)

- **Left panel** ($\rho$): Theoretical $\rho_t$ (firebrick) vs empirical
  $\bar{R}_t$ (grey30) over time at $s = 150$.
- **Right panel** ($V$): Scatter of theoretical $V_t$ vs empirical
  $s \cdot \text{Var}(R_t)$ across time points at $s = 150$. Firebrick 45°
  reference line, black semi-transparent points.

**Story**: Both moment formulas closely track their empirical counterparts at
the window size used in the end-to-end experiment.

### Row (c) — End-to-end pipeline validation ($s = 150$)

- **Left panel**: Overlaid oracle (blue) and pipeline (orange) QQ plots at
  $t = 1000$.
- **Right panel**: Pointwise 95% coverage at 5 eval points, oracle vs
  pipeline, with $\pm 1.96$ SE error bars and firebrick dashed nominal line.

**Story**: Coverage remains at or above nominal levels; the pipeline produces
somewhat conservative inference due to plug-in variance compression.

## Estimands

- **CLT experiments** (`clt-*`): The standardized statistic $Z_t$ at
  $t = 1000$. Primary check: $Z \sim N(0,1)$.
- **ρ experiments** (`rho-*`): The time-averaged empirical mean $\bar R_t$
  across replications vs theoretical $\rho_t$.
- **V experiment** (`var-s150`): The empirical $s \cdot \text{Var}(R_t)$ at
  each time point vs theoretical $V_t$.
- **End-to-end** (`e2e-s150`): Coverage of $|Z| < 1.96$ under oracle vs
  pipeline standardization at each eval point.

## Acceptance criteria

- CLT QQ plots: points fall within 95% simulation envelope for all three $s$.
- Empirical 95% coverage ≥ 0.90 at $s = 80$ and ≥ 0.93 at $s = 300$.
- $\text{Cor}(\bar R_t, \rho_t) > 0.99$ for both $s$ values.
- $\text{Cor}(V_\text{emp}, V_\text{theory}) > 0.90$.
- Pipeline coverage within ± 0.03 of oracle coverage at each eval point.

---

## File layout

```
validation/
├── design.md           ← you are here
├── template.R          ← local proof-of-concept (defines run_rep_*())
├── results/            ← output (gitignored)
│   ├── raw/            ← per-job .rds files fetched from Tide
│   ├── validation-composite.png
│   └── validation-composite.pdf
└── tide/               ← Kubernetes scaffolding
    ├── sim.R           ← container entrypoint (dispatches on SIM_EXPERIMENT)
    ├── submit.sh       ← full pipeline (PVC → jobs → wait → fetch)
    ├── submit_sweep.sh ← submits one Job per experiment
    ├── pvc.yaml        ← shared storage (create once)
    ├── accessor.yaml   ← lightweight pod for file retrieval
    ├── fetch.sh        ← copies results from PVC to local machine
    └── collect.R       ← assembles per-job files + generates figures
```

---

## Workflow

### Local development

Source `template.R` in RStudio. It runs all experiments at small scale (reduced
reps) and produces draft versions of all three figures. The `run_rep_*()`
functions it defines are the same ones used in `tide/sim.R`.

### Testing tide/sim.R locally

Before submitting to Tide, test the container entrypoint locally with a small
number of replicates:

```bash
SIM_EXPERIMENT=clt-s80 SIM_S=5 SIM_SEED=7291 \
  SIM_OUT_DIR=dev/sims/validation/results/raw \
  Rscript dev/sims/validation/tide/sim.R
```

### Shared container image

All studies share a single Docker image:

```
ghcr.io/ruizt/lomad-sims:latest
```

The image contains R and all packages (including lomad) but no simulation
scripts — `sim.R` is mounted at runtime via a Kubernetes ConfigMap. See
`../tide-example/building-containers.md` for details on building/updating.

### Running on Tide

#### One-shot pipeline

```bash
bash dev/sims/validation/tide/submit.sh
```

Creates the PVC, submits one job per experiment (6 jobs), polls until all
complete, and fetches results to `results/raw/`.

#### Step-by-step

```bash
# 1. Create PVC (once)
kubectl apply -n cal-poly-lomad -f dev/sims/validation/tide/pvc.yaml

# 2. Submit all experiments (creates ConfigMap + Jobs)
bash dev/sims/validation/tide/submit_sweep.sh

# 3. Monitor until all jobs show Complete
kubectl get jobs -n cal-poly-lomad -l app=lomad-valid

# 4. Fetch results from PVC to local
bash dev/sims/validation/tide/fetch.sh

# 5. Assemble results and generate figures
Rscript dev/sims/validation/tide/collect.R

# 6. Clean up
kubectl delete jobs -n cal-poly-lomad -l app=lomad-valid
kubectl delete -n cal-poly-lomad -f dev/sims/validation/tide/pvc.yaml
```

#### Monitoring and debugging

```bash
# Job status
kubectl get jobs -n cal-poly-lomad -l app=lomad-valid

# Container logs for a specific experiment
kubectl logs -n cal-poly-lomad job/lomad-valid-clt-s80

# Detailed events and error info
kubectl describe job -n cal-poly-lomad lomad-valid-clt-s80
```

#### Tips

- Start with small `SIM_S` (e.g., 20) in `submit_sweep.sh` to confirm
  everything works end-to-end, then delete jobs and resubmit at full scale.
- If you update `sim.R`, re-running `submit_sweep.sh` will update the
  ConfigMap automatically — no image rebuild needed.
