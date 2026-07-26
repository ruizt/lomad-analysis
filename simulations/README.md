---

editor_options: 
  markdown: 
    wrap: 72
---

# Simulation Studies

This directory contains all simulation code for the lomad paper. Each study has its own subdirectory with a self-contained design document, R scripts, and Kubernetes submission materials.

## Studies

| Directory | Purpose |
|---------------------------------------|--------------------------------|
| `validation/` | Finite-sample accuracy of CLT approximation, Proposition 1 moments, and end-to-end pipeline coverage |
| `power/` | Power curves and localization as a function of separation *d* for each trend structure |

## Shared container image

All studies share a single Docker image defined by `simulations/Dockerfile`. The image contains R 4.5 and every package needed by any study (lomad, mvtnorm, dplyr) but **no simulation scripts**.

The current public image is:

```         
ghcr.io/ruizt/lomad-sims:latest
```

### Why scripts aren't in the image

Simulation scripts (`sim.R`) are mounted into the container at runtime using Kubernetes [ConfigMaps](https://kubernetes.io/docs/concepts/configuration/configmap/). Each study's `submit_sweep.sh` creates a ConfigMap from its local `sim.R` and mounts it at `/scripts/sim.R` inside the container.

This means:

-   **Script changes don't require an image rebuild.** Just re-run `submit_sweep.sh` and the updated ConfigMap is applied.
-   **The image only needs rebuilding** when a new `lomad` version (in the `lomad-package` repo) or a new runtime dependency is needed. Rebuild manually via the `Build and push lomad-sims image` workflow (`workflow_dispatch`), optionally passing a `LOMAD_REF`.
-   **The same image works for all studies** — the MVN example, the calibration study, and any future studies.

### Rebuilding the image

Build and push from the repo root:

``` bash
# From the repo root:
docker buildx build --platform linux/amd64 \
  -f simulations/Dockerfile \
  -t ghcr.io/ruizt/lomad-sims:latest --push .
```

The image installs the `lomad` package from the `lomad-package` GitHub repo
(see the `LOMAD_REF` build arg in the Dockerfile to pin a version); it does not
copy this repository's contents. The `.dockerignore` at the repo root simply
keeps the build context small.

## Per-study conventions

-   **`design.md`** — prose description of the simulation design, parameter grid, estimands, and intended outputs.
-   **`template.R`** — defines a `run_rep()` function and runs a small local proof-of-concept. This is the source of truth for the simulation logic.
-   **`settings.R`** (if present) — visual walkthrough of single replicates.
-   **`results/`** — gitignored. Populated by local runs and Tide results.
-   **`tide/`** — Kubernetes submission materials (see below).

## Running on Tide (Kubernetes)

Each study's `tide/` directory contains:

| File | Purpose |
|-----------------------------|-------------------------------------------|
| `sim.R` | Simulation script (mounted into the container via ConfigMap) |
| `pvc.yaml` | Persistent Volume Claim for results (create once per study) |
| `submit_sweep.sh` | Creates ConfigMap from `sim.R`, submits one Job per parameter value |
| `job.yaml` | Reference Job template (not submitted directly) |
| `accessor.yaml` | Lightweight pod for copying files off the PVC |
| `fetch.sh` | Copies `.rds` results from PVC to local `results/raw/` |
| `collect.R` | Assembles per-job files into summary results and plots |

### Workflow

``` bash
# 1. Create the PVC (once per study)
kubectl apply -n cal-poly-ruiz -f simulations/<study>/tide/pvc.yaml

# 2. Submit jobs (creates ConfigMap + parallel Jobs)
bash simulations/<study>/tide/submit_sweep.sh

# 3. Monitor
kubectl get jobs -n cal-poly-ruiz -l app=<study-label>

# 4. Fetch results
bash simulations/<study>/tide/fetch.sh

# 5. Assemble
Rscript simulations/<study>/tide/collect.R

# 6. Clean up
kubectl delete jobs -n cal-poly-ruiz -l app=<study-label>
kubectl delete -n cal-poly-ruiz -f simulations/<study>/tide/pvc.yaml
```

### Local testing with Docker

You can test a simulation script locally without the cluster:

``` bash
docker run --rm \
  -e SIM_N=30 -e SIM_S=5 \
  -v $(pwd)/simulations/power/tide/sim.R:/scripts/sim.R \
  -v $(pwd)/simulations/power/results:/jobs/output \
  ghcr.io/ruizt/lomad-sims:latest
```

