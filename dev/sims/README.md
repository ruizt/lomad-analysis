# Simulation Studies

This directory contains all simulation code for the lomad paper. Each study
has its own subdirectory with a self-contained design document, R scripts, and
Kubernetes submission materials.

## Studies

| Directory | Purpose |
|-----------|---------|
| `calibration/` | Validate CLT approximation: rejection rates under H₀ across noise/window conditions |
| `power/` | Power curves as a function of separation *d* for each trend structure |
| `tide-example/`  | Minimal working example (MVN coverage) to learn the Tide workflow |

## Shared container image

All studies share a single Docker image defined by `dev/sims/Dockerfile`. The
image contains R 4.5 and every package needed by any study (lomad, mvtnorm,
dplyr) but **no simulation scripts**.

The current public image is:

```
ghcr.io/ruizt/lomad-sims:latest
```

### Why scripts aren't in the image

Simulation scripts (`sim.R`) are mounted into the container at runtime using
Kubernetes [ConfigMaps](https://kubernetes.io/docs/concepts/configuration/configmap/).
Each study's `submit_sweep.sh` creates a ConfigMap from its local `sim.R` and
mounts it at `/scripts/sim.R` inside the container.

This means:

- **Script changes don't require an image rebuild.** Just re-run
  `submit_sweep.sh` and the updated ConfigMap is applied.
- **The image only needs rebuilding** when R package dependencies change or
  the lomad package source is updated.
- **The same image works for all studies** — the MVN example, the calibration
  study, and any future studies.

### Rebuilding the image

See `tide-example/building-containers.md` for full instructions. In short:

```bash
# From the repo root:
docker buildx build --platform linux/amd64 \
  -f dev/sims/Dockerfile \
  -t ghcr.io/USERNAME/lomad-sims:latest --push .
```

The `.dockerignore` at the repo root keeps the build context lean by excluding
`.git/`, `dev/`, and other files not needed for the R package install.

## Per-study conventions

- **`design.md`** — prose description of the simulation design, parameter
  grid, estimands, and intended outputs.
- **`template.R`** — defines a `run_rep()` function and runs a small local
  proof-of-concept. This is the source of truth for the simulation logic.
- **`settings.R`** (if present) — visual walkthrough of single replicates.
- **`results/`** — gitignored. Populated by local runs and Tide results.
- **`tide/`** — Kubernetes submission materials (see below).

## Running on Tide (Kubernetes)

Each study's `tide/` directory contains:

| File | Purpose |
|------|---------|
| `sim.R` | Simulation script (mounted into the container via ConfigMap) |
| `pvc.yaml` | Persistent Volume Claim for results (create once per study) |
| `submit_sweep.sh` | Creates ConfigMap from `sim.R`, submits one Job per parameter value |
| `job.yaml` | Reference Job template (not submitted directly) |
| `accessor.yaml` | Lightweight pod for copying files off the PVC |
| `fetch.sh` | Copies `.rds` results from PVC to local `results/raw/` |
| `collect.R` | Assembles per-job files into summary results and plots |

### Workflow

```bash
# 1. Create the PVC (once per study)
kubectl apply -n cal-poly-lomad -f dev/sims/<study>/tide/pvc.yaml

# 2. Submit jobs (creates ConfigMap + parallel Jobs)
bash dev/sims/<study>/tide/submit_sweep.sh

# 3. Monitor
kubectl get jobs -n cal-poly-lomad -l app=<study-label>

# 4. Fetch results
bash dev/sims/<study>/tide/fetch.sh

# 5. Assemble
Rscript dev/sims/<study>/tide/collect.R

# 6. Clean up
kubectl delete jobs -n cal-poly-lomad -l app=<study-label>
kubectl delete -n cal-poly-lomad -f dev/sims/<study>/tide/pvc.yaml
```

### Local testing with Docker

You can test a simulation script locally without the cluster:

```bash
docker run --rm \
  -e SIM_N=30 -e SIM_S=5 \
  -v $(pwd)/dev/sims/tide-example/tide/sim.R:/scripts/sim.R \
  -v $(pwd)/dev/sims/tide-example/results:/jobs/output \
  ghcr.io/ruizt/lomad-sims:latest
```

See `tide-example/design.md` for a detailed step-by-step walkthrough
including install instructions for kubectl and kubelogin.
