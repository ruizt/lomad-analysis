# MVN Coverage Example — Tide HPC Walkthrough

This is a minimal working example of running a simulation study on [Tide](https://tide.calpoly.edu), Cal Poly's Kubernetes-based HPC cluster. The statistical content is intentionally simple so the focus stays on the infrastructure.

## What the simulation does

We draw `S` independent replicates of `n` i.i.d. samples from a 3-dimensional multivariate normal distribution with known mean `mu` and covariance `Sigma`. For each replicate we build a Hotelling T² confidence region for the mean at the 95% level and check whether `mu` is inside. The expected coverage is 95% for every sample size.

We sweep over sample sizes `n ∈ {10, 30, 100, 500}`. Each value of `n` is a separate Kubernetes **Job** that runs in parallel.

### Parameters

| Symbol   | Value                                           |
|----------|-------------------------------------------------|
| `p`      | 3 (dimension)                                   |
| `mu`     | (1, −0.5, 2)                                    |
| `Sigma`  | 3×3 positive-definite matrix (see `template.R`) |
| `alpha`  | 0.05 (nominal miscoverage)                      |
| `n_vals` | 10, 30, 100, 500                                |
| `S`      | 200 replicates per `n`                          |

### Expected result

Coverage ≈ 95% (± sampling error) for all sample sizes. The Hotelling T² region is exact for normal data, so any deviation from 95% is pure Monte Carlo noise.

------------------------------------------------------------------------

## File layout

```         
tide-example/
├── design.md                ← you are here
├── building-containers.md   ← how to build/push the Docker image (optional)
├── template.R               ← local proof-of-concept (run in RStudio)
├── results/                 ← output (gitignored)
│   ├── results.rds          ← all replicates (local or assembled from Tide)
│   ├── results_summary.rds  ← coverage rates by n
│   ├── coverage_plot.png
│   └── raw/                 ← per-job .rds files fetched from Tide
└── tide/                    ← Kubernetes scaffolding
    ├── sim.R                ← simulation script (mounted into container)
    ├── submit.sh            ← runs the full pipeline (PVC → jobs → wait → fetch)
    ├── submit_sweep.sh      ← submits one Job per n (called by submit.sh)
    ├── pvc.yaml             ← shared storage (create once)
    ├── job.yaml             ← Job template (reference only)
    ├── accessor.yaml        ← lightweight pod for file retrieval
    ├── fetch.sh             ← copies results from PVC to local machine
    └── collect.R            ← assembles per-job files into summary
```

------------------------------------------------------------------------

## How the pieces fit together

There are two phases to running a simulation on Tide:

1.  **Run** — submit Kubernetes Jobs that pull the container image and execute your simulation script. Each job writes its results to a shared persistent volume (PVC).
2.  **Collect** — copy the results from the PVC to your local machine and assemble them.

### The shared container image

All simulation studies use a common Docker image that is pre-built and publicly available at:

```         
ghcr.io/ruizt/lomad-sims:latest
```

The image contains R 4.5 and the packages needed to run simulations (lomad, mvtnorm, dplyr) but **no simulation scripts**. Instead, each study provides its own `sim.R` script, which gets mounted into the container at `/scripts/sim.R` at runtime using a Kubernetes [ConfigMap](https://kubernetes.io/docs/concepts/configuration/configmap/).

This means you can update a simulation script without rebuilding the image. The image only needs to be rebuilt when R package dependencies change. See `building-containers.md` for instructions on building and pushing the image.

### Key concepts

-   **Docker image**: a self-contained snapshot of your R environment (packages, system libraries). Shared across studies; defined by `sims/Dockerfile`.
-   **ConfigMap**: a Kubernetes object that holds small files (like `sim.R`). Created automatically by `submit_sweep.sh` from your local `sim.R`.
-   **Kubernetes Job**: a one-shot task that runs a container to completion. If it fails, Kubernetes can retry it (controlled by `backoffLimit`).
-   **PVC (Persistent Volume Claim)**: shared disk storage that persists across jobs. All jobs write their output here; we read from it when collecting.
-   **Environment variables**: how we pass parameters (like `SIM_N`) into the container without rebuilding the image.

------------------------------------------------------------------------

## Step-by-step walkthrough

### Prerequisites

Before starting, make sure you have:

-   [ ] `kubectl` installed and configured to talk to the Tide cluster
-   [ ] Access to the `cal-poly-ruiz` namespace

The Nautilus/NRP cluster documentation lives at <https://nrp.ai/documentation/>. The pages linked below are the most relevant; refer to the full docs for troubleshooting and advanced topics.

#### Installing kubectl + kubelogin

Follow the [Nautilus Getting Started](https://nrp.ai/documentation/userdocs/start/getting-started/) guide. The key steps are:

1.  **Install kubectl** — see the [official instructions](https://kubernetes.io/docs/tasks/tools/install-kubectl/) for your OS, or on macOS: `brew install kubectl`.

2.  **Install the kubelogin OIDC plugin** — this is **required** for Nautilus authentication. Follow the [kubelogin setup instructions](https://github.com/int128/kubelogin?tab=readme-ov-file#setup). On macOS with Homebrew:

    ``` bash
    brew install int128/kubelogin/kubelogin
    ```

3.  **Download the Nautilus kubeconfig** and save it as `~/.kube/config`:

    ``` bash
    mkdir -p ~/.kube
    curl -o ~/.kube/config -fSL "https://nrp.ai/config"
    ```

4.  **Authenticate** — the first `kubectl` command will open a browser window for CILogon. Select your institution, log in, and the token is cached automatically.

#### Verify your setup

``` bash
kubectl version --client               # kubectl installed?
kubectl get nodes                      # kubelogin working? (opens browser)
kubectl get ns cal-poly-ruiz          # namespace exists?
kubectl auth can-i create jobs \
  -n cal-poly-ruiz                    # permission to submit jobs?
```

------------------------------------------------------------------------

### Step 0: Run locally first

Always verify your simulation works on your own machine before sending it to the cluster. This catches bugs cheaply.

``` r
# In RStudio:
source("sims/tide-example/template.R")
```

Or from the terminal:

``` bash
Rscript sims/tide-example/template.R
```

You should see coverage rates near 95% for every `n`. If something is wrong, fix it here before proceeding.

------------------------------------------------------------------------

### Step 1: Create the PVC

The PVC is shared storage where all jobs write their output. Each simulation study needs its own PVC (with a different name), but you only create it **once** — it persists on the cluster until you explicitly delete it.

``` bash
kubectl apply -n cal-poly-ruiz -f sims/tide-example/tide/pvc.yaml
```

Verify it exists:

``` bash
kubectl get pvc -n cal-poly-ruiz mvn-example-results
```

We use `rook-cephfs-tide` (CephFS) because multiple jobs write to the PVC simultaneously, which requires `ReadWriteMany`. The default storage class (`rook-ceph-block`) is `ReadWriteOnce` and would fail when a second job tries to mount. See the [Nautilus storage docs](https://nrp.ai/documentation/userdocs/storage/ceph/) for all available storage classes.

------------------------------------------------------------------------

### Step 2: Submit the jobs

``` bash
bash sims/tide-example/tide/submit_sweep.sh
```

This does two things:

1.  **Creates a ConfigMap** from the local `tide/sim.R` and uploads it to the cluster. The ConfigMap is mounted into each container at `/scripts/sim.R`.
2.  **Submits four Jobs** (one per `n` value). Each job pulls the shared `lomad-sims` image, mounts the ConfigMap and the PVC, reads its `SIM_N` environment variable, runs `S` replicates, and saves an `.rds` file (e.g. `n30.rds`) to the PVC.

------------------------------------------------------------------------

### Step 3: Monitor progress

``` bash
# Check job status (look for "Complete" in the COMPLETIONS column)
kubectl get jobs -n cal-poly-ruiz -l app=mvn-example

# Watch a specific job's logs
kubectl logs -n cal-poly-ruiz job/mvn-example-n30

# If something went wrong, describe the job for events and error details
kubectl describe job -n cal-poly-ruiz mvn-example-n30
```

Wait until all four jobs show `1/1` in the COMPLETIONS column.

------------------------------------------------------------------------

### Step 4: Fetch results

``` bash
bash sims/tide-example/tide/fetch.sh
```

This spins up a temporary pod, copies the `.rds` files from the PVC to `sims/tide-example/results/raw/`, and cleans up the pod.

------------------------------------------------------------------------

### Step 5: Assemble and plot

``` bash
Rscript sims/tide-example/tide/collect.R
```

Or in RStudio:

``` r
source("sims/tide-example/tide/collect.R")
```

This reads the per-job `.rds` files, computes coverage rates by `n`, and saves `results/results.rds`, `results/results_summary.rds`, and `results/coverage_plot.png`.

------------------------------------------------------------------------

### Step 6: Clean up

When you're done, delete the jobs and (optionally) the PVC:

``` bash
# Delete all jobs from this example
kubectl delete jobs -n cal-poly-ruiz -l app=mvn-example

# Delete the PVC (removes the stored results from the cluster)
kubectl delete -n cal-poly-ruiz -f sims/tide-example/tide/pvc.yaml
```

------------------------------------------------------------------------

## Troubleshooting

| Symptom | Likely cause | Fix |
|-----------------------|---------------------------------|-----------------|
| Job stuck in `Pending` | Image can't be pulled | Check image name in `submit_sweep.sh`; verify the ghcr.io package is public |
| Job fails immediately | R error in `sim.R` | `kubectl logs -n cal-poly-ruiz job/<name>` to see the error |
| `fetch.sh` times out | Accessor pod can't start | Check PVC name matches in `accessor.yaml` |
| No `.rds` files after fetch | Jobs haven't finished yet | `kubectl get jobs` — wait for COMPLETIONS = 1/1 |
| Coverage far from 95% | Bug in `run_rep()` | Run `template.R` locally to debug |

------------------------------------------------------------------------

## Adapting this for your own simulation

To use this scaffolding for a different simulation:

1.  **Write your `template.R`** with a `run_rep()` function that takes parameters and a seed, and returns a one-row data frame.
2.  **Copy `tide/sim.R`** and replace the `run_rep()` function and parameter block. Decide which parameter to sweep over (the one that becomes an environment variable).
3.  **If you need new R packages**, add them to the shared `Dockerfile` at `sims/Dockerfile` and rebuild the image (see `building-containers.md`).
4.  **Create a `tide/pvc.yaml`** with a unique PVC name for your study.
5.  **Update `tide/submit_sweep.sh`** with your sweep values, ConfigMap name, and PVC name.
6.  **Update `tide/collect.R`** to summarise your specific output columns.
7.  **Update `tide/accessor.yaml`** and `tide/fetch.sh` with your PVC name and pod name.

The general pattern is always the same: local development → shared Docker image → ConfigMap + parallel Jobs → PVC → fetch → collect.

------------------------------------------------------------------------

## Nautilus documentation

These pages from the Nautilus/NRP docs are most relevant to this workflow:

-   [Getting Started](https://nrp.ai/documentation/userdocs/start/getting-started/) — account setup, kubectl, kubelogin
-   [Batch Jobs (tutorial)](https://nrp.ai/documentation/userdocs/tutorial/jobs/) — simple job example
-   [Running Batch Jobs](https://nrp.ai/documentation/userdocs/running/jobs/) — advanced patterns, git pulls, retries
-   [Storage (Ceph FS / RBD)](https://nrp.ai/documentation/userdocs/storage/ceph/) — storage classes, CephFS vs block
-   [Cluster Policies](https://nrp.ai/documentation/userdocs/start/policies/) — resource limits, acceptable use
