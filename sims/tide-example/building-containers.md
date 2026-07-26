# Building and Pushing Container Images

This document covers how to build, push, and manage the shared Docker image
used by all simulation studies on Tide. You only need this if you're setting
up the image for the first time or adding new R package dependencies.

If the image already exists and is public on ghcr.io, you can skip this
entirely and go straight to `design.md`.

## Prerequisites

- [Docker Desktop](https://www.docker.com/products/docker-desktop/) installed
  and running.
- A [GitHub Personal Access Token](https://github.com/settings/tokens) with
  `write:packages` scope.

## Authenticate with ghcr.io

```bash
echo "YOUR_GITHUB_PAT" | docker login ghcr.io -u USERNAME --password-stdin
```

Replace `USERNAME` with your GitHub username (lowercase). You only need to do
this once per machine — Docker caches the credential.

## Build and push

From the **repo root**:

```bash
docker buildx build --platform linux/amd64 \
  -f sims/Dockerfile \
  -t ghcr.io/USERNAME/lomad-sims:latest --push .
```

- `--platform linux/amd64` — ensures the image runs on the cluster, even when
  building on Apple Silicon.
- `-f sims/Dockerfile` — path to the Dockerfile.
- `-t ghcr.io/USERNAME/lomad-sims:latest` — image name and tag.
- `--push` — uploads to ghcr.io in one step.
- `.` — build context is the repo root. The Dockerfile installs the `lomad`
  package from GitHub (`ruizt/lomad-package`) rather than from this context;
  pass `--build-arg LOMAD_REF=<branch|tag|sha>` to pin a version.

Build time is ~4 minutes from scratch. Subsequent builds are faster because
Docker caches the R package install layers.

## Make the image public

New ghcr.io packages default to **private**. The cluster can't pull a private
image without an image pull secret. To avoid that complexity, make the image
public:

1. Go to `https://github.com/users/USERNAME/packages/container/lomad-sims/settings`
2. Scroll to "Danger Zone" → "Change package visibility" → set to **Public**.

## What's in the image

The shared `Dockerfile` (`sims/Dockerfile`) produces an image with:

- **Base**: `rocker/r-ver:4.5.0` (Debian + R 4.5)
- **System libraries**: libcurl, libssl, libxml2
- **R packages**: pak, dplyr, plus the lomad package (which pulls in its own
  dependencies such as mvtnorm, fda, and roll)
  (Imports only — Suggests like ggplot2 are excluded to keep the image lean)
  (installed from source)
- **No simulation scripts** — scripts are mounted at runtime via Kubernetes
  ConfigMap or Docker volume

The image expects a script at `/scripts/sim.R` and runs it via:

```
CMD ["Rscript", "/scripts/sim.R"]
```

## When to rebuild

Rebuild and push when:

- An R package dependency is added or updated
- A new `lomad` version is released in the `lomad-package` repo (rebuild to pick
  it up, optionally pinning with `--build-arg LOMAD_REF=...`)
- The base R version needs to change

You do **not** need to rebuild when:

- A simulation script (`sim.R`) changes — scripts are mounted at runtime
- Environment variable values change (like `SIM_S`)

## Updating an existing image

If you need to update the image (e.g., after a bug fix in lomad or adding a
new package to the Dockerfile), the process is the same as the initial build:

```bash
docker buildx build --platform linux/amd64 \
  -f sims/Dockerfile \
  -t ghcr.io/USERNAME/lomad-sims:latest --push .
```

Because the image is public and the job specs use `imagePullPolicy: Always`,
any new jobs submitted after the push will automatically pull the updated
image. No changes to the Kubernetes YAML are needed.

**Note for collaborators:** you need write access to the `lomad-sims` package
on ghcr.io (granted through the GitHub repo settings) and a GitHub PAT with
`write:packages` scope. Authenticate with `docker login ghcr.io` using your
own username and PAT before pushing.

### Adding a new R package

1. Edit `sims/Dockerfile` and add the package to the `pak::pkg_install()`
   call:

   ```dockerfile
   RUN Rscript -e "pak::pkg_install(c('dplyr', 'NEW_PACKAGE', 'ruizt/lomad-package@${LOMAD_REF}'))"
   ```

2. Rebuild and push:

   ```bash
   docker buildx build --platform linux/amd64 \
     -f sims/Dockerfile \
     -t ghcr.io/USERNAME/lomad-sims:latest --push .
   ```

3. Verify the package is available:

   ```bash
   docker run --rm ghcr.io/USERNAME/lomad-sims:latest \
     Rscript -e "library(NEW_PACKAGE); cat('OK\n')"
   ```

## Testing locally with Docker

Run a simulation script locally without the cluster:

```bash
docker run --rm \
  -e SIM_N=30 -e SIM_S=5 \
  -v $(pwd)/sims/tide-example/tide/sim.R:/scripts/sim.R \
  -v $(pwd)/sims/tide-example/results:/jobs/output \
  ghcr.io/USERNAME/lomad-sims:latest
```

Or open an interactive R session inside the container:

```bash
docker run --rm -it ghcr.io/USERNAME/lomad-sims:latest R
```

## Versioned tags

The `:latest` tag is mutable — it always points to the most recent push. For
reproducibility (e.g., pinning a paper submission to an exact image), push
with a versioned tag alongside `:latest`:

```bash
docker buildx build --platform linux/amd64 \
  -f sims/Dockerfile \
  -t ghcr.io/USERNAME/lomad-sims:latest \
  -t ghcr.io/USERNAME/lomad-sims:v0.1.0 \
  --push .
```

## .dockerignore

The `.dockerignore` at the repo root controls what gets sent to Docker during
the build. It excludes `.git/`, simulation `results/`, large `_`-prefixed data,
`.Renviron`, and other files. The image installs the R package from GitHub
rather than from the build context, so this mainly keeps the build context
small and prevents accidental inclusion of sensitive files.
