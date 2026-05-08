## Shared simulation image for lomad
##
## Contains R + all packages needed by any simulation study (lomad, mvtnorm,
## dplyr). Scripts are NOT baked in — mount them at /scripts/sim.R
## via a Kubernetes ConfigMap or a Docker volume.
##
## Build and push from the repo root:
##   docker buildx build --platform linux/amd64 \
##     -f dev/sims/Dockerfile \
##     -t ghcr.io/ruizt/lomad-sims:latest --push .
##
## Run locally (mount a script + output directory):
##   docker run --rm \
##     -e SIM_N=30 -e SIM_S=5 \
##     -v $(pwd)/dev/sims/tide-example/tide/sim.R:/scripts/sim.R \
##     -v $(pwd)/dev/sims/tide-example/results:/jobs/output \
##     ghcr.io/ruizt/lomad-sims:latest
##
## Interactive R session:
##   docker run --rm -it ghcr.io/ruizt/lomad-sims:latest R

FROM rocker/r-ver:4.5.0

# System libraries required by common R packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

# Install remotes and non-lomad runtime dependencies used by sim scripts
RUN Rscript -e "install.packages( \
    c('remotes', 'dplyr'), \
    repos = 'https://packagemanager.posit.co/cran/latest')"

# Copy the lomad package source and install with Imports
COPY . /build/lomad
RUN Rscript -e "remotes::install_local('/build/lomad', dependencies = NA, \
    upgrade = 'never')"

# Default: run whatever script is mounted at /scripts/sim.R
CMD ["Rscript", "/scripts/sim.R"]
