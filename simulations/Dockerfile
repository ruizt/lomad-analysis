## Shared simulation image for lomad
##
## Contains R + all packages needed by any simulation study (lomad, mvtnorm,
## dplyr). Scripts are NOT baked in — mount them at /scripts/sim.R
## via a Kubernetes ConfigMap or a Docker volume.
##
## Build and push from the repo root:
##   docker buildx build --platform linux/amd64 \
##     -f simulations/Dockerfile \
##     -t ghcr.io/ruizt/lomad-sims:latest --push .
##
## Run locally (mount a script + output directory):
##   docker run --rm \
##     -e SIM_N=30 -e SIM_S=5 \
##     -v $(pwd)/simulations/power/tide/sim.R:/scripts/sim.R \
##     -v $(pwd)/simulations/power/results:/jobs/output \
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

# Install pak, the fast dependency resolver used to install everything below
RUN Rscript -e "install.packages('pak', \
    repos = 'https://packagemanager.posit.co/cran/latest')"

# Install the lomad package (with its Imports) from GitHub, plus non-lomad
# runtime dependencies used by sim scripts. pak resolves and installs all
# dependencies automatically. The package source lives in the companion
# lomad-package repo. Override the ref at build time with
# --build-arg LOMAD_REF=<branch|tag|sha> for a reproducible pin.
ARG LOMAD_REF=main
RUN Rscript -e "pak::pkg_install(c('dplyr', 'ruizt/lomad-package@${LOMAD_REF}'))"

# Default: run whatever script is mounted at /scripts/sim.R
CMD ["Rscript", "/scripts/sim.R"]
