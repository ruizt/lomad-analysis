## Shared simulation image for lomad
##
## Contains R + all packages needed by any simulation study (lomad, mvtnorm,
## dplyr). Scripts are NOT baked in — mount them at /scripts/sim.R
## via a Kubernetes ConfigMap or a Docker volume.
##
## Build and push with simulations/build-image.sh, which builds the lomad
## tarball from the sibling lomad-package checkout first. Building this
## Dockerfile directly will fail unless that tarball is already in place.
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

# Install lomad from a locally built source tarball, plus the non-lomad
# runtime dependencies used by the sim scripts. pak resolves dependencies for
# both. build-image.sh builds the tarball from the sibling lomad-package
# checkout, unpacks it, and drops the source tree at the path below.
#
# Local rather than from GitHub because lomad-package is private: pak inside
# the container has no credentials. Once the repo is public this can go back
# to a one-line remote install:
#   ARG LOMAD_REF=main
#   RUN Rscript -e "pak::pkg_install(c('dplyr', 'ruizt/lomad-package@${LOMAD_REF}'))"
# pak's local:: ref wants a package directory, not an archive, so
# build-image.sh unpacks the tarball it builds.
COPY simulations/_build/lomad /tmp/lomad
RUN Rscript -e "pak::pkg_install(c('dplyr', 'local::/tmp/lomad'))" \
    && rm -rf /tmp/lomad

# Default: run whatever script is mounted at /scripts/sim.R
CMD ["Rscript", "/scripts/sim.R"]
