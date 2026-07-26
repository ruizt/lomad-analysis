## Shared simulation image for lomad
##
## Contains R + all packages needed by any simulation study (lomad, mvtnorm,
## dplyr). Scripts are NOT baked in — mount them at /scripts/sim.R
## via a Kubernetes ConfigMap or a Docker volume.
##
## Build and push with simulations/build-image.sh. Building this Dockerfile
## directly will fail on the COPY, which the script sets up.
##
## Run locally (mount a script + output directory):
##   docker run --rm \
##     -e SIM_N=30 -e SIM_S=5 \
##     -v $(pwd)/simulations/power/tide/sim.R:/scripts/sim.R \
##     -v $(pwd)/simulations/power/results:/jobs/output \
##     ghcr.io/ruizt/lomad-simulations:latest
##
## Interactive R session:
##   docker run --rm -it ghcr.io/ruizt/lomad-simulations:latest R

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

# Install lomad, plus the non-lomad runtime dependencies used by the sim
# scripts; pak resolves dependencies for both.
#
#   LOMAD_SOURCE=remote (default) installs from GitHub at LOMAD_REF. Requires
#     ruizt/lomad-package to be public, or pak has no credentials.
#   LOMAD_SOURCE=local installs the source tree that build-image.sh unpacks
#     into simulations/_build/lomad. Use this while the package repo is
#     private, or to build an image that matches your working tree rather than
#     a branch -- which is what you want when re-running simulations against a
#     change that has not been released.
#
# pak's local:: ref resolves a package directory, not an archive, which is why
# build-image.sh unpacks the tarball it builds.
#
# _build is always present because build-image.sh creates it; it is only read
# when LOMAD_SOURCE=local. Building this Dockerfile directly, without the
# script, will fail on the COPY.
ARG LOMAD_SOURCE=remote
ARG LOMAD_REF=main

COPY simulations/_build /tmp/build

RUN if [ "${LOMAD_SOURCE}" = "local" ]; then \
      Rscript -e "pak::pkg_install(c('dplyr', 'local::/tmp/build/lomad'))" ; \
    else \
      Rscript -e "pak::pkg_install(c('dplyr', 'ruizt/lomad-package@${LOMAD_REF}'))" ; \
    fi \
 && rm -rf /tmp/build

# Default: run whatever script is mounted at /scripts/sim.R
CMD ["Rscript", "/scripts/sim.R"]
