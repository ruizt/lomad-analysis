## Dockerfile for the lomad calibration simulation
##
## Build from the repo root (so the full package source is in context):
##   docker buildx build --platform linux/amd64 \
##     -f dev/sims/calibration/tide/Dockerfile \
##     -t ghcr.io/<org>/lomad-calib:latest --push .
##
## The container runs dev/sims/calibration/tide/sim.R and reads parameters
## from environment variables (SIM_D, SIM_S, SIM_SEED, SIM_OUT_DIR).

FROM rocker/r-ver:4.5.0

# System libraries required by common R packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

# R package dependencies (install before copying source for layer caching)
RUN Rscript -e "install.packages( \
    c('dplyr', 'remotes'), \
    repos = 'https://packagemanager.posit.co/cran/latest')"

# Copy the lomad package source and install it
COPY . /build/lomad
RUN Rscript -e "remotes::install_local('/build/lomad', dependencies = TRUE, \
    upgrade = 'never')"

# Copy the container entrypoint
COPY dev/sims/calibration/tide/sim.R /sim.R

CMD ["Rscript", "/sim.R"]
