# Dockerfile — lomad simulation runner
# Build from the repo root so the package source is available:
#   docker buildx build --platform linux/amd64 -f lomad-sim/Dockerfile -t ghcr.io/otishunt/lomad-sim:latest --push .

FROM rocker/r-ver:4.4.2

# System libraries required by lomad dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
      libcurl4-openssl-dev \
      libssl-dev \
      libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

# Install lomad's R dependencies
RUN Rscript -e " \
  install.packages( \
    c('fda', 'mvtnorm', 'roll'), \
    repos  = 'https://cloud.r-project.org', \
    Ncpus  = parallel::detectCores() \
  )"

# Copy lomad package source and install it
WORKDIR /pkg
COPY . .
RUN Rscript -e "install.packages('.', repos = NULL, type = 'source')"

# Copy simulation script
WORKDIR /jobs
COPY lomad-sim/sim.R /jobs/sim.R

# Default parameters
ENV SIM_N=500
ENV SIM_D=5
ENV SIM_SEED=32026
ENV SIM_H=30
ENV SIM_Q=30
ENV SIM_B=1000
ENV SIM_METHOD=boot
ENV SIM_NCORES=1
ENV SIM_OUT_DIR=/jobs/output

CMD ["Rscript", "/jobs/sim.R"]
