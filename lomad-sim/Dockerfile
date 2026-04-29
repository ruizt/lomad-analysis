# Dockerfile — lomad simulation runner
# Build from the repo root so the package source is available.
# Replace <github-username> with your GitHub username (e.g. otishunt):
#   docker buildx build --platform linux/amd64 -f lomad-sim/Dockerfile \
#     -t ghcr.io/<github-username>/lomad-sim:latest --push .

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

# Default parameters (SIM_H=0 and SIM_S=0 mean auto-select)
ENV SIM_N=500
ENV SIM_D=5
ENV SIM_SEED=32026
ENV SIM_H=0
ENV SIM_S=0
ENV SIM_MAX_PQ=3
ENV SIM_OUT_DIR=/jobs/output

CMD ["Rscript", "/jobs/sim.R"]
