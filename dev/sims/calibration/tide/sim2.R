## Calibration simulation — Kubernetes container entrypoint
##
## One container = all S replicates for one value of d.
## Parameters are passed as environment variables by the job spec.
##
## Your task: adapt template.R into this entrypoint. The parameter
## reading below is done for you. You need to:
##   1. Copy the run_rep() function from template.R — it can be used
##      here without modification.
##   2. Implement the simulation loop: draw S seeds, call run_rep() for each,
##      bind the results into a data frame.
##   3. Save the results to SIM_OUT_DIR as a .rds file named after d
##      (e.g. d0-0.rds for d=0, d0-5.rds for d=0.5).
##
## To test locally before submitting to Tide:
##   SIM_D=0 SIM_S=5 SIM_SEED=4853 SIM_OUT_DIR=dev/sims/calibration/results/raw \
##     Rscript dev/sims/calibration/tide/sim.R

library(lomad)
library(dplyr)

# ---- Parameters from environment -------------------------------------------

d       <- as.numeric(Sys.getenv("SIM_D",       "0"))
S       <- as.integer(Sys.getenv("SIM_S",       "200"))
seed0   <- as.integer(Sys.getenv("SIM_SEED",    "4853"))
out_dir <- Sys.getenv("SIM_OUT_DIR", "/jobs/output")

# ---- Fixed parameters (must match template.R) -----------------------------

n     <- 1000
phi   <- 0.5
snr   <- 1
alpha <- 0.05

h_win <- 10
s_win <- 50

# ---- TODO: copy run_rep() from template.R ----------------------------------
#
# Paste run_rep() here. It depends only on the parameters above and on
# the lomad package, so it can be copied verbatim. run_rep() returns a
# data frame with columns: d, seed, clt_rejected, oracle_rejected,
# identity_rejected.

# ---- TODO: simulation loop -------------------------------------------------
#
# Draw S seeds and run run_rep(d, seed) for each.
# Collect the results into a single data frame called `results`.
#
# Hint: set.seed(seed0 + as.integer(d * 100)) before drawing seeds so each
# d value gets a distinct random stream even if seed0 is the same.

# ---- TODO: save results ----------------------------------------------------
#
# Create out_dir if it does not exist.
# Save a list(d = d, S = S, seed0 = seed0, results = results) as an .rds
# file in out_dir. Name the file after d, e.g. "d0-5.rds" for d = 0.5.
# Hint: gsub("\\.", "-", format(d, nsmall = 1)) produces the d part of
# the filename.
