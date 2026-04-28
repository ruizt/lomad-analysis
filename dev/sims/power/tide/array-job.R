## Tide HPC — power study array job
##
## Called by SLURM as:
##   Rscript array-job.R <SLURM_ARRAY_TASK_ID> <out_dir>

args     <- commandArgs(trailingOnly = TRUE)
task_id  <- as.integer(args[1])
out_dir  <- args[2]

pkg_path <- Sys.getenv("LOMAD_PKG", unset = "~/lomad")
if (!requireNamespace("lomad", quietly = TRUE)) {
  devtools::load_all(pkg_path, quiet = TRUE)
} else {
  library(lomad)
}

source(file.path(dirname(sys.frame(1)$ofile), "..", "run.R"))

S_per_task <- 25L
set.seed(task_id)
master_seeds <- sample.int(1e7, nrow(grid) * S_per_task)

grid_idx <- ((task_id - 1L) %/% S_per_task) + 1L
rep_idx  <- ((task_id - 1L) %% S_per_task) + 1L

if (grid_idx > nrow(grid)) {
  message("task_id ", task_id, " out of range — nothing to do.")
  quit(status = 0)
}

params      <- as.list(grid[grid_idx, ])
params$seed <- master_seeds[(grid_idx - 1L) * S_per_task + rep_idx]

result <- run_rep(params)

out_file <- file.path(out_dir, sprintf("rep_%05d.rds", task_id))
saveRDS(result, out_file)
message("Saved: ", out_file)
