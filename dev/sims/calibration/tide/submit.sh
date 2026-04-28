#!/bin/bash
#SBATCH --job-name=lomad-calib
#SBATCH --array=1-9000           # nrow(grid) * S_per_task = 36 * 25 * 10 = 9000
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH --output=logs/calib_%A_%a.out
#SBATCH --error=logs/calib_%A_%a.err

# --- Environment -------------------------------------------------------------
module load R/4.4.0              # adjust to Tide's available R version

OUT_DIR="${SLURM_SUBMIT_DIR}/../results"
mkdir -p "${OUT_DIR}"
mkdir -p "${SLURM_SUBMIT_DIR}/logs"

# --- Run ---------------------------------------------------------------------
Rscript array-job.R "${SLURM_ARRAY_TASK_ID}" "${OUT_DIR}"
