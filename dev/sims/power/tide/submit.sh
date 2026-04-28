#!/bin/bash
#SBATCH --job-name=lomad-power
#SBATCH --array=1-TODO            # update: nrow(grid) * S_per_task once grid is final
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:45:00
#SBATCH --output=logs/power_%A_%a.out
#SBATCH --error=logs/power_%A_%a.err

module load R/4.4.0

OUT_DIR="${SLURM_SUBMIT_DIR}/../results"
mkdir -p "${OUT_DIR}"
mkdir -p "${SLURM_SUBMIT_DIR}/logs"

Rscript array-job.R "${SLURM_ARRAY_TASK_ID}" "${OUT_DIR}"
