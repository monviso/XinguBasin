#!/bin/bash
#SBATCH --job-name=xingu_gam
#SBATCH --partition=core
# Seven architecture classes x (nine weighted models + one unweighted model).
#SBATCH --array=0-69%7
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/mr3882/temp/Xingu_GAM_v3_random/logs/gam_%A_%a.out
#SBATCH --error=/scratch/mr3882/temp/Xingu_GAM_v3_random/logs/gam_%A_%a.err

set -euo pipefail
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
mkdir -p /scratch/mr3882/temp/Xingu_GAM_v3_random/logs

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate RFLD_env

PREPARED=${PREPARED:-/scratch/mr3882/temp/Xingu_GAM_v3_random/Xingu_GAM_prepared.rds}
FITS_DIR=${FITS_DIR:-/scratch/mr3882/temp/Xingu_GAM_v3_random/fits_reduced_basis}
FIT_SCRIPT=${FIT_SCRIPT:-/scratch/mr3882/temp/Xingu_GAM_package/13_fit_Xingu_GAM_task.R}

mkdir -p "$(dirname "${PREPARED}")/logs" "${FITS_DIR}"

Rscript "${FIT_SCRIPT}" \
  "${PREPARED}" \
  "${FITS_DIR}" \
  "$SLURM_ARRAY_TASK_ID"
