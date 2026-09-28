#!/bin/bash
#SBATCH --job-name=xingu_effects
#SBATCH --partition=core
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/mr3882/temp/Xingu_GAM_v3_random/logs/LSE_effects_%j.out
#SBATCH --error=/scratch/mr3882/temp/Xingu_GAM_v3_random/logs/LSE_effects_%j.err

set -euo pipefail
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate RFLD_env
export R_LIBS_USER=${R_LIBS_USER:-$HOME/R/x86_64-pc-linux-gnu-library}

PREPARED=${PREPARED:-/scratch/mr3882/temp/Xingu_GAM_v3_random/Xingu_GAM_prepared.rds}
FITS_ROOT=${FITS_ROOT:-/scratch/mr3882/temp/Xingu_GAM_v3_random/fits_reduced_basis}
OUTPUT_ROOT=${OUTPUT_ROOT:-/scratch/mr3882/temp/Xingu_GAM_v3_random/figures/LSE_all_effects}
SCRIPT_ROOT=${SCRIPT_ROOT:-/scratch/mr3882/temp/Xingu_GAM_package}

mkdir -p /scratch/mr3882/temp/Xingu_GAM_v3_random/logs
Rscript "${SCRIPT_ROOT}/22_plot_all_LSE_effects.R" \
  "${PREPARED}" "${FITS_ROOT}" "${OUTPUT_ROOT}"
