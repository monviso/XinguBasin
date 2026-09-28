#!/bin/bash
#SBATCH --job-name=xingu_land_fig
#SBATCH --partition=core
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/mr3882/temp/Xingu_GAM/logs/landscape_terms_%j.out
#SBATCH --error=/scratch/mr3882/temp/Xingu_GAM/logs/landscape_terms_%j.err

set -euo pipefail

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

mkdir -p /scratch/mr3882/temp/Xingu_GAM/logs
mkdir -p /scratch/mr3882/temp/Xingu_GAM/figures/landscape_temporal

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate RFLD_env

Rscript /scratch/mr3882/temp/Xingu_GAM_package/16_plot_landscape_temporal_terms.R \
  /scratch/mr3882/temp/Xingu_GAM/Xingu_GAM_prepared.rds \
  /scratch/mr3882/temp/Xingu_GAM/figures/landscape_temporal

