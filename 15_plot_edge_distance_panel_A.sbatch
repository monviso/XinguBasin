#!/bin/bash
#SBATCH --job-name=xingu_edge_fig
#SBATCH --partition=core
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=00:30:00
#SBATCH --output=/scratch/mr3882/temp/Xingu_GAM/logs/edge_panel_A_%j.out
#SBATCH --error=/scratch/mr3882/temp/Xingu_GAM/logs/edge_panel_A_%j.err

set -euo pipefail

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

mkdir -p /scratch/mr3882/temp/Xingu_GAM/logs
mkdir -p /scratch/mr3882/temp/Xingu_GAM/figures

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate RFLD_env

Rscript /scratch/mr3882/temp/Xingu_GAM_package/15_plot_edge_distance_panel_A.R \
  /scratch/mr3882/temp/Xingu_GAM/Xingu_GAM_prepared.rds \
  /scratch/mr3882/temp/Xingu_GAM/figures

