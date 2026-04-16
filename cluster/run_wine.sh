#!/bin/bash
# -----------------------------------------------------------------------------
# SLURM submission script for the Wine Quality case study.
#
# Runs the red-wine (two-stage) and white-wine analyses sequentially,
# then regenerates the wine-related figures. Adjust --partition,
# --time, and --mem to match the target cluster.
# -----------------------------------------------------------------------------

#SBATCH --job-name=causal_wine
#SBATCH --output=logs/wine_%j.out
#SBATCH --error=logs/wine_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --partition=default

# Load an R installation (edit to match your cluster's module layout).
module load R/4.3.0

cd $SLURM_SUBMIT_DIR

mkdir -p logs results figures data

echo "=== Red wine (two-stage selector) ==="
Rscript scripts/03b_wine_red_twostage.R

echo "=== White wine ==="
Rscript scripts/04_wine_white.R

echo "=== Wine figures ==="
Rscript scripts/06_plot_wine.R

echo "=== Done ==="
