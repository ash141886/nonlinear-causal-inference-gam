#!/bin/bash
#SBATCH --job-name=causal_wine
#SBATCH --output=logs/wine_%j.out
#SBATCH --error=logs/wine_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --partition=default

module load R/4.3.0

cd $SLURM_SUBMIT_DIR

mkdir -p logs results figures data

echo "=== Red wine (two-stage) ==="
Rscript scripts/03b_wine_red_twostage.R

echo "=== White wine ==="
Rscript scripts/04_wine_white.R

echo "=== Plotting ==="
Rscript scripts/06_plot_wine.R

echo "=== Done ==="
