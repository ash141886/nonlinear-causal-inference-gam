#!/bin/bash
#SBATCH --job-name=causal_sim
#SBATCH --output=logs/sim_%j.out
#SBATCH --error=logs/sim_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --partition=default

# Adjust --partition to your cluster's queue name
# (e.g., "compute", "general", "short", "long")

module load R/4.3.0   # adjust to your cluster's R module

cd $SLURM_SUBMIT_DIR

mkdir -p logs results figures data

echo "=== Starting non-linear simulations ==="
Rscript scripts/01_simulation_nonlinear.R

echo "=== Starting linear simulations ==="
Rscript scripts/02_simulation_linear.R

echo "=== Done ==="
