#!/bin/bash
# -----------------------------------------------------------------------------
# SLURM submission script for the simulation studies.
#
# Runs the non-linear and linear simulation pipelines sequentially on a
# single compute node. Adjust --partition, --time, and --mem to match
# the target cluster; the R module name may also need to be edited.
# -----------------------------------------------------------------------------

#SBATCH --job-name=causal_sim
#SBATCH --output=logs/sim_%j.out
#SBATCH --error=logs/sim_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --partition=default

# Load an R installation (edit to match your cluster's module layout).
module load R/4.3.0

cd $SLURM_SUBMIT_DIR

mkdir -p logs results figures data

echo "=== Non-linear simulation study ==="
Rscript scripts/01_simulation_nonlinear.R

echo "=== Linear-data residual diagnostics ==="
Rscript scripts/02_simulation_linear.R

echo "=== Done ==="
