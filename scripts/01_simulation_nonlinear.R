# =============================================================================
# Non-linear simulation study
# =============================================================================
#
# Compares the additive-HSIC discovery procedure against the ICA-based
# LiNGAM baseline on synthetic non-linear structural equation models.
# For each combination of dimension p and sample size n, a random DAG is
# drawn with a target edge density and data are simulated with a mixture
# of linear and non-linear parent effects and non-Gaussian (Laplace /
# uniform) noise. Each method is applied to every replicate and scored by
# directed F1, graph accuracy, structural Hamming distance, adjacency
# mean-squared error, count of mis-oriented edges, and wall-clock time.
#
# Grid:   p in {4, 8, 12, 16},  n in {400, 800, 1200, 1600}
# Trials: 50 independent replicates per (p, n) configuration
#
# Usage:  Rscript scripts/01_simulation_nonlinear.R
# Output: results/sim_nonlinear.rds
#         (long-format data frame; rendered into figures by
#          scripts/05_plot_simulations.R)
#
# Typical runtime: 4-8 h single-threaded on a laptop. The outermost loops
# are embarrassingly parallel and can be trivially split across cores or
# cluster nodes; an optional SLURM wrapper is provided in cluster/.
# =============================================================================

source("R/additive_hsic.R")
source("R/lingam_baseline.R")


# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

P_VALUES        <- c(4, 8, 12, 16)
N_VALUES        <- c(400, 800, 1200, 1600)
N_TRIALS        <- 50
EDGE_DENSITY    <- 0.3
ALPHA           <- 0.05
N_PERM          <- 200
NONLINEAR_FRAC  <- 0.5

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)


# -----------------------------------------------------------------------------
# Main simulation loop
# -----------------------------------------------------------------------------

cat("=== Non-linear simulation study ===\n")
cat(sprintf("p: %s\n", paste(P_VALUES, collapse = ", ")))
cat(sprintf("n: %s\n", paste(N_VALUES, collapse = ", ")))
cat(sprintf("Trials: %d, alpha: %.2f, permutations: %d\n\n",
            N_TRIALS, ALPHA, N_PERM))

results <- list()

for (p in P_VALUES) {
  for (n in N_VALUES) {
    cat(sprintf("--- p=%d, n=%d ---\n", p, n))

    trial_results <- data.frame(
      trial = integer(0), method = character(0),
      F1 = numeric(0), Accuracy = numeric(0),
      SHD = numeric(0), MSE = numeric(0),
      Misoriented = numeric(0), Time = numeric(0),
      stringsAsFactors = FALSE
    )

    for (trial in seq_len(N_TRIALS)) {
      if (trial %% 10 == 0) cat(sprintf("  Trial %d/%d\n", trial, N_TRIALS))

      # Reproducible seed built from the configuration indices.
      set.seed(1000 * p + 100 * n + trial)
      dag  <- random_dag(p, EDGE_DENSITY)
      data <- generate_sem_data(n, dag, nonlinear_frac = NONLINEAR_FRAC)

      # --- Proposed additive-HSIC method ---
      t0 <- proc.time()
      result_prop <- tryCatch(
        discover_dag(data, alpha = ALPHA, n_perm = N_PERM, verbose = FALSE),
        error = function(e) NULL
      )
      time_prop <- (proc.time() - t0)[3]

      if (!is.null(result_prop)) {
        metrics_prop <- evaluate_dag(result_prop$adjacency, dag)
        trial_results <- rbind(trial_results, data.frame(
          trial = trial, method = "Proposed",
          F1 = metrics_prop["F1"],   Accuracy = metrics_prop["Accuracy"],
          SHD = metrics_prop["SHD"], MSE = metrics_prop["MSE"],
          Misoriented = metrics_prop["Misoriented"], Time = time_prop,
          stringsAsFactors = FALSE
        ))
      }

      # --- LiNGAM baseline ---
      t0 <- proc.time()
      result_ling <- tryCatch(
        lingam_discover(data, verbose = FALSE),
        error = function(e) NULL
      )
      time_ling <- (proc.time() - t0)[3]

      if (!is.null(result_ling)) {
        metrics_ling <- evaluate_dag(result_ling$adjacency, dag)
        trial_results <- rbind(trial_results, data.frame(
          trial = trial, method = "LiNGAM",
          F1 = metrics_ling["F1"],   Accuracy = metrics_ling["Accuracy"],
          SHD = metrics_ling["SHD"], MSE = metrics_ling["MSE"],
          Misoriented = metrics_ling["Misoriented"], Time = time_ling,
          stringsAsFactors = FALSE
        ))
      }
    }

    trial_results$p <- p
    trial_results$n <- n
    results[[paste(p, n, sep = "_")]] <- trial_results
  }
}

# Combine and persist.
all_results <- do.call(rbind, results)
rownames(all_results) <- NULL
saveRDS(all_results, "results/sim_nonlinear.rds")
cat("\nResults saved to results/sim_nonlinear.rds\n")
