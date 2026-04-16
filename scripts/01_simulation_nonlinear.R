# =============================================================================
# Simulation Study: Non-linear Settings (Section 4.1)
# =============================================================================
#
# Generates Figures 1-6 and simulation results from the manuscript.
# Compares additive-HSIC against LiNGAM across:
#   p in {4, 8, 12, 15}, n in {400, 800, 1200, 1500}
#   50 Monte Carlo trials per configuration
#
# Usage: Rscript scripts/01_simulation_nonlinear.R
# Output: results/sim_nonlinear.rds, figures/f1_r34.pdf, etc.
#
# Estimated runtime: ~4-8 hours (parallelizable)
# =============================================================================

source("R/additive_hsic.R")
source("R/lingam_baseline.R")

# --- Configuration -----------------------------------------------------------

P_VALUES   <- c(4, 8, 12, 15)
N_VALUES   <- c(400, 800, 1200, 1500)
N_TRIALS   <- 50
EDGE_DENSITY <- 0.3
ALPHA      <- 0.05
N_PERM     <- 200
NONLINEAR_FRAC <- 0.5

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# --- Main Simulation Loop ----------------------------------------------------

cat("=== Non-linear Simulation Study ===\n")
cat(sprintf("p: %s\n", paste(P_VALUES, collapse = ", ")))
cat(sprintf("n: %s\n", paste(N_VALUES, collapse = ", ")))
cat(sprintf("Trials: %d, alpha: %.2f, B: %d\n\n", N_TRIALS, ALPHA, N_PERM))

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

      # Generate random DAG and data
      set.seed(1000 * p + 100 * n + trial)
      dag <- random_dag(p, EDGE_DENSITY)
      data <- generate_sem_data(n, dag, nonlinear_frac = NONLINEAR_FRAC)

      # --- Proposed method ---
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
          F1 = metrics_prop["F1"], Accuracy = metrics_prop["Accuracy"],
          SHD = metrics_prop["SHD"], MSE = metrics_prop["MSE"],
          Misoriented = metrics_prop["Misoriented"], Time = time_prop,
          stringsAsFactors = FALSE
        ))
      }

      # --- LiNGAM ---
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
          F1 = metrics_ling["F1"], Accuracy = metrics_ling["Accuracy"],
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

# Combine and save
all_results <- do.call(rbind, results)
rownames(all_results) <- NULL
saveRDS(all_results, "results/sim_nonlinear.rds")
cat("\nResults saved to results/sim_nonlinear.rds\n")
