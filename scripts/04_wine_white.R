# =============================================================================
# Wine Quality Analysis: White Wine
# =============================================================================
#
# Applies the same additive-HSIC two-stage procedure to the white wine subset
# of the UCI Wine Quality dataset (n=4,898, p=12).
#
# Usage: Rscript scripts/04_wine_white.R
# Output: results/wine_white_results.rds
# =============================================================================

source("R/additive_hsic.R")
source("R/lingam_baseline.R")

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# --- Load white wine data ----------------------------------------------------

cat("=== White Wine Analysis ===\n\n")
wine_white <- load_wine_data(color = "white", data_dir = "data")
var_names <- colnames(wine_white)

# Reference DAG (same chemical relationships apply)
ref_dag <- get_wine_reference_dag(var_names)

# --- Proposed method: two-stage ----------------------------------------------

cat("Running additive-HSIC (two-stage)...\n")
t0 <- proc.time()
result_prop <- discover_dag(
  wine_white,
  alpha = 0.05,
  n_perm = 1000,
  two_stage = TRUE,
  screening_alpha = 0.20,
  target_density = 0.15,
  n_sub = 200,         # subsample for permutation (full sample is 4898)
  seed = 42,
  verbose = TRUE
)
time_prop <- (proc.time() - t0)[3]

metrics_prop <- evaluate_dag(result_prop$adjacency, ref_dag)
cat("\nProposed method results:\n")
cat(sprintf("  F1=%.3f, Accuracy=%.3f, SHD=%d, MSE=%.3f, Time=%.1fs\n",
            metrics_prop["F1"], metrics_prop["Accuracy"],
            metrics_prop["SHD"], metrics_prop["MSE"], time_prop))
cat(sprintf("  Edges: %d (TP=%d, FP=%d, FN=%d)\n",
            sum(result_prop$adjacency),
            metrics_prop["TP"], metrics_prop["FP"], metrics_prop["FN"]))

# --- LiNGAM baseline ---------------------------------------------------------

cat("\nRunning LiNGAM...\n")
t0 <- proc.time()
result_ling <- lingam_discover(wine_white, threshold = 0.1)
time_ling <- (proc.time() - t0)[3]

metrics_ling <- evaluate_dag(result_ling$adjacency, ref_dag)
cat(sprintf("  F1=%.3f, Accuracy=%.3f, SHD=%d, MSE=%.3f, Time=%.1fs\n",
            metrics_ling["F1"], metrics_ling["Accuracy"],
            metrics_ling["SHD"], metrics_ling["MSE"], time_ling))

# --- Save results ------------------------------------------------------------

results <- list(
  proposed = list(result = result_prop, metrics = metrics_prop, time = time_prop),
  lingam = list(result = result_ling, metrics = metrics_ling, time = time_ling),
  data_info = list(color = "white", n = nrow(wine_white), p = ncol(wine_white))
)
saveRDS(results, "results/wine_white_results.rds")
cat("\nResults saved to results/wine_white_results.rds\n")

# --- Edge table (for inspection) ---------------------------------------------

cat("\nSelected edges (proposed method):\n")
edges <- result_prop$edges
edges <- edges[order(-edges$hsic_score), ]
print(edges, row.names = FALSE)
