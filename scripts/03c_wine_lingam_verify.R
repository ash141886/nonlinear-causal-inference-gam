# =============================================================================
# LiNGAM baseline: selection-rule comparison on the red-wine data
# =============================================================================
#
# Compares two edge-selection rules for the ICA-based LiNGAM estimator
# on the UCI Wine Quality (red) dataset:
#
#   Method 1.  Hard threshold on |B_{ij}| at 0.1 (the value used in the
#              LiNGAM literature as a simple, deterministic selector).
#   Method 2.  Element-wise bootstrap p-values followed by
#              Benjamini-Hochberg FDR control; p-values are obtained
#              from the sign-stability of B_{ij} across bootstrap
#              replicates of the ICA decomposition.
#
# The script also prints a consistency check on the adjacency MSE for
# binary matrices, exploiting the identity MSE = SHD / p^2.
#
# Usage:  Rscript scripts/03c_wine_lingam_verify.R
# =============================================================================

suppressPackageStartupMessages({
  library(fastICA)
  library(MASS)
})

# -----------------------------------------------------------------------------
# Data and reference graph
# -----------------------------------------------------------------------------

if (!file.exists("data/winequality-red.csv")) {
  dir.create("data", showWarnings = FALSE, recursive = TRUE)
  download.file("https://archive.ics.uci.edu/ml/machine-learning-databases/wine-quality/winequality-red.csv",
                "data/winequality-red.csv")
}
wine_data <- as.data.frame(scale(read.csv("data/winequality-red.csv", sep = ";")))
var_names <- colnames(wine_data)
n_vars    <- ncol(wine_data)

ref_dag <- matrix(0, n_vars, n_vars, dimnames = list(var_names, var_names))
ref_dag["fixed.acidity",    "pH"] <- 1
ref_dag["volatile.acidity", "pH"] <- 1
ref_dag["citric.acid",      "pH"] <- 1
ref_dag["citric.acid",      "fixed.acidity"] <- 1
ref_dag["residual.sugar",   "density"] <- 1
ref_dag["alcohol",          "density"] <- 1
ref_dag["free.sulfur.dioxide", "total.sulfur.dioxide"] <- 1
ref_dag["alcohol",          "quality"] <- 1
ref_dag["volatile.acidity", "quality"] <- 1
ref_dag["sulphates",        "quality"] <- 1

cat("Reference edges:", sum(ref_dag), "\n\n")

#' Evaluation metrics, including a redundant \code{MSE_check = SHD / p^2}
#' used to verify the algebraic identity for binary matrices.
evaluate_dag <- function(est, ref) {
  n <- nrow(ref)
  tp <- sum(est == 1 & ref == 1)
  fp <- sum(est == 1 & ref == 0)
  fn <- sum(est == 0 & ref == 1)
  tn <- sum(est == 0 & ref == 0) - n  # subtract diagonal matches

  prec <- if ((tp + fp) > 0) tp / (tp + fp) else 0
  rec  <- if ((tp + fn) > 0) tp / (tp + fn) else 0
  f1   <- if ((prec + rec) > 0) 2 * prec * rec / (prec + rec) else 0
  shd  <- sum(abs(est - ref))
  mse  <- mean((est - ref)^2)

  c(F1 = f1, Accuracy = (tp + tn) / (n * (n - 1)), SHD = shd, MSE = mse,
    TP = tp, FP = fp, FN = fn, Precision = prec, Recall = rec,
    Edges = sum(est), MSE_check = shd / (n^2))
}

# =============================================================================
# Method 1: hard threshold on |B_{ij}|
# =============================================================================

cat(paste(rep("=", 60), collapse = ""), "\n")
cat("Method 1: LiNGAM with hard threshold |B| < 0.1\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

set.seed(42)
t0 <- Sys.time()

ica_result <- fastICA(as.matrix(wine_data), n.comp = n_vars,
                      alg.typ = "parallel", fun = "logcosh",
                      method = "C", verbose = FALSE)
# Combine the pre-whitening factor K with the rotation W so that W
# corresponds to the full demixing of the standardised data.
W      <- ica_result$K %*% ica_result$W
B      <- diag(n_vars) - W           # LiNGAM parameterisation: X = B X + e
diag(B) <- 0
B_orig <- B                          # retained for the bootstrap variant

cat(sprintf("Non-zero entries in B (pre-threshold): %d / %d\n",
            sum(abs(B) > 0), n_vars * (n_vars - 1)))
cat(sprintf("Entries with |B| >= 0.1:  %d\n", sum(abs(B) >= 0.1)))
cat(sprintf("Entries with |B| >= 0.05: %d\n", sum(abs(B) >= 0.05)))

B[abs(B) < 0.1] <- 0
dag1 <- (abs(B) > 0) * 1
dimnames(dag1) <- list(var_names, var_names)

time1 <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

m1 <- evaluate_dag(dag1, ref_dag)
cat(sprintf("\nResults:\n"))
cat(sprintf("  Edges: %d\n", m1["Edges"]))
cat(sprintf("  TP = %d  FP = %d  FN = %d\n", m1["TP"], m1["FP"], m1["FN"]))
cat(sprintf("  F1 = %.3f  Accuracy = %.3f\n", m1["F1"], m1["Accuracy"]))
cat(sprintf("  SHD = %d  MSE = %.3f  (cross-check SHD / p^2 = %.3f)\n",
            m1["SHD"], m1["MSE"], m1["MSE_check"]))
cat(sprintf("  Wall-clock time: %.1fs\n", time1))

cat("\nEstimated edges:\n")
for (i in 1:n_vars) for (j in 1:n_vars) {
  if (dag1[i, j] == 1) {
    in_ref <- ref_dag[i, j] == 1
    cat(sprintf("  %s -> %s %s\n",
                var_names[i], var_names[j],
                if (in_ref) "(TP)" else ""))
  }
}

# =============================================================================
# Method 2: bootstrap element-wise BH-FDR
# =============================================================================

cat(paste(rep("\n=", 60), collapse = ""), "\n")
cat("Method 2: LiNGAM with bootstrap BH-FDR selection\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

N_BOOT <- 500
cat(sprintf("Running %d bootstrap replicates...\n", N_BOOT))
t0 <- Sys.time()

n_samples <- nrow(wine_data)
B_boot    <- array(0, dim = c(N_BOOT, n_vars, n_vars))

set.seed(123)
for (b in 1:N_BOOT) {
  boot_idx  <- sample.int(n_samples, replace = TRUE)
  boot_data <- wine_data[boot_idx, ]

  ica_b <- tryCatch(
    fastICA(as.matrix(boot_data), n.comp = n_vars,
            alg.typ = "parallel", fun = "logcosh",
            method = "C", verbose = FALSE),
    error = function(e) NULL
  )

  if (!is.null(ica_b)) {
    W_b <- ica_b$W
    A_b <- tryCatch(solve(W_b), error = function(e) MASS::ginv(W_b))
    diag(A_b) <- 0
    B_boot[b, , ] <- A_b
  }

  if (b %% 100 == 0) cat(sprintf("  %d/%d done\n", b, N_BOOT))
}

boot_time <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("  Bootstrap done in %.0fs.\n\n", boot_time))

# Element-wise two-sided sign-stability p-value: for each off-diagonal
# entry B[i, j], estimate the probability that its bootstrap
# distribution straddles zero. Large values indicate unstable sign and
# are treated as evidence against the corresponding edge.
pvals_boot <- matrix(1, n_vars, n_vars, dimnames = list(var_names, var_names))
for (i in 1:n_vars) {
  for (j in 1:n_vars) {
    if (i == j) next
    boot_vals <- B_boot[, i, j]
    boot_vals <- boot_vals[boot_vals != 0]   # drop failed replicates
    if (length(boot_vals) < 10) next
    orig_val  <- B_orig[i, j]
    if (abs(orig_val) < 1e-10) { pvals_boot[i, j] <- 1; next }
    pvals_boot[i, j] <- mean(boot_vals * sign(orig_val) <= 0)
  }
}

# Apply Benjamini-Hochberg FDR control on the off-diagonal p-values.
off_diag <- which(row(pvals_boot) != col(pvals_boot))
pv_vec   <- pvals_boot[off_diag]
qv_vec   <- p.adjust(pv_vec, method = "BH")
qv_mat   <- matrix(1, n_vars, n_vars, dimnames = list(var_names, var_names))
qv_mat[off_diag] <- qv_vec

for (alpha in c(0.01, 0.05, 0.10, 0.20)) {
  dag2 <- (qv_mat < alpha & abs(B_orig) > 0) * 1
  diag(dag2) <- 0
  dimnames(dag2) <- list(var_names, var_names)

  m2 <- evaluate_dag(dag2, ref_dag)
  cat(sprintf("alpha = %.2f: Edges = %d  F1 = %.3f  Acc = %.3f  SHD = %d  MSE = %.3f  TP = %d  FP = %d  FN = %d\n",
              alpha, m2["Edges"], m2["F1"], m2["Accuracy"],
              m2["SHD"], m2["MSE"], m2["TP"], m2["FP"], m2["FN"]))
}

# =============================================================================
# Algebraic consistency check: MSE = SHD / p^2 for binary matrices
# =============================================================================

cat(paste(rep("\n=", 60), collapse = ""), "\n")
cat("Binary-matrix identity: MSE = SHD / p^2\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat(sprintf("For p = %d, the identity gives MSE = SHD / %d.\n",
            n_vars, n_vars^2))
cat(sprintf("Method 1 SHD = %d -> MSE check = %.4f  (observed MSE = %.4f)\n",
            m1["SHD"], m1["SHD"] / n_vars^2, m1["MSE"]))

# =============================================================================
# LaTeX row for the aggregate-performance table (Method 1)
# =============================================================================

cat(paste(rep("\n=", 60), collapse = ""), "\n")
cat("LaTeX: aggregate-performance row (Method 1, hard threshold)\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat(sprintf("LiNGAM & %.3f & %.3f & %d & %.3f & %.1f \\\\\n",
            m1["F1"], m1["Accuracy"], m1["SHD"], m1["MSE"], time1))

cat("\nDone.\n")
