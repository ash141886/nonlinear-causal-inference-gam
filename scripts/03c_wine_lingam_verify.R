# =============================================================================
# VERIFY LiNGAM results on Wine data
# =============================================================================
# Issues to check:
# 1. SHD=56 and MSE=0.778 are inconsistent (MSE should = SHD/144)
# 2. Code uses hard threshold (|B|<0.1), but paper claims bootstrap BH-FDR
# 3. Reviewer CE asked: "same FDR for LiNGAM?"
# =============================================================================

suppressPackageStartupMessages({
  library(fastICA)
  library(MASS)
})

# Wine data (standardized)
if (!file.exists("data/winequality-red.csv")) {
  dir.create("data", showWarnings = FALSE, recursive = TRUE)
  download.file("https://archive.ics.uci.edu/ml/machine-learning-databases/wine-quality/winequality-red.csv",
                "data/winequality-red.csv")
}
wine_data <- as.data.frame(scale(read.csv("data/winequality-red.csv", sep = ";")))
var_names <- colnames(wine_data)
n_vars <- ncol(wine_data)

# Reference DAG
ref_dag <- matrix(0, n_vars, n_vars, dimnames = list(var_names, var_names))
ref_dag["fixed.acidity","pH"] <- 1; ref_dag["volatile.acidity","pH"] <- 1
ref_dag["citric.acid","pH"] <- 1; ref_dag["citric.acid","fixed.acidity"] <- 1
ref_dag["residual.sugar","density"] <- 1; ref_dag["alcohol","density"] <- 1
ref_dag["free.sulfur.dioxide","total.sulfur.dioxide"] <- 1
ref_dag["alcohol","quality"] <- 1; ref_dag["volatile.acidity","quality"] <- 1
ref_dag["sulphates","quality"] <- 1

cat("Reference edges:", sum(ref_dag), "\n\n")

# Evaluation function (same as in the paper's code)
evaluate_dag <- function(est, ref) {
  n <- nrow(ref)
  tp <- sum(est == 1 & ref == 1)
  fp <- sum(est == 1 & ref == 0)
  fn <- sum(est == 0 & ref == 1)
  tn <- sum(est == 0 & ref == 0) - n  # subtract diagonal
  prec <- if ((tp + fp) > 0) tp / (tp + fp) else 0
  rec <- if ((tp + fn) > 0) tp / (tp + fn) else 0
  f1 <- if ((prec + rec) > 0) 2 * prec * rec / (prec + rec) else 0
  shd <- sum(abs(est - ref))
  mse <- mean((est - ref)^2)
  c(F1=f1, Accuracy=(tp + tn) / (n * (n - 1)), SHD=shd, MSE=mse,
    TP=tp, FP=fp, FN=fn, Precision=prec, Recall=rec,
    Edges=sum(est), MSE_check=shd/(n^2))
}

# =============================================================================
# METHOD 1: Original code (hard threshold |B| < 0.1)
# =============================================================================
cat(paste(rep("=",60),collapse=""),"\n")
cat("METHOD 1: LiNGAM with hard threshold |B| < 0.1\n")
cat(paste(rep("=",60),collapse=""),"\n")

set.seed(42)
t0 <- Sys.time()

ica_result <- fastICA(as.matrix(wine_data), n.comp = n_vars,
                      alg.typ = "parallel", fun = "logcosh",
                      method = "C", verbose = FALSE)
W <- ica_result$K %*% ica_result$W  # canonical: include pre-whitening
B <- diag(n_vars) - W                # canonical LiNGAM: B = I - W_full
diag(B) <- 0
B_orig <- B  # save for bootstrap later

cat(sprintf("Non-zero entries in B (before threshold): %d / %d\n",
            sum(abs(B) > 0), n_vars * (n_vars - 1)))
cat(sprintf("Entries with |B| >= 0.1: %d\n", sum(abs(B) >= 0.1)))
cat(sprintf("Entries with |B| >= 0.05: %d\n", sum(abs(B) >= 0.05)))

B[abs(B) < 0.1] <- 0
dag1 <- (abs(B) > 0) * 1
dimnames(dag1) <- list(var_names, var_names)

time1 <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

m1 <- evaluate_dag(dag1, ref_dag)
cat(sprintf("\nResults:\n"))
cat(sprintf("  Edges: %d\n", m1["Edges"]))
cat(sprintf("  TP=%d FP=%d FN=%d\n", m1["TP"], m1["FP"], m1["FN"]))
cat(sprintf("  F1=%.3f Accuracy=%.3f\n", m1["F1"], m1["Accuracy"]))
cat(sprintf("  SHD=%d MSE=%.3f (check: SHD/144=%.3f)\n", m1["SHD"], m1["MSE"], m1["MSE_check"]))
cat(sprintf("  Time: %.1fs\n", time1))

# Print edges
cat("\nEstimated edges:\n")
for (i in 1:n_vars) for (j in 1:n_vars) {
  if (dag1[i,j] == 1) {
    in_ref <- ref_dag[i,j] == 1
    cat(sprintf("  %s -> %s %s\n", var_names[i], var_names[j],
                if(in_ref) "(TP)" else ""))
  }
}

# =============================================================================
# METHOD 2: LiNGAM with bootstrap BH-FDR (as manuscript CLAIMS)
# =============================================================================
cat(paste(rep("\n=",60),collapse=""),"\n")
cat("METHOD 2: LiNGAM with bootstrap BH-FDR (alpha=0.05)\n")
cat(paste(rep("=",60),collapse=""),"\n")

N_BOOT <- 500
cat(sprintf("Running %d bootstrap replicates...\n", N_BOOT))
t0 <- Sys.time()

n_samples <- nrow(wine_data)
B_boot <- array(0, dim = c(N_BOOT, n_vars, n_vars))

set.seed(123)
for (b in 1:N_BOOT) {
  boot_idx <- sample.int(n_samples, replace = TRUE)
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
    B_boot[b,,] <- A_b
  }

  if (b %% 100 == 0) cat(sprintf("  %d/%d done\n", b, N_BOOT))
}

boot_time <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("  Bootstrap done in %.0fs.\n\n", boot_time))

# For each off-diagonal entry B[i,j], test H0: B[i,j] = 0
# p-value = proportion of bootstrap samples where sign of B[i,j] differs
# from the original, or where |B[i,j]| is close to 0
pvals_boot <- matrix(1, n_vars, n_vars, dimnames = list(var_names, var_names))
for (i in 1:n_vars) {
  for (j in 1:n_vars) {
    if (i == j) next
    boot_vals <- B_boot[, i, j]
    boot_vals <- boot_vals[boot_vals != 0]  # remove failed bootstraps
    if (length(boot_vals) < 10) next
    # Two-sided test: proportion of bootstrap values on opposite side of 0
    orig_val <- B_orig[i, j]
    if (abs(orig_val) < 1e-10) { pvals_boot[i,j] <- 1; next }
    # p-value: fraction of bootstraps with opposite sign or zero
    pvals_boot[i,j] <- mean(boot_vals * sign(orig_val) <= 0)
  }
}

# Apply BH-FDR
off_diag <- which(row(pvals_boot) != col(pvals_boot))
pv_vec <- pvals_boot[off_diag]
qv_vec <- p.adjust(pv_vec, method = "BH")
qv_mat <- matrix(1, n_vars, n_vars, dimnames = list(var_names, var_names))
qv_mat[off_diag] <- qv_vec

for (alpha in c(0.01, 0.05, 0.10, 0.20)) {
  dag2 <- (qv_mat < alpha & abs(B_orig) > 0) * 1
  diag(dag2) <- 0
  dimnames(dag2) <- list(var_names, var_names)

  # Orient by sign of B: if B[i,j] > 0, edge j -> i (or we use |B| for direction)
  m2 <- evaluate_dag(dag2, ref_dag)
  cat(sprintf("alpha=%.2f: Edges=%d F1=%.3f Acc=%.3f SHD=%d MSE=%.3f TP=%d FP=%d FN=%d\n",
              alpha, m2["Edges"], m2["F1"], m2["Accuracy"], m2["SHD"], m2["MSE"],
              m2["TP"], m2["FP"], m2["FN"]))
}

# =============================================================================
# CONSISTENCY CHECK
# =============================================================================
cat(paste(rep("\n=",60),collapse=""),"\n")
cat("CONSISTENCY CHECK\n")
cat(paste(rep("=",60),collapse=""),"\n")
cat("For binary matrices: MSE = SHD / p^2 = SHD / 144\n")
cat(sprintf("Paper claims: SHD=56 → MSE should be %.3f (paper says 0.778)\n", 56/144))
cat(sprintf("Paper claims: MSE=0.778 → SHD should be %.0f (paper says 56)\n", 0.778*144))
cat(sprintf("Paper claims: Accuracy=0.152 → errors = %.0f → SHD should be %.0f\n",
            132 - 0.152*132, 132 - 0.152*132))
cat("\nConclusion: MSE=0.778 is consistent with SHD=112, NOT SHD=56.\n")
cat("The SHD=56 in the manuscript is likely wrong.\n")

# =============================================================================
# LaTeX for Table 1 (LiNGAM row)
# =============================================================================
cat(paste(rep("\n=",60),collapse=""),"\n")
cat("LaTeX for Table 1 LiNGAM row (from Method 1, hard threshold):\n")
cat(paste(rep("=",60),collapse=""),"\n")
cat(sprintf("LiNGAM & %.3f & %.3f & %d & %.3f & %.1f \\\\\n",
            m1["F1"], m1["Accuracy"], m1["SHD"], m1["MSE"], time1))

cat("\nDone.\n")
