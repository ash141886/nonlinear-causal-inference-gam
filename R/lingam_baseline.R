# =============================================================================
# LiNGAM Baseline for Comparison
# =============================================================================
#
# ICA-based LiNGAM (Shimizu et al., 2006) with hard threshold |B_ij| >= 0.1.
# Bootstrap BH-FDR is also attempted but yields zero edges due to ICA
# permutation/sign indeterminacy (documented in manuscript Section 4).
#
# Repository: https://github.com/ash141886/nonlinear-causal-inference-gam
# =============================================================================

source("R/utils.R")

#' LiNGAM causal discovery via ICA
#'
#' Estimates the mixing matrix using fastICA, then thresholds to obtain
#' the adjacency matrix. Uses standard threshold |B_ij| >= 0.1.
#'
#' @param data Data frame (n x p), standardized
#' @param threshold Hard threshold for mixing-matrix coefficients (default 0.1)
#' @param n_comp Number of ICA components (default = p)
#' @param verbose Print progress? (default TRUE)
#' @return List with: adjacency, mixing_matrix
lingam_discover <- function(data, threshold = 0.1, n_comp = NULL, verbose = TRUE) {
  data <- as.data.frame(scale(data))
  p <- ncol(data)
  var_names <- colnames(data)

  if (is.null(n_comp)) n_comp <- p

  # Run fastICA
  ica_result <- fastICA::fastICA(as.matrix(data), n.comp = n_comp, method = "C")

  # Estimate mixing matrix W = A^{-1}
  W <- ica_result$K %*% ica_result$W
  A <- solve(W)

  # Causal order: permute rows of W to maximize diagonal
  # Use the LiNGAM ordering heuristic
  B <- diag(p) - solve(A)  # B matrix: X = B X + e

  # Threshold
  adjacency <- matrix(0, p, p, dimnames = list(var_names, var_names))
  adjacency[abs(B) >= threshold] <- 1
  diag(adjacency) <- 0

  if (verbose) {
    cat(sprintf("LiNGAM: %d edges selected (threshold=%.2f)\n", sum(adjacency), threshold))
  }

  list(adjacency = adjacency, mixing_matrix = B)
}

#' LiNGAM with bootstrap BH-FDR (for fairness comparison)
#'
#' Attempts element-wise bootstrap inference on ICA coefficients.
#' Typically yields zero significant edges due to ICA indeterminacy.
#'
#' @param data Data frame (n x p)
#' @param n_boot Number of bootstrap replicates (default 500)
#' @param alpha FDR level (default 0.05)
#' @param verbose Print progress?
#' @return List with: adjacency, n_significant
lingam_bootstrap_bh <- function(data, n_boot = 500, alpha = 0.05, verbose = TRUE) {
  data <- as.data.frame(scale(data))
  n <- nrow(data)
  p <- ncol(data)
  var_names <- colnames(data)

  if (verbose) cat(sprintf("LiNGAM bootstrap: %d replicates...\n", n_boot))

  # Collect bootstrap B matrices
  B_boots <- array(0, dim = c(p, p, n_boot))
  for (b in seq_len(n_boot)) {
    idx <- sample(n, replace = TRUE)
    tryCatch({
      ica <- fastICA::fastICA(as.matrix(data[idx, ]), n.comp = p, method = "C")
      W <- ica$K %*% ica$W
      A <- solve(W)
      B_boots[, , b] <- diag(p) - solve(A)
    }, error = function(e) {
      B_boots[, , b] <<- matrix(NA, p, p)
    })
  }

  # Element-wise p-values (two-sided test: H0: B_ij = 0)
  pvalues <- matrix(1, p, p, dimnames = list(var_names, var_names))
  for (i in seq_len(p)) {
    for (j in seq_len(p)) {
      if (i == j) next
      vals <- B_boots[i, j, ]
      vals <- vals[!is.na(vals)]
      if (length(vals) > 10) {
        # Proportion of sign changes indicates instability
        pvalues[i, j] <- 2 * min(mean(vals > 0), mean(vals < 0))
      }
    }
  }

  # BH-FDR
  off_diag <- pvalues[row(pvalues) != col(pvalues)]
  m <- length(off_diag)
  ordered <- sort(off_diag)
  thresholds <- seq_len(m) * alpha / m
  k_star <- max(c(0, which(ordered <= thresholds)))

  adjacency <- matrix(0, p, p, dimnames = list(var_names, var_names))
  if (k_star > 0) {
    cutoff <- ordered[k_star]
    adjacency[pvalues <= cutoff & row(pvalues) != col(pvalues)] <- 1
  }

  n_sig <- sum(adjacency)
  if (verbose) {
    cat(sprintf("LiNGAM bootstrap BH-FDR (alpha=%.2f): %d significant edges\n", alpha, n_sig))
  }

  list(adjacency = adjacency, n_significant = n_sig)
}
