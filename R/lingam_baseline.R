# =============================================================================
# Linear non-Gaussian acyclic model (LiNGAM) baseline
# =============================================================================
#
# Implements the ICA-based LiNGAM estimator of Shimizu et al. (2006) for use
# as a linear benchmark against the additive-HSIC procedure. The estimator
# fits independent components on standardised data with \code{fastICA},
# converts the mixing matrix to the structural coefficient matrix
# \eqn{B = I - A^{-1}}, and applies a hard threshold to obtain an
# adjacency matrix.
#
# A bootstrap Benjamini-Hochberg variant is also provided for comparison;
# in practice it tends to return very few edges because the ICA solution
# is identified only up to permutation and sign of its components, which
# inflates the apparent element-wise variability of the bootstrapped
# \eqn{B} entries.
# =============================================================================

source("R/utils.R")

#' LiNGAM causal discovery via independent component analysis
#'
#' Standardises the input, fits a p-component FastICA, converts the
#' resulting mixing matrix to the LiNGAM structural coefficient matrix
#' \eqn{B = I - A^{-1}}, and retains the entries whose absolute value
#' exceeds \code{threshold} as directed edges.
#'
#' @param data Data frame with n rows and p numeric columns.
#' @param threshold Hard cut-off on \eqn{|B_{ij}|}; default 0.1.
#' @param n_comp Number of ICA components. Defaults to \code{ncol(data)}.
#' @param verbose Logical; whether to print a summary line.
#' @return A list with elements \code{adjacency} (p-by-p binary matrix)
#'   and \code{mixing_matrix} (the estimated \eqn{B}).
lingam_discover <- function(data, threshold = 0.1, n_comp = NULL, verbose = TRUE) {
  data <- as.data.frame(scale(data))
  p          <- ncol(data)
  var_names  <- colnames(data)

  if (is.null(n_comp)) n_comp <- p

  ica_result <- fastICA::fastICA(as.matrix(data), n.comp = n_comp, method = "C")

  # W = K %*% W_rot is the whitened ICA demixing matrix; invert to obtain A.
  W <- ica_result$K %*% ica_result$W
  A <- solve(W)

  # Structural coefficient matrix B in the LiNGAM parameterisation X = B X + e.
  B <- diag(p) - solve(A)

  adjacency <- matrix(0, p, p, dimnames = list(var_names, var_names))
  adjacency[abs(B) >= threshold] <- 1
  diag(adjacency) <- 0

  if (verbose) {
    cat(sprintf("LiNGAM: %d edges selected (threshold=%.2f)\n",
                sum(adjacency), threshold))
  }

  list(adjacency = adjacency, mixing_matrix = B)
}

#' LiNGAM with bootstrap-based Benjamini-Hochberg inference
#'
#' Repeatedly resamples the rows of \code{data}, re-estimates the LiNGAM
#' coefficient matrix on each replicate, converts the sign-stability of
#' each off-diagonal entry into a two-sided element-wise p-value, and
#' finally selects edges by Benjamini-Hochberg FDR control at level
#' \code{alpha}.
#'
#' The procedure is offered only for completeness: because the ICA
#' decomposition is identified up to signed permutations of its
#' components, the bootstrapped \eqn{B} entries exhibit high apparent
#' variance and the BH selection is therefore typically empty, even when
#' the underlying linear structure is well identified.
#'
#' @param data Data frame with n rows and p numeric columns.
#' @param n_boot Number of bootstrap replicates.
#' @param alpha FDR level.
#' @param verbose Logical; whether to print a summary line.
#' @return A list with elements \code{adjacency} (p-by-p binary matrix)
#'   and \code{n_significant} (integer).
lingam_bootstrap_bh <- function(data, n_boot = 500, alpha = 0.05, verbose = TRUE) {
  data <- as.data.frame(scale(data))
  n          <- nrow(data)
  p          <- ncol(data)
  var_names  <- colnames(data)

  if (verbose) cat(sprintf("LiNGAM bootstrap: %d replicates...\n", n_boot))

  # Store one p-by-p coefficient matrix per bootstrap replicate.
  B_boots <- array(0, dim = c(p, p, n_boot))
  for (b in seq_len(n_boot)) {
    idx <- sample(n, replace = TRUE)
    tryCatch({
      ica <- fastICA::fastICA(as.matrix(data[idx, ]), n.comp = p, method = "C")
      W   <- ica$K %*% ica$W
      A   <- solve(W)
      B_boots[, , b] <- diag(p) - solve(A)
    }, error = function(e) {
      B_boots[, , b] <<- matrix(NA, p, p)
    })
  }

  # Two-sided sign-stability p-values for each off-diagonal entry.
  pvalues <- matrix(1, p, p, dimnames = list(var_names, var_names))
  for (i in seq_len(p)) {
    for (j in seq_len(p)) {
      if (i == j) next
      vals <- B_boots[i, j, ]
      vals <- vals[!is.na(vals)]
      if (length(vals) > 10) {
        pvalues[i, j] <- 2 * min(mean(vals > 0), mean(vals < 0))
      }
    }
  }

  # Benjamini-Hochberg FDR selection on the off-diagonal p-values.
  off_diag   <- pvalues[row(pvalues) != col(pvalues)]
  m          <- length(off_diag)
  ordered    <- sort(off_diag)
  thresholds <- seq_len(m) * alpha / m
  k_star     <- max(c(0, which(ordered <= thresholds)))

  adjacency <- matrix(0, p, p, dimnames = list(var_names, var_names))
  if (k_star > 0) {
    cutoff <- ordered[k_star]
    adjacency[pvalues <= cutoff & row(pvalues) != col(pvalues)] <- 1
  }

  n_sig <- sum(adjacency)
  if (verbose) {
    cat(sprintf("LiNGAM bootstrap BH-FDR (alpha=%.2f): %d significant edges\n",
                alpha, n_sig))
  }

  list(adjacency = adjacency, n_significant = n_sig)
}
