# =============================================================================
# Utility functions for additive-HSIC causal discovery
# =============================================================================
#
# Provides the shared building blocks used by the main discovery procedure
# and by the analysis scripts:
#   * Gaussian-RBF kernel construction and median-heuristic bandwidth
#   * Biased empirical Hilbert-Schmidt Independence Criterion (HSIC)
#   * Graph-level evaluation metrics (directed F1, accuracy, SHD, MSE)
#   * Data-generating mechanisms for linear and non-linear structural
#     equation models with non-Gaussian (Laplace/uniform) noise
#   * Loaders for the UCI Wine Quality datasets and a literature-based
#     reference DAG for the red-wine variables
#
# All routines are written to be callable stand-alone; no global state is
# assumed beyond the supplied arguments.
# =============================================================================


# -----------------------------------------------------------------------------
# Kernel construction
# -----------------------------------------------------------------------------

#' Gaussian radial-basis-function kernel matrix
#'
#' Computes the n-by-n matrix with entries
#' \eqn{K_{ij} = \exp(-\|x_i - x_j\|_2^2 / (2 \sigma^2))}.
#'
#' @param x Numeric vector or matrix of n rows (samples).
#' @param sigma Positive scalar bandwidth.
#' @return An n-by-n symmetric positive-semidefinite matrix.
rbf_kernel <- function(x, sigma) {
  x <- as.matrix(x)
  D <- as.matrix(dist(x))
  exp(-D^2 / (2 * sigma^2))
}

#' Median-heuristic bandwidth for Gaussian kernels
#'
#' Returns the median of all pairwise Euclidean distances between samples,
#' a standard data-driven choice for the Gaussian-RBF bandwidth.
#'
#' @param x Numeric vector or matrix (n samples).
#' @return Positive scalar; falls back to 1 if the median is non-positive
#'   or non-finite (degenerate inputs).
median_bandwidth <- function(x) {
  x <- as.matrix(x)
  D <- as.matrix(dist(x))
  med <- median(D[lower.tri(D)])
  if (!is.finite(med) || med <= 0) med <- 1
  med
}

#' Centre a kernel matrix
#'
#' Applies the two-sided projection
#' \eqn{H K H} with \eqn{H = I_n - n^{-1} \mathbf{1} \mathbf{1}^\top},
#' which removes the mean from each row and each column.
#'
#' @param K An n-by-n kernel matrix.
#' @return The centred kernel matrix.
center_kernel <- function(K) {
  n <- nrow(K)
  rm <- rowMeans(K)
  K - outer(rm, rep(1, n)) - outer(rep(1, n), rm) + mean(rm)
}


# -----------------------------------------------------------------------------
# HSIC estimation
# -----------------------------------------------------------------------------

#' Biased empirical HSIC
#'
#' Returns \eqn{\widehat{\mathrm{HSIC}}(X, Y) = n^{-2} \mathrm{tr}(K H L H)},
#' the standard biased V-statistic estimator, with Gaussian-RBF kernels
#' whose bandwidths are set by the median heuristic separately for each
#' variable.
#'
#' @param x Numeric vector or matrix (n samples).
#' @param y Numeric vector or matrix (n samples).
#' @return Non-negative scalar.
hsic_biased <- function(x, y) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  n <- nrow(x)
  K <- rbf_kernel(x, median_bandwidth(x))
  L <- rbf_kernel(y, median_bandwidth(y))
  Kc <- center_kernel(K)
  as.numeric(sum(Kc * L) / n^2)
}

#' HSIC from pre-computed kernel matrices
#'
#' Convenience wrapper that avoids recomputing K and L in tight permutation
#' loops. The first argument must already be centred; the second is used
#' as-is.
#'
#' @param Kc Centred kernel matrix for the first variable.
#' @param L Kernel matrix for the second variable.
#' @param n Sample size.
#' @return Non-negative scalar.
hsic_from_kernels <- function(Kc, L, n) {
  as.numeric(sum(Kc * L) / n^2)
}


# -----------------------------------------------------------------------------
# Graph-level evaluation metrics
# -----------------------------------------------------------------------------

#' Compare an estimated DAG against a reference DAG
#'
#' Reports the standard directed-edge metrics together with a scalar
#' structural Hamming distance and an adjacency mean-squared error.
#' A reversed edge is counted as one false positive and one false negative
#' and therefore contributes 2 to the SHD.
#'
#' @param est_dag p-by-p binary adjacency matrix (estimated).
#' @param ref_dag p-by-p binary adjacency matrix (reference).
#' @return A named numeric vector: \code{F1}, \code{Accuracy}, \code{SHD},
#'   \code{MSE}, \code{TP}, \code{FP}, \code{FN}, \code{Misoriented}.
#'   \code{Misoriented} counts edges that are present in the estimate but
#'   reversed relative to the reference; it is reported for diagnostic use
#'   and is already accounted for inside \code{SHD}.
evaluate_dag <- function(est_dag, ref_dag) {
  stopifnot(all(dim(est_dag) == dim(ref_dag)))
  p <- nrow(est_dag)

  tp <- sum(est_dag == 1 & ref_dag == 1)
  fp <- sum(est_dag == 1 & ref_dag == 0)
  fn <- sum(est_dag == 0 & ref_dag == 1)

  precision <- if (tp + fp > 0) tp / (tp + fp) else 0
  recall    <- if (tp + fn > 0) tp / (tp + fn) else 0
  f1        <- if (precision + recall > 0) 2 * precision * recall / (precision + recall) else 0

  # Accuracy: fraction of correctly classified ordered off-diagonal pairs.
  n_pairs  <- p * (p - 1)
  correct  <- sum(est_dag == ref_dag) - p      # subtract the p diagonal matches
  accuracy <- correct / n_pairs

  # Structural Hamming distance: cell-count convention.
  shd <- sum(abs(est_dag - ref_dag))

  # Misoriented edges: present in the estimate, reversed in the reference.
  misoriented <- 0
  for (i in seq_len(p)) {
    for (j in seq_len(p)) {
      if (i != j && est_dag[i, j] == 1 && ref_dag[i, j] == 0 && ref_dag[j, i] == 1) {
        misoriented <- misoriented + 1
      }
    }
  }

  mse <- mean((est_dag - ref_dag)^2)

  c(F1 = f1, Accuracy = accuracy, SHD = shd, MSE = mse,
    TP = tp, FP = fp, FN = fn, Misoriented = misoriented)
}


# -----------------------------------------------------------------------------
# Data loading
# -----------------------------------------------------------------------------

#' Load and standardise a UCI Wine Quality dataset
#'
#' Reads the semicolon-separated CSV from \code{data_dir}, downloading it
#' from the UCI Machine Learning Repository on first use, and standardises
#' every column to mean zero and unit variance. The resulting frame is used
#' by the wine-analysis scripts.
#'
#' @param color Either "red" or "white".
#' @param data_dir Destination directory for cached data files.
#' @return A data frame with n rows and p standardised numeric columns.
load_wine_data <- function(color = "red", data_dir = "data") {
  filename <- paste0("winequality-", color, ".csv")
  filepath <- file.path(data_dir, filename)

  if (!file.exists(filepath)) {
    dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
    url <- paste0(
      "https://archive.ics.uci.edu/ml/machine-learning-databases/wine-quality/",
      filename
    )
    download.file(url, filepath)
    cat("Downloaded:", filepath, "\n")
  }

  wine <- read.csv(filepath, sep = ";")
  wine <- as.data.frame(scale(wine))

  cat(sprintf("Wine (%s) loaded: n=%d, p=%d\n", color, nrow(wine), ncol(wine)))
  cat("Variables:", paste(colnames(wine), collapse = ", "), "\n\n")
  wine
}

#' Literature-based reference DAG for the red-wine variables
#'
#' Encodes well-established chemical relationships: titratable acids drive
#' pH; citric acid contributes to fixed acidity; residual sugar and alcohol
#' jointly determine density; free sulphur dioxide is a component of total
#' sulphur dioxide; and alcohol, volatile acidity, and sulphates influence
#' the sensory quality score. The returned matrix is intended only as a
#' post-hoc comparator and not as ground truth.
#'
#' @param var_names Character vector of variable names (order preserved).
#' @return A p-by-p binary adjacency matrix with row/column names.
get_wine_reference_dag <- function(var_names) {
  p <- length(var_names)
  dag <- matrix(0, p, p, dimnames = list(var_names, var_names))

  # Acids -> pH
  dag["fixed.acidity",   "pH"] <- 1
  dag["volatile.acidity","pH"] <- 1
  dag["citric.acid",     "pH"] <- 1

  # Citric acid -> fixed acidity
  dag["citric.acid", "fixed.acidity"] <- 1

  # Mass-balance relationships on density
  dag["residual.sugar", "density"] <- 1
  dag["alcohol",        "density"] <- 1

  # Sulphur dioxide: free is a component of total
  dag["free.sulfur.dioxide", "total.sulfur.dioxide"] <- 1

  # Sensory quality
  dag["alcohol",         "quality"] <- 1
  dag["volatile.acidity","quality"] <- 1
  dag["sulphates",       "quality"] <- 1

  dag
}


# -----------------------------------------------------------------------------
# Random DAG generation
# -----------------------------------------------------------------------------

#' Generate a random DAG
#'
#' Samples a random permutation of the p nodes (treated as a topological
#' order) and then uniformly retains a target fraction of the upper-
#' triangular entries as directed edges.
#'
#' @param p Number of variables.
#' @param edge_density Target fraction of the p(p-1)/2 possible edges.
#' @return A p-by-p binary adjacency matrix encoding a DAG.
random_dag <- function(p, edge_density = 0.3) {
  n_possible <- p * (p - 1) / 2
  n_edges    <- round(n_possible * edge_density)

  order <- sample(p)
  dag   <- matrix(0, p, p)

  possible <- which(upper.tri(dag), arr.ind = TRUE)
  if (nrow(possible) > 0 && n_edges > 0) {
    selected <- possible[sample(nrow(possible), min(n_edges, nrow(possible))), , drop = FALSE]
    for (r in seq_len(nrow(selected))) {
      i <- order[selected[r, 1]]
      j <- order[selected[r, 2]]
      dag[i, j] <- 1
    }
  }
  dag
}


# -----------------------------------------------------------------------------
# Structural equation model data generation
# -----------------------------------------------------------------------------

#' Simulate from a non-linear additive SEM with non-Gaussian noise
#'
#' Traverses the DAG in topological order. For each node, a random subset
#' of its parent effects is routed through a non-linear link drawn from a
#' fixed library (sin, cos, exponential, quadratic-sign, log-sign, tanh);
#' the remaining parent effects are linear. Each edge receives an
#' independent sign and a uniformly sampled magnitude. The noise term for
#' every node is drawn, with equal probability, from a centred Laplace
#' distribution or a centred uniform distribution, each rescaled to the
#' target standard deviation.
#'
#' @param n Sample size.
#' @param dag p-by-p binary adjacency matrix of the generating DAG.
#' @param nonlinear_frac Fraction of parent effects made non-linear.
#' @param noise_sd Target noise standard deviation.
#' @return A data frame with n rows and p columns, named \code{X1, ..., Xp}.
generate_sem_data <- function(n, dag, nonlinear_frac = 0.5, noise_sd = 1.0) {
  p <- nrow(dag)
  X <- matrix(0, n, p)

  topo <- topological_sort(dag)

  nl_funs <- list(
    function(x) sin(x),
    function(x) cos(x),
    function(x) exp(-abs(x)) * sign(x),
    function(x) x^2 * sign(x) * 0.5,
    function(x) log(abs(x) + 1) * sign(x),
    function(x) tanh(x)
  )

  for (i in topo) {
    parents <- which(dag[, i] == 1)

    # Non-Gaussian noise: equal-probability mixture of Laplace and uniform.
    noise <- if (runif(1) > 0.5) {
      rlaplace(n, scale = noise_sd / sqrt(2))
    } else {
      runif(n, -noise_sd * sqrt(3), noise_sd * sqrt(3))
    }

    if (length(parents) == 0) {
      X[, i] <- noise
    } else {
      effect <- rep(0, n)
      for (idx in seq_along(parents)) {
        j    <- parents[idx]
        coef <- runif(1, 0.5, 1.5) * sample(c(-1, 1), 1)
        if (idx <= ceiling(length(parents) * nonlinear_frac)) {
          fn     <- nl_funs[[sample(length(nl_funs), 1)]]
          effect <- effect + coef * fn(X[, j])
        } else {
          effect <- effect + coef * X[, j]
        }
      }
      X[, i] <- effect + noise
    }
  }

  colnames(X) <- paste0("X", seq_len(p))
  as.data.frame(X)
}

#' Simulate from a linear SEM with non-Gaussian noise
#'
#' Thin wrapper around \code{generate_sem_data} with
#' \code{nonlinear_frac = 0}, used by the linear-data diagnostics.
#'
#' @param n Sample size.
#' @param dag Adjacency matrix of the generating DAG.
#' @param noise_sd Target noise standard deviation.
#' @return A data frame with n rows and p columns.
generate_linear_sem_data <- function(n, dag, noise_sd = 1.0) {
  generate_sem_data(n, dag, nonlinear_frac = 0.0, noise_sd = noise_sd)
}

#' Topological sort of a DAG
#'
#' Classical in-degree-zero peeling (Kahn's algorithm). Throws an error if
#' the graph contains a directed cycle.
#'
#' @param dag p-by-p binary adjacency matrix.
#' @return An integer vector of node indices in a valid topological order.
topological_sort <- function(dag) {
  p         <- nrow(dag)
  in_degree <- colSums(dag)
  order     <- integer(0)
  remaining <- seq_len(p)

  while (length(remaining) > 0) {
    roots <- remaining[in_degree[remaining] == 0]
    if (length(roots) == 0) stop("Graph contains a cycle")
    node      <- roots[1]
    order     <- c(order, node)
    remaining <- setdiff(remaining, node)
    children  <- which(dag[node, ] == 1)
    in_degree[children] <- in_degree[children] - 1
  }
  order
}

#' Draw from a Laplace (double-exponential) distribution
#'
#' Inverse-CDF sampling; avoids any external-package dependency.
#'
#' @param n Sample size.
#' @param location Location parameter.
#' @param scale Positive scale parameter.
#' @return A numeric vector of length n.
rlaplace <- function(n, location = 0, scale = 1) {
  u <- runif(n, -0.5, 0.5)
  location - scale * sign(u) * log(1 - 2 * abs(u))
}
