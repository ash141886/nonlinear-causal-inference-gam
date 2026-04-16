# =============================================================================
# Utility Functions for Additive-HSIC Causal Discovery
# =============================================================================
#
# Shared helper functions: kernel computation, HSIC estimation, evaluation
# metrics, and data loading.
#
# Reference:
#   Islam, M.A. and Suzuki, J. "Non-linear Causal Inference in Observational
#   Data using Additive Models and Kernel-based Independence Testing"
#
# Repository: https://github.com/ash141886/nonlinear-causal-inference-gam
# =============================================================================

# --- Kernel Functions --------------------------------------------------------

#' Gaussian RBF kernel matrix
#'
#' @param x Numeric vector or matrix (n x d)
#' @param sigma Bandwidth parameter
#' @return n x n kernel matrix
rbf_kernel <- function(x, sigma) {
  x <- as.matrix(x)
  D <- as.matrix(dist(x))
  exp(-D^2 / (2 * sigma^2))
}

#' Median heuristic bandwidth
#'
#' Sets bandwidth to the median of all pairwise Euclidean distances.
#' See Gretton et al. (2008) for justification.
#'
#' @param x Numeric vector or matrix
#' @return Scalar bandwidth
median_bandwidth <- function(x) {
  x <- as.matrix(x)
  D <- as.matrix(dist(x))
  med <- median(D[lower.tri(D)])
  if (!is.finite(med) || med <= 0) med <- 1
  med
}

#' Center a kernel matrix
#'
#' Applies H K H where H = I_n - (1/n) 1 1^T
#'
#' @param K n x n kernel matrix
#' @return Centered kernel matrix
center_kernel <- function(K) {
  n <- nrow(K)
  rm <- rowMeans(K)
  K - outer(rm, rep(1, n)) - outer(rep(1, n), rm) + mean(rm)
}

# --- HSIC Estimation ---------------------------------------------------------

#' Biased empirical HSIC
#'
#' Computes (1/n^2) * tr(K H L H) using Gaussian RBF kernels
#' with median-heuristic bandwidth.
#'
#' @param x Numeric vector or matrix (n samples)
#' @param y Numeric vector or matrix (n samples)
#' @return Scalar HSIC value
hsic_biased <- function(x, y) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  n <- nrow(x)
  K <- rbf_kernel(x, median_bandwidth(x))
  L <- rbf_kernel(y, median_bandwidth(y))
  Kc <- center_kernel(K)
  as.numeric(sum(Kc * L) / n^2)
}

#' HSIC from precomputed kernel matrices
#'
#' @param Kc Centered kernel matrix for first variable
#' @param L  Kernel matrix for second variable (uncentered)
#' @param n  Sample size
#' @return Scalar HSIC value
hsic_from_kernels <- function(Kc, L, n) {
  as.numeric(sum(Kc * L) / n^2)
}

# --- Evaluation Metrics ------------------------------------------------------

#' Evaluate estimated DAG against a reference
#'
#' Computes directed F1, graph accuracy, SHD, and adjacency MSE.
#'
#' @param est_dag Estimated adjacency matrix (p x p, binary)
#' @param ref_dag Reference adjacency matrix (p x p, binary)
#' @return Named vector: F1, Accuracy, SHD, MSE
evaluate_dag <- function(est_dag, ref_dag) {
  stopifnot(all(dim(est_dag) == dim(ref_dag)))
  p <- nrow(est_dag)

  # Directed: TP, FP, FN
  tp <- sum(est_dag == 1 & ref_dag == 1)
  fp <- sum(est_dag == 1 & ref_dag == 0)
  fn <- sum(est_dag == 0 & ref_dag == 1)

  precision <- if (tp + fp > 0) tp / (tp + fp) else 0
  recall    <- if (tp + fn > 0) tp / (tp + fn) else 0
  f1        <- if (precision + recall > 0) 2 * precision * recall / (precision + recall) else 0

  # Graph accuracy: fraction of correctly classified ordered pairs
  n_pairs <- p * (p - 1)
  correct <- sum(est_dag == ref_dag) - p  # exclude diagonal
  accuracy <- correct / n_pairs

  # SHD: cell-count convention (sum of per-cell differences)
  # This matches the formula used for the wine analysis (evaluate_dag_wine)
  # and the manuscript-reported values. A reversed edge contributes 2.
  extra   <- sum(est_dag == 1 & ref_dag == 0)
  missing <- sum(est_dag == 0 & ref_dag == 1)
  shd     <- sum(abs(est_dag - ref_dag))
  # Diagnostic: count reversed-orientation edges (for reporting only,
  # NOT added to SHD; already captured via extra + missing).
  misoriented <- 0
  for (i in seq_len(p)) {
    for (j in seq_len(p)) {
      if (i != j && est_dag[i, j] == 1 && ref_dag[i, j] == 0 && ref_dag[j, i] == 1) {
        misoriented <- misoriented + 1
      }
    }
  }

  # MSE: mean squared difference
  mse <- mean((est_dag - ref_dag)^2)

  c(F1 = f1, Accuracy = accuracy, SHD = shd, MSE = mse,
    TP = tp, FP = fp, FN = fn, Misoriented = misoriented)
}

# --- Data Loading ------------------------------------------------------------

#' Load and standardize UCI Wine Quality data
#'
#' Downloads from UCI repository if not present locally.
#'
#' @param color "red" or "white"
#' @param data_dir Directory for data files
#' @return Standardized data frame
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

#' Construct literature-based reference DAG for wine data
#'
#' Based on Jackson (2020) and Waterhouse et al. (2024).
#'
#' @param var_names Variable names from data
#' @return Binary adjacency matrix
get_wine_reference_dag <- function(var_names) {
  p <- length(var_names)
  dag <- matrix(0, p, p, dimnames = list(var_names, var_names))

  # Acid -> pH
  dag["fixed.acidity", "pH"] <- 1
  dag["volatile.acidity", "pH"] <- 1
  dag["citric.acid", "pH"] <- 1

  # Citric acid -> Fixed acidity
  dag["citric.acid", "fixed.acidity"] <- 1

  # Density relationships
  dag["residual.sugar", "density"] <- 1
  dag["alcohol", "density"] <- 1

  # SO2
  dag["free.sulfur.dioxide", "total.sulfur.dioxide"] <- 1

  # Quality
  dag["alcohol", "quality"] <- 1
  dag["volatile.acidity", "quality"] <- 1
  dag["sulphates", "quality"] <- 1

  dag
}

# --- Random DAG Generation ---------------------------------------------------

#' Generate a random DAG adjacency matrix
#'
#' @param p Number of variables
#' @param edge_density Target edge density (fraction of possible edges)
#' @return Binary adjacency matrix (lower-triangular in topological order)
random_dag <- function(p, edge_density = 0.3) {
  n_possible <- p * (p - 1) / 2
  n_edges <- round(n_possible * edge_density)

  # Random topological ordering
  order <- sample(p)
  dag <- matrix(0, p, p)

  # Add edges from earlier to later in ordering
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

# --- Simulation Data Generation ----------------------------------------------

#' Generate data from a non-linear additive SEM
#'
#' @param n Sample size
#' @param dag Adjacency matrix (p x p)
#' @param nonlinear_frac Fraction of parent effects that are non-linear
#' @param noise_sd Noise standard deviation
#' @return Data frame with p columns
generate_sem_data <- function(n, dag, nonlinear_frac = 0.5, noise_sd = 1.0) {
  p <- nrow(dag)
  X <- matrix(0, n, p)

  # Topological sort
  topo <- topological_sort(dag)

  # Non-linear transformations
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
    # Non-Gaussian noise: mixture of Laplace and uniform
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
        j <- parents[idx]
        coef <- runif(1, 0.5, 1.5) * sample(c(-1, 1), 1)
        if (idx <= ceiling(length(parents) * nonlinear_frac)) {
          # Non-linear effect
          fn <- nl_funs[[sample(length(nl_funs), 1)]]
          effect <- effect + coef * fn(X[, j])
        } else {
          # Linear effect
          effect <- effect + coef * X[, j]
        }
      }
      X[, i] <- effect + noise
    }
  }

  colnames(X) <- paste0("X", seq_len(p))
  as.data.frame(X)
}

#' Generate data from a linear SEM (for Section 4.2)
#'
#' @param n Sample size
#' @param dag Adjacency matrix
#' @param noise_sd Noise standard deviation
#' @return Data frame
generate_linear_sem_data <- function(n, dag, noise_sd = 1.0) {
  generate_sem_data(n, dag, nonlinear_frac = 0.0, noise_sd = noise_sd)
}

#' Topological sort of a DAG
#'
#' @param dag Adjacency matrix
#' @return Integer vector of node indices in topological order
topological_sort <- function(dag) {
  p <- nrow(dag)
  in_degree <- colSums(dag)
  order <- integer(0)
  remaining <- seq_len(p)

  while (length(remaining) > 0) {
    roots <- remaining[in_degree[remaining] == 0]
    if (length(roots) == 0) stop("Graph contains a cycle")
    node <- roots[1]
    order <- c(order, node)
    remaining <- setdiff(remaining, node)
    children <- which(dag[node, ] == 1)
    in_degree[children] <- in_degree[children] - 1
  }
  order
}

#' Laplace random numbers
#'
#' @param n Number of samples
#' @param location Location parameter
#' @param scale Scale parameter
#' @return Numeric vector
rlaplace <- function(n, location = 0, scale = 1) {
  u <- runif(n, -0.5, 0.5)
  location - scale * sign(u) * log(1 - 2 * abs(u))
}
