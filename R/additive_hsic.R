# =============================================================================
# Core Method: Additive-HSIC Causal Discovery
# =============================================================================
#
# Implementation of Algorithm 1 from:
#   Islam, M.A. and Suzuki, J. "Non-linear Causal Inference in Observational
#   Data using Additive Models and Kernel-based Independence Testing"
#
# Repository: https://github.com/ash141886/nonlinear-causal-inference-gam
# =============================================================================

source("R/utils.R")

# --- Leave-One-Out Additive Model Fitting ------------------------------------

#' Fit reduced additive model excluding one variable
#'
#' Fits X_i ~ sum_{k != i, k != j} f_ik(X_k) using mgcv::gam with REML.
#' Returns residuals r_i^{(-j)}.
#'
#' @param data Data frame (n x p)
#' @param i Response variable index
#' @param j Excluded variable index
#' @param k_max Maximum basis dimension for each smooth
#' @return Numeric vector of residuals (length n)
fit_reduced_model <- function(data, i, j, k_max = NULL) {
  n <- nrow(data)
  p <- ncol(data)
  var_names <- colnames(data)

  if (is.null(k_max)) k_max <- min(floor(n / 4), 40)

  # Predictors: all except i and j
  pred_idx <- setdiff(seq_len(p), c(i, j))

  if (length(pred_idx) == 0) {
    return(data[[i]] - mean(data[[i]]))
  }

  # Build formula: X_i ~ s(X_k1) + s(X_k2) + ...
  smooth_terms <- paste0("s(", var_names[pred_idx], ", k=", k_max, ")")
  formula_str <- paste(var_names[i], "~", paste(smooth_terms, collapse = " + "))
  formula_obj <- as.formula(formula_str)

  # Fit GAM with REML and gamma=1.4 for extra smoothing
  fit <- mgcv::gam(formula_obj, data = data, method = "REML",
                    gamma = 1.4, select = TRUE)

  as.numeric(residuals(fit))
}

# --- Main Discovery Algorithm ------------------------------------------------

#' Additive-HSIC Causal Discovery (Algorithm 1)
#'
#' Two modes:
#'   - Single-stage: BH-FDR at level alpha (default for simulations)
#'   - Two-stage: BH-FDR screening at liberal alpha + top-k by HSIC
#'                (recommended for real data with large n/p ratio)
#'
#' @param data Data frame (n x p), standardized
#' @param alpha FDR level (default 0.05)
#' @param n_perm Number of permutations (default 200)
#' @param k_max Maximum GAM basis dimension (NULL = auto)
#' @param two_stage Use two-stage selection? (default FALSE)
#' @param screening_alpha FDR level for screening stage (default 0.20)
#' @param target_density Target edge density for top-k (default 0.15)
#' @param n_sub Subsample size for permutation test (NULL = use full data)
#' @param seed Random seed for reproducibility
#' @param verbose Print progress? (default TRUE)
#' @return List with: adjacency, scores, pvalues, edges, selected_edges
discover_dag <- function(data, alpha = 0.05, n_perm = 200, k_max = NULL,
                         two_stage = FALSE, screening_alpha = 0.20,
                         target_density = 0.15, n_sub = NULL,
                         seed = 42, verbose = TRUE) {
  set.seed(seed)
  n <- nrow(data)
  p <- ncol(data)
  var_names <- colnames(data)

  if (is.null(k_max)) k_max <- min(floor(n / 4), 40)

  # Step 1-2: Standardize and initialize
  data <- as.data.frame(scale(data))
  M <- matrix(0, p, p, dimnames = list(var_names, var_names))

  # Steps 3-9: Compute HSIC scores for all ordered pairs
  if (verbose) cat("Computing leave-one-out residuals and HSIC scores...\n")
  resid_cache <- list()

  for (i in seq_len(p)) {
    if (verbose) cat(sprintf("  Variable %d/%d (%s)\n", i, p, var_names[i]))
    for (j in seq_len(p)) {
      if (i == j) next

      # Fit reduced model and cache residuals
      key <- paste(i, j, sep = "_")
      r_ij <- fit_reduced_model(data, i, j, k_max)
      resid_cache[[key]] <- r_ij

      # Compute HSIC score on full data
      M[i, j] <- hsic_biased(r_ij, data[[j]])
    }
  }

  # Step 10: Compute T_ij and permutation p-values
  if (verbose) cat("Computing permutation p-values...\n")
  m_tests <- p * (p - 1) / 2
  pairs <- combn(p, 2)
  T_obs <- numeric(ncol(pairs))
  pvalues <- numeric(ncol(pairs))
  pair_names <- character(ncol(pairs))

  for (idx in seq_len(ncol(pairs))) {
    i <- pairs[1, idx]
    j <- pairs[2, idx]
    T_obs[idx] <- max(M[i, j], M[j, i])
    pair_names[idx] <- paste(var_names[i], var_names[j], sep = "--")

    # Permutation test
    r_ij <- resid_cache[[paste(i, j, sep = "_")]]
    r_ji <- resid_cache[[paste(j, i, sep = "_")]]

    # Subsample for permutation if requested
    if (!is.null(n_sub) && n_sub < n) {
      sub_idx <- sample(n, n_sub)
      x_j_perm <- data[[j]][sub_idx]
      x_i_perm <- data[[i]][sub_idx]
      r_ij_perm <- r_ij[sub_idx]
      r_ji_perm <- r_ji[sub_idx]
    } else {
      x_j_perm <- data[[j]]
      x_i_perm <- data[[i]]
      r_ij_perm <- r_ij
      r_ji_perm <- r_ji
    }

    T_null <- numeric(n_perm)
    for (b in seq_len(n_perm)) {
      perm <- sample(length(x_j_perm))
      h1 <- hsic_biased(r_ij_perm, x_j_perm[perm])
      h2 <- hsic_biased(r_ji_perm, x_i_perm[perm])
      T_null[b] <- max(h1, h2)
    }

    pvalues[idx] <- (1 + sum(T_null >= T_obs[idx])) / (n_perm + 1)
  }

  names(pvalues) <- pair_names
  names(T_obs) <- pair_names

  # Step 11: Edge selection
  if (two_stage) {
    # Two-stage: BH screening + top-k
    effective_alpha <- screening_alpha
    k_target <- floor(m_tests * target_density)
    if (verbose) cat(sprintf("Two-stage: BH at alpha=%.2f, then top-%d by HSIC\n",
                             effective_alpha, k_target))
  } else {
    effective_alpha <- alpha
  }

  # BH procedure
  ordered_idx <- order(pvalues)
  thresholds <- seq_len(length(pvalues)) * effective_alpha / length(pvalues)
  rejected <- pvalues[ordered_idx] <= thresholds
  k_star <- max(c(0, which(rejected)))

  if (k_star > 0) {
    selected_idx <- ordered_idx[seq_len(k_star)]
  } else {
    selected_idx <- integer(0)
  }

  # Two-stage: further filter by top-k HSIC
  if (two_stage && length(selected_idx) > k_target) {
    hsic_ranks <- order(T_obs[selected_idx], decreasing = TRUE)
    selected_idx <- selected_idx[hsic_ranks[seq_len(k_target)]]
  }

  # Steps 12-15: Orient edges
  adjacency <- matrix(0, p, p, dimnames = list(var_names, var_names))
  selected_edges <- data.frame(
    parent = character(0), child = character(0),
    hsic_score = numeric(0), pvalue = numeric(0),
    stringsAsFactors = FALSE
  )

  for (idx in selected_idx) {
    i <- pairs[1, idx]
    j <- pairs[2, idx]
    if (M[i, j] > M[j, i]) {
      # j -> i (j is parent of i)
      adjacency[j, i] <- 1
      selected_edges <- rbind(selected_edges, data.frame(
        parent = var_names[j], child = var_names[i],
        hsic_score = M[i, j], pvalue = pvalues[idx],
        stringsAsFactors = FALSE
      ))
    } else {
      # i -> j
      adjacency[i, j] <- 1
      selected_edges <- rbind(selected_edges, data.frame(
        parent = var_names[i], child = var_names[j],
        hsic_score = M[j, i], pvalue = pvalues[idx],
        stringsAsFactors = FALSE
      ))
    }
  }

  # Steps 16-18: Break cycles
  adjacency <- break_cycles(adjacency, M)

  if (verbose) {
    cat(sprintf("Edges selected: %d\n", sum(adjacency)))
    cat(sprintf("Pairs tested: %d, BH-rejected: %d\n", length(pvalues), k_star))
  }

  list(
    adjacency = adjacency,
    scores = M,
    pvalues = pvalues,
    T_obs = T_obs,
    edges = selected_edges,
    resid_cache = resid_cache
  )
}

# --- Cycle Breaking ----------------------------------------------------------

#' Break directed cycles by removing the weakest edge
#'
#' Iteratively finds cycles and removes the edge with smallest HSIC score.
#'
#' @param adj Adjacency matrix
#' @param M HSIC score matrix
#' @return Acyclic adjacency matrix
break_cycles <- function(adj, M) {
  max_iter <- sum(adj)
  for (iter in seq_len(max_iter)) {
    cycle <- find_cycle(adj)
    if (is.null(cycle)) break

    # Find weakest edge in cycle
    min_score <- Inf
    min_edge <- NULL
    for (k in seq_along(cycle)) {
      from <- cycle[k]
      to <- cycle[if (k < length(cycle)) k + 1 else 1]
      if (adj[from, to] == 1 && M[to, from] < min_score) {
        min_score <- M[to, from]
        min_edge <- c(from, to)
      }
    }

    if (!is.null(min_edge)) {
      adj[min_edge[1], min_edge[2]] <- 0
    } else {
      break
    }
  }
  adj
}

#' Find a directed cycle in an adjacency matrix
#'
#' @param adj Adjacency matrix
#' @return Vector of node indices forming a cycle, or NULL
find_cycle <- function(adj) {
  p <- nrow(adj)
  visited <- rep(FALSE, p)
  in_stack <- rep(FALSE, p)
  parent <- rep(0, p)

  for (start in seq_len(p)) {
    if (visited[start]) next
    result <- dfs_cycle(adj, start, visited, in_stack, parent)
    if (!is.null(result$cycle)) return(result$cycle)
    visited <- result$visited
    in_stack <- result$in_stack
  }
  NULL
}

dfs_cycle <- function(adj, node, visited, in_stack, parent) {
  visited[node] <- TRUE
  in_stack[node] <- TRUE

  children <- which(adj[node, ] == 1)
  for (child in children) {
    if (!visited[child]) {
      parent[child] <- node
      result <- dfs_cycle(adj, child, visited, in_stack, parent)
      if (!is.null(result$cycle)) return(result)
      visited <- result$visited
      in_stack <- result$in_stack
    } else if (in_stack[child]) {
      # Found cycle: trace back
      cycle <- c(child)
      curr <- node
      while (curr != child) {
        cycle <- c(curr, cycle)
        curr <- parent[curr]
      }
      return(list(cycle = cycle, visited = visited, in_stack = in_stack))
    }
  }

  in_stack[node] <- FALSE
  list(cycle = NULL, visited = visited, in_stack = in_stack)
}
