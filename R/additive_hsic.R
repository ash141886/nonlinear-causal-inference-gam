# =============================================================================
# Core procedure: additive-HSIC causal discovery
# =============================================================================
#
# This file implements the non-linear causal-discovery routine used
# throughout the package. For every ordered pair of variables (i, j) the
# procedure:
#
#   (1) fits a penalised additive model of X_i on every variable except
#       X_i itself and the candidate X_j, using cubic B-splines and REML
#       smoothing (via mgcv::gam);
#   (2) evaluates the residual dependence between r_i^{(-j)} and X_j with
#       the biased empirical HSIC under a Gaussian kernel whose bandwidth
#       is set by the median heuristic;
#   (3) builds a pairwise test statistic T_{ij} = max(M_{ij}, M_{ji})
#       whose null distribution is approximated by permuting sample
#       indices.
#
# Edges are selected either by Benjamini-Hochberg FDR control at level
# alpha (single-stage mode) or by liberal BH screening followed by an
# effect-size filter that retains the top-k pairs by full-sample HSIC
# (two-stage mode). Orientation is determined by the asymmetry
# M_{ij} vs. M_{ji}, and the resulting graph is post-processed to be
# acyclic by iteratively removing the weakest edge of any remaining
# cycle.
# =============================================================================

source("R/utils.R")


# -----------------------------------------------------------------------------
# Leave-one-out additive-model fit
# -----------------------------------------------------------------------------

#' Fit a reduced additive model excluding one variable
#'
#' Fits \deqn{X_i = \sum_{k \neq i, j} f_k(X_k) + r_i^{(-j)}} by penalised
#' additive regression (cubic B-splines, REML smoothing, extra smoothing
#' penalty \eqn{\gamma = 1.4}, and term selection via \code{select = TRUE}).
#' The residual vector is returned and used downstream as an estimate of
#' the component of X_i that cannot be explained by \{X_k : k \neq i, j\}.
#'
#' @param data Data frame with n rows and p numeric columns.
#' @param i Integer index of the response variable.
#' @param j Integer index of the variable to be excluded from the predictors.
#' @param k_max Upper bound on the per-smooth basis dimension; defaults to
#'   \eqn{\min(\lfloor n/4 \rfloor, 40)} when \code{NULL}.
#' @return Numeric vector of length n containing the fitted residuals.
fit_reduced_model <- function(data, i, j, k_max = NULL) {
  n          <- nrow(data)
  p          <- ncol(data)
  var_names  <- colnames(data)

  if (is.null(k_max)) k_max <- min(floor(n / 4), 40)

  # All predictors except the response i and the excluded variable j.
  pred_idx <- setdiff(seq_len(p), c(i, j))

  if (length(pred_idx) == 0) {
    return(data[[i]] - mean(data[[i]]))
  }

  smooth_terms <- paste0("s(", var_names[pred_idx], ", k=", k_max, ")")
  formula_str  <- paste(var_names[i], "~", paste(smooth_terms, collapse = " + "))
  formula_obj  <- as.formula(formula_str)

  fit <- mgcv::gam(formula_obj, data = data, method = "REML",
                   gamma = 1.4, select = TRUE)

  as.numeric(residuals(fit))
}


# -----------------------------------------------------------------------------
# Main discovery routine
# -----------------------------------------------------------------------------

#' Additive-HSIC causal discovery
#'
#' Implements the procedure described at the top of this file. Two edge-
#' selection modes are supported:
#' \itemize{
#'   \item Single-stage: Benjamini-Hochberg control of the FDR at level
#'     \code{alpha} on the pairwise permutation p-values.
#'   \item Two-stage: liberal BH screening at \code{screening_alpha},
#'     followed by retention of the top-\eqn{k} pairs ranked by full-
#'     sample HSIC effect size, with
#'     \eqn{k = \lfloor p(p-1)\,\rho / 2 \rfloor} for a target edge
#'     density \eqn{\rho = } \code{target_density}.
#' }
#'
#' When \code{n_sub} is supplied, permutations are evaluated on a random
#' sub-sample of size \code{n_sub} while the GAM residuals continue to be
#' fitted on the full sample; this yields tractable kernel matrices for
#' large \eqn{n} without discarding regression precision.
#'
#' @param data Data frame with n rows and p numeric columns; will be
#'   internally standardised.
#' @param alpha FDR level (single-stage mode).
#' @param n_perm Number of permutations used to approximate the null.
#' @param k_max Upper bound on the per-smooth basis dimension.
#' @param two_stage Logical; if \code{TRUE}, use the two-stage selector.
#' @param screening_alpha Liberal BH level for the screening stage.
#' @param target_density Target edge density for the top-\eqn{k} filter.
#' @param n_sub Optional sub-sample size for the permutation kernels.
#' @param seed Random seed used to initialise the permutation loop.
#' @param verbose Logical; whether to print progress messages.
#' @return A list with components:
#'   \item{adjacency}{p-by-p binary adjacency matrix (the estimated DAG).}
#'   \item{scores}{p-by-p matrix of directional HSIC scores M_{ij}.}
#'   \item{pvalues}{Named numeric vector of permutation p-values, one per
#'     unordered pair.}
#'   \item{T_obs}{Named numeric vector of observed T_{ij} = max(M_{ij},
#'     M_{ji}).}
#'   \item{edges}{Data frame of selected edges with columns
#'     \code{parent}, \code{child}, \code{hsic_score}, \code{pvalue}.}
#'   \item{resid_cache}{List of cached leave-one-out residual vectors,
#'     keyed by "i_j".}
discover_dag <- function(data, alpha = 0.05, n_perm = 200, k_max = NULL,
                         two_stage = FALSE, screening_alpha = 0.20,
                         target_density = 0.15, n_sub = NULL,
                         seed = 42, verbose = TRUE) {
  set.seed(seed)
  n          <- nrow(data)
  p          <- ncol(data)
  var_names  <- colnames(data)

  if (is.null(k_max)) k_max <- min(floor(n / 4), 40)

  # Standardise and allocate the directional HSIC score matrix.
  data <- as.data.frame(scale(data))
  M    <- matrix(0, p, p, dimnames = list(var_names, var_names))

  # Leave-one-out residuals and directional HSIC scores for every ordered pair.
  if (verbose) cat("Computing leave-one-out residuals and HSIC scores...\n")
  resid_cache <- list()

  for (i in seq_len(p)) {
    if (verbose) cat(sprintf("  Variable %d/%d (%s)\n", i, p, var_names[i]))
    for (j in seq_len(p)) {
      if (i == j) next
      key                 <- paste(i, j, sep = "_")
      r_ij                <- fit_reduced_model(data, i, j, k_max)
      resid_cache[[key]]  <- r_ij
      M[i, j]             <- hsic_biased(r_ij, data[[j]])
    }
  }

  # Permutation p-values for the symmetric pair statistic T_{ij}.
  if (verbose) cat("Computing permutation p-values...\n")
  m_tests    <- p * (p - 1) / 2
  pairs      <- combn(p, 2)
  T_obs      <- numeric(ncol(pairs))
  pvalues    <- numeric(ncol(pairs))
  pair_names <- character(ncol(pairs))

  for (idx in seq_len(ncol(pairs))) {
    i <- pairs[1, idx]
    j <- pairs[2, idx]
    T_obs[idx]      <- max(M[i, j], M[j, i])
    pair_names[idx] <- paste(var_names[i], var_names[j], sep = "--")

    r_ij <- resid_cache[[paste(i, j, sep = "_")]]
    r_ji <- resid_cache[[paste(j, i, sep = "_")]]

    # Optionally evaluate permutations on a sub-sample to keep kernel
    # matrices of manageable size while retaining full-data residuals.
    if (!is.null(n_sub) && n_sub < n) {
      sub_idx   <- sample(n, n_sub)
      x_j_perm  <- data[[j]][sub_idx]
      x_i_perm  <- data[[i]][sub_idx]
      r_ij_perm <- r_ij[sub_idx]
      r_ji_perm <- r_ji[sub_idx]
    } else {
      x_j_perm  <- data[[j]]
      x_i_perm  <- data[[i]]
      r_ij_perm <- r_ij
      r_ji_perm <- r_ji
    }

    T_null <- numeric(n_perm)
    for (b in seq_len(n_perm)) {
      perm       <- sample(length(x_j_perm))
      h1         <- hsic_biased(r_ij_perm, x_j_perm[perm])
      h2         <- hsic_biased(r_ji_perm, x_i_perm[perm])
      T_null[b]  <- max(h1, h2)
    }

    pvalues[idx] <- (1 + sum(T_null >= T_obs[idx])) / (n_perm + 1)
  }

  names(pvalues) <- pair_names
  names(T_obs)   <- pair_names

  # Stage 1: Benjamini-Hochberg selection at the chosen level.
  if (two_stage) {
    effective_alpha <- screening_alpha
    k_target        <- floor(m_tests * target_density)
    if (verbose) cat(sprintf("Two-stage: BH at alpha=%.2f, then top-%d by HSIC\n",
                             effective_alpha, k_target))
  } else {
    effective_alpha <- alpha
  }

  ordered_idx <- order(pvalues)
  thresholds  <- seq_len(length(pvalues)) * effective_alpha / length(pvalues)
  rejected    <- pvalues[ordered_idx] <= thresholds
  k_star      <- max(c(0, which(rejected)))

  if (k_star > 0) {
    selected_idx <- ordered_idx[seq_len(k_star)]
  } else {
    selected_idx <- integer(0)
  }

  # Stage 2: effect-size filter (two-stage mode only).
  if (two_stage && length(selected_idx) > k_target) {
    hsic_ranks   <- order(T_obs[selected_idx], decreasing = TRUE)
    selected_idx <- selected_idx[hsic_ranks[seq_len(k_target)]]
  }

  # Orient each selected edge by the direction of larger HSIC asymmetry.
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
      adjacency[j, i] <- 1
      selected_edges <- rbind(selected_edges, data.frame(
        parent = var_names[j], child = var_names[i],
        hsic_score = M[i, j], pvalue = pvalues[idx],
        stringsAsFactors = FALSE
      ))
    } else {
      adjacency[i, j] <- 1
      selected_edges <- rbind(selected_edges, data.frame(
        parent = var_names[i], child = var_names[j],
        hsic_score = M[j, i], pvalue = pvalues[idx],
        stringsAsFactors = FALSE
      ))
    }
  }

  # Enforce acyclicity as a post-processing step.
  adjacency <- break_cycles(adjacency, M)

  if (verbose) {
    cat(sprintf("Edges selected: %d\n", sum(adjacency)))
    cat(sprintf("Pairs tested: %d, BH-rejected: %d\n", length(pvalues), k_star))
  }

  list(
    adjacency   = adjacency,
    scores      = M,
    pvalues     = pvalues,
    T_obs       = T_obs,
    edges       = selected_edges,
    resid_cache = resid_cache
  )
}


# -----------------------------------------------------------------------------
# Acyclicity post-processing
# -----------------------------------------------------------------------------

#' Remove directed cycles by deleting the weakest edge
#'
#' Repeatedly searches for a directed cycle and, when one is found, deletes
#' the edge in the cycle with the smallest HSIC score (using the reverse-
#' direction score as a proxy for edge strength in the selected orientation).
#' Terminates once the graph is acyclic or no edge can be safely removed.
#'
#' @param adj Adjacency matrix (will be modified).
#' @param M Matrix of directional HSIC scores (same shape as \code{adj}).
#' @return The updated, acyclic adjacency matrix.
break_cycles <- function(adj, M) {
  max_iter <- sum(adj)
  for (iter in seq_len(max_iter)) {
    cycle <- find_cycle(adj)
    if (is.null(cycle)) break

    min_score <- Inf
    min_edge  <- NULL
    for (k in seq_along(cycle)) {
      from <- cycle[k]
      to   <- cycle[if (k < length(cycle)) k + 1 else 1]
      if (adj[from, to] == 1 && M[to, from] < min_score) {
        min_score <- M[to, from]
        min_edge  <- c(from, to)
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

#' Locate a directed cycle in an adjacency matrix
#'
#' Depth-first search with a recursion stack to detect back-edges.
#'
#' @param adj Adjacency matrix.
#' @return An integer vector of node indices forming a cycle, or
#'   \code{NULL} if the graph is acyclic.
find_cycle <- function(adj) {
  p        <- nrow(adj)
  visited  <- rep(FALSE, p)
  in_stack <- rep(FALSE, p)
  parent   <- rep(0, p)

  for (start in seq_len(p)) {
    if (visited[start]) next
    result   <- dfs_cycle(adj, start, visited, in_stack, parent)
    if (!is.null(result$cycle)) return(result$cycle)
    visited  <- result$visited
    in_stack <- result$in_stack
  }
  NULL
}

#' Recursive DFS helper used by \code{find_cycle}
#'
#' @param adj Adjacency matrix.
#' @param node Current node index.
#' @param visited Logical vector: nodes already fully explored.
#' @param in_stack Logical vector: nodes currently on the recursion stack.
#' @param parent Integer vector used for cycle reconstruction.
#' @return A list with elements \code{cycle}, \code{visited}, \code{in_stack}.
dfs_cycle <- function(adj, node, visited, in_stack, parent) {
  visited[node]  <- TRUE
  in_stack[node] <- TRUE

  children <- which(adj[node, ] == 1)
  for (child in children) {
    if (!visited[child]) {
      parent[child] <- node
      result        <- dfs_cycle(adj, child, visited, in_stack, parent)
      if (!is.null(result$cycle)) return(result)
      visited  <- result$visited
      in_stack <- result$in_stack
    } else if (in_stack[child]) {
      # Reconstruct the cycle by walking the parent pointers back to child.
      cycle <- c(child)
      curr  <- node
      while (curr != child) {
        cycle <- c(curr, cycle)
        curr  <- parent[curr]
      }
      return(list(cycle = cycle, visited = visited, in_stack = in_stack))
    }
  }

  in_stack[node] <- FALSE
  list(cycle = NULL, visited = visited, in_stack = in_stack)
}
