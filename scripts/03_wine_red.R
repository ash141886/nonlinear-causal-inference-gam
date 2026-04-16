# =============================================================================
# Wine Quality (red) case study: single-stage discovery
# =============================================================================
#
# Applies the additive-HSIC causal-discovery procedure to the UCI Wine
# Quality (red) dataset of Cortez et al. (2009) and compares its output
# against a literature-informed reference DAG built from well-established
# oenological relationships. A fast implementation of DirectLiNGAM is
# included as a linear, non-Gaussian baseline.
#
# The script produces:
#   * a side-by-side table of directed-F1, accuracy, SHD, adjacency MSE,
#     count of mis-oriented edges, and wall-clock time for both methods;
#   * partial-effect plots for three representative GAM components that
#     are used as qualitative evidence of non-linearity.
#
# Data source:
#   Cortez, P., Cerdeira, A., Almeida, F., Matos, T., Reis, J. (2009).
#   Modeling wine preferences by data mining from physicochemical
#   properties. Decision Support Systems, 47(4), 547-553.
#
# Usage:  Rscript scripts/03_wine_red.R
# Output: figures/wine_partial_effects.{pdf,png}
# =============================================================================

# -----------------------------------------------------------------------------
# Dependencies
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(mgcv)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(fastICA)
  library(gridExtra)
  library(MASS)
})

# -----------------------------------------------------------------------------
# Data loading
# -----------------------------------------------------------------------------

#' Download (if necessary) and standardise the UCI red-wine table.
#'
#' The semicolon-separated CSV is cached in \code{data/} on first use.
#' Each column is rescaled to mean zero and unit variance.
#'
#' @return A standardised data frame with n = 1599 rows and p = 12 columns.
load_wine_data <- function() {
  if (!file.exists("data/winequality-red.csv")) {
    dir.create("data", showWarnings = FALSE, recursive = TRUE)
    download.file(
      "https://archive.ics.uci.edu/ml/machine-learning-databases/wine-quality/winequality-red.csv",
      "data/winequality-red.csv"
    )
  }
  wine_data <- read.csv("data/winequality-red.csv", sep = ";")
  wine_data <- as.data.frame(scale(wine_data))

  cat("Wine data loaded and standardised.\n")
  cat("Dimensions:", nrow(wine_data), "x", ncol(wine_data), "\n")
  cat("Variables:", paste(colnames(wine_data), collapse = ", "), "\n\n")

  wine_data
}

# -----------------------------------------------------------------------------
# Literature-based reference graph
# -----------------------------------------------------------------------------

#' Literature-informed reference DAG for the red-wine variables.
#'
#' Encodes chemical relationships that are well established in the
#' oenology literature (for example, Jackson (2020), Waterhouse et al.
#' (2024)): titratable acids drive pH; citric acid contributes to fixed
#' acidity; residual sugar and alcohol jointly determine density; free
#' sulphur dioxide is a component of total sulphur dioxide; and alcohol,
#' volatile acidity, and sulphates each influence the sensory quality
#' score.
#'
#' @param var_names Column names in the wine table.
#' @return A p-by-p binary adjacency matrix with named rows and columns.
get_reference_dag <- function(var_names) {
  n_vars <- length(var_names)
  dag <- matrix(0, n_vars, n_vars, dimnames = list(var_names, var_names))

  # Titratable acids -> pH.
  dag["fixed.acidity",    "pH"] <- 1
  dag["volatile.acidity", "pH"] <- 1
  dag["citric.acid",      "pH"] <- 1

  # Citric acid contributes to fixed acidity.
  dag["citric.acid", "fixed.acidity"] <- 1

  # Mass-balance relationships on density.
  dag["residual.sugar", "density"] <- 1
  dag["alcohol",        "density"] <- 1

  # Free sulphur dioxide is a component of total sulphur dioxide.
  dag["free.sulfur.dioxide", "total.sulfur.dioxide"] <- 1

  # Sensory quality determinants.
  dag["alcohol",         "quality"] <- 1
  dag["volatile.acidity","quality"] <- 1
  dag["sulphates",       "quality"] <- 1

  dag
}

# =============================================================================
# Gaussian-RBF kernel utilities
# =============================================================================

rbf_kernel <- function(x, sigma) {
  x <- as.matrix(x)
  dist_sq <- as.matrix(dist(x))^2
  exp(-dist_sq / (2 * sigma^2))
}

median_bandwidth <- function(x) {
  x <- as.matrix(x)
  dists <- as.matrix(dist(x))
  med <- median(dists[lower.tri(dists)])
  if (!is.finite(med) || med <= 0) med <- 1
  med
}

# =============================================================================
# HSIC computation
# =============================================================================

#' Centre a kernel matrix in O(n^2) memory-bound operations.
#'
#' Uses the algebraic identity
#'   (H K H)_{ab} = K_{ab} - mean_a(K) - mean_b(K) + mean(K),
#' avoiding the two dense O(n^3) matrix multiplications implied by
#' forming H K H directly.
center_kernel <- function(K) {
  n <- nrow(K)
  row_m   <- rowMeans(K)
  grand_m <- mean(row_m)
  K - outer(row_m, rep(1, n)) - outer(rep(1, n), row_m) + grand_m
}

#' Biased empirical HSIC with Gaussian-RBF kernels.
hsic_biased <- function(x, y, sigma_x = NULL, sigma_y = NULL) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  n <- nrow(x)

  if (is.null(sigma_x)) sigma_x <- median_bandwidth(x)
  if (is.null(sigma_y)) sigma_y <- median_bandwidth(y)

  K  <- rbf_kernel(x, sigma_x)
  L  <- rbf_kernel(y, sigma_y)
  Kc <- center_kernel(K)

  as.numeric(sum(Kc * L) / n^2)
}

#' HSIC from pre-computed kernel matrices (one already centred).
#'
#' Avoids recomputing distances and the centring projection inside tight
#' permutation loops.
hsic_from_kernels <- function(Kc, L, n) {
  as.numeric(sum(Kc * L) / n^2)
}

# =============================================================================
# Pairwise max-statistic permutation test
# =============================================================================

#' Permutation test for the symmetric pair statistic
#' \eqn{T_{ij} = \max(M_{ij}, M_{ji})}.
#'
#' Kernel matrices are computed once per pair; each permutation
#' corresponds to a symmetric row/column reordering of the raw kernel,
#' which is an O(n^2) indexing operation and avoids repeated kernel
#' construction.
#'
#' @param resid_i Leave-one-out residuals for response X_i excluding X_j.
#' @param Xj Vector of values of X_j.
#' @param resid_j Leave-one-out residuals for response X_j excluding X_i.
#' @param Xi Vector of values of X_i.
#' @param B Number of permutations.
#' @param seed Optional random seed.
#' @return A list with the Monte Carlo p-value, the two directional
#'   HSIC scores, and the observed pair statistic.
hsic_pairwise_test <- function(resid_i, Xj, resid_j, Xi, B = 200, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  resid_i <- as.matrix(resid_i)
  resid_j <- as.matrix(resid_j)
  Xi <- as.matrix(Xi)
  Xj <- as.matrix(Xj)
  n  <- nrow(Xi)

  # Pre-compute and cache the four relevant kernel matrices.
  sig_ri <- median_bandwidth(resid_i)
  sig_rj <- median_bandwidth(resid_j)
  sig_xj <- median_bandwidth(Xj)
  sig_xi <- median_bandwidth(Xi)

  Kc_ri <- center_kernel(rbf_kernel(resid_i, sig_ri))  # centred, fixed
  Kc_rj <- center_kernel(rbf_kernel(resid_j, sig_rj))  # centred, fixed
  L_xj  <- rbf_kernel(Xj, sig_xj)                      # raw, will be permuted
  L_xi  <- rbf_kernel(Xi, sig_xi)                      # raw, will be permuted

  score_ij <- hsic_from_kernels(Kc_ri, L_xj, n)
  score_ji <- hsic_from_kernels(Kc_rj, L_xi, n)
  t_obs    <- max(score_ij, score_ji)

  # Null distribution via symmetric row/column permutations of L.
  t_perm <- numeric(B)
  for (b in seq_len(B)) {
    perm     <- sample.int(n)
    L_xj_p   <- L_xj[perm, perm]
    L_xi_p   <- L_xi[perm, perm]
    s_ij_p   <- hsic_from_kernels(Kc_ri, L_xj_p, n)
    s_ji_p   <- hsic_from_kernels(Kc_rj, L_xi_p, n)
    t_perm[b] <- max(s_ij_p, s_ji_p)
  }

  pvalue <- (sum(t_perm >= t_obs) + 1) / (B + 1)

  list(
    pvalue   = pvalue,
    score_ij = score_ij,
    score_ji = score_ji,
    t_obs    = t_obs
  )
}

# =============================================================================
# Leave-one-out residuals
# =============================================================================

#' Residuals from a reduced additive model excluding one variable.
#'
#' Fits X_{target} on the remaining variables except X_{exclude} via a
#' penalised additive model with cubic B-splines and REML smoothing. The
#' per-smooth basis dimension is set adaptively to
#' \eqn{\min(k_{max}, \lfloor n/40 \rfloor)} and is further capped by the
#' number of unique values of each predictor; predictors with fewer than
#' three unique values are entered linearly. If the GAM fit fails, the
#' routine falls back to ordinary least squares.
compute_residuals_wine <- function(data, target, exclude, var_names, gam_k_max = 10) {
  n_vars    <- ncol(data)
  n_samples <- nrow(data)
  predictors <- setdiff(seq_len(n_vars), c(target, exclude))

  target_name <- var_names[target]

  if (length(predictors) > 0) {
    k_basis <- min(gam_k_max, max(3, floor(n_samples / 40)))

    terms <- sapply(predictors, function(idx) {
      pred_name <- var_names[idx]
      n_unique  <- length(unique(data[[pred_name]]))
      k_use     <- min(k_basis, n_unique - 1)
      if (k_use >= 3) {
        paste0("s(`", pred_name, "`, k=", k_use, ")")
      } else {
        paste0("`", pred_name, "`")
      }
    })

    formula_str <- paste0("`", target_name, "` ~ ", paste(terms, collapse = " + "))

    model <- tryCatch(
      suppressWarnings(
        gam(as.formula(formula_str), data = data, method = "REML", gamma = 1.4)
      ),
      error = function(e) {
        lm_terms   <- paste0("`", var_names[predictors], "`", collapse = " + ")
        lm_formula <- paste0("`", target_name, "` ~ ", lm_terms)
        lm(as.formula(lm_formula), data = data)
      }
    )
    residuals(model)
  } else {
    data[[target_name]] - mean(data[[target_name]])
  }
}

# =============================================================================
# Causal-discovery algorithm
# =============================================================================

#' Additive-HSIC causal discovery with BH-FDR selection.
#'
#' Steps:
#'   1. Cache leave-one-out residuals for every ordered pair (i, j).
#'   2. For every unordered pair, compute the directional HSIC scores and
#'      obtain a Monte Carlo p-value via the max-statistic permutation
#'      test.
#'   3. Apply Benjamini-Hochberg FDR control at level \code{fdr_level}.
#'   4. Orient each selected edge by argmax of the directional scores.
#'   5. Remove any directed cycles by iteratively deleting the weakest
#'      edge on the cycle.
#'
#' @param data Standardised data frame.
#' @param fdr_level FDR level.
#' @param n_perm Number of permutations per pair.
#' @param gam_k_max Upper bound on the per-smooth basis dimension.
#' @param perm_seed Optional base seed used to seed each pair-wise test
#'   reproducibly.
#' @return A list with the estimated adjacency matrix and the matrix of
#'   directional HSIC scores.
discover_dag_wine <- function(data, fdr_level = 0.05, n_perm = 200,
                              gam_k_max = 10, perm_seed = NULL) {
  n_vars    <- ncol(data)
  var_names <- colnames(data)

  cat("Computing leave-one-out residuals...\n")
  resid_cache <- list()
  for (i in seq_len(n_vars)) {
    for (j in seq_len(n_vars)) {
      if (i == j) next
      key <- paste(i, j, sep = "_")
      resid_cache[[key]] <- compute_residuals_wine(data, i, j, var_names, gam_k_max)
    }
  }

  cat("Testing pair-wise independence...\n")
  pair_results <- list()
  idx <- 0

  for (i in seq_len(n_vars - 1)) {
    for (j in (i + 1):n_vars) {
      idx <- idx + 1

      resid_i <- resid_cache[[paste(i, j, sep = "_")]]
      resid_j <- resid_cache[[paste(j, i, sep = "_")]]
      Xi      <- data[[var_names[i]]]
      Xj      <- data[[var_names[j]]]

      seed_pair <- if (!is.null(perm_seed)) perm_seed + i * 100 + j else NULL

      test <- hsic_pairwise_test(
        resid_i = resid_i, Xj = Xj,
        resid_j = resid_j, Xi = Xi,
        B = n_perm, seed = seed_pair
      )

      # Orient the edge by the direction of larger HSIC asymmetry.
      if (test$score_ij >= test$score_ji) {
        from_node <- j
        to_node   <- i
        score     <- test$score_ij
      } else {
        from_node <- i
        to_node   <- j
        score     <- test$score_ji
      }

      pair_results[[idx]] <- list(
        i = i, j = j,
        from = from_node, to = to_node,
        score = score, pvalue = test$pvalue,
        score_ij = test$score_ij, score_ji = test$score_ji
      )
    }
  }

  # Benjamini-Hochberg FDR selection on the permutation p-values.
  pvalues <- sapply(pair_results, `[[`, "pvalue")
  qvalues <- p.adjust(pvalues, method = "BH")

  adjacency    <- matrix(0, n_vars, n_vars, dimnames = list(var_names, var_names))
  score_matrix <- matrix(0, n_vars, n_vars, dimnames = list(var_names, var_names))

  for (k in seq_along(pair_results)) {
    pr <- pair_results[[k]]
    score_matrix[var_names[pr$j], var_names[pr$i]] <- pr$score_ij
    score_matrix[var_names[pr$i], var_names[pr$j]] <- pr$score_ji

    if (qvalues[k] < fdr_level) {
      adjacency[var_names[pr$from], var_names[pr$to]] <- 1
    }
  }

  cat("Enforcing acyclicity...\n")
  adjacency <- enforce_acyclicity_wine(adjacency, score_matrix)

  list(adjacency = adjacency, scores = score_matrix)
}

# =============================================================================
# Acyclicity post-processing
# =============================================================================

#' Repeatedly delete the weakest edge of any directed cycle until the
#' graph is acyclic. Edge strength is measured by the reverse-direction
#' HSIC score.
enforce_acyclicity_wine <- function(adjacency, scores) {
  n        <- nrow(adjacency)
  max_iter <- n * n
  iter     <- 0

  while (iter < max_iter) {
    cycle <- find_cycle_wine(adjacency)
    if (is.null(cycle)) break

    min_score <- Inf
    weakest   <- NULL

    for (k in seq_along(cycle)) {
      from <- cycle[k]
      to   <- cycle[if (k == length(cycle)) 1 else k + 1]

      if (adjacency[from, to] == 1) {
        edge_score <- scores[to, from]
        if (is.na(edge_score)) edge_score <- 0

        if (edge_score < min_score) {
          min_score <- edge_score
          weakest   <- c(from, to)
        }
      }
    }

    if (!is.null(weakest)) {
      adjacency[weakest[1], weakest[2]] <- 0
    }
    iter <- iter + 1
  }

  adjacency
}

#' Locate a directed cycle in an adjacency matrix via DFS.
find_cycle_wine <- function(adjacency) {
  n      <- nrow(adjacency)
  color  <- rep("white", n)
  parent <- rep(NA_integer_, n)

  dfs <- function(v) {
    color[v] <<- "gray"
    children <- which(adjacency[v, ] == 1)

    for (u in children) {
      if (color[u] == "gray") {
        path    <- u
        current <- v
        while (current != u && !is.na(parent[current])) {
          path    <- c(path, current)
          current <- parent[current]
        }
        return(path)
      }
      if (color[u] == "white") {
        parent[u] <<- v
        result     <- dfs(u)
        if (!is.null(result)) return(result)
      }
    }
    color[v] <<- "black"
    NULL
  }

  for (v in seq_len(n)) {
    if (color[v] == "white") {
      result <- dfs(v)
      if (!is.null(result)) return(result)
    }
  }
  NULL
}

# =============================================================================
# LiNGAM baseline (ICA-based)
# =============================================================================

#' ICA-based LiNGAM estimator with a hard threshold at |B_{ij}| = 0.1.
run_lingam_wine <- function(data) {
  start_time <- Sys.time()
  var_names  <- colnames(data)

  dag_result <- tryCatch({
    ica_result <- fastICA(as.matrix(data), n.comp = ncol(data),
                          alg.typ = "parallel", fun = "logcosh",
                          method = "C", verbose = FALSE)
    # Include the pre-whitening factor so that W corresponds to the full
    # demixing of the standardised data.
    W <- ica_result$K %*% ica_result$W
    # Canonical LiNGAM parameterisation X = B X + e, so B = I - W.
    B <- diag(ncol(data)) - W
    diag(B) <- 0
    B[abs(B) < 0.1] <- 0
    dag <- (abs(B) > 0) * 1
    dag
  }, error = function(e) {
    matrix(0, ncol(data), ncol(data))
  })

  time_taken <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  dimnames(dag_result) <- list(var_names, var_names)

  list(adjacency = dag_result, time = time_taken)
}

# =============================================================================
# Evaluation metrics
# =============================================================================

#' Compare an estimated DAG against a reference DAG.
#'
#' Reports the standard directed-edge metrics together with a scalar
#' structural Hamming distance and an adjacency mean-squared error. The
#' \code{Misoriented} entry counts edges that are present in the
#' estimate but reversed relative to the reference and is already
#' accounted for inside \code{SHD}.
evaluate_dag_wine <- function(estimated, reference) {
  n <- nrow(reference)

  tp <- sum(estimated == 1 & reference == 1)
  fp <- sum(estimated == 1 & reference == 0)
  fn <- sum(estimated == 0 & reference == 1)
  tn <- sum(estimated == 0 & reference == 0) - n  # subtract diagonal matches

  misoriented <- 0
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i < j) {
        if (reference[i, j] == 1 && estimated[j, i] == 1 && estimated[i, j] == 0) {
          misoriented <- misoriented + 1
        } else if (reference[j, i] == 1 && estimated[i, j] == 1 && estimated[j, i] == 0) {
          misoriented <- misoriented + 1
        }
      }
    }
  }

  precision <- if ((tp + fp) > 0) tp / (tp + fp) else 0
  recall    <- if ((tp + fn) > 0) tp / (tp + fn) else 0
  f1        <- if ((precision + recall) > 0) 2 * precision * recall / (precision + recall) else 0
  accuracy  <- (tp + tn) / (n * (n - 1))
  shd       <- sum(abs(estimated - reference))
  mse       <- mean((estimated - reference)^2)

  c(F1 = f1, Accuracy = accuracy, SHD = shd, MSE = mse,
    Misoriented = misoriented, Precision = precision, Recall = recall)
}

# =============================================================================
# Partial-effect visualisation
# =============================================================================

#' Extract and plot the partial effect of one predictor from a GAM fit of
#' the target on all other variables.
#'
#' The per-smooth basis dimension is set adaptively; predictors with too
#' few unique values are entered linearly. Pointwise 95% confidence bands
#' are drawn from the GAM standard errors, and the effective degrees of
#' freedom of the smooth are reported as the plot title.
plot_partial_effect <- function(data, target, predictor, gam_k_max = 10) {
  var_names <- colnames(data)
  n_samples <- nrow(data)

  other_vars <- setdiff(var_names, target)
  k_basis    <- min(gam_k_max, max(3, floor(n_samples / 40)))

  terms <- sapply(other_vars, function(v) {
    n_unique <- length(unique(data[[v]]))
    k_use    <- min(k_basis, n_unique - 1)
    if (k_use >= 3) {
      paste0("s(`", v, "`, k=", k_use, ")")
    } else {
      paste0("`", v, "`")
    }
  })

  formula_str <- paste0("`", target, "` ~ ", paste(terms, collapse = " + "))
  model       <- gam(as.formula(formula_str), data = data, method = "REML")

  # Evaluate the predictor's smooth on an equispaced grid.
  pred_seq <- seq(min(data[[predictor]]), max(data[[predictor]]), length.out = 100)
  newdata  <- data.frame(matrix(0, nrow = 100, ncol = length(var_names)))
  colnames(newdata) <- var_names
  newdata[[predictor]] <- pred_seq

  pred_result <- predict(model, newdata = newdata, type = "terms", se.fit = TRUE)

  term_names <- colnames(pred_result$fit)
  pred_col   <- grep(predictor, term_names, value = TRUE)[1]

  if (!is.na(pred_col)) {
    partial_effect <- pred_result$fit[, pred_col]
    se             <- pred_result$se.fit[, pred_col]

    smooth_idx <- grep(predictor, names(model$smooth))
    edf        <- if (length(smooth_idx) > 0) round(model$edf[smooth_idx[1]], 1) else NA

    df <- data.frame(
      x     = pred_seq,
      y     = partial_effect,
      lower = partial_effect - 1.96 * se,
      upper = partial_effect + 1.96 * se
    )

    p <- ggplot(df, aes(x = x, y = y)) +
      geom_ribbon(aes(ymin = lower, ymax = upper), fill = "gray80", alpha = 0.5) +
      geom_line(linewidth = 1) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      labs(
        x     = paste0(gsub("\\.", " ", predictor), " (std.)"),
        y     = paste0(gsub("\\.", " ", target),    " (partial)"),
        title = if (!is.na(edf)) paste0("EDF = ", edf) else NULL
      ) +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(size = 11, hjust = 0))

    return(p)
  }

  NULL
}

# =============================================================================
# Main analysis
# =============================================================================

#' End-to-end red-wine analysis.
#'
#' Loads the data, runs both methods, evaluates each against the
#' literature-informed reference, and writes three partial-effect plots
#' to \code{figures/}.
run_wine_analysis <- function() {
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("Red-wine causal-discovery analysis\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")

  wine_data <- load_wine_data()
  var_names <- colnames(wine_data)

  reference_dag <- get_reference_dag(var_names)
  cat("Reference graph constructed from the oenology literature.\n")
  cat("Number of reference edges:", sum(reference_dag), "\n\n")

  cat("Running the proposed additive-HSIC method...\n")
  t0 <- Sys.time()
  proposed_result <- discover_dag_wine(
    wine_data,
    fdr_level  = 0.05,
    n_perm     = 200,
    gam_k_max  = 10,
    perm_seed  = 42
  )
  proposed_time <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("Completed in %.1f seconds.\n\n", proposed_time))

  cat("Running the LiNGAM baseline...\n")
  lingam_result <- run_lingam_wine(wine_data)
  cat(sprintf("Completed in %.3f seconds.\n\n", lingam_result$time))

  proposed_metrics <- evaluate_dag_wine(proposed_result$adjacency, reference_dag)
  lingam_metrics   <- evaluate_dag_wine(lingam_result$adjacency,   reference_dag)

  results_df <- data.frame(
    Method   = c("Proposed method", "LiNGAM"),
    F1       = c(proposed_metrics["F1"],       lingam_metrics["F1"]),
    Accuracy = c(proposed_metrics["Accuracy"], lingam_metrics["Accuracy"]),
    SHD      = c(proposed_metrics["SHD"],      lingam_metrics["SHD"]),
    MSE      = c(proposed_metrics["MSE"],      lingam_metrics["MSE"]),
    Time     = c(proposed_time,                lingam_result$time)
  )

  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("Performance comparison\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  print(results_df, row.names = FALSE, digits = 3)

  cat("\nGenerating partial-effect plots...\n")

  if (!dir.exists("figures")) dir.create("figures", recursive = TRUE)

  p1 <- plot_partial_effect(wine_data, "quality", "volatile.acidity")
  p2 <- plot_partial_effect(wine_data, "density", "alcohol")
  p3 <- plot_partial_effect(wine_data, "pH",      "citric.acid")

  if (!is.null(p1) && !is.null(p2) && !is.null(p3)) {
    combined <- gridExtra::grid.arrange(p1, p2, p3, nrow = 1)
    ggsave("figures/wine_partial_effects.pdf", combined, width = 12, height = 4, dpi = 300)
    ggsave("figures/wine_partial_effects.png", combined, width = 12, height = 4, dpi = 300)
    cat("Partial-effect plots written to figures/.\n")
  }

  cat("\nAnalysis complete.\n")

  list(
    proposed  = proposed_result,
    lingam    = lingam_result,
    reference = reference_dag,
    metrics   = results_df
  )
}

# =============================================================================
# Entry point
# =============================================================================

if (interactive() || !exists("SKIP_MAIN")) {
  results <- run_wine_analysis()
}
