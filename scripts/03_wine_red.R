# =============================================================================
# Causal Discovery on Wine Quality Data
# =============================================================================
#
# Application of the additive-HSIC method to the UCI Wine Quality dataset.
#
# Reference:
#   Islam, M.A. and Suzuki, J. "Non-linear Causal Inference in Observational
#   Data using Additive Models and Kernel-based Independence Testing"
#
# Data source:
#   Cortez et al. (2009). Modeling wine preferences by data mining from
#   physicochemical properties. Decision Support Systems, 47(4), 547-553.
#
# Repository: https://github.com/ash141886/nonlinear-causal-inference-gam
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
# Data Loading
# -----------------------------------------------------------------------------

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
  
  cat("Wine data loaded and standardized.\n")
  cat("Dimensions:", nrow(wine_data), "x", ncol(wine_data), "\n")
  cat("Variables:", paste(colnames(wine_data), collapse = ", "), "\n\n")
  
  wine_data
}

# -----------------------------------------------------------------------------
# Literature-Based Reference Graph
# -----------------------------------------------------------------------------

#' Construct reference DAG from oenological literature
#'
#' Based on established chemical relationships from:
#'   - Jackson (2020). Wine Science: Principles and Applications.
#'   - Waterhouse et al. (2024). Understanding Wine Chemistry.
#'
#' @param var_names Variable names from data
#' @return Adjacency matrix
get_reference_dag <- function(var_names) {
  n_vars <- length(var_names)
  dag <- matrix(0, n_vars, n_vars, dimnames = list(var_names, var_names))
  
  # Acid-pH relationships
  dag["fixed.acidity", "pH"] <- 1
  dag["volatile.acidity", "pH"] <- 1
  dag["citric.acid", "pH"] <- 1
  
  # Citric acid contributes to fixed acidity
  dag["citric.acid", "fixed.acidity"] <- 1
  
  # Density relationships
  dag["residual.sugar", "density"] <- 1
  dag["alcohol", "density"] <- 1
  
  # Sulfur dioxide relationship
  dag["free.sulfur.dioxide", "total.sulfur.dioxide"] <- 1
  
  # Quality determinants
  dag["alcohol", "quality"] <- 1
  dag["volatile.acidity", "quality"] <- 1
  dag["sulphates", "quality"] <- 1
  
  dag
}

# =============================================================================
# Kernel Functions
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
# HSIC Computation
# =============================================================================

#' Center a kernel matrix in O(n^2) instead of O(n^3)
#'
#' Replaces  H %*% K %*% H  (two dense matrix multiplies)
#' with the algebraic identity:
#'   (HKH)_{ab} = K_{ab} - mean_a(K) - mean_b(K) + mean(K)
#' which only requires row-means, column-means, and the grand mean.
center_kernel <- function(K) {
  n <- nrow(K)
  row_m <- rowMeans(K)               # O(n^2)
  grand_m <- mean(row_m)             # O(n)
  K - outer(row_m, rep(1, n)) -
      outer(rep(1, n), row_m) + grand_m
}

hsic_biased <- function(x, y, sigma_x = NULL, sigma_y = NULL) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  n <- nrow(x)

  if (is.null(sigma_x)) sigma_x <- median_bandwidth(x)
  if (is.null(sigma_y)) sigma_y <- median_bandwidth(y)

  K <- rbf_kernel(x, sigma_x)
  L <- rbf_kernel(y, sigma_y)
  Kc <- center_kernel(K)             # O(n^2) centering

  as.numeric(sum(Kc * L) / n^2)
}

#' Compute HSIC from precomputed centered kernel + raw kernel
#' Used inside the permutation loop to avoid recomputing kernels
hsic_from_kernels <- function(Kc, L, n) {
  as.numeric(sum(Kc * L) / n^2)
}

# =============================================================================
# Max-Statistic Permutation Test
# =============================================================================

#' Pairwise HSIC test using max-statistic permutation null
#'
#' Optimized: kernel matrices are computed ONCE per pair; each permutation
#' only shuffles rows/columns of the precomputed kernel (O(n^2) index
#' operation) instead of recomputing dist + exp + centering from scratch.
#'
#' @param resid_i Residuals from model for X_i excluding X_j
#' @param Xj Variable X_j
#' @param resid_j Residuals from model for X_j excluding X_i
#' @param Xi Variable X_i
#' @param B Number of permutations
#' @param seed Random seed
#' @return List with pvalue, scores, and test statistic
hsic_pairwise_test <- function(resid_i, Xj, resid_j, Xi, B = 200, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  resid_i <- as.matrix(resid_i)
  resid_j <- as.matrix(resid_j)
  Xi <- as.matrix(Xi)
  Xj <- as.matrix(Xj)
  n <- nrow(Xi)

  # Precompute all four kernel matrices ONCE
  sig_ri <- median_bandwidth(resid_i)
  sig_rj <- median_bandwidth(resid_j)
  sig_xj <- median_bandwidth(Xj)
  sig_xi <- median_bandwidth(Xi)

  Kc_ri <- center_kernel(rbf_kernel(resid_i, sig_ri))   # centered, fixed
  Kc_rj <- center_kernel(rbf_kernel(resid_j, sig_rj))   # centered, fixed
  L_xj  <- rbf_kernel(Xj, sig_xj)                       # raw (will be permuted)
  L_xi  <- rbf_kernel(Xi, sig_xi)                        # raw (will be permuted)

  # Observed scores  (O(n^2) each)
  score_ij <- hsic_from_kernels(Kc_ri, L_xj, n)
  score_ji <- hsic_from_kernels(Kc_rj, L_xi, n)
  t_obs <- max(score_ij, score_ji)

  # Permutation null: permute rows & cols of L  (O(n^2) per draw)
  t_perm <- numeric(B)
  for (b in seq_len(B)) {
    perm <- sample.int(n)
    L_xj_p <- L_xj[perm, perm]
    L_xi_p <- L_xi[perm, perm]
    s_ij_p <- hsic_from_kernels(Kc_ri, L_xj_p, n)
    s_ji_p <- hsic_from_kernels(Kc_rj, L_xi_p, n)
    t_perm[b] <- max(s_ij_p, s_ji_p)
  }

  # Monte Carlo p-value
  pvalue <- (sum(t_perm >= t_obs) + 1) / (B + 1)

  list(
    pvalue = pvalue,
    score_ij = score_ij,
    score_ji = score_ji,
    t_obs = t_obs
  )
}

# =============================================================================
# Leave-One-Out Residual Computation
# =============================================================================

#' Compute residuals from reduced additive model
#'
#' @param data Data frame
#' @param target Target variable index
#' @param exclude Variable index to exclude
#' @param var_names Variable names
#' @param gam_k_max Maximum basis dimension
#' @return Residual vector
compute_residuals_wine <- function(data, target, exclude, var_names, gam_k_max = 10) {
  n_vars <- ncol(data)
  n_samples <- nrow(data)
  predictors <- setdiff(seq_len(n_vars), c(target, exclude))
  
  target_name <- var_names[target]
  
  if (length(predictors) > 0) {
    k_basis <- min(gam_k_max, max(3, floor(n_samples / 40)))
    
    terms <- sapply(predictors, function(idx) {
      pred_name <- var_names[idx]
      n_unique <- length(unique(data[[pred_name]]))
      k_use <- min(k_basis, n_unique - 1)
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
        lm_terms <- paste0("`", var_names[predictors], "`", collapse = " + ")
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
# Causal Discovery Algorithm
# =============================================================================

#' Discover causal DAG using additive-HSIC method
#'
#' @param data Data frame
#' @param fdr_level FDR control level
#' @param n_perm Number of permutations
#' @param gam_k_max Maximum GAM basis dimension
#' @param perm_seed Random seed
#' @return List with adjacency matrix and score matrix
discover_dag_wine <- function(data, fdr_level = 0.05, n_perm = 200, 
                               gam_k_max = 10, perm_seed = NULL) {
  n_vars <- ncol(data)
  var_names <- colnames(data)
  
  cat("Computing leave-one-out residuals...\n")
  
  # Compute all leave-one-out residuals
  resid_cache <- list()
  for (i in seq_len(n_vars)) {
    for (j in seq_len(n_vars)) {
      if (i == j) next
      key <- paste(i, j, sep = "_")
      resid_cache[[key]] <- compute_residuals_wine(data, i, j, var_names, gam_k_max)
    }
  }
  
  cat("Testing pairwise independence...\n")
  
  # Test each unordered pair
  pair_results <- list()
  idx <- 0
  
  for (i in seq_len(n_vars - 1)) {
    for (j in (i + 1):n_vars) {
      idx <- idx + 1
      
      resid_i <- resid_cache[[paste(i, j, sep = "_")]]
      resid_j <- resid_cache[[paste(j, i, sep = "_")]]
      Xi <- data[[var_names[i]]]
      Xj <- data[[var_names[j]]]
      
      seed_pair <- if (!is.null(perm_seed)) perm_seed + i * 100 + j else NULL
      
      test <- hsic_pairwise_test(
        resid_i = resid_i, Xj = Xj,
        resid_j = resid_j, Xi = Xi,
        B = n_perm, seed = seed_pair
      )
      
      # Orient by argmax
      if (test$score_ij >= test$score_ji) {
        from_node <- j
        to_node <- i
        score <- test$score_ij
      } else {
        from_node <- i
        to_node <- j
        score <- test$score_ji
      }
      
      pair_results[[idx]] <- list(
        i = i, j = j,
        from = from_node, to = to_node,
        score = score, pvalue = test$pvalue,
        score_ij = test$score_ij, score_ji = test$score_ji
      )
    }
  }
  
  # BH-FDR control
  pvalues <- sapply(pair_results, `[[`, "pvalue")
  qvalues <- p.adjust(pvalues, method = "BH")
  
  # Build adjacency matrix
  adjacency <- matrix(0, n_vars, n_vars, dimnames = list(var_names, var_names))
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
  
  # Enforce acyclicity
  adjacency <- enforce_acyclicity_wine(adjacency, score_matrix)
  
  list(adjacency = adjacency, scores = score_matrix)
}

# =============================================================================
# Cycle Removal
# =============================================================================

enforce_acyclicity_wine <- function(adjacency, scores) {
  n <- nrow(adjacency)
  max_iter <- n * n
  iter <- 0
  
  while (iter < max_iter) {
    cycle <- find_cycle_wine(adjacency)
    if (is.null(cycle)) break
    
    min_score <- Inf
    weakest <- NULL
    
    for (k in seq_along(cycle)) {
      from <- cycle[k]
      to <- cycle[if (k == length(cycle)) 1 else k + 1]
      
      if (adjacency[from, to] == 1) {
        edge_score <- scores[to, from]
        if (is.na(edge_score)) edge_score <- 0
        
        if (edge_score < min_score) {
          min_score <- edge_score
          weakest <- c(from, to)
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

find_cycle_wine <- function(adjacency) {
  n <- nrow(adjacency)
  color <- rep("white", n)
  parent <- rep(NA_integer_, n)
  
  dfs <- function(v) {
    color[v] <<- "gray"
    children <- which(adjacency[v, ] == 1)
    
    for (u in children) {
      if (color[u] == "gray") {
        path <- u
        current <- v
        while (current != u && !is.na(parent[current])) {
          path <- c(path, current)
          current <- parent[current]
        }
        return(path)
      }
      if (color[u] == "white") {
        parent[u] <<- v
        result <- dfs(u)
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
# LiNGAM Baseline
# =============================================================================

run_lingam_wine <- function(data) {
  start_time <- Sys.time()
  var_names <- colnames(data)
  
  dag_result <- tryCatch({
    ica_result <- fastICA(as.matrix(data), n.comp = ncol(data),
                          alg.typ = "parallel", fun = "logcosh",
                          method = "C", verbose = FALSE)
    W <- ica_result$K %*% ica_result$W  # canonical: include pre-whitening
    B <- diag(ncol(data)) - W            # canonical LiNGAM: B = I - W_full
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
# Evaluation Metrics
# =============================================================================

#' Evaluate estimated DAG against reference
#'
#' @param estimated Estimated adjacency matrix
#' @param reference Reference adjacency matrix
#' @return Named vector of metrics
evaluate_dag_wine <- function(estimated, reference) {
  n <- nrow(reference)
  
  tp <- sum(estimated == 1 & reference == 1)
  fp <- sum(estimated == 1 & reference == 0)
  fn <- sum(estimated == 0 & reference == 1)
  tn <- sum(estimated == 0 & reference == 0) - n
  
  # Misoriented edges
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
  recall <- if ((tp + fn) > 0) tp / (tp + fn) else 0
  f1 <- if ((precision + recall) > 0) 2 * precision * recall / (precision + recall) else 0
  accuracy <- (tp + tn) / (n * (n - 1))
  shd <- sum(abs(estimated - reference))
  mse <- mean((estimated - reference)^2)
  
  c(F1 = f1, Accuracy = accuracy, SHD = shd, MSE = mse, 
    Misoriented = misoriented, Precision = precision, Recall = recall)
}

# =============================================================================
# Partial Effect Visualization
# =============================================================================

#' Extract and plot partial effects from GAM
#'
#' @param data Data frame
#' @param target Target variable name
#' @param predictor Predictor variable name
#' @param gam_k_max Maximum basis dimension
#' @return ggplot object
plot_partial_effect <- function(data, target, predictor, gam_k_max = 10) {
  var_names <- colnames(data)
  n_samples <- nrow(data)
  
  # Fit GAM with all predictors
  other_vars <- setdiff(var_names, target)
  k_basis <- min(gam_k_max, max(3, floor(n_samples / 40)))
  
  terms <- sapply(other_vars, function(v) {
    n_unique <- length(unique(data[[v]]))
    k_use <- min(k_basis, n_unique - 1)
    if (k_use >= 3) {
      paste0("s(`", v, "`, k=", k_use, ")")
    } else {
      paste0("`", v, "`")
    }
  })
  
  formula_str <- paste0("`", target, "` ~ ", paste(terms, collapse = " + "))
  model <- gam(as.formula(formula_str), data = data, method = "REML")
  
  # Extract smooth for predictor
  pred_seq <- seq(min(data[[predictor]]), max(data[[predictor]]), length.out = 100)
  newdata <- data.frame(matrix(0, nrow = 100, ncol = length(var_names)))
  colnames(newdata) <- var_names
  newdata[[predictor]] <- pred_seq
  
  # Get partial effect
  pred_result <- predict(model, newdata = newdata, type = "terms", se.fit = TRUE)
  
  # Find the column for the predictor
  term_names <- colnames(pred_result$fit)
  pred_col <- grep(predictor, term_names, value = TRUE)[1]
  
  if (!is.na(pred_col)) {
    partial_effect <- pred_result$fit[, pred_col]
    se <- pred_result$se.fit[, pred_col]
    
    # Get EDF
    smooth_idx <- grep(predictor, names(model$smooth))
    edf <- if (length(smooth_idx) > 0) round(model$edf[smooth_idx[1]], 1) else NA
    
    df <- data.frame(
      x = pred_seq,
      y = partial_effect,
      lower = partial_effect - 1.96 * se,
      upper = partial_effect + 1.96 * se
    )
    
    p <- ggplot(df, aes(x = x, y = y)) +
      geom_ribbon(aes(ymin = lower, ymax = upper), fill = "gray80", alpha = 0.5) +
      geom_line(linewidth = 1) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      labs(
        x = paste0(gsub("\\.", " ", predictor), " (std.)"),
        y = paste0(gsub("\\.", " ", target), " (partial)"),
        title = if (!is.na(edf)) paste0("EDF = ", edf) else NULL
      ) +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(size = 11, hjust = 0))
    
    return(p)
  }
  
  NULL
}

# =============================================================================
# Main Analysis
# =============================================================================

run_wine_analysis <- function() {
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("Wine Quality Causal Discovery Analysis\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
  
  # Load data
  wine_data <- load_wine_data()
  var_names <- colnames(wine_data)
  
  # Reference graph
  reference_dag <- get_reference_dag(var_names)
  cat("Reference graph constructed from literature.\n")
  cat("Number of reference edges:", sum(reference_dag), "\n\n")
  
  # Run proposed method
  cat("Running proposed method...\n")
  t0 <- Sys.time()
  proposed_result <- discover_dag_wine(
    wine_data,
    fdr_level = 0.05,
    n_perm = 200,
    gam_k_max = 10,
    perm_seed = 42
  )
  proposed_time <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("Completed in %.1f seconds.\n\n", proposed_time))
  
  # Run LiNGAM
  cat("Running LiNGAM...\n")
  lingam_result <- run_lingam_wine(wine_data)
  cat(sprintf("Completed in %.3f seconds.\n\n", lingam_result$time))
  
  # Evaluate
  proposed_metrics <- evaluate_dag_wine(proposed_result$adjacency, reference_dag)
  lingam_metrics <- evaluate_dag_wine(lingam_result$adjacency, reference_dag)
  
  # Results table
  results_df <- data.frame(
    Method = c("Proposed method", "LiNGAM"),
    F1 = c(proposed_metrics["F1"], lingam_metrics["F1"]),
    Accuracy = c(proposed_metrics["Accuracy"], lingam_metrics["Accuracy"]),
    SHD = c(proposed_metrics["SHD"], lingam_metrics["SHD"]),
    MSE = c(proposed_metrics["MSE"], lingam_metrics["MSE"]),
    Time = c(proposed_time, lingam_result$time)
  )
  
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("Performance Comparison\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  print(results_df, row.names = FALSE, digits = 3)
  
  # Generate partial effect plots
  cat("\nGenerating partial effect plots...\n")
  
  if (!dir.exists("figures")) dir.create("figures", recursive = TRUE)
  
  p1 <- plot_partial_effect(wine_data, "quality", "volatile.acidity")
  p2 <- plot_partial_effect(wine_data, "density", "alcohol")
  p3 <- plot_partial_effect(wine_data, "pH", "citric.acid")
  
  if (!is.null(p1) && !is.null(p2) && !is.null(p3)) {
    combined <- gridExtra::grid.arrange(p1, p2, p3, nrow = 1)
    ggsave("figures/wine_partial_effects.pdf", combined, width = 12, height = 4, dpi = 300)
    ggsave("figures/wine_partial_effects.png", combined, width = 12, height = 4, dpi = 300)
    cat("Partial effect plots saved to figures/\n")
  }
  
  cat("\nAnalysis complete.\n")
  
  list(
    proposed = proposed_result,
    lingam = lingam_result,
    reference = reference_dag,
    metrics = results_df
  )
}

# =============================================================================
# Execute
# =============================================================================

if (interactive() || !exists("SKIP_MAIN")) {
  results <- run_wine_analysis()
}
