# =============================================================================
# Wine Quality (red) case study: two-stage hybrid selector
# =============================================================================
#
# Runs the additive-HSIC procedure on the UCI red-wine data under a
# two-stage edge-selection rule:
#
#   Stage 1.  Fit leave-one-out GAMs on the full sample (n = 1599).
#   Stage 2.  Compute permutation p-values for the max-statistic pair
#             test on a random sub-sample (n_sub = 200, B = 1000
#             permutations) to keep each kernel matrix tractable.
#   Stage 3.  Liberal BH-FDR screening at alpha = 0.20, followed by
#             retention of the top k = 10 pairs by full-sample HSIC
#             effect size. Directions are assigned from the
#             full-sample asymmetry M_{ij} vs. M_{ji}, and the result
#             is made acyclic by iteratively deleting the weakest edge
#             of any remaining cycle.
#
# The script prints a summary performance block and LaTeX-ready rows
# for the two tables that summarise (i) overall accuracy versus the
# baseline and (ii) the ranked list of retained edges and their
# overlap with the literature-informed reference graph.
#
# Usage:  Rscript scripts/03b_wine_red_twostage.R
# =============================================================================

suppressPackageStartupMessages({ library(mgcv); library(MASS) })

# -----------------------------------------------------------------------------
# Data loading and reference graph
# -----------------------------------------------------------------------------

load_wine_data <- function() {
  if (!file.exists("data/winequality-red.csv")) {
    dir.create("data", showWarnings = FALSE, recursive = TRUE)
    download.file("https://archive.ics.uci.edu/ml/machine-learning-databases/wine-quality/winequality-red.csv",
                  "data/winequality-red.csv")
  }
  as.data.frame(scale(read.csv("data/winequality-red.csv", sep = ";")))
}

#' Literature-informed reference DAG for the red-wine variables.
#' Encodes titratable-acid -> pH, citric-acid -> fixed-acidity, sugar
#' and alcohol -> density, free SO2 -> total SO2, and the standard set
#' of sensory-quality determinants.
get_reference_dag <- function(var_names) {
  n <- length(var_names)
  dag <- matrix(0, n, n, dimnames = list(var_names, var_names))
  dag["fixed.acidity",    "pH"] <- 1
  dag["volatile.acidity", "pH"] <- 1
  dag["citric.acid",      "pH"] <- 1
  dag["citric.acid",      "fixed.acidity"] <- 1
  dag["residual.sugar",   "density"] <- 1
  dag["alcohol",          "density"] <- 1
  dag["free.sulfur.dioxide", "total.sulfur.dioxide"] <- 1
  dag["alcohol",          "quality"] <- 1
  dag["volatile.acidity", "quality"] <- 1
  dag["sulphates",        "quality"] <- 1
  dag
}

# -----------------------------------------------------------------------------
# Kernel and HSIC utilities
# -----------------------------------------------------------------------------

rbf_kernel <- function(x, sigma) {
  x <- as.matrix(x)
  exp(-as.matrix(dist(x))^2 / (2 * sigma^2))
}

median_bandwidth <- function(x) {
  x <- as.matrix(x)
  d <- as.matrix(dist(x))
  med <- median(d[lower.tri(d)])
  if (!is.finite(med) || med <= 0) med <- 1
  med
}

#' Double-centring of a kernel matrix via the O(n^2) algebraic identity
#' (H K H)_{ab} = K_{ab} - mean_a(K) - mean_b(K) + mean(K).
center_kernel <- function(K) {
  rm <- rowMeans(K)
  K - outer(rm, rep(1, nrow(K))) - outer(rep(1, nrow(K)), rm) + mean(rm)
}

hsic_from_kernels <- function(Kc, L, n) as.numeric(sum(Kc * L) / n^2)

#' Biased empirical HSIC with Gaussian-RBF kernels and median-heuristic
#' bandwidths.
hsic_biased <- function(x, y) {
  x <- as.matrix(x); y <- as.matrix(y); n <- nrow(x)
  K <- rbf_kernel(x, median_bandwidth(x))
  L <- rbf_kernel(y, median_bandwidth(y))
  hsic_from_kernels(center_kernel(K), L, n)
}

# -----------------------------------------------------------------------------
# Leave-one-out residuals
# -----------------------------------------------------------------------------

#' Fit a reduced additive model for \code{target} excluding \code{exclude}
#' and return the residuals. Falls back to ordinary least squares if the
#' GAM fit fails; predictors with fewer than three unique values are
#' entered linearly.
compute_residuals_wine <- function(data, target, exclude, var_names, gam_k_max = 10) {
  n_vars    <- ncol(data)
  n_samples <- nrow(data)
  predictors <- setdiff(seq_len(n_vars), c(target, exclude))
  target_name <- var_names[target]

  if (length(predictors) > 0) {
    k_basis <- min(gam_k_max, max(3, floor(n_samples / 40)))
    terms <- sapply(predictors, function(idx) {
      pn <- var_names[idx]
      nu <- length(unique(data[[pn]]))
      ku <- min(k_basis, nu - 1)
      if (ku >= 3) paste0("s(`", pn, "`, k=", ku, ")") else paste0("`", pn, "`")
    })
    fs <- paste0("`", target_name, "` ~ ", paste(terms, collapse = " + "))
    model <- tryCatch(
      suppressWarnings(gam(as.formula(fs), data = data, method = "REML", gamma = 1.4)),
      error = function(e) {
        lt <- paste0("`", var_names[predictors], "`", collapse = " + ")
        lm(as.formula(paste0("`", target_name, "` ~ ", lt)), data = data)
      }
    )
    residuals(model)
  } else {
    data[[target_name]] - mean(data[[target_name]])
  }
}

# -----------------------------------------------------------------------------
# Acyclicity post-processing
# -----------------------------------------------------------------------------

#' Remove directed cycles by repeatedly deleting the edge with the
#' smallest reverse-direction HSIC score on the detected cycle.
enforce_acyclicity_wine <- function(adj, scores) {
  n <- nrow(adj); iter <- 0
  while (iter < n * n) {
    cyc <- find_cycle_wine(adj); if (is.null(cyc)) break
    ms <- Inf; w <- NULL
    for (k in seq_along(cyc)) {
      fr <- cyc[k]
      to <- cyc[if (k == length(cyc)) 1 else k + 1]
      if (adj[fr, to] == 1) {
        sc <- scores[to, fr]; if (is.na(sc)) sc <- 0
        if (sc < ms) { ms <- sc; w <- c(fr, to) }
      }
    }
    if (!is.null(w)) adj[w[1], w[2]] <- 0
    iter <- iter + 1
  }
  adj
}

#' Locate a directed cycle via depth-first search with colouring.
find_cycle_wine <- function(adj) {
  n <- nrow(adj); col <- rep("w", n); par <- rep(NA_integer_, n)
  dfs <- function(v) {
    col[v] <<- "g"
    for (u in which(adj[v, ] == 1)) {
      if (col[u] == "g") {
        p <- u; c <- v
        while (c != u && !is.na(par[c])) { p <- c(p, c); c <- par[c] }
        return(p)
      }
      if (col[u] == "w") {
        par[u] <<- v
        r <- dfs(u); if (!is.null(r)) return(r)
      }
    }
    col[v] <<- "b"; NULL
  }
  for (v in 1:n) if (col[v] == "w") { r <- dfs(v); if (!is.null(r)) return(r) }
  NULL
}

# -----------------------------------------------------------------------------
# Evaluation
# -----------------------------------------------------------------------------

evaluate_dag_wine <- function(est, ref) {
  n <- nrow(ref)
  tp <- sum(est == 1 & ref == 1)
  fp <- sum(est == 1 & ref == 0)
  fn <- sum(est == 0 & ref == 1)
  tn <- sum(est == 0 & ref == 0) - n   # subtract diagonal matches
  prec <- if ((tp + fp) > 0) tp / (tp + fp) else 0
  rec  <- if ((tp + fn) > 0) tp / (tp + fn) else 0
  f1   <- if ((prec + rec) > 0) 2 * prec * rec / (prec + rec) else 0
  c(F1 = f1, Accuracy = (tp + tn) / (n * (n - 1)),
    SHD = sum(abs(est - ref)), MSE = mean((est - ref)^2),
    TP = tp, FP = fp, FN = fn, Precision = prec, Recall = rec,
    Edges = tp + fp)
}

#' Pretty LaTeX label for each wine variable.
latex_name <- function(v) {
  map <- c(
    "fixed.acidity"        = "Fixed acidity",
    "volatile.acidity"     = "Volatile acidity",
    "citric.acid"          = "Citric acid",
    "residual.sugar"       = "Residual sugar",
    "chlorides"            = "Chlorides",
    "free.sulfur.dioxide"  = "Free SO$_2$",
    "total.sulfur.dioxide" = "Total SO$_2$",
    "density"              = "Density",
    "pH"                   = "pH",
    "sulphates"            = "Sulphates",
    "alcohol"              = "Alcohol",
    "quality"              = "Quality"
  )
  if (v %in% names(map)) map[v] else v
}

# =============================================================================
# Main analysis
# =============================================================================

cat(paste(rep("=", 70), collapse = ""), "\n")
cat("Two-stage hybrid selector (BH-FDR screening + HSIC top-k)\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

wine_data <- load_wine_data()
var_names <- colnames(wine_data)
n_vars    <- ncol(wine_data)
ref_dag   <- get_reference_dag(var_names)
cat("Reference edges:", sum(ref_dag), "\n\n")

# -----------------------------------------------------------------------------
# Stage 1: leave-one-out GAM fits on the full sample
# -----------------------------------------------------------------------------

cat("Stage 1: fitting GAMs on the full sample (n = 1599)...\n")
t0 <- Sys.time()
rc <- list()
for (i in 1:n_vars) {
  cat(sprintf("  Variable %d/%d (%s)\n", i, n_vars, var_names[i]))
  for (j in 1:n_vars) {
    if (i == j) next
    rc[[paste(i, j, sep = "_")]] <-
      compute_residuals_wine(wine_data, i, j, var_names, 10)
  }
}
gam_time <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("  Done in %.0fs.\n\n", gam_time))

# -----------------------------------------------------------------------------
# Stage 2: permutation p-values on a sub-sample
# -----------------------------------------------------------------------------

N_SUB        <- 200
B_PERM       <- 1000
SCREEN_ALPHA <- 0.20
K_TARGET     <- 10

cat(sprintf("Stage 2: HSIC + permutations (n_sub = %d, B = %d)...\n", N_SUB, B_PERM))
t1 <- Sys.time()
set.seed(42)
sub_idx <- sample.int(nrow(wine_data), N_SUB)

pair_results <- list(); idx <- 0
for (i in 1:(n_vars - 1)) {
  for (j in (i + 1):n_vars) {
    idx <- idx + 1
    ri <- rc[[paste(i, j, sep = "_")]]; rj <- rc[[paste(j, i, sep = "_")]]
    Xi <- wine_data[[var_names[i]]];    Xj <- wine_data[[var_names[j]]]

    ri_s <- as.matrix(ri[sub_idx]); rj_s <- as.matrix(rj[sub_idx])
    xi_s <- as.matrix(Xi[sub_idx]); xj_s <- as.matrix(Xj[sub_idx])
    n <- N_SUB

    # Use a reproducible per-pair seed for the permutation loop.
    set.seed(42 + i * 100 + j)
    Kc_ri <- center_kernel(rbf_kernel(ri_s, median_bandwidth(ri_s)))
    Kc_rj <- center_kernel(rbf_kernel(rj_s, median_bandwidth(rj_s)))
    L_xj  <- rbf_kernel(xj_s, median_bandwidth(xj_s))
    L_xi  <- rbf_kernel(xi_s, median_bandwidth(xi_s))

    sij  <- hsic_from_kernels(Kc_ri, L_xj, n)
    sji  <- hsic_from_kernels(Kc_rj, L_xi, n)
    tobs <- max(sij, sji)

    tp <- numeric(B_PERM)
    for (b in 1:B_PERM) {
      pm    <- sample.int(n)
      tp[b] <- max(hsic_from_kernels(Kc_ri, L_xj[pm, pm], n),
                   hsic_from_kernels(Kc_rj, L_xi[pm, pm], n))
    }
    pv <- (sum(tp >= tobs) + 1) / (B_PERM + 1)

    # Full-sample HSIC scores drive both the effect-size ranking and the
    # final orientation, independent of the sub-sampled p-value.
    sij_f <- hsic_biased(as.matrix(ri), as.matrix(Xj))
    sji_f <- hsic_biased(as.matrix(rj), as.matrix(Xi))

    if (sij_f >= sji_f) { fr <- j; to2 <- i; sc <- sij_f }
    else                { fr <- i; to2 <- j; sc <- sji_f }

    pair_results[[idx]] <- list(i = i, j = j, from = fr, to = to2,
                                score = sc, pvalue = pv,
                                score_sub = tobs,
                                score_ij = sij_f, score_ji = sji_f)

    if (idx %% 11 == 0 || idx == 66)
      cat(sprintf("  %d/66 pairs done\n", idx))
  }
}
hsic_time  <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
total_time <- gam_time + hsic_time
cat(sprintf("  Done in %.0fs. Total: %.0fs.\n\n", hsic_time, total_time))

pvals  <- sapply(pair_results, `[[`, "pvalue")
scores <- sapply(pair_results, `[[`, "score")

cat("P-value distribution:\n")
cat(sprintf("  min = %.4f  median = %.4f  max = %.4f\n",
            min(pvals), median(pvals), max(pvals)))
cat(sprintf("  p < 0.001: %d | p < 0.01: %d | p < 0.05: %d | p >= 0.05: %d\n\n",
            sum(pvals <= 0.001), sum(pvals < 0.01),
            sum(pvals  < 0.05),  sum(pvals >= 0.05)))

# -----------------------------------------------------------------------------
# Stage 3: BH-FDR screening followed by top-k effect-size filter
# -----------------------------------------------------------------------------

cat(sprintf("Stage 3: BH-FDR screening (alpha = %.2f) + top-%d by HSIC...\n",
            SCREEN_ALPHA, K_TARGET))
qvals   <- p.adjust(pvals, method = "BH")
sig_idx <- which(qvals < SCREEN_ALPHA)
cat(sprintf("  BH at alpha = %.2f: %d pairs pass screening\n",
            SCREEN_ALPHA, length(sig_idx)))

sig_scores  <- scores[sig_idx]
top_k_order <- order(-sig_scores)
top_k_idx   <- sig_idx[top_k_order[1:min(K_TARGET, length(sig_idx))]]

cat(sprintf("  Retaining top-%d by HSIC effect size\n\n", length(top_k_idx)))

# Build the adjacency matrix from the retained, oriented edges.
adj  <- matrix(0, n_vars, n_vars, dimnames = list(var_names, var_names))
smat <- matrix(0, n_vars, n_vars, dimnames = list(var_names, var_names))
for (k in seq_along(pair_results)) {
  pr <- pair_results[[k]]
  smat[var_names[pr$j], var_names[pr$i]] <- pr$score_ij
  smat[var_names[pr$i], var_names[pr$j]] <- pr$score_ji
}
for (k in top_k_idx) {
  pr <- pair_results[[k]]
  adj[var_names[pr$from], var_names[pr$to]] <- 1
}
adj <- enforce_acyclicity_wine(adj, smat)
m   <- evaluate_dag_wine(adj, ref_dag)

cat(paste(rep("=", 70), collapse = ""), "\n")
cat("Results\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")
cat(sprintf("F1 = %.3f  Accuracy = %.3f  SHD = %d  MSE = %.3f  Edges = %d  TP = %d  FP = %d  FN = %d\n",
            m["F1"], m["Accuracy"], m["SHD"], m["MSE"],
            m["Edges"], m["TP"], m["FP"], m["FN"]))
cat(sprintf("Precision = %.3f  Recall = %.3f\n", m["Precision"], m["Recall"]))
cat(sprintf("Wall-clock time: %.1f seconds\n\n", total_time))

# -----------------------------------------------------------------------------
# LaTeX output blocks
# -----------------------------------------------------------------------------

cat(paste(rep("=", 70), collapse = ""), "\n")
cat("LaTeX: aggregate performance row\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat(sprintf("Proposed method & %.3f & %.3f & %d  & %.3f & %.1f \\\\\n",
            m["F1"], m["Accuracy"], m["SHD"], m["MSE"], total_time))
cat("\n")

cat(paste(rep("=", 70), collapse = ""), "\n")
cat("LaTeX: ranked list of retained edges\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

edge_list <- list()
for (k in top_k_idx) {
  pr <- pair_results[[k]]
  fn <- var_names[pr$from]; tn <- var_names[pr$to]
  # Only list edges that survive cycle-breaking.
  if (adj[fn, tn] == 1) {
    in_ref <- ref_dag[fn, tn] == 1
    edge_list[[length(edge_list) + 1]] <- list(
      from = fn, to = tn, score = pr$score, pvalue = pr$pvalue, in_ref = in_ref
    )
  }
}
edge_list <- edge_list[order(-sapply(edge_list, `[[`, "score"))]

for (e in edge_list) {
  ref_mark <- if (e$in_ref) "$\\checkmark$" else ""
  cat(sprintf("%s $\\to$ %s & %.4f & %.3f & %s \\\\\n",
              latex_name(e$from), latex_name(e$to),
              e$score, e$pvalue, ref_mark))
}

cat(sprintf("\nEdges after cycle-breaking: %d\n", sum(adj)))
cat(sprintf("Edges before cycle-breaking: %d\n", length(top_k_idx)))

cat("\nEdges removed by cycle-breaking:\n")
for (k in top_k_idx) {
  pr <- pair_results[[k]]
  fn <- var_names[pr$from]; tn <- var_names[pr$to]
  if (adj[fn, tn] == 0) {
    cat(sprintf("  %s -> %s (HSIC = %.4f)\n", fn, tn, pr$score))
  }
}

# -----------------------------------------------------------------------------
# Sensitivity to the BH level (no top-k filter)
# -----------------------------------------------------------------------------

cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("Sensitivity: pure BH-FDR selection at several levels\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
for (alpha in c(0.01, 0.02, 0.05, 0.10, 0.20)) {
  sel <- which(qvals < alpha)
  if (length(sel) == 0) {
    cat(sprintf("  alpha = %.2f: 0 edges\n", alpha))
    next
  }
  adj2 <- matrix(0, n_vars, n_vars, dimnames = list(var_names, var_names))
  for (k in sel) {
    pr <- pair_results[[k]]
    adj2[var_names[pr$from], var_names[pr$to]] <- 1
  }
  adj2 <- enforce_acyclicity_wine(adj2, smat)
  m2   <- evaluate_dag_wine(adj2, ref_dag)
  cat(sprintf("  alpha = %.2f: %d edges, F1 = %.3f, SHD = %d, TP = %d, FP = %d\n",
              alpha, m2["Edges"], m2["F1"], m2["SHD"], m2["TP"], m2["FP"]))
}

cat("\nDone.\n")
