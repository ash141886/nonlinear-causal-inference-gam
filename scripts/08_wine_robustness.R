# =============================================================================
# Robustness Sensitivity Analysis for the Wine Two-Stage Procedure
# =============================================================================
# Re-runs the red-wine two-stage additive-HSIC procedure under a local grid
# of (a) GAM basis-dimension caps (+/-20% around the default k=10) and
# (b) kernel bandwidths (+/-50% around the median-heuristic default).
#
# For each configuration we record:
#   - final F1 / SHD / MSE / Accuracy against the literature reference
#   - the ranked list of the top-10 pairs by full-data HSIC and its
#     Spearman correlation with the default ranking
#   - EDFs of the three partial-effect GAMs used in Figure 9
#
# Addresses the robustness claim in Section 5 (manuscript):
#   "EDF values varied by less than 15% ... Spearman correlation exceeding 0.90".
#
# Output:
#   - results/wine_robustness.csv  (one row per configuration)
#
# Usage: Rscript scripts/08_wine_robustness.R
# =============================================================================

suppressPackageStartupMessages({ library(mgcv) })

dir.create("results", showWarnings = FALSE)

# --- Wine data + reference ---------------------------------------------------

load_wine_data <- function() {
  if (!file.exists("data/winequality-red.csv")) {
    dir.create("data", showWarnings = FALSE, recursive = TRUE)
    download.file(
      "https://archive.ics.uci.edu/ml/machine-learning-databases/wine-quality/winequality-red.csv",
      "data/winequality-red.csv"
    )
  }
  as.data.frame(scale(read.csv("data/winequality-red.csv", sep = ";")))
}

get_reference_dag <- function(var_names) {
  n <- length(var_names)
  dag <- matrix(0, n, n, dimnames = list(var_names, var_names))
  dag["fixed.acidity", "pH"] <- 1
  dag["volatile.acidity", "pH"] <- 1
  dag["citric.acid", "pH"] <- 1
  dag["citric.acid", "fixed.acidity"] <- 1
  dag["residual.sugar", "density"] <- 1
  dag["alcohol", "density"] <- 1
  dag["free.sulfur.dioxide", "total.sulfur.dioxide"] <- 1
  dag["alcohol", "quality"] <- 1
  dag["volatile.acidity", "quality"] <- 1
  dag["sulphates", "quality"] <- 1
  dag
}

# --- Kernel / HSIC helpers (with a bandwidth multiplier) --------------------

rbf_kernel <- function(x, sigma) {
  x <- as.matrix(x)
  exp(-as.matrix(dist(x))^2 / (2 * sigma^2))
}
median_bandwidth <- function(x) {
  x <- as.matrix(x); d <- as.matrix(dist(x))
  med <- median(d[lower.tri(d)]); if (!is.finite(med) || med <= 0) med <- 1; med
}
center_kernel <- function(K) {
  rm <- rowMeans(K)
  K - outer(rm, rep(1, nrow(K))) - outer(rep(1, nrow(K)), rm) + mean(rm)
}
hsic_from_kernels <- function(Kc, L, n) as.numeric(sum(Kc * L) / n^2)

# HSIC with a bandwidth multiplier applied to the median heuristic.
hsic_biased_bw <- function(x, y, bw_mult = 1.0) {
  x <- as.matrix(x); y <- as.matrix(y); n <- nrow(x)
  sx <- bw_mult * median_bandwidth(x)
  sy <- bw_mult * median_bandwidth(y)
  K <- rbf_kernel(x, sx); L <- rbf_kernel(y, sy)
  hsic_from_kernels(center_kernel(K), L, n)
}

# --- Reduced GAM with configurable basis cap --------------------------------

compute_residuals_wine <- function(data, target, exclude, var_names, gam_k_max) {
  n_vars <- ncol(data); n_samples <- nrow(data)
  predictors <- setdiff(seq_len(n_vars), c(target, exclude))
  tn <- var_names[target]
  if (length(predictors) > 0) {
    k_basis <- min(gam_k_max, max(3, floor(n_samples / 40)))
    terms <- sapply(predictors, function(idx) {
      pn <- var_names[idx]
      nu <- length(unique(data[[pn]]))
      ku <- min(k_basis, nu - 1)
      if (ku >= 3) paste0("s(`", pn, "`, k=", ku, ")") else paste0("`", pn, "`")
    })
    fs <- paste0("`", tn, "` ~ ", paste(terms, collapse = " + "))
    model <- tryCatch(
      suppressWarnings(gam(as.formula(fs), data = data, method = "REML", gamma = 1.4)),
      error = function(e) {
        lt <- paste0("`", var_names[predictors], "`", collapse = " + ")
        lm(as.formula(paste0("`", tn, "` ~ ", lt)), data = data)
      })
    residuals(model)
  } else data[[tn]] - mean(data[[tn]])
}

# --- Evaluation --------------------------------------------------------------

evaluate_dag_wine <- function(est, ref) {
  n <- nrow(ref)
  tp <- sum(est == 1 & ref == 1); fp <- sum(est == 1 & ref == 0)
  fn <- sum(est == 0 & ref == 1); tn <- sum(est == 0 & ref == 0) - n
  prec <- if ((tp + fp) > 0) tp / (tp + fp) else 0
  rec  <- if ((tp + fn) > 0) tp / (tp + fn) else 0
  f1 <- if ((prec + rec) > 0) 2 * prec * rec / (prec + rec) else 0
  c(F1 = f1, Accuracy = (tp + tn) / (n * (n - 1)),
    SHD = sum(abs(est - ref)), MSE = mean((est - ref)^2))
}

# --- Cycle breaking ----------------------------------------------------------

find_cycle_wine <- function(adj) {
  n <- nrow(adj); col <- rep("w", n); par <- rep(NA_integer_, n)
  dfs <- function(v) {
    col[v] <<- "g"
    for (u in which(adj[v, ] == 1)) {
      if (col[u] == "g") { p <- u; c <- v
        while (c != u && !is.na(par[c])) { p <- c(p, c); c <- par[c] }
        return(p) }
      if (col[u] == "w") { par[u] <<- v; r <- dfs(u); if (!is.null(r)) return(r) }
    }
    col[v] <<- "b"; NULL
  }
  for (v in 1:n) if (col[v] == "w") { r <- dfs(v); if (!is.null(r)) return(r) }
  NULL
}
enforce_acyclicity_wine <- function(adj, scores) {
  n <- nrow(adj); iter <- 0
  while (iter < n * n) {
    cyc <- find_cycle_wine(adj); if (is.null(cyc)) break
    ms <- Inf; w <- NULL
    for (k in seq_along(cyc)) {
      fr <- cyc[k]; to <- cyc[if (k == length(cyc)) 1 else k + 1]
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

# --- Core two-stage pipeline ------------------------------------------------

run_two_stage <- function(wine_data, gam_k_max = 10, bw_mult = 1.0,
                          screen_alpha = 0.20, k_target = 10,
                          n_sub = 200, b_perm = 1000, seed = 42) {
  var_names <- colnames(wine_data); n_vars <- ncol(wine_data)

  # Stage 1: leave-one-out residuals (full data)
  rc <- list()
  for (i in 1:n_vars) for (j in 1:n_vars) if (i != j)
    rc[[paste(i, j, sep = "_")]] <-
      compute_residuals_wine(wine_data, i, j, var_names, gam_k_max)

  # Stage 2: HSIC + permutation on a subsample
  set.seed(seed)
  sub_idx <- sample.int(nrow(wine_data), n_sub)
  pair_results <- list(); idx <- 0
  for (i in 1:(n_vars - 1)) {
    for (j in (i + 1):n_vars) {
      idx <- idx + 1
      ri <- rc[[paste(i, j, sep = "_")]]; rj <- rc[[paste(j, i, sep = "_")]]
      Xi <- wine_data[[var_names[i]]]; Xj <- wine_data[[var_names[j]]]

      ri_s <- as.matrix(ri[sub_idx]); rj_s <- as.matrix(rj[sub_idx])
      xi_s <- as.matrix(Xi[sub_idx]); xj_s <- as.matrix(Xj[sub_idx])
      n <- n_sub

      set.seed(seed + i * 100 + j)
      Kc_ri <- center_kernel(rbf_kernel(ri_s, bw_mult * median_bandwidth(ri_s)))
      Kc_rj <- center_kernel(rbf_kernel(rj_s, bw_mult * median_bandwidth(rj_s)))
      L_xj  <- rbf_kernel(xj_s, bw_mult * median_bandwidth(xj_s))
      L_xi  <- rbf_kernel(xi_s, bw_mult * median_bandwidth(xi_s))

      sij <- hsic_from_kernels(Kc_ri, L_xj, n)
      sji <- hsic_from_kernels(Kc_rj, L_xi, n)
      tobs <- max(sij, sji)

      tp <- numeric(b_perm)
      for (b in 1:b_perm) {
        pm <- sample.int(n)
        tp[b] <- max(hsic_from_kernels(Kc_ri, L_xj[pm, pm], n),
                      hsic_from_kernels(Kc_rj, L_xi[pm, pm], n))
      }
      pv <- (sum(tp >= tobs) + 1) / (b_perm + 1)

      sij_f <- hsic_biased_bw(as.matrix(ri), as.matrix(Xj), bw_mult)
      sji_f <- hsic_biased_bw(as.matrix(rj), as.matrix(Xi), bw_mult)
      if (sij_f >= sji_f) { fr <- j; to <- i; sc <- sij_f }
      else                { fr <- i; to <- j; sc <- sji_f }

      pair_results[[idx]] <- list(i = i, j = j, from = fr, to = to,
                                   score = sc, pvalue = pv,
                                   score_ij = sij_f, score_ji = sji_f,
                                   pair_key = paste(sort(c(var_names[i], var_names[j])),
                                                     collapse = "--"))
    }
  }

  pvals  <- sapply(pair_results, `[[`, "pvalue")
  scores <- sapply(pair_results, `[[`, "score")

  qvals   <- p.adjust(pvals, method = "BH")
  sig_idx <- which(qvals < screen_alpha)
  top_k   <- sig_idx[order(-scores[sig_idx])[1:min(k_target, length(sig_idx))]]

  adj  <- matrix(0, n_vars, n_vars, dimnames = list(var_names, var_names))
  smat <- matrix(0, n_vars, n_vars, dimnames = list(var_names, var_names))
  for (k in seq_along(pair_results)) {
    pr <- pair_results[[k]]
    smat[var_names[pr$j], var_names[pr$i]] <- pr$score_ij
    smat[var_names[pr$i], var_names[pr$j]] <- pr$score_ji
  }
  for (k in top_k) {
    pr <- pair_results[[k]]
    adj[var_names[pr$from], var_names[pr$to]] <- 1
  }
  adj <- enforce_acyclicity_wine(adj, smat)

  ref_dag <- get_reference_dag(var_names)
  m <- evaluate_dag_wine(adj, ref_dag)

  list(metrics = m,
       scores = scores,
       pair_keys = sapply(pair_results, `[[`, "pair_key"))
}

# --- Partial-effect EDFs (Figure 9) -----------------------------------------

fig9_edfs <- function(wine, gam_k_max = 10) {
  f1 <- gam(density       ~ s(alcohol,     k = gam_k_max),
            data = wine, method = "REML", gamma = 1.4)
  f2 <- gam(fixed.acidity ~ s(citric.acid, k = gam_k_max),
            data = wine, method = "REML", gamma = 1.4)
  f3 <- gam(fixed.acidity ~ s(pH,          k = gam_k_max),
            data = wine, method = "REML", gamma = 1.4)
  c(edf_alc_den = sum(f1$edf),
    edf_cit_fix = sum(f2$edf),
    edf_pH_fix  = sum(f3$edf))
}

# --- Run the grid ------------------------------------------------------------

wine <- load_wine_data()

# Baseline (k = 10, bw_mult = 1.0) -----------------------------------------
cat("Baseline run (k=10, bw_mult=1.0)...\n")
base <- run_two_stage(wine, gam_k_max = 10, bw_mult = 1.0)
base_edfs <- fig9_edfs(wine, 10)

# Grid: basis dim +/- 20%, bandwidth +/- 50% --------------------------------
grid <- expand.grid(
  gam_k_max = c(8, 10, 12),                     # +/-20% around 10
  bw_mult   = c(0.5, 1.0, 1.5),                 # +/-50% around median
  stringsAsFactors = FALSE
)

rows <- list()
for (r in seq_len(nrow(grid))) {
  k_cfg  <- grid$gam_k_max[r]
  bw_cfg <- grid$bw_mult[r]
  cat(sprintf("Config %d/%d: k=%d bw_mult=%.1f\n",
              r, nrow(grid), k_cfg, bw_cfg))
  out <- run_two_stage(wine, gam_k_max = k_cfg, bw_mult = bw_cfg)
  edfs <- fig9_edfs(wine, k_cfg)

  # Spearman rank correlation between this and the baseline score ordering
  rho <- suppressWarnings(cor(out$scores, base$scores, method = "spearman"))

  # EDF deviations relative to baseline (as fractions)
  rel_edf <- abs(edfs - base_edfs) / base_edfs

  rows[[r]] <- data.frame(
    gam_k_max       = k_cfg,
    bw_mult         = bw_cfg,
    F1              = unname(out$metrics["F1"]),
    Accuracy        = unname(out$metrics["Accuracy"]),
    SHD             = unname(out$metrics["SHD"]),
    MSE             = unname(out$metrics["MSE"]),
    spearman_vs_base = rho,
    edf_alc_den     = unname(edfs["edf_alc_den"]),
    edf_cit_fix     = unname(edfs["edf_cit_fix"]),
    edf_pH_fix      = unname(edfs["edf_pH_fix"]),
    rel_edf_max     = max(rel_edf, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}
df <- do.call(rbind, rows)

write.csv(df, "results/wine_robustness.csv", row.names = FALSE)

# --- Summary -----------------------------------------------------------------

cat("\n-- Robustness summary (wine, two-stage) --\n")
cat(sprintf("Baseline EDFs: alc->den=%.2f, citric->fix=%.2f, pH->fix=%.2f\n",
            base_edfs["edf_alc_den"], base_edfs["edf_cit_fix"], base_edfs["edf_pH_fix"]))

non_base <- df$gam_k_max != 10 | df$bw_mult != 1.0
cat(sprintf("Max relative EDF deviation across grid (excluding baseline): %.1f%%\n",
            100 * max(df$rel_edf_max[non_base])))
cat(sprintf("Min Spearman rank-correlation of pair scores vs baseline  : %.3f\n",
            min(df$spearman_vs_base[non_base])))
cat(sprintf("F1 range across grid                                       : [%.3f, %.3f]\n",
            min(df$F1), max(df$F1)))
cat(sprintf("SHD range across grid                                      : [%d, %d]\n",
            min(df$SHD), max(df$SHD)))

cat("\nWrote results/wine_robustness.csv\n")
cat("Done.\n")
