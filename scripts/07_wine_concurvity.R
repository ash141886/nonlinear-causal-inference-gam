# =============================================================================
# Concurvity Diagnostic for the Wine Leave-One-Out Additive Models
# =============================================================================
# Addresses Remark 1 in the manuscript (Section 3.1): reports mgcv::concurvity()
# for every (target, excluded-predictor) pair fitted on the red-wine data.
#
# Output:
#   - results/wine_concurvity.csv    (one row per reduced model)
#   - results/wine_concurvity_summary.txt
#
# Usage:  Rscript scripts/07_wine_concurvity.R
# =============================================================================

suppressPackageStartupMessages({ library(mgcv) })

dir.create("results", showWarnings = FALSE)

# --- Data --------------------------------------------------------------------

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

wine_data <- load_wine_data()
var_names <- colnames(wine_data)
n_vars    <- ncol(wine_data)
n_samples <- nrow(wine_data)

cat(sprintf("Wine red data: %d x %d\n\n", n_samples, n_vars))

# --- Leave-one-out reduced model + concurvity --------------------------------

fit_reduced_gam <- function(target, exclude, gam_k_max = 10) {
  predictors <- setdiff(seq_len(n_vars), c(target, exclude))
  target_name <- var_names[target]
  k_basis <- min(gam_k_max, max(3, floor(n_samples / 40)))

  terms <- sapply(predictors, function(idx) {
    pn <- var_names[idx]
    nu <- length(unique(wine_data[[pn]]))
    ku <- min(k_basis, nu - 1)
    if (ku >= 3) paste0("s(`", pn, "`, k=", ku, ")") else paste0("`", pn, "`")
  })

  fs <- paste0("`", target_name, "` ~ ", paste(terms, collapse = " + "))
  suppressWarnings(
    gam(as.formula(fs), data = wine_data, method = "REML", gamma = 1.4)
  )
}

# mgcv::concurvity() returns a 3 x (n_smooth + 1) matrix with rows
# "worst", "observed", "estimate".  We record the maximum across non-"para"
# smooths (the intercept column is excluded).
max_concurvity <- function(m) {
  cc <- tryCatch(mgcv::concurvity(m, full = TRUE), error = function(e) NULL)
  if (is.null(cc)) return(c(worst = NA_real_, observed = NA_real_, estimate = NA_real_))
  # Drop the "para" column if present
  cc <- cc[, setdiff(colnames(cc), "para"), drop = FALSE]
  if (ncol(cc) == 0) return(c(worst = 0, observed = 0, estimate = 0))
  apply(cc, 1, max, na.rm = TRUE)
}

# --- Loop over all (target, excluded) pairs ----------------------------------

cat("Fitting 12 x 11 = 132 reduced additive models and computing concurvity...\n")
rows <- list(); k <- 0
for (i in seq_len(n_vars)) {
  for (j in seq_len(n_vars)) {
    if (i == j) next
    k <- k + 1
    m <- fit_reduced_gam(target = i, exclude = j)
    mc <- max_concurvity(m)
    rows[[k]] <- data.frame(
      target   = var_names[i],
      excluded = var_names[j],
      conc_worst    = unname(mc["worst"]),
      conc_observed = unname(mc["observed"]),
      conc_estimate = unname(mc["estimate"]),
      stringsAsFactors = FALSE
    )
    if (k %% 22 == 0) cat(sprintf("  %d / 132 fits done\n", k))
  }
}
df <- do.call(rbind, rows)

write.csv(df, "results/wine_concurvity.csv", row.names = FALSE)

# --- Summary ------------------------------------------------------------------

summ <- function(x) c(min = min(x, na.rm = TRUE),
                      median = median(x, na.rm = TRUE),
                      mean = mean(x, na.rm = TRUE),
                      max = max(x, na.rm = TRUE),
                      frac_above_0.6 = mean(x >= 0.6, na.rm = TRUE),
                      frac_above_0.8 = mean(x >= 0.8, na.rm = TRUE))

sink("results/wine_concurvity_summary.txt")
cat("Wine red-wine leave-one-out reduced-GAM concurvity diagnostics\n")
cat("Data: n = 1599, p = 12. 132 reduced models (one per (target, excluded) pair).\n\n")
cat("'worst'    : worst-case (full-basis) concurvity index of smooth terms\n")
cat("'observed' : concurvity among the smooths that are actually used\n")
cat("'estimate' : data-driven concurvity estimate\n\n")

for (nm in c("conc_worst", "conc_observed", "conc_estimate")) {
  cat(sprintf("%s:\n", nm))
  s <- summ(df[[nm]])
  cat(sprintf("  min=%.3f  median=%.3f  mean=%.3f  max=%.3f\n",
              s["min"], s["median"], s["mean"], s["max"]))
  cat(sprintf("  fraction of reduced models with index >= 0.6 : %.3f\n",
              s["frac_above_0.6"]))
  cat(sprintf("  fraction of reduced models with index >= 0.8 : %.3f\n\n",
              s["frac_above_0.8"]))
}

cat("Models with conc_observed >= 0.6 (if any):\n")
bad <- df[!is.na(df$conc_observed) & df$conc_observed >= 0.6, ]
if (nrow(bad) == 0) {
  cat("  (none)\n")
} else {
  print(bad, row.names = FALSE)
}
sink()

cat("\nWrote results/wine_concurvity.csv and results/wine_concurvity_summary.txt\n")
cat("Done.\n")
