# =============================================================================
# Wine-analysis figures
# =============================================================================
#
# Produces two wine-related figures for the red-wine case study:
#
#   * \code{figures/wine_dag.pdf} -- side-by-side comparison of the
#     estimated and literature-informed reference DAGs, rendered on a
#     shared hand-tuned layout so that the two panels can be read
#     against each other without visual realignment.
#
#   * \code{figures/wine_smooths.pdf} -- partial-effect panels for three
#     representative smooths (alcohol -> density, citric acid -> fixed
#     acidity, pH -> fixed acidity). Each panel reports the effective
#     degrees of freedom of the fitted smooth, a compact indicator of
#     how non-linear the component is.
#
# If \code{results/wine_red_results.rds} is unavailable, the estimated
# adjacency is reconstructed from the ten-edge two-stage solution used
# in the rest of the analysis.
#
# Usage:  Rscript scripts/06_plot_wine.R
# Output: figures/wine_dag.pdf, figures/wine_smooths.pdf
# =============================================================================

library(ggplot2)
library(igraph)
library(mgcv)
library(gridExtra)

dir.create("figures", showWarnings = FALSE)

source("R/utils.R")

# -----------------------------------------------------------------------------
# Load data and reference graph
# -----------------------------------------------------------------------------

wine      <- load_wine_data(color = "red", data_dir = "data")
var_names <- colnames(wine)
ref_dag   <- get_wine_reference_dag(var_names)

# =============================================================================
# Figure: estimated vs. reference DAG
# =============================================================================

cat("Generating DAG comparison figure...\n")

# Use the cached adjacency if available; otherwise fall back to the
# ten-edge two-stage solution.
if (file.exists("results/wine_red_results.rds")) {
  wine_res <- readRDS("results/wine_red_results.rds")
  est_dag  <- wine_res$proposed$result$adjacency
} else {
  cat("  No saved .rds found. Using the ten-edge two-stage solution.\n")
  est_dag <- matrix(0, length(var_names), length(var_names),
                    dimnames = list(var_names, var_names))
  est_dag["free.sulfur.dioxide", "total.sulfur.dioxide"] <- 1
  est_dag["alcohol",              "density"]             <- 1
  est_dag["volatile.acidity",     "citric.acid"]         <- 1
  est_dag["pH",                   "fixed.acidity"]       <- 1
  est_dag["residual.sugar",       "density"]             <- 1
  est_dag["alcohol",              "residual.sugar"]      <- 1
  est_dag["density",              "fixed.acidity"]       <- 1
  est_dag["pH",                   "density"]             <- 1
  est_dag["alcohol",              "fixed.acidity"]       <- 1
  est_dag["citric.acid",          "fixed.acidity"]       <- 1
}

# Short node labels for display.
short_names <- c(
  "fixed.acidity"        = "Fixed\nAcidity",
  "volatile.acidity"     = "Volatile\nAcidity",
  "citric.acid"          = "Citric\nAcid",
  "residual.sugar"       = "Residual\nSugar",
  "chlorides"            = "Chlorides",
  "free.sulfur.dioxide"  = "Free\nSO2",
  "total.sulfur.dioxide" = "Total\nSO2",
  "density"              = "Density",
  "pH"                   = "pH",
  "sulphates"            = "Sulphates",
  "alcohol"              = "Alcohol",
  "quality"              = "Quality"
)

# Hand-tuned layout chosen to minimise edge-node overlap. The rows
# follow the same order as \code{var_names}.
manual_layout <- matrix(c(
  -0.3,  1.0,   # fixed.acidity       (top-centre-left)
  -0.9,  0.6,   # volatile.acidity    (upper-left)
   0.5,  0.8,   # citric.acid         (upper-right)
   0.7,  0.2,   # residual.sugar      (mid-right)
  -1.2,  0.1,   # chlorides           (far-left)
   1.0, -0.4,   # free.sulfur.dioxide (lower-right)
   1.0, -0.9,   # total.sulfur.dioxide (bottom-right)
   0.2,  0.0,   # density             (centre; spaced from edges)
  -0.8, -0.1,   # pH                  (mid-left)
  -0.6, -0.6,   # sulphates           (lower-left)
  -0.2, -0.5,   # alcohol             (lower-centre)
  -0.2, -1.0    # quality             (bottom-centre)
), ncol = 2, byrow = TRUE)

#' Render a DAG with edges coloured by membership in a reference graph.
#'
#' Edges that also appear in \code{ref} are drawn in green; additional
#' edges (false positives with respect to \code{ref}) are drawn in red.
plot_dag <- function(adj, ref, title, layout) {
  g <- graph_from_adjacency_matrix(adj, mode = "directed")
  V(g)$label <- short_names[V(g)$name]

  edge_colors <- sapply(seq_len(ecount(g)), function(e) {
    ends <- ends(g, e)
    from <- ends[1]; to <- ends[2]
    if (ref[from, to] == 1) "#2ca02c" else "#d62728"
  })

  plot(g, layout = layout,
       edge.color   = edge_colors, edge.width = 2.5, edge.arrow.size = 0.6,
       edge.curved  = autocurve.edges(g, start = 0.2),
       vertex.color = "#dceefb", vertex.frame.color = "#4a90d9",
       vertex.size  = 28, vertex.label.cex = 0.55,
       vertex.label.color = "black", vertex.label.font = 2,
       main         = title)
}

pdf("figures/wine_dag.pdf", width = 14, height = 7)
par(mfrow = c(1, 2), mar = c(1, 1, 3, 1))
plot_dag(est_dag, ref_dag, "(a) Estimated DAG (proposed method)", manual_layout)
legend("bottomleft",
       legend = c("True positive (matches reference)",
                  "False positive (additional edge)"),
       col    = c("#2ca02c", "#d62728"),
       lwd    = 2.5, cex = 0.8, bty = "n")

ref_g <- graph_from_adjacency_matrix(ref_dag, mode = "directed")
V(ref_g)$label <- short_names[V(ref_g)$name]
plot(ref_g, layout = manual_layout,
     edge.color   = "#4a90d9", edge.width = 2.5, edge.arrow.size = 0.6,
     edge.curved  = autocurve.edges(ref_g, start = 0.2),
     vertex.color = "#dceefb", vertex.frame.color = "#4a90d9",
     vertex.size  = 28, vertex.label.cex = 0.55,
     vertex.label.color = "black", vertex.label.font = 2,
     main         = "(b) Literature-informed reference DAG")
legend("bottomleft",
       legend = c("Reference edge", "Variable node"),
       col    = c("#4a90d9", "#dceefb"),
       lwd    = c(2.5, NA), pch = c(NA, 21),
       pt.bg  = "#dceefb", cex = 0.8, bty = "n")
dev.off()
cat("Saved: figures/wine_dag.pdf\n")

# =============================================================================
# Figure: partial-effect smooths for three representative edges
# =============================================================================

cat("Generating partial-effect figure...\n")

# Fit three univariate-smooth GAMs for the representative components.
fit1 <- gam(density       ~ s(alcohol,     k = 10), data = wine, method = "REML", gamma = 1.4)
fit2 <- gam(fixed.acidity ~ s(citric.acid, k = 10), data = wine, method = "REML", gamma = 1.4)
fit3 <- gam(fixed.acidity ~ s(pH,          k = 10), data = wine, method = "REML", gamma = 1.4)

edf1 <- round(sum(fit1$edf), 1)
edf2 <- round(sum(fit2$edf), 1)
edf3 <- round(sum(fit3$edf), 1)

pdf("figures/wine_smooths.pdf", width = 12, height = 4)
par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3, 1))

plot(fit1, shade = TRUE, shade.col = "gray85",
     xlab = "Alcohol (std.)",      ylab = "Density (partial)",
     main = paste0("EDF = ", edf1))

plot(fit2, shade = TRUE, shade.col = "gray85",
     xlab = "Citric acid (std.)",  ylab = "Fixed acidity (partial)",
     main = paste0("EDF = ", edf2))

plot(fit3, shade = TRUE, shade.col = "gray85",
     xlab = "pH (std.)",           ylab = "Fixed acidity (partial)",
     main = paste0("EDF = ", edf3))

dev.off()
cat("Saved: figures/wine_smooths.pdf\n")

cat(sprintf("\nEDF values: alcohol -> density = %.1f, citric -> fixed = %.1f, pH -> fixed = %.1f\n",
            edf1, edf2, edf3))
cat("\nAll wine figures generated.\n")
