# =============================================================================
# Plotting: Wine Figures (Figures 8-9)
# =============================================================================
#
# Figure 8: Estimated DAG vs. reference DAG (shared layout)
# Figure 9: Partial effect functions (alcohol->density, citric->fixed, pH->fixed)
#
# Usage: Rscript scripts/06_plot_wine.R
# Output: figures/wine_dag_R5_shared.pdf, figures/wine_smooths_R5.pdf
# =============================================================================

library(ggplot2)
library(igraph)
library(mgcv)
library(gridExtra)

dir.create("figures", showWarnings = FALSE)

source("R/utils.R")

# --- Load wine data ----------------------------------------------------------

wine <- load_wine_data(color = "red", data_dir = "data")
var_names <- colnames(wine)
ref_dag <- get_wine_reference_dag(var_names)

# --- Figure 8: DAG Comparison ------------------------------------------------

cat("Generating Figure 8: DAG comparison...\n")

# Load estimated adjacency (from 03b results or hard-coded from Table 3)
if (file.exists("results/wine_red_results.rds")) {
  wine_res <- readRDS("results/wine_red_results.rds")
  est_dag <- wine_res$proposed$result$adjacency
} else {
  cat("  No saved .rds found. Using Table 3 edges from manuscript.\n")
  # 10 edges from Table 3 (two-stage: BH-FDR alpha=0.20 + top-10 HSIC)
  est_dag <- matrix(0, length(var_names), length(var_names),
                    dimnames = list(var_names, var_names))
  est_dag["free.sulfur.dioxide", "total.sulfur.dioxide"] <- 1
  est_dag["alcohol", "density"]                          <- 1
  est_dag["volatile.acidity", "citric.acid"]             <- 1
  est_dag["pH", "fixed.acidity"]                         <- 1
  est_dag["residual.sugar", "density"]                   <- 1
  est_dag["alcohol", "residual.sugar"]                   <- 1
  est_dag["density", "fixed.acidity"]                    <- 1
  est_dag["pH", "density"]                               <- 1
  est_dag["alcohol", "fixed.acidity"]                    <- 1
  est_dag["citric.acid", "fixed.acidity"]                <- 1
}

# Short variable names for plotting
short_names <- c(
  "fixed.acidity" = "Fixed\nAcidity",
  "volatile.acidity" = "Volatile\nAcidity",
  "citric.acid" = "Citric\nAcid",
  "residual.sugar" = "Residual\nSugar",
  "chlorides" = "Chlorides",
  "free.sulfur.dioxide" = "Free\nSO2",
  "total.sulfur.dioxide" = "Total\nSO2",
  "density" = "Density",
  "pH" = "pH",
  "sulphates" = "Sulphates",
  "alcohol" = "Alcohol",
  "quality" = "Quality"
)

# Manual layout: hand-tuned coordinates to avoid edge-node crowding
# Order: fixed.acidity, volatile.acidity, citric.acid, residual.sugar,
#        chlorides, free.sulfur.dioxide, total.sulfur.dioxide, density,
#        pH, sulphates, alcohol, quality
manual_layout <- matrix(c(
  -0.3,  1.0,   # fixed.acidity      (top-center-left)
  -0.9,  0.6,   # volatile.acidity   (upper-left)
   0.5,  0.8,   # citric.acid        (upper-right)
   0.7,  0.2,   # residual.sugar     (mid-right)
  -1.2,  0.1,   # chlorides          (far-left)
   1.0, -0.4,   # free.sulfur.dioxide (lower-right)
   1.0, -0.9,   # total.sulfur.dioxide (bottom-right)
   0.2,  0.0,   # density            (center — spaced from edges)
  -0.8, -0.1,   # pH                 (mid-left)
  -0.6, -0.6,   # sulphates          (lower-left)
  -0.2, -0.5,   # alcohol            (lower-center)
  -0.2, -1.0    # quality            (bottom-center)
), ncol = 2, byrow = TRUE)

plot_dag <- function(adj, ref, title, layout) {
  g <- graph_from_adjacency_matrix(adj, mode = "directed")
  V(g)$label <- short_names[V(g)$name]

  # Color edges: green if in reference, red if not
  edge_colors <- sapply(seq_len(ecount(g)), function(e) {
    ends <- ends(g, e)
    from <- ends[1]; to <- ends[2]
    if (ref[from, to] == 1) "#2ca02c" else "#d62728"
  })

  # Auto-curve edges where multiple edges share endpoints to avoid overlap
  plot(g, layout = layout,
       edge.color = edge_colors, edge.width = 2.5, edge.arrow.size = 0.6,
       edge.curved = autocurve.edges(g, start = 0.2),
       vertex.color = "#dceefb", vertex.frame.color = "#4a90d9",
       vertex.size = 28, vertex.label.cex = 0.55,
       vertex.label.color = "black", vertex.label.font = 2,
       main = title)
}

pdf("figures/wine_dag_R5_shared.pdf", width = 14, height = 7)
par(mfrow = c(1, 2), mar = c(1, 1, 3, 1))
plot_dag(est_dag, ref_dag, "(a) Estimated DAG (Proposed Method)", manual_layout)
legend("bottomleft", legend = c("True positive (matches reference)",
                                 "False positive (additional edge)"),
       col = c("#2ca02c", "#d62728"), lwd = 2.5, cex = 0.8, bty = "n")

# Reference DAG
ref_g <- graph_from_adjacency_matrix(ref_dag, mode = "directed")
V(ref_g)$label <- short_names[V(ref_g)$name]
plot(ref_g, layout = manual_layout,
     edge.color = "#4a90d9", edge.width = 2.5, edge.arrow.size = 0.6,
     edge.curved = autocurve.edges(ref_g, start = 0.2),
     vertex.color = "#dceefb", vertex.frame.color = "#4a90d9",
     vertex.size = 28, vertex.label.cex = 0.55,
     vertex.label.color = "black", vertex.label.font = 2,
     main = "(b) Literature-Based Reference DAG")
legend("bottomleft", legend = c("Reference edge", "Variable node"),
       col = c("#4a90d9", "#dceefb"), lwd = c(2.5, NA), pch = c(NA, 21),
       pt.bg = "#dceefb", cex = 0.8, bty = "n")
dev.off()
cat("Saved: wine_dag_R5_shared.pdf\n")

# --- Figure 9: Partial Effect Functions --------------------------------------

cat("Generating Figure 9: Partial effects...\n")

# Fit the three GAMs for selected edges
# 1. Alcohol -> Density
fit1 <- gam(density ~ s(alcohol, k = 10), data = wine, method = "REML", gamma = 1.4)
# 2. Citric acid -> Fixed acidity
fit2 <- gam(fixed.acidity ~ s(citric.acid, k = 10), data = wine, method = "REML", gamma = 1.4)
# 3. pH -> Fixed acidity
fit3 <- gam(fixed.acidity ~ s(pH, k = 10), data = wine, method = "REML", gamma = 1.4)

edf1 <- round(sum(fit1$edf), 1)
edf2 <- round(sum(fit2$edf), 1)
edf3 <- round(sum(fit3$edf), 1)

pdf("figures/wine_smooths_R5.pdf", width = 12, height = 4)
par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3, 1))

plot(fit1, shade = TRUE, shade.col = "gray85",
     xlab = "Alcohol (std.)", ylab = "Density (partial)",
     main = paste0("EDF = ", edf1))

plot(fit2, shade = TRUE, shade.col = "gray85",
     xlab = "Citric acid (std.)", ylab = "Fixed acidity (partial)",
     main = paste0("EDF = ", edf2))

plot(fit3, shade = TRUE, shade.col = "gray85",
     xlab = "pH (std.)", ylab = "Fixed acidity (partial)",
     main = paste0("EDF = ", edf3))

dev.off()
cat("Saved: wine_smooths_R5.pdf\n")

cat(sprintf("\nEDF values: Alcohol->Density=%.1f, Citric->Fixed=%.1f, pH->Fixed=%.1f\n",
            edf1, edf2, edf3))
cat("\nAll wine figures generated.\n")
