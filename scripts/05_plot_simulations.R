# =============================================================================
# Plotting: Non-linear Simulation Figures (Figures 1-6) -- Publication-grade
# =============================================================================
#
# Reads results from scripts/01_simulation_nonlinear.R and generates the
# non-linear simulation figures in the journal-ready style: italic math symbols
# for n and p (via plotmath), en-dash in "Kolmogorov-Smirnov", explicit x-axis
# ticks matching the simulation grid, clean captions, and a unified theme with
# adequate base font size for single-column / 180 mm PDF reproduction.
#
# Figure 7 (linear-case diagnostics, Linear_R5.pdf) is produced directly by
# scripts/02_simulation_linear.R and is NOT drawn here.
#
# Usage:  Rscript scripts/05_plot_simulations.R
# Output: figures/f1_r34.pdf, accuracy_r34.pdf, misoriented_r34.pdf,
#         shd_r34.pdf, mse_r34.pdf, time_r34.pdf
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

dir.create("figures", showWarnings = FALSE)

# --- Load Results ------------------------------------------------------------

sim_nl <- readRDS("results/sim_nonlinear.rds")

# --- Publication Theme -------------------------------------------------------

theme_paper <- theme_bw(base_size = 12) +
  theme(
    text              = element_text(size = 12, colour = "black"),
    axis.title        = element_text(size = 12),
    axis.text         = element_text(size = 10, colour = "black"),
    axis.ticks        = element_line(colour = "black", linewidth = 0.35),
    panel.border      = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    panel.grid.minor  = element_blank(),
    panel.grid.major  = element_line(colour = "grey90", linewidth = 0.3),
    plot.title        = element_text(size = 12, face = "plain", hjust = 0),
    strip.background  = element_rect(fill = "white", colour = "black", linewidth = 0.4),
    strip.text        = element_text(size = 11, face = "plain"),
    legend.position   = "bottom",
    legend.title      = element_blank(),
    legend.key.width  = unit(1.2, "cm"),
    legend.text       = element_text(size = 11),
    plot.margin       = margin(6, 8, 4, 6)
  )

method_colors    <- c("Proposed" = "#1f77b4", "LiNGAM" = "#d62728")
method_linetypes <- c("Proposed" = "dashed",  "LiNGAM" = "solid")
method_shapes    <- c("Proposed" = 17,        "LiNGAM" = 16)

# Italic math labels (plotmath): x-axis uses italic n / p, facet strips use
# "n = <value>" and "p = <value>" with italic symbol.
xlab_n <- expression("Sample size," ~ italic(n))
xlab_p <- expression("Number of variables," ~ italic(p))

facet_n_labeller <- as_labeller(function(v) paste0("italic(n)==", v), default = label_parsed)
facet_p_labeller <- as_labeller(function(v) paste0("italic(p)==", v), default = label_parsed)

# --- Helper: Two-panel plot (metric vs. n, metric vs. p) ---------------------

make_dual_plot <- function(data, metric, ylabel, title_prefix,
                           filename, width = 10, height = 6) {

  medians <- data %>%
    group_by(p, n, method) %>%
    summarise(value = median(.data[[metric]], na.rm = TRUE), .groups = "drop")

  # Panel 1: metric vs. sample size, faceted by p (italic p = <value>)
  p1 <- ggplot(medians, aes(x = n, y = value, color = method,
                            linetype = method, shape = method)) +
    geom_line(linewidth = 0.7) + geom_point(size = 2.2) +
    facet_wrap(~p, nrow = 1, labeller = facet_p_labeller) +
    scale_color_manual(values = method_colors) +
    scale_linetype_manual(values = method_linetypes) +
    scale_shape_manual(values = method_shapes) +
    scale_x_continuous(breaks = c(400, 800, 1200, 1600)) +
    labs(x = xlab_n, y = ylabel,
         title = bquote(.(title_prefix) ~ "vs. sample size" ~ italic(n))) +
    theme_paper

  # Panel 2: metric vs. number of variables, faceted by n (italic n = <value>)
  p2 <- ggplot(medians, aes(x = p, y = value, color = method,
                            linetype = method, shape = method)) +
    geom_line(linewidth = 0.7) + geom_point(size = 2.2) +
    facet_wrap(~n, nrow = 1, labeller = facet_n_labeller) +
    scale_color_manual(values = method_colors) +
    scale_linetype_manual(values = method_linetypes) +
    scale_shape_manual(values = method_shapes) +
    scale_x_continuous(breaks = c(4, 8, 12, 16)) +
    labs(x = xlab_p, y = ylabel,
         title = bquote(.(title_prefix) ~ "vs. number of variables" ~ italic(p))) +
    theme_paper

  pdf(file.path("figures", filename), width = width, height = height, useDingbats = FALSE)
  grid.arrange(p1, p2, nrow = 2)
  dev.off()
  cat("Saved:", filename, "\n")
}

# --- Generate Figures 1-6 ----------------------------------------------------

cat("Generating simulation figures...\n\n")

# Figure 1: Directed F1 score
make_dual_plot(sim_nl, "F1", "Directed F1", "Directed F1", "f1_r34.pdf")

# Figure 2: Graph accuracy
make_dual_plot(sim_nl, "Accuracy", "Graph accuracy", "Graph accuracy",
               "accuracy_r34.pdf")

# Figure 3: Misoriented edges
make_dual_plot(sim_nl, "Misoriented", "Misoriented edges", "Misoriented edges",
               "misoriented_r34.pdf")

# Figure 4: Structural Hamming distance
make_dual_plot(sim_nl, "SHD", "Structural Hamming distance",
               "Structural Hamming distance", "shd_r34.pdf")

# Figure 5: Adjacency MSE
make_dual_plot(sim_nl, "MSE", "Adjacency MSE", "Adjacency MSE", "mse_r34.pdf")

# Figure 6: Runtime
make_dual_plot(sim_nl, "Time", "Runtime (seconds)", "Runtime", "time_r34.pdf")

# Note: Figure 7 (linear-case residual diagnostics, Linear_R5.pdf) is produced
# directly by scripts/02_simulation_linear.R, which runs its own simulation
# and plotting in a single pass. No action is needed here.

cat("\nAll simulation figures generated.\n")
