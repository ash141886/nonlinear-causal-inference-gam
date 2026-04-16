# =============================================================================
# Simulation figures from the non-linear study
# =============================================================================
#
# Reads the long-format results saved by
# \code{scripts/01_simulation_nonlinear.R} and produces six two-panel
# figures, one per performance metric. Each figure shows the median of
# the metric across Monte Carlo replicates, first as a function of the
# sample size (panels faceted by the number of variables p) and then as
# a function of the number of variables (panels faceted by the sample
# size n). The styling aims at a clean, journal-quality look: italic
# math symbols for n and p via plotmath, en-dash in compound names,
# explicit x-axis ticks that match the simulation grid, and a base font
# size that remains legible under typical two-column reductions.
#
# The linear-case residual diagnostics
# (\code{figures/linear_diagnostics.pdf}) are produced directly by
# \code{scripts/02_simulation_linear.R}, which runs its own simulation
# and plotting in a single pass; nothing further is done for it here.
#
# Usage:  Rscript scripts/05_plot_simulations.R
# Output: figures/f1.pdf, accuracy.pdf, misoriented.pdf, shd.pdf,
#         mse.pdf, time.pdf
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

dir.create("figures", showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Load simulation results
# -----------------------------------------------------------------------------

sim_nl <- readRDS("results/sim_nonlinear.rds")

# -----------------------------------------------------------------------------
# Global plotting theme and palette
# -----------------------------------------------------------------------------

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

# Italic math for the axis and facet labels.
xlab_n <- expression("Sample size," ~ italic(n))
xlab_p <- expression("Number of variables," ~ italic(p))

facet_n_labeller <- as_labeller(function(v) paste0("italic(n)==", v), default = label_parsed)
facet_p_labeller <- as_labeller(function(v) paste0("italic(p)==", v), default = label_parsed)

# -----------------------------------------------------------------------------
# Helper: two-panel plot of one metric (vs. n and vs. p)
# -----------------------------------------------------------------------------

make_dual_plot <- function(data, metric, ylabel, title_prefix,
                           filename, width = 10, height = 6) {

  medians <- data %>%
    group_by(p, n, method) %>%
    summarise(value = median(.data[[metric]], na.rm = TRUE), .groups = "drop")

  # Upper panel: metric vs. sample size, faceted by p.
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

  # Lower panel: metric vs. number of variables, faceted by n.
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

# -----------------------------------------------------------------------------
# Produce the six metric figures
# -----------------------------------------------------------------------------

cat("Generating simulation figures...\n\n")

make_dual_plot(sim_nl, "F1",          "Directed F1",                "Directed F1",
               "f1.pdf")
make_dual_plot(sim_nl, "Accuracy",    "Graph accuracy",             "Graph accuracy",
               "accuracy.pdf")
make_dual_plot(sim_nl, "Misoriented", "Misoriented edges",          "Misoriented edges",
               "misoriented.pdf")
make_dual_plot(sim_nl, "SHD",         "Structural Hamming distance","Structural Hamming distance",
               "shd.pdf")
make_dual_plot(sim_nl, "MSE",         "Adjacency MSE",              "Adjacency MSE",
               "mse.pdf")
make_dual_plot(sim_nl, "Time",        "Runtime (seconds)",          "Runtime",
               "time.pdf")

cat("\nAll simulation figures generated.\n")
