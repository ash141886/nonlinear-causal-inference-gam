# =============================================================================
# Simulation Study: Linear Data (Section 4.2, Figure 7)
# =============================================================================
#
# Reproduces Figure 7 of Islam & Suzuki (JJSD, 2026): cross-method agreement
# diagnostics on a 4-node linear SEM with Laplace noise.
#
# Usage:  Rscript scripts/02_simulation_linear.R
# Output: figures/Linear_R5.pdf, results/linear_diag_R.csv
#
# Structural model (Laplace noise, scale 1/sqrt(2), unit variance):
#   x1 =              e1
#   x2 = 2*x1       + e2
#   x4 =     x2     + e4
#   x3 =         x4 + e3
#
# For each variable at each n:
#   LiNGAM residuals       : e_hat = (I - B_ica) x,  B_ica from ICA whitening
#   Proposed GAM residuals : leave-one-out GAM of x on remaining 3 via mgcv
#
# Metrics reported (mean over N_REP replications, with Monte Carlo SE):
#   (a) Fitted-component norm ||x - residual||_2 / sqrt(n)      <-- Panel 1
#   (b) True residual RMSE  sqrt(mean(residual^2))              <-- bonus (CSV only)
#   (c) Pearson correlation between LiNGAM and Proposed residuals <-- Panel 2
#   (d) Two-sample KS statistic + p-value (cross-method agreement) <-- Panels 3-4
#   (e) Paired t statistic + p-value                            <-- CSV only
#       (dropped from the figure: averaging t/p across reps is not a valid
#        summary, and the KS panels already test distributional agreement.)
#
# NOTE on terminology: the "RMSE" axis in the original Python code and the
# manuscript's Fig 7 caption is computed as sqrt(mean((x - residual)^2)), which
# equals ||fitted||_2 / sqrt(n), NOT the standard residual RMSE. This script
# reports BOTH quantities so the reader can verify numbers either way.
#
# Dependencies: mgcv (only).  DirectLiNGAM is implemented from scratch below,
# so no `lingam` / `fastICA` packages are needed.
# Run time: ~3-5 min for 50 replications x 4 sample sizes on a laptop.
# =============================================================================

suppressPackageStartupMessages({
  library(mgcv)
})

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
SAMPLE_SIZES <- c(400, 800, 1200, 1600)
N_REP        <- 50            # Monte Carlo replications per sample size
BASE_SEED    <- 42
OUT_DIR_RES  <- "results"
OUT_DIR_FIG  <- "figures"
dir.create(OUT_DIR_RES, showWarnings = FALSE)
dir.create(OUT_DIR_FIG, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Laplace sampler (base R has no rlaplace; inverse-CDF trick)
# -----------------------------------------------------------------------------
rlaplace <- function(n, location = 0, scale = 1) {
  u <- runif(n, -0.5, 0.5)
  location - scale * sign(u) * log(1 - 2 * abs(u))
}

# -----------------------------------------------------------------------------
# Data generation
# -----------------------------------------------------------------------------
generate_data <- function(n_samples) {
  scale <- 1 / sqrt(2)
  x1 <- rlaplace(n_samples, 0, scale)
  x2 <- 2 * x1 + rlaplace(n_samples, 0, scale)
  x4 <- x2     + rlaplace(n_samples, 0, scale)
  x3 <- x4     + rlaplace(n_samples, 0, scale)
  list(x1 = x1, x2 = x2, x3 = x3, x4 = x4)
}

# -----------------------------------------------------------------------------
# DirectLiNGAM (Shimizu et al. 2011; Hyvarinen & Smith 2013)
# -----------------------------------------------------------------------------
# Returns (i) the inferred causal order as a permutation of 1:p, and
#         (ii) the strictly lower-triangular-in-order B matrix estimated by
#              ordinary least squares on the inferred order.
#
# Pairwise score: the asymmetric entropy-based likelihood ratio of Hyvarinen &
# Smith (2013), using Hyvarinen's differential-entropy approximation for
# standardized variables.  If x -> y, H(x) + H(r_{y|x}) < H(y) + H(r_{x|y}),
# so a NEGATIVE pw_score(x, y) favours "x causes y".
# -----------------------------------------------------------------------------
H_hat <- function(u) {
  # Hyvarinen differential entropy approximation for standardized u.
  # Returns H(u) up to a constant (constants cancel in pairwise comparisons).
  sdu <- sd(u); if (sdu <= 0) return(0)
  u <- u / sdu
  k1 <- 79.047; k2 <- 7.4129; gamma_lc <- 0.37457
  -k1 * (mean(log(cosh(u))) - gamma_lc)^2 - k2 * mean(u * exp(-u^2 / 2))^2
}

pw_score <- function(x, y) {
  x <- x - mean(x); y <- y - mean(y)
  sx2 <- sum(x * x); sy2 <- sum(y * y)
  if (sx2 <= 0 || sy2 <= 0) return(0)
  beta_yx <- sum(x * y) / sx2               # y = beta*x + r_{y|x}
  beta_xy <- sum(x * y) / sy2               # x = beta*y + r_{x|y}
  r_yx <- y - beta_yx * x
  r_xy <- x - beta_xy * y
  # (H(x) + H(r_yx)) - (H(y) + H(r_xy)) ; negative => x -> y
  (H_hat(x) + H_hat(r_yx)) - (H_hat(y) + H_hat(r_xy))
}

direct_lingam <- function(X) {
  p <- ncol(X)
  Xres <- scale(X, center = TRUE, scale = FALSE)   # mean-centred, scale unchanged
  remaining <- seq_len(p)
  order <- integer(0)
  while (length(remaining) > 1) {
    # T(j) = sum over i != j of max(0, pw_score(x_j, x_i))^2.
    # Small T(j) means x_j is "most exogenous among remaining".
    Tj <- sapply(remaining, function(j) {
      others <- setdiff(remaining, j)
      sum(vapply(others, function(i) {
        s <- pw_score(Xres[, j], Xres[, i])
        max(0, s)^2
      }, numeric(1)))
    })
    j_star <- remaining[which.min(Tj)]
    order <- c(order, j_star)
    # Residualize remaining on x_{j_star}
    vj <- Xres[, j_star]
    vj2 <- sum(vj * vj)
    if (vj2 > 0) {
      for (k in setdiff(remaining, j_star)) {
        b <- sum(Xres[, k] * vj) / vj2
        Xres[, k] <- Xres[, k] - b * vj
      }
    }
    remaining <- setdiff(remaining, j_star)
  }
  order <- c(order, remaining)

  # Estimate B by OLS using the inferred causal order (all on original X).
  Xc <- scale(X, center = TRUE, scale = FALSE)
  B <- matrix(0, p, p)
  for (idx in 2:p) {
    j <- order[idx]
    parents <- order[1:(idx - 1)]
    Xp <- Xc[, parents, drop = FALSE]
    yj <- Xc[, j]
    # OLS with ridge fallback for numerical safety
    G <- crossprod(Xp)
    beta <- tryCatch(solve(G, crossprod(Xp, yj)),
                     error = function(e) solve(G + 1e-8 * diag(nrow(G)),
                                               crossprod(Xp, yj)))
    B[j, parents] <- as.numeric(beta)
  }
  list(B = B, order = order)
}

# Thin wrapper so the rest of the script just calls lingam_B(X)
lingam_B <- function(X) direct_lingam(X)$B

# -----------------------------------------------------------------------------
# One replication at sample size n
# -----------------------------------------------------------------------------
one_rep <- function(n) {
  d <- generate_data(n)
  keys <- sort(names(d))                         # x1, x2, x3, x4
  X <- do.call(cbind, d[keys])                   # n x 4 matrix

  # --- LiNGAM residuals -----------------------------------------------------
  B <- lingam_B(X)
  lingam_resid <- X - X %*% t(B)                 # (I - B) x

  # --- Proposed GAM residuals (leave-one-out) -------------------------------
  prop_resid <- matrix(NA_real_, nrow = n, ncol = length(keys),
                       dimnames = list(NULL, keys))
  for (i in seq_along(keys)) {
    y   <- d[[keys[i]]]
    Xmi <- X[, -i, drop = FALSE]
    df  <- as.data.frame(Xmi); colnames(df) <- paste0("v", 1:3); df$y <- y
    fit <- mgcv::gam(y ~ s(v1, k = 10) + s(v2, k = 10) + s(v3, k = 10),
                     data = df, method = "REML", gamma = 1.4)
    prop_resid[, i] <- y - predict(fit, newdata = df)
  }

  # --- Per-variable metrics -------------------------------------------------
  out <- list(fitted_lin = numeric(length(keys)),
              fitted_prop = numeric(length(keys)),
              resid_lin   = numeric(length(keys)),
              resid_prop  = numeric(length(keys)),
              corr        = numeric(length(keys)),
              ks_stat     = numeric(length(keys)),
              ks_p        = numeric(length(keys)),
              tt_stat     = numeric(length(keys)),
              tt_p        = numeric(length(keys)))
  for (i in seq_along(keys)) {
    xi <- d[[keys[i]]]
    l  <- lingam_resid[, i]
    p  <- prop_resid[, i]
    out$fitted_lin[i]  <- sqrt(mean((xi - l)^2))      # ||fitted||/sqrt(n)
    out$fitted_prop[i] <- sqrt(mean((xi - p)^2))
    out$resid_lin[i]   <- sqrt(mean(l^2))             # true residual RMSE
    out$resid_prop[i]  <- sqrt(mean(p^2))
    out$corr[i]        <- cor(l, p)
    ks                 <- suppressWarnings(ks.test(l, p))
    out$ks_stat[i]     <- unname(ks$statistic)
    out$ks_p[i]        <- ks$p.value
    tt                 <- t.test(l, p, paired = TRUE)
    out$tt_stat[i]     <- unname(tt$statistic)
    out$tt_p[i]        <- tt$p.value
  }
  out
}

# -----------------------------------------------------------------------------
# Main loop
# -----------------------------------------------------------------------------
set.seed(BASE_SEED)
agg <- data.frame(n                = SAMPLE_SIZES,
                  fitted_lin       = NA_real_, se_fitted_lin   = NA_real_,
                  fitted_prop      = NA_real_, se_fitted_prop  = NA_real_,
                  resid_lin        = NA_real_, se_resid_lin    = NA_real_,
                  resid_prop       = NA_real_, se_resid_prop   = NA_real_,
                  corr             = NA_real_, se_corr         = NA_real_,
                  ks_stat          = NA_real_, se_ks           = NA_real_,
                  ks_p             = NA_real_,
                  tt_stat          = NA_real_, se_tt           = NA_real_,
                  tt_p             = NA_real_)

for (row in seq_along(SAMPLE_SIZES)) {
  n <- SAMPLE_SIZES[row]
  cat(sprintf("n = %4d: ", n))
  means_per_rep <- data.frame(fl = numeric(N_REP), fp = numeric(N_REP),
                              rl = numeric(N_REP), rp = numeric(N_REP),
                              co = numeric(N_REP), ks = numeric(N_REP),
                              kp = numeric(N_REP), tt = numeric(N_REP),
                              tp = numeric(N_REP))
  for (r in 1:N_REP) {
    set.seed(BASE_SEED + 1000L * n + r)           # matches Python per-rep seeding
    o <- one_rep(n)
    means_per_rep$fl[r] <- mean(o$fitted_lin)
    means_per_rep$fp[r] <- mean(o$fitted_prop)
    means_per_rep$rl[r] <- mean(o$resid_lin)
    means_per_rep$rp[r] <- mean(o$resid_prop)
    means_per_rep$co[r] <- mean(o$corr)
    means_per_rep$ks[r] <- mean(o$ks_stat)
    means_per_rep$kp[r] <- mean(o$ks_p)
    means_per_rep$tt[r] <- mean(o$tt_stat)
    means_per_rep$tp[r] <- mean(o$tt_p)
    cat(".")
  }
  cat(" done.\n")
  se <- function(v) sd(v) / sqrt(length(v))
  agg$fitted_lin[row]     <- mean(means_per_rep$fl);  agg$se_fitted_lin[row]  <- se(means_per_rep$fl)
  agg$fitted_prop[row]    <- mean(means_per_rep$fp);  agg$se_fitted_prop[row] <- se(means_per_rep$fp)
  agg$resid_lin[row]      <- mean(means_per_rep$rl);  agg$se_resid_lin[row]   <- se(means_per_rep$rl)
  agg$resid_prop[row]     <- mean(means_per_rep$rp);  agg$se_resid_prop[row]  <- se(means_per_rep$rp)
  agg$corr[row]           <- mean(means_per_rep$co);  agg$se_corr[row]        <- se(means_per_rep$co)
  agg$ks_stat[row]        <- mean(means_per_rep$ks);  agg$se_ks[row]          <- se(means_per_rep$ks)
  agg$ks_p[row]           <- mean(means_per_rep$kp)
  agg$tt_stat[row]        <- mean(means_per_rep$tt);  agg$se_tt[row]          <- se(means_per_rep$tt)
  agg$tt_p[row]           <- mean(means_per_rep$tp)
}

cat("\n=== Aggregated results (mean over", N_REP, "replications) ===\n")
print(round(agg, 4), row.names = FALSE)

# Focused side-by-side view of Panel 1 so the reader can verify at a glance
# whether the LiNGAM and Proposed fitted-component norms actually differ.
cat("\n--- Panel 1 focus: fitted-component norm ||x - e_hat||/sqrt(n) ---\n")
panel1 <- data.frame(n        = agg$n,
                     LiNGAM   = round(agg$fitted_lin,  4),
                     Proposed = round(agg$fitted_prop, 4),
                     diff     = round(agg$fitted_lin - agg$fitted_prop, 4))
print(panel1, row.names = FALSE)
write.csv(agg, file.path(OUT_DIR_RES, "linear_diag_R.csv"), row.names = FALSE)
cat(sprintf("\nSaved: %s\n", file.path(OUT_DIR_RES, "linear_diag_R.csv")))

# -----------------------------------------------------------------------------
# Figure: four-panel aggregate diagnostics (base R, 2 x 2 layout)
# -----------------------------------------------------------------------------
# Helper: draw all 4 panels on whatever device is currently active.
draw_panels <- function(agg) {
  par(mfrow = c(2, 2), mar = c(4.5, 4.8, 3.2, 1.2), cex.lab = 1.1,
      cex.main = 1.15, mgp = c(2.8, 0.8, 0))

  # Panel 1: fitted-component norm (both methods with error bars).
  # We pad ylim generously so the horizontal legend at the top has room, and
  # we draw Proposed first and LiNGAM last so the blue curve is never hidden
  # under the green one when the two values nearly coincide.  A small x-jitter
  # keeps the markers from exactly overlapping.
  all_y   <- c(agg$fitted_lin  + agg$se_fitted_lin,
               agg$fitted_lin  - agg$se_fitted_lin,
               agg$fitted_prop + agg$se_fitted_prop,
               agg$fitted_prop - agg$se_fitted_prop)
  y_lo <- min(all_y, na.rm = TRUE) - 0.35    # bottom padding for legend
  y_hi <- max(all_y, na.rm = TRUE) + 0.10
  x_rng <- diff(range(agg$n))
  dx    <- 0.008 * x_rng                     # ~ +/- 1% of x-range
  x_lin <- agg$n + dx
  x_pro <- agg$n - dx

  # frame first (empty), then grid, then curves on top -----------------------
  plot(agg$n, agg$fitted_lin, type = "n",
       xlab = "Sample size n",
       ylab = expression("||x -" ~ hat(e) ~ "||"[2] / sqrt(n)),
       main = "Fitted-component norm", ylim = c(y_lo, y_hi))
  grid(col = "gray85", lty = 3)

  # Proposed (drawn first => ends up underneath)
  arrows(x_pro, agg$fitted_prop - agg$se_fitted_prop,
         x_pro, agg$fitted_prop + agg$se_fitted_prop,
         angle = 90, code = 3, length = 0.04, col = "#2ca02c", lwd = 1.5)
  lines(x_pro, agg$fitted_prop, type = "o", pch = 15, lty = 2, lwd = 2.2,
        col = "#2ca02c", cex = 1.3)

  # LiNGAM (drawn last => sits on top)
  arrows(x_lin, agg$fitted_lin - agg$se_fitted_lin,
         x_lin, agg$fitted_lin + agg$se_fitted_lin,
         angle = 90, code = 3, length = 0.04, col = "#1f77b4", lwd = 1.5)
  lines(x_lin, agg$fitted_lin, type = "o", pch = 21, bg = "white",
        lty = 1, lwd = 2.2, col = "#1f77b4", cex = 1.5)

  # One-line legend in the bottom-right corner, white background.
  # ncol = 2 keeps it on a single row without the column-padding that
  # horiz = TRUE introduces; seg.len/x.intersp shrink the preview-line and
  # inter-entry gaps so the legend is as compact as its content.
  legend("bottomright", legend = c("LiNGAM", "Proposed"),
         pch = c(21, 15), pt.bg = c("white", NA),
         lty = c(1, 2), lwd = 2.2,
         col = c("#1f77b4", "#2ca02c"),
         ncol = 2,
         seg.len = 1.6,     # default ~2; shorter line sample
         x.intersp = 0.6,   # default 1.0; tighter marker-to-text gap
         y.intersp = 0.9,
         bg = "white", box.col = "gray70",
         inset = c(0.03, 0.03), cex = 0.95)

  # ---- helper: pad ylim so the legend always has headroom ------------------
  pad_ylim <- function(y, se = NULL, bottom_pad = 0.0, top_pad = 0.15) {
    if (is.null(se)) { lo <- min(y, na.rm = TRUE); hi <- max(y, na.rm = TRUE) }
    else             { lo <- min(y - se, na.rm = TRUE); hi <- max(y + se, na.rm = TRUE) }
    rng <- hi - lo; if (rng == 0) rng <- abs(hi) + 1
    c(lo - bottom_pad * rng, hi + top_pad * rng)
  }

  # Panel 2: Pearson correlation between LiNGAM and Proposed residuals -------
  plot(agg$n, agg$corr, type = "o", pch = 19, lwd = 2, cex = 1.3,
       col = "#d62728",
       xlab = "Sample size n", ylab = "Mean Pearson correlation",
       main = "Correlation between residuals",
       ylim = pad_ylim(agg$corr, agg$se_corr, 0.20, 0.20))
  grid(col = "gray85", lty = 3)
  arrows(agg$n, agg$corr - agg$se_corr, agg$n, agg$corr + agg$se_corr,
         angle = 90, code = 3, length = 0.04, col = "#d62728", lwd = 1.5)
  legend("bottomright",
         legend = expression("Mean" %+-% "SE"),
         pch = 19, lty = 1, lwd = 2, col = "#d62728",
         bg = "white", box.col = "gray70",
         inset = c(0.03, 0.03), cex = 0.9)

  # Panel 3: Two-sample KS statistic between the two sets of residuals -------
  plot(agg$n, agg$ks_stat, type = "o", pch = 19, lwd = 2, cex = 1.3,
       col = "#9467bd",
       xlab = "Sample size n", ylab = "Mean KS statistic",
       main = "KS distance between residual distributions",
       ylim = pad_ylim(agg$ks_stat, agg$se_ks, 0.20, 0.20))
  grid(col = "gray85", lty = 3)
  arrows(agg$n, agg$ks_stat - agg$se_ks, agg$n, agg$ks_stat + agg$se_ks,
         angle = 90, code = 3, length = 0.04, col = "#9467bd", lwd = 1.5)
  legend("topright",
         legend = expression("Mean" %+-% "SE"),
         pch = 19, lty = 1, lwd = 2, col = "#9467bd",
         bg = "white", box.col = "gray70",
         inset = c(0.03, 0.03), cex = 0.9)

  # Panel 4: KS p-value (with alpha = 0.05 reference line) -------------------
  plot(agg$n, agg$ks_p, type = "o", pch = 19, lwd = 2, cex = 1.3,
       col = "#17becf",
       xlab = "Sample size n", ylab = "Mean KS p-value",
       main = "KS test p-value (H0: same distribution)",
       ylim = c(0, 1))
  grid(col = "gray85", lty = 3)
  abline(h = 0.05, col = "red", lty = 2, lwd = 2)
  legend("topright",
         legend = c("Mean p-value", expression(alpha == 0.05)),
         pch = c(19, NA), lty = c(1, 2), lwd = c(2, 2),
         col = c("#17becf", "red"),
         bg = "white", box.col = "gray70",
         inset = c(0.03, 0.03), cex = 0.9)

  # Panels 5-6 (paired-t statistic and p-value) were removed on purpose:
  # averaging t-stats and p-values across replications/variables is not a
  # valid summary, and the paired-t here is hyper-powered by large n so it
  # rejects H0 of equal means at microscopic offsets that KS (Panel 3-4)
  # already shows are distributionally irrelevant.  If a replication-level
  # summary is wanted in a future revision, replace with median |t| and the
  # rejection rate (fraction of reps with p < alpha), not the means.
}

# --- (1) Always save the PDF for archiving ----------------------------------
pdf_path <- file.path(OUT_DIR_FIG, "linear_diag_R.pdf")
pdf(pdf_path, width = 12, height = 9)
draw_panels(agg)
dev.off()
cat(sprintf("Saved: %s\n", pdf_path))

# --- (2) Also show it on-screen if running interactively (RStudio, R GUI) ---
if (interactive()) {
  # Use the platform's default screen device (RStudio's Plots pane, X11, quartz)
  dev.new(width = 12, height = 9, noRStudioGD = FALSE)
  draw_panels(agg)
  cat("Plot displayed in the active graphics device.\n")
} else {
  cat("Non-interactive session (Rscript): open",
      pdf_path, "to view the figure,\n",
      "or run this script inside RStudio via Source to see the plot on-screen.\n")
}

cat("\nAll done.\n")
