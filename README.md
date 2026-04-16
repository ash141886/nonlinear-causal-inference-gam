# Non-linear Causal Inference using Additive Models and Kernel-based Independence Testing

R implementation of the additive-HSIC method for non-linear causal discovery from observational data.

**Paper:** Islam, M.A. and Suzuki, J. "Non-linear Causal Inference in Observational Data using Additive Models and Kernel-based Independence Testing." *Japanese Journal of Statistics and Data Science* (under review).

## Quick Start

```bash
git clone https://github.com/ash141886/nonlinear-causal-inference-gam.git
cd nonlinear-causal-inference-gam
Rscript -e 'install.packages(c("mgcv","ggplot2","dplyr","tidyr","gridExtra","fastICA","MASS","igraph"))'
Rscript scripts/03b_wine_red_twostage.R
```

## Method

The method pairs penalized additive models (GAMs via `mgcv`) with the Hilbert-Schmidt Independence Criterion (HSIC) in a leave-one-out testing procedure. For each potential causal link, a reduced additive model is fitted excluding one variable, and the remaining dependence between the prediction errors and the excluded variable is measured with HSIC to determine edge direction.

Edge selection uses either:
- **Single-stage:** Benjamini-Hochberg FDR control at level α (default for simulations)
- **Two-stage:** BH-FDR screening at a liberal level (α = 0.20) followed by HSIC effect-size filtering (recommended for real data with large n relative to p)

## Repository Structure

```
nonlinear-causal-inference-gam/
├── R/                          # Core functions
│   ├── utils.R                 # Kernels, HSIC, metrics, data loading
│   ├── additive_hsic.R         # Main method (Algorithm 1)
│   └── lingam_baseline.R       # LiNGAM comparison
├── scripts/                    # Reproducible analysis scripts
│   ├── 01_simulation_nonlinear.R   # Section 4.1: non-linear benchmarks (Figs 1-6)
│   ├── 02_simulation_linear.R      # Section 4.2: linear-data diagnostics (Fig 7)
│   ├── 03_wine_red.R               # Section 5: red wine analysis
│   ├── 03b_wine_red_twostage.R     # Two-stage variant for red wine
│   ├── 03c_wine_lingam_verify.R    # LiNGAM verification on wine
│   ├── 04_wine_white.R             # White wine analysis
│   ├── 05_plot_simulations.R       # Renders Figures 1-6 from 01_* output
│   └── 06_plot_wine.R              # Figures 8-9
├── data/                       # Data files (auto-downloaded)
├── figures/                    # Generated figures
├── results/                    # Saved results (.rds)
└── README.md
```

## Requirements

R (>= 4.0) with packages:

```r
install.packages(c("mgcv", "ggplot2", "dplyr", "tidyr", "gridExtra",
                    "fastICA", "MASS", "igraph"))
```

## Reproducing the Paper

Run scripts in order from the repository root:

```bash
# 1. Simulations (Section 4) — ~4-8 hours
Rscript scripts/01_simulation_nonlinear.R
Rscript scripts/02_simulation_linear.R

# 2. Wine analysis (Section 5) — ~5-10 minutes each
Rscript scripts/03b_wine_red_twostage.R
Rscript scripts/04_wine_white.R

# 3. Generate figures
# Figures 1-6 (non-linear simulation summaries):
Rscript scripts/05_plot_simulations.R
# Figure 7 is produced directly by scripts/02_simulation_linear.R above.
# Figures 8-9 (wine):
Rscript scripts/06_plot_wine.R
```

Wine data is automatically downloaded from the [UCI Machine Learning Repository](https://archive.ics.uci.edu/ml/datasets/Wine+Quality) on first run.

## Key Parameters

| Parameter | Simulation | Wine (Red) | Description |
|-----------|-----------|------------|-------------|
| `alpha` | 0.05 | 0.20 (screening) | FDR level |
| `n_perm` | 200 | 1,000 | Permutation count |
| `n_sub` | — | 200 | Subsample size for permutation test |
| `target_density` | — | 0.15 | Target edge density for top-k filter |
| `gamma` | 1.4 | 1.4 | Extra smoothing penalty in mgcv |

## Citation

```bibtex
@article{islam2026nonlinear,
  title={Non-linear Causal Inference in Observational Data using Additive
         Models and Kernel-based Independence Testing},
  author={Islam, Md. Ashraful and Suzuki, Joe},
  journal={Japanese Journal of Statistics and Data Science},
  year={2026},
  note={Under review},
  url={https://github.com/ash141886/nonlinear-causal-inference-gam}
}
```

## License

MIT
