# Non-linear Causal Discovery via Additive Models and HSIC

An R implementation of a non-linear causal-discovery procedure for
continuous observational data. The method combines penalised additive
models, fitted by restricted maximum likelihood (REML), with the
Hilbert–Schmidt Independence Criterion (HSIC) in a permutation test for
residual dependence. Edges are selected by Benjamini–Hochberg FDR
control on the pairwise HSIC permutation p-values (optionally followed
by an effect-size filter for high-sample, low-dimensional data), and
orientation is determined by the asymmetry of the directional HSIC
scores. A post-processing step enforces acyclicity.

## Quick start

```bash
git clone https://github.com/ash141886/nonlinear-causal-inference-gam.git
cd nonlinear-causal-inference-gam
Rscript -e 'install.packages(c("mgcv","ggplot2","dplyr","tidyr","gridExtra","fastICA","MASS","igraph"))'
Rscript scripts/03b_wine_red_twostage.R
```

The quick-start command runs the two-stage procedure on the UCI Wine
Quality (red) dataset and writes the estimated causal graph, per-edge
HSIC scores, and permutation p-values to `results/` and `figures/`.

## Method

For every ordered pair of variables $(i, j)$, the procedure fits a
penalised additive model of $X_i$ on all remaining variables *except*
$X_j$,
$$
  X_i \;=\; \sum_{k \neq i, j} f_k(X_k) \;+\; r_i^{(-j)},
$$
using cubic B-splines with REML-selected smoothness penalties, and then
measures the residual dependence $\widehat{\mathrm{HSIC}}(r_i^{(-j)},
X_j)$ with a Gaussian kernel whose bandwidth is set by the median
heuristic. The directional asymmetry $M_{ij} - M_{ji}$ is used to
orient the edge. Significance is assessed by permuting the sample
indices of $X_j$ (with a matched permutation of $X_i$), yielding a
calibrated finite-sample null.

Two edge-selection modes are provided:

- **Single stage.** Benjamini–Hochberg control of the false discovery
  rate at level $\alpha$ on the pairwise p-values
  $T_{ij} = \max(M_{ij}, M_{ji})$. Appropriate for simulation studies
  and for datasets where the sample size is moderate relative to the
  dimension.
- **Two stage.** Liberal BH screening at $\alpha = 0.20$ followed by
  retention of the top-$k$ pairs ranked by full-sample HSIC effect
  size, where $k = \lceil p(p-1)\rho/2 \rceil$ for a target density
  $\rho$. Recommended when $n$ is large relative to $p$, so that pure
  FDR control tends to retain statistically significant but
  substantively weak edges.

## Repository layout

```
nonlinear-causal-inference-gam/
├── R/                             # Core library
│   ├── utils.R                    # Kernels, HSIC, metrics, I/O
│   ├── additive_hsic.R            # Main discovery procedure
│   └── lingam_baseline.R          # Linear-non-Gaussian baseline
├── scripts/                       # Reproducible analyses
│   ├── 01_simulation_nonlinear.R  # Non-linear benchmark grid
│   ├── 02_simulation_linear.R     # Linear-data residual diagnostics
│   ├── 03_wine_red.R              # Red-wine analysis (single-stage)
│   ├── 03b_wine_red_twostage.R    # Red-wine analysis (two-stage)
│   ├── 03c_wine_lingam_verify.R   # Cross-check against LiNGAM
│   ├── 04_wine_white.R            # White-wine analysis
│   ├── 05_plot_simulations.R      # Figures from 01_* results
│   ├── 06_plot_wine.R             # DAG and partial-effect figures
│   ├── 07_wine_concurvity.R       # Concurvity diagnostics (wine)
│   └── 08_wine_robustness.R       # Robustness/perturbation study
├── cluster/                       # Optional SLURM wrappers
├── data/                          # Bundled + auto-downloaded data
├── figures/                       # Generated figures (git-ignored)
├── results/                       # Cached results (git-ignored)
├── LICENSE
└── README.md
```

## Requirements

R (≥ 4.0) with the following packages:

```r
install.packages(c("mgcv", "ggplot2", "dplyr", "tidyr", "gridExtra",
                   "fastICA", "MASS", "igraph"))
```

No Python or cluster dependency is required: every script can be run
locally as a single `Rscript` invocation. The `cluster/` directory
contains optional SLURM wrappers for users with cluster access.

## Reproducing the analyses

Run scripts from the repository root:

```bash
# 1. Simulations
Rscript scripts/01_simulation_nonlinear.R   # ~4-8 h on a laptop
Rscript scripts/02_simulation_linear.R      # ~3-5 min

# 2. Wine analyses
Rscript scripts/03b_wine_red_twostage.R     # ~5-10 min
Rscript scripts/04_wine_white.R             # ~5-10 min

# 3. Figures
Rscript scripts/05_plot_simulations.R       # non-linear summary figures
Rscript scripts/06_plot_wine.R              # wine DAG + partial effects
```

The linear-data diagnostic figure is produced directly by
`scripts/02_simulation_linear.R`; it does not go through
`scripts/05_plot_simulations.R`. Wine data are downloaded automatically
from the UCI Machine Learning Repository on first run and cached in
`data/`.

## Key parameters

| Parameter | Simulation | Wine (red) | Role |
|---|---|---|---|
| `alpha` | 0.05 | 0.20 (screening) | BH false-discovery-rate level |
| `n_perm` | 200 | 1 000 | Permutation count for HSIC null |
| `n_sub` | – | 200 | Sub-sample size for permutation test |
| `target_density` | – | 0.15 | Target edge density in the two-stage filter |
| `gamma` | 1.4 | 1.4 | Extra smoothing penalty in `mgcv::gam` |
| Kernel | Gaussian RBF | Gaussian RBF | Median-heuristic bandwidth |
| Spline basis | cubic B | cubic B | `s(..., bs = "cr")`, REML smoothing |

## Data

The red-wine CSV (`data/winequality-red.csv`, $n = 1\,599$, $p = 12$)
is bundled with the repository for immediate reproducibility. The
white-wine file is downloaded on demand. Both datasets are from
[Cortez et al., 2009](https://archive.ics.uci.edu/ml/datasets/Wine+Quality).

## Citation

If you use this code, please cite:

```bibtex
@article{islam_suzuki_additive_hsic,
  author  = {Islam, Md. Ashraful and Suzuki, Joe},
  title   = {Non-linear Causal Inference in Observational Data using
             Additive Models and Kernel-based Independence Testing},
  journal = {Japanese Journal of Statistics and Data Science},
  year    = {2026},
  url     = {https://github.com/ash141886/nonlinear-causal-inference-gam}
}
```

## License

Released under the MIT License; see `LICENSE`.
