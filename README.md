# P-value SPC manuscript code

R code and simulation outputs for the numerical studies in **Statistical process control via p-values**.

The repository contains the processed CSV files used in the manuscript, the scripts used to regenerate the figures, and the main simulation scripts.  The code uses base R only.

## Layout

```text
R/
  utils.R                 shared simulation utilities
  sim_ic.R                IC ARL and EWMA sensitivity studies
  sim_ks_benchmark.R      matched-IC-ARL KS benchmark study
  sim_localisation.R      normal and Cauchy localisation studies
  plot_uniform_ewma.R     uniform-EWMA PDF/CDF figures
  plot_sensitivity.R      EWMA conservativeness heatmap
  plot_benchmark.R        matched-ARL delay-ratio boxplot
  make_figures.R          regenerates all figures
  tables.R                compact CSV summaries for benchmark tables

scripts/                  short wrappers for common runs
data/                     processed simulation outputs used in the paper
figures/                  manuscript figures in PDF and PNG formats
output/                   default destination for rerun simulations
tables/                   generated table summaries
```

## Figures

From the repository root:

```r
source("scripts/make_figures.R")
```

This recreates the four figures in `figures/`:

```text
Unif_EWMA_Figs.pdf
Super_Unif_CDF.pdf
ewma_qtilde_sensitivity_heatmap.pdf
ks_benchmark_relative_delay_boxplot.pdf
```

PNG copies are written at the same time.

## Table summaries

```r
source("scripts/make_tables.R")
```

The summaries are written to `tables/`.

## Simulations

The processed outputs in `data/` are the results used in the manuscript. To rerun the simulations, use:

```r
source("scripts/run_ic_study.R")
source("scripts/run_ks_benchmark_smoke_test.R")  # quick check
source("scripts/run_ks_benchmark.R")             # full benchmark
source("scripts/run_localisation.R")
```

New simulation results are written to `output/`. The plotting and table scripts look in `output/` first when a freshly generated file is present; otherwise they use the packaged files in `data/`.

The KS benchmark is the longest run. The smoke test uses a small grid and is useful for checking that R and the local file paths are set up correctly.
