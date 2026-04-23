# Regenerate all figures used in the manuscript.

source("R/plot_uniform_ewma.R")
source("R/plot_sensitivity.R")
source("R/plot_benchmark.R")

make_all_figures <- function(figure_dir = "figures") {
  make_uniform_ewma_figures(figure_dir)
  make_ewma_sensitivity_heatmap(figure_dir)
  make_ks_benchmark_delay_plot(figure_dir)
  invisible(TRUE)
}

if (sys.nframe() == 0L) make_all_figures()
