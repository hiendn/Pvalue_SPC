# Regenerate the uniform-EWMA PDF and CDF figures used in the manuscript.

source("R/plot_uniform_ewma.R")

make_all_figures <- function(figure_dir = "generated_figures") {
  make_uniform_ewma_figures(figure_dir)
  invisible(TRUE)
}

if (sys.nframe() == 0L) make_all_figures()
