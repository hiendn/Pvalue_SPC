# Matched-IC-ARL benchmark delay plot.
#
# Reports OOC delay ratios relative to the raw p-value chart.

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

find_csv <- function(filename) {
  candidates <- c(file.path("output", filename), file.path("data", filename), filename)
  hit <- candidates[file.exists(candidates)]
  if (!length(hit)) stop("Could not find ", filename)
  hit[[1L]]
}

make_ratio_data <- function(csv = find_csv("ks_benchmark_ooc_fast.csv")) {
  ooc <- utils::read.csv(csv, stringsAsFactors = FALSE, check.names = FALSE)
  persistent <- c("N_0.5_1", "N_1_1", "N_0_2", "Cauchy")
  ooc$scenario <- ifelse(ooc$ooc %in% persistent, "persistent", "dynamic")
  key_cols <- c("n0", "target_arl", "ooc", "k")
  raw <- ooc[ooc$chart == "p-value Shewhart", c(key_cols, "mean_Rk"), drop = FALSE]
  names(raw)[names(raw) == "mean_Rk"] <- "raw_mean_Rk"
  m <- merge(ooc[ooc$chart != "p-value Shewhart", , drop = FALSE], raw, by = key_cols, all.x = TRUE)
  m$ratio <- m$mean_Rk / m$raw_mean_Rk
  m <- m[is.finite(m$ratio) & m$ratio > 0, , drop = FALSE]
  m$panel <- paste0("R", m$k, ", ", m$scenario)
  m
}

plot_ratio_boxplot <- function(outfile, device = c("pdf", "png")) {
  device <- match.arg(device)
  m <- make_ratio_data()
  charts <- c("KS Shewhart statistic", "KS-EWMA statistic", "KS-CUSUM statistic",
              "Qbar lambda=0.95", "Qtilde r=-0.8", "Qtilde r=-0.9")
  panels <- c("R1, persistent", "R1, dynamic", "R5, persistent", "R5, dynamic")
  box_data <- list(); at <- numeric(); panel_centres <- numeric(); pos <- 1
  for (pp in panels) {
    start <- pos
    for (cc in charts) {
      vals <- m$ratio[m$panel == pp & m$chart == cc]
      box_data[[length(box_data) + 1L]] <- vals
      at <- c(at, pos)
      pos <- pos + 1
    }
    panel_centres <- c(panel_centres, mean(start:(pos - 1)))
    pos <- pos + 1.2
  }
  if (device == "pdf") grDevices::pdf(outfile, width = 11.5, height = 3.2)
  if (device == "png") grDevices::png(outfile, width = 1800, height = 500, res = 160)
  on.exit(grDevices::dev.off(), add = TRUE)
  old <- graphics::par(mar = c(8.8, 4.5, 2.7, 1.2))
  on.exit(graphics::par(old), add = TRUE)
  graphics::boxplot(box_data, at = at, log = "y", xaxt = "n", outline = FALSE,
                    xlab = "", ylab = expression("OOC delay relative to raw " * P[t] * " chart"),
                    main = "Matched-IC-ARL OOC delays relative to raw p-value chart")
  graphics::abline(h = 1, lty = 2)
  lab <- rep(c("KS-Shewhart", "KS-EWMA", "KS-CUSUM", "Qbar(0.95)", "Qtilde(-0.8)", "Qtilde(-0.9)"), length(panels))
  graphics::axis(1, at = at, labels = FALSE)
  graphics::text(x = at, y = graphics::par("usr")[3], labels = lab, srt = 65, adj = 1, xpd = TRUE, cex = 0.78)
  graphics::mtext(c(expression(R[1] * ", persistent"), expression(R[1] * ", dynamic"),
                   expression(R[5] * ", persistent"), expression(R[5] * ", dynamic")),
                  side = 1, at = panel_centres, line = 5.9)
}

make_ks_benchmark_delay_plot <- function(figure_dir = "figures") {
  ensure_dir(figure_dir)
  plot_ratio_boxplot(file.path(figure_dir, "ks_benchmark_relative_delay_boxplot.pdf"), "pdf")
  plot_ratio_boxplot(file.path(figure_dir, "ks_benchmark_relative_delay_boxplot.png"), "png")
  message("KS benchmark relative-delay figure written to ", figure_dir, "/.")
}

if (sys.nframe() == 0L) make_ks_benchmark_delay_plot()
