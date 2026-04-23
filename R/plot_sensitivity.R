# EWMA sensitivity heatmap.
#
# Plots log10(conservativeness ratio) for the \tilde Q chart at
# n0 = 50, alpha = 0.05, and k = 1.

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

plot_heatmap <- function(outfile, device = c("pdf", "png"),
                         csv = find_csv("ks_ic_ewma_sensitivity_fast.csv"),
                         n0_value = 50, alpha_value = 0.05, k_value = 1) {
  device <- match.arg(device)
  dat <- utils::read.csv(csv, stringsAsFactors = FALSE, check.names = FALSE)
  dat <- dat[dat$chart_type == "qtilde" & dat$n0 == n0_value &
               abs(dat$alpha_nominal - alpha_value) < 1e-12 & dat$k == k_value, , drop = FALSE]
  if (nrow(dat) == 0L) stop("No matching qtilde rows found in ", csv)
  lambdas <- sort(unique(dat$lambda))
  rs <- sort(unique(dat$r))
  z <- matrix(NA_real_, nrow = length(lambdas), ncol = length(rs),
              dimnames = list(sprintf("%.2f", lambdas), sprintf("%.2f", rs)))
  ratio <- matrix(NA_real_, nrow = length(lambdas), ncol = length(rs), dimnames = dimnames(z))
  for (i in seq_len(nrow(dat))) {
    ii <- match(dat$lambda[i], lambdas)
    jj <- match(dat$r[i], rs)
    z[ii, jj] <- log10(dat$ratio[i])
    ratio[ii, jj] <- dat$ratio[i]
  }
  rng <- range(z, na.rm = TRUE)
  if (device == "pdf") grDevices::pdf(outfile, width = 8.3, height = 5.3)
  if (device == "png") grDevices::png(outfile, width = 1890, height = 1205, res = 220)
  on.exit(grDevices::dev.off(), add = TRUE)
  old <- graphics::par(mar = c(4.5, 4.5, 4.7, 5.7), xpd = FALSE)
  on.exit(graphics::par(old), add = TRUE)
  graphics::image(x = seq_along(lambdas), y = seq_along(rs), z = z,
                  axes = FALSE, xlab = expression(lambda), ylab = expression(r),
                  main = bquote(log[10] ~ " conservativeness ratio for " ~ tilde(Q)[lambda,t]^{(r)} * "\n" *
                                  n[0] == .(n0_value) * ", " * alpha == .(alpha_value) * ", " * k == .(k_value)),
                  col = grDevices::hcl.colors(64, "viridis"))
  graphics::axis(1, at = seq_along(lambdas), labels = sprintf("%.2f", lambdas))
  graphics::axis(2, at = seq_along(rs), labels = sprintf("%.2f", rs), las = 1)
  graphics::box()
  for (ii in seq_along(lambdas)) for (jj in seq_along(rs)) {
    if (is.finite(ratio[ii, jj])) graphics::text(ii, jj, labels = format(round(ratio[ii, jj]), big.mark = ","), cex = 0.85)
  }
  # Simple colour legend.
  usr <- graphics::par("usr")
  x0 <- usr[2] + 0.28; x1 <- usr[2] + 0.43
  yseq <- seq(usr[3], usr[4], length.out = 65)
  graphics::rect(x0, yseq[-length(yseq)], x1, yseq[-1], col = grDevices::hcl.colors(64, "viridis"), border = NA, xpd = TRUE)
  graphics::rect(x0, usr[3], x1, usr[4], border = "black", xpd = TRUE)
  zticks <- pretty(rng, n = 5)
  ypos <- usr[3] + (zticks - rng[1]) / diff(rng) * diff(usr[3:4])
  graphics::axis(4, at = ypos, labels = sprintf("%.2f", zticks), las = 1, line = 4.6)
  graphics::mtext(expression(log[10](hat(L)[1] / B[1])), side = 4, line = 3.1)
}

make_ewma_sensitivity_heatmap <- function(figure_dir = "figures") {
  ensure_dir(figure_dir)
  plot_heatmap(file.path(figure_dir, "ewma_qtilde_sensitivity_heatmap.pdf"), "pdf")
  plot_heatmap(file.path(figure_dir, "ewma_qtilde_sensitivity_heatmap.png"), "png")
  message("EWMA sensitivity heatmap written to ", figure_dir, "/.")
}

if (sys.nframe() == 0L) make_ewma_sensitivity_heatmap()
