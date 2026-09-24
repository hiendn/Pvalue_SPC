# Uniform EWMA PDF and CDF plots.
#
# Recreates the two uniform-EWMA figures used in the manuscript.

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

pdf_tildeU <- function(u, lambda, t, u0 = 0.5) {
  if (t < 1) stop("t must be >= 1")
  u <- as.numeric(u)
  if (t == 1) {
    offset <- (1 - lambda) * u0
    return(stats::dunif(u, min = offset, max = offset + lambda))
  }
  a_ts <- lambda * (1 - lambda)^(t - seq_len(t))
  offset <- (1 - lambda)^t * u0
  denom <- factorial(t - 1) * prod(a_ts)
  out <- numeric(length(u))
  for (m in 0:(2^t - 1)) {
    bits <- as.integer(intToBits(m))[seq_len(t)]
    idx <- which(bits == 1L)
    a_S <- if (length(idx) == 0L) 0 else sum(a_ts[idx])
    sign <- if ((length(idx) %% 2L) == 0L) 1 else -1
    out <- out + sign * pmax(u - offset - a_S, 0)^(t - 1)
  }
  out / denom
}

cdf_tildeU <- function(u, lambda, t, u0 = 0.5) {
  if (t < 1) stop("t must be >= 1")
  u <- as.numeric(u)
  if (t == 1) {
    offset <- (1 - lambda) * u0
    return(stats::punif(u, min = offset, max = offset + lambda))
  }
  a_ts <- lambda * (1 - lambda)^(t - seq_len(t))
  offset <- (1 - lambda)^t * u0
  denom <- factorial(t) * prod(a_ts)
  out <- numeric(length(u))
  for (m in 0:(2^t - 1)) {
    bits <- as.integer(intToBits(m))[seq_len(t)]
    idx <- which(bits == 1L)
    a_S <- if (length(idx) == 0L) 0 else sum(a_ts[idx])
    sign <- if ((length(idx) %% 2L) == 0L) 1 else -1
    out <- out + sign * pmax(u - offset - a_S, 0)^t
  }
  out / denom
}

simulate_tildeU <- function(n, lambda, t, u0 = 0.5) {
  U <- matrix(stats::runif(n * t), nrow = n, ncol = t)
  a_ts <- lambda * (1 - lambda)^(t - seq_len(t))
  offset <- (1 - lambda)^t * u0
  as.numeric(offset + U %*% a_ts)
}

plot_pdf_grid <- function(outfile, device = c("pdf", "png"), seed = 123L,
                          n_sim = 1e5, lambdas = c(0.3, 0.5, 0.7),
                          ts = c(2L, 3L, 4L), u0 = 0.5) {
  device <- match.arg(device)
  if (device == "pdf") grDevices::pdf(outfile, width = 15 / 2.54, height = 15 / 2.54)
  if (device == "png") grDevices::png(outfile, width = 1800, height = 1800, res = 240)
  on.exit(grDevices::dev.off(), add = TRUE)
  set.seed(seed)
  old <- graphics::par(mfrow = c(3, 3), mar = c(4, 4, 2, 1))
  on.exit(graphics::par(old), add = TRUE)
  for (lambda in lambdas) {
    for (tt in ts) {
      sim <- simulate_tildeU(n_sim, lambda, tt, u0 = u0)
      s_min <- (1 - lambda)^tt * u0
      s_max <- 1 - (1 - lambda)^tt * (1 - u0)
      x_grid <- seq(s_min, s_max, length.out = 400)
      graphics::hist(sim, breaks = 60, freq = FALSE, xlim = c(s_min, s_max),
                     main = bquote(lambda == .(lambda) ~ "," ~ t == .(tt)),
                     xlab = expression(tilde(U)[lambda,t]), ylab = "Density",
                     col = "grey80", border = "grey40")
      graphics::lines(x_grid, pdf_tildeU(x_grid, lambda, tt, u0 = u0), lwd = 2)
    }
  }
}

plot_cdf_grid <- function(outfile, device = c("pdf", "png"),
                          lambdas = c(0.3, 0.5, 0.7), ts = c(2L, 3L, 4L),
                          u0 = 0.5) {
  device <- match.arg(device)
  if (device == "pdf") grDevices::pdf(outfile, width = 15 / 2.54, height = 15 / 2.54)
  if (device == "png") grDevices::png(outfile, width = 1800, height = 1800, res = 240)
  on.exit(grDevices::dev.off(), add = TRUE)
  old <- graphics::par(mfrow = c(3, 3), mar = c(4, 4, 2, 1))
  on.exit(graphics::par(old), add = TRUE)
  alphas <- seq(0, 0.5, length.out = 400)
  for (lambda in lambdas) {
    for (tt in ts) {
      F_vals <- cdf_tildeU(alphas, lambda, tt, u0 = u0)
      graphics::plot(alphas, F_vals, type = "l", xlab = expression(alpha),
                     ylab = expression(P(tilde(U)[lambda,t] <= alpha)),
                     xlim = c(0, 0.5), ylim = c(0, 0.5),
                     main = bquote(lambda == .(lambda) ~ "," ~ t == .(tt)))
      graphics::abline(0, 1, lty = 2)
    }
  }
}

make_uniform_ewma_figures <- function(figure_dir = "generated_figures") {
  root <- getwd()
  ensure_dir(file.path(root, figure_dir))
  plot_pdf_grid(file.path(root, figure_dir, "Unif_EWMA_Figs.pdf"), "pdf")
  plot_pdf_grid(file.path(root, figure_dir, "Unif_EWMA_Figs.png"), "png")
  plot_cdf_grid(file.path(root, figure_dir, "Super_Unif_CDF.pdf"), "pdf")
  plot_cdf_grid(file.path(root, figure_dir, "Super_Unif_CDF.png"), "png")
  message("Uniform-EWMA figures written to ", figure_dir, "/.")
}

if (sys.nframe() == 0L) make_uniform_ewma_figures()
