# Directional and coordinate localisation simulations.
#
# Run from the repository root with:
#   source("R/sim_localisation.R")
#
# The normal example uses Z tests. The two-phase Cauchy example uses exact
# Mann--Whitney p-values with cached null distributions.

source("R/utils.R")
ensure_dir("output")

# -----------------------------------------------------------------------------
# User-adjustable parameters
# -----------------------------------------------------------------------------

N_REP_RUN <- 1000L
N_REP_FWER <- 10000L
MAX_T_RUN <- 200000L
BATCH_SIZE_RUN <- 50L
BATCH_SIZE_FWER <- 1000L
VERBOSE <- TRUE

DELTA_GRID <- c(0.5, 1)
RHO_GRID <- c(0, 0.5, 0.9)
ALPHA_GRID <- c(0.01, 0.05)
N0_CAUCHY_GRID <- c(20L, 50L, 100L)
D <- 3L

# -----------------------------------------------------------------------------
# Shared localisation helpers
# -----------------------------------------------------------------------------

make_sigma <- function(d = 3L, rho = 0) {
  S <- matrix(rho, nrow = d, ncol = d)
  diag(S) <- 1
  S
}

rmvnorm_chol <- function(n, mean, cholSigma) {
  d <- length(mean)
  Z <- matrix(stats::rnorm(n * d), nrow = n, ncol = d) %*% cholSigma
  sweep(Z, 2, mean, `+`)
}

rmvcauchy_chol <- function(n, location, cholSigma) {
  d <- length(location)
  Z <- matrix(stats::rnorm(n * d), nrow = n, ncol = d) %*% cholSigma
  U <- stats::rchisq(n, df = 1)
  X <- sweep(Z, 1, sqrt(U), "/")
  sweep(X, 2, location, `+`)
}

localise_from_p <- function(p_up, p_down, alpha, agg_type = c("bonferroni", "mean")) {
  agg_type <- match.arg(agg_type)
  d <- length(p_up)
  p_coord <- pmin(1, 2 * pmin(p_up, p_down))
  p_global <- if (agg_type == "bonferroni") {
    min(1, d * min(p_coord))
  } else {
    min(1, min(2, d) * mean(p_coord))
  }

  # Holm step-down, explicit implementation to avoid p.adjust overhead in tight loops.
  ord <- order(p_coord)
  J <- integer(0)
  for (k in seq_len(d)) {
    if (p_coord[ord[k]] <= alpha / (d - k + 1L)) {
      J <- c(J, ord[k])
    } else {
      break
    }
  }
  direction <- if (length(J) > 0L) ifelse(p_up[J] <= p_down[J], "up", "down") else character(0)
  list(p_coord = p_coord, p_global = p_global, J = J, direction = direction)
}

count_true_directions <- function(J, direction, true_direction) {
  if (length(J) == 0L) return(0L)
  sum(mapply(function(j, dir) identical(true_direction[j], dir), J, direction))
}

count_false_directions <- function(J, direction, true_direction) {
  if (length(J) == 0L) return(0L)
  sum(mapply(function(j, dir) !identical(true_direction[j], dir), J, direction))
}

summarise_localisation <- function(run_df, fwer_vals, model, n0 = NA_integer_, delta, rho, alpha) {
  data.frame(
    model = model,
    n0 = n0,
    delta = delta,
    rho = rho,
    alpha = alpha,
    n_rep_run = nrow(run_df),
    mean_R = mean(run_df$R),
    se_R = se(run_df$R),
    censored_rate = mean(run_df$censored),
    mean_true_ooc = mean(run_df$true_ooc),
    se_true_ooc = se(run_df$true_ooc),
    mean_false_ooc_at_alarm = mean(run_df$false_ooc),
    n_rep_fwer = length(fwer_vals),
    fwer_hat = mean(fwer_vals),
    fwer_se = se(fwer_vals),
    stringsAsFactors = FALSE
  )
}

# -----------------------------------------------------------------------------
# Normal one-phase localisation
# -----------------------------------------------------------------------------

normal_pvals <- function(x) {
  z <- as.numeric(x)
  list(p_up = 1 - stats::pnorm(z), p_down = stats::pnorm(z))
}

sim_one_normal_localisation <- function(delta, rho, alpha, max_t = MAX_T_RUN,
                                        agg_type = "bonferroni") {
  Sigma <- make_sigma(D, rho)
  cholSigma <- chol(Sigma)
  mu <- c(delta, rep(0, D - 2L), -delta)
  true_direction <- c("up", rep("none", D - 2L), "down")

  for (t in seq_len(max_t)) {
    x <- rmvnorm_chol(1L, mu, cholSigma)[1, ]
    pv <- normal_pvals(x)
    loc <- localise_from_p(pv$p_up, pv$p_down, alpha, agg_type = agg_type)
    if (loc$p_global <= alpha) {
      return(data.frame(R = t,
                        true_ooc = count_true_directions(loc$J, loc$direction, true_direction),
                        false_ooc = count_false_directions(loc$J, loc$direction, true_direction),
                        censored = FALSE))
    }
  }
  data.frame(R = max_t + 1L, true_ooc = 0L, false_ooc = 0L, censored = TRUE)
}

fwer_one_normal <- function(delta, rho, alpha, agg_type = "bonferroni") {
  # Fixed-time directional FWER under the same alternative used in the
  # localisation power study.  An error is counted only after the global
  # p-value chart would alarm at this fixed time.
  Sigma <- make_sigma(D, rho)
  cholSigma <- chol(Sigma)
  mu <- c(delta, rep(0, D - 2L), -delta)
  true_direction <- c("up", rep("none", D - 2L), "down")
  x <- rmvnorm_chol(1L, mu, cholSigma)[1, ]
  pv <- normal_pvals(x)
  loc <- localise_from_p(pv$p_up, pv$p_down, alpha, agg_type = agg_type)
  if (loc$p_global > alpha) return(0L)
  as.integer(count_false_directions(loc$J, loc$direction, true_direction) > 0L)
}


run_normal_config <- function(delta, rho, alpha, seed, checkpoint_file = NULL,
                              verbose = TRUE) {
  log_msg("Normal localisation config: delta=%.3g rho=%.3g alpha=%.3g", delta, rho, alpha, verbose = verbose)
  run_rows <- vector("list", N_REP_RUN)
  t0 <- Sys.time()
  for (i in seq_len(N_REP_RUN)) {
    set.seed(seed + i)
    run_rows[[i]] <- sim_one_normal_localisation(delta, rho, alpha)
    run_rows[[i]]$replicate <- i
    if (verbose && (i == 1L || i %% BATCH_SIZE_RUN == 0L || i == N_REP_RUN)) {
      part <- do.call(rbind, run_rows[seq_len(i)])
      log_msg("Normal run reps %d/%d: mean_R=%.3f mean_true=%.3f elapsed=%s",
              i, N_REP_RUN, mean(part$R), mean(part$true_ooc),
              format_elapsed(as.numeric(difftime(Sys.time(), t0, units = "secs"))), verbose = verbose)
      if (!is.null(checkpoint_file)) utils::write.csv(part, checkpoint_file, row.names = FALSE)
    }
  }
  run_df <- do.call(rbind, run_rows)

  fwer_vals <- integer(N_REP_FWER)
  for (i in seq_len(N_REP_FWER)) {
    set.seed(seed + 10000000L + i)
    fwer_vals[i] <- fwer_one_normal(delta, rho, alpha)
    if (verbose && (i == 1L || i %% BATCH_SIZE_FWER == 0L || i == N_REP_FWER)) {
      log_msg("Normal FWER reps %d/%d: fwer_hat=%.4f", i, N_REP_FWER, mean(fwer_vals[seq_len(i)]), verbose = verbose)
    }
  }

  summarise_localisation(run_df, fwer_vals, model = "normal_Z", delta = delta, rho = rho, alpha = alpha)
}

normal_grid <- expand.grid(delta = DELTA_GRID, rho = RHO_GRID, alpha = ALPHA_GRID, stringsAsFactors = FALSE)
normal_rows <- vector("list", nrow(normal_grid))
for (i in seq_len(nrow(normal_grid))) {
  g <- normal_grid[i, ]
  normal_rows[[i]] <- run_normal_config(g$delta, g$rho, g$alpha,
                                        seed = 80000L + i * 100000L,
                                        checkpoint_file = NULL,
                                        verbose = VERBOSE)
}
normal_out <- do.call(rbind, normal_rows)
utils::write.csv(normal_out, "output/localisation_normal_fast.csv", row.names = FALSE)
log_msg("Saved output/localisation_normal_fast.csv", verbose = VERBOSE)

# -----------------------------------------------------------------------------
# Cauchy two-phase localisation using fast exact Mann--Whitney p-values
# -----------------------------------------------------------------------------

precompute_wilcox_null <- function(n0, n1) {
  u_vals <- 0:(n0 * n1)
  pmf <- stats::dwilcox(u_vals, m = n1, n = n0)
  cdf <- cumsum(pmf)
  list(n0 = n0, n1 = n1, cdf = cdf)
}

make_wilcox_null_cache <- function(n0_vals, width = 10L, verbose = TRUE) {
  cache <- list()
  for (n0 in n0_vals) {
    n1_candidates <- (n0 - width):(n0 + width)
    n1_candidates <- unique(n1_candidates[n1_candidates > 0L])
    for (n1 in n1_candidates) {
      key <- paste0(n0, "_", n1)
      if (is.null(cache[[key]])) {
        cache[[key]] <- precompute_wilcox_null(n0, n1)
      }
    }
  }
  log_msg("Precomputed %d Wilcoxon null distributions.", length(cache), verbose = verbose)
  cache
}

wilcox_U_stat <- function(x0, xt) {
  n1 <- length(xt)
  r <- rank(c(xt, x0))
  W <- sum(r[seq_len(n1)])
  U <- W - n1 * (n1 + 1) / 2
  as.integer(round(U))
}

mw_onesided_exact_fast <- function(x0, xt, null) {
  n0 <- null$n0
  n1 <- null$n1
  U <- wilcox_U_stat(x0, xt)
  U <- max(0L, min(U, n0 * n1))
  cdf <- null$cdf
  p_down <- cdf[U + 1L]                # Phase-II shifted downward: P(U <= Uobs)
  p_up <- if (U == 0L) 1.0 else 1.0 - cdf[U]  # Phase-II shifted upward: P(U >= Uobs)
  c(p_up = p_up, p_down = p_down)
}

mw_pvals_fast <- function(x0, xt, null_cache) {
  d <- ncol(x0)
  n0 <- nrow(x0)
  nt <- nrow(xt)
  key <- paste0(n0, "_", nt)
  null <- null_cache[[key]]
  if (is.null(null)) stop("Missing Wilcoxon null cache for key ", key)
  p_up <- numeric(d)
  p_down <- numeric(d)
  for (j in seq_len(d)) {
    tmp <- mw_onesided_exact_fast(x0[, j], xt[, j], null)
    p_up[j] <- tmp["p_up"]
    p_down[j] <- tmp["p_down"]
  }
  list(p_up = p_up, p_down = p_down)
}

sim_one_cauchy_localisation <- function(n0, delta, rho, alpha, null_cache,
                                        max_t = MAX_T_RUN, agg_type = "bonferroni") {
  Sigma <- make_sigma(D, rho)
  cholSigma <- chol(Sigma)
  mu0 <- rep(0, D)
  mu <- c(delta, rep(0, D - 2L), -delta)
  true_direction <- c("up", rep("none", D - 2L), "down")

  x0 <- rmvcauchy_chol(n0, mu0, cholSigma)
  for (t in seq_len(max_t)) {
    nt <- sample_vss(n0, 10L)
    xt <- rmvcauchy_chol(nt, mu, cholSigma)
    pv <- mw_pvals_fast(x0, xt, null_cache)
    loc <- localise_from_p(pv$p_up, pv$p_down, alpha, agg_type = agg_type)
    if (loc$p_global <= alpha) {
      return(data.frame(R = t,
                        true_ooc = count_true_directions(loc$J, loc$direction, true_direction),
                        false_ooc = count_false_directions(loc$J, loc$direction, true_direction),
                        censored = FALSE))
    }
  }
  data.frame(R = max_t + 1L, true_ooc = 0L, false_ooc = 0L, censored = TRUE)
}

fwer_one_cauchy <- function(n0, delta, rho, alpha, null_cache, agg_type = "bonferroni") {
  # Fixed-time directional FWER under the same Cauchy alternative used in the
  # localisation power study.  An error is counted only after the global
  # p-value chart would alarm at this fixed time.
  Sigma <- make_sigma(D, rho)
  cholSigma <- chol(Sigma)
  mu0 <- rep(0, D)
  mu <- c(delta, rep(0, D - 2L), -delta)
  true_direction <- c("up", rep("none", D - 2L), "down")
  x0 <- rmvcauchy_chol(n0, mu0, cholSigma)
  nt <- sample_vss(n0, 10L)
  xt <- rmvcauchy_chol(nt, mu, cholSigma)
  pv <- mw_pvals_fast(x0, xt, null_cache)
  loc <- localise_from_p(pv$p_up, pv$p_down, alpha, agg_type = agg_type)
  if (loc$p_global > alpha) return(0L)
  as.integer(count_false_directions(loc$J, loc$direction, true_direction) > 0L)
}


run_cauchy_config <- function(n0, delta, rho, alpha, null_cache, seed,
                              checkpoint_file = NULL, verbose = TRUE) {
  log_msg("Cauchy localisation config: n0=%d delta=%.3g rho=%.3g alpha=%.3g", n0, delta, rho, alpha, verbose = verbose)
  run_rows <- vector("list", N_REP_RUN)
  t0 <- Sys.time()
  for (i in seq_len(N_REP_RUN)) {
    set.seed(seed + i)
    run_rows[[i]] <- sim_one_cauchy_localisation(n0, delta, rho, alpha, null_cache)
    run_rows[[i]]$replicate <- i
    if (verbose && (i == 1L || i %% BATCH_SIZE_RUN == 0L || i == N_REP_RUN)) {
      part <- do.call(rbind, run_rows[seq_len(i)])
      log_msg("Cauchy run reps %d/%d: mean_R=%.3f mean_true=%.3f elapsed=%s",
              i, N_REP_RUN, mean(part$R), mean(part$true_ooc),
              format_elapsed(as.numeric(difftime(Sys.time(), t0, units = "secs"))), verbose = verbose)
      if (!is.null(checkpoint_file)) utils::write.csv(part, checkpoint_file, row.names = FALSE)
    }
  }
  run_df <- do.call(rbind, run_rows)

  fwer_vals <- integer(N_REP_FWER)
  for (i in seq_len(N_REP_FWER)) {
    set.seed(seed + 10000000L + i)
    fwer_vals[i] <- fwer_one_cauchy(n0, delta, rho, alpha, null_cache)
    if (verbose && (i == 1L || i %% BATCH_SIZE_FWER == 0L || i == N_REP_FWER)) {
      log_msg("Cauchy FWER reps %d/%d: fwer_hat=%.4f", i, N_REP_FWER, mean(fwer_vals[seq_len(i)]), verbose = verbose)
    }
  }

  summarise_localisation(run_df, fwer_vals, model = "cauchy_MW_exact_fast",
                         n0 = n0, delta = delta, rho = rho, alpha = alpha)
}

null_cache <- make_wilcox_null_cache(N0_CAUCHY_GRID, width = 10L, verbose = VERBOSE)
cauchy_grid <- expand.grid(n0 = N0_CAUCHY_GRID, delta = DELTA_GRID,
                           rho = RHO_GRID, alpha = ALPHA_GRID,
                           stringsAsFactors = FALSE)
cauchy_rows <- vector("list", nrow(cauchy_grid))
for (i in seq_len(nrow(cauchy_grid))) {
  g <- cauchy_grid[i, ]
  cauchy_rows[[i]] <- run_cauchy_config(g$n0, g$delta, g$rho, g$alpha, null_cache,
                                        seed = 200000L + i * 100000L,
                                        checkpoint_file = NULL,
                                        verbose = VERBOSE)
}

cauchy_out <- do.call(rbind, cauchy_rows)
utils::write.csv(cauchy_out, "output/localisation_cauchy_fast.csv", row.names = FALSE)
log_msg("Saved output/localisation_cauchy_fast.csv", verbose = VERBOSE)
log_msg("Finished localisation simulations.", verbose = VERBOSE)
