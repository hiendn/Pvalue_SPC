# IC calibration and EWMA sensitivity simulations.
#
# Run from the repository root with:
#   source("R/sim_ic.R")
#
# The script writes the elementary IC examples, the AR(1) example, the raw KS
# chart study, and the KS EWMA sensitivity study to output/.

source("R/utils.R")
ensure_dir("output")

# -----------------------------------------------------------------------------
# User-adjustable parameters
# -----------------------------------------------------------------------------

# Use "asymp" for fast screening and sensitivity analysis.  For final exact
# finite-sample KS p-values, try "psmirnov_exact".  If psmirnov is unavailable,
# that engine falls back to stats::ks.test() and will be slower.
KS_ENGINE <- "asymp"
JMAX <- 50L
EXACT_MAX_PRODUCT <- 15000L
FINITE_SAMPLE_CORRECTION <- FALSE

N_CORES <- max(1L, parallel::detectCores(logical = FALSE) - 1L)
BATCH_SIZE <- 25L
VERBOSE <- TRUE

N_REP_ELEMENTARY <- 1000L
N_REP_AR1 <- 1000L
N_REP_KS_RAW <- 1000L
N_REP_EWMA_GRID <- 500L

MAX_T_AR1 <- 200000L
MAX_T_KS <- 200000L

# KS grids.
KS_RAW_N0 <- c(20L, 50L, 100L)
EWMA_N0 <- c(50L, 100L)
ALPHA_RAW <- c(0.01, 0.05)
ALPHA_EWMA <- c(0.05, 0.10)
K_VALUES <- c(1L, 5L)

# Sensitivity grid.  The multi-chart simulator computes all charts on the same
# KS p-value stream.
QTILDE_LAMBDAS <- c(0.5, 0.7, 0.9)
QTILDE_RS <- c(-0.9, -0.8, -0.5, 0.5)
QBAR_LAMBDAS <- c(0.8, 0.9, 0.95)
QBAR_RS <- c(1, 2)

# -----------------------------------------------------------------------------
# Elementary two-phase normal example: negative-binomial speedup
# -----------------------------------------------------------------------------

karl_two_phase_nb_fast <- function(alpha = 0.01, k_values = c(1L, 5L), reps = 1000L,
                                   mu0 = 0, mu1 = 0, seed = 1L,
                                   verbose = TRUE, batch_size = 100L) {
  stopifnot(alpha > 0, alpha < 1, reps >= 1L)
  set.seed(seed)
  zthr <- stats::qnorm(1 - alpha / 2)
  k_values <- sort(unique(as.integer(k_values)))
  rows <- vector("list", length(k_values))
  Tk <- matrix(NA_real_, nrow = reps, ncol = length(k_values))
  colnames(Tk) <- paste0("R", k_values)
  t0 <- Sys.time()

  for (i in seq_len(reps)) {
    x0 <- stats::rnorm(1, mean = mu0, sd = 1)
    m <- (mu1 - x0) / sqrt(2)
    s <- 1 / sqrt(2)
    q <- stats::pnorm(-zthr, mean = m, sd = s) + (1 - stats::pnorm(zthr, mean = m, sd = s))
    q <- min(max(q, .Machine$double.xmin), 1)
    for (j in seq_along(k_values)) {
      kk <- k_values[j]
      Tk[i, j] <- stats::rnbinom(1, size = kk, prob = q) + kk
    }
    if (verbose && (i == 1L || i %% batch_size == 0L || i == reps)) {
      log_msg("Normal two-phase NB: alpha=%.3g rep=%d/%d elapsed=%s",
              alpha, i, reps,
              format_elapsed(as.numeric(difftime(Sys.time(), t0, units = "secs"))),
              verbose = verbose)
    }
  }

  for (j in seq_along(k_values)) {
    kk <- k_values[j]
    vals <- Tk[, j]
    ci <- mc_ci(vals)
    lower <- bound_prop3_marginal(alpha, kk)
    rows[[j]] <- data.frame(
      chart = "two-phase normal p-value",
      alpha = alpha,
      validity = "marginal",
      k = kk,
      n_rep = reps,
      mean_Rk = mean(vals),
      sd_Rk = stats::sd(vals),
      se_Rk = se(vals),
      ci_lower = ci[1],
      ci_upper = ci[2],
      mcse_ratio = se(vals) / mean(vals),
      censored_rate = 0,
      lower_bound = lower,
      ratio = mean(vals) / lower,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

log_msg("Running elementary two-phase normal simulations using negative-binomial speedup.", verbose = VERBOSE)
normal_out <- do.call(rbind, lapply(seq_along(ALPHA_RAW), function(i) {
  karl_two_phase_nb_fast(alpha = ALPHA_RAW[i], k_values = K_VALUES,
                         reps = N_REP_ELEMENTARY, seed = 11000L + 1000L * i,
                         verbose = VERBOSE, batch_size = max(1L, floor(N_REP_ELEMENTARY / 10L)))
}))
utils::write.csv(normal_out, "output/elementary_normal_ic_fast.csv", row.names = FALSE)
log_msg("Saved output/elementary_normal_ic_fast.csv", verbose = VERBOSE)

# -----------------------------------------------------------------------------
# AR(1) example: both Pprime and Pstar in one stream
# -----------------------------------------------------------------------------

p_marginal_ar1 <- function(x_t) 2 * (1 - stats::pnorm(abs(x_t)))

p_star_ar1 <- function(x_t, x_tm1) {
  s <- sqrt(max(x_t^2 - x_tm1^2, 0))
  2 * (1 - stats::pnorm(s))
}

simulate_one_ar1_both <- function(beta, alpha, k_values = c(1L, 5L),
                                  delta = 0, max_t = 200000L) {
  sigma_eps <- sqrt(1 - beta^2)
  x_prev <- stats::rnorm(1)
  k_values <- sort(unique(as.integer(k_values)))
  k_max <- max(k_values)

  counts <- c(Pprime = 0L, Pstar = 0L)
  R <- matrix(NA_integer_, nrow = 2L, ncol = length(k_values),
              dimnames = list(c("Pprime", "Pstar"), paste0("R", k_values)))
  finished <- c(Pprime = FALSE, Pstar = FALSE)

  for (t in seq_len(max_t)) {
    x_t <- delta + beta * x_prev + stats::rnorm(1, 0, sigma_eps)
    pm <- p_marginal_ar1(x_t)
    ps <- p_star_ar1(x_t, x_prev)

    if (!finished["Pprime"] && pm <= alpha) {
      counts["Pprime"] <- counts["Pprime"] + 1L
      hit <- which(k_values == counts["Pprime"])
      if (length(hit) == 1L) R["Pprime", hit] <- t
      if (counts["Pprime"] >= k_max) finished["Pprime"] <- TRUE
    }
    if (!finished["Pstar"] && ps <= alpha) {
      counts["Pstar"] <- counts["Pstar"] + 1L
      hit <- which(k_values == counts["Pstar"])
      if (length(hit) == 1L) R["Pstar", hit] <- t
      if (counts["Pstar"] >= k_max) finished["Pstar"] <- TRUE
    }
    if (all(finished)) break
    x_prev <- x_t
  }

  rows <- list()
  idx <- 1L
  for (method in rownames(R)) {
    for (j in seq_along(k_values)) {
      cens <- is.na(R[method, j])
      rows[[idx]] <- data.frame(
        method = method,
        k = k_values[j],
        Rk = if (cens) max_t + 1L else R[method, j],
        censored = cens,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }
  do.call(rbind, rows)
}

run_ar1_grid_fast <- function(betas = c(0.1, 0.5), alphas = c(0.01, 0.05),
                              k_values = c(1L, 5L), reps = 1000L,
                              max_t = 200000L, seed = 1L,
                              verbose = TRUE, batch_size = 100L) {
  grid <- expand.grid(beta = betas, alpha = alphas, stringsAsFactors = FALSE)
  out_rows <- list()
  out_idx <- 1L
  for (gidx in seq_len(nrow(grid))) {
    g <- grid[gidx, ]
    log_msg("AR(1) grid start: beta=%.2f alpha=%.3g reps=%d", g$beta, g$alpha, reps, verbose = verbose)
    ans <- vector("list", reps)
    t0 <- Sys.time()
    for (i in seq_len(reps)) {
      set.seed(seed + 100000L * gidx + i)
      ans[[i]] <- simulate_one_ar1_both(beta = g$beta, alpha = g$alpha,
                                        k_values = k_values, max_t = max_t)
      ans[[i]]$replicate <- i
      if (verbose && (i == 1L || i %% batch_size == 0L || i == reps)) {
        part <- do.call(rbind, ans[seq_len(i)])
        sm <- aggregate(Rk ~ method + k, data = part, FUN = mean)
        log_msg("AR(1): beta=%.2f alpha=%.3g rep=%d/%d elapsed=%s",
                g$beta, g$alpha, i, reps,
                format_elapsed(as.numeric(difftime(Sys.time(), t0, units = "secs"))),
                verbose = verbose)
        print(sm, row.names = FALSE)
      }
    }
    df <- do.call(rbind, ans)
    for (method in unique(df$method)) {
      for (kk in k_values) {
        z <- df[df$method == method & df$k == kk, , drop = FALSE]
        validity <- if (method == "Pstar") "conditional" else "marginal"
        lower <- bound_for_chart(g$alpha, kk, validity)
        ci <- mc_ci(z$Rk)
        out_rows[[out_idx]] <- data.frame(
          chart = method,
          beta = g$beta,
          alpha = g$alpha,
          validity = validity,
          k = kk,
          n_rep = nrow(z),
          mean_Rk = mean(z$Rk),
          sd_Rk = stats::sd(z$Rk),
          se_Rk = se(z$Rk),
          ci_lower = ci[1],
          ci_upper = ci[2],
          mcse_ratio = se(z$Rk) / mean(z$Rk),
          censored_rate = mean(z$censored),
          lower_bound = lower,
          ratio = mean(z$Rk) / lower,
          stringsAsFactors = FALSE
        )
        out_idx <- out_idx + 1L
      }
    }
  }
  do.call(rbind, out_rows)
}

log_msg("Running AR(1) simulations.", verbose = VERBOSE)
ar1_out <- run_ar1_grid_fast(reps = N_REP_AR1, max_t = MAX_T_AR1, seed = 12000L,
                             verbose = VERBOSE, batch_size = max(1L, floor(N_REP_AR1 / 10L)))
utils::write.csv(ar1_out, "output/ar1_ic_fast.csv", row.names = FALSE)
log_msg("Saved output/ar1_ic_fast.csv", verbose = VERBOSE)

# -----------------------------------------------------------------------------
# Raw KS p-value chart under VSS
# -----------------------------------------------------------------------------

log_msg("Running raw KS IC simulations with engine=%s.", KS_ENGINE, verbose = VERBOSE)
ks_raw_rows <- list()
idx <- 1L
for (n0 in KS_RAW_N0) {
  for (alpha in ALPHA_RAW) {
    charts <- list(make_p_chart(alpha, validity = "marginal",
                                label = sprintf("KS raw P_t alpha=%.3g", alpha)))
    cp <- sprintf("output/ks_raw_n0%d_alpha%s_raw_long.csv", n0, gsub("\\.", "p", as.character(alpha)))
    sm <- estimate_arl_ks_multi(
      charts = charts, n_rep = N_REP_KS_RAW, seed = 13000L + idx * 1000L,
      n0 = n0, max_t = MAX_T_KS, k_values = K_VALUES, n_cores = N_CORES,
      phase2_name = "IC", ks_engine = KS_ENGINE,
      exact_max_product = EXACT_MAX_PRODUCT, jmax = JMAX,
      finite_sample_correction = FINITE_SAMPLE_CORRECTION,
      batch_size = BATCH_SIZE, verbose = VERBOSE, checkpoint_file = NULL
    )
    sm$n0 <- n0
    sm$ks_engine <- KS_ENGINE
    ks_raw_rows[[idx]] <- sm
    idx <- idx + 1L
  }
}
ks_raw_out <- do.call(rbind, ks_raw_rows)
utils::write.csv(ks_raw_out, "output/ks_ic_raw_fast.csv", row.names = FALSE)
log_msg("Saved output/ks_ic_raw_fast.csv", verbose = VERBOSE)

# -----------------------------------------------------------------------------
# EWMA-like p-value chart sensitivity grid, computed multi-chart per path
# -----------------------------------------------------------------------------

make_sensitivity_charts <- function(alpha) {
  charts <- list(make_p_chart(alpha, validity = "marginal",
                              label = sprintf("raw P_t alpha=%.3g", alpha)))
  for (lam in QTILDE_LAMBDAS) {
    for (rr in QTILDE_RS) {
      charts[[length(charts) + 1L]] <- make_qtilde_chart(alpha, lambda = lam, r = rr,
                                                         validity = "marginal")
    }
  }
  for (lam in QBAR_LAMBDAS) {
    for (rr in QBAR_RS) {
      charts[[length(charts) + 1L]] <- make_qbar_chart(alpha, lambda = lam, r = rr,
                                                       validity = "marginal")
    }
  }
  charts
}

log_msg("Running KS EWMA sensitivity simulations with multi-chart stream evaluation.", verbose = VERBOSE)
ewma_summary_rows <- list()
ewma_raw_rows <- list()
idx <- 1L
for (n0 in EWMA_N0) {
  for (alpha in ALPHA_EWMA) {
    charts <- make_sensitivity_charts(alpha)
    cp <- sprintf("output/ks_ewma_sensitivity_n0%d_alpha%s_raw_long.csv",
                  n0, gsub("\\.", "p", as.character(alpha)))
    raw_df <- run_replicates_ks_multi(
      n_rep = N_REP_EWMA_GRID, seed = 14000L + idx * 1000L,
      charts = charts, n0 = n0, n_sampler = function() sample_vss(n0, 10L),
      phase1_gen = phase1_normal, phase2_gen = make_phase2_generator("IC"),
      k_values = K_VALUES, max_t = MAX_T_KS,
      ks_engine = KS_ENGINE, exact_max_product = EXACT_MAX_PRODUCT, jmax = JMAX,
      finite_sample_correction = FINITE_SAMPLE_CORRECTION,
      n_cores = N_CORES, batch_size = BATCH_SIZE,
      verbose = VERBOSE, checkpoint_file = NULL
    )
    raw_df$n0 <- n0
    raw_df$alpha <- alpha
    raw_df$ks_engine <- KS_ENGINE
    sm <- summarise_run_df(raw_df)
    sm$n0 <- n0
    sm$alpha_nominal <- alpha
    sm$ks_engine <- KS_ENGINE
    ewma_raw_rows[[idx]] <- raw_df
    ewma_summary_rows[[idx]] <- sm
    idx <- idx + 1L
  }
}

ewma_raw_out <- do.call(rbind, ewma_raw_rows)
ewma_out <- do.call(rbind, ewma_summary_rows)
utils::write.csv(ewma_raw_out, "output/ks_ic_ewma_sensitivity_fast_raw_long.csv", row.names = FALSE)
utils::write.csv(ewma_out, "output/ks_ic_ewma_sensitivity_fast.csv", row.names = FALSE)
log_msg("Saved EWMA sensitivity outputs.", verbose = VERBOSE)

# Deterministic multiplier table used to explain conservativeness in the paper.
multiplier_out <- data.frame(
  chart = c("Qtilde r=-0.9", "Qtilde r=-0.8", "Qtilde r=-0.5",
            "Qbar lambda=0.90 r=1", "Qbar lambda=0.95 r=1"),
  lambda = c(0.5, 0.5, 0.5, 0.90, 0.95),
  r = c(-0.9, -0.8, -0.5, 1, 1),
  multiplier = c(ewma_multiplier_qtilde(0.5, -0.9),
                 ewma_multiplier_qtilde(0.5, -0.8),
                 ewma_multiplier_qtilde(0.5, -0.5),
                 ewma_multiplier_qbar(0.90, 1),
                 ewma_multiplier_qbar(0.95, 1)),
  effective_alpha_0.05 = c(0.05 / ewma_multiplier_qtilde(0.5, -0.9),
                           0.05 / ewma_multiplier_qtilde(0.5, -0.8),
                           0.05 / ewma_multiplier_qtilde(0.5, -0.5),
                           0.05 / ewma_multiplier_qbar(0.90, 1),
                           0.05 / ewma_multiplier_qbar(0.95, 1)),
  effective_alpha_0.01 = c(0.01 / ewma_multiplier_qtilde(0.5, -0.9),
                           0.01 / ewma_multiplier_qtilde(0.5, -0.8),
                           0.01 / ewma_multiplier_qtilde(0.5, -0.5),
                           0.01 / ewma_multiplier_qbar(0.90, 1),
                           0.01 / ewma_multiplier_qbar(0.95, 1)),
  stringsAsFactors = FALSE
)
utils::write.csv(multiplier_out, "output/ewma_multiplier_table.csv", row.names = FALSE)
log_msg("Finished fast IC calibration and sensitivity simulations.", verbose = VERBOSE)
