# Shared functions for the p-value SPC simulation studies.
#
# The KS routines compute the two-sample statistic directly from sorted samples.
# This keeps the two-phase KS simulations practical while retaining the option
# to use exact or hybrid p-value calculations for smaller sample sizes.

`%||%` <- function(x, y) if (is.null(x)) y else x

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

.now_string <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

log_msg <- function(..., verbose = TRUE) {
  if (isTRUE(verbose)) message(sprintf("[%s] %s", .now_string(), sprintf(...)))
  invisible(NULL)
}

format_elapsed <- function(seconds) {
  seconds <- as.numeric(seconds)
  if (!is.finite(seconds)) return("NA")
  if (seconds < 60) return(sprintf("%.1fs", seconds))
  if (seconds < 3600) return(sprintf("%.1fmin", seconds / 60))
  sprintf("%.2fh", seconds / 3600)
}

se <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) <= 1L) return(NA_real_)
  stats::sd(x) / sqrt(length(x))
}

mc_ci <- function(x, level = 0.95) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(c(NA_real_, NA_real_))
  z <- stats::qnorm(1 - (1 - level) / 2)
  m <- mean(x)
  s <- se(x)
  c(m - z * s, m + z * s)
}

bound_prop3_marginal <- function(alpha, k) {
  nu <- floor(k / alpha)
  (nu + 1) * (1 - (alpha * nu) / (2 * k))
}

bound_prop4_conditional <- function(alpha, k) k / alpha

bound_for_chart <- function(alpha, k, validity) {
  if (is.null(alpha) || is.na(alpha) || is.null(validity) || is.na(validity)) return(NA_real_)
  if (validity == "conditional") return(bound_prop4_conditional(alpha, k))
  if (validity == "marginal") return(bound_prop3_marginal(alpha, k))
  NA_real_
}

sample_vss <- function(n0, width = 10L) {
  as.integer(sample.int(2L * width + 1L, size = 1L) + n0 - width - 1L)
}

safe_p <- function(p) {
  p <- as.numeric(p)
  p[!is.finite(p)] <- 1
  pmin(pmax(p, .Machine$double.xmin), 1)
}

clamp01 <- function(x) {
  x <- as.numeric(x)
  x[!is.finite(x)] <- 1
  pmin(pmax(x, 0), 1)
}

# -----------------------------------------------------------------------------
# Fast two-sample Kolmogorov--Smirnov statistic and p-value engines
# -----------------------------------------------------------------------------

# Fast vectorised KS statistic for continuous data.  Ties are unlikely in all
# simulations in the manuscript.  For discrete/tied data, use ks_two_sample_D_tie_safe().
ks_two_sample_D_order <- function(x0_sorted, xt) {
  y_sorted <- sort(xt)
  n <- length(x0_sorted)
  m <- length(y_sorted)
  v <- c(x0_sorted, y_sorted)
  g <- c(rep.int(1L, n), rep.int(2L, m))
  o <- order(v)
  g <- g[o]
  c1 <- cumsum(g == 1L)
  c2 <- cumsum(g == 2L)
  max(abs(c1 / n - c2 / m))
}

# Tie-safe merge implementation.  It is more robust but can be slower in pure R.
ks_two_sample_D_tie_safe <- function(x0_sorted, xt) {
  y <- sort(xt)
  n <- length(x0_sorted)
  m <- length(y)
  i <- 1L
  j <- 1L
  c1 <- 0L
  c2 <- 0L
  D <- 0.0
  while (i <= n || j <= m) {
    if (j > m || (i <= n && x0_sorted[i] <= y[j])) {
      z <- x0_sorted[i]
    } else {
      z <- y[j]
    }
    while (i <= n && x0_sorted[i] <= z) {
      c1 <- c1 + 1L
      i <- i + 1L
    }
    while (j <= m && y[j] <= z) {
      c2 <- c2 + 1L
      j <- j + 1L
    }
    d <- abs(c1 / n - c2 / m)
    if (d > D) D <- d
  }
  D
}

kolmogorov_sf_series <- function(x, jmax = 50L) {
  x <- as.numeric(x)
  out <- numeric(length(x))
  jj <- seq_len(jmax)
  signs <- (-1)^(jj - 1L)
  for (i in seq_along(x)) {
    if (!is.finite(x[i]) || x[i] <= 0) {
      out[i] <- 1
    } else {
      terms <- signs * exp(-2 * (jj^2) * (x[i]^2))
      out[i] <- 2 * sum(terms)
    }
  }
  clamp01(out)
}

ks_p_asymp_from_D <- function(D, n0, nt, jmax = 50L, finite_sample_correction = FALSE) {
  ne <- (n0 * nt) / (n0 + nt)
  z <- sqrt(ne) * D
  if (isTRUE(finite_sample_correction)) {
    # The Stephens-type correction is commonly used as a finite-sample adjustment.
    # Set finite_sample_correction=FALSE to reproduce the simpler asymptotic
    # engine used in the earlier fast scripts.
    z <- (sqrt(ne) + 0.12 + 0.11 / sqrt(ne)) * D
  }
  kolmogorov_sf_series(z, jmax = jmax)
}

.psmirnov_fun <- function() {
  f <- get0("psmirnov", envir = asNamespace("stats"), mode = "function", inherits = FALSE)
  if (is.null(f)) f <- get0("psmirnov", envir = as.environment("package:stats"), mode = "function", inherits = FALSE)
  f
}

has_psmirnov <- function() !is.null(.psmirnov_fun())

ks_p_psmirnov_from_D <- function(D, n0, nt, exact = TRUE) {
  f <- .psmirnov_fun()
  if (is.null(f)) stop("stats::psmirnov is not available in this R installation")

  attempts <- list(
    function() f(q = D, sizes = c(n0, nt), alternative = "two.sided", exact = exact, lower.tail = FALSE),
    function() f(D, sizes = c(n0, nt), alternative = "two.sided", exact = exact, lower.tail = FALSE),
    function() f(q = D, sizes = c(n0, nt), exact = exact, lower.tail = FALSE),
    function() f(D, c(n0, nt), exact = exact, lower.tail = FALSE)
  )

  last_error <- NULL
  for (a in attempts) {
    val <- tryCatch(a(), error = function(e) { last_error <<- e; NA_real_ })
    if (length(val) == 1L && is.finite(val)) return(clamp01(val))
  }
  stop("All psmirnov() call signatures failed. Last error: ", conditionMessage(last_error))
}

ks_p_stats_test <- function(x0_sorted, xt, exact = TRUE) {
  kt <- suppressWarnings(stats::ks.test(x0_sorted, xt, exact = exact))
  p <- kt$p.value
  D <- unname(as.numeric(kt$statistic))
  list(D = D, p = clamp01(p))
}

# Main KS engine.  Recommended engines:
#   "asymp"          : fastest; approximate p-values, suitable for screening.
#   "hybrid"         : exact via psmirnov/ks.test for small n0*nt, asymptotic otherwise.
#   "psmirnov_exact" : exact if stats::psmirnov is available; otherwise falls back to ks.test.
#   "stats_exact"    : slow but close to the original implementation.
ks_p_and_d_fast <- function(x0_sorted, xt,
                            engine = c("asymp", "asymp_tie_safe", "hybrid",
                                       "psmirnov_exact", "psmirnov_asymp",
                                       "stats_exact", "stats_asymp"),
                            exact_max_product = 15000L,
                            jmax = 50L,
                            finite_sample_correction = FALSE) {
  engine <- match.arg(engine)
  n0 <- length(x0_sorted)
  nt <- length(xt)

  if (engine == "stats_exact") return(ks_p_stats_test(x0_sorted, xt, exact = TRUE))
  if (engine == "stats_asymp") return(ks_p_stats_test(x0_sorted, xt, exact = FALSE))

  D <- if (engine == "asymp_tie_safe") {
    ks_two_sample_D_tie_safe(x0_sorted, xt)
  } else {
    ks_two_sample_D_order(x0_sorted, xt)
  }

  if (engine == "asymp" || engine == "asymp_tie_safe") {
    return(list(D = D, p = ks_p_asymp_from_D(D, n0, nt, jmax, finite_sample_correction)))
  }

  if (engine == "psmirnov_asymp") {
    p <- tryCatch(ks_p_psmirnov_from_D(D, n0, nt, exact = FALSE),
                  error = function(e) ks_p_asymp_from_D(D, n0, nt, jmax, finite_sample_correction))
    return(list(D = D, p = clamp01(p)))
  }

  if (engine == "psmirnov_exact") {
    p <- tryCatch(ks_p_psmirnov_from_D(D, n0, nt, exact = TRUE),
                  error = function(e) ks_p_stats_test(x0_sorted, xt, exact = TRUE)$p)
    return(list(D = D, p = clamp01(p)))
  }

  if (engine == "hybrid") {
    if (n0 * nt <= exact_max_product) {
      p <- tryCatch(ks_p_psmirnov_from_D(D, n0, nt, exact = TRUE),
                    error = function(e) ks_p_stats_test(x0_sorted, xt, exact = TRUE)$p)
    } else {
      p <- tryCatch(ks_p_psmirnov_from_D(D, n0, nt, exact = FALSE),
                    error = function(e) ks_p_asymp_from_D(D, n0, nt, jmax, finite_sample_correction))
    }
    return(list(D = D, p = clamp01(p)))
  }

  stop("Unknown KS engine: ", engine)
}

# -----------------------------------------------------------------------------
# Chart definitions and updates
# -----------------------------------------------------------------------------

ewma_multiplier_qtilde <- function(lambda, r) {
  stopifnot(lambda > 0, lambda < 1, r > -1, r != 0)
  if (r >= 1) {
    min(1 + r, 1 / lambda)^(1 / r)
  } else {
    (1 + r)^(1 / r)
  }
}

ewma_multiplier_qbar <- function(lambda, r) {
  stopifnot(lambda > 0, lambda < 1, r >= 1)
  lambda^(-1 / r)
}

make_chart_label <- function(chart) {
  if (!is.null(chart$label)) return(chart$label)
  switch(chart$type,
         p_shewhart = sprintf("P_t <= %.4g", chart$threshold),
         qtilde = sprintf("Qtilde(lambda=%.3g,r=%.3g,alpha=%.4g)", chart$lambda, chart$r, chart$threshold),
         qbar = sprintf("Qbar(lambda=%.3g,r=%.3g,alpha=%.4g)", chart$lambda, chart$r, chart$threshold),
         ks_shewhart = sprintf("KS Shewhart(h=%.4g)", chart$threshold),
         ks_ewma = sprintf("KS-EWMA(lambda=%.3g,h=%.4g)", chart$lambda, chart$threshold),
         ks_cusum = sprintf("KS-CUSUM(kappa=%.4g,h=%.4g)", chart$kappa, chart$threshold),
         chart$type)
}

chart_param_frame <- function(charts) {
  rows <- lapply(seq_along(charts), function(i) {
    ch <- charts[[i]]
    data.frame(
      chart_index = i,
      chart = make_chart_label(ch),
      chart_type = ch$type %||% NA_character_,
      threshold = ch$threshold %||% NA_real_,
      validity = ch$validity %||% NA_character_,
      lambda = ch$lambda %||% NA_real_,
      r = ch$r %||% NA_real_,
      kappa = ch$kappa %||% NA_real_,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

update_chart <- function(state, chart, p_value, ks_stat) {
  p_value <- safe_p(p_value)
  type <- chart$type

  if (type == "p_shewhart") {
    stat <- p_value
    alarm <- stat <= chart$threshold
  } else if (type == "qtilde") {
    r <- chart$r
    lambda <- chart$lambda
    prev <- state$S %||% NULL
    S <- if (is.null(prev)) p_value^r else lambda * p_value^r + (1 - lambda) * prev
    state$S <- S
    stat <- ewma_multiplier_qtilde(lambda, r) * S^(1 / r)
    alarm <- stat <= chart$threshold
  } else if (type == "qbar") {
    r <- chart$r
    lambda <- chart$lambda
    prev <- state$S %||% NULL
    S <- if (is.null(prev)) p_value^r else lambda * p_value^r + (1 - lambda) * prev
    state$S <- S
    stat <- ewma_multiplier_qbar(lambda, r) * S^(1 / r)
    alarm <- stat <= chart$threshold
  } else if (type == "ks_shewhart") {
    stat <- ks_stat
    alarm <- stat >= chart$threshold
  } else if (type == "ks_ewma") {
    lambda <- chart$lambda
    prev <- state$W %||% NULL
    W <- if (is.null(prev)) ks_stat else lambda * ks_stat + (1 - lambda) * prev
    state$W <- W
    stat <- W
    alarm <- stat >= chart$threshold
  } else if (type == "ks_cusum") {
    Cprev <- state$C %||% 0
    C <- max(0, Cprev + ks_stat - chart$kappa)
    state$C <- C
    stat <- C
    alarm <- stat >= chart$threshold
  } else {
    stop("Unknown chart type: ", type)
  }

  list(state = state, statistic = stat, alarm = isTRUE(alarm))
}

# -----------------------------------------------------------------------------
# Phase-I and Phase-II data generators
# -----------------------------------------------------------------------------

phase1_normal <- function(n) stats::rnorm(n, mean = 0, sd = 1)

make_phase2_generator <- function(name) {
  force(name)
  function(n, t) {
    switch(name,
           IC = stats::rnorm(n, 0, 1),
           N_0.5_1 = stats::rnorm(n, 0.5, 1),
           N_1_1 = stats::rnorm(n, 1, 1),
           N_0_2 = stats::rnorm(n, 0, sqrt(2)),
           Cauchy = stats::rcauchy(n, 0, 1),
           dyn_mu_half = {
             mu <- stats::rnorm(1, 0, sqrt(1 / 2))
             stats::rnorm(n, mu, 1)
           },
           dyn_mu_quarter = {
             mu <- stats::rnorm(1, 0, sqrt(1 / 4))
             stats::rnorm(n, mu, 1)
           },
           dyn_var_chi1 = {
             sig2 <- stats::rchisq(1, df = 1)
             stats::rnorm(n, 0, sqrt(sig2))
           },
           dyn_var_chi2 = {
             sig2 <- stats::rchisq(1, df = 2)
             stats::rnorm(n, 0, sqrt(sig2))
           },
           stop("Unknown phase-II generator name: ", name))
  }
}

phase2_label <- function(name) {
  switch(name,
         IC = "IC N(0,1)",
         N_0.5_1 = "N(1/2,1)",
         N_1_1 = "N(1,1)",
         N_0_2 = "N(0,2)",
         Cauchy = "Cauchy(0,1)",
         dyn_mu_half = "N(mu_t,1), mu_t~N(0,1/2)",
         dyn_mu_quarter = "N(mu_t,1), mu_t~N(0,1/4)",
         dyn_var_chi1 = "N(0,sigma_t^2), sigma_t^2~chi^2_1",
         dyn_var_chi2 = "N(0,sigma_t^2), sigma_t^2~chi^2_2",
         name)
}

# -----------------------------------------------------------------------------
# Multi-chart KS stream simulation
# -----------------------------------------------------------------------------

run_stream_ks_multi <- function(charts, n0, n_sampler,
                                phase1_gen = phase1_normal,
                                phase2_gen = make_phase2_generator("IC"),
                                k_values = c(1L, 5L),
                                max_t = 100000L,
                                ks_engine = c("asymp", "hybrid", "psmirnov_exact", "psmirnov_asymp", "stats_exact", "stats_asymp", "asymp_tie_safe"),
                                exact_max_product = 15000L,
                                jmax = 50L,
                                finite_sample_correction = FALSE,
                                return_last_state = FALSE) {
  ks_engine <- match.arg(ks_engine)
  stopifnot(length(charts) >= 1L, n0 >= 1L, max_t >= 1L)

  x0_sorted <- sort(phase1_gen(n0))
  n_charts <- length(charts)
  states <- vector("list", n_charts)
  states[] <- list(list())

  k_values <- sort(unique(as.integer(k_values)))
  k_max <- max(k_values)
  n_k <- length(k_values)

  alarm_counts <- integer(n_charts)
  finished <- rep(FALSE, n_charts)
  Rmat <- matrix(NA_integer_, nrow = n_charts, ncol = n_k)
  colnames(Rmat) <- paste0("R", k_values)

  t <- 0L
  while (t < max_t && !all(finished)) {
    t <- t + 1L
    nt <- n_sampler()
    xt <- phase2_gen(nt, t)
    kp <- ks_p_and_d_fast(x0_sorted, xt, engine = ks_engine,
                          exact_max_product = exact_max_product, jmax = jmax,
                          finite_sample_correction = finite_sample_correction)

    for (i in seq_len(n_charts)) {
      if (finished[i]) next
      upd <- update_chart(states[[i]], charts[[i]], p_value = kp$p, ks_stat = kp$D)
      states[[i]] <- upd$state
      if (upd$alarm) {
        alarm_counts[i] <- alarm_counts[i] + 1L
        hit <- which(k_values == alarm_counts[i])
        if (length(hit) == 1L) Rmat[i, hit] <- t
        if (alarm_counts[i] >= k_max) finished[i] <- TRUE
      }
    }
  }

  params <- chart_param_frame(charts)
  rows <- vector("list", n_charts * n_k)
  idx <- 1L
  for (i in seq_len(n_charts)) {
    for (j in seq_len(n_k)) {
      cens <- is.na(Rmat[i, j])
      rows[[idx]] <- data.frame(
        chart_index = i,
        chart = params$chart[i],
        chart_type = params$chart_type[i],
        threshold = params$threshold[i],
        validity = params$validity[i],
        lambda = params$lambda[i],
        r = params$r[i],
        kappa = params$kappa[i],
        k = k_values[j],
        Rk = if (cens) max_t + 1L else Rmat[i, j],
        censored = cens,
        stopped_at = t,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }
  out <- do.call(rbind, rows)
  if (isTRUE(return_last_state)) attr(out, "states") <- states
  out
}

summarise_run_df <- function(df) {
  split_key <- interaction(df$chart_index, df$k, drop = TRUE)
  pieces <- split(df, split_key)
  rows <- lapply(pieces, function(z) {
    k <- unique(z$k)
    alpha <- unique(z$threshold)
    validity <- unique(z$validity)
    alpha <- if (length(alpha) == 1L) alpha else NA_real_
    validity <- if (length(validity) == 1L) validity else NA_character_
    finite_R <- z$Rk[is.finite(z$Rk)]
    mean_R <- mean(finite_R)
    se_R <- se(finite_R)
    ci <- mc_ci(finite_R)
    lower <- bound_for_chart(alpha, k, validity)
    data.frame(
      chart_index = unique(z$chart_index),
      chart = unique(z$chart),
      chart_type = unique(z$chart_type),
      threshold = alpha,
      validity = validity,
      lambda = unique(z$lambda),
      r = unique(z$r),
      kappa = unique(z$kappa),
      k = k,
      n_rep = nrow(z),
      mean_Rk = mean_R,
      sd_Rk = stats::sd(finite_R),
      se_Rk = se_R,
      ci_lower = ci[1],
      ci_upper = ci[2],
      mcse_ratio = ifelse(mean_R > 0, se_R / mean_R, NA_real_),
      censored_rate = mean(z$censored),
      lower_bound = lower,
      ratio = ifelse(is.finite(lower) && lower > 0, mean_R / lower, NA_real_),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out[order(out$chart_index, out$k), ]
}

.print_partial_summary <- function(df, batch_label = "", verbose = TRUE) {
  if (!isTRUE(verbose) || nrow(df) == 0L) return(invisible(NULL))
  k_min <- min(df$k)
  z <- df[df$k == k_min, , drop = FALSE]
  sm <- aggregate(cbind(Rk, censored) ~ chart, data = z,
                  FUN = function(x) if (is.logical(x)) mean(x) else mean(x, na.rm = TRUE))
  names(sm)[names(sm) == "Rk"] <- paste0("mean_R", k_min)
  names(sm)[names(sm) == "censored"] <- "censor_rate"
  log_msg("Partial summary %s for k=%d:", batch_label, k_min, verbose = verbose)
  print(utils::head(sm, 12L), row.names = FALSE)
  invisible(NULL)
}

run_replicates_ks_multi <- function(n_rep, seed, charts, n0, n_sampler,
                                    phase1_gen = phase1_normal,
                                    phase2_gen = make_phase2_generator("IC"),
                                    k_values = c(1L, 5L), max_t = 100000L,
                                    ks_engine = "asymp",
                                    exact_max_product = 15000L,
                                    jmax = 50L,
                                    finite_sample_correction = FALSE,
                                    n_cores = 1L,
                                    batch_size = 25L,
                                    verbose = TRUE,
                                    checkpoint_file = NULL) {
  stopifnot(n_rep >= 1L, batch_size >= 1L)
  ensure_dir(dirname(checkpoint_file %||% "output/checkpoints/dummy.csv"))

  all_ans <- vector("list", 0L)
  batches <- split(seq_len(n_rep), ceiling(seq_len(n_rep) / batch_size))
  t0 <- Sys.time()

  log_msg("Starting %d KS replications: n0=%d, charts=%d, max_t=%d, engine=%s, cores=%d, batch_size=%d",
          n_rep, n0, length(charts), max_t, ks_engine, n_cores, batch_size, verbose = verbose)

  for (b in seq_along(batches)) {
    ids <- batches[[b]]
    worker <- function(i) {
      set.seed(seed + i)
      out <- run_stream_ks_multi(
        charts = charts, n0 = n0, n_sampler = n_sampler,
        phase1_gen = phase1_gen, phase2_gen = phase2_gen,
        k_values = k_values, max_t = max_t, ks_engine = ks_engine,
        exact_max_product = exact_max_product, jmax = jmax,
        finite_sample_correction = finite_sample_correction
      )
      out$replicate <- i
      out
    }

    if (n_cores > 1L && .Platform$OS.type != "windows") {
      batch_ans <- parallel::mclapply(ids, worker, mc.cores = n_cores)
    } else {
      batch_ans <- lapply(ids, worker)
    }
    all_ans <- c(all_ans, batch_ans)

    partial <- do.call(rbind, all_ans)
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    log_msg("Finished batch %d/%d; reps done=%d/%d; elapsed=%s",
            b, length(batches), length(all_ans), n_rep, format_elapsed(elapsed), verbose = verbose)
    .print_partial_summary(partial, sprintf("after %d reps", length(all_ans)), verbose = verbose)

    if (!is.null(checkpoint_file)) {
      checkpoint_rds <- sub("\\.csv$", "_summary.rds", checkpoint_file)
      saveRDS(summarise_run_df(partial), checkpoint_rds, compress = "xz")
      log_msg("Checkpoint summary written to %s", checkpoint_rds, verbose = verbose)
    }
  }

  out <- do.call(rbind, all_ans)
  rownames(out) <- NULL
  out
}

estimate_arl_ks_multi <- function(charts, n_rep, seed, n0, max_t,
                                  k_values = c(1L, 5L), n_cores = 1L,
                                  phase2_name = "IC", ks_engine = "asymp",
                                  exact_max_product = 15000L,
                                  jmax = 50L,
                                  finite_sample_correction = FALSE,
                                  vss_width = 10L,
                                  batch_size = 25L,
                                  verbose = TRUE,
                                  checkpoint_file = NULL) {
  df <- run_replicates_ks_multi(
    n_rep = n_rep, seed = seed, charts = charts, n0 = n0,
    n_sampler = function() sample_vss(n0, vss_width),
    phase1_gen = phase1_normal,
    phase2_gen = make_phase2_generator(phase2_name),
    k_values = k_values, max_t = max_t, ks_engine = ks_engine,
    exact_max_product = exact_max_product, jmax = jmax,
    finite_sample_correction = finite_sample_correction,
    n_cores = n_cores, batch_size = batch_size,
    verbose = verbose, checkpoint_file = checkpoint_file
  )
  summarise_run_df(df)
}

estimate_ic_ks_stat_mean <- function(n_rep = 1000L, n0 = 50L, seed = 1L,
                                     vss_width = 10L, ks_engine = "asymp",
                                     jmax = 50L, verbose = TRUE) {
  D <- numeric(n_rep)
  t0 <- Sys.time()
  for (i in seq_len(n_rep)) {
    set.seed(seed + i)
    x0 <- sort(phase1_normal(n0))
    nt <- sample_vss(n0, vss_width)
    xt <- phase1_normal(nt)
    D[i] <- ks_p_and_d_fast(x0, xt, engine = ks_engine, jmax = jmax)$D
    if (verbose && (i == 1L || i %% max(1L, floor(n_rep / 10L)) == 0L)) {
      log_msg("Estimating IC mean KS D: n0=%d rep=%d/%d mean_D_so_far=%.5f elapsed=%s",
              n0, i, n_rep, mean(D[seq_len(i)]),
              format_elapsed(as.numeric(difftime(Sys.time(), t0, units = "secs"))),
              verbose = verbose)
    }
  }
  mean(D)
}

# -----------------------------------------------------------------------------
# Threshold calibration for benchmark charts
# -----------------------------------------------------------------------------

calibrate_threshold_ks_fast <- function(chart_template, target_arl,
                                        direction = c("lower", "upper"),
                                        lower, upper, log_scale = FALSE,
                                        n_iter = 12L, n_rep_cal = 200L,
                                        seed = 1000L, n0 = 50L, max_t = 100000L,
                                        n_cores = 1L, phase2_name = "IC",
                                        ks_engine = "asymp",
                                        exact_max_product = 15000L,
                                        jmax = 50L,
                                        vss_width = 10L,
                                        batch_size = 25L,
                                        verbose = TRUE,
                                        checkpoint_file = NULL) {
  direction <- match.arg(direction)
  stopifnot(lower > 0, upper > lower, target_arl > 0)

  history <- data.frame(iter = integer(), threshold = numeric(), mean_R1 = numeric(),
                        se_R1 = numeric(), censored_rate = numeric(), stringsAsFactors = FALSE)
  lo <- lower
  hi <- upper

  log_msg("Calibration start: target ARL=%g, chart=%s, n0=%d, direction=%s, bracket=[%.6g, %.6g]",
          target_arl, make_chart_label(chart_template), n0, direction, lo, hi, verbose = verbose)

  for (iter in seq_len(n_iter)) {
    mid <- if (isTRUE(log_scale)) exp((log(lo) + log(hi)) / 2) else (lo + hi) / 2
    chart <- chart_template
    chart$threshold <- mid

    iter_checkpoint <- NULL
    if (!is.null(checkpoint_file)) {
      iter_checkpoint <- sub("\\.csv$", sprintf("_iter%02d_raw.csv", iter), checkpoint_file)
    }

    sm <- estimate_arl_ks_multi(
      charts = list(chart), n_rep = n_rep_cal, seed = seed + 10000L * iter,
      n0 = n0, max_t = max_t, k_values = 1L, n_cores = n_cores,
      phase2_name = phase2_name, ks_engine = ks_engine,
      exact_max_product = exact_max_product, jmax = jmax,
      vss_width = vss_width, batch_size = batch_size,
      verbose = FALSE, checkpoint_file = iter_checkpoint
    )
    m <- sm$mean_Rk[1]
    s <- sm$se_Rk[1]
    c_rate <- sm$censored_rate[1]
    history <- rbind(history, data.frame(iter = iter, threshold = mid, mean_R1 = m,
                                         se_R1 = s, censored_rate = c_rate))
    log_msg("Calibration iter %02d/%02d: h=%.8g, mean_R1=%.3f, se=%.3f, censor=%.3f, bracket=[%.8g, %.8g]",
            iter, n_iter, mid, m, s, c_rate, lo, hi, verbose = verbose)

    # For lower-tail charts, increasing h increases alarm probability and decreases ARL.
    # For upper-tail charts, increasing h decreases alarm probability and increases ARL.
    if (direction == "lower") {
      if (m > target_arl) lo <- mid else hi <- mid
    } else {
      if (m > target_arl) hi <- mid else lo <- mid
    }

    if (!is.null(checkpoint_file)) {
      hist_out <- history
      hist_out$target_arl <- target_arl
      hist_out$chart <- make_chart_label(chart_template)
      utils::write.csv(hist_out, checkpoint_file, row.names = FALSE)
    }
  }

  score <- abs(log(history$mean_R1 / target_arl))
  best <- history[which.min(score), , drop = FALSE]
  log_msg("Calibration done: selected h=%.8g with mean_R1=%.3f for target=%g",
          best$threshold, best$mean_R1, target_arl, verbose = verbose)

  list(threshold = best$threshold, history = history, selected = best)
}

# -----------------------------------------------------------------------------
# Common chart factories
# -----------------------------------------------------------------------------

make_p_chart <- function(alpha, validity = "marginal", label = NULL) {
  list(type = "p_shewhart", threshold = alpha, validity = validity,
       label = label %||% sprintf("P_t (alpha=%.3g)", alpha))
}

make_qtilde_chart <- function(alpha, lambda, r, validity = "marginal", label = NULL) {
  list(type = "qtilde", threshold = alpha, lambda = lambda, r = r, validity = validity,
       label = label %||% sprintf("Qtilde(lambda=%.2f,r=%.2f,alpha=%.3g)", lambda, r, alpha))
}

make_qbar_chart <- function(alpha, lambda, r = 1, validity = "marginal", label = NULL) {
  list(type = "qbar", threshold = alpha, lambda = lambda, r = r, validity = validity,
       label = label %||% sprintf("Qbar(lambda=%.2f,r=%.2f,alpha=%.3g)", lambda, r, alpha))
}

make_ks_shewhart_chart <- function(threshold = NA_real_, label = NULL) {
  list(type = "ks_shewhart", threshold = threshold,
       label = label %||% sprintf("KS Shewhart(h=%.4g)", threshold))
}

make_ks_ewma_chart <- function(lambda = 0.2, threshold = NA_real_, label = NULL) {
  list(type = "ks_ewma", lambda = lambda, threshold = threshold,
       label = label %||% sprintf("KS-EWMA(lambda=%.2f,h=%.4g)", lambda, threshold))
}

make_ks_cusum_chart <- function(kappa, threshold = NA_real_, label = NULL) {
  list(type = "ks_cusum", kappa = kappa, threshold = threshold,
       label = label %||% sprintf("KS-CUSUM(kappa=%.4g,h=%.4g)", kappa, threshold))
}
