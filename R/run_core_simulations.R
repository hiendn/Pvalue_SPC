#!/usr/bin/env Rscript

# Base-R simulation for the five charts reported in the revised manuscript.
#
# Five charts are compared on exactly the same Phase-I/Phase-II streams:
#   1. raw exact two-sample KS p-value chart;
#   2. marginally valid Qtilde(lambda=.95,r=1);
#   3. ordinary lower p-EWMA(lambda=.2);
#   4. upper KS-statistic EWMA(lambda=.2); and
#   5. Bakir's published two-sided Shew-KS chart.
#
# Usage:
#   Rscript R/run_core_simulations.R pilot
#   Rscript R/run_core_simulations.R full
#
# Only base/recommended R packages are used.  The exact KS implementation is
# for continuous, tie-free data; a unit check against ks.test(exact=TRUE) is
# run before any simulation.

options(stringsAsFactors = FALSE, warn = 1)
RNGkind("L'Ecuyer-CMRG")

args <- commandArgs(trailingOnly = TRUE)
MODE <- if (length(args)) tolower(args[[1L]]) else "pilot"
if (!MODE %in% c("pilot", "full")) stop("Mode must be 'pilot' or 'full'.")

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
SCRIPT_FILE <- if (length(script_arg)) {
  normalizePath(sub("^--file=", "", script_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath("R/run_core_simulations.R", mustWork = TRUE)
}
ROOT <- dirname(SCRIPT_FILE)
OUT <- file.path(ROOT, "output", MODE)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

CFG <- if (MODE == "pilot") {
  list(
    mode = MODE, seed = 26081800L, m = 100L,
    n_values = c(3L, 5L, 8L, 20L),
    target_arl = 200, lambda_q = 0.95, lambda_ewma = 0.20,
    n_unit = 200L, n_centre = 1500L,
    n_cal = 500L, h_cal = 700L,
    n_verify = 800L, h_verify = 1800L,
    n_ooc = 500L, n_delayed = 700L, max_delay = 900L,
    change_times = c(1L, 50L), chunk = 96L,
    scenarios = c("normal_location", "normal_scale", "cauchy"),
    observation_budgets = c(50L, 100L, 200L),
    vsi_n = 5L, n_vsi = 600L, max_vsi = 1200L,
    worked_example_batches = 12L
  )
} else {
  list(
    mode = MODE, seed = 26081800L, m = 100L,
    n_values = c(3L, 5L, 8L, 20L),
    target_arl = 200, lambda_q = 0.95, lambda_ewma = 0.20,
    n_unit = 1000L, n_centre = 20000L,
    n_cal = 10000L, h_cal = 1600L,
    n_verify = 15000L, h_verify = 5000L,
    n_ooc = 8000L, n_delayed = 10000L, max_delay = 2000L,
    change_times = c(1L, 50L), chunk = 128L,
    scenarios = c("normal_location", "normal_scale", "cauchy"),
    observation_budgets = c(50L, 100L, 200L),
    vsi_n = 5L, n_vsi = 6000L, max_vsi = 3000L,
    worked_example_batches = 15L
  )
}

CHARTS <- c("Raw exact KS p", "Qtilde(0.95,1)",
            "p-EWMA(0.2)", "KS-EWMA(0.2)", "Bakir Shew-KS")

stage_seed <- function(n, stage) {
  as.integer(CFG$seed + 1000L * as.integer(n) + as.integer(stage))
}

write_csv <- function(x, name) {
  write.csv(x, file.path(OUT, name), row.names = FALSE, na = "")
}

log_line <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ",
                paste0(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = file.path(OUT, "run.log"), append = TRUE)
}

# ---- Exact two-sample KS calculation ---------------------------------------

make_exact_ks_table <- function(m, n) {
  mn <- as.integer(m * n)
  p <- getFromNamespace("psmirnov", "stats")(
    (0:mn) / mn, sizes = c(m, n), alternative = "two.sided",
    exact = TRUE, lower.tail = FALSE
  )
  list(m = as.integer(m), n = as.integer(n), mn = mn, p = p)
}

# Vectorised over the rows of y. x0_sorted is one reused Phase-I sample.
ks_exact_rows <- function(x0_sorted, y, tab) {
  if (is.null(dim(y))) y <- matrix(y, nrow = 1L)
  m <- tab$m
  n <- tab$n
  stopifnot(length(x0_sorted) == m, ncol(y) == n)
  nr <- nrow(y)

  # Within-row ranks of each Phase-II observation (ties have probability zero
  # in all simulations below).
  y_rank <- matrix(1L, nr, n)
  if (n > 1L) {
    for (j in seq_len(n)) {
      for (k in seq_len(n)) {
        if (k != j) y_rank[, j] <- y_rank[, j] + (y[, k] < y[, j])
      }
    }
  }

  k_max <- integer(nr)
  placements <- matrix(0L, nr, n)
  for (j in seq_len(n)) {
    rj <- findInterval(y[, j], x0_sorted)
    placements[cbind(seq_len(nr), y_rank[, j])] <- rj
    kj <- pmax.int(n * rj - m * (y_rank[, j] - 1L),
                   m * y_rank[, j] - n * rj)
    k_max <- pmax.int(k_max, kj)
  }

  # Bakir's two-sided Shew-KS statistic maximises only over the ordered
  # Phase-I observations.  If k_j Phase-I observations do not exceed the
  # jth ordered Phase-II observation, then St(x_(i))=j/n for
  # i=k_j+1,...,k_(j+1).  A linear absolute difference attains its maximum
  # at an endpoint, so only these O(n) endpoints need to be checked.
  bakir_key <- n * placements[, 1L]
  if (n > 1L) {
    for (j in seq_len(n - 1L)) {
      lo <- placements[, j] + 1L
      hi <- placements[, j + 1L]
      valid <- lo <= hi
      candidate <- integer(nr)
      candidate[valid] <- pmax.int(
        abs(n * lo[valid] - m * j),
        abs(n * hi[valid] - m * j)
      )
      bakir_key <- pmax.int(bakir_key, candidate)
    }
  }
  lo <- placements[, n] + 1L
  valid <- lo <= m
  candidate <- integer(nr)
  candidate[valid] <- abs(n * lo[valid] - m * n)
  bakir_key <- pmax.int(bakir_key, candidate)

  list(D = k_max / tab$mn, p = tab$p[k_max + 1L], key = k_max,
       bakir = bakir_key / tab$mn, bakir_key = bakir_key)
}

unit_check_exact_ks <- function(m, n, tab, B) {
  d_diff <- numeric(B)
  p_diff <- numeric(B)
  bakir_diff <- numeric(B)
  for (b in seq_len(B)) {
    x <- rnorm(m)
    y <- rnorm(n)
    fast <- ks_exact_rows(sort(x), matrix(y, nrow = 1L), tab)
    ref <- suppressWarnings(ks.test(x, y, exact = TRUE))
    x_sorted <- sort(x)
    bakir_ref <- max(abs(seq_len(m) / m -
                         vapply(x_sorted, function(z) mean(y <= z),
                                numeric(1))))
    d_diff[b] <- abs(fast$D - unname(ref$statistic))
    p_diff[b] <- abs(fast$p - ref$p.value)
    bakir_diff[b] <- abs(fast$bakir - bakir_ref)
  }
  data.frame(n = n, comparisons = B,
             max_abs_D_difference = max(d_diff),
             max_abs_p_difference = max(p_diff),
             max_abs_Bakir_difference = max(bakir_diff),
             passed = max(d_diff) < 1e-12 && max(p_diff) < 1e-12 &&
               max(bakir_diff) < 1e-12)
}

generate_ic_bank <- function(n, B, H, tab) {
  p <- matrix(NA_real_, B, H)
  d <- matrix(NA_real_, B, H)
  bakir <- matrix(NA_real_, B, H)
  for (b in seq_len(B)) {
    x0 <- sort(rnorm(CFG$m))
    y <- matrix(rnorm(H * n), nrow = H, ncol = n)
    z <- ks_exact_rows(x0, y, tab)
    p[b, ] <- z$p
    d[b, ] <- z$D
    bakir[b, ] <- z$bakir
  }
  list(p = p, d = d, bakir = bakir)
}

estimate_ic_centres <- function(n, tab, B) {
  p <- d <- bakir <- numeric(B)
  for (b in seq_len(B)) {
    z <- ks_exact_rows(sort(rnorm(CFG$m)), matrix(rnorm(n), 1L, n), tab)
    p[b] <- z$p
    d[b] <- z$D
    bakir[b] <- z$bakir
  }
  data.frame(n = n, n_draws = B, mean_p = mean(p), se_mean_p = sd(p) / sqrt(B),
             mean_D = mean(d), se_mean_D = sd(d) / sqrt(B),
             mean_Bakir = mean(bakir),
             se_mean_Bakir = sd(bakir) / sqrt(B))
}

# ---- Chart paths and threshold calibration --------------------------------

ewma_sequence <- function(x, lambda, initial) {
  as.numeric(stats::filter(lambda * x, filter = 1 - lambda,
                           method = "recursive", init = initial))
}

chart_bank <- function(p, d, bakir, mean_p, mean_d) {
  B <- nrow(p)
  H <- ncol(p)
  q <- matrix(NA_real_, B, H)
  pe <- matrix(NA_real_, B, H)
  ke <- matrix(NA_real_, B, H)

  s <- p[, 1L]
  q[, 1L] <- s / CFG$lambda_q
  w <- rep(mean_p, B)
  z <- rep(mean_d, B)
  for (tt in seq_len(H)) {
    if (tt > 1L) {
      s <- CFG$lambda_q * p[, tt] + (1 - CFG$lambda_q) * s
      q[, tt] <- s / CFG$lambda_q
    }
    w <- CFG$lambda_ewma * p[, tt] + (1 - CFG$lambda_ewma) * w
    z <- CFG$lambda_ewma * d[, tt] + (1 - CFG$lambda_ewma) * z
    pe[, tt] <- w
    ke[, tt] <- z
  }
  list(raw = p, qtilde = q, p_ewma = pe, ks_ewma = ke,
       bakir = bakir)
}

first_crossing <- function(stat, threshold, direction = c("lower", "upper"),
                           boundary_u = NULL, boundary_gamma = NA_real_) {
  direction <- match.arg(direction)
  if (!is.null(boundary_u)) {
    hit <- if (direction == "lower") {
      stat < threshold | (stat == threshold & boundary_u <= boundary_gamma)
    } else {
      stat > threshold | (stat == threshold & boundary_u <= boundary_gamma)
    }
  } else if (direction == "lower") {
    hit <- stat <= threshold
  } else {
    hit <- stat >= threshold
  }
  has <- rowSums(hit) > 0L
  ans <- max.col(hit, ties.method = "first")
  ans[!has] <- ncol(stat) + 1L
  as.integer(ans)
}

summarise_rl <- function(r, horizon) {
  se <- sd(r) / sqrt(length(r))
  data.frame(n_rep = length(r), mean_RL = mean(r), sd_RL = sd(r), se_RL = se,
             ci_lower = mean(r) - 1.96 * se,
             ci_upper = mean(r) + 1.96 * se,
             relative_MCSE = se / mean(r),
             censored_rate = mean(r == horizon + 1L))
}

calibrate_continuous <- function(stat, direction, target, iterations = 30L) {
  H <- ncol(stat)
  eps <- .Machine$double.eps^0.5
  if (direction == "lower") {
    lo <- min(stat) - eps
    hi <- max(stat) + eps
    for (i in seq_len(iterations)) {
      mid <- (lo + hi) / 2
      mr <- mean(first_crossing(stat, mid, "lower"))
      if (mr > target) lo <- mid else hi <- mid
    }
    threshold <- (lo + hi) / 2
  } else {
    lo <- min(stat) - eps
    hi <- max(stat) + eps
    for (i in seq_len(iterations)) {
      mid <- (lo + hi) / 2
      mr <- mean(first_crossing(stat, mid, "upper"))
      if (mr > target) hi <- mid else lo <- mid
    }
    threshold <- (lo + hi) / 2
  }
  r <- first_crossing(stat, threshold, direction)
  list(threshold = threshold, gamma = NA_real_, r = r,
       summary = summarise_rl(r, H))
}

calibrate_randomised_boundary <- function(stat, u, target,
                                           direction = c("lower", "upper")) {
  direction <- match.arg(direction)
  H <- ncol(stat)
  vals <- sort(unique(as.vector(stat)))
  left <- 1L
  right <- length(vals)
  if (direction == "lower") {
    while (left < right) {
      mid <- floor((left + right) / 2)
      mr <- mean(first_crossing(stat, vals[mid], direction))
      if (mr <= target) right <- mid else left <- mid + 1L
    }
  } else {
    while (left < right) {
      mid <- ceiling((left + right) / 2)
      mr <- mean(first_crossing(stat, vals[mid], direction))
      if (mr <= target) left <- mid else right <- mid - 1L
    }
  }
  boundary <- vals[left]
  r_gamma0 <- first_crossing(stat, boundary, direction, u, 0)
  r_gamma1 <- first_crossing(stat, boundary, direction, u, 1)
  lo <- 0
  hi <- 1
  for (i in seq_len(30L)) {
    gamma <- (lo + hi) / 2
    mr <- mean(first_crossing(stat, boundary, direction, u, gamma))
    if (mr > target) lo <- gamma else hi <- gamma
  }
  gamma <- (lo + hi) / 2
  r <- first_crossing(stat, boundary, direction, u, gamma)
  list(threshold = boundary, gamma = gamma, r = r,
       deterministic_gamma0_mean = mean(r_gamma0),
       deterministic_gamma1_mean = mean(r_gamma1),
       summary = summarise_rl(r, H))
}

calibrate_all <- function(bank, mean_p, mean_d) {
  cb <- chart_bank(bank$p, bank$d, bank$bakir, mean_p, mean_d)
  u_raw <- matrix(runif(length(bank$p)), nrow(bank$p), ncol(bank$p))
  u_bakir <- matrix(runif(length(bank$p)), nrow(bank$p), ncol(bank$p))
  out <- list(
    raw = calibrate_randomised_boundary(cb$raw, u_raw, CFG$target_arl,
                                        "lower"),
    qtilde = calibrate_continuous(cb$qtilde, "lower", CFG$target_arl),
    p_ewma = calibrate_continuous(cb$p_ewma, "lower", CFG$target_arl),
    ks_ewma = calibrate_continuous(cb$ks_ewma, "upper", CFG$target_arl),
    bakir = calibrate_randomised_boundary(cb$bakir, u_bakir,
                                          CFG$target_arl, "upper")
  )
  names(out) <- CHARTS
  out
}

threshold_frame <- function(n, cal) {
  direction <- c("lower", "lower", "lower", "upper", "upper")
  do.call(rbind, lapply(seq_along(cal), function(j) {
    z <- cal[[j]]
    data.frame(n = n, chart = names(cal)[j], direction = direction[j],
               threshold = z$threshold, boundary_gamma = z$gamma,
               reported_limit = z$threshold,
               threshold_scale = "chart scale",
               lambda = c(NA, CFG$lambda_q, CFG$lambda_ewma,
                          CFG$lambda_ewma, NA)[j],
               r = c(NA, 1, NA, NA, NA)[j],
               deterministic_gamma0_mean =
                 if (!is.null(z$deterministic_gamma0_mean))
                   z$deterministic_gamma0_mean else NA_real_,
               deterministic_gamma1_mean =
                 if (!is.null(z$deterministic_gamma1_mean))
                   z$deterministic_gamma1_mean else NA_real_,
               calibration_horizon = CFG$h_cal,
               z$summary, check.names = FALSE)
  }))
}

# ---- Independent path simulation ------------------------------------------

draw_phase2 <- function(number, n, scenario) {
  switch(scenario,
    ic = matrix(rnorm(number * n), number, n),
    normal_location = matrix(rnorm(number * n, mean = 0.5), number, n),
    normal_scale = matrix(rnorm(number * n, sd = sqrt(2)), number, n),
    cauchy = matrix(rcauchy(number * n), number, n),
    stop("Unknown scenario: ", scenario)
  )
}

simulate_first_alarm_paths <- function(n, tab, cal, mean_p, mean_d, N,
                                       scenario = "ic", change_time = 1L,
                                       max_delay = CFG$max_delay) {
  stopifnot(change_time >= 1L)
  max_total <- as.integer(change_time - 1L + max_delay)
  Tmat <- matrix(max_total + 1L, N, length(CHARTS),
                 dimnames = list(NULL, CHARTS))

  thresholds <- vapply(cal, `[[`, numeric(1), "threshold")
  gamma_raw <- cal[[1L]]$gamma
  gamma_bakir <- cal[[5L]]$gamma

  for (b in seq_len(N)) {
    x0 <- sort(rnorm(CFG$m))
    first <- rep(NA_integer_, length(CHARTS))
    s_q <- NA_real_
    w_p <- mean_p
    w_d <- mean_d
    tt <- 0L

    while (anyNA(first) && tt < max_total) {
      len <- min(CFG$chunk, max_total - tt)
      times <- tt + seq_len(len)
      y <- matrix(NA_real_, len, n)
      before <- times < change_time
      if (any(before)) y[before, ] <- draw_phase2(sum(before), n, "ic")
      if (any(!before)) y[!before, ] <- draw_phase2(sum(!before), n, scenario)
      zd <- ks_exact_rows(x0, y, tab)
      p <- zd$p
      d <- zd$D
      bakir <- zd$bakir

      if (is.na(s_q)) {
        sq <- c(p[1L], if (len > 1L) ewma_sequence(
          p[-1L], CFG$lambda_q, p[1L]) else numeric())
      } else {
        sq <- ewma_sequence(p, CFG$lambda_q, s_q)
      }
      q <- sq / CFG$lambda_q
      pe <- ewma_sequence(p, CFG$lambda_ewma, w_p)
      ke <- ewma_sequence(d, CFG$lambda_ewma, w_d)
      u <- runif(len)

      hit <- list(
        p < thresholds[1L] | (p == thresholds[1L] & u <= gamma_raw),
        q <= thresholds[2L],
        pe <= thresholds[3L],
        ke >= thresholds[4L],
        bakir > thresholds[5L] |
          (bakir == thresholds[5L] & u <= gamma_bakir)
      )
      for (j in seq_along(first)) {
        if (is.na(first[j])) {
          jj <- which(hit[[j]])
          if (length(jj)) first[j] <- tt + jj[1L]
        }
      }

      s_q <- tail(sq, 1L)
      w_p <- tail(pe, 1L)
      w_d <- tail(ke, 1L)
      tt <- tt + len
    }
    first[is.na(first)] <- max_total + 1L
    Tmat[b, ] <- first
  }
  Tmat
}

wilson_interval <- function(x, n, z = 1.96) {
  if (!n) return(c(NA_real_, NA_real_))
  phat <- x / n
  den <- 1 + z^2 / n
  centre <- (phat + z^2 / (2 * n)) / den
  half <- z * sqrt(phat * (1 - phat) / n + z^2 / (4 * n^2)) / den
  c(max(0, centre - half), min(1, centre + half))
}

summarise_ic_validation <- function(n, Tmat, horizon) {
  do.call(rbind, lapply(seq_along(CHARTS), function(j) {
    r <- Tmat[, j]
    s <- summarise_rl(r, horizon)
    data.frame(n = n, chart = CHARTS[j], target_ARL = CFG$target_arl,
               ARL_ratio = s$mean_RL / CFG$target_arl,
               within_5_percent = abs(s$mean_RL / CFG$target_arl - 1) <= .05,
               CI_contains_target = s$ci_lower <= CFG$target_arl &
                 s$ci_upper >= CFG$target_arl,
               s, check.names = FALSE)
  }))
}

summarise_ooc <- function(n, scenario, change_time, Tmat, max_delay) {
  do.call(rbind, lapply(seq_along(CHARTS), function(j) {
    T <- Tmat[, j]
    keep <- T >= change_time
    delay <- T[keep] - change_time + 1L
    cens <- delay > max_delay
    restricted <- pmin(delay, max_delay + 1L)
    se <- sd(restricted) / sqrt(length(restricted))
    p10 <- mean(delay <= 10L)
    p25 <- mean(delay <= 25L)
    ci10 <- wilson_interval(sum(delay <= 10L), length(delay))
    ci25 <- wilson_interval(sum(delay <= 25L), length(delay))
    data.frame(
      n = n, scenario = scenario, change_time = change_time,
      chart = CHARTS[j], n_started = nrow(Tmat), n_conditioned = length(delay),
      prechange_signal_rate = mean(!keep),
      mean_delay = mean(restricted), sd_delay = sd(restricted), se_delay = se,
      mean_postchange_observations = n * mean(restricted),
      delay_ci_lower = mean(restricted) - 1.96 * se,
      delay_ci_upper = mean(restricted) + 1.96 * se,
      delay_censored_rate = mean(cens),
      delay_estimand = if (any(cens)) paste0("restricted/capped at ", max_delay + 1L)
                       else "ECD",
      PTS_10 = p10, PTS_10_lower = ci10[1L], PTS_10_upper = ci10[2L],
      PTS_25 = p25, PTS_25_lower = ci25[1L], PTS_25_upper = ci25[2L],
      PTS_obs_50 = mean(n * delay <= 50L),
      PTS_obs_100 = mean(n * delay <= 100L),
      PTS_obs_200 = mean(n * delay <= 200L)
    )
  }))
}

# ---- Main exact-p, calibration, IC and OOC run -----------------------------

unlink(file.path(OUT, "run.log"), force = TRUE)
log_line("Starting ", MODE, " run in ", OUT)
saveRDS(CFG, file.path(OUT, "config.rds"))
writeLines(capture.output(str(CFG)), file.path(OUT, "config.txt"))
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))

tables <- setNames(lapply(CFG$n_values, function(n) make_exact_ks_table(CFG$m, n)),
                   CFG$n_values)

unit_rows <- list()
for (n in CFG$n_values) {
  set.seed(stage_seed(n, 1L))
  unit_rows[[as.character(n)]] <- unit_check_exact_ks(
    CFG$m, n, tables[[as.character(n)]], CFG$n_unit)
}
unit_df <- do.call(rbind, unit_rows)
write_csv(unit_df, "exact_ks_unit_check.csv")
if (!all(unit_df$passed)) stop("Exact KS unit check failed; simulation aborted.")
log_line("Exact-KS and Bakir-statistic unit checks passed.")

centre_rows <- threshold_rows <- ic_rows <- ooc_rows <- list()
replicates <- list()

for (n in CFG$n_values) {
  key <- as.character(n)
  tab <- tables[[key]]
  log_line("n=", n, ": estimating IC centres")
  set.seed(stage_seed(n, 10L))
  centres <- estimate_ic_centres(n, tab, CFG$n_centre)
  centre_rows[[key]] <- centres

  log_line("n=", n, ": generating common calibration bank (", CFG$n_cal,
           " x ", CFG$h_cal, ")")
  set.seed(stage_seed(n, 20L))
  bank <- generate_ic_bank(n, CFG$n_cal, CFG$h_cal, tab)
  set.seed(stage_seed(n, 21L))
  cal <- calibrate_all(bank, centres$mean_p, centres$mean_D)
  threshold_rows[[key]] <- threshold_frame(n, cal)
  rm(bank); invisible(gc())

  log_line("n=", n, ": independent IC validation (N=", CFG$n_verify, ")")
  set.seed(stage_seed(n, 30L))
  T_ic <- simulate_first_alarm_paths(
    n, tab, cal, centres$mean_p, centres$mean_D, CFG$n_verify,
    scenario = "ic", change_time = 1L, max_delay = CFG$h_verify)
  ic_rows[[key]] <- summarise_ic_validation(n, T_ic, CFG$h_verify)
  replicates[[paste0("n", n, "_IC")]] <- T_ic

  for (sc_i in seq_along(CFG$scenarios)) {
    scenario <- CFG$scenarios[sc_i]
    for (nu in CFG$change_times) {
      N <- if (nu == 1L) CFG$n_ooc else CFG$n_delayed
      log_line("n=", n, ": ", scenario, ", change time ", nu,
               ", N=", N)
      set.seed(stage_seed(n, 100L + 10L * sc_i + nu))
      T_ooc <- simulate_first_alarm_paths(
        n, tab, cal, centres$mean_p, centres$mean_D, N,
        scenario = scenario, change_time = nu, max_delay = CFG$max_delay)
      nm <- paste(n, scenario, nu, sep = "_")
      ooc_rows[[nm]] <- summarise_ooc(n, scenario, nu, T_ooc,
                                      CFG$max_delay)
      replicates[[paste0("n", nm)]] <- T_ooc
    }
  }
  saveRDS(list(calibration = cal, centres = centres),
          file.path(OUT, paste0("calibration_n", n, ".rds")))
}

write_csv(do.call(rbind, centre_rows), "ic_centres.csv")
write_csv(do.call(rbind, threshold_rows), "calibrated_thresholds.csv")
write_csv(do.call(rbind, ic_rows), "ic_validation.csv")
write_csv(do.call(rbind, ooc_rows), "ooc_delay_metrics.csv")
saveRDS(replicates, file.path(OUT, "replicate_run_lengths.rds"), compress = TRUE)

log_line("Core exact-p simulation completed.")
