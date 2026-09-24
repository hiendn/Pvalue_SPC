#!/usr/bin/env Rscript

# Independent VSI illustration and numerical implementation example for the
# five-chart simulation. Run only after run_core_simulations.R has saved
# calibration_n5.rds for the selected mode.
#
# Usage:
#   Rscript R/postprocess_vsi_example.R pilot
#   Rscript R/postprocess_vsi_example.R full

options(stringsAsFactors = FALSE, warn = 1)
RNGkind("L'Ecuyer-CMRG")

args <- commandArgs(trailingOnly = TRUE)
MODE <- if (length(args)) tolower(args[[1L]]) else "pilot"
if (!MODE %in% c("pilot", "full")) stop("Mode must be 'pilot' or 'full'.")
TARGET_ONLY <- identical(tolower(Sys.getenv("SPC_TARGETED_N8_ONLY", "false")),
                         "true")

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
SCRIPT_FILE <- normalizePath(sub("^--file=", "", script_arg[[1L]]),
                             mustWork = TRUE)
ROOT <- dirname(SCRIPT_FILE)
OUT <- file.path(ROOT, "output", MODE)
CORE_CONFIG_FILE <- file.path(OUT, "config.rds")
CAL_FILE <- file.path(OUT, "calibration_n5.rds")
if (!file.exists(CORE_CONFIG_FILE) || !file.exists(CAL_FILE)) {
  stop("Core config/calibration missing in ", OUT,
       "; run run_core_simulations.R first.")
}

CORE <- readRDS(CORE_CONFIG_FILE)
SAVED <- readRDS(CAL_FILE)
CAL <- SAVED$calibration
CENTRES <- SAVED$centres

PP <- if (MODE == "pilot") {
  list(mode = MODE, seed = 26081800L, n = 5L, n_cal_vsi = 500L,
       n_verify_vsi = 800L, max_inspections = 1800L, chunk = 128L,
       n8_cal_per_bank = 800L, n8_cal_horizon = 1200L,
       n8_final_validation = 1500L, n8_final_horizon = 2500L,
       d_short = 0.1, d_long_initial = 1.9, warning_nominal = 0.5,
       first_interval = 1,
       target_ATS = 200, worked_example_batches = 15L,
       worked_change_after = 5L, worked_shift = 1.5)
} else {
  list(mode = MODE, seed = 26081800L, n = 5L, n_cal_vsi = 4000L,
       n_verify_vsi = 10000L, max_inspections = 5000L, chunk = 128L,
       n8_cal_per_bank = 6000L, n8_cal_horizon = 2500L,
       n8_final_validation = 15000L, n8_final_horizon = 5000L,
       d_short = 0.1, d_long_initial = 1.9, warning_nominal = 0.5,
       first_interval = 1,
       target_ATS = 200, worked_example_batches = 15L,
       worked_change_after = 5L, worked_shift = 1.5)
}

CHARTS <- c("Raw exact KS p", "Qtilde(0.95,1)",
            "p-EWMA(0.2)", "KS-EWMA(0.2)", "Bakir Shew-KS")

log_line <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ",
                paste0(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = file.path(OUT, "postprocess.log"), append = TRUE)
}

write_csv <- function(x, name) {
  write.csv(x, file.path(OUT, name), row.names = FALSE, na = "")
}

make_exact_ks_table <- function(m, n) {
  mn <- as.integer(m * n)
  p <- getFromNamespace("psmirnov", "stats")(
    (0:mn) / mn, sizes = c(m, n), alternative = "two.sided",
    exact = TRUE, lower.tail = FALSE
  )
  list(m = as.integer(m), n = as.integer(n), mn = mn, p = p)
}

ks_exact_rows <- function(x0_sorted, y, tab) {
  if (is.null(dim(y))) y <- matrix(y, nrow = 1L)
  m <- tab$m
  n <- tab$n
  nr <- nrow(y)
  y_rank <- matrix(1L, nr, n)
  if (n > 1L) {
    for (j in seq_len(n)) for (k in seq_len(n)) if (k != j) {
      y_rank[, j] <- y_rank[, j] + (y[, k] < y[, j])
    }
  }
  k_max <- integer(nr)
  for (j in seq_len(n)) {
    rj <- findInterval(y[, j], x0_sorted)
    kj <- pmax.int(n * rj - m * (y_rank[, j] - 1L),
                   m * y_rank[, j] - n * rj)
    k_max <- pmax.int(k_max, kj)
  }
  list(D = k_max / tab$mn, p = tab$p[k_max + 1L])
}

# ---- Targeted n=8 recalibration of the two limits that missed validation ---

TARGET_CHARTS <- c("Qtilde(0.95,1)", "p-EWMA(0.2)")

generate_p_bank <- function(n, B, H, tab) {
  p <- matrix(NA_real_, B, H)
  for (b in seq_len(B)) {
    x0 <- sort(rnorm(tab$m))
    y <- matrix(rnorm(H * n), H, n)
    p[b, ] <- ks_exact_rows(x0, y, tab)$p
  }
  p
}

make_qtilde_bank <- function(p, lambda) {
  ans <- matrix(NA_real_, nrow(p), ncol(p))
  state <- p[, 1L]
  ans[, 1L] <- state / lambda
  if (ncol(p) > 1L) for (tt in 2:ncol(p)) {
    state <- lambda * p[, tt] + (1 - lambda) * state
    ans[, tt] <- state / lambda
  }
  ans
}

make_ewma_bank <- function(p, lambda, centre) {
  ans <- matrix(NA_real_, nrow(p), ncol(p))
  state <- rep(centre, nrow(p))
  for (tt in seq_len(ncol(p))) {
    state <- lambda * p[, tt] + (1 - lambda) * state
    ans[, tt] <- state
  }
  ans
}

first_lower_crossing <- function(stat, threshold) {
  hit <- stat <= threshold
  has <- rowSums(hit) > 0L
  ans <- max.col(hit, ties.method = "first")
  ans[!has] <- ncol(stat) + 1L
  as.integer(ans)
}

rl_summary <- function(r, horizon) {
  se <- sd(r) / sqrt(length(r))
  data.frame(n_rep = length(r), mean_RL = mean(r), sd_RL = sd(r),
             se_RL = se, ci_lower = mean(r) - 1.96 * se,
             ci_upper = mean(r) + 1.96 * se,
             relative_MCSE = se / mean(r),
             censored_rate = mean(r == horizon + 1L))
}

calibrate_lower <- function(stat, target, iterations = 24L) {
  lo <- min(stat) - sqrt(.Machine$double.eps)
  hi <- max(stat) + sqrt(.Machine$double.eps)
  for (iter in seq_len(iterations)) {
    mid <- (lo + hi) / 2
    if (mean(first_lower_crossing(stat, mid)) > target) lo <- mid else hi <- mid
  }
  threshold <- (lo + hi) / 2
  r <- first_lower_crossing(stat, threshold)
  list(threshold = threshold, r = r,
       summary = rl_summary(r, ncol(stat)))
}

targeted_calibration_bank <- function(bank_id, tab8, centre8) {
  set.seed(PP$seed + 8000L + bank_id)
  p <- generate_p_bank(8L, PP$n8_cal_per_bank, PP$n8_cal_horizon, tab8)
  q <- make_qtilde_bank(p, CORE$lambda_q)
  cq <- calibrate_lower(q, CORE$target_arl)
  rm(q); invisible(gc())
  pe <- make_ewma_bank(p, CORE$lambda_ewma, centre8$mean_p)
  cp <- calibrate_lower(pe, CORE$target_arl)
  rm(pe, p); invisible(gc())
  list(q = cq, p = cp)
}

simulate_qp_paths <- function(N, scenario, change_time, max_delay, tab8,
                              centre8, thresholds) {
  max_total <- as.integer(change_time - 1L + max_delay)
  T <- matrix(max_total + 1L, N, 2L,
              dimnames = list(NULL, TARGET_CHARTS))
  for (b in seq_len(N)) {
    x0 <- sort(rnorm(tab8$m))
    first <- rep(NA_integer_, 2L)
    s_q <- NA_real_
    w_p <- centre8$mean_p
    tt <- 0L
    while (anyNA(first) && tt < max_total) {
      len <- min(PP$chunk, max_total - tt)
      times <- tt + seq_len(len)
      y <- matrix(NA_real_, len, 8L)
      before <- times < change_time
      if (any(before)) y[before, ] <- draw_phase2(sum(before), 8L, "ic")
      if (any(!before)) y[!before, ] <- draw_phase2(sum(!before), 8L, scenario)
      p <- ks_exact_rows(x0, y, tab8)$p
      if (is.na(s_q)) {
        sq <- c(p[1L], if (len > 1L) ewma_sequence(
          p[-1L], CORE$lambda_q, p[1L]) else numeric())
      } else {
        sq <- ewma_sequence(p, CORE$lambda_q, s_q)
      }
      q <- sq / CORE$lambda_q
      pe <- ewma_sequence(p, CORE$lambda_ewma, w_p)
      hits <- list(q <= thresholds[1L], pe <= thresholds[2L])
      for (j in 1:2) if (is.na(first[j])) {
        jj <- which(hits[[j]])
        if (length(jj)) first[j] <- tt + jj[1L]
      }
      s_q <- tail(sq, 1L)
      w_p <- tail(pe, 1L)
      tt <- tt + len
    }
    first[is.na(first)] <- max_total + 1L
    T[b, ] <- first
  }
  T
}

summarise_target_ic <- function(T, horizon) {
  do.call(rbind, lapply(1:2, function(j) {
    s <- rl_summary(T[, j], horizon)
    data.frame(n = 8L, chart = TARGET_CHARTS[j],
               target_ARL = CORE$target_arl,
               ARL_ratio = s$mean_RL / CORE$target_arl,
               within_5_percent = abs(s$mean_RL / CORE$target_arl - 1) <= .05,
               CI_contains_target = s$ci_lower <= CORE$target_arl &
                 s$ci_upper >= CORE$target_arl,
               s, check.names = FALSE)
  }))
}

summarise_target_ooc <- function(scenario, change_time, T, max_delay) {
  do.call(rbind, lapply(1:2, function(j) {
    tj <- T[, j]
    keep <- tj >= change_time
    delay <- tj[keep] - change_time + 1L
    cens <- delay > max_delay
    restricted <- pmin(delay, max_delay + 1L)
    se <- sd(restricted) / sqrt(length(restricted))
    p10 <- mean(delay <= 10L); p25 <- mean(delay <= 25L)
    ci10 <- wilson_interval(sum(delay <= 10L), length(delay))
    ci25 <- wilson_interval(sum(delay <= 25L), length(delay))
    data.frame(
      n = 8L, scenario = scenario, change_time = change_time,
      chart = TARGET_CHARTS[j], n_started = nrow(T),
      n_conditioned = length(delay), prechange_signal_rate = mean(!keep),
      mean_delay = mean(restricted), sd_delay = sd(restricted), se_delay = se,
      mean_postchange_observations = 8L * mean(restricted),
      delay_ci_lower = mean(restricted) - 1.96 * se,
      delay_ci_upper = mean(restricted) + 1.96 * se,
      delay_censored_rate = mean(cens),
      delay_estimand = if (any(cens)) paste0("restricted/capped at ",
                                             max_delay + 1L) else "ECD",
      PTS_10 = p10, PTS_10_lower = ci10[1L], PTS_10_upper = ci10[2L],
      PTS_25 = p25, PTS_25_lower = ci25[1L], PTS_25_upper = ci25[2L],
      PTS_obs_50 = mean(8L * delay <= 50L),
      PTS_obs_100 = mean(8L * delay <= 100L),
      PTS_obs_200 = mean(8L * delay <= 200L)
    )
  }))
}

draw_phase2 <- function(number, n, scenario) {
  switch(scenario,
    ic = matrix(rnorm(number * n), number, n),
    normal_location = matrix(rnorm(number * n, mean = 0.5), number, n),
    normal_scale = matrix(rnorm(number * n, sd = sqrt(2)), number, n),
    cauchy = matrix(rcauchy(number * n), number, n),
    stop("Unknown scenario: ", scenario)
  )
}

# Simulate the raw action rule, retaining only p-values preceding its signal.
# The retained histories are sufficient to evaluate any two-interval VSI rule
# without resimulating or changing the action event.
simulate_raw_histories <- function(N, scenario, tab, threshold, gamma,
                                   max_inspections, chunk) {
  T <- integer(N)
  pre_p <- vector("list", N)
  for (b in seq_len(N)) {
    x0 <- sort(rnorm(tab$m))
    seen <- numeric()
    signalled <- FALSE
    tt <- 0L
    while (!signalled && tt < max_inspections) {
      len <- min(chunk, max_inspections - tt)
      z <- ks_exact_rows(x0, draw_phase2(len, tab$n, scenario), tab)
      u <- runif(len)
      hit <- z$p < threshold | (z$p == threshold & u <= gamma)
      jj <- which(hit)
      if (length(jj)) {
        j1 <- jj[1L]
        if (j1 > 1L) seen <- c(seen, z$p[seq_len(j1 - 1L)])
        T[b] <- tt + j1
        signalled <- TRUE
      } else {
        seen <- c(seen, z$p)
        tt <- tt + len
      }
    }
    if (!signalled) T[b] <- max_inspections + 1L
    pre_p[[b]] <- seen
  }
  list(T = T, pre_p = pre_p,
       censored = T == max_inspections + 1L)
}

calendar_times <- function(histories, warning_threshold, d_short, d_long,
                           first_interval = 1) {
  first_interval + vapply(histories$pre_p, function(p) {
    sum(ifelse(p <= warning_threshold, d_short, d_long))
  }, numeric(1))
}

calibrate_warning_threshold <- function(histories, action_threshold) {
  candidates <- sort(unique(unlist(histories$pre_p, use.names = FALSE)))
  candidates <- candidates[candidates > action_threshold & candidates < 1]
  if (!length(candidates)) stop("No attainable warning thresholds.")
  w <- candidates[which.min(abs(candidates - PP$warning_nominal))]
  n_short <- vapply(histories$pre_p, function(p) sum(p <= w), numeric(1))
  n_long <- lengths(histories$pre_p) - n_short
  target_fsi <- mean(histories$T)
  d_long <- (target_fsi - PP$first_interval -
             PP$d_short * mean(n_short)) / mean(n_long)
  if (!is.finite(d_long) || d_long <= PP$d_short) {
    stop("Calibrated long interval is not larger than the short interval.")
  }
  times <- calendar_times(histories, w, PP$d_short, d_long,
                          PP$first_interval)
  se <- sd(times) / sqrt(length(times))
  data.frame(
    action_threshold = action_threshold,
    action_boundary_gamma = CAL[[1L]]$gamma,
    warning_nominal = PP$warning_nominal, warning_threshold = w,
    d_short = PP$d_short, d_long_initial = PP$d_long_initial,
    d_long = d_long,
    first_interval = PP$first_interval, n_calibration = length(times),
    target_FSI_ATS = target_fsi, calibrated_VSI_ATS = mean(times),
    se_ATS = se, ci_lower = mean(times) - 1.96 * se,
    ci_upper = mean(times) + 1.96 * se,
    inspection_ARL = mean(histories$T),
    censoring_rate = mean(histories$censored),
    matched_by_linear_interval_solution =
      abs(mean(times) - target_fsi) < 1e-10
  )
}

wilson_interval <- function(x, n, z = 1.96) {
  phat <- x / n
  den <- 1 + z^2 / n
  centre <- (phat + z^2 / (2 * n)) / den
  half <- z * sqrt(phat * (1 - phat) / n + z^2 / (4 * n^2)) / den
  c(max(0, centre - half), min(1, centre + half))
}

summarise_vsi <- function(histories, scenario, warning_threshold, d_long) {
  vsi <- calendar_times(histories, warning_threshold, PP$d_short,
                        d_long, PP$first_interval)
  fsi <- as.numeric(histories$T)
  do.call(rbind, lapply(list(FSI = fsi, VSI = vsi), function(times) {
    se <- sd(times) / sqrt(length(times))
    ci10 <- wilson_interval(sum(times <= 10), length(times))
    ci25 <- wilson_interval(sum(times <= 25), length(times))
    data.frame(
      n = PP$n, scenario = scenario,
      method = if (identical(times, fsi)) "fixed interval (1.0)" else
        paste0("VSI (", PP$d_short, "/",
               formatC(d_long, digits = 4L, format = "f"), ")"),
      estimand = if (scenario == "ic") "IC ATS" else "zero-state OOC ATS",
      n_rep = length(times), mean_inspections = mean(histories$T),
      mean_elapsed_time = mean(times), sd_elapsed_time = sd(times),
      se_elapsed_time = se,
      mean_interval_to_signal = mean(times) / mean(histories$T),
      inspections_per_time_unit = mean(histories$T) / mean(times),
      ci_lower = mean(times) - 1.96 * se,
      ci_upper = mean(times) + 1.96 * se,
      censoring_rate = mean(histories$censored),
      P_signal_by_time_10 = mean(times <= 10),
      P10_lower = ci10[1L], P10_upper = ci10[2L],
      P_signal_by_time_25 = mean(times <= 25),
      P25_lower = ci25[1L], P25_upper = ci25[2L]
    )
  }))
}

ewma_sequence <- function(x, lambda, initial) {
  as.numeric(stats::filter(lambda * x, filter = 1 - lambda,
                           method = "recursive", init = initial))
}

make_worked_example <- function(tab) {
  set.seed(PP$seed + 9005L)
  phase1 <- rnorm(tab$m)
  x0 <- sort(phase1)
  B <- PP$worked_example_batches
  y <- matrix(NA_real_, B, tab$n)
  for (tt in seq_len(B)) {
    y[tt, ] <- if (tt <= PP$worked_change_after) rnorm(tab$n) else
      rnorm(tab$n, mean = PP$worked_shift)
  }
  kd <- ks_exact_rows(x0, y, tab)
  p <- kd$p
  d <- kd$D
  sq <- c(p[1L], ewma_sequence(p[-1L], CORE$lambda_q, p[1L]))
  q <- pmin(1, sq / CORE$lambda_q)
  pe <- ewma_sequence(p, CORE$lambda_ewma, CENTRES$mean_p)
  ke <- ewma_sequence(d, CORE$lambda_ewma, CENTRES$mean_D)
  u <- runif(B)
  thresholds <- vapply(CAL, `[[`, numeric(1), "threshold")
  out <- data.frame(
    inspection = seq_len(B),
    regime = ifelse(seq_len(B) <= PP$worked_change_after, "IC N(0,1)",
                    paste0("OOC N(", PP$worked_shift, ",1)")),
    phase2_values = apply(y, 1L, function(z) paste(formatC(z, digits = 3,
                                                          format = "f"),
                                                   collapse = "; ")),
    phase2_mean = rowMeans(y), KS_D = d, exact_KS_p = p,
    raw_boundary_u = u,
    raw_alarm = p < thresholds[1L] |
      (p == thresholds[1L] & u <= CAL[[1L]]$gamma),
    Qtilde = q, Qtilde_alarm = q <= thresholds[2L],
    p_EWMA = pe, p_EWMA_alarm = pe <= thresholds[3L],
    KS_EWMA = ke, KS_EWMA_alarm = ke >= thresholds[4L]
  )
  list(phase1 = data.frame(index = seq_along(phase1), value = phase1),
       monitoring = out)
}

unlink(file.path(OUT, "postprocess.log"), force = TRUE)
saveRDS(PP, file.path(OUT, "postprocess_config.rds"))
writeLines(capture.output(str(PP)), file.path(OUT, "postprocess_config.txt"))
TAB <- make_exact_ks_table(CORE$m, PP$n)
action_threshold <- CAL[[1L]]$threshold
action_gamma <- CAL[[1L]]$gamma

if (TARGET_ONLY) {
# This legacy recovery entry point is retained for auditing an n=8 limit if a
# future common-core run fails its frozen validation criterion. It is not part
# of the normal publication workflow, which keeps every chart on the common
# core streams.
log_line("Targeted n=8 recalibration: fresh bank 1 of 2")
saved8 <- readRDS(file.path(OUT, "calibration_n8.rds"))
tab8 <- make_exact_ks_table(CORE$m, 8L)
c8_1 <- targeted_calibration_bank(1L, tab8, saved8$centres)
log_line("Targeted n=8 recalibration: fresh bank 2 of 2")
c8_2 <- targeted_calibration_bank(2L, tab8, saved8$centres)
selected8 <- c(mean(c(c8_1$q$threshold, c8_2$q$threshold)),
               mean(c(c8_1$p$threshold, c8_2$p$threshold)))
names(selected8) <- TARGET_CHARTS

calibration_bank_rows <- do.call(rbind, list(
  data.frame(bank = 1L, chart = TARGET_CHARTS[1L],
             threshold = c8_1$q$threshold, c8_1$q$summary),
  data.frame(bank = 1L, chart = TARGET_CHARTS[2L],
             threshold = c8_1$p$threshold, c8_1$p$summary),
  data.frame(bank = 2L, chart = TARGET_CHARTS[1L],
             threshold = c8_2$q$threshold, c8_2$q$summary),
  data.frame(bank = 2L, chart = TARGET_CHARTS[2L],
             threshold = c8_2$p$threshold, c8_2$p$summary)
))
write_csv(calibration_bank_rows, "targeted_n8_calibration_banks.csv")
write_csv(data.frame(chart = TARGET_CHARTS, threshold = selected8,
                     selection = "arithmetic mean of two independent-bank limits",
                     n_per_bank = PP$n8_cal_per_bank,
                     horizon_per_bank = PP$n8_cal_horizon),
          "targeted_n8_selected_thresholds.csv")

log_line("Targeted n=8 third-bank IC validation (N=",
         PP$n8_final_validation, ")")
set.seed(PP$seed + 8903L)
T8_ic <- simulate_qp_paths(
  PP$n8_final_validation, "ic", 1L, PP$n8_final_horizon,
  tab8, saved8$centres, selected8)
ic8 <- summarise_target_ic(T8_ic, PP$n8_final_horizon)
write_csv(ic8, "targeted_n8_ic_validation.csv")

log_line("Targeted n=8 replacement OOC runs")
ooc8_rows <- list()
ooc8_replicates <- list(IC = T8_ic)
for (sc_i in seq_along(CORE$scenarios)) {
  sc <- CORE$scenarios[sc_i]
  for (nu in CORE$change_times) {
    N <- if (nu == 1L) CORE$n_ooc else CORE$n_delayed
    set.seed(PP$seed + 9000L + 100L * sc_i + nu)
    z <- simulate_qp_paths(N, sc, nu, CORE$max_delay, tab8,
                           saved8$centres, selected8)
    nm <- paste(sc, nu, sep = "_")
    ooc8_rows[[nm]] <- summarise_target_ooc(sc, nu, z, CORE$max_delay)
    ooc8_replicates[[nm]] <- z
  }
}
ooc8 <- do.call(rbind, ooc8_rows)
write_csv(ooc8, "targeted_n8_ooc_delay_metrics.csv")
saveRDS(ooc8_replicates, file.path(OUT, "targeted_n8_replicates.rds"),
        compress = TRUE)

# Preserve the original core files and create explicitly named final files in
# which only the two failed n=8 rows (and their OOC counterparts) are replaced.
core_thresholds <- read.csv(file.path(OUT, "calibrated_thresholds.csv"),
                            check.names = FALSE)
core_thresholds$integer_boundary <- ifelse(
  core_thresholds$chart == "Bakir Shew-KS",
  round(core_thresholds$threshold * CORE$m * core_thresholds$n), NA_real_)
new_threshold_rows <- list()
for (j in 1:2) {
  row <- core_thresholds[core_thresholds$n == 8L &
                         core_thresholds$chart == TARGET_CHARTS[j], , drop = FALSE]
  row$threshold <- selected8[j]
  row$calibration_horizon <- PP$n8_cal_horizon
  # The selected limit is the average of two independently calibrated limits.
  # Bank-specific run lengths were evaluated at their own limits, not at this
  # average, so reporting their pooled summary here would be misleading.  The
  # selected limit's operating performance is reported only from the untouched
  # third-bank validation below.
  summary_names <- intersect(names(rl_summary(rep(1, 2), 1L)), names(row))
  for (nm in summary_names) row[[nm]] <- NA_real_
  new_threshold_rows[[j]] <- row
}
new_threshold_rows <- do.call(rbind, new_threshold_rows)
final_thresholds <- rbind(
  transform(core_thresholds[!(core_thresholds$n == 8L &
    core_thresholds$chart %in% TARGET_CHARTS), , drop = FALSE],
    result_source = "original full core"),
  transform(new_threshold_rows, result_source = "two fresh calibration banks")
)
final_thresholds <- final_thresholds[order(final_thresholds$n,
  match(final_thresholds$chart, CHARTS)), ]
write_csv(final_thresholds, "final_calibrated_thresholds.csv")

core_ic <- read.csv(file.path(OUT, "ic_validation.csv"), check.names = FALSE)
stopifnot(setequal(names(core_ic), names(ic8)))
ic8 <- ic8[, names(core_ic), drop = FALSE]
final_ic <- rbind(
  transform(core_ic[!(core_ic$n == 8L & core_ic$chart %in% TARGET_CHARTS),
                    , drop = FALSE], result_source = "original full core"),
  transform(ic8, result_source = "independent third bank")
)
final_ic <- final_ic[order(final_ic$n, match(final_ic$chart, CHARTS)), ]
write_csv(final_ic, "final_ic_validation.csv")

core_ooc <- read.csv(file.path(OUT, "ooc_delay_metrics.csv"),
                     check.names = FALSE)
stopifnot(setequal(names(core_ooc), names(ooc8)))
ooc8 <- ooc8[, names(core_ooc), drop = FALSE]
final_ooc <- rbind(
  transform(core_ooc[!(core_ooc$n == 8L &
    core_ooc$chart %in% TARGET_CHARTS), , drop = FALSE],
    result_source = "original full core"),
  transform(ooc8, result_source = "rerun at targeted final limit")
)
final_ooc <- final_ooc[order(final_ooc$n,
  match(final_ooc$scenario, CORE$scenarios), final_ooc$change_time,
  match(final_ooc$chart, CHARTS)), ]
write_csv(final_ooc, "final_ooc_delay_metrics.csv")

acceptance8 <- transform(ic8,
  accepted = within_5_percent & CI_contains_target & censored_rate <= .002)
write_csv(acceptance8, "targeted_n8_acceptance_audit.csv")
if (!all(acceptance8$accepted)) {
  warning("At least one targeted n=8 final limit failed the prespecified audit; ",
          "inspect targeted_n8_acceptance_audit.csv before manuscript use.")
}
log_line("Targeted n=8 final files written")

  log_line("Target-only mode completed")
  quit(save = "no", status = 0L)
}

# In the normal workflow the large common-core run is already the frozen final
# analysis.  Copy its three summary tables without replacing individual cells.
legacy_targeted <- list.files(
  OUT, pattern = "^(targeted_n8_|TARGETED_N8_)", full.names = TRUE)
if (length(legacy_targeted)) unlink(legacy_targeted, force = TRUE)
core_thresholds <- read.csv(file.path(OUT, "calibrated_thresholds.csv"),
                            check.names = FALSE)
core_ic <- read.csv(file.path(OUT, "ic_validation.csv"), check.names = FALSE)
core_ooc <- read.csv(file.path(OUT, "ooc_delay_metrics.csv"),
                     check.names = FALSE)
write_csv(core_thresholds, "final_calibrated_thresholds.csv")
write_csv(core_ic, "final_ic_validation.csv")
write_csv(core_ooc, "final_ooc_delay_metrics.csv")
log_line("Common-core final summary files written")

log_line("Calibrating the n=5 two-interval VSI warning threshold")
set.seed(PP$seed + 5101L)
vsi_cal_paths <- simulate_raw_histories(
  PP$n_cal_vsi, "ic", TAB, action_threshold, action_gamma,
  PP$max_inspections, PP$chunk)
vsi_cal <- calibrate_warning_threshold(vsi_cal_paths, action_threshold)
write_csv(vsi_cal, "vsi_calibration.csv")

log_line("Independent VSI validation and OOC comparison")
vsi_rows <- list()
vsi_replicates <- list()
scenarios <- c("ic", "normal_location", "normal_scale", "cauchy")
for (j in seq_along(scenarios)) {
  sc <- scenarios[j]
  set.seed(PP$seed + 5200L + j)
  z <- simulate_raw_histories(
    PP$n_verify_vsi, sc, TAB, action_threshold, action_gamma,
    PP$max_inspections, PP$chunk)
  vsi_rows[[sc]] <- summarise_vsi(z, sc, vsi_cal$warning_threshold,
                                  vsi_cal$d_long)
  vsi_replicates[[sc]] <- list(
    T = z$T,
    VSI_time = calendar_times(z, vsi_cal$warning_threshold,
                              PP$d_short, vsi_cal$d_long,
                              PP$first_interval),
    censored = z$censored)
}
write_csv(do.call(rbind, vsi_rows), "vsi_metrics.csv")
saveRDS(vsi_replicates, file.path(OUT, "vsi_replicates.rds"), compress = TRUE)

log_line("Creating deterministic implementation example")
example <- make_worked_example(TAB)
write_csv(example$phase1, "worked_example_phaseI.csv")
write_csv(example$monitoring, "worked_example_monitoring_full.csv")

design <- data.frame(
  field = c("Phase I law", "Phase I m", "Phase II n", "exact p-value",
            "IC target", "location alternative", "scale alternative",
            "distribution alternative", "VSI action", "VSI warning",
            "VSI intervals", "elapsed-time convention"),
  value = c("N(0,1)", CORE$m, PP$n, "exact two-sample KS, continuous/tie-free",
            PP$target_ATS, "N(0.5,1)", "N(0,2) (variance 2)",
            "standard Cauchy", "raw exact-KS action boundary",
            format(vsi_cal$warning_threshold, digits = 12),
            paste(PP$d_short, format(vsi_cal$d_long, digits = 12), sep = "/"),
            "first inspection at time 1; next interval chosen from current p"))
write_csv(design, "design_specification.csv")

log_line("Postprocess completed.")
