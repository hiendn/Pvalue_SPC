#!/usr/bin/env Rscript

# Build the supplementary tables cited in the revised manuscript from the
# committed numerical summaries. Only base R is required.

options(stringsAsFactors = FALSE)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- normalizePath(sub("^--file=", "", script_arg[[1L]]),
                             mustWork = TRUE)
root <- dirname(dirname(script_file))
data_dir <- file.path(root, "data")
out_dir <- file.path(root, "supplementary_tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_data <- function(name) {
  read.csv(file.path(data_dir, name), check.names = FALSE,
           stringsAsFactors = FALSE)
}

write_table <- function(x, name) {
  write.csv(x, file.path(out_dir, name), row.names = FALSE, na = "")
}

charts <- c("Raw exact KS p", "Qtilde(0.95,1)", "p-EWMA(0.2)",
            "KS-EWMA(0.2)", "Bakir Shew-KS")
chart_order <- function(x) match(x, charts)

# S1a: IC calibration and independent validation.
cal <- read_data("final_calibrated_thresholds.csv")
val <- read_data("final_ic_validation.csv")
cal <- cal[cal$chart %in% charts, ]
val <- val[val$chart %in% charts, ]
z <- merge(cal, val, by = c("n", "chart"), suffixes = c("_cal", "_val"))
z <- z[order(z$n, chart_order(z$chart)), ]
s1a <- data.frame(
  n = z$n,
  chart = z$chart,
  direction = z$direction,
  limit = z$reported_limit,
  boundary_gamma = z$boundary_gamma,
  calibration_n_rep = z$n_rep_cal,
  calibration_horizon = z$calibration_horizon,
  calibration_mean_recorded_RL = z$mean_RL_cal,
  calibration_se_RL = z$se_RL_cal,
  calibration_ci_lower = z$ci_lower_cal,
  calibration_ci_upper = z$ci_upper_cal,
  calibration_relative_MCSE = z$relative_MCSE_cal,
  calibration_censored_rate = z$censored_rate_cal,
  deterministic_gamma0_mean_recorded_RL = z$deterministic_gamma0_mean,
  deterministic_gamma1_mean_recorded_RL = z$deterministic_gamma1_mean,
  target_ARL = z$target_ARL,
  validation_n_rep = z$n_rep_val,
  validation_mean_recorded_RL = z$mean_RL_val,
  validation_se_RL = z$se_RL_val,
  validation_ci_lower = z$ci_lower_val,
  validation_ci_upper = z$ci_upper_val,
  validation_relative_MCSE = z$relative_MCSE_val,
  validation_censored_rate = z$censored_rate_val
)
stopifnot(nrow(s1a) == 20L,
          sum(s1a$validation_censored_rate > 0) == 1L,
          s1a$n[s1a$validation_censored_rate > 0] == 8L,
          s1a$chart[s1a$validation_censored_rate > 0] == "Raw exact KS p")
write_table(s1a, "Table_S1a_IC_calibration_and_validation.csv")

# S1b: complete five-chart OOC comparison.
ooc <- read_data("final_ooc_delay_metrics.csv")
ooc <- ooc[ooc$chart %in% charts, ]
ooc <- ooc[order(ooc$n, ooc$scenario, chart_order(ooc$chart),
                 ooc$change_time), ]
s1b_cols <- c(
  "n", "scenario", "change_time", "chart", "n_started",
  "n_conditioned", "prechange_signal_rate", "mean_delay", "sd_delay",
  "se_delay", "delay_ci_lower", "delay_ci_upper", "delay_censored_rate",
  "delay_estimand", "mean_postchange_observations", "PTS_10",
  "PTS_10_lower", "PTS_10_upper", "PTS_25", "PTS_25_lower",
  "PTS_25_upper", "PTS_obs_50", "PTS_obs_100", "PTS_obs_200"
)
s1b <- ooc[, s1b_cols]
stopifnot(nrow(s1b) == 120L, all(s1b$delay_censored_rate == 0))
write_table(s1b, "Table_S1b_full_OOC_performance.csv")

# S2: VSI calibration and independent performance evaluation.
vcal <- read_data("vsi_calibration.csv")
vsi <- read_data("vsi_metrics.csv")
stopifnot(nrow(vcal) == 1L, nrow(vsi) == 8L)
s2 <- data.frame(
  n = vsi$n,
  scenario = vsi$scenario,
  method = vsi$method,
  estimand = vsi$estimand,
  n_rep = vsi$n_rep,
  calibration_n_rep = vcal$n_calibration,
  calibration_action_threshold = vcal$action_threshold,
  calibration_boundary_gamma = vcal$action_boundary_gamma,
  calibration_warning_threshold = vcal$warning_threshold,
  d_short = vcal$d_short,
  d_long = vcal$d_long,
  first_interval = vcal$first_interval,
  mean_inspections = vsi$mean_inspections,
  mean_elapsed_time = vsi$mean_elapsed_time,
  sd_elapsed_time = vsi$sd_elapsed_time,
  se_elapsed_time = vsi$se_elapsed_time,
  ci_lower = vsi$ci_lower,
  ci_upper = vsi$ci_upper,
  mean_interval_to_signal = vsi$mean_interval_to_signal,
  inspections_per_time_unit = vsi$inspections_per_time_unit,
  censoring_rate = vsi$censoring_rate,
  P_signal_by_time_10 = vsi$P_signal_by_time_10,
  P10_lower = vsi$P10_lower,
  P10_upper = vsi$P10_upper,
  P_signal_by_time_25 = vsi$P_signal_by_time_25,
  P25_lower = vsi$P25_lower,
  P25_upper = vsi$P25_upper
)
write_table(s2, "Table_S2_VSI_performance.csv")

# S3: the nine inspections displayed in the implementation example.
trace <- read_data("worked_example_monitoring_full.csv")
s3 <- trace[trace$inspection <= 9L,
            c("inspection", "regime", "phase2_values", "phase2_mean",
              "KS_D", "exact_KS_p", "Qtilde")]
stopifnot(nrow(s3) == 9L)
write_table(s3, "Table_S3_implementation_example.csv")

# S4: elementary theorem-bound studies.
e1 <- read_data("elementary_normal_ic.csv")
ar1 <- read_data("ar1_ic.csv")
e1_out <- data.frame(
  setting = "E1", study = "two-phase reusable baseline", chart = e1$chart,
  beta = NA_real_, alpha = e1$alpha, validity = e1$validity, k = e1$k,
  n_rep = e1$n_rep, mean_Rk = e1$mean_Rk, se_Rk = e1$se_Rk,
  ci_lower = e1$ci_lower, ci_upper = e1$ci_upper,
  relative_MCSE = e1$mcse_ratio, censored_rate = e1$censored_rate,
  marginal_or_conditional_lower_bound = e1$lower_bound,
  ratio_to_marginal_or_conditional_bound = e1$ratio,
  reusable_baseline_lower_bound = e1$reuse_lower_bound,
  ratio_to_reusable_baseline_bound = e1$reuse_ratio,
  applicable_lower_bound = e1$reuse_lower_bound,
  ratio_to_applicable_bound = e1$reuse_ratio
)
ar1_out <- data.frame(
  setting = ifelse(ar1$chart == "Pprime", "E2", "E3"),
  study = "stationary AR(1)", chart = ar1$chart, beta = ar1$beta,
  alpha = ar1$alpha, validity = ar1$validity, k = ar1$k,
  n_rep = ar1$n_rep, mean_Rk = ar1$mean_Rk, se_Rk = ar1$se_Rk,
  ci_lower = ar1$ci_lower, ci_upper = ar1$ci_upper,
  relative_MCSE = ar1$mcse_ratio, censored_rate = ar1$censored_rate,
  marginal_or_conditional_lower_bound = ar1$lower_bound,
  ratio_to_marginal_or_conditional_bound = ar1$ratio,
  reusable_baseline_lower_bound = NA_real_,
  ratio_to_reusable_baseline_bound = NA_real_,
  applicable_lower_bound = ar1$lower_bound,
  ratio_to_applicable_bound = ar1$ratio
)
s4 <- rbind(e1_out, ar1_out)
s4 <- s4[order(s4$setting, s4$beta, s4$alpha, s4$k), ]
stopifnot(nrow(s4) == 20L, all(s4$censored_rate == 0))
write_table(s4, "Table_S4_elementary_bound_results.csv")

# S5: localisation performance and stopped-run diagnostics.
loc <- rbind(read_data("localisation_normal.csv"),
             read_data("localisation_cauchy.csv"))
diag <- read_data("localisation_stopped_diagnostics.csv")
make_key <- function(x) {
  paste(x$model, ifelse(is.na(x$n0), "NA", x$n0), x$delta, x$rho,
        x$alpha, sep = "|")
}
idx <- match(make_key(loc), make_key(diag))
stopifnot(!anyNA(idx))
s5 <- data.frame(
  setting = ifelse(loc$model == "normal_Z", "L1", "L2"),
  model = loc$model, n0 = loc$n0, delta = loc$delta, rho = loc$rho,
  alpha = loc$alpha, n_rep_run = loc$n_rep_run, mean_R = loc$mean_R,
  se_R = loc$se_R, relative_MCSE = diag$relative_MCSE[idx],
  median_R = diag$median_R[idx], q90_R = diag$q90_R[idx],
  PTS_25 = diag$PTS_25[idx], PTS_100 = diag$PTS_100[idx],
  censored_rate = loc$censored_rate,
  mean_true_ooc = loc$mean_true_ooc, se_true_ooc = loc$se_true_ooc,
  mean_false_ooc_at_alarm = loc$mean_false_ooc_at_alarm,
  precision_proxy = loc$mean_true_ooc /
    (loc$mean_true_ooc + loc$mean_false_ooc_at_alarm),
  n_rep_fwer = loc$n_rep_fwer, fwer_hat = loc$fwer_hat,
  fwer_se = loc$fwer_se
)
s5 <- s5[order(s5$setting, s5$n0, s5$delta, s5$rho, s5$alpha), ]
stopifnot(nrow(s5) == 48L,
          sum(s5$censored_rate > 0) == 3L,
          max(s5$censored_rate) == 0.005)
write_table(s5, "Table_S5_localisation_results.csv")

message("Wrote Supplementary Tables S1a, S1b, and S2--S5 to ", out_dir)
