# KS matched-IC-ARL benchmark simulations.
#
# The script calibrates each chart to a common realised IC ARL and then evaluates
# IC verification and OOC delays. It is self-contained so that the benchmark can
# be run independently of the other simulation scripts.

message("Running KS matched-ARL benchmark.")

# ---- Defaults; override before sourcing if needed -----------------------------
set_default <- function(name, value, envir = .GlobalEnv) {
  if (!exists(name, envir = envir, inherits = FALSE)) assign(name, value, envir = envir)
  invisible(get(name, envir = envir, inherits = FALSE))
}

set_default("SMOKE_TEST", FALSE)
set_default("RESUME", TRUE)
set_default("FORCE_RECALIBRATE", FALSE)
set_default("FORCE_FINAL_IC", FALSE)
set_default("FORCE_OOC", FALSE)
set_default("START_AT_CHART", NULL)  # for example, "KS-CUSUM statistic"
set_default("RUN_STAGES", c("calibration", "ic", "ooc"))
set_default("OUTPUT_DIR", "output")
set_default("KS_ENGINE", "asymp")    # "asymp", "asymp_tie_safe", "hybrid", "stats_asymp", "stats_exact"
set_default("JMAX", 50L)
set_default("EXACT_MAX_PRODUCT", 15000L)
set_default("FINITE_SAMPLE_CORRECTION", FALSE)
set_default("N0_SET", c(50L, 100L))
set_default("TARGET_ARL", c(200, 500))
set_default("OOC_NAMES", c("N_0.5_1", "N_1_1", "N_0_2", "Cauchy",
                           "dyn_mu_half", "dyn_mu_quarter", "dyn_var_chi1", "dyn_var_chi2"))
set_default("N_REP_CAL", 200L)
set_default("N_REP_FINAL_IC", 1000L)
set_default("N_REP_OOC", 1000L)
set_default("N_REP_MEAN_D", 1000L)
set_default("MAX_T_CAL", 100000L)
set_default("MAX_T_FINAL", 200000L)
set_default("MAX_T_OOC", 100000L)
set_default("N_ITER_CAL", 15L)
set_default("CAL_HORIZON_MULTIPLIER", 12L)
set_default("CAL_HORIZON_MIN", 3000L)
set_default("N_CORES", max(1L, parallel::detectCores(logical = FALSE) - 1L))
set_default("BATCH_SIZE", 25L)
set_default("VSS_WIDTH", 10L)
set_default("VERBOSE", TRUE)
set_default("SAVE_RAW_RDS", FALSE)

if (isTRUE(SMOKE_TEST)) {
  N0_SET <- c(30L); TARGET_ARL <- c(30); OOC_NAMES <- c("N_0.5_1", "dyn_mu_half")
  N_REP_CAL <- 20L; N_REP_FINAL_IC <- 30L; N_REP_OOC <- 30L; N_REP_MEAN_D <- 40L
  MAX_T_CAL <- 1000L; MAX_T_FINAL <- 3000L; MAX_T_OOC <- 3000L; N_ITER_CAL <- 6L
  CAL_HORIZON_MIN <- 300L; CAL_HORIZON_MULTIPLIER <- 8L; N_CORES <- 1L; BATCH_SIZE <- 10L
  OUTPUT_DIR <- file.path(OUTPUT_DIR, "smoke_test")
}

# ---- Small utilities ----------------------------------------------------------
`%||%` <- function(x, y) if (is.null(x)) y else x
ensure_dir <- function(path) { if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE); invisible(path) }
ensure_dir(OUTPUT_DIR); ensure_dir(file.path(OUTPUT_DIR, "raw_rds"))
now_string <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
log_msg <- function(..., verbose = TRUE) { if (isTRUE(verbose)) { cat(sprintf("[%s] %s\n", now_string(), sprintf(...))); flush.console() }; invisible(NULL) }
format_elapsed <- function(seconds) { seconds <- as.numeric(seconds); if (!is.finite(seconds)) "NA" else if (seconds < 60) sprintf("%.1fs", seconds) else if (seconds < 3600) sprintf("%.1fmin", seconds / 60) else sprintf("%.2fh", seconds / 3600) }
se <- function(x) { x <- x[is.finite(x)]; if (length(x) <= 1L) NA_real_ else stats::sd(x) / sqrt(length(x)) }
mc_ci <- function(x, level = 0.95) { x <- x[is.finite(x)]; if (!length(x)) return(c(NA_real_, NA_real_)); z <- stats::qnorm(1 - (1 - level) / 2); m <- mean(x); s <- se(x); c(m - z*s, m + z*s) }
sample_vss <- function(n0, width = 10L) max(1L, as.integer(sample.int(2L * width + 1L, 1L) + n0 - width - 1L))
clamp01 <- function(x) { x <- as.numeric(x); x[!is.finite(x)] <- 1; pmin(pmax(x, 0), 1) }
safe_p <- function(p) { p <- as.numeric(p); p[!is.finite(p)] <- 1; pmin(pmax(p, .Machine$double.xmin), 1) }
sanitize_key <- function(x) gsub("[^A-Za-z0-9]+", "_", as.character(x))

# Small CSV writer for compact summary files.
csv_field <- function(x) {
  if (length(x) == 0L || is.na(x)) return("")
  if (is.numeric(x)) { if (!is.finite(x)) return(""); return(format(x, scientific = FALSE, digits = 15, trim = TRUE)) }
  if (is.logical(x)) return(ifelse(isTRUE(x), "TRUE", "FALSE"))
  s <- as.character(x); s <- gsub('"', '""', s, fixed = TRUE)
  if (grepl('[",\n\r]', s)) s <- paste0('"', s, '"')
  s
}
write_csv_manual <- function(df, file) {
  ensure_dir(dirname(file)); tmp <- paste0(file, ".tmp.", Sys.getpid(), ".", sample.int(1e9, 1))
  if (is.null(df) || !nrow(df)) {
    lines <- paste(names(df %||% data.frame()), collapse = ",")
  } else {
    header <- paste(vapply(names(df), csv_field, character(1L)), collapse = ",")
    body <- apply(df, 1L, function(row) paste(vapply(as.list(row), csv_field, character(1L)), collapse = ","))
    lines <- c(header, body)
  }
  writeLines(lines, tmp, useBytes = TRUE)
  if (file.exists(file)) unlink(file)
  if (!file.rename(tmp, file)) { file.copy(tmp, file, overwrite = TRUE); unlink(tmp) }
  invisible(file)
}
safe_save_rds <- function(object, file, compress = "xz") {
  ensure_dir(dirname(file)); tmp <- paste0(file, ".tmp.", Sys.getpid(), ".", sample.int(1e9, 1))
  saveRDS(object, tmp, compress = compress)
  if (file.exists(file)) unlink(file)
  if (!file.rename(tmp, file)) { file.copy(tmp, file, overwrite = TRUE); unlink(tmp) }
  invisible(file)
}

# ---- Fast KS statistic and p-value engines -----------------------------------
ks_D_findinterval <- function(x0_sorted, xt) {
  y <- sort(xt); n <- length(x0_sorted); m <- length(y)
  fy_at_x <- findInterval(x0_sorted, y, rightmost.closed = TRUE) / m
  fx_at_y <- findInterval(y, x0_sorted, rightmost.closed = TRUE) / n
  max(abs(seq_len(n) / n - fy_at_x), abs(fx_at_y - seq_len(m) / m))
}
ks_D_merge_tie_safe <- function(x0_sorted, xt) {
  y <- sort(xt); n <- length(x0_sorted); m <- length(y); i <- 1L; j <- 1L; c0 <- 0L; c1 <- 0L; D <- 0
  while (i <= n || j <= m) {
    if (j > m || (i <= n && x0_sorted[i] <= y[j])) z <- x0_sorted[i] else z <- y[j]
    while (i <= n && x0_sorted[i] <= z) { c0 <- c0 + 1L; i <- i + 1L }
    while (j <= m && y[j] <= z) { c1 <- c1 + 1L; j <- j + 1L }
    d <- abs(c0 / n - c1 / m); if (d > D) D <- d
  }
  D
}
kolmogorov_sf_series <- function(z, jmax = 50L) {
  z <- as.numeric(z); if (!is.finite(z) || z <= 0) return(1)
  j <- seq_len(jmax); clamp01(2 * sum(((-1)^(j - 1L)) * exp(-2 * (j^2) * (z^2))))
}
ks_p_asymp_from_D <- function(D, n0, nt, jmax = 50L, finite_sample_correction = FALSE) {
  ne <- (n0 * nt) / (n0 + nt); z <- sqrt(ne) * D
  if (isTRUE(finite_sample_correction)) z <- (sqrt(ne) + 0.12 + 0.11 / sqrt(ne)) * D
  kolmogorov_sf_series(z, jmax)
}
ks_p_and_D <- function(x0_sorted, xt, engine = c("asymp", "asymp_tie_safe", "hybrid", "stats_asymp", "stats_exact"), exact_max_product = 15000L, jmax = 50L, finite_sample_correction = FALSE) {
  engine <- match.arg(engine); n0 <- length(x0_sorted); nt <- length(xt)
  if (engine %in% c("stats_asymp", "stats_exact")) {
    kt <- suppressWarnings(stats::ks.test(x0_sorted, xt, exact = identical(engine, "stats_exact")))
    return(list(D = unname(as.numeric(kt$statistic)), p = clamp01(kt$p.value)))
  }
  D <- if (identical(engine, "asymp_tie_safe")) ks_D_merge_tie_safe(x0_sorted, xt) else ks_D_findinterval(x0_sorted, xt)
  if (identical(engine, "hybrid") && n0 * nt <= exact_max_product) {
    kt <- tryCatch(suppressWarnings(stats::ks.test(x0_sorted, xt, exact = TRUE)), error = function(e) suppressWarnings(stats::ks.test(x0_sorted, xt, exact = FALSE)))
    return(list(D = unname(as.numeric(kt$statistic)), p = clamp01(kt$p.value)))
  }
  list(D = D, p = ks_p_asymp_from_D(D, n0, nt, jmax, finite_sample_correction))
}

# ---- Data generators ----------------------------------------------------------
phase1_normal <- function(n) stats::rnorm(n, 0, 1)
make_phase2_generator <- function(name) {
  force(name)
  function(n, t) switch(name,
    IC = stats::rnorm(n, 0, 1),
    N_0.5_1 = stats::rnorm(n, 0.5, 1),
    N_1_1 = stats::rnorm(n, 1, 1),
    N_0_2 = stats::rnorm(n, 0, sqrt(2)),
    Cauchy = stats::rcauchy(n, 0, 1),
    dyn_mu_half = { mu <- stats::rnorm(1, 0, sqrt(1/2)); stats::rnorm(n, mu, 1) },
    dyn_mu_quarter = { mu <- stats::rnorm(1, 0, sqrt(1/4)); stats::rnorm(n, mu, 1) },
    dyn_var_chi1 = { sig2 <- stats::rchisq(1, 1); stats::rnorm(n, 0, sqrt(sig2)) },
    dyn_var_chi2 = { sig2 <- stats::rchisq(1, 2); stats::rnorm(n, 0, sqrt(sig2)) },
    stop("Unknown phase-II generator name: ", name))
}
phase2_label <- function(name) switch(name,
  IC = "IC N(0,1)", N_0.5_1 = "N(1/2,1)", N_1_1 = "N(1,1)", N_0_2 = "N(0,2)", Cauchy = "Cauchy(0,1)",
  dyn_mu_half = "N(mu_t,1), mu_t~N(0,1/2)", dyn_mu_quarter = "N(mu_t,1), mu_t~N(0,1/4)",
  dyn_var_chi1 = "N(0,sigma_t^2), sigma_t^2~chi^2_1", dyn_var_chi2 = "N(0,sigma_t^2), sigma_t^2~chi^2_2", name)

# ---- Chart update rules -------------------------------------------------------
ewma_multiplier_qtilde <- function(lambda, r) if (r >= 1) min(1 + r, 1 / lambda)^(1 / r) else (1 + r)^(1 / r)
ewma_multiplier_qbar <- function(lambda, r) lambda^(-1 / r)
make_chart_label <- function(ch) ch$label %||% switch(ch$type, p_shewhart = "p-value Shewhart", qtilde = "Qtilde", qbar = "Qbar", ks_shewhart = "KS Shewhart statistic", ks_ewma = "KS-EWMA statistic", ks_cusum = "KS-CUSUM statistic", ch$type)
chart_stat_update <- function(state, ch, p_value, ks_stat) {
  p_value <- safe_p(p_value)
  if (ch$type == "p_shewhart") stat <- p_value
  else if (ch$type == "qtilde") { S <- if (is.null(state$S)) p_value^ch$r else ch$lambda*p_value^ch$r + (1-ch$lambda)*state$S; state$S <- S; stat <- ewma_multiplier_qtilde(ch$lambda, ch$r) * S^(1/ch$r) }
  else if (ch$type == "qbar") { S <- if (is.null(state$S)) p_value^ch$r else ch$lambda*p_value^ch$r + (1-ch$lambda)*state$S; state$S <- S; stat <- ewma_multiplier_qbar(ch$lambda, ch$r) * S^(1/ch$r) }
  else if (ch$type == "ks_shewhart") stat <- ks_stat
  else if (ch$type == "ks_ewma") { W <- if (is.null(state$W)) ks_stat else ch$lambda*ks_stat + (1-ch$lambda)*state$W; state$W <- W; stat <- W }
  else if (ch$type == "ks_cusum") { C <- max(0, (state$C %||% 0) + ks_stat - ch$kappa); state$C <- C; stat <- C }
  else stop("Unknown chart type: ", ch$type)
  list(state = state, statistic = stat)
}
chart_param_frame <- function(charts) {
  do.call(rbind, lapply(seq_along(charts), function(i) {
    ch <- charts[[i]]
    data.frame(chart_index=i, chart=make_chart_label(ch), chart_type=ch$type %||% NA_character_, threshold=ch$threshold %||% NA_real_, direction=ch$direction %||% NA_character_, validity=ch$validity %||% NA_character_, lambda=ch$lambda %||% NA_real_, r=ch$r %||% NA_real_, kappa=ch$kappa %||% NA_real_, stringsAsFactors=FALSE)
  }))
}
# ---- Compact state and output -------------------------------------------------
state_file <- file.path(OUTPUT_DIR, "ks_benchmark_state.rds")
empty_thresholds <- data.frame(n0=integer(), target_arl=numeric(), chart=character(), chart_label=character(), chart_type=character(), threshold=numeric(), direction=character(), lambda=numeric(), r=numeric(), kappa=numeric(), ks_engine=character(), selected_mean_R1=numeric(), selected_se_R1=numeric(), selected_censored_rate=numeric(), stringsAsFactors=FALSE)
empty_cal_history <- data.frame(n0=integer(), target_arl=numeric(), chart=character(), iter=integer(), threshold=numeric(), mean_R1=numeric(), se_R1=numeric(), censored_rate=numeric(), stringsAsFactors=FALSE)
empty_ooc <- data.frame(n0=integer(), target_arl=numeric(), ooc=character(), ooc_label=character(), chart_index=integer(), chart=character(), chart_type=character(), threshold=numeric(), direction=character(), validity=character(), lambda=numeric(), r=numeric(), kappa=numeric(), k=integer(), n_rep=integer(), mean_Rk=numeric(), sd_Rk=numeric(), se_Rk=numeric(), ci_lower=numeric(), ci_upper=numeric(), mcse_ratio=numeric(), censored_rate=numeric(), ks_engine=character(), stringsAsFactors=FALSE)
empty_ic <- empty_ooc[0, setdiff(names(empty_ooc), c("ooc", "ooc_label")), drop=FALSE]
init_state <- function() list(version="ks-benchmark", thresholds=empty_thresholds, calibration_history=empty_cal_history, ic=empty_ic, ooc=empty_ooc, created_at=now_string(), updated_at=now_string())
coerce_to_template <- function(df, template) { if (is.null(df) || !nrow(df)) return(template[0,,drop=FALSE]); for (nm in setdiff(names(template), names(df))) df[[nm]] <- NA; rbind(template[0,,drop=FALSE], df[, names(template), drop=FALSE]) }
import_existing_csv_state <- function(state, output_dir) {
  f <- list(thresholds=file.path(output_dir,"ks_benchmark_thresholds_fast.csv"), calibration_history=file.path(output_dir,"ks_benchmark_calibration_history_fast.csv"), ic=file.path(output_dir,"ks_benchmark_ic_fast.csv"), ooc=file.path(output_dir,"ks_benchmark_ooc_fast.csv"))
  if (file.exists(f$thresholds)) { log_msg("Importing existing threshold CSV: %s", f$thresholds, verbose=VERBOSE); z <- tryCatch(utils::read.csv(f$thresholds, stringsAsFactors=FALSE), error=function(e) NULL); if (!is.null(z)) state$thresholds <- coerce_to_template(z, empty_thresholds) }
  if (file.exists(f$calibration_history)) { log_msg("Importing existing calibration CSV: %s", f$calibration_history, verbose=VERBOSE); z <- tryCatch(utils::read.csv(f$calibration_history, stringsAsFactors=FALSE), error=function(e) NULL); if (!is.null(z)) state$calibration_history <- coerce_to_template(z, empty_cal_history) }
  if (file.exists(f$ic)) { log_msg("Importing existing IC CSV: %s", f$ic, verbose=VERBOSE); z <- tryCatch(utils::read.csv(f$ic, stringsAsFactors=FALSE), error=function(e) NULL); if (!is.null(z)) state$ic <- coerce_to_template(z, empty_ic) }
  if (file.exists(f$ooc)) { log_msg("Importing existing OOC CSV: %s", f$ooc, verbose=VERBOSE); z <- tryCatch(utils::read.csv(f$ooc, stringsAsFactors=FALSE), error=function(e) NULL); if (!is.null(z)) state$ooc <- coerce_to_template(z, empty_ooc) }
  state
}
state <- if (isTRUE(RESUME) && file.exists(state_file)) { log_msg("Loading compact state: %s", state_file, verbose=VERBOSE); readRDS(state_file) } else { st <- init_state(); if (isTRUE(RESUME)) st <- import_existing_csv_state(st, OUTPUT_DIR); st }
save_state <- function(state) { state$updated_at <- now_string(); safe_save_rds(state, state_file) }
write_public_outputs <- function(state) { write_csv_manual(state$thresholds, file.path(OUTPUT_DIR,"ks_benchmark_thresholds_fast.csv")); write_csv_manual(state$calibration_history, file.path(OUTPUT_DIR,"ks_benchmark_calibration_history_fast.csv")); write_csv_manual(state$ic, file.path(OUTPUT_DIR,"ks_benchmark_ic_fast.csv")); write_csv_manual(state$ooc, file.path(OUTPUT_DIR,"ks_benchmark_ooc_fast.csv")); invisible(NULL) }
state_has_threshold <- function(state, n0, target, chart) { z <- state$thresholds; nrow(z)>0 && any(z$n0==n0 & z$target_arl==target & z$chart==chart) }
state_has_ic <- function(state, n0, target) { z <- state$ic; nrow(z)>0 && any(z$n0==n0 & z$target_arl==target) }
state_has_ooc <- function(state, n0, target, ooc) { z <- state$ooc; nrow(z)>0 && any(z$n0==n0 & z$target_arl==target & z$ooc==ooc) }
remove_threshold <- function(state, n0, target, chart) { state$thresholds <- state$thresholds[!(state$thresholds$n0==n0 & state$thresholds$target_arl==target & state$thresholds$chart==chart),,drop=FALSE]; state$calibration_history <- state$calibration_history[!(state$calibration_history$n0==n0 & state$calibration_history$target_arl==target & state$calibration_history$chart==chart),,drop=FALSE]; state }
remove_ic <- function(state, n0, target) { state$ic <- state$ic[!(state$ic$n0==n0 & state$ic$target_arl==target),,drop=FALSE]; state }
remove_ooc <- function(state, n0, target, ooc) { state$ooc <- state$ooc[!(state$ooc$n0==n0 & state$ooc$target_arl==target & state$ooc$ooc==ooc),,drop=FALSE]; state }

# ---- Templates and CUSUM reference -------------------------------------------
estimate_ic_ks_stat_mean <- function(n_rep, n0, seed, verbose=TRUE) {
  D <- numeric(n_rep); t0 <- Sys.time()
  for (i in seq_len(n_rep)) {
    set.seed(seed+i); x0 <- sort(phase1_normal(n0)); nt <- sample_vss(n0, VSS_WIDTH); xt <- phase1_normal(nt)
    D[i] <- ks_p_and_D(x0, xt, engine=KS_ENGINE, exact_max_product=EXACT_MAX_PRODUCT, jmax=JMAX, finite_sample_correction=FINITE_SAMPLE_CORRECTION)$D
    if (isTRUE(verbose) && (i==1L || i %% max(1L, floor(n_rep/10L)) == 0L)) log_msg("Estimating mean KS D: n0=%d rep=%d/%d mean=%.5f elapsed=%s", n0, i, n_rep, mean(D[seq_len(i)]), format_elapsed(difftime(Sys.time(), t0, units="secs")))
  }
  mean(D)
}
get_existing_cusum_kappa <- function(state, n0) { z <- state$thresholds; zz <- z[z$n0==n0 & z$chart=="KS-CUSUM statistic" & is.finite(z$kappa),,drop=FALSE]; if (nrow(zz)) zz$kappa[1] else NULL }
make_templates <- function(n0, state=NULL) {
  kappa <- if (!is.null(state)) get_existing_cusum_kappa(state, n0) else NULL
  if (is.null(kappa)) { log_msg("Estimating CUSUM reference kappa for n0=%d", n0, verbose=VERBOSE); kappa <- 1.05 * estimate_ic_ks_stat_mean(N_REP_MEAN_D, n0, 2026L+n0, VERBOSE) } else log_msg("Using existing CUSUM kappa=%.6f for n0=%d", kappa, n0, verbose=VERBOSE)
  list(
    list(name="p-value Shewhart", direction="lower", lower=1e-7, upper=0.50, log_scale=TRUE, chart=list(type="p_shewhart", validity="marginal", label="p-value Shewhart")),
    list(name="Qtilde r=-0.8", direction="lower", lower=1e-9, upper=0.50, log_scale=TRUE, chart=list(type="qtilde", lambda=0.5, r=-0.8, validity="marginal", label="Qtilde r=-0.8")),
    list(name="Qtilde r=-0.9", direction="lower", lower=1e-10, upper=0.50, log_scale=TRUE, chart=list(type="qtilde", lambda=0.5, r=-0.9, validity="marginal", label="Qtilde r=-0.9")),
    list(name="Qbar lambda=0.95", direction="lower", lower=1e-7, upper=0.50, log_scale=TRUE, chart=list(type="qbar", lambda=0.95, r=1, validity="marginal", label="Qbar lambda=0.95")),
    list(name="KS Shewhart statistic", direction="upper", lower=0.001, upper=1.0, log_scale=FALSE, chart=list(type="ks_shewhart", label="KS Shewhart statistic")),
    list(name="KS-EWMA statistic", direction="upper", lower=0.001, upper=1.0, log_scale=FALSE, chart=list(type="ks_ewma", lambda=0.2, label="KS-EWMA statistic")),
    list(name="KS-CUSUM statistic", direction="upper", lower=0.001, upper=5.0, log_scale=FALSE, chart=list(type="ks_cusum", kappa=kappa, label="KS-CUSUM statistic"))
  )
}

# ---- Calibration from reusable paths -----------------------------------------
calibration_horizon <- function() as.integer(min(MAX_T_CAL, max(CAL_HORIZON_MIN, ceiling(CAL_HORIZON_MULTIPLIER * max(TARGET_ARL)))))
build_calibration_paths <- function(templates, n_rep, n0, T_len, seed, verbose=TRUE) {
  charts <- lapply(templates, `[[`, "chart"); n_charts <- length(charts)
  paths <- lapply(seq_len(n_charts), function(i) matrix(NA_real_, n_rep, T_len)); names(paths) <- vapply(templates, function(x) x$name, character(1L))
  t0 <- Sys.time(); log_msg("Building reusable IC calibration paths: n0=%d reps=%d T=%d charts=%d", n0, n_rep, T_len, n_charts, verbose=verbose)
  for (rep in seq_len(n_rep)) {
    set.seed(seed+rep); x0 <- sort(phase1_normal(n0)); states <- vector("list", n_charts); states[] <- list(list())
    for (tt in seq_len(T_len)) {
      nt <- sample_vss(n0, VSS_WIDTH); xt <- phase1_normal(nt)
      kp <- ks_p_and_D(x0, xt, engine=KS_ENGINE, exact_max_product=EXACT_MAX_PRODUCT, jmax=JMAX, finite_sample_correction=FINITE_SAMPLE_CORRECTION)
      for (j in seq_len(n_charts)) { upd <- chart_stat_update(states[[j]], charts[[j]], kp$p, kp$D); states[[j]] <- upd$state; paths[[j]][rep, tt] <- upd$statistic }
    }
    if (isTRUE(verbose) && (rep==1L || rep %% max(1L, floor(n_rep/10L)) == 0L)) log_msg("Calibration path progress: n0=%d rep=%d/%d elapsed=%s", n0, rep, n_rep, format_elapsed(difftime(Sys.time(), t0, units="secs")))
  }
  paths
}
first_crossing_vec <- function(path_mat, threshold, direction) {
  hit <- if (direction == "lower") path_mat <= threshold else path_mat >= threshold
  has <- rowSums(hit, na.rm=TRUE) > 0L; idx <- max.col(hit, ties.method="first"); idx[!has] <- ncol(path_mat) + 1L; as.integer(idx)
}
calibrate_from_path <- function(path_mat, target_arl, direction, lower, upper, log_scale, n_iter, chart_name) {
  lo <- lower; hi <- upper; T_len <- ncol(path_mat)
  hist <- data.frame(iter=integer(), threshold=numeric(), mean_R1=numeric(), se_R1=numeric(), censored_rate=numeric())
  for (iter in seq_len(n_iter)) {
    mid <- if (isTRUE(log_scale)) exp((log(lo)+log(hi))/2) else (lo+hi)/2
    R1 <- first_crossing_vec(path_mat, mid, direction); cens <- R1 > T_len
    m <- mean(R1); s <- se(R1); cr <- mean(cens)
    hist <- rbind(hist, data.frame(iter=iter, threshold=mid, mean_R1=m, se_R1=s, censored_rate=cr))
    log_msg("Calibration %s target=%g iter=%02d/%02d h=%.8g mean_R1=%.2f se=%.2f censor=%.3f bracket=[%.6g, %.6g]", chart_name, target_arl, iter, n_iter, mid, m, s, cr, lo, hi, verbose=VERBOSE)
    if (direction == "lower") { if (m > target_arl) lo <- mid else hi <- mid } else { if (m > target_arl) hi <- mid else lo <- mid }
  }
  best <- hist[which.min(abs(log(hist$mean_R1 / target_arl))), , drop=FALSE]
  list(threshold=best$threshold[1], history=hist, selected=best)
}
# ---- Final IC/OOC stream simulation ------------------------------------------
run_stream_final <- function(charts, n0, phase2_name, k_values, max_t, seed) {
  set.seed(seed); x0 <- sort(phase1_normal(n0)); phase2_gen <- make_phase2_generator(phase2_name)
  n_charts <- length(charts); states <- vector("list", n_charts); states[] <- list(list())
  params <- chart_param_frame(charts); k_values <- sort(unique(as.integer(k_values))); k_max <- max(k_values)
  Rmat <- matrix(NA_integer_, n_charts, length(k_values)); alarm_counts <- integer(n_charts); finished <- rep(FALSE, n_charts)
  tt <- 0L
  while (tt < max_t && !all(finished)) {
    tt <- tt + 1L; nt <- sample_vss(n0, VSS_WIDTH); xt <- phase2_gen(nt, tt)
    kp <- ks_p_and_D(x0, xt, engine=KS_ENGINE, exact_max_product=EXACT_MAX_PRODUCT, jmax=JMAX, finite_sample_correction=FINITE_SAMPLE_CORRECTION)
    for (j in seq_len(n_charts)) {
      if (finished[j]) next
      upd <- chart_stat_update(states[[j]], charts[[j]], kp$p, kp$D); states[[j]] <- upd$state; ch <- charts[[j]]
      alarm <- if (ch$direction == "lower") upd$statistic <= ch$threshold else upd$statistic >= ch$threshold
      if (isTRUE(alarm)) {
        alarm_counts[j] <- alarm_counts[j] + 1L; hit <- which(k_values == alarm_counts[j])
        if (length(hit) == 1L) Rmat[j, hit] <- tt
        if (alarm_counts[j] >= k_max) finished[j] <- TRUE
      }
    }
  }
  rows <- vector("list", n_charts * length(k_values)); idx <- 1L
  for (j in seq_len(n_charts)) for (kk_i in seq_along(k_values)) {
    cens <- is.na(Rmat[j, kk_i])
    rows[[idx]] <- data.frame(chart_index=j, chart=params$chart[j], chart_type=params$chart_type[j], threshold=params$threshold[j], direction=params$direction[j], validity=params$validity[j], lambda=params$lambda[j], r=params$r[j], kappa=params$kappa[j], k=k_values[kk_i], Rk=if (cens) max_t+1L else Rmat[j, kk_i], censored=cens, stopped_at=tt, stringsAsFactors=FALSE)
    idx <- idx + 1L
  }
  do.call(rbind, rows)
}
summarise_run_df <- function(df) {
  if (!nrow(df)) return(data.frame())
  pieces <- split(df, interaction(df$chart_index, df$k, drop=TRUE))
  out <- do.call(rbind, lapply(pieces, function(z) {
    R <- z$Rk[is.finite(z$Rk)]; ci <- mc_ci(R)
    data.frame(chart_index=unique(z$chart_index), chart=unique(z$chart), chart_type=unique(z$chart_type), threshold=unique(z$threshold), direction=unique(z$direction), validity=unique(z$validity), lambda=unique(z$lambda), r=unique(z$r), kappa=unique(z$kappa), k=unique(z$k), n_rep=nrow(z), mean_Rk=mean(R), sd_Rk=stats::sd(R), se_Rk=se(R), ci_lower=ci[1], ci_upper=ci[2], mcse_ratio=ifelse(mean(R)>0, se(R)/mean(R), NA_real_), censored_rate=mean(z$censored), stringsAsFactors=FALSE)
  }))
  rownames(out) <- NULL; out[order(out$chart_index, out$k), , drop=FALSE]
}
print_partial_summary <- function(raw_df, label) {
  if (!isTRUE(VERBOSE) || !nrow(raw_df)) return(invisible(NULL))
  sm <- summarise_run_df(raw_df); z <- sm[sm$k == min(sm$k), c("chart", "mean_Rk", "se_Rk", "censored_rate"), drop=FALSE]
  log_msg("Partial summary %s, k=%d", label, min(sm$k)); print(z, row.names=FALSE); flush.console()
}
run_replicates_final <- function(n_rep, seed, charts, n0, phase2_name, k_values, max_t) {
  batches <- split(seq_len(n_rep), ceiling(seq_len(n_rep) / BATCH_SIZE)); raw_list <- vector("list", 0L); t0 <- Sys.time()
  log_msg("Starting final runs: phase2=%s n0=%d reps=%d charts=%d max_t=%d cores=%d", phase2_name, n0, n_rep, length(charts), max_t, N_CORES, verbose=VERBOSE)
  for (b in seq_along(batches)) {
    ids <- batches[[b]]
    worker <- function(i) { out <- run_stream_final(charts, n0, phase2_name, k_values, max_t, seed+i); out$replicate <- i; out }
    ans <- if (N_CORES > 1L && .Platform$OS.type != "windows") parallel::mclapply(ids, worker, mc.cores=N_CORES) else lapply(ids, worker)
    raw_list <- c(raw_list, ans); raw_df <- do.call(rbind, raw_list)
    log_msg("Finished batch %d/%d; reps done=%d/%d; elapsed=%s", b, length(batches), length(raw_list), n_rep, format_elapsed(difftime(Sys.time(), t0, units="secs")), verbose=VERBOSE)
    print_partial_summary(raw_df, sprintf("after %d reps", length(raw_list)))
  }
  raw_df <- do.call(rbind, raw_list); rownames(raw_df) <- NULL; raw_df
}

# ---- Reconstruct chart objects from threshold rows ----------------------------
threshold_row_to_chart <- function(row) {
  ch <- list(type=as.character(row$chart_type), threshold=as.numeric(row$threshold), direction=as.character(row$direction), label=as.character(row$chart), validity=as.character(row$validity %||% NA_character_))
  if (is.finite(as.numeric(row$lambda))) ch$lambda <- as.numeric(row$lambda)
  if (is.finite(as.numeric(row$r))) ch$r <- as.numeric(row$r)
  if (is.finite(as.numeric(row$kappa))) ch$kappa <- as.numeric(row$kappa)
  if (ch$type %in% c("p_shewhart", "qtilde", "qbar")) ch$validity <- "marginal"
  ch
}
get_charts_from_state <- function(state, n0, target) {
  z <- state$thresholds[state$thresholds$n0==n0 & state$thresholds$target_arl==target,,drop=FALSE]
  if (!nrow(z)) return(list())
  desired <- c("p-value Shewhart", "Qtilde r=-0.8", "Qtilde r=-0.9", "Qbar lambda=0.95", "KS Shewhart statistic", "KS-EWMA statistic", "KS-CUSUM statistic")
  z$order <- match(z$chart, desired); z <- z[order(z$order),,drop=FALSE]
  lapply(seq_len(nrow(z)), function(i) threshold_row_to_chart(z[i,]))
}

# ---- Main execution -----------------------------------------------------------
log_msg("Settings: RESUME=%s SMOKE_TEST=%s KS_ENGINE=%s OUTPUT_DIR=%s N_CORES=%d", RESUME, SMOKE_TEST, KS_ENGINE, OUTPUT_DIR, N_CORES, verbose=VERBOSE)
log_msg("Compact state file: %s", state_file, verbose=VERBOSE)

if ("calibration" %in% RUN_STAGES) {
  for (n0 in N0_SET) {
    templates <- make_templates(n0, state=state)
    missing <- FALSE; started <- is.null(START_AT_CHART)
    for (target in TARGET_ARL) for (tpl in templates) {
      if (!started) { if (identical(tpl$name, START_AT_CHART)) started <- TRUE else next }
      if (isTRUE(FORCE_RECALIBRATE) || !state_has_threshold(state, n0, target, tpl$name)) missing <- TRUE
    }
    if (!missing) { log_msg("All calibration thresholds already present for n0=%d; skipping.", n0, verbose=VERBOSE); next }
    paths <- build_calibration_paths(templates, N_REP_CAL, n0, calibration_horizon(), 30000L+n0, VERBOSE)
    started <- is.null(START_AT_CHART)
    for (target in TARGET_ARL) for (tpl_i in seq_along(templates)) {
      tpl <- templates[[tpl_i]]
      if (!started) { if (identical(tpl$name, START_AT_CHART)) started <- TRUE else { log_msg("START_AT_CHART=%s: skipping %s", START_AT_CHART, tpl$name, verbose=VERBOSE); next } }
      if (!isTRUE(FORCE_RECALIBRATE) && state_has_threshold(state, n0, target, tpl$name)) { log_msg("Threshold already present: n0=%d target=%g chart=%s", n0, target, tpl$name, verbose=VERBOSE); next }
      state <- remove_threshold(state, n0, target, tpl$name)
      cal <- calibrate_from_path(paths[[tpl_i]], target, tpl$direction, tpl$lower, tpl$upper, tpl$log_scale, N_ITER_CAL, tpl$name)
      chart <- tpl$chart; chart$threshold <- cal$threshold; chart$direction <- tpl$direction
      thr <- data.frame(n0=n0, target_arl=target, chart=tpl$name, chart_label=make_chart_label(chart), chart_type=chart$type, threshold=cal$threshold, direction=tpl$direction, lambda=chart$lambda %||% NA_real_, r=chart$r %||% NA_real_, kappa=chart$kappa %||% NA_real_, ks_engine=KS_ENGINE, selected_mean_R1=cal$selected$mean_R1, selected_se_R1=cal$selected$se_R1, selected_censored_rate=cal$selected$censored_rate, stringsAsFactors=FALSE)
      hist <- cal$history; hist$n0 <- n0; hist$target_arl <- target; hist$chart <- tpl$name; hist <- hist[, names(empty_cal_history), drop=FALSE]
      state$thresholds <- rbind(state$thresholds, thr[, names(empty_thresholds), drop=FALSE]); state$calibration_history <- rbind(state$calibration_history, hist)
      save_state(state); write_public_outputs(state)
    }
    rm(paths); invisible(gc())
  }
}

if ("ic" %in% RUN_STAGES) {
  for (n0 in N0_SET) for (target in TARGET_ARL) {
    if (!isTRUE(FORCE_FINAL_IC) && state_has_ic(state, n0, target)) { log_msg("Final IC already present: n0=%d target=%g; skipping.", n0, target, verbose=VERBOSE); next }
    charts <- get_charts_from_state(state, n0, target); if (!length(charts)) { warning(sprintf("No thresholds for n0=%d target=%g; skipping IC", n0, target)); next }
    state <- remove_ic(state, n0, target); log_msg("Final IC verification: n0=%d target=%g charts=%d", n0, target, length(charts), verbose=VERBOSE)
    raw_df <- run_replicates_final(N_REP_FINAL_IC, 50000L + 1000L*as.integer(n0) + as.integer(target), charts, n0, "IC", c(1L,5L), MAX_T_FINAL)
    if (isTRUE(SAVE_RAW_RDS)) safe_save_rds(raw_df, file.path(OUTPUT_DIR,"raw_rds",sprintf("final_ic_n0%d_target%d_raw.rds", n0, as.integer(target))))
    sm <- summarise_run_df(raw_df); sm$n0 <- n0; sm$target_arl <- target; sm$ks_engine <- KS_ENGINE; sm <- sm[, names(empty_ic), drop=FALSE]
    state$ic <- rbind(state$ic, sm); save_state(state); write_public_outputs(state); rm(raw_df); invisible(gc())
  }
}

if ("ooc" %in% RUN_STAGES) {
  for (n0 in N0_SET) for (target in TARGET_ARL) {
    charts <- get_charts_from_state(state, n0, target); if (!length(charts)) { warning(sprintf("No thresholds for n0=%d target=%g; skipping OOC", n0, target)); next }
    for (ooc in OOC_NAMES) {
      if (!isTRUE(FORCE_OOC) && state_has_ooc(state, n0, target, ooc)) { log_msg("OOC already present: n0=%d target=%g ooc=%s; skipping.", n0, target, ooc, verbose=VERBOSE); next }
      state <- remove_ooc(state, n0, target, ooc); log_msg("OOC evaluation: n0=%d target=%g OOC=%s charts=%d", n0, target, phase2_label(ooc), length(charts), verbose=VERBOSE)
      raw_df <- run_replicates_final(N_REP_OOC, 70000L + 100000L*as.integer(match(ooc, OOC_NAMES)) + 1000L*as.integer(n0) + as.integer(target), charts, n0, ooc, c(1L,5L), MAX_T_OOC)
      if (isTRUE(SAVE_RAW_RDS)) safe_save_rds(raw_df, file.path(OUTPUT_DIR,"raw_rds",sprintf("ooc_n0%d_target%d_%s_raw.rds", n0, as.integer(target), sanitize_key(ooc))))
      sm <- summarise_run_df(raw_df); sm$n0 <- n0; sm$target_arl <- target; sm$ooc <- ooc; sm$ooc_label <- phase2_label(ooc); sm$ks_engine <- KS_ENGINE; sm <- sm[, names(empty_ooc), drop=FALSE]
      state$ooc <- rbind(state$ooc, sm); save_state(state); write_public_outputs(state); rm(raw_df); invisible(gc())
    }
  }
}

write_public_outputs(state); save_state(state)
log_msg("Finished matched-ARL benchmark simulations.", verbose=VERBOSE)
log_msg("Outputs: %s", normalizePath(OUTPUT_DIR, mustWork=FALSE), verbose=VERBOSE)
