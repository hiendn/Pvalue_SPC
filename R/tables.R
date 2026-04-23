# Table summaries for the manuscript numerics.
#
# The LaTeX tables are typeset in the paper; this script regenerates the compact
# CSV summaries used to check the benchmark calibration and delay comparisons.

set_default <- function(name, value, envir = .GlobalEnv) {
  if (!exists(name, envir = envir, inherits = FALSE)) assign(name, value, envir = envir)
  invisible(get(name, envir = envir, inherits = FALSE))
}
set_default("DATA_DIR", "data")
set_default("OUTPUT_DIR", "output")
set_default("TABLE_DIR", "tables")

ensure_dir <- function(path) { if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE); invisible(path) }
ensure_dir(TABLE_DIR)

find_csv <- function(filename) {
  candidates <- c(file.path(OUTPUT_DIR, filename), file.path(DATA_DIR, filename), filename)
  hit <- candidates[file.exists(candidates)]
  if (!length(hit)) stop("Could not find required CSV: ", filename)
  hit[[1L]]
}
read_csv <- function(filename) utils::read.csv(find_csv(filename), stringsAsFactors = FALSE, check.names = FALSE)
write_out <- function(df, filename) {
  utils::write.csv(df, file.path(TABLE_DIR, filename), row.names = FALSE)
  message("Wrote ", file.path(TABLE_DIR, filename))
  invisible(df)
}
range_string <- function(x, digits = 2) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_character_)
  sprintf(paste0("%.", digits, "f--%.", digits, "f"), min(x), max(x))
}
chart_math_label <- function(x) {
  map <- c(
    "p-value Shewhart" = "Raw P_t",
    "Qbar lambda=0.95" = "Qbar(lambda=0.95,r=1)",
    "Qtilde r=-0.8" = "Qtilde(lambda=0.5,r=-0.8)",
    "Qtilde r=-0.9" = "Qtilde(lambda=0.5,r=-0.9)",
    "KS Shewhart statistic" = "KS Shewhart statistic",
    "KS-EWMA statistic" = "KS-EWMA statistic",
    "KS-CUSUM statistic" = "KS-CUSUM statistic"
  )
  unname(ifelse(x %in% names(map), map[x], x))
}

make_benchmark_ic_summary <- function() {
  ic <- read_csv("ks_benchmark_ic_fast.csv")
  ic$target_for_k <- ifelse(ic$k == 1, ic$target_arl, 5 * ic$target_arl)
  ic$ratio_to_target <- ic$mean_Rk / ic$target_for_k
  charts <- c("p-value Shewhart", "Qbar lambda=0.95", "Qtilde r=-0.8", "Qtilde r=-0.9",
              "KS Shewhart statistic", "KS-EWMA statistic", "KS-CUSUM statistic")
  rows <- lapply(charts, function(ch) {
    z <- ic[ic$chart == ch, , drop = FALSE]
    data.frame(
      chart = chart_math_label(ch),
      R1_ratio_range = range_string(z$ratio_to_target[z$k == 1], 2),
      R5_ratio_range = range_string(z$ratio_to_target[z$k == 5], 2),
      max_relative_MCSE = max(z$mcse_ratio, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  write_out(out, "table_ks_benchmark_ic_summary.csv")
}

make_benchmark_delay_summary <- function() {
  ooc <- read_csv("ks_benchmark_ooc_fast.csv")
  raw <- ooc[ooc$chart == "p-value Shewhart", c("n0", "target_arl", "ooc", "k", "mean_Rk"), drop = FALSE]
  names(raw)[names(raw) == "mean_Rk"] <- "raw_mean_Rk"
  m <- merge(ooc, raw, by = c("n0", "target_arl", "ooc", "k"), all.x = TRUE)
  m$delay_ratio <- m$mean_Rk / m$raw_mean_Rk
  m$scenario <- ifelse(grepl("^dyn", m$ooc), "dynamic", "persistent")
  charts <- c("KS Shewhart statistic", "KS-EWMA statistic", "KS-CUSUM statistic",
              "Qbar lambda=0.95", "Qtilde r=-0.8", "Qtilde r=-0.9")
  rows <- lapply(charts, function(ch) {
    z <- m[m$chart == ch, , drop = FALSE]
    data.frame(
      chart = chart_math_label(ch),
      persistent_R1 = stats::median(z$delay_ratio[z$scenario == "persistent" & z$k == 1], na.rm = TRUE),
      persistent_R5 = stats::median(z$delay_ratio[z$scenario == "persistent" & z$k == 5], na.rm = TRUE),
      dynamic_R1 = stats::median(z$delay_ratio[z$scenario == "dynamic" & z$k == 1], na.rm = TRUE),
      dynamic_R5 = stats::median(z$delay_ratio[z$scenario == "dynamic" & z$k == 5], na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  write_out(out, "table_ks_benchmark_relative_delay_summary.csv")
}

make_benchmark_threshold_tables <- function() {
  thr <- read_csv("ks_benchmark_thresholds_fast.csv")
  ic <- read_csv("ks_benchmark_ic_fast.csv")
  ic_wide <- reshape(ic[, c("n0", "target_arl", "chart", "k", "mean_Rk", "se_Rk", "mcse_ratio", "censored_rate")],
                     idvar = c("n0", "target_arl", "chart"), timevar = "k", direction = "wide")
  merged <- merge(thr, ic_wide, by = c("n0", "target_arl", "chart"), all.x = TRUE)
  for (n0 in sort(unique(merged$n0))) {
    write_out(merged[merged$n0 == n0, , drop = FALSE], sprintf("table_ks_benchmark_thresholds_n0_%s.csv", n0))
  }
}

make_ewma_multiplier_table <- function() {
  mult <- read_csv("ewma_multiplier_table.csv")
  write_out(mult, "table_ewma_effective_thresholds.csv")
}

make_all_table_summaries <- function() {
  make_benchmark_ic_summary()
  make_benchmark_delay_summary()
  make_benchmark_threshold_tables()
  make_ewma_multiplier_table()
}

if (sys.nframe() == 0L) make_all_table_summaries()
