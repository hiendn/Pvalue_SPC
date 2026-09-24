#!/usr/bin/env Rscript

# Derive robust stopped-run diagnostics from the retained per-path files.
# Run from R/auxiliary.

options(stringsAsFactors = FALSE)
input_dir <- file.path("output", "checkpoints")
files <- sort(list.files(input_dir,
                         pattern = "^localisation_.*_raw\\.csv$",
                         full.names = TRUE))
if (length(files) != 48L) {
  stop("Expected 48 localisation checkpoint files; found ", length(files))
}

decode <- function(x, prefix) {
  as.numeric(gsub("p", ".", sub(prefix, "", x, fixed = TRUE), fixed = TRUE))
}

rows <- lapply(files, function(path) {
  stem <- sub("_raw\\.csv$", "", basename(path))
  bits <- strsplit(stem, "_", fixed = TRUE)[[1L]]
  is_cauchy <- bits[2L] == "cauchy"
  offset <- if (is_cauchy) 1L else 0L
  n0 <- if (is_cauchy) as.integer(sub("^n0", "", bits[3L])) else NA_integer_
  delta <- decode(bits[3L + offset], "delta")
  rho <- decode(bits[4L + offset], "rho")
  alpha <- decode(bits[5L + offset], "alpha")
  z <- read.csv(path)
  se_r <- stats::sd(z$R) / sqrt(nrow(z))
  data.frame(
    model = if (is_cauchy) "cauchy_MW_exact_fast" else "normal_Z",
    n0 = n0, delta = delta, rho = rho, alpha = alpha,
    n_rep = nrow(z), mean_R = mean(z$R), se_R = se_r,
    relative_MCSE = se_r / mean(z$R), median_R = stats::median(z$R),
    q90_R = unname(stats::quantile(z$R, 0.90, type = 1)),
    PTS_25 = mean(z$R <= 25L), PTS_100 = mean(z$R <= 100L),
    censored_rate = mean(z$censored),
    stringsAsFactors = FALSE
  )
})

out <- do.call(rbind, rows)
write.csv(out, file.path("output", "localisation_stopped_diagnostics.csv"),
          row.names = FALSE)
cat("Wrote", nrow(out), "rows to output/localisation_stopped_diagnostics.csv\n")
