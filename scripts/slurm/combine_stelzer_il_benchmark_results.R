#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) == 0L) return(default)
  sub(paste0("^--", name, "="), "", hit[[1]])
}

out_dir <- path.expand(get_arg("out_dir", "~/Mimesys_validation/stelzer_integratedlearner_benchmark"))
task_dir <- file.path(out_dir, "tasks")
raw_out <- file.path(out_dir, "stelzer_integratedlearner_benchmark_raw.csv")
summary_out <- file.path(out_dir, "stelzer_integratedlearner_benchmark_summary.csv")
failures_out <- file.path(out_dir, "stelzer_integratedlearner_benchmark_failures.csv")

files <- sort(Sys.glob(file.path(task_dir, "task_*.csv")))
if (length(files) == 0L) stop("No task CSV files found in ", task_dir, call. = FALSE)

read_one <- function(path) {
  x <- read.csv(path, stringsAsFactors = FALSE)
  x$task_file <- basename(path)
  x
}

raw <- do.call(rbind, lapply(files, read_one))
write.csv(raw, raw_out, row.names = FALSE)

ok <- raw[raw$learner != "fit_failed", , drop = FALSE]
groups <- c(
  "dataset", "scenario", "feature_setting", "feature_counts", "feature_n", "nsample", "crosstalk", "outcome_type",
  "effect_setting", "dependence", "snr", "de_prob", "de_fac_loc",
  "de_fac_scale", "signal_share", "learner"
)
groups <- groups[groups %in% names(ok)]

if (nrow(ok) == 0L) {
  summary <- ok[, groups, drop = FALSE]
  summary$performance_mean <- numeric()
  summary$performance_sd <- numeric()
  summary$gain_over_best_single_mean <- numeric()
  summary$gain_over_best_single_sd <- numeric()
  summary$win_rate <- numeric()
  summary$n_success <- integer()
  summary$n_total <- integer()
  write.csv(summary, summary_out, row.names = FALSE)

  failures <- raw[!is.na(raw$error) & nzchar(raw$error), , drop = FALSE]
  write.csv(failures, failures_out, row.names = FALSE)

  cat("Combined task files:", length(files), "\n")
  cat("Raw:", raw_out, "\n")
  cat("Summary:", summary_out, "\n")
  cat("Failures:", failures_out, "\n")
  quit(save = "no", status = 0)
}

suppressPackageStartupMessages({
  library(dplyr)
})

summary <- ok |>
  dplyr::group_by(dplyr::across(dplyr::all_of(groups))) |>
  dplyr::summarise(
    performance_mean = mean(performance, na.rm = TRUE),
    performance_sd = stats::sd(performance, na.rm = TRUE),
    gain_over_best_single_mean = mean(gain_over_best_single, na.rm = TRUE),
    gain_over_best_single_sd = stats::sd(gain_over_best_single, na.rm = TRUE),
    win_rate = mean(gain_over_best_single > 0, na.rm = TRUE),
    n_success = sum(!is.na(performance)),
    n_total = dplyr::n(),
    .groups = "drop"
  ) |>
  dplyr::arrange(feature_setting, feature_n, nsample, effect_setting, crosstalk, learner)

write.csv(summary, summary_out, row.names = FALSE)

failures <- raw[!is.na(raw$error) & nzchar(raw$error), , drop = FALSE]
write.csv(failures, failures_out, row.names = FALSE)

cat("Combined task files:", length(files), "\n")
cat("Raw:", raw_out, "\n")
cat("Summary:", summary_out, "\n")
cat("Failures:", failures_out, "\n")
