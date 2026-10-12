#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) == 0L) return(default)
  sub(paste0("^--", name, "="), "", hit[[1]])
}

out_dir <- path.expand(get_arg("out_dir", "~/Mimesys_validation/stelzer_multiview_cooperative_benchmark"))
task_dir <- file.path(out_dir, "tasks")
raw_out <- file.path(out_dir, "stelzer_multiview_cooperative_benchmark_raw.csv")
summary_out <- file.path(out_dir, "stelzer_multiview_cooperative_benchmark_summary.csv")
failures_out <- file.path(out_dir, "stelzer_multiview_cooperative_benchmark_failures.csv")

files <- sort(Sys.glob(file.path(task_dir, "task_*.csv")))
if (length(files) == 0L) stop("No task CSV files found in ", task_dir, call. = FALSE)

read_one <- function(path) {
  x <- read.csv(path, stringsAsFactors = FALSE)
  x$task_file <- basename(path)
  x
}

raw <- do.call(rbind, lapply(files, read_one))
write.csv(raw, raw_out, row.names = FALSE)

suppressPackageStartupMessages(library(dplyr))

groups <- c(
  "dataset", "scenario", "ygen_mode", "feature_n", "nsample", "crosstalk",
  "rho_grid", "folds", "snr", "de_prob", "de_fac_loc", "de_fac_scale",
  "dependence", "friedman_share", "method"
)
groups <- groups[groups %in% names(raw)]

summary <- raw |>
  dplyr::group_by(dplyr::across(dplyr::all_of(groups))) |>
  dplyr::summarise(
    performance_mean = mean(test_r2, na.rm = TRUE),
    performance_sd = stats::sd(test_r2, na.rm = TRUE),
    selected_rho_mean = mean(selected_rho, na.rm = TRUE),
    selected_rho_sd = stats::sd(selected_rho, na.rm = TRUE),
    cv_mse_mean = mean(cv_mse, na.rm = TRUE),
    cv_mse_sd = stats::sd(cv_mse, na.rm = TRUE),
    n_success = sum(is.finite(test_r2)),
    n_total = dplyr::n(),
    .groups = "drop"
  ) |>
  dplyr::arrange(ygen_mode, feature_n, nsample, scenario, crosstalk, method)

write.csv(summary, summary_out, row.names = FALSE)

failures <- raw[!is.na(raw$error) & nzchar(raw$error), , drop = FALSE]
write.csv(failures, failures_out, row.names = FALSE)

cat("Combined task files:", length(files), "\n")
cat("Raw:", raw_out, "\n")
cat("Summary:", summary_out, "\n")
cat("Failures:", failures_out, "\n")
