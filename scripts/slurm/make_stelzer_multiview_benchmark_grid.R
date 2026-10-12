#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) == 0L) return(default)
  sub(paste0("^--", name, "="), "", hit[[1]])
}

split_arg <- function(x, sep = ",") {
  if (is.null(x) || !nzchar(x)) return(character())
  out <- trimws(strsplit(x, sep, fixed = TRUE)[[1]])
  out[nzchar(out)]
}

base_dir <- path.expand(get_arg("base_dir", "~/Mimesys_validation"))
out_dir <- path.expand(get_arg("out_dir", file.path(base_dir, "stelzer_multiview_cooperative_benchmark")))
feature_sizes <- as.integer(split_arg(get_arg("feature_sizes", "100,200,500")))
sample_sizes <- as.integer(split_arg(get_arg("sample_sizes", "53,100,200,500")))
crosstalk_grid <- as.numeric(split_arg(get_arg("crosstalk", "0,0.1,0.25,0.5,0.75,0.9,1")))
rho_grid <- as.numeric(split_arg(get_arg("rho_grid", "0.1,0.25,0.5,0.75")))
scenarios <- split_arg(get_arg("scenarios", "balanced_aligned,balanced_random"))
ygen_modes <- split_arg(get_arg("ygen_modes", "LM"))
friedman_share <- get_arg("friedman_share", "")
nrep <- as.integer(get_arg("nrep", "20"))
folds <- as.integer(get_arg("folds", "5"))
seed_base <- as.integer(get_arg("seed_base", "30000"))
snr <- as.numeric(get_arg("snr", "0.5"))
de_prob <- as.numeric(get_arg("de_prob", "0.02"))
de_fac_loc <- as.numeric(get_arg("de_fac_loc", "0.8"))
de_fac_scale <- as.numeric(get_arg("de_fac_scale", "0.2"))
dependence <- get_arg("dependence", "factor")

if (any(is.na(feature_sizes)) || any(feature_sizes <= 0)) stop("--feature_sizes must be positive integers.", call. = FALSE)
if (any(is.na(sample_sizes)) || any(sample_sizes <= 0)) stop("--sample_sizes must be positive integers.", call. = FALSE)
if (any(is.na(crosstalk_grid))) stop("--crosstalk must be numeric.", call. = FALSE)
if (any(is.na(rho_grid))) stop("--rho_grid must be numeric.", call. = FALSE)
if (is.na(nrep) || nrep <= 0) stop("--nrep must be positive.", call. = FALSE)
if (folds < 3) stop("--folds must be at least 3.", call. = FALSE)

allowed_scenarios <- c("balanced_aligned", "balanced_random")
unknown <- setdiff(scenarios, allowed_scenarios)
if (length(unknown) > 0L) stop("Unknown scenario(s): ", paste(unknown, collapse = ", "), call. = FALSE)
allowed_ygen_modes <- c("LM", "Friedman", "Friedman2", "SharedLatent")
unknown_ygen <- setdiff(ygen_modes, allowed_ygen_modes)
if (length(unknown_ygen) > 0L) stop("Unknown ygen mode(s): ", paste(unknown_ygen, collapse = ", "), call. = FALSE)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

grid <- expand.grid(
  feature_n = feature_sizes,
  nsample = sample_sizes,
  crosstalk = crosstalk_grid,
  scenario = scenarios,
  ygen_mode = ygen_modes,
  rep = seq_len(nrep),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

grid$job_id <- seq_len(nrow(grid))
grid$dataset <- "StelzerDOS"
grid$seed <- seed_base + grid$job_id
grid$rho_grid <- paste(rho_grid, collapse = ",")
grid$folds <- folds
grid$snr <- snr
grid$de_prob <- de_prob
grid$de_fac_loc <- de_fac_loc
grid$de_fac_scale <- de_fac_scale
grid$dependence <- dependence
grid$friedman_share <- friedman_share

grid <- grid[, c(
  "job_id", "dataset", "scenario", "ygen_mode", "feature_n", "nsample", "crosstalk",
  "rep", "seed", "rho_grid", "folds", "snr", "de_prob", "de_fac_loc",
  "de_fac_scale", "dependence", "friedman_share"
)]

grid_file <- file.path(out_dir, "stelzer_multiview_benchmark_grid.csv")
write.csv(grid, grid_file, row.names = FALSE)

cat("Wrote grid:\n")
cat(grid_file, "\n")
cat("Rows:", nrow(grid), "\n")
cat("Feature sizes:", paste(feature_sizes, collapse = ", "), "\n")
cat("Sample sizes:", paste(sample_sizes, collapse = ", "), "\n")
cat("Crosstalk:", paste(crosstalk_grid, collapse = ", "), "\n")
cat("Intermediate rho grid:", paste(rho_grid, collapse = ", "), "\n")
cat("Scenarios:", paste(scenarios, collapse = ", "), "\n")
cat("Outcome generators:", paste(ygen_modes, collapse = ", "), "\n")
if (nzchar(friedman_share)) cat("Friedman share:", friedman_share, "\n")
