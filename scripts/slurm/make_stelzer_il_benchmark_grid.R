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

parse_key_values <- function(x, label, value_fun = identity) {
  if (is.null(x) || !nzchar(x)) return(character())
  pieces <- split_arg(x, ",")
  out <- vector("list", length(pieces))
  names(out) <- character(length(pieces))
  for (i in seq_along(pieces)) {
    kv <- strsplit(pieces[[i]], "=", fixed = TRUE)[[1]]
    if (length(kv) != 2L || !nzchar(kv[[1]]) || !nzchar(kv[[2]])) {
      stop(label, " entries must look like name=value.", call. = FALSE)
    }
    names(out)[[i]] <- trimws(kv[[1]])
    out[[i]] <- value_fun(trimws(kv[[2]]))
  }
  unlist(out, use.names = TRUE)
}

parse_feature_sets <- function(x) {
  pieces <- split_arg(x, ";")
  if (length(pieces) == 0L) stop("--feature_sets cannot be empty.", call. = FALSE)
  out <- vector("list", length(pieces))
  names(out) <- character(length(pieces))
  for (i in seq_along(pieces)) {
    parts <- strsplit(pieces[[i]], ":", fixed = TRUE)[[1]]
    if (length(parts) != 2L || !nzchar(parts[[1]]) || !nzchar(parts[[2]])) {
      stop("--feature_sets entries must look like setting:view=count,view=count.", call. = FALSE)
    }
    names(out)[[i]] <- trimws(parts[[1]])
    out[[i]] <- parse_key_values(trimws(parts[[2]]), paste0("feature set ", names(out)[[i]]), as.integer)
  }
  out
}

format_feature_counts <- function(x) {
  paste(paste(names(x), as.integer(x), sep = "="), collapse = ",")
}

base_dir <- path.expand(get_arg("base_dir", "~/Mimesys_validation"))
out_dir <- path.expand(get_arg("out_dir", file.path(base_dir, "stelzer_integratedlearner_benchmark")))

feature_grid <- as.integer(split_arg(get_arg("feature_grid", "100,200")))
feature_sets_arg <- get_arg("feature_sets", NULL)
sample_sizes <- as.integer(split_arg(get_arg("sample_sizes", "100,200,500")))
crosstalk_grid <- as.numeric(split_arg(get_arg("crosstalk", "0,0.5,0.75,1")))
outcome_types <- split_arg(get_arg("outcomes", "continuous"))
effect_settings <- split_arg(get_arg("effects", "random,aligned"))
nrep <- as.integer(get_arg("nrep", "20"))

dependence <- get_arg("dependence", "factor")
base_learner <- get_arg("base_learner", "SL.randomForest")
meta_learner <- get_arg("meta_learner", "sl_lasso")
folds <- as.integer(get_arg("folds", "2"))
snr <- as.numeric(get_arg("snr", "0.5"))
de_prob <- as.numeric(get_arg("de_prob", "0.02"))
de_fac_loc <- as.numeric(get_arg("de_fac_loc", "0.8"))
de_fac_scale <- as.numeric(get_arg("de_fac_scale", "0.2"))
signal_share <- get_arg("signal_share", "")
seed_base <- as.integer(get_arg("seed_base", "1234"))

if (!is.null(feature_sets_arg)) {
  feature_sets <- parse_feature_sets(feature_sets_arg)
  if (any(vapply(feature_sets, function(x) any(is.na(x)) || any(x <= 0), logical(1)))) {
    stop("--feature_sets counts must be positive integers.", call. = FALSE)
  }
} else {
  if (any(is.na(feature_grid)) || any(feature_grid <= 0)) stop("--feature_grid must be positive integers.", call. = FALSE)
  feature_sets <- stats::setNames(
    lapply(feature_grid, function(n) c(Immunome = n, Metabolome = n, Proteome = n)),
    paste0("features_", feature_grid)
  )
}
if (any(is.na(sample_sizes)) || any(sample_sizes <= 0)) stop("--sample_sizes must be positive integers.", call. = FALSE)
if (any(is.na(crosstalk_grid))) stop("--crosstalk must be numeric.", call. = FALSE)
if (is.na(nrep) || nrep <= 0) stop("--nrep must be a positive integer.", call. = FALSE)
if (is.na(folds) || folds <= 1) stop("--folds must be an integer > 1.", call. = FALSE)
if (any(!outcome_types %in% c("continuous", "binary"))) stop("--outcomes must contain continuous and/or binary.", call. = FALSE)
if (any(!effect_settings %in% c("random", "aligned"))) stop("--effects must contain random and/or aligned.", call. = FALSE)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

grid <- expand.grid(
  dataset = "StelzerDOS",
  feature_setting = names(feature_sets),
  nsample = sample_sizes,
  crosstalk = crosstalk_grid,
  outcome_type = outcome_types,
  effect_setting = effect_settings,
  rep = seq_len(nrep),
  stringsAsFactors = FALSE
)

grid <- grid[order(
  grid$feature_setting,
  grid$nsample,
  grid$effect_setting,
  grid$outcome_type,
  grid$crosstalk,
  grid$rep
), ]

grid$job_id <- seq_len(nrow(grid))
grid$seed <- seed_base + grid$job_id
grid$feature_counts <- vapply(feature_sets[grid$feature_setting], format_feature_counts, character(1))
grid$feature_n <- vapply(feature_sets[grid$feature_setting], function(x) {
  if (length(unique(as.integer(x))) == 1L) as.integer(x[[1]]) else NA_integer_
}, integer(1))
grid$dependence <- dependence
grid$base_learner <- base_learner
grid$meta_learner <- meta_learner
grid$folds <- folds
grid$snr <- snr
grid$de_prob <- de_prob
grid$de_fac_loc <- de_fac_loc
grid$de_fac_scale <- de_fac_scale
grid$signal_share <- signal_share

grid <- grid[
  ,
  c(
    "job_id", "dataset", "feature_setting", "feature_counts", "feature_n", "nsample", "crosstalk", "outcome_type",
    "effect_setting", "rep", "seed", "dependence", "base_learner",
    "meta_learner", "folds", "snr", "de_prob", "de_fac_loc", "de_fac_scale", "signal_share"
  )
]

grid_file <- file.path(out_dir, "stelzer_integratedlearner_benchmark_grid.csv")
write.csv(grid, grid_file, row.names = FALSE)

cat("Wrote grid:\n")
cat(grid_file, "\n")
cat("Rows:", nrow(grid), "\n")
cat("Feature sets:\n")
for (nm in names(feature_sets)) {
  cat(" ", nm, ":", format_feature_counts(feature_sets[[nm]]), "\n")
}
cat("Sample sizes:", paste(sample_sizes, collapse = ", "), "\n")
cat("Crosstalk:", paste(crosstalk_grid, collapse = ", "), "\n")
cat("Outcomes:", paste(outcome_types, collapse = ", "), "\n")
cat("Effects:", paste(effect_settings, collapse = ", "), "\n")
cat("Replicates per setting:", nrep, "\n")
