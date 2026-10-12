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

format_feature_counts <- function(x) {
  paste(paste(names(x), as.integer(x), sep = "="), collapse = ",")
}

make_balanced_feature_counts <- function(view_names, n) {
  stats::setNames(rep(as.integer(n), length(view_names)), view_names)
}

base_dir <- path.expand(get_arg("base_dir", "~/Mimesys_validation"))
out_dir <- path.expand(get_arg("out_dir", file.path(base_dir, "stelzer_il_outcome_scenarios")))

view_names <- split_arg(get_arg("views", "Immunome,Metabolome,Proteome"))
feature_sizes_arg <- get_arg("feature_sizes", NULL)
feature_counts_arg <- get_arg("feature_counts", NULL)
if (!is.null(feature_sizes_arg) && !is.null(feature_counts_arg)) {
  stop("Use either --feature_sizes for balanced layer sizes or --feature_counts for one named setting, not both.", call. = FALSE)
}
if (!is.null(feature_sizes_arg)) {
  feature_sizes <- as.integer(split_arg(feature_sizes_arg))
  feature_count_sets <- lapply(feature_sizes, function(n) make_balanced_feature_counts(view_names, n))
} else {
  feature_counts <- parse_key_values(
    if (is.null(feature_counts_arg)) "Immunome=200,Metabolome=200,Proteome=200" else feature_counts_arg,
    "--feature_counts",
    as.integer
  )
  feature_count_sets <- list(feature_counts)
}
sample_sizes <- as.integer(split_arg(get_arg("sample_sizes", "53,100,200,500")))
crosstalk_dense <- as.numeric(split_arg(get_arg("crosstalk", "0,0.1,0.25,0.5,0.75,0.9,1")))
nrep <- as.integer(get_arg("nrep", "20"))
outcome_types <- split_arg(get_arg("outcomes", "continuous"))

dependence <- get_arg("dependence", "factor")
base_learner <- get_arg("base_learner", "SL.xgboost")
meta_learner <- get_arg("meta_learner", "sl_lasso")
folds <- as.integer(get_arg("folds", "2"))
snr <- as.numeric(get_arg("snr", "0.5"))
de_prob <- as.numeric(get_arg("de_prob", "0.02"))
de_fac_loc <- as.numeric(get_arg("de_fac_loc", "0.8"))
de_fac_scale <- as.numeric(get_arg("de_fac_scale", "0.2"))
seed_base <- as.integer(get_arg("seed_base", "10000"))

if (length(view_names) == 0L) {
  stop("--views must contain at least one view name.", call. = FALSE)
}
if (length(feature_count_sets) == 0L) {
  stop("No feature-count settings were requested.", call. = FALSE)
}
for (feature_counts in feature_count_sets) {
  if (length(feature_counts) == 0L || any(is.na(feature_counts)) || any(feature_counts <= 0)) {
    stop("Feature counts must be positive integers.", call. = FALSE)
  }
  if (length(setdiff(view_names, names(feature_counts))) > 0L) {
    stop("Feature counts are missing view(s): ", paste(setdiff(view_names, names(feature_counts)), collapse = ", "), call. = FALSE)
  }
}
if (any(is.na(sample_sizes)) || any(sample_sizes <= 0)) stop("--sample_sizes must be positive integers.", call. = FALSE)
if (any(is.na(crosstalk_dense))) stop("Crosstalk values must be numeric.", call. = FALSE)
if (is.na(nrep) || nrep <= 0) stop("--nrep must be positive.", call. = FALSE)
if (any(!outcome_types %in% c("continuous", "binary"))) stop("--outcomes must contain continuous and/or binary.", call. = FALSE)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

scenario_defs <- list(
  single_immunome = list(
    signal_share = c(Immunome = 1, Metabolome = 0, Proteome = 0),
    effect_setting = "random",
    crosstalk = crosstalk_dense
  ),
  single_metabolome = list(
    signal_share = c(Immunome = 0, Metabolome = 1, Proteome = 0),
    effect_setting = "random",
    crosstalk = crosstalk_dense
  ),
  single_proteome = list(
    signal_share = c(Immunome = 0, Metabolome = 0, Proteome = 1),
    effect_setting = "random",
    crosstalk = crosstalk_dense
  ),
  balanced_aligned = list(
    signal_share = c(Immunome = 1, Metabolome = 1, Proteome = 1),
    effect_setting = "aligned",
    crosstalk = crosstalk_dense
  ),
  balanced_random = list(
    signal_share = c(Immunome = 1, Metabolome = 1, Proteome = 1),
    effect_setting = "random",
    crosstalk = crosstalk_dense
  )
)

format_signal_share <- function(x) {
  x <- x[view_names]
  paste(paste(names(x), x, sep = "="), collapse = ",")
}

rows <- list()
counter <- 1L
for (feature_counts in feature_count_sets) {
  feature_counts <- feature_counts[view_names]
  feature_setting <- paste0("fixed_", paste(as.integer(feature_counts), collapse = "_"))
  feature_n <- if (length(unique(as.integer(feature_counts))) == 1L) as.integer(feature_counts[[1]]) else NA_integer_
  feature_counts_string <- format_feature_counts(feature_counts)

  for (scenario in names(scenario_defs)) {
    def <- scenario_defs[[scenario]]
    for (nsample in sample_sizes) {
      for (crosstalk in def$crosstalk) {
        for (outcome_type in outcome_types) {
          for (rep in seq_len(nrep)) {
            rows[[counter]] <- data.frame(
              job_id = counter,
              dataset = "StelzerDOS",
              scenario = scenario,
              feature_setting = feature_setting,
              feature_counts = feature_counts_string,
              feature_n = feature_n,
              nsample = nsample,
              crosstalk = crosstalk,
              outcome_type = outcome_type,
              effect_setting = def$effect_setting,
              rep = rep,
              seed = seed_base + counter,
              dependence = dependence,
              base_learner = base_learner,
              meta_learner = meta_learner,
              folds = folds,
              snr = snr,
              de_prob = de_prob,
              de_fac_loc = de_fac_loc,
              de_fac_scale = de_fac_scale,
              signal_share = format_signal_share(def$signal_share),
              stringsAsFactors = FALSE
            )
            counter <- counter + 1L
          }
        }
      }
    }
  }
}

grid <- do.call(rbind, rows)
grid_file <- file.path(out_dir, "stelzer_integratedlearner_benchmark_grid.csv")
write.csv(grid, grid_file, row.names = FALSE)

cat("Wrote grid:\n")
cat(grid_file, "\n")
cat("Rows:", nrow(grid), "\n")
cat("Feature settings:", paste(unique(grid$feature_setting), collapse = ", "), "\n")
cat("Sample sizes:", paste(sample_sizes, collapse = ", "), "\n")
cat("Crosstalk:", paste(crosstalk_dense, collapse = ", "), "\n")
cat("Scenarios:", paste(names(scenario_defs), collapse = ", "), "\n")
cat("Base learner:", base_learner, "\n")
cat("Meta learner:", meta_learner, "\n")
