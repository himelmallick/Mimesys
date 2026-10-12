#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) == 0L) return(default)
  sub(paste0("^--", name, "="), "", hit[[1]])
}

keep_top_variable <- function(x, n_features) {
  vars <- apply(x, 2, stats::var, na.rm = TRUE)
  vars <- vars[is.finite(vars)]
  vars <- sort(vars, decreasing = TRUE)
  n_keep <- min(as.integer(n_features), length(vars))
  if (n_keep <= 0L) stop("No variable features available.", call. = FALSE)
  x[, names(vars)[seq_len(n_keep)], drop = FALSE]
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

load_stelzer_views <- function(data_path) {
  env <- new.env(parent = emptyenv())
  load(data_path, envir = env)
  if (!exists("pcl", envir = env, inherits = FALSE)) {
    stop("Expected object `pcl` in ", data_path, call. = FALSE)
  }
  pcl <- get("pcl", envir = env, inherits = FALSE)
  feature_table <- pcl$feature_table
  sample_metadata <- pcl$sample_metadata
  feature_metadata <- pcl$feature_metadata

  subject_id <- sample_metadata$subjectID
  subject_feature <- t(apply(feature_table, 1, function(x) {
    tapply(x, subject_id, mean, na.rm = TRUE)
  }))

  feature_type <- feature_metadata$featureType
  names(feature_type) <- rownames(feature_metadata)

  list(
    Immunome = t(subject_feature[feature_type == "CYTOF", , drop = FALSE]),
    Metabolome = t(subject_feature[feature_type == "Metabolomics", , drop = FALSE]),
    Proteome = t(subject_feature[feature_type == "Proteomics", , drop = FALSE])
  )
}

grid_file <- path.expand(get_arg("grid", "~/Mimesys_validation/stelzer_integratedlearner_benchmark/stelzer_integratedlearner_benchmark_grid.csv"))
out_dir <- path.expand(get_arg("out_dir", "~/Mimesys_validation/stelzer_integratedlearner_benchmark"))
data_path <- path.expand(get_arg("data_path", "~/Mimesys_2.0/inst/extdata/StelzerDOS.RData"))
task_id_arg <- get_arg("task_id", Sys.getenv("SLURM_ARRAY_TASK_ID", unset = NA_character_))
task_id <- as.integer(task_id_arg)
if (is.na(task_id)) stop("No task id supplied. Use --task_id=... or Slurm array.", call. = FALSE)

load_mimesys <- function() {
  project_dir <- Sys.getenv("MIMESYS_PROJECT_DIR", unset = getwd())
  if (requireNamespace("devtools", quietly = TRUE) && file.exists(file.path(project_dir, "DESCRIPTION"))) {
    devtools::load_all(project_dir, quiet = TRUE)
  } else {
    suppressPackageStartupMessages(library(MimeSys))
  }
}

load_mimesys()

grid <- read.csv(grid_file, stringsAsFactors = FALSE)
if (task_id < 1L || task_id > nrow(grid)) stop("task_id outside grid range: ", task_id, call. = FALSE)
row <- grid[task_id, , drop = FALSE]

task_dir <- file.path(out_dir, "tasks")
dir.create(task_dir, recursive = TRUE, showWarnings = FALSE)
out_csv <- file.path(task_dir, sprintf("task_%05d.csv", task_id))
out_rds <- file.path(task_dir, sprintf("task_%05d.rds", task_id))

if (file.exists(out_csv) && file.exists(out_rds)) {
  message("Task already complete; skipping: ", task_id)
  quit(save = "no", status = 0)
}

views_full <- load_stelzer_views(data_path)
feature_counts <- if ("feature_counts" %in% names(row) && nzchar(row$feature_counts[[1]])) {
  parse_key_values(row$feature_counts[[1]], "feature_counts", as.integer)
} else {
  stats::setNames(rep(row$feature_n[[1]], length(views_full)), names(views_full))
}
missing_counts <- setdiff(names(views_full), names(feature_counts))
if (length(missing_counts) > 0L) {
  stop("feature_counts missing view(s): ", paste(missing_counts, collapse = ", "), call. = FALSE)
}
views <- lapply(names(views_full), function(view) {
  keep_top_variable(views_full[[view]], feature_counts[[view]])
})
names(views) <- names(views_full)
template <- fit_template(views)

effect_setting <- row$effect_setting[[1]]
de_align <- if (identical(effect_setting, "aligned")) 1 else NULL
signal_share <- if ("signal_share" %in% names(row) && nzchar(row$signal_share[[1]])) {
  parse_key_values(row$signal_share[[1]], "signal_share", as.numeric)
} else {
  NULL
}

sim <- gen_simmba(
  nsample = row$nsample[[1]],
  template = template,
  crosstalk = row$crosstalk[[1]],
  outcome.type = row$outcome_type[[1]],
  de.align = de_align,
  dependence = row$dependence[[1]],
  de.prob = row$de_prob[[1]],
  de.facLoc = rep(row$de_fac_loc[[1]], length(template$views)),
  de.facScale = rep(row$de_fac_scale[[1]], length(template$views)),
  signal.share = signal_share,
  snr = row$snr[[1]],
  nrep = 1,
  seed = row$seed[[1]]
)

bench <- benchmark_integratedlearner_simulation(
  sim = sim,
  outcome_type = row$outcome_type[[1]],
  folds = row$folds[[1]],
  base_learner = row$base_learner[[1]],
  meta_learner = row$meta_learner[[1]],
  run_concat = TRUE,
  run_stacked = TRUE,
  extra_metadata = list(
    job_id = row$job_id[[1]],
    dataset = row$dataset[[1]],
    scenario = if ("scenario" %in% names(row)) row$scenario[[1]] else "",
    feature_setting = if ("feature_setting" %in% names(row)) row$feature_setting[[1]] else paste0("features_", row$feature_n[[1]]),
    feature_counts = paste(paste(names(feature_counts), feature_counts, sep = "="), collapse = ","),
    feature_n = row$feature_n[[1]],
    nsample = row$nsample[[1]],
    crosstalk = row$crosstalk[[1]],
    outcome_type = row$outcome_type[[1]],
    effect_setting = row$effect_setting[[1]],
    grid_rep = row$rep[[1]],
    seed = row$seed[[1]],
    dependence = row$dependence[[1]],
    base_learner = row$base_learner[[1]],
    meta_learner = row$meta_learner[[1]],
    folds = row$folds[[1]],
    snr = row$snr[[1]],
    de_prob = row$de_prob[[1]],
    de_fac_loc = row$de_fac_loc[[1]],
    de_fac_scale = row$de_fac_scale[[1]],
    signal_share = if (is.null(signal_share)) "" else paste(paste(names(signal_share), signal_share, sep = "="), collapse = ",")
  )
)

saveRDS(list(grid_row = row, benchmark = bench), out_rds)
write.csv(bench, out_csv, row.names = FALSE)
message("Wrote ", out_csv)
