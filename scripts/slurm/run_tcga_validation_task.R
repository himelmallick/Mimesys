#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) == 0) {
    return(default)
  }
  sub(paste0("^--", name, "="), "", hit[[1]])
}

split_arg <- function(x, sep = ",") {
  if (is.null(x) || !nzchar(x)) {
    return(character())
  }
  out <- trimws(strsplit(x, sep, fixed = TRUE)[[1]])
  out[nzchar(out)]
}

parse_counts <- function(x, label) {
  pieces <- split_arg(x, ",")
  out <- vector("list", length(pieces))
  names(out) <- character(length(pieces))
  for (i in seq_along(pieces)) {
    kv <- strsplit(pieces[[i]], "=", fixed = TRUE)[[1]]
    if (length(kv) != 2L || !nzchar(kv[[1]]) || !nzchar(kv[[2]])) {
      stop(label, " entries must look like view=count.", call. = FALSE)
    }
    names(out)[[i]] <- trimws(kv[[1]])
    out[[i]] <- as.integer(trimws(kv[[2]]))
  }
  values <- unlist(out, use.names = TRUE)
  if (any(is.na(values)) || any(values <= 0)) {
    stop(label, " counts must be positive integers.", call. = FALSE)
  }
  values
}

parse_view_transforms <- function(x) {
  pieces <- split_arg(x, ",")
  if (length(pieces) == 0L) {
    return(character())
  }
  out <- vector("list", length(pieces))
  names(out) <- character(length(pieces))
  for (i in seq_along(pieces)) {
    kv <- strsplit(pieces[[i]], "=", fixed = TRUE)[[1]]
    if (length(kv) != 2L || !nzchar(kv[[1]]) || !nzchar(kv[[2]])) {
      stop("--view_transforms entries must look like view=transform.", call. = FALSE)
    }
    names(out)[[i]] <- trimws(kv[[1]])
    out[[i]] <- tolower(trimws(kv[[2]]))
  }
  values <- unlist(out, use.names = TRUE)
  allowed <- c("none", "identity", "log", "log1p", "clr")
  bad <- setdiff(unique(values), allowed)
  if (length(bad) > 0L) {
    stop("Unsupported transform(s): ", paste(bad, collapse = ", "), ". Allowed: ", paste(allowed, collapse = ", "))
  }
  values
}

is_present <- function(x) {
  length(x) == 1L && !is.na(x) && nzchar(trimws(as.character(x))) &&
    !identical(toupper(trimws(as.character(x))), "NA")
}

extract_template_input <- function(x) {
  if (is.list(x) && !is.null(x$template)) {
    return(x$template)
  }
  x
}

required_block_metrics <- c("mean_cor", "median_cor", "mean_abs_cor", "median_abs_cor")

task_output_current <- function(path) {
  if (!file.exists(path)) {
    return(FALSE)
  }
  x <- tryCatch(readRDS(path), error = function(e) NULL)
  if (is.null(x) || !is.list(x) || is.null(x$block_comparison)) {
    return(FALSE)
  }
  bc <- x$block_comparison
  bs <- if (!is.null(x$block_structure_agreement)) x$block_structure_agreement else x$correlation_structure_diagnostics
  bp <- x$block_structure_pairs
  is.data.frame(bc) && is.data.frame(bs) && is.data.frame(bp) &&
    "value" %in% names(bc) &&
    all(required_block_metrics %in% unique(as.character(bc$value))) &&
    all(c("cor_pearson", "cor_spearman", "mae_cor", "rmse_cor") %in% names(bs)) &&
    all(c("observed_cor", "simulated_cor") %in% names(bp))
}

keep_top_variable <- function(x, n_features) {
  if (is.null(n_features) || is.na(n_features)) {
    return(x)
  }
  vars <- apply(x, 2, stats::var, na.rm = TRUE)
  vars <- vars[is.finite(vars)]
  vars <- sort(vars, decreasing = TRUE)
  n_keep <- min(as.integer(n_features), length(vars))
  if (n_keep == 0L) {
    stop("No variable features available.")
  }
  x[, names(vars)[seq_len(n_keep)], drop = FALSE]
}

clr_transform <- function(x, pseudocount = NULL) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (any(x < 0, na.rm = TRUE)) {
    stop("CLR transform requires non-negative values.", call. = FALSE)
  }
  positive <- x[is.finite(x) & x > 0]
  if (is.null(pseudocount)) {
    pseudocount <- if (length(positive) > 0L) min(positive, na.rm = TRUE) / 2 else 1e-6
  }
  if (!is.finite(pseudocount) || pseudocount <= 0) {
    pseudocount <- 1e-6
  }
  x[!is.finite(x)] <- NA_real_
  x <- x + pseudocount
  log_x <- log(x)
  row_center <- rowMeans(log_x, na.rm = TRUE)
  out <- sweep(log_x, 1, row_center, "-")
  out[!is.finite(out)] <- NA_real_
  rownames(out) <- rownames(x)
  colnames(out) <- colnames(x)
  out
}

log_transform <- function(x, pseudocount = NULL, log1p = FALSE) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  min_x <- suppressWarnings(min(x, na.rm = TRUE))
  if (!is.finite(min_x)) {
    stop("Cannot log-transform a matrix with no finite values.", call. = FALSE)
  }
  if (isTRUE(log1p)) {
    if (min_x < -1) {
      shift <- abs(min_x) + 1
      x <- x + shift
    }
    out <- log1p(x)
  } else {
    if (is.null(pseudocount)) {
      positive <- x[is.finite(x) & x > 0]
      pseudocount <- if (length(positive) > 0L) min(positive, na.rm = TRUE) / 2 else 1e-6
    }
    if (min_x <= 0) {
      x <- x + abs(min_x) + pseudocount
    }
    out <- log(x)
  }
  out[!is.finite(out)] <- NA_real_
  rownames(out) <- rownames(x)
  colnames(out) <- colnames(x)
  out
}

apply_view_transforms <- function(mats, transforms) {
  if (length(transforms) == 0L) {
    return(mats)
  }
  missing_views <- setdiff(names(transforms), names(mats))
  if (length(missing_views) > 0L) {
    stop("Transform requested views not present in template: ", paste(missing_views, collapse = ", "))
  }
  for (view in names(transforms)) {
    transform <- transforms[[view]]
    mats[[view]] <- switch(
      transform,
      none = mats[[view]],
      identity = mats[[view]],
      log = log_transform(mats[[view]], log1p = FALSE),
      log1p = log_transform(mats[[view]], log1p = TRUE),
      clr = clr_transform(mats[[view]]),
      stop("Unsupported transform: ", transform)
    )
  }
  mats
}

add_task_metadata <- function(df, row, dataset, feature_counts) {
  df$job_id <- row$job_id[[1]]
  df$dataset <- dataset
  df$cancer <- dataset
  df$feature_setting <- row$feature_setting[[1]]
  df$feature_counts <- paste(paste(names(feature_counts), feature_counts, sep = "="), collapse = ",")
  if ("view_transforms" %in% names(row)) {
    df$view_transforms <- row$view_transforms[[1]]
  }
  df$n_subjects <- row$n_subjects[[1]]
  df$de_align_label <- row$de_align_label[[1]]
  if (!"crosstalk" %in% names(df)) {
    df$crosstalk <- row$crosstalk[[1]]
  }
  df$grid_rep <- row$rep[[1]]
  df$seed <- row$seed[[1]]
  df
}

grid_file <- path.expand(get_arg("grid", "~/Mimesys_validation/tcga_validation_grid/tcga_mimesys_validation_grid.csv"))
base_dir <- path.expand(get_arg("base_dir", "~/Mimesys_validation"))
out_dir <- path.expand(get_arg("out_dir", file.path(base_dir, "tcga_validation_grid")))
snr <- as.numeric(get_arg("snr", "2"))
figure_caps_arg <- get_arg("figure_caps", "")

task_id_arg <- get_arg("task_id", Sys.getenv("SLURM_ARRAY_TASK_ID", unset = NA_character_))
task_id <- as.integer(task_id_arg)
if (is.na(task_id)) {
  stop("No task id supplied. Use --task_id=... or run as a Slurm array job.")
}

suppressPackageStartupMessages({
  library(MimeSys)
})

grid <- read.csv(grid_file, stringsAsFactors = FALSE)
if (task_id < 1 || task_id > nrow(grid)) {
  stop("task_id is outside the grid range: ", task_id)
}

row <- grid[task_id, , drop = FALSE]
dataset <- if ("dataset" %in% names(row) && is_present(row$dataset[[1]])) {
  row$dataset[[1]]
} else {
  row$cancer[[1]]
}

task_dir <- file.path(out_dir, "tasks")
dir.create(task_dir, recursive = TRUE, showWarnings = FALSE)

out_rds <- file.path(task_dir, sprintf("task_%05d.rds", task_id))
out_csv <- file.path(task_dir, sprintf("task_%05d.csv", task_id))

if (file.exists(out_rds) && file.exists(out_csv) && task_output_current(out_rds)) {
  message("Task already complete; skipping: ", task_id)
  quit(save = "no", status = 0)
} else if (file.exists(out_rds) || file.exists(out_csv)) {
  message("Task output exists but is missing current block metrics; regenerating: ", task_id)
}

template_path <- if ("template_path" %in% names(row) && is_present(row$template_path[[1]])) {
  path.expand(row$template_path[[1]])
} else {
  file.path(base_dir, dataset, paste0("prepared_", dataset), paste0(dataset, "_multiomics_template_clean.rds"))
}
if (!file.exists(template_path)) {
  stop(
    "Missing template: ", template_path, "\n",
    "For TCGA, run run_tcga_mimesys_validation('", dataset, "', nrep=1, crosstalk_grid=0) once. ",
    "For other datasets, save a MultiAssayExperiment or named matrix-list RDS and pass --template_paths."
  )
}

template_raw <- extract_template_input(readRDS(template_path))
template_clean <- as_mimesys_view_list(template_raw)

view_transforms <- if ("view_transforms" %in% names(row) && is_present(row$view_transforms[[1]])) {
  parse_view_transforms(row$view_transforms[[1]])
} else {
  character()
}
template_clean <- apply_view_transforms(template_clean, view_transforms)

feature_counts <- if ("feature_counts" %in% names(row) && is_present(row$feature_counts[[1]])) {
  parse_counts(row$feature_counts[[1]], "feature_counts")
} else {
  old_cols <- intersect(c("mrna_features", "mirna_features", "protein_features"), names(row))
  old_values <- unlist(row[old_cols], use.names = FALSE)
  names(old_values) <- sub("_features$", "", old_cols)
  old_values[!is.na(old_values)]
}

if (length(feature_counts) == 0L) {
  feature_counts <- stats::setNames(vapply(template_clean, ncol, integer(1)), names(template_clean))
}
missing_views <- setdiff(names(feature_counts), names(template_clean))
if (length(missing_views) > 0L) {
  stop("Feature counts requested views not present in template: ", paste(missing_views, collapse = ", "))
}

template_input <- lapply(names(feature_counts), function(view) {
  keep_top_variable(template_clean[[view]], feature_counts[[view]])
})
names(template_input) <- names(feature_counts)

de_align <- switch(
  row$de_align_label[[1]],
  random_effects = NULL,
  aligned_effects = 1,
  stop("Unknown de_align_label: ", row$de_align_label[[1]])
)

message(sprintf(
  "Task %s | %s | %s | n=%s | %s | crosstalk=%s | rep=%s",
  task_id,
  dataset,
  row$feature_setting[[1]],
  row$n_subjects[[1]],
  row$de_align_label[[1]],
  row$crosstalk[[1]],
  row$rep[[1]]
))

tmpl <- fit_template(template_input)

sim <- gen_simmba(
  nsample = row$n_subjects[[1]],
  template = tmpl,
  crosstalk = row$crosstalk[[1]],
  de.align = de_align,
  nrep = 1,
  snr = snr,
  seed = row$seed[[1]]
)

diag <- validate_template_simulation(sim)
diag <- add_task_metadata(diag, row, dataset, feature_counts)

sim_mats <- simulation_view_matrices(
  sim,
  rep = 1,
  split = "train",
  feature_names = "template_feature"
)

figure_caps <- if (nzchar(figure_caps_arg)) {
  parse_counts(figure_caps_arg, "--figure_caps")
} else {
  stats::setNames(vapply(template_input, ncol, integer(1)), names(template_input))
}
figure_caps <- stats::setNames(
  vapply(names(template_input), function(view) {
    cap <- if (view %in% names(figure_caps)) figure_caps[[view]] else ncol(template_input[[view]])
    min(as.integer(cap), ncol(template_input[[view]]))
  }, integer(1)),
  names(template_input)
)

marginal_diag <- marginal_distribution_diagnostic(
  observed_matrices = template_input,
  simulated_matrices = sim_mats,
  max_features_per_view = figure_caps
)
marginal_diag <- add_task_metadata(marginal_diag, row, dataset, feature_counts)

pca_diag <- tryCatch(
  {
    out <- multiomics_pca_diagnostic(
      observed_matrices = template_input,
      simulated_matrices = sim_mats,
      max_features_per_view = figure_caps
    )
    out$scores$pc1_var_explained <- out$var_explained[[1]]
    out$scores$pc2_var_explained <- out$var_explained[[2]]
    out$scores
  },
  error = function(e) {
    warning("PCA diagnostic failed for task ", task_id, ": ", conditionMessage(e))
    data.frame(
      sample_id = character(),
      data_type = character(),
      PC1 = numeric(),
      PC2 = numeric(),
      pc1_var_explained = numeric(),
      pc2_var_explained = numeric(),
      stringsAsFactors = FALSE
    )
  }
)
pca_diag <- add_task_metadata(pca_diag, row, dataset, feature_counts)

modality_pca_diag <- tryCatch(
  {
    modality_pca_diagnostic(
      observed_matrices = template_input,
      simulated_matrices = sim_mats,
      max_features_per_view = figure_caps
    )
  },
  error = function(e) {
    warning("Modality-wise PCA diagnostic failed for task ", task_id, ": ", conditionMessage(e))
    data.frame(
      view = character(),
      sample_id = character(),
      data_type = character(),
      PC1 = numeric(),
      PC2 = numeric(),
      pc1_var_explained = numeric(),
      pc2_var_explained = numeric(),
      stringsAsFactors = FALSE
    )
  }
)
modality_pca_diag <- add_task_metadata(modality_pca_diag, row, dataset, feature_counts)

obs_corr <- blockwise_feature_correlation(
  template_input,
  max_features_per_view = figure_caps
)
sim_corr <- blockwise_feature_correlation(
  sim_mats,
  max_features_per_view = figure_caps
)
obs_block_summary <- summarize_blockwise_feature_correlation(obs_corr)
sim_block_summary <- summarize_blockwise_feature_correlation(sim_corr)
block_comparison <- do.call(rbind, lapply(required_block_metrics, function(metric) {
  compare_blockwise_correlation_summary(
    obs_block_summary,
    sim_block_summary,
    value = metric
  )
}))
block_comparison <- add_task_metadata(block_comparison, row, dataset, feature_counts)

block_structure <- compare_feature_correlation_structure(obs_corr, sim_corr, return_pairs = TRUE)
block_structure_agreement <- add_task_metadata(block_structure$summary, row, dataset, feature_counts)
block_structure_pairs <- add_task_metadata(block_structure$pairs, row, dataset, feature_counts)

out_bundle <- list(
  diagnostics = diag,
  marginal_diagnostics = marginal_diag,
  pca_diagnostics = pca_diag,
  modality_pca_diagnostics = modality_pca_diag,
  block_comparison = block_comparison,
  block_structure_agreement = block_structure_agreement,
  block_structure_pairs = block_structure_pairs
)

saveRDS(out_bundle, out_rds)
write.csv(diag, out_csv, row.names = FALSE)

message("Wrote: ", out_rds)
