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
  out <- trimws(strsplit(x, sep, fixed = TRUE)[[1]])
  out[nzchar(out)]
}

parse_key_values <- function(x, label, value_fun = identity) {
  if (is.null(x) || !nzchar(x)) {
    return(character())
  }
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
  if (anyDuplicated(names(out))) {
    stop(label, " contains duplicate names: ", paste(names(out)[duplicated(names(out))], collapse = ", "))
  }
  unlist(out, use.names = TRUE)
}

parse_feature_counts <- function(x, views, label) {
  if (grepl("=", x, fixed = TRUE)) {
    values <- parse_key_values(x, label, value_fun = as.integer)
    if (any(is.na(values)) || any(values <= 0)) {
      stop(label, " feature counts must be positive integers.", call. = FALSE)
    }
    return(values)
  }

  values <- as.integer(split_arg(x, ","))
  if (length(values) != length(views) || any(is.na(values)) || any(values <= 0)) {
    stop(
      label, " must contain either named counts like view1=100,view2=50 ",
      "or ", length(views), " positive integers in --views order.",
      call. = FALSE
    )
  }
  stats::setNames(values, views)
}

parse_feature_sets <- function(x, views) {
  pieces <- split_arg(x, ";")
  if (length(pieces) == 0L) {
    stop("--feature_sets cannot be empty.", call. = FALSE)
  }

  out <- vector("list", length(pieces))
  names(out) <- character(length(pieces))
  for (i in seq_along(pieces)) {
    parts <- strsplit(pieces[[i]], ":", fixed = TRUE)[[1]]
    if (length(parts) != 2L || !nzchar(parts[[1]]) || !nzchar(parts[[2]])) {
      stop(
        "--feature_sets entries must look like set:view=count,view=count. ",
        "Example: --feature_sets='small:mrna=500,protein=50;large:mrna=2000,protein=200'",
        call. = FALSE
      )
    }
    names(out)[[i]] <- trimws(parts[[1]])
    out[[i]] <- parse_feature_counts(trimws(parts[[2]]), views, paste0("feature set ", names(out)[[i]]))
  }

  if (anyDuplicated(names(out))) {
    stop("--feature_sets contains duplicate names: ", paste(names(out)[duplicated(names(out))], collapse = ", "))
  }
  out
}

format_feature_counts <- function(x) {
  paste(paste(names(x), as.integer(x), sep = "="), collapse = ",")
}

format_view_transforms <- function(x) {
  if (length(x) == 0L) {
    return("")
  }
  paste(paste(names(x), x, sep = "="), collapse = ",")
}

base_dir <- path.expand(get_arg("base_dir", "~/Mimesys_validation"))
out_dir <- path.expand(get_arg("out_dir", file.path(base_dir, "tcga_validation_grid")))
nrep <- as.integer(get_arg("nrep", "25"))

datasets_arg <- get_arg("datasets", get_arg("cancers", "HNSC,LUAD"))
datasets <- split_arg(datasets_arg, ",")
sample_sizes <- as.integer(split_arg(get_arg("sample_sizes", "100,200,400"), ","))
crosstalk_grid <- as.numeric(split_arg(get_arg("crosstalk", "0,0.25,0.5,0.75,1"), ","))
views <- split_arg(get_arg("views", "mrna,mirna,protein"), ",")

if (length(datasets) == 0L) stop("No datasets supplied.", call. = FALSE)
if (any(is.na(sample_sizes)) || any(sample_sizes <= 0)) stop("--sample_sizes must be positive integers.", call. = FALSE)
if (any(is.na(crosstalk_grid))) stop("--crosstalk must be numeric.", call. = FALSE)
if (is.na(nrep) || nrep <= 0) stop("--nrep must be a positive integer.", call. = FALSE)

template_paths <- parse_key_values(get_arg("template_paths", ""), "--template_paths", value_fun = path.expand)
single_template_path <- get_arg("template_path", NULL)
if (!is.null(single_template_path) && length(datasets) == 1L) {
  template_paths <- stats::setNames(path.expand(single_template_path), datasets)
}
view_transforms <- parse_key_values(get_arg("view_transforms", ""), "--view_transforms")

feature_sets_arg <- get_arg("feature_sets", NULL)
if (!is.null(feature_sets_arg)) {
  feature_sets <- parse_feature_sets(feature_sets_arg, views)
} else {
  feature_sets <- list(
    small = parse_feature_counts(get_arg("small_features", "500,100,50"), views, "--small_features"),
    medium = parse_feature_counts(get_arg("medium_features", "1000,250,100"), views, "--medium_features"),
    large = parse_feature_counts(get_arg("large_features", "2000,500,200"), views, "--large_features")
  )
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

grid <- expand.grid(
  dataset = datasets,
  feature_setting = names(feature_sets),
  n_subjects = sample_sizes,
  de_align_label = c("random_effects", "aligned_effects"),
  crosstalk = crosstalk_grid,
  rep = seq_len(nrep),
  stringsAsFactors = FALSE
)

grid <- grid[order(
  grid$dataset,
  grid$feature_setting,
  grid$n_subjects,
  grid$de_align_label,
  grid$crosstalk,
  grid$rep
), ]

grid$job_id <- seq_len(nrow(grid))
grid$seed <- 100000 + grid$job_id
grid$feature_counts <- vapply(feature_sets[grid$feature_setting], format_feature_counts, character(1))
grid$view_transforms <- format_view_transforms(view_transforms)
grid$template_path <- ifelse(
  grid$dataset %in% names(template_paths),
  unname(template_paths[grid$dataset]),
  ""
)

# Backward-compatible alias for older TCGA analysis notebooks.
grid$cancer <- grid$dataset

grid <- grid[
  ,
  c(
    "job_id", "dataset", "cancer", "template_path", "feature_setting",
    "feature_counts", "view_transforms", "n_subjects", "de_align_label", "crosstalk", "rep", "seed"
  )
]

grid_file <- file.path(out_dir, "tcga_mimesys_validation_grid.csv")
write.csv(grid, grid_file, row.names = FALSE)

cat("Wrote grid:\n")
cat(grid_file, "\n")
cat("Rows:", nrow(grid), "\n")
cat("Datasets:", paste(datasets, collapse = ", "), "\n")
cat("Feature sets:\n")
for (nm in names(feature_sets)) {
  cat("  ", nm, ": ", format_feature_counts(feature_sets[[nm]]), "\n", sep = "")
}
if (length(view_transforms) > 0L) {
  cat("View transforms:", format_view_transforms(view_transforms), "\n")
}
cat("Suggested batches of 1000:", ceiling(nrow(grid) / 1000), "\n")
