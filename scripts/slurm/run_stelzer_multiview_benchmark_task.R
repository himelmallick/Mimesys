#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) == 0L) return(default)
  sub(paste0("^--", name, "="), "", hit[[1]])
}

split_arg <- function(x, sep = ",") {
  if (is.null(x) || length(x) == 0L || is.na(x[[1]])) return(character())
  x <- as.character(x[[1]])
  if (!nzchar(x)) return(character())
  out <- trimws(strsplit(x, sep, fixed = TRUE)[[1]])
  out[nzchar(out)]
}

parse_key_values <- function(x, label, value_fun = as.numeric) {
  if (is.null(x) || length(x) == 0L || is.na(x[[1]])) return(NULL)
  x <- as.character(x[[1]])
  if (!nzchar(x)) return(NULL)
  parts <- split_arg(x)
  out <- vapply(parts, function(part) {
    kv <- strsplit(part, "=", fixed = TRUE)[[1]]
    if (length(kv) != 2L || !nzchar(kv[[1]])) {
      stop("Invalid ", label, " entry: ", part, call. = FALSE)
    }
    value_fun(kv[[2]])
  }, numeric(1))
  names(out) <- vapply(parts, function(part) strsplit(part, "=", fixed = TRUE)[[1]][[1]], character(1))
  out
}

has_value <- function(x) {
  !is.null(x) && length(x) > 0L && !is.na(x[[1]]) && nzchar(as.character(x[[1]]))
}

load_mimesys <- function() {
  project_dir <- Sys.getenv("MIMESYS_PROJECT_DIR", unset = getwd())
  if (requireNamespace("devtools", quietly = TRUE) && file.exists(file.path(project_dir, "DESCRIPTION"))) {
    devtools::load_all(project_dir, quiet = TRUE)
  } else {
    suppressPackageStartupMessages(library(MimeSys))
  }
}

keep_top_variable <- function(x, n_features) {
  vars <- apply(x, 2, stats::var, na.rm = TRUE)
  vars <- vars[is.finite(vars)]
  vars <- sort(vars, decreasing = TRUE)
  n_keep <- min(as.integer(n_features), length(vars))
  if (n_keep <= 0L) stop("No variable features available.", call. = FALSE)
  x[, names(vars)[seq_len(n_keep)], drop = FALSE]
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

prep_view_matrices <- function(train_x, test_x) {
  means <- colMeans(train_x, na.rm = TRUE)
  sds <- apply(train_x, 2, stats::sd, na.rm = TRUE)
  sds[!is.finite(sds) | sds == 0] <- 1

  train_scaled <- sweep(train_x, 2, means, "-")
  train_scaled <- sweep(train_scaled, 2, sds, "/")
  test_scaled <- sweep(test_x, 2, means, "-")
  test_scaled <- sweep(test_scaled, 2, sds, "/")
  train_scaled[!is.finite(train_scaled)] <- 0
  test_scaled[!is.finite(test_scaled)] <- 0
  list(train = train_scaled, test = test_scaled)
}

prep_views <- function(train_views, test_views) {
  out <- lapply(names(train_views), function(view) {
    prep_view_matrices(train_views[[view]], test_views[[view]])
  })
  names(out) <- names(train_views)
  list(train = lapply(out, `[[`, "train"), test = lapply(out, `[[`, "test"))
}

r2_score <- function(y, pred) {
  ok <- is.finite(y) & is.finite(pred)
  y <- y[ok]
  pred <- pred[ok]
  if (length(y) < 2L) return(NA_real_)
  denom <- sum((y - mean(y))^2)
  if (!is.finite(denom) || denom <= 0) return(NA_real_)
  1 - sum((y - pred)^2) / denom
}

make_foldid <- function(n, folds) {
  sample(rep(seq_len(folds), length.out = n))
}

fit_cv_multiview <- function(x_train, y_train, rho, foldid) {
  multiview::cv.multiview(
    x_list = x_train,
    y = y_train,
    rho = rho,
    family = stats::gaussian(),
    type.measure = "mse",
    nfolds = length(unique(foldid)),
    foldid = foldid
  )
}

cv_error_min <- function(cvfit) {
  suppressWarnings(min(cvfit$cvm, na.rm = TRUE))
}

fit_eval_fixed_rho <- function(x_train, y_train, x_test, y_test, rho, foldid, method_label) {
  fit <- fit_cv_multiview(x_train, y_train, rho = rho, foldid = foldid)
  pred <- as.numeric(stats::predict(fit, newx = x_test, s = "lambda.min"))
  data.frame(
    method = method_label,
    rho = rho,
    selected_rho = rho,
    lambda_min = fit$lambda.min,
    cv_mse = cv_error_min(fit),
    test_r2 = r2_score(y_test, pred),
    stringsAsFactors = FALSE
  )
}

fit_eval_cv_rho <- function(x_train, y_train, x_test, y_test, rho_grid, foldid) {
  fits <- lapply(rho_grid, function(rho) {
    tryCatch(fit_cv_multiview(x_train, y_train, rho = rho, foldid = foldid), error = function(e) e)
  })
  cv_mse <- vapply(fits, function(fit) {
    if (inherits(fit, "error")) return(Inf)
    cv_error_min(fit)
  }, numeric(1))
  if (all(!is.finite(cv_mse))) {
    return(data.frame(
      method = "rho_cv", rho = NA_real_, selected_rho = NA_real_,
      lambda_min = NA_real_, cv_mse = NA_real_, test_r2 = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  best_idx <- which.min(cv_mse)
  best_fit <- fits[[best_idx]]
  pred <- as.numeric(stats::predict(best_fit, newx = x_test, s = "lambda.min"))
  data.frame(
    method = "rho_cv",
    rho = NA_real_,
    selected_rho = rho_grid[[best_idx]],
    lambda_min = best_fit$lambda.min,
    cv_mse = cv_mse[[best_idx]],
    test_r2 = r2_score(y_test, pred),
    stringsAsFactors = FALSE
  )
}

grid_file <- path.expand(get_arg("grid", "~/Mimesys_validation/stelzer_multiview_cooperative_benchmark/stelzer_multiview_benchmark_grid.csv"))
out_dir <- path.expand(get_arg("out_dir", "~/Mimesys_validation/stelzer_multiview_cooperative_benchmark"))
data_path <- path.expand(get_arg("data_path", "~/Mimesys_2.0/inst/extdata/StelzerDOS.RData"))
task_id_arg <- get_arg("task_id", Sys.getenv("SLURM_ARRAY_TASK_ID", unset = NA_character_))
task_id <- as.integer(task_id_arg)
if (is.na(task_id)) stop("No task id supplied. Use --task_id=... or Slurm array.", call. = FALSE)

load_mimesys()
if (!requireNamespace("multiview", quietly = TRUE)) {
  stop("This benchmark requires the multiview package.", call. = FALSE)
}

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
views <- lapply(views_full, keep_top_variable, n_features = row$feature_n[[1]])
tmpl <- fit_template(views)

scenario_defs <- list(
  balanced_aligned = list(signal_share = c(Immunome = 1, Metabolome = 1, Proteome = 1), de_align = 1),
  balanced_random = list(signal_share = c(Immunome = 1, Metabolome = 1, Proteome = 1), de_align = NULL)
)
scenario <- row$scenario[[1]]
if (!scenario %in% names(scenario_defs)) stop("Unknown scenario: ", scenario, call. = FALSE)
scen <- scenario_defs[[scenario]]
rho_grid <- as.numeric(split_arg(row$rho_grid[[1]]))
friedman_share <- if ("friedman_share" %in% names(row) && has_value(row$friedman_share)) {
  parse_key_values(row$friedman_share[[1]], "friedman_share", as.numeric)
} else {
  NULL
}

set.seed(row$seed[[1]])
sim <- gen_simmba(
  nsample = row$nsample[[1]],
  template = tmpl,
  crosstalk = row$crosstalk[[1]],
  ygen.mode = if ("ygen_mode" %in% names(row)) row$ygen_mode[[1]] else "LM",
  outcome.type = "continuous",
  de.align = scen$de_align,
  signal.share = scen$signal_share,
  friedman.share = friedman_share,
  dependence = row$dependence[[1]],
  de.prob = row$de_prob[[1]],
  de.facLoc = rep(row$de_fac_loc[[1]], length(tmpl$views)),
  de.facScale = rep(row$de_fac_scale[[1]], length(tmpl$views)),
  snr = row$snr[[1]],
  nrep = 1,
  seed = row$seed[[1]]
)

x_train <- simulation_view_matrices(sim, rep = 1, split = "train", feature_names = "template_feature")
x_test <- simulation_view_matrices(sim, rep = 1, split = "test", feature_names = "template_feature")
y_train <- sim$trainDat[[1]]$sample_metadata$Y
y_test <- sim$testDat[[1]]$sample_metadata$Y
x_prepped <- prep_views(x_train, x_test)
foldid <- make_foldid(length(y_train), row$folds[[1]])

bench <- tryCatch(
  rbind(
    fit_eval_fixed_rho(x_prepped$train, y_train, x_prepped$test, y_test, rho = 0, foldid = foldid, method_label = "rho_0"),
    fit_eval_fixed_rho(x_prepped$train, y_train, x_prepped$test, y_test, rho = 1, foldid = foldid, method_label = "rho_1"),
    fit_eval_cv_rho(x_prepped$train, y_train, x_prepped$test, y_test, rho_grid = rho_grid, foldid = foldid)
  ),
  error = function(e) {
    data.frame(
      method = c("rho_0", "rho_1", "rho_cv"),
      rho = c(0, 1, NA_real_),
      selected_rho = c(0, 1, NA_real_),
      lambda_min = NA_real_,
      cv_mse = NA_real_,
      test_r2 = NA_real_,
      error = conditionMessage(e),
      stringsAsFactors = FALSE
    )
  }
)
if (!"error" %in% names(bench)) bench$error <- NA_character_

meta <- data.frame(
  job_id = row$job_id[[1]],
  dataset = row$dataset[[1]],
  scenario = row$scenario[[1]],
  ygen_mode = if ("ygen_mode" %in% names(row)) row$ygen_mode[[1]] else "LM",
  feature_n = row$feature_n[[1]],
  nsample = row$nsample[[1]],
  crosstalk = row$crosstalk[[1]],
  grid_rep = row$rep[[1]],
  seed = row$seed[[1]],
  rho_grid = row$rho_grid[[1]],
  folds = row$folds[[1]],
  snr = row$snr[[1]],
  de_prob = row$de_prob[[1]],
  de_fac_loc = row$de_fac_loc[[1]],
  de_fac_scale = row$de_fac_scale[[1]],
  dependence = row$dependence[[1]],
  friedman_share = if (is.null(friedman_share)) "" else paste(paste(names(friedman_share), friedman_share, sep = "="), collapse = ","),
  stringsAsFactors = FALSE
)
out <- cbind(meta[rep(1, nrow(bench)), , drop = FALSE], bench)
saveRDS(list(grid_row = row, benchmark = out), out_rds)
write.csv(out, out_csv, row.names = FALSE)
message("Wrote ", out_csv)
