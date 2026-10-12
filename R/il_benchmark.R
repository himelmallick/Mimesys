#' Benchmark IntegratedLearner on a MimeSys simulation
#'
#' Converts each train/test replicate from `gen_simmba()` into
#' `MultiAssayExperiment` objects, runs `IntegratedLearner()`, and returns
#' predictive performance for single-view and fusion learners.
#'
#' @param sim A list returned by template-based `gen_simmba()`.
#' @param outcome_type Either `"continuous"` or `"binary"`.
#' @param folds Number of cross-validation folds used inside
#'   `IntegratedLearner()`.
#' @param base_learner Base learner passed to `IntegratedLearner()`.
#' @param meta_learner Meta learner passed to `IntegratedLearner()` for stacked
#'   late fusion.
#' @param run_concat,run_stacked Logical flags passed to `IntegratedLearner()`.
#' @param extra_metadata Optional one-row data frame of simulation setting
#'   metadata to attach to every output row.
#'
#' @return A data frame with one row per replicate and learner.
#' @export
benchmark_integratedlearner_simulation <- function(
  sim,
  outcome_type = c("continuous", "binary"),
  folds = 2,
  base_learner = "SL.randomForest",
  meta_learner = "SL.nnls.auc",
  run_concat = TRUE,
  run_stacked = TRUE,
  extra_metadata = NULL
) {
  outcome_type <- match.arg(outcome_type)

  for (pkg in c("IntegratedLearner", "MultiAssayExperiment", "S4Vectors", "SummarizedExperiment")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("benchmark_integratedlearner_simulation() requires package '%s'.", pkg), call. = FALSE)
    }
  }
  if (!is.list(sim) || is.null(sim$trainDat) || is.null(sim$testDat)) {
    stop("sim must be a list returned by gen_simmba().", call. = FALSE)
  }

  family <- switch(
    outcome_type,
    continuous = stats::gaussian(),
    binary = stats::binomial()
  )
  metric_name <- switch(
    outcome_type,
    continuous = "R2.test",
    binary = "AUC.test"
  )

  pcl_to_mae <- function(dat) {
    fm <- dat$feature_metadata
    sm <- dat$sample_metadata
    ft <- as.matrix(dat$feature_table)
    if (is.null(sm$subjectID)) {
      sm$subjectID <- rownames(sm)
    }
    sm <- S4Vectors::DataFrame(sm)

    views <- unique(fm$featureType)
    experiments <- lapply(views, function(view) {
      sel <- fm$featureType == view
      assay <- ft[sel, , drop = FALSE]
      SummarizedExperiment::SummarizedExperiment(
        assays = list(data = assay),
        colData = sm
      )
    })
    names(experiments) <- views
    MultiAssayExperiment::MultiAssayExperiment(
      experiments = experiments,
      colData = sm
    )
  }

  extract_metric <- function(fit) {
    values <- switch(
      outcome_type,
      continuous = fit$R2.test,
      binary = fit$AUC.test
    )
    if (is.null(values)) {
      return(data.frame(
        metric = metric_name,
        learner = NA_character_,
        performance = NA_real_,
        best_single = NA_real_,
        gain_over_best_single = NA_real_,
        stringsAsFactors = FALSE
      ))
    }

    values <- unlist(values)
    single_idx <- !grepl(
      "stack|concat|early|late|intermediate|integrat|ensemble",
      names(values),
      ignore.case = TRUE
    )
    best_single <- max(values[single_idx], na.rm = TRUE)

    learner <- names(values)
    learner[learner == "concatenated"] <- "early_fusion"
    learner[learner == "stacked"] <- "late_fusion"

    data.frame(
      metric = metric_name,
      learner = learner,
      performance = as.numeric(values),
      best_single = best_single,
      gain_over_best_single = as.numeric(values) - best_single,
      stringsAsFactors = FALSE
    )
  }

  rows <- vector("list", length(sim$trainDat))
  for (rep_idx in seq_along(sim$trainDat)) {
    fit <- tryCatch(
      suppressWarnings(IntegratedLearner::IntegratedLearner(
        MAE_train = pcl_to_mae(sim$trainDat[[rep_idx]]),
        MAE_valid = pcl_to_mae(sim$testDat[[rep_idx]]),
        experiment = sim$template$views,
        assay.type = rep("data", length(sim$template$views)),
        outcome_col = "Y",
        subject_id_col = "subjectID",
        folds = folds,
        base_learner = base_learner,
        meta_learner = meta_learner,
        run_concat = run_concat,
        run_stacked = run_stacked,
        print_learner = FALSE,
        family = family
      )),
      error = function(e) e
    )

    if (inherits(fit, "error")) {
      out <- data.frame(
        metric = metric_name,
        learner = "fit_failed",
        performance = NA_real_,
        best_single = NA_real_,
        gain_over_best_single = NA_real_,
        error = conditionMessage(fit),
        stringsAsFactors = FALSE
      )
    } else {
      out <- extract_metric(fit)
      out$error <- NA_character_
    }
    out$rep <- rep_idx
    rows[[rep_idx]] <- out
  }

  out <- do.call(rbind, rows)
  if (!is.null(extra_metadata)) {
    extra_metadata <- as.data.frame(extra_metadata, stringsAsFactors = FALSE)
    extra_metadata <- extra_metadata[rep(1, nrow(out)), , drop = FALSE]
    rownames(extra_metadata) <- NULL
    out <- cbind(extra_metadata, out)
  }
  rownames(out) <- NULL
  out
}

#' Run IntegratedLearner benchmark using cached TCGA MimeSys templates
#'
#' Loads TCGA templates created by `run_tcga_mimesys_validation()`, simulates
#' continuous and/or binary outcomes under a crosstalk/effect-alignment grid,
#' runs IntegratedLearner, and saves raw and summarized benchmark results.
#'
#' @param cancers Character vector of TCGA cancer codes.
#' @param base_dir Base directory containing cancer-specific MimeSys validation
#'   folders.
#' @param outcome_types Outcomes to benchmark: `"continuous"` and/or `"binary"`.
#' @param feature_setting Which feature-size setting to use.
#' @param feature_sets Named list of feature counts per view.
#' @param crosstalk_grid Numeric crosstalk values.
#' @param nrep Number of simulated replicates per setting.
#' @param snr Signal-to-noise ratio for `gen_simmba()`.
#' @param folds Number of folds used inside `IntegratedLearner()`.
#' @param base_learner Base learner passed to `IntegratedLearner()`.
#' @param meta_learner Meta learner passed to `IntegratedLearner()` for stacked
#'   late fusion.
#' @param seed Base seed.
#'
#' @return A list containing raw results, summary table, and output paths.
#' @export
run_tcga_integratedlearner_benchmark <- function(
  cancers = c("HNSC", "LUAD"),
  base_dir = path.expand("~/Mimesys_validation"),
  outcome_types = c("continuous", "binary"),
  feature_setting = "medium",
  feature_sets = list(
    small = c(mrna = 500, mirna = 100, protein = 50),
    medium = c(mrna = 1000, mirna = 250, protein = 100),
    large = c(mrna = 2000, mirna = 500, protein = 200)
  ),
  crosstalk_grid = c(0, 0.25, 0.5, 0.75, 1),
  nrep = 5,
  snr = 2,
  folds = 2,
  base_learner = "SL.randomForest",
  meta_learner = "SL.nnls.auc",
  seed = 5000
) {
  outcome_types <- match.arg(outcome_types, choices = c("continuous", "binary"), several.ok = TRUE)
  if (!feature_setting %in% names(feature_sets)) {
    stop("feature_setting must be one of: ", paste(names(feature_sets), collapse = ", "), call. = FALSE)
  }

  base_dir <- path.expand(base_dir)
  de_align_grid <- list(random_effects = NULL, aligned_effects = 1)
  all_results <- list()
  counter <- 1L

  for (cancer in toupper(cancers)) {
    cancer_dir <- file.path(base_dir, cancer)
    cache_dir <- file.path(cancer_dir, paste0("prepared_", cancer))
    out_dir <- file.path(cancer_dir, "integratedlearner_benchmark")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    template_path <- file.path(cache_dir, paste0(cancer, "_multiomics_template_clean.rds"))
    if (!file.exists(template_path)) {
      stop(
        "Missing clean template for ", cancer, ": ", template_path, "\n",
        "Run run_tcga_mimesys_validation('", cancer, "') first.",
        call. = FALSE
      )
    }
    template_clean <- readRDS(template_path)
    counts <- feature_sets[[feature_setting]]
    template_input <- list(
      mrna = .keep_top_variable_for_il(template_clean$mrna, counts[["mrna"]]),
      mirna = .keep_top_variable_for_il(template_clean$mirna, counts[["mirna"]]),
      protein = .keep_top_variable_for_il(template_clean$protein, counts[["protein"]])
    )
    tmpl <- fit_template(template_input)

    for (outcome_type in outcome_types) {
      for (de_name in names(de_align_grid)) {
        for (ct in crosstalk_grid) {
          message(sprintf(
            "IL benchmark | %s | %s | %s | crosstalk=%s",
            cancer, outcome_type, de_name, ct
          ))
          sim <- gen_simmba(
            nsample = nrow(template_input[[1]]),
            template = tmpl,
            crosstalk = ct,
            de.align = de_align_grid[[de_name]],
            nrep = nrep,
            snr = snr,
            outcome.type = outcome_type,
            seed = seed + counter
          )
          diagnostics <- validate_template_simulation(sim)

          extra <- data.frame(
            cancer = cancer,
            outcome_type = outcome_type,
            feature_setting = feature_setting,
            crosstalk = ct,
            de_align_label = de_name,
            mean_signal_agreement = mean(diagnostics$signal_agreement, na.rm = TRUE),
            mean_feature_cor = mean(diagnostics$mean_abs_feature_cor, na.rm = TRUE),
            y_signal_cor = mean(diagnostics$y_signal_cor, na.rm = TRUE),
            stringsAsFactors = FALSE
          )
          result <- benchmark_integratedlearner_simulation(
            sim = sim,
            outcome_type = outcome_type,
            folds = folds,
            base_learner = base_learner,
            meta_learner = meta_learner,
            extra_metadata = extra
          )
          result$setting_id <- counter
          all_results[[counter]] <- result
          counter <- counter + 1L
        }
      }
    }
  }

  raw <- do.call(rbind, all_results)
  out_dir <- file.path(base_dir, "integratedlearner_benchmark_all")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  raw_csv <- file.path(out_dir, "integratedlearner_benchmark_raw.csv")
  summary_csv <- file.path(out_dir, "integratedlearner_benchmark_summary.csv")
  rds_file <- file.path(out_dir, "integratedlearner_benchmark.rds")
  utils::write.csv(raw, raw_csv, row.names = FALSE)

  raw_summary <- raw[
    !is.na(raw$learner) &
      raw$learner != "fit_failed" &
      !is.na(raw$metric) &
      is.finite(raw$performance),
    ,
    drop = FALSE
  ]
  if (nrow(raw_summary) == 0) {
    summary <- data.frame()
  } else {
    summary <- stats::aggregate(
      cbind(performance, gain_over_best_single) ~
        cancer + outcome_type + feature_setting + de_align_label + crosstalk + learner + metric,
      data = raw_summary,
      FUN = function(x) mean(x, na.rm = TRUE)
    )
  }

  utils::write.csv(summary, summary_csv, row.names = FALSE)
  saveRDS(list(raw = raw, summary = summary), rds_file)

  list(
    raw = raw,
    summary = summary,
    files = list(raw_csv = raw_csv, summary_csv = summary_csv, rds = rds_file)
  )
}

.keep_top_variable_for_il <- function(x, n_features) {
  vars <- apply(x, 2, stats::var, na.rm = TRUE)
  vars <- vars[is.finite(vars)]
  vars <- sort(vars, decreasing = TRUE)
  n_keep <- min(n_features, length(vars))
  if (n_keep == 0) {
    stop("No variable features available.", call. = FALSE)
  }
  x[, names(vars)[seq_len(n_keep)], drop = FALSE]
}
