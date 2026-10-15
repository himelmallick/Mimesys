#' Convert a MimeSys simulation replicate to a MultiAssayExperiment
#'
#' Converts one replicate from [gen_simmba()] to a
#' [MultiAssayExperiment::MultiAssayExperiment()]. Simulated assays follow the
#' Bioconductor convention: features in rows and samples in columns. Outcome,
#' train/test split, and simulated sample identifiers are stored in
#' `colData()`. Feature-level simulation truth, including spike-in coefficients,
#' is stored in the `rowData()` of each assay.
#'
#' @param x A simulation returned by [gen_simmba()].
#' @param rep Replicate to convert.
#' @param split Which samples to include: `"all"`, `"train"`, or `"test"`.
#'
#' @return A `MultiAssayExperiment`.
#' @export
as_mimesys_mae <- function(x, rep = 1, split = c("all", "train", "test")) {
  if (!requireNamespace("MultiAssayExperiment", quietly = TRUE) ||
      !requireNamespace("SummarizedExperiment", quietly = TRUE) ||
      !requireNamespace("S4Vectors", quietly = TRUE)) {
    stop(
      "as_mimesys_mae() requires MultiAssayExperiment, SummarizedExperiment, and S4Vectors.",
      call. = FALSE
    )
  }
  split <- match.arg(split)
  if (!is.list(x$trainDat) || !is.list(x$testDat) || length(x$trainDat) < rep || rep < 1L) {
    stop("x must be a gen_simmba() result containing the requested replicate.", call. = FALSE)
  }

  train <- x$trainDat[[rep]]
  test <- x$testDat[[rep]]
  dat <- switch(split,
    train = train,
    test = test,
    all = .mimesys_bind_simulation_splits(train, test)
  )
  .mimesys_replicate_to_mae(dat, simulation = .mimesys_simulation_settings(x), split = split)
}

.mimesys_bind_simulation_splits <- function(train, test) {
  out <- train
  out$feature_table <- cbind(train$feature_table, test$feature_table)
  out$sample_metadata <- rbind(train$sample_metadata, test$sample_metadata)
  out$sample_metadata$split <- c(
    rep("train", nrow(train$sample_metadata)),
    rep("test", nrow(test$sample_metadata))
  )
  out$truth <- list(
    scores = rbind(train$truth$scores, test$truth$scores),
    shared = rbind(train$truth$shared, test$truth$shared)
  )
  out
}

.mimesys_replicate_to_mae <- function(dat, simulation = list(), split = "all") {
  feature_table <- as.matrix(dat$feature_table)
  feature_metadata <- as.data.frame(dat$feature_metadata, stringsAsFactors = FALSE)
  sample_metadata <- as.data.frame(dat$sample_metadata, stringsAsFactors = FALSE)
  if (nrow(feature_metadata) != nrow(feature_table) ||
      nrow(sample_metadata) != ncol(feature_table)) {
    stop("Simulation metadata dimensions do not match the simulated feature table.", call. = FALSE)
  }
  if (is.null(feature_metadata$featureType)) {
    stop("Simulation feature metadata must contain featureType.", call. = FALSE)
  }

  sample_ids <- sample_metadata$subjectID
  if (is.null(sample_ids)) sample_ids <- colnames(feature_table)
  sample_ids <- make.unique(as.character(sample_ids), sep = "_")
  rownames(sample_metadata) <- sample_ids
  colnames(feature_table) <- sample_ids
  if (identical(split, "train") || identical(split, "test")) {
    sample_metadata$split <- split
  } else if (is.null(sample_metadata$split)) {
    sample_metadata$split <- NA_character_
  }

  view_index <- split(seq_len(nrow(feature_metadata)), as.character(feature_metadata$featureType))
  experiments <- lapply(names(view_index), function(view) {
    idx <- view_index[[view]]
    assay <- feature_table[idx, , drop = FALSE]
    row_data <- feature_metadata[idx, , drop = FALSE]
    feature_ids <- row_data$featureID
    if (is.null(feature_ids)) feature_ids <- rownames(assay)
    feature_ids <- make.unique(as.character(feature_ids), sep = "_")
    rownames(assay) <- feature_ids
    rownames(row_data) <- feature_ids
    SummarizedExperiment::SummarizedExperiment(
      assays = list(simulated = assay),
      rowData = S4Vectors::DataFrame(row_data),
      colData = S4Vectors::DataFrame(sample_metadata)
    )
  })
  names(experiments) <- names(view_index)

  sample_map <- do.call(rbind, lapply(names(experiments), function(view) {
    S4Vectors::DataFrame(
      assay = factor(view, levels = names(experiments)),
      primary = sample_ids,
      colname = sample_ids
    )
  }))
  truth <- dat$truth
  if (!is.null(truth$shared)) rownames(truth$shared) <- sample_ids
  if (!is.null(truth$scores)) rownames(truth$scores) <- sample_ids

  MultiAssayExperiment::MultiAssayExperiment(
    experiments = MultiAssayExperiment::ExperimentList(experiments),
    colData = S4Vectors::DataFrame(sample_metadata),
    sampleMap = sample_map,
    metadata = list(
      simulation = simulation,
      truth = truth
    )
  )
}

.mimesys_simulation_settings <- function(x) {
  drop <- c("trainDat", "testDat", "template", "agreement.emp")
  x[setdiff(names(x), drop)]
}

.mimesys_merge_options <- function(values, defaults, name) {
  if (is.null(values)) return(defaults)
  if (!is.list(values)) {
    stop(name, " must be a named list.", call. = FALSE)
  }
  unknown <- setdiff(names(values), names(defaults))
  if (length(unknown) > 0L) {
    stop(name, " contains unsupported option(s): ", paste(unknown, collapse = ", "), ".", call. = FALSE)
  }
  utils::modifyList(defaults, values)
}

.mimesys_outcome_mode <- function(model) {
  lookup <- c(
    linear = "LM", lm = "LM", friedman = "Friedman", friedman2 = "Friedman2",
    shared_latent = "SharedLatent", sharedlatent = "SharedLatent"
  )
  model <- tolower(as.character(model))
  if (length(model) != 1L || !model %in% names(lookup)) {
    stop("outcome$model must be one of: linear, friedman, friedman2, shared_latent.", call. = FALSE)
  }
  unname(lookup[[model]])
}

.mimesys_preset_options <- function(preset) {
  presets <- list(
    default = list(
      description = "Aligned linear outcome with factor-based null simulation.",
      outcome = list(model = "linear"), signal = list(alignment = 1),
      null_model = list(dependence = "factor")
    ),
    linear = list(
      description = "Aligned additive linear outcome.",
      outcome = list(model = "linear"), signal = list(alignment = 1),
      null_model = list(dependence = "factor")
    ),
    friedman = list(
      description = "Nonlinear Friedman outcome from spiked observed features.",
      outcome = list(model = "friedman"), signal = list(alignment = 1),
      null_model = list(dependence = "factor")
    ),
    friedman2 = list(
      description = "Additive linear plus nonlinear Friedman outcome.",
      outcome = list(model = "friedman2"), signal = list(alignment = 1),
      null_model = list(dependence = "factor")
    ),
    shared_latent = list(
      description = "Linear outcome from all selected shared latent factors.",
      outcome = list(model = "shared_latent", shared_latent_function = "linear"),
      signal = list(alignment = 1), null_model = list(dependence = "factor")
    ),
    shared_latent_friedman = list(
      description = "Friedman outcome from the first five selected shared latent factors.",
      outcome = list(model = "shared_latent", shared_latent_function = "friedman"),
      signal = list(alignment = 1), null_model = list(dependence = "factor")
    ),
    shared_latent_friedman2 = list(
      description = "Linear plus Friedman outcome from selected shared latent factors.",
      outcome = list(model = "shared_latent", shared_latent_function = "friedman2"),
      signal = list(alignment = 1), null_model = list(dependence = "factor")
    )
  )
  if (!preset %in% names(presets)) {
    stop("Unknown preset: ", preset, ". Use mimesys_presets() to list available presets.", call. = FALSE)
  }
  presets[[preset]]
}

#' List built-in MimeSys simulation presets
#'
#' @return A data frame describing built-in presets for [simulate_mimesys()].
#' @export
mimesys_presets <- function() {
  preset_names <- c(
    "default", "linear", "friedman", "friedman2", "shared_latent",
    "shared_latent_friedman", "shared_latent_friedman2"
  )
  settings <- lapply(preset_names, .mimesys_preset_options)
  data.frame(
    preset = preset_names,
    outcome_model = vapply(settings, function(x) x$outcome$model, character(1)),
    shared_latent_function = vapply(
      settings,
      function(x) if (is.null(x$outcome$shared_latent_function)) NA_character_ else x$outcome$shared_latent_function,
      character(1)
    ),
    description = vapply(settings, `[[`, character(1), "description"),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

.mimesys_resolve_crosstalk <- function(tmpl, n, p, crosstalk, calibration,
                                       null_model, seed) {
  if (is.numeric(crosstalk)) {
    return(list(value = crosstalk, calibration = NULL))
  }
  if (!is.character(crosstalk) || length(crosstalk) != 1L ||
      !identical(tolower(crosstalk), "calibrate")) {
    stop("crosstalk must be numeric or \"calibrate\".", call. = FALSE)
  }
  calibration_defaults <- list(
    grid = seq(0, 1, by = 0.25), nrep = 5, nsample = NULL,
    max_features_per_view = NULL
  )
  calibration <- .mimesys_merge_options(calibration, calibration_defaults, "calibration")
  fit <- calibrate_crosstalk(
    template = tmpl,
    crosstalk_grid = calibration$grid,
    nsample = calibration$nsample,
    nrep = calibration$nrep,
    p = p,
    feature.select = null_model$feature_select,
    feature.jitter = null_model$feature_jitter,
    max_features_per_view = calibration$max_features_per_view,
    dependence = null_model$dependence,
    cor.shrinkage = null_model$cor_shrinkage,
    module.count = null_model$module_count,
    module.between.shrinkage = null_model$module_between_shrinkage,
    view.models = null_model$view_models,
    seed = seed
  )
  list(value = fit$best_crosstalk, calibration = fit)
}

#' Simulate null multi-omics data as a MultiAssayExperiment
#'
#' Fits or accepts a MimeSys template, calibrates crosstalk by default using
#' empirical cross-view correlation blocks, and returns simulated null feature
#' data in a [MultiAssayExperiment::MultiAssayExperiment()]. No spike-in
#' coefficients or outcomes are generated.
#'
#' @inheritParams simulate_mimesys
#' @param calibration Named list controlling automatic crosstalk calibration.
#'   Options are `grid`, `nrep`, `nsample`, and `max_features_per_view`.
#'
#' @return A null-data `MultiAssayExperiment`. The selected crosstalk value and
#'   full calibration result are stored in `metadata(x)$simulation`.
#' @export
simulate_mimesys_null <- function(template,
                                  n,
                                  p = NULL,
                                  crosstalk = "calibrate",
                                  null_model = list(),
                                  template_args = list(),
                                  calibration = list(),
                                  seed = 1234) {
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 2L) {
    stop("n must be a single integer of at least 2.", call. = FALSE)
  }
  null_defaults <- list(
    dependence = "factor", feature_select = "random", feature_jitter = 0.2,
    cor_shrinkage = 0.05, module_count = NULL,
    module_between_shrinkage = 0.75, view_models = NULL
  )
  null_model <- .mimesys_merge_options(null_model, null_defaults, "null_model")
  if (!is.list(template_args)) stop("template_args must be a named list.", call. = FALSE)
  tmpl <- if (inherits(template, "mimesys_template")) {
    template
  } else {
    do.call(fit_template, c(list(template = template), template_args))
  }
  resolved <- .mimesys_resolve_crosstalk(
    tmpl, n = as.integer(n), p = p, crosstalk = crosstalk,
    calibration = calibration, null_model = null_model, seed = seed
  )
  pcl <- sim_from_template(
    tmpl = tmpl, n = as.integer(n), p = p, crosstalk = resolved$value,
    feature.select = null_model$feature_select,
    feature.jitter = null_model$feature_jitter,
    dependence = null_model$dependence,
    cor.shrinkage = null_model$cor_shrinkage,
    module.count = null_model$module_count,
    module.between.shrinkage = null_model$module_between_shrinkage,
    view.models = null_model$view_models
  )
  .mimesys_replicate_to_mae(
    pcl,
    simulation = list(
      purpose = "null", n = as.integer(n), p = p,
      crosstalk = resolved$value, crosstalk_calibration = resolved$calibration,
      null_model = null_model, seed = seed,
      provenance = list(
        package_version = as.character(utils::packageVersion("MimeSys")),
        created_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
        selected_crosstalk = resolved$value,
        calibration_grid = if (is.null(resolved$calibration)) NULL else resolved$calibration$provenance$crosstalk_grid,
        template = tmpl$provenance
      )
    )
  )
}

#' Simulate multi-omics data with a compact MAE-first interface
#'
#' This is the recommended user-facing simulation interface. It accepts a
#' `MultiAssayExperiment` or a fitted `mimesys_template`, uses practical
#' defaults, and returns simulated data as a `MultiAssayExperiment`. The
#' lower-level [gen_simmba()] remains available when every simulation control
#' must be supplied separately.
#'
#' @param template A `MultiAssayExperiment`, a `SummarizedExperiment`, a named
#'   list of matrices, or a fitted `mimesys_template`.
#' @param n Number of samples per replicate.
#' @param p Number of features per view. See [gen_simmba()].
#' @param crosstalk Cross-view dependence strength. The default `"calibrate"`
#'   selects a value by matching empirical and simulated cross-view correlation
#'   blocks. Supply a numeric value to bypass calibration; `0` removes template
#'   shared structure and `1` uses the fitted template level.
#' @param preset Built-in simulation setting. Use [mimesys_presets()] to list
#'   choices. Explicit `outcome`, `signal`, and `null_model` values override
#'   the preset.
#' @param outcome Either one of `"linear"`, `"friedman"`, `"friedman2"`, or
#'   `"shared_latent"`, or a named list. List options are `model`, `type`,
#'   `snr`, `friedman_share`, `shared_latent_function`,
#'   `shared_latent_n_factors`, `surv_hscale`, `cens_lower`, and `cens_upper`.
#' @param signal Named list controlling spike-ins. Options are `de_prob`,
#'   `down_probability`, `effect_location`, `effect_scale`, `alignment`, and
#'   `share`. By default, `alignment = 1`; supply `share` to enforce a desired
#'   relative contribution among views with nonzero spike-in signal.
#' @param null_model Named list controlling null-data generation. Options are
#'   `dependence`, `feature_select`, `feature_jitter`, `cor_shrinkage`,
#'   `module_count`, `module_between_shrinkage`, and `view_models`.
#' @param template_args Named list passed to [fit_template()] only when
#'   `template` is not already fitted.
#' @param calibration Named list controlling automatic crosstalk calibration.
#'   Options are `grid`, `nrep`, `nsample`, and `max_features_per_view`.
#' @param train_fraction Fraction of samples labelled `"train"` in `colData()`.
#' @param nrep Number of independent replicates. One replicate returns a
#'   `MultiAssayExperiment`; multiple replicates return a
#'   `mimesys_mae_simulation` list with one MAE per replicate.
#' @param seed Random seed.
#' @param return_type `"mae"` (default) returns MAE-oriented output. Use
#'   `"simulation"` to return the underlying [gen_simmba()] result.
#'
#' @return A `MultiAssayExperiment`, a `mimesys_mae_simulation` list, or the
#'   lower-level simulation list when `return_type = "simulation"`.
#' @export
simulate_mimesys <- function(template,
                             n,
                             p = NULL,
                             crosstalk = "calibrate",
                             preset = "default",
                             outcome = NULL,
                             signal = list(),
                             null_model = list(),
                             template_args = list(),
                             calibration = list(),
                             train_fraction = 0.7,
                             nrep = 1,
                             seed = 1234,
                             return_type = c("mae", "simulation")) {
  return_type <- match.arg(return_type)
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 2L) {
    stop("n must be a single integer of at least 2.", call. = FALSE)
  }
  if (!is.numeric(train_fraction) || length(train_fraction) != 1L ||
      !is.finite(train_fraction) || train_fraction <= 0 || train_fraction >= 1) {
    stop("train_fraction must be strictly between 0 and 1.", call. = FALSE)
  }
  if (!is.numeric(nrep) || length(nrep) != 1L || nrep < 1L) {
    stop("nrep must be a positive integer.", call. = FALSE)
  }
  if (!is.character(preset) || length(preset) != 1L) {
    stop("preset must be a single preset name.", call. = FALSE)
  }
  preset_settings <- .mimesys_preset_options(preset)

  outcome_defaults <- list(
    model = "linear", type = "continuous", snr = 1,
    friedman_share = NULL, shared_latent_function = "linear",
    shared_latent_n_factors = NULL, surv_hscale = 1,
    cens_lower = 1, cens_upper = 3
  )
  outcome_defaults <- utils::modifyList(outcome_defaults, preset_settings$outcome)
  if (is.character(outcome)) outcome <- list(model = outcome)
  outcome <- .mimesys_merge_options(outcome, outcome_defaults, "outcome")
  outcome$model <- .mimesys_outcome_mode(outcome$model)

  signal_defaults <- list(
    de_prob = 0.1, down_probability = 0.5, effect_location = 1,
    effect_scale = 0.4, alignment = 1, share = NULL
  )
  signal_defaults <- utils::modifyList(signal_defaults, preset_settings$signal)
  signal <- .mimesys_merge_options(signal, signal_defaults, "signal")

  null_defaults <- list(
    dependence = "factor", feature_select = "random", feature_jitter = 0.2,
    cor_shrinkage = 0.05, module_count = NULL,
    module_between_shrinkage = 0.75, view_models = NULL
  )
  null_defaults <- utils::modifyList(null_defaults, preset_settings$null_model)
  null_model <- .mimesys_merge_options(null_model, null_defaults, "null_model")
  if (!is.list(template_args)) stop("template_args must be a named list.", call. = FALSE)

  fitted_template <- if (inherits(template, "mimesys_template")) {
    template
  } else {
    do.call(fit_template, c(list(template = template), template_args))
  }
  resolved <- .mimesys_resolve_crosstalk(
    fitted_template, n = as.integer(n), p = p, crosstalk = crosstalk,
    calibration = calibration, null_model = null_model, seed = seed
  )
  sim <- gen_simmba(
    nsample = as.integer(n), template = fitted_template, p = p, crosstalk = resolved$value,
    feature.select = null_model$feature_select,
    feature.jitter = null_model$feature_jitter,
    dependence = null_model$dependence,
    cor.shrinkage = null_model$cor_shrinkage,
    module.count = null_model$module_count,
    module.between.shrinkage = null_model$module_between_shrinkage,
    view.models = null_model$view_models,
    snr = outcome$snr, p.train = train_fraction,
    de.prob = signal$de_prob, de.downProb = signal$down_probability,
    de.facLoc = signal$effect_location, de.facScale = signal$effect_scale,
    de.align = signal$alignment, signal.share = signal$share,
    friedman.share = outcome$friedman_share, ygen.mode = outcome$model,
    shared.latent.fun = outcome$shared_latent_function,
    shared.latent.n_factors = outcome$shared_latent_n_factors,
    outcome.type = outcome$type, surv.hscale = outcome$surv_hscale,
    cens.lower = outcome$cens_lower, cens.upper = outcome$cens_upper,
    nrep = as.integer(nrep), seed = seed
  )
  sim$crosstalk_calibration <- resolved$calibration
  sim$provenance <- list(
    package_version = as.character(utils::packageVersion("MimeSys")),
    created_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    purpose = "outcome",
    preset = preset,
    seed = seed,
    selected_crosstalk = resolved$value,
    calibration_grid = if (is.null(resolved$calibration)) NULL else resolved$calibration$provenance$crosstalk_grid,
    template = fitted_template$provenance
  )
  if (identical(return_type, "simulation")) return(sim)

  maes <- lapply(seq_len(sim$nrep), function(i) as_mimesys_mae(sim, rep = i))
  if (length(maes) == 1L) return(maes[[1]])
  structure(
    list(replicates = maes, settings = .mimesys_simulation_settings(sim)),
    class = "mimesys_mae_simulation"
  )
}

#' @export
print.mimesys_mae_simulation <- function(x, ...) {
  cat("MimeSys MAE simulation collection\n")
  cat("  Replicates: ", length(x$replicates), "\n", sep = "")
  if (!is.null(x$settings$provenance$selected_crosstalk)) {
    cat("  Crosstalk: ", x$settings$provenance$selected_crosstalk, "\n", sep = "")
  }
  if (!is.null(x$settings$ygen.mode)) {
    cat("  Outcome: ", x$settings$ygen.mode, "\n", sep = "")
  }
  invisible(x)
}
