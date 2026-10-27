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
#' @examples
#' data(stelzer_mae)
#' tmpl <- fit_template(stelzer_mae, shared.rank = 2)
#' sim <- gen_simmba(n = 40, template = tmpl, p = 10, nrep = 1, seed = 1)
#' as_mimesys_mae(sim, rep = 1, split = "train")
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
        stop("x must be a gen_simmba() result containing the requested replicate.",
            call. = FALSE)
    }

    train <- x$trainDat[[rep]]
    test <- x$testDat[[rep]]
    dat <- switch(split,
        train = train,
        test = test,
        all = .mimesys_bind_simulation_splits(train, test)
    )
    .mimesys_replicate_to_mae(dat, simulation = .mimesys_simulation_settings(x),
        split = split)
}

.mimesys_bind_simulation_splits <- function(train, test) {
    out <- train
    out$feature_table <- cbind(train$feature_table, test$feature_table)
    out$sample_metadata <- rbind(train$sample_metadata, test$sample_metadata)
    out$sample_metadata$split <- c(
        rep("train", nrow(train$sample_metadata)),
        rep("test", nrow(test$sample_metadata))
    )
    out$truth <- train$truth
    out$truth$scores <- rbind(train$truth$scores, test$truth$scores)
    out$truth$shared <- rbind(train$truth$shared, test$truth$shared)
    out
}

.mimesys_replicate_to_mae <- function(dat, simulation = list(), split = "all") {
    feature_table <- as.matrix(dat$feature_table)
    feature_metadata <- as.data.frame(dat$feature_metadata,
        stringsAsFactors = FALSE)
    sample_metadata <- as.data.frame(dat$sample_metadata,
        stringsAsFactors = FALSE)
    if (nrow(feature_metadata) != nrow(feature_table) ||
        nrow(sample_metadata) != ncol(feature_table)) {
        stop("Simulation metadata dimensions do not match the simulated feature table.",
            call. = FALSE)
    }
    if (is.null(feature_metadata$featureType)) {
        stop("Simulation feature metadata must contain featureType.",
            call. = FALSE)
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

    view_index <- split(seq_len(nrow(feature_metadata)),
        factor(feature_metadata$featureType, levels = unique(feature_metadata$featureType)))
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
        stop(name, " contains unsupported option(s): ", paste(unknown,
            collapse = ", "), ".", call. = FALSE)
    }
    utils::modifyList(defaults, values)
}

.mimesys_outcome_mode <- function(model) {
    lookup <- c(
        linear = "linear", lm = "linear", friedman = "friedman", friedman2 = "friedman2",
        shared_latent = "shared_latent", sharedlatent = "shared_latent",
        interaction = "interaction"
    )
    model <- tolower(as.character(model))
    if (length(model) != 1L || !model %in% names(lookup)) {
        stop("outcome$model must be one of: linear, friedman, friedman2, shared_latent, interaction.",
            call. = FALSE)
    }
    unname(lookup[[model]])
}

.mimesys_preset_options <- function(preset) {
    presets <- list(
        default = list(
            description = "Aligned linear outcome with factor-based null simulation.",
            outcome = list(model = "linear"), signal = list(align = 1), null_model = list()
        ),
        linear = list(
            description = "Aligned additive linear outcome.",
            outcome = list(model = "linear"), signal = list(align = 1), null_model = list()
        ),
        friedman = list(
            description = "Nonlinear Friedman outcome from spiked observed features.",
            outcome = list(model = "friedman"), signal = list(align = 1), null_model = list()
        ),
        friedman2 = list(
            description = "Additive linear plus nonlinear Friedman outcome.",
            outcome = list(model = "friedman2"), signal = list(align = 1), null_model = list()
        ),
        shared_latent = list(
            description = "Linear outcome from all selected shared latent factors.",
            outcome = list(model = "shared_latent",
                shared_latent_function = "linear"),
            signal = list(align = 1), null_model = list()
        ),
        shared_latent_friedman = list(
            description = "Friedman outcome from the first five selected shared latent factors.",
            outcome = list(model = "shared_latent",
                shared_latent_function = "friedman"),
            signal = list(align = 1), null_model = list()
        ),
        shared_latent_friedman2 = list(
            description = "Linear plus Friedman outcome from selected shared latent factors.",
            outcome = list(model = "shared_latent",
                shared_latent_function = "friedman2"),
            signal = list(align = 1), null_model = list()
        ),
        interaction = list(
            description = "Linear plus cross-view feature-interaction outcome.",
            outcome = list(model = "interaction"),
            signal = list(align = 1), null_model = list()
        )
    )
    if (!preset %in% names(presets)) {
        stop("Unknown preset: ", preset,
            ". Use mimesys_presets() to list available presets.", call. = FALSE)
    }
    presets[[preset]]
}

#' List built-in MimeSys simulation presets
#'
#' @return A data frame describing built-in presets for [simulate_mimesys()].
#' @examples
#' mimesys_presets()
#' @export
mimesys_presets <- function() {
    preset_names <- c(
        "default", "linear", "friedman", "friedman2", "shared_latent",
        "shared_latent_friedman", "shared_latent_friedman2", "interaction"
    )
    settings <- lapply(preset_names, .mimesys_preset_options)
    data.frame(
        preset = preset_names,
        outcome_model = vapply(settings, function(x) x$outcome$model,
            character(1)),
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
    if (length(tmpl$views) < 2L || identical(as.integer(tmpl$shared.rank), 0L)) {
        message("Skipping crosstalk calibration: the template has ",
            if (length(tmpl$views) < 2L) "a single view" else "no shared latent factors",
            ", so crosstalk has no effect. Using crosstalk = 0.")
        return(list(value = 0, calibration = NULL))
    }
    p_check <- .normalize_feature_counts(p, tmpl$views,
        vapply(tmpl$fits, function(f) nrow(f$Lambda), integer(1)))
    if (any(p_check < 2L)) {
        message("Skipping crosstalk calibration: it needs at least two features per view, ",
            "but p is below 2 for ", paste(tmpl$views[p_check < 2L], collapse = ", "),
            ". Using the template level, crosstalk = 1.")
        return(list(value = 1, calibration = NULL))
    }
    calibration_defaults <- list(
        grid = seq(0, 1, by = 0.25), nrep = 5, nsample = NULL,
        max_features_per_view = NULL
    )
    calibration <- .mimesys_merge_options(calibration, calibration_defaults,
        "calibration")
    fit <- calibrate_crosstalk(
        template = tmpl,
        crosstalk.grid = calibration$grid,
        n = calibration$nsample,
        nrep = calibration$nrep,
        p = p,
        feature.select = null_model$feature_select,
        feature.jitter = null_model$feature_jitter,
        max.features.per.view = calibration$max_features_per_view,
        view.models = null_model$view_models,
        ## Offset the seed so the calibration replicates do not share random
        ## numbers with the simulation that follows
        seed = if (is.null(seed)) NULL else seed + 1L
    )
    list(value = fit$best_crosstalk, calibration = fit)
}

## design = "null" of simulate_mimesys(): null feature data, no outcome
.simulate_null_mae <- function(template,
    n,
    p = NULL,
    crosstalk = "calibrate",
    null_model = list(),
    template_args = list(),
    calibration = list(),
    nrep = 1,
    seed = 1234,
    return_type = c("mae", "simulation")) {
    return_type <- match.arg(return_type)
    ## Seed everything (template fitting, calibration, simulation) and restore
    ## the caller's random-number state on exit.
    if (!is.null(seed)) withr::local_seed(seed)
    setup <- .mimesys_wrapper_setup(template, n = n, p = p, nrep = nrep,
        crosstalk = crosstalk, null_model = null_model,
        template_args = template_args, calibration = calibration, seed = seed)
    tmpl <- setup$template
    null_model <- setup$null_model
    resolved <- setup$crosstalk
    reps <- lapply(seq_len(as.integer(nrep)), function(i) {
        .sim_from_template_impl(
            tmpl = tmpl, n = as.integer(n), p = p, crosstalk = resolved$value,
            feature.select = null_model$feature_select,
            feature.jitter = null_model$feature_jitter,
            view.models = null_model$view_models,
            back.transform = isTRUE(null_model$back_transform)
        )
    })
    settings <- list(
        purpose = "null", n = as.integer(n), p = p, nrep = as.integer(nrep),
        crosstalk = resolved$value, crosstalk_calibration = resolved$calibration,
        null_model = null_model, seed = seed,
        provenance = .mimesys_wrapper_provenance(tmpl, resolved)
    )
    if (identical(return_type, "simulation")) {
        return(list(replicates = reps, settings = settings, template = tmpl))
    }
    maes <- lapply(reps, .mimesys_replicate_to_mae, simulation = settings)
    .mimesys_collect_maes(maes, settings)
}


#' Simulate multi-omics data: the main MimeSys interface
#'
#' One function for every MimeSys simulation. It accepts a
#' `MultiAssayExperiment`, a list of matrices or a fitted template, fits the
#' template if needed, calibrates crosstalk by default, and returns the
#' simulated data as a `MultiAssayExperiment`. `design` chooses what is
#' simulated:
#'
#' * `"outcome"` (default): multi-omics features with planted outcome signal
#'   (continuous, binary or survival), split into training and test samples
#'   (engine: [gen_simmba()]);
#' * `"null"`: feature data only, with no outcome or signal (engine:
#'   [sim_from_template()]);
#' * `"multistudy"`: several studies with between-study heterogeneity, pooled
#'   in one MAE with a `study` column (engine: [gen_multistudy()]);
#' * `"mediation"`: an exposure, mediators in one or more views and an outcome,
#'   with the exact natural indirect, direct and total effects, including
#'   sequential (multi-stage) mediation (engine: [gen_mediation()]).
#'
#' Each design reads its own option lists (see Details); supplying a list that
#' a design does not use is an error, so a setting is never silently ignored.
#' The lower-level engines remain exported for full control.
#'
#' @param template A `MultiAssayExperiment`, a `SummarizedExperiment`, a named
#'   list of matrices, or a fitted `mimesys_template`. With
#'   `design = "multistudy"` it may also be one template (or data set) per
#'   study: a list of templates from [fit_study_templates()] or a list of data
#'   sets (fitted with `fit_study_templates(template, ...)`, using
#'   `template_args`); see [gen_multistudy()].
#' @param n Number of samples per replicate. With `design = "multistudy"`, the
#'   samples per study: one value for every study (with `nstudy`) or one value
#'   per study.
#' @param p Number of features per view. See [gen_simmba()].
#' @param crosstalk Cross-view dependence strength. The default `"calibrate"`
#'   selects a value by matching empirical and simulated cross-view correlation
#'   blocks. Supply a numeric value to bypass calibration; `0` removes template
#'   shared structure and `1` uses the fitted template level. With one
#'   template per study, calibration is done for each study on its own
#'   template and size.
#' @param design What to simulate: `"outcome"` (default), `"null"`,
#'   `"multistudy"` or `"mediation"`.
#' @param preset Built-in outcome setting (designs `"outcome"` and
#'   `"multistudy"`). Use [mimesys_presets()] to list choices. Explicit
#'   `outcome`, `signal`, and `null_model` values override the preset.
#' @param outcome Outcome model (designs `"outcome"` and `"multistudy"`):
#'   either one of `"linear"`, `"friedman"`, `"friedman2"`, `"shared_latent"`,
#'   or `"interaction"`, or a named list. List options are
#'   `model`, `type`, `snr`, `friedman_share`, `shared_latent_function`,
#'   `shared_latent_n_factors`, `shared_latent_weight` (share of signal from
#'   the shared factors; the rest is view-specific), `nonlinear_weight`,
#'   `interaction_weight`, `n_interactions`, `prevalence` (binary target case
#'   proportion), `baseline_hazard`, `cens_rate` (survival target censoring
#'   proportion; with `"multistudy"` one value or one per study), `cens_lower`,
#'   `cens_upper`, `cens_type` (`"uniform"`, `"exponential"` or
#'   `"administrative"`), `cens_times` (user censoring times), `cens_dependence`
#'   (informative censoring), `cens_accrual`, `cens_dropout`, `baseline_shape`
#'   (Weibull shape) and `keep_true_times`. See [gen_simmba()].
#' @param signal Named list controlling spike-ins (designs `"outcome"` and
#'   `"multistudy"`). Options are `prob`, `neg_prob`, `effect_meanlog`,
#'   `effect_sdlog`, `align`, and `share` (see `signal.prob`,
#'   `signal.neg.prob`, `effect.meanlog`, `effect.sdlog`, `signal.align` and
#'   `signal.share` in [gen_simmba()]). By default, `align = 1`; supply `share`
#'   to enforce a desired relative contribution among views with nonzero
#'   spike-in signal.
#' @param heterogeneity Named list of between-study differences (design
#'   `"multistudy"`). Options are `effect_cor`, `support_overlap`,
#'   `intercept_sd`, `snr_sd`, `shift_sd`, `batch_loc_sd`, `batch_scale_sd`
#'   and `holdout_studies` (see `effect.cor`, `support.overlap`,
#'   `intercept.sd`, `snr.sd`, `shift.sd`, `batch.loc.sd`, `batch.scale.sd`
#'   and `holdout.studies` in [gen_multistudy()]). The defaults give
#'   homogeneous studies.
#' @param nstudy Number of studies (design `"multistudy"`). Defaults to
#'   `length(n)`, or to the number of study templates.
#' @param mediation Named list of mediation settings (design `"mediation"`):
#'   `n_mediators`, `tie`, `de`, `dir_a`, `dir_b`, `effect_type`,
#'   `exposure_type`, `exposure_prevalence`, `n_interactions`, `sigma`,
#'   `outcome_type` (`"continuous"` or `"binary"`), `outcome_prevalence`,
#'   `min_variability`, `n_mc`, and for sequential mediation `stages`,
#'   `stage_links`, `stage_direct` and `stage_share`. See [gen_mediation()] for
#'   their meaning and defaults.
#' @param null_model Named list controlling null-data generation (every
#'   design). Options are `feature_select`, `feature_jitter`, `view_models`,
#'   and `back_transform` (default `FALSE`: views fitted with a
#'   `view.transform` are returned on the fitted, transformed scale, which is
#'   the scale outcomes, mediation paths and true effects are defined on; set
#'   `TRUE` for the original measurement scale, see [sim_from_template()]).
#' @param template_args Named list passed to [fit_template()] (or
#'   [fit_study_templates()]) only when `template` is not already fitted.
#' @param calibration Named list controlling automatic crosstalk calibration.
#'   Options are `grid`, `nrep`, `nsample`, and `max_features_per_view`. With
#'   one template for several studies, calibration uses the average study
#'   size.
#' @param train_fraction Fraction of samples labelled `"train"` in `colData()`
#'   (designs `"outcome"` and `"multistudy"`; with `"multistudy"`, within each
#'   non-holdout study).
#' @param nrep Number of independent replicates. One replicate returns a
#'   `MultiAssayExperiment`; multiple replicates return a
#'   `mimesys_mae_simulation` list with one MAE per replicate.
#' @param seed Random seed. The caller's random-number state is restored on
#'   exit.
#' @param return_type `"mae"` (default) returns MAE output. `"simulation"`
#'   returns the underlying engine result: the [gen_simmba()],
#'   [gen_multistudy()] or [gen_mediation()] list, or for `design = "null"` a
#'   list with the pcl-style `replicates` from [sim_from_template()],
#'   `settings` and the fitted `template`.
#'
#' @details
#' Options used by each design (`template`, `n`, `p`, `crosstalk`,
#' `null_model`, `template_args`, `calibration`, `nrep`, `seed` and
#' `return_type` apply to all):
#'
#' | `design` | Option lists |
#' |---|---|
#' | `"outcome"` | `preset`, `outcome`, `signal`, `train_fraction` |
#' | `"null"` | none |
#' | `"multistudy"` | `preset`, `outcome`, `signal`, `heterogeneity`, `nstudy`, `train_fraction` |
#' | `"mediation"` | `mediation` |
#'
#' Results are those of the engine with the same settings and seed. The
#' selected crosstalk and the calibration result are stored in
#' `metadata(x)$simulation`; the simulation truth (coefficients, study
#' parameters, mediation effects) in `rowData()`, `colData()` and
#' `metadata(x)$truth`.
#'
#' @return A `MultiAssayExperiment`, a `mimesys_mae_simulation` list (several
#'   replicates), or the engine result when `return_type = "simulation"`.
#' @examples
#' data(stelzer_mae)
#' tmpl <- fit_template(stelzer_mae, shared.rank = 2, seed = 1)
#'
#' ## Outcome (default design)
#' sim_mae <- simulate_mimesys(tmpl, n = 60, p = 10, crosstalk = 0.5, seed = 1)
#' head(MultiAssayExperiment::colData(sim_mae))
#'
#' ## Survival outcome with administrative censoring
#' surv_mae <- simulate_mimesys(tmpl, n = 60, p = 10, crosstalk = 0.5,
#'     outcome = list(type = "survival", cens_rate = 0.3,
#'         cens_type = "administrative"), seed = 1)
#'
#' ## Null data
#' null_mae <- simulate_mimesys(tmpl, n = 40, p = 10, design = "null",
#'     calibration = list(grid = c(0, 0.5, 1), nrep = 1))
#' S4Vectors::metadata(null_mae)$simulation$crosstalk
#'
#' ## Several studies
#' ms_mae <- simulate_mimesys(tmpl, n = c(40, 50, 30), p = 10, crosstalk = 1,
#'     design = "multistudy",
#'     heterogeneity = list(effect_cor = 0.5, shift_sd = 0.3,
#'         holdout_studies = 3), seed = 1)
#' table(MultiAssayExperiment::colData(ms_mae)$study,
#'     MultiAssayExperiment::colData(ms_mae)$split)
#'
#' ## Mediation
#' med_mae <- simulate_mimesys(tmpl, n = 150, p = 20, crosstalk = 1,
#'     design = "mediation",
#'     mediation = list(n_mediators = c(CYTOF = 2, Metabolomics = 2, Proteomics = 0),
#'         tie = 0.5, n_mc = 2000))
#' S4Vectors::metadata(med_mae)$truth$mediation$NIE
#' @export
simulate_mimesys <- function(template,
    n,
    p = NULL,
    crosstalk = "calibrate",
    design = c("outcome", "null", "multistudy", "mediation"),
    preset = "default",
    outcome = NULL,
    signal = list(),
    heterogeneity = list(),
    nstudy = NULL,
    mediation = list(),
    null_model = list(),
    template_args = list(),
    calibration = list(),
    train_fraction = 0.7,
    nrep = 1,
    seed = 1234,
    return_type = c("mae", "simulation")) {
    design <- match.arg(design)
    return_type <- match.arg(return_type)
    .check_design_arguments(design, names(match.call())[-1L], n, template)
    switch(design,
        outcome = .simulate_outcome_mae(template = template, n = n, p = p,
            crosstalk = crosstalk, preset = preset, outcome = outcome,
            signal = signal, null_model = null_model,
            template_args = template_args, calibration = calibration,
            train_fraction = train_fraction, nrep = nrep, seed = seed,
            return_type = return_type),
        null = .simulate_null_mae(template = template, n = n, p = p,
            crosstalk = crosstalk, null_model = null_model,
            template_args = template_args, calibration = calibration,
            nrep = nrep, seed = seed, return_type = return_type),
        multistudy = .simulate_multistudy_mae(template = template, n = n,
            nstudy = nstudy, p = p, crosstalk = crosstalk, preset = preset,
            outcome = outcome, signal = signal, heterogeneity = heterogeneity,
            null_model = null_model, template_args = template_args,
            calibration = calibration, train_fraction = train_fraction,
            nrep = nrep, seed = seed, return_type = return_type),
        mediation = .simulate_mediation_mae(template = template, n = n, p = p,
            crosstalk = crosstalk, mediation = mediation,
            null_model = null_model, template_args = template_args,
            calibration = calibration, nrep = nrep, seed = seed,
            return_type = return_type))
}

## Option lists each design reads; a supplied option that the design does not
## use is an error rather than silently ignored.
.design_options <- list(
    outcome = c("preset", "outcome", "signal", "train_fraction"),
    null = character(0),
    multistudy = c("preset", "outcome", "signal", "heterogeneity", "nstudy",
        "train_fraction"),
    mediation = "mediation")

.check_design_arguments <- function(design, supplied, n, template) {
    specific <- unique(unlist(.design_options))
    unused <- setdiff(intersect(supplied, specific), .design_options[[design]])
    if (length(unused) > 0L) {
        users <- vapply(unused, function(a) {
            d <- names(.design_options)[vapply(.design_options,
                function(o) a %in% o, logical(1))]
            paste0("`", a, "` (design ", paste0("\"", d, "\"", collapse = " or "), ")")
        }, character(1))
        stop("design = \"", design, "\" does not use ", paste(users, collapse = ", "),
            ".", call. = FALSE)
    }
    if (!identical(design, "multistudy")) {
        if (is.numeric(n) && length(n) > 1L) {
            stop("n has ", length(n), " values; use design = \"multistudy\" for ",
                "several studies.", call. = FALSE)
        }
        if (inherits(template, "mimesys_template_list")) {
            stop("A list of study templates needs design = \"multistudy\".", call. = FALSE)
        }
    }
    invisible(TRUE)
}

## design = "outcome" of simulate_mimesys()
.simulate_outcome_mae <- function(template,
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
    .check_count(n, "n", min = 2L)
    .check_count(nrep, "nrep", min = 1L)
    .check_prevalence(train_fraction, "train_fraction")
    if (!is.null(seed)) withr::local_seed(seed)
    opts <- .mimesys_outcome_signal_options(preset, outcome, signal)
    preset_settings <- opts$preset
    outcome <- opts$outcome
    signal <- opts$signal

    setup <- .mimesys_wrapper_setup(template, n = n, p = p, nrep = nrep,
        crosstalk = crosstalk, null_model = null_model,
        template_args = template_args, calibration = calibration, seed = seed,
        null_preset = preset_settings$null_model)
    fitted_template <- setup$template
    null_model <- setup$null_model
    resolved <- setup$crosstalk
    sim <- gen_simmba(
        n = as.integer(n), template = fitted_template, p = p, crosstalk = resolved$value,
        feature.select = null_model$feature_select,
        feature.jitter = null_model$feature_jitter,
        view.models = null_model$view_models,
        back.transform = isTRUE(null_model$back_transform),
        snr = outcome$snr, train.fraction = train_fraction,
        signal.prob = signal$prob, signal.neg.prob = signal$neg_prob,
        effect.meanlog = signal$effect_meanlog, effect.sdlog = signal$effect_sdlog,
        signal.align = signal$align, signal.share = signal$share,
        friedman.share = outcome$friedman_share, outcome.model = outcome$model,
        shared.latent.function = outcome$shared_latent_function,
        shared.latent.n.factors = outcome$shared_latent_n_factors,
        shared.latent.weight = outcome$shared_latent_weight,
        nonlinear.weight = outcome$nonlinear_weight,
        interaction.weight = outcome$interaction_weight,
        n.interactions = outcome$n_interactions,
        outcome.type = outcome$type, binary.prevalence = outcome$prevalence,
        baseline.hazard = outcome$baseline_hazard, cens.rate = outcome$cens_rate,
        cens.lower = outcome$cens_lower, cens.upper = outcome$cens_upper,
        cens.type = outcome$cens_type, cens.times = outcome$cens_times,
        cens.dependence = outcome$cens_dependence,
        cens.accrual = outcome$cens_accrual, cens.dropout = outcome$cens_dropout,
        baseline.shape = outcome$baseline_shape,
        keep.true.times = isTRUE(outcome$keep_true_times),
        nrep = as.integer(nrep), seed = seed
    )
    sim$crosstalk_calibration <- resolved$calibration
    sim$provenance <- .mimesys_wrapper_provenance(fitted_template, resolved,
        purpose = "outcome", preset = preset, seed = seed)
    if (identical(return_type, "simulation")) return(sim)

    maes <- lapply(seq_len(sim$nrep), function(i) as_mimesys_mae(sim, rep = i))
    .mimesys_collect_maes(maes, .mimesys_simulation_settings(sim))
}

## Preset, outcome and signal options of the outcome wrappers
## (designs "outcome" and "multistudy" of simulate_mimesys())
.mimesys_outcome_signal_options <- function(preset, outcome, signal) {
    if (!is.character(preset) || length(preset) != 1L) {
        stop("preset must be a single preset name.", call. = FALSE)
    }
    preset_settings <- .mimesys_preset_options(preset)

    outcome_defaults <- list(
        model = "linear", type = "continuous", snr = 1,
        friedman_share = NULL, shared_latent_function = "linear",
        shared_latent_n_factors = NULL, nonlinear_weight = 0.5,
        interaction_weight = 0.5, n_interactions = NULL,
        prevalence = 0.5, baseline_hazard = 1,
        cens_rate = NULL, cens_lower = 1, cens_upper = 3,
        shared_latent_weight = 1, cens_type = "uniform", cens_times = NULL,
        cens_dependence = 0, cens_accrual = 0.5, cens_dropout = 0,
        baseline_shape = 1, keep_true_times = FALSE
    )
    outcome_defaults <- utils::modifyList(outcome_defaults,
        preset_settings$outcome)
    if (is.character(outcome)) outcome <- list(model = outcome)
    outcome <- .mimesys_merge_options(outcome, outcome_defaults, "outcome")
    outcome$model <- .mimesys_outcome_mode(outcome$model)

    signal_defaults <- list(
        prob = 0.1, neg_prob = 0.5, effect_meanlog = 1,
        effect_sdlog = 0.4, align = 1, share = NULL
    )
    signal_defaults <- utils::modifyList(signal_defaults,
        preset_settings$signal)
    signal <- .mimesys_merge_options(signal, signal_defaults, "signal")
    list(preset = preset_settings, outcome = outcome, signal = signal)
}

#' @export
print.mimesys_mae_simulation <- function(x, ...) {
    cat("MimeSys MAE simulation collection\n")
    cat("  Replicates: ", length(x$replicates), "\n", sep = "")
    if (!is.null(x$settings$provenance$selected_crosstalk)) {
        cat("  Crosstalk: ", x$settings$provenance$selected_crosstalk, "\n",
            sep = "")
    }
    if (!is.null(x$settings$outcome.model)) {
        cat("  Outcome: ", x$settings$outcome.model, "\n", sep = "")
    }
    invisible(x)
}
