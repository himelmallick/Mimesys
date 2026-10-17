#' Validate a MimeSys input container
#'
#' Checks a `MultiAssayExperiment`, `SummarizedExperiment`, or named list for
#' the minimum structural requirements of MimeSys and returns a per-view input
#' summary. When sample identifiers are available, the summary includes the
#' number shared across all views.
#'
#' @param x Input accepted by [fit_template()].
#' @param assay Assay name or index passed to [as_mimesys_view_list()].
#' @param min_common_samples Minimum number of samples required across views.
#'
#' @return A data frame with one row per view and a `common_samples` attribute.
#' @examples
#' data(stelzer_mae)
#' validate_mimesys_input(stelzer_mae)
#' @export
validate_mimesys_input <- function(x, assay = NULL, min_common_samples = 10) {
    if (!is.numeric(min_common_samples) || length(min_common_samples) != 1L || min_common_samples < 2L) {
        stop("min_common_samples must be a single number of at least 2.",
            call. = FALSE)
    }
    views <- as_mimesys_view_list(x, assay = assay)
    common <- Reduce(intersect, lapply(views, rownames))
    if (length(views) > 1L && length(common) < min_common_samples) {
        stop(
            "Only ", length(common), " samples are shared across views; at least ",
            min_common_samples, " are required.", call. = FALSE
        )
    }
    out <- data.frame(
        view = names(views),
        samples = vapply(views, nrow, integer(1)),
        features = vapply(views, ncol, integer(1)),
        missing_values = vapply(views, function(x) sum(!is.finite(x)),
            integer(1)),
        constant_features = vapply(views, function(x) sum(apply(x, 2,
            function(v) length(unique(v[is.finite(v)])) <= 1L)), integer(1)),
        stringsAsFactors = FALSE,
        row.names = NULL
    )
    attr(out, "common_samples") <- length(common)
    out
}

.mimesys_marginal_diagnostics <- function(observed, simulated,
    max_features_per_view = NULL) {
    aligned <- .mimesys_align_observed_simulated(
        observed, simulated, max_features_per_view = max_features_per_view
    )
    probs <- seq(0.05, 0.95, length.out = 19)
    do.call(rbind, lapply(names(aligned$observed), function(view) {
        obs <- aligned$observed[[view]]
        sim <- aligned$simulated[[view]]
        stats_j <- vapply(seq_len(ncol(obs)), function(j) {
            o <- obs[is.finite(obs[, j]), j]
            s <- sim[is.finite(sim[, j]), j]
            qo <- stats::quantile(o, probs = probs, names = FALSE)
            qs <- stats::quantile(s, probs = probs, names = FALSE)
            sd_o <- stats::sd(o)
            ## Kolmogorov-Smirnov distance between empirical CDFs (scale-free)
            grid <- sort(unique(c(o, s)))
            ks <- max(abs(stats::ecdf(o)(grid) - stats::ecdf(s)(grid)))
            c(
                quantile_rmse = sqrt(mean((qo - qs)^2)),
                std_quantile_rmse = if (is.finite(sd_o) && sd_o > 0) sqrt(mean((qo - qs)^2)) / sd_o else NA_real_,
                ks = ks
            )
            }, numeric(3))
        data.frame(
            view = view,
            features = ncol(obs),
            mean_ks = mean(stats_j["ks", ], na.rm = TRUE),
            median_ks = stats::median(stats_j["ks", ], na.rm = TRUE),
            mean_std_quantile_rmse = mean(stats_j["std_quantile_rmse", ],
                na.rm = TRUE),
            mean_quantile_rmse = mean(stats_j["quantile_rmse", ], na.rm = TRUE),
            median_quantile_rmse = stats::median(stats_j["quantile_rmse", ],
                na.rm = TRUE),
            stringsAsFactors = FALSE
        )
    }))
}

#' Validate a fitted template against its empirical data
#'
#' Calibrates crosstalk using null simulations and returns diagnostics for
#' cross-view reconstruction, correlation-pattern preservation, and marginal
#' distribution preservation. The report is designed to be inspected directly
#' or plotted with [plot_template_validation()].
#'
#' @param template A fitted `mimesys_template` or input accepted by
#'   [fit_template()].
#' @param p Feature counts for validation simulation. See [sim_from_template()].
#' @param nsample Number of samples per null validation replicate. Defaults to
#'   the template sample size.
#' @param nrep Number of null replicates used in calibration.
#' @param crosstalk_grid Candidate crosstalk values.
#' @param max_features_per_view Optional cap for detailed correlation and
#'   marginal diagnostics: the most variable observed features of each view are
#'   kept, and the same features are used for the simulated data. This does
#'   not alter the calibration, which uses all matched features by default.
#' @param view.models Optional named vector of per-view observation models
#'   passed to [sim_from_template()].
#' @param template_args Named list passed to [fit_template()] for raw input.
#' @param seed Random seed.
#'
#' @return An object of class `mimesys_template_validation`. Its
#'   `marginal_diagnostics` table reports, per view, scale-free summaries
#'   (`mean_ks`, `median_ks`: Kolmogorov-Smirnov distance between observed and
#'   simulated feature distributions; `mean_std_quantile_rmse`: quantile RMSE
#'   divided by the observed feature SD) alongside the raw-unit quantile RMSE.
#' @examples
#' data(stelzer_mae)
#' tmpl <- fit_template(stelzer_mae, shared.rank = 2)
#' report <- validate_template(
#'     tmpl, nrep = 1, crosstalk_grid = c(0, 0.5, 1),
#'     max_features_per_view = 20
#' )
#' report
#' report$marginal_diagnostics
#' @export
validate_template <- function(template,
    p = NULL,
    nsample = NULL,
    nrep = 5,
    crosstalk_grid = seq(0, 1, by = 0.25),
    max_features_per_view = 100,
    view.models = NULL,
    template_args = list(),
    seed = 1234) {
    if (!is.null(seed)) withr::local_seed(seed)
    if (!is.list(template_args)) stop("template_args must be a named list.",
        call. = FALSE)
    tmpl <- if (inherits(template, "mimesys_template")) template else {
        do.call(fit_template, c(list(template = template), template_args))
    }
    if (is.null(tmpl$observed)) {
        stop("Template must retain observed matrices; refit with fit_template().",
            call. = FALSE)
    }
    if (is.null(nsample)) nsample <- tmpl$n0

    calibrate <- length(tmpl$views) >= 2L && !identical(as.integer(tmpl$shared.rank), 0L)
    if (calibrate) {
        calibration <- calibrate_crosstalk(
            tmpl,
            crosstalk_grid = crosstalk_grid,
            nsample = nsample,
            nrep = nrep,
            p = p,
            view.models = view.models,
            seed = seed
        )
        selected <- calibration$best_crosstalk
    } else {
        message("Skipping crosstalk calibration: the template has ",
            if (length(tmpl$views) < 2L) "a single view" else "no shared latent factors",
            ", so crosstalk has no effect.")
        calibration <- NULL
        selected <- 0
    }
    if (!is.null(seed)) withr::local_seed(seed + 1L)
    pcl <- .sim_from_template_impl(
        tmpl, n = nsample, p = p, crosstalk = selected,
        view.models = view.models
    )
    simulated <- .mimesys_pcl_view_matrices(pcl,
        feature_names = "template_feature")
    observed <- tmpl$observed
    ## Choose the top-variable features on the observed data and use the same
    ## features in the simulated data, so the two correlation structures are
    ## compared feature for feature.
    aligned <- .mimesys_align_observed_simulated(observed, simulated,
        max_features_per_view = max_features_per_view)
    obs_corr <- blockwise_feature_correlation(aligned$observed)
    sim_corr <- blockwise_feature_correlation(aligned$simulated)

    out <- list(
        template = tmpl,
        calibration = calibration,
        selected_crosstalk = selected,
        reconstruction_error = feature_correlation_frobenius(observed,
            simulated),
        correlation_structure = compare_feature_correlation_structure(obs_corr,
            sim_corr),
        marginal_diagnostics = .mimesys_marginal_diagnostics(
            observed, simulated, max_features_per_view = max_features_per_view
        ),
        settings = list(
            p = p, nsample = nsample, nrep = nrep,
            max_features_per_view = max_features_per_view, view.models = view.models, seed = seed,
            package_version = as.character(utils::packageVersion("MimeSys")),
            created_at = format(Sys.time(), tz = "UTC", usetz = TRUE)
        )
    )
    class(out) <- "mimesys_template_validation"
    out
}

#' @export
print.mimesys_template_validation <- function(x, ...) {
    cat("MimeSys template validation\n")
    cat("  Selected crosstalk: ", x$selected_crosstalk, "\n", sep = "")
    if (is.null(x$calibration)) {
        cat("  Calibration: skipped (single view or no shared factors)\n")
    } else {
        cat("  Near-optimal values: ", paste(signif(x$calibration$near_optimal$crosstalk, 3),
            collapse = ", "), "\n", sep = "")
    }
    cross <- x$reconstruction_error[!x$reconstruction_error$diagonal_block, ,
        drop = FALSE]
    if (nrow(cross) > 0L) {
        cat("  Mean cross-block RMSE: ", round(mean(cross$rescaled_frobenius), 4),
            "\n", sep = "")
    }
    invisible(x)
}

#' Plot a MimeSys template-validation report
#'
#' @param x A `mimesys_template_validation` object.
#' @param type Plot type: `"calibration"`, `"correlation_structure"`, or
#'   `"marginal"`.
#'
#' @return A ggplot object.
#' @examples
#' data(stelzer_mae)
#' tmpl <- fit_template(stelzer_mae, shared.rank = 2)
#' report <- validate_template(
#'     tmpl, nrep = 1, crosstalk_grid = c(0, 0.5, 1),
#'     max_features_per_view = 20
#' )
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'     plot_template_validation(report, type = "marginal")
#' }
#' @export
plot_template_validation <- function(x,
    type = c("calibration", "correlation_structure", "marginal")) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("plot_template_validation() requires ggplot2.", call. = FALSE)
    }
    if (!inherits(x, "mimesys_template_validation")) {
        stop("x must be returned by validate_template().", call. = FALSE)
    }
    type <- match.arg(type)
    if (identical(type, "calibration")) {
        if (is.null(x$calibration)) {
            stop("No calibration to plot: calibration was skipped for this template.",
                call. = FALSE)
        }
        return(plot_crosstalk_calibration(x$calibration))
    }
    if (identical(type, "correlation_structure")) {
        structure_df <- x$correlation_structure
        structure_df$block <- paste(
            structure_df$view_x, structure_df$view_y, sep = " vs "
        )
        return(
            ggplot2::ggplot(structure_df, ggplot2::aes(x = block,
                y = pearson)) +
                ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                    color = "grey60") +
                ggplot2::geom_col(fill = "#2C7FB8") +
                ggplot2::coord_cartesian(ylim = c(-1, 1)) +
                ggplot2::labs(x = NULL, y = "Correlation-pattern preservation",
                    title = "Observed versus simulated correlation patterns") +
                ggplot2::theme_bw(base_size = 12) +
                ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30,
                    hjust = 1))
        )
    }
    ggplot2::ggplot(x$marginal_diagnostics, ggplot2::aes(x = view,
        y = mean_ks)) +
        ggplot2::geom_col(fill = "#41AB5D") +
        ggplot2::labs(x = NULL, y = "Mean feature KS distance",
            title = "Marginal distribution reconstruction") +
        ggplot2::theme_bw(base_size = 12)
}
