#' Calibrate crosstalk to empirical cross-view correlation structure
#'
#' Sweeps candidate crosstalk values, simulates null data from the supplied
#' template, and chooses the value that minimizes the rescaled Frobenius
#' distance
#' between empirical and simulated cross-view feature-correlation blocks.
#'
#' @param template A fitted \code{mimesys_template} or any object accepted by
#'   \code{fit_template()}.
#' @param crosstalk.grid Numeric vector of crosstalk values to evaluate.
#' @param n Number of simulated samples per replicate.
#' @param nrep Number of replicates per crosstalk value.
#' @param p Optional named integer vector of simulated features per view.
#' @param feature.select,feature.jitter Feature-generation settings passed to
#'   \code{sim_from_template()}.
#' @param max.features.per.view Optional named cap used when computing
#'   correlation distances. If `NULL`, every matched empirical feature is used.
#' @param view.models Optional named vector of per-view observation models
#'   passed
#'   to \code{sim_from_template()}.
#' @param seed Random seed.
#' @param ... Passed to \code{fit_template()} when \code{template} is raw data.
#'
#' @details Simulated data are always compared with the matrices stored in the
#'   fitted template (after any \code{view.transform}), so raw input fitted with
#'   a transform is calibrated on the transformed scale. Grid values above
#'   \code{crosstalk_max(template)} are dropped with a warning, and a warning is
#'   given when the best value sits at the upper end of the grid.
#'
#' @return A list. `best_crosstalk` is the grid value with the smallest mean
#'   cross-view distance, where each replicate's distance is averaged over all
#'   cross-view blocks. `cross_view_summary` gives that mean, its SD and SE per
#'   grid value; `near_optimal` lists values within one SE of the minimum;
#'   `selection_frequency` gives how often each value is best in individual
#'   replicates. `cross_block_summary` and `summary` give per-block distances
#'   and `raw` the replicate-level distances.
#' @examples
#' data(stelzer_mae)
#' tmpl <- fit_template(stelzer_mae, shared.rank = 2)
#' cal <- calibrate_crosstalk(tmpl, crosstalk.grid = c(0, 0.5, 1), nrep = 2)
#' cal$best_crosstalk
#' cal$selection_frequency
#' @export
calibrate_crosstalk <- function(
    template,
    crosstalk.grid = seq(0, 1, by = 0.25),
    n = NULL,
    nrep = 5,
    p = NULL,
    feature.select = c("random", "variance"),
    feature.jitter = 0.2,
    max.features.per.view = NULL,
    view.models = NULL,
    seed = 1234,
    ...
) {
    feature.select <- match.arg(feature.select)
    if (!is.null(seed)) withr::local_seed(seed)
    tmpl <- .resolve_template(template, list(...))
    if (is.null(tmpl$observed)) {
        stop(
            "This fitted template does not retain its empirical matrices. Refit it with fit_template() before calibrating crosstalk.",
            call. = FALSE
        )
    }
    ## Always compare against the matrices the template was fitted on, i.e. after
    ## any view.transform and feature filtering, so observed and simulated data
    ## are on the same scale.
    observed_matrices <- tmpl$observed
    if (length(tmpl$views) < 2L) {
        stop("calibrate_crosstalk() needs a template with at least two views.",
            call. = FALSE)
    }
    if (is.null(n)) {
        n <- tmpl$n0
    }
    if (!is.numeric(crosstalk.grid) || length(crosstalk.grid) == 0L ||
        any(!is.finite(crosstalk.grid)) || any(crosstalk.grid < 0)) {
        stop("crosstalk.grid must be a vector of finite, non-negative numbers.",
            call. = FALSE)
    }
    crosstalk.grid <- sort(unique(crosstalk.grid))
    ct_limit <- crosstalk_max(tmpl)
    if (identical(as.integer(tmpl$shared.rank), 0L)) {
        warning("The template has no shared latent factors, so crosstalk has no effect ",
            "on the simulated data; calibration cannot distinguish grid values.",
            call. = FALSE)
    }
    too_big <- crosstalk.grid > ct_limit + 1e-10
    if (any(too_big)) {
        warning(
            "Dropping crosstalk.grid value(s) above crosstalk_max(template) = ",
            sprintf("%.3f", ct_limit), ": ", paste(crosstalk.grid[too_big],
                collapse = ", "),
            call. = FALSE
        )
        crosstalk.grid <- crosstalk.grid[!too_big]
        if (length(crosstalk.grid) == 0L) {
            stop("No crosstalk.grid values are within the valid range.",
                call. = FALSE)
        }
    }

    ## Feature counts: each view needs at least two features for correlation
    ## blocks.
    p_check <- .normalize_feature_counts(p, tmpl$views,
        vapply(tmpl$fits, function(f) nrow(f$Lambda), integer(1)))
    if (any(p_check < 2L)) {
        stop("calibrate_crosstalk() needs at least two features per view to compare ",
            "correlation blocks, but p is below 2 for: ",
            paste(tmpl$views[p_check < 2L], collapse = ", "),
            ". Use a larger p or supply crosstalk as a number.", call. = FALSE)
    }
    .check_count(nrep, "nrep", min = 1L)
    .check_count(n, "n", min = 3L)
    obs_cache <- .mimesys_observed_cor_cache(observed_matrices)

    rows <- list()
    idx <- 1L
    for (ct in crosstalk.grid) {
        pcl <- .sim_from_template_impl(
            tmpl = tmpl,
            n = n,
            p = p,
            crosstalk = ct,
            feature.select = feature.select,
            feature.jitter = feature.jitter,
            view.models = view.models
        )
        for (rep_idx in seq_len(nrep)) {
            ## Use a fresh null replicate for every distance calculation. No outcome
            ## generation or spike-in parameters enter crosstalk calibration.
            if (rep_idx > 1L) {
                pcl <- .sim_from_template_impl(
                    tmpl = tmpl, n = n, p = p, crosstalk = ct,
                    feature.select = feature.select, feature.jitter = feature.jitter,
                    view.models = view.models
                )
            }
            sim_mats <- .mimesys_pcl_view_matrices(pcl,
                feature_names = "template_feature")
            aligned <- .mimesys_align_observed_simulated(observed_matrices,
                sim_mats, max_features_per_view = max.features.per.view)
            distances <- .mimesys_frobenius_blocks(aligned$observed,
                aligned$simulated, observed_cache = obs_cache)
            distances$crosstalk <- ct
            distances$rep <- rep_idx
            rows[[idx]] <- distances
            idx <- idx + 1L
        }
    }

    raw <- do.call(rbind, rows)
    summary <- .mimesys_summarize_frobenius(raw)
    cross <- summary[!summary$diagonal_block, , drop = FALSE]
    ## Selection criterion: for each replicate, the mean distance over all
    ## cross-view blocks; then average over replicates for each crosstalk value.
    cross_by_replicate <- stats::aggregate(
        rescaled_frobenius ~ crosstalk + rep,
        data = raw[!raw$diagonal_block, , drop = FALSE],
        FUN = mean
    )
    overall <- stats::aggregate(
        rescaled_frobenius ~ crosstalk,
        data = cross_by_replicate,
        FUN = function(v) c(mean = mean(v), sd = if (length(v) > 1L) stats::sd(v) else NA_real_)
    )
    overall <- do.call(data.frame, overall)
    names(overall) <- c("crosstalk", "mean", "sd")
    overall$se <- overall$sd / sqrt(nrep)
    best <- overall[which.min(overall$mean), , drop = FALSE]
    chosen_by_replicate <- do.call(rbind, lapply(
        split(cross_by_replicate, cross_by_replicate$rep),
        function(x) x[which.min(x$rescaled_frobenius), , drop = FALSE]
    ))
    selection_frequency <- as.data.frame(table(
        factor(chosen_by_replicate$crosstalk,
            levels = sort(unique(raw$crosstalk)))
    ))
    names(selection_frequency) <- c("crosstalk", "n_selected")
    selection_frequency$crosstalk <- as.numeric(as.character(selection_frequency$crosstalk))
    selection_frequency$selection_frequency <- selection_frequency$n_selected / nrep
    best_ct <- best$crosstalk[[1]]
    if (length(crosstalk.grid) > 1L && best_ct >= max(crosstalk.grid) &&
        is.finite(ct_limit) && ct_limit - max(crosstalk.grid) > 0.1) {
        warning(
            "Best crosstalk (", best_ct, ") is at the upper end of crosstalk.grid; ",
            "the optimum may be larger. Consider extending the grid up to crosstalk_max(template) = ",
            sprintf("%.3f", ct_limit), ".",
            call. = FALSE
        )
    }
    best_se <- best$se[[1]]
    if (!is.finite(best_se)) best_se <- 0
    overall$near_optimal <- overall$mean <= best$mean[[1]] + best_se
    rownames(best) <- rownames(overall) <- NULL

    list(
        best_crosstalk = best$crosstalk[[1]],
        best = best,
        summary = summary,
        cross_view_summary = overall,
        cross_block_summary = cross,
        raw = raw,
        cross_block_by_replicate = cross_by_replicate,
        selected_by_replicate = chosen_by_replicate,
        selection_frequency = selection_frequency,
        near_optimal = overall[overall$near_optimal, , drop = FALSE],
        view.models = view.models,
        provenance = list(
            package_version = as.character(utils::packageVersion("MimeSys")),
            created_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
            seed = seed,
            crosstalk.grid = crosstalk.grid,
            n = n,
            nrep = nrep,
            p = p,
            max.features.per.view = max.features.per.view,
            feature_select = feature.select,
            feature_jitter = feature.jitter
        )
    )
}

.mimesys_pcl_view_matrices <- function(pcl,
    feature_names = "template_feature") {
    x <- t(as.matrix(pcl$feature_table))
    feature_type <- as.character(pcl$feature_metadata$featureType)
    idx <- split(seq_len(ncol(x)), factor(feature_type, levels = unique(feature_type)))
    lapply(idx, function(j) {
        out <- x[, j, drop = FALSE]
        feature_ids <- pcl$feature_metadata[[feature_names]][j]
        if (!is.null(feature_ids) && !anyNA(feature_ids)) {
            colnames(out) <- make.unique(as.character(feature_ids), sep = "_")
        }
        out
    })
}

#' Plot crosstalk calibration distances
#'
#' Visualizes the output of \code{calibrate_crosstalk()} by plotting the mean
#' rescaled Frobenius distance between observed and simulated correlation blocks
#' across candidate crosstalk values. Lower values indicate better
#' reconstruction
#' of the observed correlation structure.
#'
#' @param calibration A list returned by \code{calibrate_crosstalk()}.
#' @param include_diagonal Logical; if \code{FALSE}, only cross-view blocks are
#'   shown because those are the blocks used to choose the best crosstalk.
#' @param show_sd Logical; if \code{TRUE}, add \code{mean +/- sd} error bars
#'   when
#'   replicate-level standard deviations are available.
#'
#' @return A \code{ggplot} object.
#' @examples
#' data(stelzer_mae)
#' tmpl <- fit_template(stelzer_mae, shared.rank = 2)
#' cal <- calibrate_crosstalk(tmpl, crosstalk.grid = c(0, 0.5, 1), nrep = 2)
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'     plot_crosstalk_calibration(cal)
#' }
#' @export
plot_crosstalk_calibration <- function(
    calibration,
    include_diagonal = FALSE,
    show_sd = TRUE
) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("plot_crosstalk_calibration() requires ggplot2.", call. = FALSE)
    }
    if (!is.list(calibration) || is.null(calibration$summary)) {
        stop("calibration must be a list returned by calibrate_crosstalk().",
            call. = FALSE)
    }
    df <- calibration$summary
    if (!include_diagonal) {
        df <- df[!df$diagonal_block, , drop = FALSE]
    }
    if (nrow(df) == 0L) {
        stop("No calibration rows available to plot.", call. = FALSE)
    }
    df$crosstalk <- as.numeric(as.character(df$crosstalk))
    df$block <- factor(df$block, levels = unique(df$block))
    best <- calibration$best_crosstalk

    p <- ggplot2::ggplot(
        df,
        ggplot2::aes(x = crosstalk, y = mean, color = block, group = block)
    ) +
        ggplot2::geom_line(linewidth = 0.8) +
        ggplot2::geom_point(size = 2.2) +
        ggplot2::geom_vline(
        xintercept = best,
        linetype = "dashed",
        linewidth = 0.5,
        color = "grey35"
    ) +
        ggplot2::labs(
        x = "Crosstalk",
        y = "Rescaled Frobenius distance",
        color = "Correlation block"
    ) +
        ggplot2::theme_classic(base_size = 12) +
        ggplot2::theme(
        legend.position = "bottom",
        strip.background = ggplot2::element_blank()
    )

    if (show_sd && "sd" %in% names(df) && any(is.finite(df$sd))) {
        p <- p +
            ggplot2::geom_errorbar(
            ggplot2::aes(
                ymin = pmax(mean - sd, 0),
                ymax = mean + sd
            ),
            width = 0.02,
            alpha = 0.55
        )
    }
    p
}

#' Summarize a MimeSys simulation object
#'
#' @param sim A list returned by \code{gen_simmba()}.
#'
#' @return A list containing sample, feature, crosstalk, outcome, and signal
#'   summaries: the number of outcome-driving features, the outcome
#'   distribution (event times for survival), case prevalence for binary and
#'   event rate for survival outcomes, and the per-replicate agreement between
#'   the outcome and the true signal.
#' @examples
#' data(stelzer_mae)
#' tmpl <- fit_template(stelzer_mae, shared.rank = 2)
#' sim <- gen_simmba(n = 40, template = tmpl, p = 10, nrep = 2, seed = 1)
#' summarize_simulation(sim)
#' @export
summarize_simulation <- function(sim) {
    if (!is.list(sim) || is.null(sim$trainDat) || is.null(sim$testDat)) {
        stop("sim must be a list returned by gen_simmba().", call. = FALSE)
    }
    first <- sim$trainDat[[1]]
    if (is.null(first$feature_metadata) || is.null(first$sample_metadata)) {
        stop("sim does not look like a MimeSys simulation object.",
            call. = FALSE)
    }
    ft <- first$feature_metadata$featureType
    view_order <- unique(ft)
    feature_counts <- stats::setNames(
        as.integer(table(factor(ft, levels = view_order))), view_order
    )
    train_n <- vapply(sim$trainDat, function(x) nrow(x$sample_metadata),
        integer(1))
    test_n <- vapply(sim$testDat, function(x) nrow(x$sample_metadata),
        integer(1))
    outcome_type <- if (is.null(sim$outcome.type)) "continuous" else sim$outcome.type
    all_meta <- lapply(seq_along(sim$trainDat), function(k) {
        rbind(sim$trainDat[[k]]$sample_metadata, sim$testDat[[k]]$sample_metadata)
    })
    fm <- first$feature_metadata
    active_features <- if (!is.null(fm$outcome_feature)) {
        sum(fm$outcome_feature)
    } else if (!is.null(fm$beta)) {
        sum(fm$beta != 0, na.rm = TRUE)
    } else {
        NA_integer_
    }

    ## Outcome summary and agreement between the outcome and the true signal
    if (identical(outcome_type, "survival")) {
        outcome_values <- unlist(lapply(all_meta, `[[`, "time"), use.names = FALSE)
        event_rate <- vapply(all_meta, function(sm) mean(sm$status), numeric(1))
        signal_cor <- vapply(all_meta, function(sm) {
            ev <- sm$status == 1
            if (sum(ev) < 3L || stats::sd(sm$Xbeta[ev]) == 0) return(NA_real_)
            stats::cor(log(sm$time[ev]), -sm$Xbeta[ev], method = "spearman")
        }, numeric(1))
        prevalence <- NULL
    } else {
        outcome_values <- unlist(lapply(all_meta, `[[`, "Y"), use.names = FALSE)
        event_rate <- NULL
        prevalence <- if (identical(outcome_type, "binary")) {
            vapply(all_meta, function(sm) mean(sm$Y), numeric(1))
        } else {
            NULL
        }
        signal_cor <- vapply(all_meta, function(sm) {
            if (is.null(sm$Y) || is.null(sm$Xbeta) ||
                stats::sd(sm$Y) == 0 || stats::sd(sm$Xbeta) == 0) {
                return(NA_real_)
            }
            stats::cor(sm$Y, sm$Xbeta)
        }, numeric(1))
    }

    out <- list(
        nrep = length(sim$trainDat),
        train_n = stats::setNames(train_n, names(sim$trainDat)),
        test_n = stats::setNames(test_n, names(sim$testDat)),
        feature_counts = feature_counts,
        crosstalk = sim$crosstalk,
        outcome.model = sim$outcome.model,
        outcome.type = sim$outcome.type,
        shared.latent.function = sim$shared.latent.function,
        shared.latent.n.factors = sim$shared.latent.n.factors,
        signal.share = sim$signal.share,
        friedman.share = sim$friedman.share,
        signal.align = sim$signal.align,
        active_features_first_rep = active_features,
        outcome_summary = summary(outcome_values),
        prevalence = prevalence,
        event_rate = event_rate,
        y_signal_cor_mean = mean(signal_cor, na.rm = TRUE),
        y_signal_cor_sd = if (sum(is.finite(signal_cor)) > 1L) stats::sd(signal_cor,
            na.rm = TRUE) else NA_real_
    )
    class(out) <- "mimesys_simulation_summary"
    out
}

#' @export
print.mimesys_simulation_summary <- function(x, ...) {
    cat("MimeSys simulation\n")
    cat("  Replicates: ", x$nrep, "\n", sep = "")
    cat("  Train n: ", paste(range(x$train_n), collapse = "-"), "\n", sep = "")
    cat("  Test n: ", paste(range(x$test_n), collapse = "-"), "\n", sep = "")
    cat("  Features: ", paste(names(x$feature_counts), x$feature_counts,
        sep = "=", collapse = ", "), "\n", sep = "")
    cat("  Crosstalk: ", paste(x$crosstalk, collapse = ", "), "\n", sep = "")
    cat("  Outcome: ", x$outcome.type, " | ", x$outcome.model, "\n", sep = "")
    if (identical(x$outcome.model, "shared_latent")) {
        shared_fun <- if (is.null(x$shared.latent.function)) "linear" else x$shared.latent.function
        shared_n <- if (is.null(x$shared.latent.n.factors)) "all" else x$shared.latent.n.factors
        cat("  Shared latent function: ", shared_fun, "\n", sep = "")
        cat("  Shared latent factors: ", shared_n, "\n", sep = "")
    }
    if (!is.null(x$signal.share)) {
        cat("  Signal share: ", paste(names(x$signal.share),
            round(x$signal.share,
                3), sep = "=", collapse = ", "), "\n", sep = "")
    }
    cat("  Outcome features in first replicate: ", x$active_features_first_rep,
        "\n", sep = "")
    if (!is.null(x$prevalence)) {
        cat("  Case prevalence: ", round(mean(x$prevalence), 3), " (range ",
            paste(round(range(x$prevalence), 3), collapse = "-"), ")\n", sep = "")
    }
    if (!is.null(x$event_rate)) {
        cat("  Event rate: ", round(mean(x$event_rate), 3), " (censored ",
            round(1 - mean(x$event_rate), 3), ")\n", sep = "")
        cat("  Mean Spearman cor(log time, -signal) among events: ",
            round(x$y_signal_cor_mean, 3), "\n", sep = "")
    } else {
        cat("  Mean cor(Y, signal): ", round(x$y_signal_cor_mean, 3), "\n", sep = "")
    }
    invisible(x)
}

#' Compute rescaled Frobenius correlation-block distances
#'
#' @param observed_matrices,simulated_matrices Named sample-by-feature matrix
#'   lists, MAE/SE objects, or compatible containers.
#' @param max_features_per_view Optional feature cap.
#'
#' @return Data frame with one row per upper-triangle view block.
#' @examples
#' data(stelzer_mae)
#' tmpl <- fit_template(stelzer_mae, shared.rank = 2)
#' sim <- gen_simmba(n = 53, template = tmpl, nrep = 1, seed = 1)
#' simulated <- simulation_view_matrices(sim, feature_names =
#' "template_feature")
#' feature_correlation_frobenius(tmpl$observed, simulated)
#' @export
feature_correlation_frobenius <- function(
    observed_matrices,
    simulated_matrices,
    max_features_per_view = NULL
) {
    aligned <- .mimesys_align_observed_simulated(
        observed_matrices,
        simulated_matrices,
        max_features_per_view = max_features_per_view
    )
    .mimesys_frobenius_blocks(aligned$observed, aligned$simulated)
}

## Correlation of the columns of x with those of y; uses the fast complete-data
## path when there are no missing values. Constant columns (possible for sparse
## or discrete features in small simulations) give NA and are skipped.
.mimesys_cor <- function(x, y = x) {
    use <- if (anyNA(x) || anyNA(y)) "pairwise.complete.obs" else "everything"
    suppressWarnings(stats::cor(x, y, use = use))
}

## Observed within- and cross-view correlation blocks for all template
## features, computed once and subset by feature name for each replicate.
.mimesys_observed_cor_cache <- function(observed) {
    views <- names(observed)
    out <- list()
    for (i in seq_along(views)) {
        for (j in seq.int(i, length(views))) {
            out[[paste(views[i], views[j], sep = "\r")]] <-
                .mimesys_cor(observed[[i]], observed[[j]])
        }
    }
    out
}

## Rescaled Frobenius distance between simulated and observed correlation
## blocks. Only the upper-triangular view pairs are computed. When
## observed_cache is given (from .mimesys_observed_cor_cache), observed blocks
## are subset from it instead of being recomputed.
.mimesys_frobenius_blocks <- function(obs, sim, observed_cache = NULL) {
    views <- names(obs)
    rows <- vector("list", length(views) * (length(views) + 1L) / 2L)
    k <- 1L
    for (i in seq_along(views)) {
        for (j in seq.int(i, length(views))) {
            view_x <- views[[i]]
            view_y <- views[[j]]
            key <- paste(view_x, view_y, sep = "\r")
            block_obs <- if (!is.null(observed_cache[[key]])) {
                observed_cache[[key]][colnames(obs[[view_x]]),
                    colnames(obs[[view_y]]), drop = FALSE]
            } else {
                .mimesys_cor(obs[[view_x]], obs[[view_y]])
            }
            block_sim <- .mimesys_cor(sim[[view_x]], sim[[view_y]])
            if (identical(view_x, view_y)) {
                diag(block_obs) <- NA_real_
                diag(block_sim) <- NA_real_
            }
            diff <- block_sim - block_obs
            rows[[k]] <- data.frame(
                block = paste0("cor(", view_x, if (!identical(view_x,
                    view_y)) paste0(",", view_y) else "", ")"),
                view_x = view_x,
                view_y = view_y,
                diagonal_block = identical(view_x, view_y),
                rescaled_frobenius = sqrt(sum(diff^2,
                    na.rm = TRUE) / sum(is.finite(diff))),
                stringsAsFactors = FALSE
            )
            k <- k + 1L
        }
    }
    do.call(rbind, rows)
}

.mimesys_summarize_frobenius <- function(x) {
    out <- stats::aggregate(
        rescaled_frobenius ~ crosstalk + block + view_x + view_y + diagonal_block,
        data = x,
        FUN = function(v) c(mean = mean(v),
            sd = if (length(v) > 1L) stats::sd(v) else NA_real_)
    )
    out <- do.call(data.frame, out)
    names(out) <- sub("rescaled_frobenius\\.", "", names(out))
    out
}
