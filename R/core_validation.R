#' Extract view-specific matrices from a simulated MimeSys replicate
#'
#' @param sim A list returned by template-based `gen_simmba()`.
#' @param rep Replicate index to extract.
#' @param split Either `"train"` or `"test"`.
#' @param feature_names Feature metadata column used to name matrix columns.
#'   Use `"template_feature"` for template-aligned feature names or
#'   `"featureID"` for generated feature IDs.
#'
#' @return A named list of sample-by-feature matrices, one per view.
#' @examples
#' data(stelzer_mae)
#' tmpl <- fit_template(stelzer_mae, shared.rank = 2)
#' sim <- gen_simmba(n = 40, template = tmpl, nrep = 1, seed = 1)
#' views <- simulation_view_matrices(sim, rep = 1, split = "train")
#' vapply(views, dim, integer(2))
#' @export
simulation_view_matrices <- function(
    sim,
    rep = 1,
    split = c("train", "test"),
    feature_names = c("template_feature", "featureID")
) {
    split <- match.arg(split)
    feature_names <- match.arg(feature_names)
    dat_list <- if (identical(split, "train")) sim$trainDat else sim$testDat
    if (!is.list(dat_list) || length(dat_list) < rep) {
        stop("sim does not contain the requested replicate in the requested split.",
            call. = FALSE)
    }

    dat <- dat_list[[rep]]
    if (is.null(dat$feature_table) || is.null(dat$feature_metadata$featureType)) {
        stop("Simulation replicate must contain feature_table and feature_metadata$featureType.",
            call. = FALSE)
    }

    x <- t(as.matrix(dat$feature_table))
    feature_type <- dat$feature_metadata$featureType
    if (length(feature_type) != ncol(x)) {
        stop("feature_metadata$featureType length does not match the feature table.",
            call. = FALSE)
    }

    idx <- split(seq_len(ncol(x)), factor(feature_type, levels = unique(feature_type)))
    out <- lapply(idx, function(j) {
        mat <- x[, j, drop = FALSE]
        names_j <- dat$feature_metadata[[feature_names]][j]
        if (!is.null(names_j) && length(names_j) == ncol(mat) && all(!is.na(names_j))) {
            colnames(mat) <- make.unique(as.character(names_j), sep = "_")
        }
        mat
    })
    out
}

.mimesys_keep_top_variable <- function(x, n_features) {
    if (is.null(n_features) || is.na(n_features)) {
        return(x)
    }
    n_features <- as.integer(n_features)
    if (n_features <= 0) {
        stop("max_features_per_view values must be positive.", call. = FALSE)
    }

    vars <- if (anyNA(x)) apply(x, 2, stats::var,
        na.rm = TRUE) else .col_vars(x)
    names(vars) <- colnames(x)
    vars <- vars[is.finite(vars)]
    vars <- sort(vars, decreasing = TRUE)
    n_keep <- min(n_features, length(vars))
    if (n_keep == 0) {
        stop("No variable features available for one of the views.",
            call. = FALSE)
    }

    x[, names(vars)[seq_len(n_keep)], drop = FALSE]
}

.mimesys_align_observed_simulated <- function(
    observed_matrices,
    simulated_matrices,
    max_features_per_view = NULL
) {
    observed_matrices <- as_mimesys_view_list(observed_matrices)
    simulated_matrices <- as_mimesys_view_list(simulated_matrices)
    views <- intersect(names(observed_matrices), names(simulated_matrices))
    if (length(views) == 0L) {
        stop("No shared view names between observed and simulated matrices.",
            call. = FALSE)
    }

    max_features_per_view <- if (is.null(max_features_per_view)) {
        NULL
    } else {
        unlist(max_features_per_view)
    }
    if (!is.null(max_features_per_view) &&
        (is.null(names(max_features_per_view)) || any(!nzchar(names(max_features_per_view))))) {
        if (length(max_features_per_view) != 1L) {
            stop("max_features_per_view must be named or a single shared value.",
                call. = FALSE)
        }
        max_features_per_view <- stats::setNames(rep(max_features_per_view[[1]],
            length(views)), views)
    }

    obs_out <- sim_out <- vector("list", length(views))
    names(obs_out) <- names(sim_out) <- views
    for (view in views) {
        obs <- as.matrix(observed_matrices[[view]])
        sim <- as.matrix(simulated_matrices[[view]])
        storage.mode(obs) <- "double"
        storage.mode(sim) <- "double"
        common_features <- intersect(colnames(obs), colnames(sim))
        if (length(common_features) < 2L) {
            stop(
                "Fewer than two shared features for view ", view,
                ". Use feature_names = 'template_feature' when extracting simulated matrices, ",
                "or provide matrices with matching column names.",
                call. = FALSE
            )
        }
        obs <- obs[, common_features, drop = FALSE]
        sim <- sim[, common_features, drop = FALSE]
        cap <- if (!is.null(max_features_per_view)) {
            if (view %in% names(max_features_per_view)) max_features_per_view[[view]] else max_features_per_view[[1]]
        } else {
            NULL
        }
        if (!is.null(cap)) {
            obs <- .mimesys_keep_top_variable(obs, cap)
            sim <- sim[, colnames(obs), drop = FALSE]
        }
        obs_out[[view]] <- obs
        sim_out[[view]] <- sim
    }

    list(observed = obs_out, simulated = sim_out)
}

#' Compute blockwise feature-feature correlations across omics views
#'
#' @param matrices Named list of sample-by-feature matrices, one matrix per
#'   omics view.
#' @param max_features_per_view Optional named integer vector/list giving the
#'   number of top-variable features to keep per view. If a single unnamed
#'   value is supplied, the same cap is used for every view.
#' @param use_abs Logical; if `TRUE`, store absolute correlations.
#' @param method Correlation method passed to `stats::cor()`.
#' @param feature_order Feature ordering within each view. `"input"` keeps the
#'   current column order; `"cluster"` orders features by hierarchical
#'   clustering
#'   on within-view correlation distance.
#'
#' @return A data frame with one row per feature-feature correlation and
#'   columns describing the block membership.
#' @examples
#' data(stelzer_mae)
#' views <- as_mimesys_view_list(stelzer_mae)
#' corr <- blockwise_feature_correlation(views, max_features_per_view = 10)
#' head(corr)
#' @export
blockwise_feature_correlation <- function(
    matrices,
    max_features_per_view = NULL,
    use_abs = FALSE,
    method = c("pearson", "spearman"),
    feature_order = c("input", "cluster")
) {
    method <- match.arg(method)
    feature_order <- match.arg(feature_order)
    matrices <- as_mimesys_view_list(matrices)
    if (length(matrices) == 0L) {
        stop("matrices must contain at least one view.", call. = FALSE)
    }

    max_features_per_view <- if (is.null(max_features_per_view)) {
        NULL
    } else {
        unlist(max_features_per_view)
    }
    if (!is.null(max_features_per_view) &&
        (is.null(names(max_features_per_view)) || any(!nzchar(names(max_features_per_view))))) {
        if (length(max_features_per_view) != 1L) {
            stop("max_features_per_view must be named by view or a single shared value.",
                call. = FALSE)
        }
        max_features_per_view <- stats::setNames(rep(max_features_per_view[[1]],
            length(matrices)), names(matrices))
    }

    view_names_in <- names(matrices)
    matrices <- lapply(view_names_in, function(view) {
        x <- as.matrix(matrices[[view]])
        storage.mode(x) <- "double"
        if (is.null(colnames(x))) {
            colnames(x) <- paste0(view, "_feature_", seq_len(ncol(x)))
        }
        cap <- if (!is.null(max_features_per_view)) {
            if (view %in% names(max_features_per_view)) max_features_per_view[[view]] else max_features_per_view[[1]]
        } else {
            NULL
        }
        x <- .mimesys_keep_top_variable(x, cap)
        if (identical(feature_order, "cluster") && ncol(x) > 2L) {
            ord <- .mimesys_cluster_feature_order(x)
            x <- x[, ord, drop = FALSE]
        }
        x
    })
    matrices <- stats::setNames(matrices, view_names_in)

    all_x <- do.call(cbind, matrices)
    cor_mat <- suppressWarnings(stats::cor(all_x, use = "pairwise.complete.obs",
        method = method))
    if (use_abs) {
        cor_mat <- abs(cor_mat)
    }

    view_sizes <- vapply(matrices, ncol, integer(1))
    view_names <- names(matrices)
    starts <- cumsum(c(1L, head(view_sizes, -1L)))
    ends <- cumsum(view_sizes)
    feature_view <- rep(view_names, view_sizes)
    feature_order_index <- unlist(lapply(view_sizes, seq_len),
        use.names = FALSE)

    ## Build the long table in one vectorized step, ordered row by row as before.
    p_all <- nrow(cor_mat)
    i <- rep.int(seq_len(p_all), p_all)
    j <- rep(seq_len(p_all), each = p_all)
    out <- data.frame(
        feature_x = rownames(cor_mat)[i],
        feature_y = colnames(cor_mat)[j],
        view_x = feature_view[i],
        view_y = feature_view[j],
        order_x = feature_order_index[i],
        order_y = feature_order_index[j],
        block = paste(feature_view[i], feature_view[j], sep = " vs "),
        diagonal_block = feature_view[i] == feature_view[j],
        feature_diagonal = i == j,
        correlation = as.vector(cor_mat),
        use_abs = use_abs,
        method = method,
        stringsAsFactors = FALSE
    )
    out <- out[order(i, j), , drop = FALSE]
    rownames(out) <- NULL
    out$view_x <- factor(out$view_x, levels = view_names)
    out$view_y <- factor(out$view_y, levels = view_names)
    out$block <- factor(
        out$block,
        levels = as.vector(outer(view_names, view_names, paste, sep = " vs "))
    )
    attr(out, "view_boundaries") <- data.frame(
        view = view_names,
        start = starts,
        end = ends,
        stringsAsFactors = FALSE
    )
    out
}

.mimesys_cluster_feature_order <- function(x) {
    cmat <- stats::cor(x, use = "pairwise.complete.obs")
    cmat[!is.finite(cmat)] <- 0
    dist_mat <- stats::as.dist(1 - abs(cmat))
    stats::hclust(dist_mat)$order
}

#' Plot blockwise feature-feature correlations across omics views
#'
#' @param corr_df Data frame returned by `blockwise_feature_correlation()`.
#' @param title Plot title.
#' @param use_abs Optional logical controlling the color scale. If `NULL`, this
#'   is read from `corr_df`.
#' @param layout Either `"matrix"` for one full matrix with block boundaries or
#'   `"faceted"` for one facet per view-pair block.
#' @param show_feature_diagonal Logical; if `FALSE`, the feature
#'   self-correlation
#'   diagonal is hidden so block structure is easier to see.
#' @param color_limit Optional positive number used as the symmetric color limit
#'   for signed correlations or upper color limit for absolute correlations.
#'
#' @return A `ggplot` object.
#' @examples
#' data(stelzer_mae)
#' views <- as_mimesys_view_list(stelzer_mae)
#' corr <- blockwise_feature_correlation(views, max_features_per_view = 10)
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'     plot_blockwise_feature_correlation(corr)
#' }
#' @export
plot_blockwise_feature_correlation <- function(
    corr_df,
    title = "Blockwise feature correlation",
    use_abs = NULL,
    layout = c("matrix", "faceted"),
    show_feature_diagonal = TRUE,
    color_limit = NULL
) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("plot_blockwise_feature_correlation() requires ggplot2.",
            call. = FALSE)
    }
    layout <- match.arg(layout)
    needed <- c("order_x", "order_y", "view_x", "view_y", "correlation",
        "feature_diagonal")
    if (!is.data.frame(corr_df) || !all(needed %in% names(corr_df))) {
        stop("corr_df must come from blockwise_feature_correlation().",
            call. = FALSE)
    }
    plot_df <- corr_df
    if (!show_feature_diagonal) {
        plot_df$correlation[plot_df$feature_diagonal] <- NA_real_
    }
    if (is.null(use_abs)) {
        use_abs <- isTRUE(plot_df$use_abs[[1]])
    }

    if (identical(layout, "faceted")) {
        p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = order_y, y = order_x,
            fill = correlation)) +
            ggplot2::geom_tile() +
            ggplot2::facet_grid(ggplot2::vars(view_x), ggplot2::vars(view_y),
                scales = "free", space = "free") +
            ggplot2::scale_y_reverse(expand = c(0, 0)) +
            ggplot2::scale_x_continuous(expand = c(0, 0)) +
            ggplot2::labs(title = title, x = "Feature", y = "Feature",
                fill = "Correlation") +
            ggplot2::theme_bw(base_size = 11) +
            ggplot2::theme(
            plot.title = ggplot2::element_text(face = "bold"),
            axis.text = ggplot2::element_blank(),
            axis.ticks = ggplot2::element_blank(),
            panel.grid = ggplot2::element_blank(),
            strip.background = ggplot2::element_rect(fill = "white",
                color = "grey75")
        )
    } else {
        boundaries <- attr(corr_df, "view_boundaries")
        if (is.null(boundaries)) {
            boundaries <- data.frame(view = unique(as.character(corr_df$view_x)),
                start = NA, end = NA)
        }
        plot_df$feature_y_index <- as.numeric(factor(plot_df$feature_y,
            levels = unique(plot_df$feature_y)))
        plot_df$feature_x_index <- as.numeric(factor(plot_df$feature_x,
            levels = unique(plot_df$feature_x)))
        p <- ggplot2::ggplot(
            plot_df,
            ggplot2::aes(
                x = feature_y_index,
                y = feature_x_index,
                fill = correlation
            )
        ) +
            ggplot2::geom_tile() +
            ggplot2::scale_y_reverse(expand = c(0, 0)) +
            ggplot2::scale_x_continuous(expand = c(0, 0)) +
            ggplot2::labs(title = title, x = "Feature", y = "Feature",
                fill = "Correlation") +
            ggplot2::theme_bw(base_size = 11) +
            ggplot2::theme(
            plot.title = ggplot2::element_text(face = "bold"),
            axis.text = ggplot2::element_blank(),
            axis.ticks = ggplot2::element_blank(),
            panel.grid = ggplot2::element_blank()
        )
        if (all(is.finite(boundaries$end))) {
            p <- p +
                ggplot2::geom_vline(xintercept = boundaries$end + 0.5,
                    color = "grey30", linewidth = 0.25) +
                ggplot2::geom_hline(yintercept = boundaries$end + 0.5,
                    color = "grey30", linewidth = 0.25)
        }
    }

    if (isTRUE(use_abs)) {
        p + ggplot2::scale_fill_gradient(
            low = "white",
            high = "#2166ac",
            limits = c(0, if (is.null(color_limit)) 1 else color_limit),
            na.value = "grey95"
        )
    } else {
        p + ggplot2::scale_fill_gradient2(
            low = "#2166ac",
            mid = "white",
            high = "#b2182b",
            midpoint = 0,
            limits = c(-if (is.null(color_limit)) 1 else color_limit,
                if (is.null(color_limit)) 1 else color_limit),
            na.value = "grey95"
        )
    }
}

#' Summarize blockwise feature correlations
#'
#' @param corr_df Data frame returned by `blockwise_feature_correlation()`.
#' @param unique_pairs Logical; if `TRUE`, symmetric off-diagonal blocks are
#'   collapsed so each view pair appears once.
#'
#' @return A data frame with one row per view-pair block.
#' @examples
#' data(stelzer_mae)
#' views <- as_mimesys_view_list(stelzer_mae)
#' corr <- blockwise_feature_correlation(views, max_features_per_view = 10)
#' summarize_blockwise_feature_correlation(corr)
#' @export
summarize_blockwise_feature_correlation <- function(corr_df,
    unique_pairs = TRUE) {
    needed <- c("view_x", "view_y", "correlation", "feature_diagonal")
    if (!is.data.frame(corr_df) || !all(needed %in% names(corr_df))) {
        stop("corr_df must come from blockwise_feature_correlation().",
            call. = FALSE)
    }
    x <- corr_df[!corr_df$feature_diagonal, , drop = FALSE]
    x$view_x_chr <- as.character(x$view_x)
    x$view_y_chr <- as.character(x$view_y)
    if (unique_pairs) {
        x$view_a <- pmin(x$view_x_chr, x$view_y_chr)
        x$view_b <- pmax(x$view_x_chr, x$view_y_chr)
    } else {
        x$view_a <- x$view_x_chr
        x$view_b <- x$view_y_chr
    }
    x$diagonal_block <- x$view_a == x$view_b
    out <- stats::aggregate(
        correlation ~ view_a + view_b + diagonal_block,
        data = x,
        FUN = function(v) {
            c(
                mean = mean(v, na.rm = TRUE),
                mean_abs = mean(abs(v), na.rm = TRUE),
                sd = stats::sd(v, na.rm = TRUE),
                n = sum(is.finite(v))
            )
        }
    )
    out <- do.call(data.frame, out)
    names(out) <- sub("correlation\\.", "", names(out))
    names(out)[names(out) == "view_a"] <- "view_x"
    names(out)[names(out) == "view_b"] <- "view_y"
    rownames(out) <- NULL
    out
}

#' Plot a compact heatmap of block-level correlations
#'
#' @param block_summary Data frame returned by
#'   `summarize_blockwise_feature_correlation()`.
#' @param value Column to plot.
#' @param title Plot title.
#'
#' @return A `ggplot` object.
#' @examples
#' data(stelzer_mae)
#' views <- as_mimesys_view_list(stelzer_mae)
#' corr <- blockwise_feature_correlation(views, max_features_per_view = 10)
#' block <- summarize_blockwise_feature_correlation(corr)
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'     plot_blockwise_correlation_summary(block)
#' }
#' @export
plot_blockwise_correlation_summary <- function(
    block_summary,
    value = c("mean_abs", "mean", "sd"),
    title = "Block-level feature correlation summary"
) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("plot_blockwise_correlation_summary() requires ggplot2.",
            call. = FALSE)
    }
    value <- match.arg(value)
    if (!is.data.frame(block_summary) || !all(c("view_x", "view_y",
        value) %in% names(block_summary))) {
        stop("block_summary must come from summarize_blockwise_feature_correlation().",
            call. = FALSE)
    }
    plot_df <- block_summary
    plot_df$value <- plot_df[[value]]
    views <- unique(c(as.character(plot_df$view_x),
        as.character(plot_df$view_y)))
    plot_df$view_x <- factor(plot_df$view_x, levels = views)
    plot_df$view_y <- factor(plot_df$view_y, levels = views)

    ggplot2::ggplot(plot_df, ggplot2::aes(x = view_y, y = view_x,
        fill = value)) +
        ggplot2::geom_tile(color = "white", linewidth = 0.5) +
        ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f", value)),
            size = 3) +
        ggplot2::scale_fill_gradient(low = "white", high = "#2166ac",
            na.value = "grey95") +
        ggplot2::labs(title = title, x = NULL, y = NULL, fill = value) +
        ggplot2::theme_bw(base_size = 12) +
        ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold"),
        panel.grid = ggplot2::element_blank()
    )
}

#' Compare observed and simulated block-level correlations
#'
#' @param observed_summary,simulated_summary Data frames from
#'   `summarize_blockwise_feature_correlation()`.
#' @param value Column to compare.
#'
#' @return A data frame with observed, simulated, and difference values.
#' @examples
#' data(stelzer_mae)
#' tmpl <- fit_template(stelzer_mae, shared.rank = 2)
#' sim <- gen_simmba(n = 53, template = tmpl, nrep = 1, seed = 1)
#' obs <- summarize_blockwise_feature_correlation(
#'     blockwise_feature_correlation(tmpl$observed)
#' )
#' simulated <- summarize_blockwise_feature_correlation(
#'     blockwise_feature_correlation(simulation_view_matrices(sim))
#' )
#' compare_blockwise_correlation_summary(obs, simulated)
#' @export
compare_blockwise_correlation_summary <- function(
    observed_summary,
    simulated_summary,
    value = c("mean_abs", "mean", "sd")
) {
    value <- match.arg(value)
    needed <- c("view_x", "view_y", value)
    if (!all(needed %in% names(observed_summary)) || !all(needed %in% names(simulated_summary))) {
        stop("Both summaries must contain view_x, view_y, and the requested value column.",
            call. = FALSE)
    }
    obs <- observed_summary[, needed, drop = FALSE]
    sim <- simulated_summary[, needed, drop = FALSE]
    names(obs)[names(obs) == value] <- "observed"
    names(sim)[names(sim) == value] <- "simulated"
    out <- merge(obs, sim, by = c("view_x", "view_y"), all = FALSE)
    out$difference <- out$simulated - out$observed
    out
}

#' Compare observed and simulated feature-correlation structure
#'
#' Correlates the vector of observed feature-feature correlations with the
#' matching vector of simulated feature-feature correlations within each view
#' block. This asks whether the simulated data preserve the pattern of the
#' correlation matrix, not only its average magnitude. Each pair of distinct
#' features is counted once: self-correlations are left out, and each
#' cross-view block appears once (views in alphabetical order).
#'
#' @param observed_corr,simulated_corr Data frames returned by
#'   `blockwise_feature_correlation()` for observed and simulated matrices.
#' @param return_pairs Logical; if `TRUE`, also return the matched
#'   feature-pair correlations used to compute the block-level summaries.
#'
#' @return If `return_pairs = FALSE`, a data frame with one row per view-pair
#'   block and columns for Pearson/Spearman corr-of-corr, MAE, RMSE, and number
#'   of feature pairs. If `return_pairs = TRUE`, a list with `summary` and
#'   `pairs`.
#' @examples
#' data(stelzer_mae)
#' tmpl <- fit_template(stelzer_mae, shared.rank = 2)
#' sim <- gen_simmba(n = 53, template = tmpl, nrep = 1, seed = 1)
#' obs <- blockwise_feature_correlation(tmpl$observed)
#' simulated <- blockwise_feature_correlation(
#'     simulation_view_matrices(sim, feature_names = "template_feature")
#' )
#' compare_feature_correlation_structure(obs, simulated)
#' @export
compare_feature_correlation_structure <- function(observed_corr, simulated_corr,
    return_pairs = FALSE) {
    key_cols <- c("feature_x", "feature_y", "view_x", "view_y")
    needed <- c(key_cols, "correlation")
    if (!all(needed %in% names(observed_corr)) || !all(needed %in% names(simulated_corr))) {
        stop("observed_corr and simulated_corr must come from blockwise_feature_correlation().",
            call. = FALSE)
    }
    ## Each feature pair once: drop self-correlations and keep one orientation
    ## of every pair (blockwise_feature_correlation() lists both), chosen by
    ## view and feature name so that the two tables agree
    one_orientation <- function(x) {
        vx <- as.character(x$view_x); vy <- as.character(x$view_y)
        fx <- as.character(x$feature_x); fy <- as.character(x$feature_y)
        keep <- (vx < vy) | (vx == vy & fx < fy)
        x[keep, needed, drop = FALSE]
    }
    obs <- one_orientation(observed_corr)
    sim <- one_orientation(simulated_corr)
    names(obs)[names(obs) == "correlation"] <- "observed"
    names(sim)[names(sim) == "correlation"] <- "simulated"
    pairs <- merge(obs, sim, by = key_cols, all = FALSE)
    pairs <- pairs[is.finite(pairs$observed) & is.finite(pairs$simulated), ,
        drop = FALSE]
    pairs$abs_error <- abs(pairs$simulated - pairs$observed)
    pairs$sq_error <- (pairs$simulated - pairs$observed)^2

    split_idx <- split(seq_len(nrow(pairs)), paste(pairs$view_x, pairs$view_y,
        sep = "\r"))
    summary <- do.call(rbind, lapply(split_idx, function(idx) {
        x <- pairs[idx, , drop = FALSE]
        data.frame(
            view_x = x$view_x[[1]],
            view_y = x$view_y[[1]],
            diagonal_block = identical(as.character(x$view_x[[1]]),
                as.character(x$view_y[[1]])),
            pearson = suppressWarnings(stats::cor(x$observed, x$simulated,
                method = "pearson")),
            spearman = suppressWarnings(stats::cor(x$observed, x$simulated,
                method = "spearman")),
            mae = mean(x$abs_error),
            rmse = sqrt(mean(x$sq_error)),
            n = nrow(x),
            stringsAsFactors = FALSE
        )
    }))
    rownames(summary) <- NULL
    if (return_pairs) {
        list(summary = summary, pairs = pairs)
    } else {
        summary
    }
}

#' Plot observed-vs-simulated block correlation differences
#'
#' @param comparison_df Data frame from
#'   `compare_blockwise_correlation_summary()`.
#' @param title Plot title.
#' @param difference_limit Symmetric color limit for simulated-minus-observed
#'   differences.
#'
#' @return A `ggplot` object.
#' @examples
#' data(stelzer_mae)
#' tmpl <- fit_template(stelzer_mae, shared.rank = 2)
#' sim <- gen_simmba(n = 53, template = tmpl, nrep = 1, seed = 1)
#' obs <- summarize_blockwise_feature_correlation(
#'     blockwise_feature_correlation(tmpl$observed)
#' )
#' simulated <- summarize_blockwise_feature_correlation(
#'     blockwise_feature_correlation(simulation_view_matrices(sim))
#' )
#' cmp <- compare_blockwise_correlation_summary(obs, simulated)
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'     plot_blockwise_correlation_difference(cmp)
#' }
#' @export
plot_blockwise_correlation_difference <- function(
    comparison_df,
    title = "Simulated minus observed block correlation",
    difference_limit = NULL
) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("plot_blockwise_correlation_difference() requires ggplot2.",
            call. = FALSE)
    }
    if (!is.data.frame(comparison_df) || !all(c("view_x", "view_y",
        "difference") %in% names(comparison_df))) {
        stop("comparison_df must come from compare_blockwise_correlation_summary().",
            call. = FALSE)
    }
    views <- unique(c(as.character(comparison_df$view_x),
        as.character(comparison_df$view_y)))
    comparison_df$view_x <- factor(comparison_df$view_x, levels = views)
    comparison_df$view_y <- factor(comparison_df$view_y, levels = views)
    lim <- if (is.null(difference_limit)) {
        max(abs(comparison_df$difference), na.rm = TRUE)
    } else {
        difference_limit
    }
    if (!is.finite(lim) || lim == 0) {
        lim <- 1
    }

    ggplot2::ggplot(comparison_df, ggplot2::aes(x = view_y, y = view_x,
        fill = difference)) +
        ggplot2::geom_tile(color = "white", linewidth = 0.5) +
        ggplot2::geom_text(ggplot2::aes(label = sprintf("%+.3f", difference)),
            size = 3) +
        ggplot2::scale_fill_gradient2(
        low = "#2166ac",
        mid = "white",
        high = "#b2182b",
        midpoint = 0,
        limits = c(-lim, lim),
        na.value = "grey95"
    ) +
        ggplot2::labs(title = title, x = NULL, y = NULL, fill = "Difference") +
        ggplot2::theme_bw(base_size = 12) +
        ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold"),
        panel.grid = ggplot2::element_blank()
    )
}

#' Save a blockwise feature-correlation plot
#'
#' @param matrices Named list of sample-by-feature matrices.
#' @param plot_dir Directory where the PDF should be saved.
#' @param prefix File prefix for the saved PDF.
#' @param title Plot title.
#' @param max_features_per_view,use_abs,method Passed to
#'   `blockwise_feature_correlation()`.
#' @param layout,show_feature_diagonal,color_limit Passed to
#'   `plot_blockwise_feature_correlation()`.
#' @param width,height PDF dimensions in inches.
#'
#' @return Path to the saved PDF.
#' @examples
#' data(stelzer_mae)
#' views <- as_mimesys_view_list(stelzer_mae)
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'     path <- save_blockwise_feature_correlation_plot(
#'         views, plot_dir = tempdir(), max_features_per_view = 10
#'     )
#'     file.exists(path)
#' }
#' @export
save_blockwise_feature_correlation_plot <- function(
    matrices,
    plot_dir,
    prefix = "blockwise_feature_correlation",
    title = "Blockwise feature correlation",
    max_features_per_view = NULL,
    use_abs = FALSE,
    method = c("pearson", "spearman"),
    layout = c("matrix", "faceted"),
    show_feature_diagonal = TRUE,
    color_limit = NULL,
    width = 8,
    height = 7
) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("save_blockwise_feature_correlation_plot() requires ggplot2.",
            call. = FALSE)
    }
    method <- match.arg(method)
    layout <- match.arg(layout)
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
    corr_df <- blockwise_feature_correlation(
        matrices = matrices,
        max_features_per_view = max_features_per_view,
        use_abs = use_abs,
        method = method
    )
    p <- plot_blockwise_feature_correlation(
        corr_df,
        title = title,
        use_abs = use_abs,
        layout = layout,
        show_feature_diagonal = show_feature_diagonal,
        color_limit = color_limit
    )
    path <- file.path(plot_dir, paste0(prefix, ".pdf"))
    ggplot2::ggsave(path, p, width = width, height = height)
    path
}
