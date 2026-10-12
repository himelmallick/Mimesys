#' Summarize template-based simulation diagnostics
#'
#' Computes the core diagnostics needed to validate whether a template-based
#' `gen_simmba()` run produced the expected cross-view dependence, aligned
#' outcome signal, and recoverable response.
#'
#' @param sim A list returned by the template-based `gen_simmba()`.
#' @param view_pairs Optional two-column character matrix/data frame of view
#'   pairs. If `NULL`, all pairwise view combinations are used.
#'
#' @return A data frame with one row per replicate and view pair.
#' @export
validate_template_simulation <- function(sim, view_pairs = NULL) {
  if (!is.list(sim) || is.null(sim$trainDat) || is.null(sim$testDat)) {
    stop("sim must be a list returned by template-based gen_simmba().", call. = FALSE)
  }
  if (is.null(sim$template) || !inherits(sim$template, "mimesys_template")) {
    stop("sim$template must be a mimesys_template object.", call. = FALSE)
  }

  views <- sim$template$views
  if (length(views) < 2) {
    stop("At least two views are required for crosstalk diagnostics.", call. = FALSE)
  }
  if (is.null(view_pairs)) {
    view_pairs <- utils::combn(views, 2)
  } else {
    view_pairs <- as.matrix(view_pairs)
    if (nrow(view_pairs) != 2 && ncol(view_pairs) == 2) {
      view_pairs <- t(view_pairs)
    }
    if (nrow(view_pairs) != 2) {
      stop("view_pairs must have two rows or two columns.", call. = FALSE)
    }
  }

  rows <- vector("list", length(sim$trainDat) * ncol(view_pairs))
  idx <- 1L
  for (rep_idx in seq_along(sim$trainDat)) {
    dat <- sim$trainDat[[rep_idx]]
    fm <- dat$feature_metadata
    x <- t(as.matrix(dat$feature_table))
    sm <- dat$sample_metadata

    if (is.null(fm$featureType) || is.null(fm$beta)) {
      stop("feature_metadata must contain featureType and beta.", call. = FALSE)
    }

    for (pair_idx in seq_len(ncol(view_pairs))) {
      a <- view_pairs[1, pair_idx]
      b <- view_pairs[2, pair_idx]
      sel_a <- fm$featureType == a
      sel_b <- fm$featureType == b
      if (!any(sel_a) || !any(sel_b)) {
        stop(sprintf("Missing features for view pair %s:%s.", a, b), call. = FALSE)
      }

      pred_a <- as.numeric(x[, sel_a, drop = FALSE] %*% fm$beta[sel_a])
      pred_b <- as.numeric(x[, sel_b, drop = FALSE] %*% fm$beta[sel_b])
      y_signal <- if (!is.null(sm$Xbeta)) sm$Xbeta else rep(NA_real_, nrow(x))

      feature_cor <- suppressWarnings(stats::cor(
        x[, sel_a, drop = FALSE],
        x[, sel_b, drop = FALSE],
        use = "pairwise.complete.obs"
      ))
      feature_cor <- as.numeric(feature_cor)
      feature_cor <- feature_cor[is.finite(feature_cor)]
      mean_abs_feature_cor <- if (length(feature_cor) == 0L) {
        NA_real_
      } else {
        mean(abs(feature_cor))
      }

      rows[[idx]] <- data.frame(
        rep = rep_idx,
        view_a = a,
        view_b = b,
        crosstalk = paste(sim$crosstalk, collapse = ","),
        de_align = if (is.null(sim$de.align)) NA_real_ else sim$de.align,
        mean_abs_feature_cor = mean_abs_feature_cor,
        signal_agreement = suppressWarnings(stats::cor(pred_a, pred_b)),
        y_signal_cor = if ("Y" %in% names(sm)) {
          suppressWarnings(stats::cor(sm$Y, y_signal))
        } else {
          NA_real_
        },
        event_rate = if ("Y" %in% names(sm) && all(sm$Y %in% c(0, 1))) mean(sm$Y) else NA_real_,
        observed_event_rate = if ("status" %in% names(sm)) mean(sm$status) else NA_real_,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Aggregate template simulation diagnostics
#'
#' @param diagnostics A data frame returned by `validate_template_simulation()`.
#'
#' @return A compact data frame grouped by view pair.
#' @export
summarize_template_diagnostics <- function(diagnostics) {
  needed <- c("view_a", "view_b", "mean_abs_feature_cor", "signal_agreement", "y_signal_cor")
  if (!is.data.frame(diagnostics) || !all(needed %in% names(diagnostics))) {
    stop("diagnostics must come from validate_template_simulation().", call. = FALSE)
  }

  diagnostics$effect_alignment <- .mimesys_effect_alignment_label(diagnostics)

  grouping <- c("crosstalk", "effect_alignment", "view_a", "view_b")
  grouping <- grouping[grouping %in% names(diagnostics)]
  key <- interaction(diagnostics[, grouping, drop = FALSE], drop = TRUE, sep = "\r")
  split_idx <- split(seq_len(nrow(diagnostics)), key)

  rows <- lapply(split_idx, function(idx) {
    first <- diagnostics[idx[1], grouping, drop = FALSE]
    data.frame(
      first,
      mean_abs_feature_cor = {
        x <- diagnostics$mean_abs_feature_cor[idx]
        if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
      },
      signal_agreement = {
        x <- diagnostics$signal_agreement[idx]
        if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
      },
      y_signal_cor = {
        x <- diagnostics$y_signal_cor[idx]
        if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
      },
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.mimesys_effect_alignment_label <- function(x) {
  if ("de_align_label" %in% names(x)) {
    return(gsub("_", " ", as.character(x$de_align_label)))
  }
  if ("effect_alignment" %in% names(x)) {
    return(as.character(x$effect_alignment))
  }
  if ("de_align" %in% names(x)) {
    return(ifelse(is.na(x$de_align), "random effects", "aligned effects"))
  }
  rep("not specified", nrow(x))
}

#' Extract view-specific matrices from a simulated MimeSys replicate
#'
#' @param sim A list returned by template-based `gen_simmba()`.
#' @param rep Replicate index to extract.
#' @param split Either `"train"` or `"test"`.
#' @param feature_names Feature metadata column used to name matrix columns.
#'   Use `"template_feature"` for simulated feature identifiers or
#'   `"featureID"` for the original generated feature IDs.
#'
#' @return A named list of sample-by-feature matrices, one per view.
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
    stop("sim does not contain the requested replicate in the requested split.", call. = FALSE)
  }

  dat <- dat_list[[rep]]
  if (is.null(dat$feature_table) || is.null(dat$feature_metadata$featureType)) {
    stop("Simulation replicate must contain feature_table and feature_metadata$featureType.", call. = FALSE)
  }

  x <- t(as.matrix(dat$feature_table))
  feature_type <- dat$feature_metadata$featureType
  if (length(feature_type) != ncol(x)) {
    stop("feature_metadata$featureType length does not match the feature table.", call. = FALSE)
  }

  idx <- split(seq_len(ncol(x)), feature_type)
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

  vars <- apply(x, 2, stats::var, na.rm = TRUE)
  vars <- vars[is.finite(vars)]
  vars <- sort(vars, decreasing = TRUE)
  n_keep <- min(n_features, length(vars))
  if (n_keep == 0) {
    stop("No variable features available for one of the views.", call. = FALSE)
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
    stop("No shared view names between observed and simulated matrices.", call. = FALSE)
  }

  max_features_per_view <- if (is.null(max_features_per_view)) {
    NULL
  } else {
    unlist(max_features_per_view)
  }
  if (!is.null(max_features_per_view) &&
      (is.null(names(max_features_per_view)) || any(!nzchar(names(max_features_per_view))))) {
    if (length(max_features_per_view) != 1L) {
      stop("max_features_per_view must be named or a single shared value.", call. = FALSE)
    }
    max_features_per_view <- stats::setNames(rep(max_features_per_view[[1]], length(views)), views)
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
        ". For simulated data, call simulation_view_matrices(..., feature_names = 'template_feature').",
        call. = FALSE
      )
    }
    obs <- obs[, common_features, drop = FALSE]
    sim <- sim[, common_features, drop = FALSE]
    finite_enough <- apply(obs, 2, function(x) sum(is.finite(x)) >= 2) &
      apply(sim, 2, function(x) sum(is.finite(x)) >= 2)
    obs <- obs[, finite_enough, drop = FALSE]
    sim <- sim[, finite_enough, drop = FALSE]
    if (ncol(obs) < 2L) {
      stop("Fewer than two finite shared features remain for view ", view, ".", call. = FALSE)
    }

    cap <- if (!is.null(max_features_per_view) && view %in% names(max_features_per_view)) {
      max_features_per_view[[view]]
    } else if (!is.null(max_features_per_view) && length(max_features_per_view) == 1L) {
      max_features_per_view[[1]]
    } else {
      NA_integer_
    }

    if (!is.na(cap)) {
      vars <- apply(obs, 2, stats::var, na.rm = TRUE)
      vars <- sort(vars[is.finite(vars)], decreasing = TRUE)
      keep <- names(vars)[seq_len(min(as.integer(cap), length(vars)))]
      obs <- obs[, keep, drop = FALSE]
      sim <- sim[, keep, drop = FALSE]
    }

    obs_out[[view]] <- obs
    sim_out[[view]] <- sim
  }

  list(observed = obs_out, simulated = sim_out)
}

.mimesys_scale_using_observed <- function(obs, sim) {
  obs[!is.finite(obs)] <- NA_real_
  sim[!is.finite(sim)] <- NA_real_
  keep <- colSums(!is.na(obs)) >= 2 & colSums(!is.na(sim)) >= 2
  obs <- obs[, keep, drop = FALSE]
  sim <- sim[, keep, drop = FALSE]
  if (ncol(obs) == 0L) {
    stop("No usable features remain after removing non-finite values.", call. = FALSE)
  }

  center <- apply(obs, 2, mean, na.rm = TRUE)
  scale <- apply(obs, 2, stats::sd, na.rm = TRUE)
  scale[!is.finite(scale) | scale == 0] <- 1

  for (j in seq_len(ncol(obs))) {
    obs[is.na(obs[, j]), j] <- center[[j]]
    sim[is.na(sim[, j]), j] <- center[[j]]
  }

  obs_z <- sweep(sweep(obs, 2, center, "-"), 2, scale, "/")
  sim_z <- sweep(sweep(sim, 2, center, "-"), 2, scale, "/")
  obs_z[!is.finite(obs_z)] <- 0
  sim_z[!is.finite(sim_z)] <- 0
  list(observed = obs_z, simulated = sim_z)
}

#' Compute a real-vs-simulated multi-omics PCA diagnostic
#'
#' Features are aligned by column name within each view, standardized using the
#' observed data, optionally weighted so views contribute equally, concatenated,
#' and then projected into a joint PCA.
#'
#' @param observed_matrices Named list of observed sample-by-feature matrices.
#' @param simulated_matrices Named list of simulated sample-by-feature matrices.
#' @param max_features_per_view Optional cap on top-variable features per view.
#' @param equal_view_weight Logical; if `TRUE`, each view is divided by the
#'   square root of its retained feature count before concatenation.
#'
#' @return A list with PCA scores, percent variance explained, loadings, and the
#'   aligned matrices used.
#' @export
multiomics_pca_diagnostic <- function(
  observed_matrices,
  simulated_matrices,
  max_features_per_view = NULL,
  equal_view_weight = TRUE
) {
  aligned <- .mimesys_align_observed_simulated(
    observed_matrices,
    simulated_matrices,
    max_features_per_view = max_features_per_view
  )

  scaled <- lapply(names(aligned$observed), function(view) {
    mats <- .mimesys_scale_using_observed(aligned$observed[[view]], aligned$simulated[[view]])
    if (isTRUE(equal_view_weight)) {
      weight <- sqrt(ncol(mats$observed))
      mats$observed <- mats$observed / weight
      mats$simulated <- mats$simulated / weight
    }
    colnames(mats$observed) <- paste(view, colnames(mats$observed), sep = "::")
    colnames(mats$simulated) <- paste(view, colnames(mats$simulated), sep = "::")
    mats
  })
  names(scaled) <- names(aligned$observed)

  x_obs <- do.call(cbind, lapply(scaled, `[[`, "observed"))
  x_sim <- do.call(cbind, lapply(scaled, `[[`, "simulated"))
  x_all <- rbind(x_obs, x_sim)
  keep <- apply(x_all, 2, function(x) all(is.finite(x)) && stats::var(x) > 0)
  x_all <- x_all[, keep, drop = FALSE]
  if (ncol(x_all) < 2L) {
    stop("Fewer than two finite, variable features remain for PCA.", call. = FALSE)
  }
  type <- c(rep("observed", nrow(x_obs)), rep("simulated", nrow(x_sim)))
  sample_id <- c(
    if (is.null(rownames(x_obs))) paste0("observed_", seq_len(nrow(x_obs))) else rownames(x_obs),
    if (is.null(rownames(x_sim))) paste0("simulated_", seq_len(nrow(x_sim))) else rownames(x_sim)
  )

  pca <- stats::prcomp(x_all, center = FALSE, scale. = FALSE)
  var_explained <- 100 * (pca$sdev^2) / sum(pca$sdev^2)
  scores <- data.frame(
    sample_id = sample_id,
    data_type = type,
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    stringsAsFactors = FALSE
  )

  list(
    scores = scores,
    var_explained = var_explained,
    rotation = pca$rotation,
    aligned = aligned
  )
}

#' Plot a real-vs-simulated multi-omics PCA diagnostic
#'
#' @param pca_result Result from `multiomics_pca_diagnostic()`.
#' @param title Plot title.
#'
#' @return A `ggplot` object.
#' @export
plot_multiomics_pca_diagnostic <- function(
  pca_result,
  title = "Observed versus simulated multi-omics PCA"
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot_multiomics_pca_diagnostic() requires ggplot2.", call. = FALSE)
  }
  scores <- pca_result$scores
  ve <- pca_result$var_explained
  ggplot2::ggplot(
    scores,
    ggplot2::aes(x = PC1, y = PC2, color = data_type)
  ) +
    ggplot2::geom_point(alpha = 0.65, size = 2) +
    ggplot2::stat_ellipse(linewidth = 0.7, alpha = 0.8, show.legend = FALSE) +
    ggplot2::scale_color_manual(values = c(observed = "#B2182B", simulated = "#2166AC")) +
    ggplot2::labs(
      title = title,
      x = sprintf("PC1 (%.1f%%)", ve[[1]]),
      y = sprintf("PC2 (%.1f%%)", ve[[2]]),
      color = "Data"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold")
    )
}

#' Compute modality-wise real-vs-simulated PCA diagnostics
#'
#' Runs one PCA per omics view after aligning observed and simulated features
#' and standardizing each feature using the observed data. This complements the
#' joint multi-omics PCA by checking whether each view's internal sample
#' structure is preserved.
#'
#' @param observed_matrices Named list of observed sample-by-feature matrices.
#' @param simulated_matrices Named list of simulated sample-by-feature matrices.
#' @param max_features_per_view Optional cap on top-variable features per view.
#'
#' @return A data frame with PCA scores and variance explained for each view.
#' @export
modality_pca_diagnostic <- function(
  observed_matrices,
  simulated_matrices,
  max_features_per_view = NULL
) {
  aligned <- .mimesys_align_observed_simulated(
    observed_matrices,
    simulated_matrices,
    max_features_per_view = max_features_per_view
  )

  rows <- lapply(names(aligned$observed), function(view) {
    mats <- .mimesys_scale_using_observed(aligned$observed[[view]], aligned$simulated[[view]])
    x_obs <- mats$observed
    x_sim <- mats$simulated
    x_all <- rbind(x_obs, x_sim)
    keep <- apply(x_all, 2, function(x) all(is.finite(x)) && stats::var(x) > 0)
    x_all <- x_all[, keep, drop = FALSE]
    if (ncol(x_all) < 2L) {
      return(data.frame(
        view = character(),
        sample_id = character(),
        data_type = character(),
        PC1 = numeric(),
        PC2 = numeric(),
        pc1_var_explained = numeric(),
        pc2_var_explained = numeric(),
        stringsAsFactors = FALSE
      ))
    }

    type <- c(rep("observed", nrow(x_obs)), rep("simulated", nrow(x_sim)))
    sample_id <- c(
      if (is.null(rownames(x_obs))) paste0("observed_", seq_len(nrow(x_obs))) else rownames(x_obs),
      if (is.null(rownames(x_sim))) paste0("simulated_", seq_len(nrow(x_sim))) else rownames(x_sim)
    )

    pca <- stats::prcomp(x_all, center = FALSE, scale. = FALSE)
    var_explained <- 100 * (pca$sdev^2) / sum(pca$sdev^2)
    data.frame(
      view = view,
      sample_id = sample_id,
      data_type = type,
      PC1 = pca$x[, 1],
      PC2 = pca$x[, 2],
      pc1_var_explained = var_explained[[1]],
      pc2_var_explained = var_explained[[2]],
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

#' Plot modality-wise real-vs-simulated PCA diagnostics
#'
#' @param pca_df Data frame from `modality_pca_diagnostic()`.
#' @param title Plot title.
#'
#' @return A `ggplot` object.
#' @export
plot_modality_pca_diagnostic <- function(
  pca_df,
  title = "Observed versus simulated modality-wise PCA"
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot_modality_pca_diagnostic() requires ggplot2.", call. = FALSE)
  }
  needed <- c("view", "data_type", "PC1", "PC2")
  if (!is.data.frame(pca_df) || !all(needed %in% names(pca_df))) {
    stop("pca_df must come from modality_pca_diagnostic().", call. = FALSE)
  }

  ggplot2::ggplot(
    pca_df,
    ggplot2::aes(x = PC1, y = PC2, color = data_type)
  ) +
    ggplot2::geom_point(alpha = 0.6, size = 1.5) +
    ggplot2::stat_ellipse(linewidth = 0.6, alpha = 0.8, show.legend = FALSE) +
    ggplot2::facet_wrap(~ view, scales = "free") +
    ggplot2::scale_color_manual(values = c(observed = "#B2182B", simulated = "#2166AC")) +
    ggplot2::labs(
      title = title,
      x = "PC1",
      y = "PC2",
      color = "Data"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "white", color = "grey70"),
      plot.title = ggplot2::element_text(face = "bold")
    )
}

#' Compare observed and simulated marginal feature distributions
#'
#' @param observed_matrices,simulated_matrices Named lists of sample-by-feature
#'   matrices.
#' @param max_features_per_view Optional cap on top-variable features per view.
#'
#' @return A data frame with feature-level marginal distribution diagnostics.
#' @export
marginal_distribution_diagnostic <- function(
  observed_matrices,
  simulated_matrices,
  max_features_per_view = NULL
) {
  aligned <- .mimesys_align_observed_simulated(
    observed_matrices,
    simulated_matrices,
    max_features_per_view = max_features_per_view
  )

  rows <- list()
  idx <- 1L
  for (view in names(aligned$observed)) {
    obs <- aligned$observed[[view]]
    sim <- aligned$simulated[[view]]
    for (feature in colnames(obs)) {
      obs_v <- obs[, feature]
      sim_v <- sim[, feature]
      obs_v <- obs_v[is.finite(obs_v)]
      sim_v <- sim_v[is.finite(sim_v)]
      if (length(obs_v) < 2L || length(sim_v) < 2L) {
        next
      }
      obs_sd <- stats::sd(obs_v, na.rm = TRUE)
      sim_sd <- stats::sd(sim_v, na.rm = TRUE)
      rows[[idx]] <- data.frame(
        view = view,
        feature = feature,
        ks_stat = suppressWarnings(stats::ks.test(obs_v, sim_v)$statistic[[1]]),
        mean_diff = mean(sim_v, na.rm = TRUE) - mean(obs_v, na.rm = TRUE),
        sd_ratio = if (is.finite(obs_sd) && obs_sd > 0) sim_sd / obs_sd else NA_real_,
        q25_diff = stats::quantile(sim_v, 0.25, na.rm = TRUE, names = FALSE) -
          stats::quantile(obs_v, 0.25, na.rm = TRUE, names = FALSE),
        median_diff = stats::median(sim_v, na.rm = TRUE) - stats::median(obs_v, na.rm = TRUE),
        q75_diff = stats::quantile(sim_v, 0.75, na.rm = TRUE, names = FALSE) -
          stats::quantile(obs_v, 0.75, na.rm = TRUE, names = FALSE),
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }
  do.call(rbind, rows)
}

#' Plot marginal distribution preservation by view
#'
#' @param marginal_df Data frame from `marginal_distribution_diagnostic()`.
#' @param metric Metric to plot.
#' @param title Plot title.
#'
#' @return A `ggplot` object.
#' @export
plot_marginal_distribution_diagnostic <- function(
  marginal_df,
  metric = c("ks_stat", "mean_diff", "median_diff", "sd_ratio"),
  title = "Observed versus simulated marginal distributions"
) {
  metric <- match.arg(metric)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot_marginal_distribution_diagnostic() requires ggplot2.", call. = FALSE)
  }
  if (!is.data.frame(marginal_df) || !all(c("view", metric) %in% names(marginal_df))) {
    stop("marginal_df must come from marginal_distribution_diagnostic().", call. = FALSE)
  }
  plot_df <- marginal_df
  plot_df$value <- plot_df[[metric]]

  ggplot2::ggplot(plot_df, ggplot2::aes(x = view, y = value, fill = view)) +
    ggplot2::geom_boxplot(width = 0.55, outlier.alpha = 0.25) +
    ggplot2::geom_jitter(width = 0.12, alpha = 0.25, size = 0.7) +
    ggplot2::guides(fill = "none") +
    ggplot2::labs(
      title = title,
      x = "View",
      y = metric
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold")
    )
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
#'   current column order; `"cluster"` orders features by hierarchical clustering
#'   on within-view correlation distance.
#'
#' @return A data frame with one row per feature-feature correlation and
#'   columns describing the block membership.
#' @export
blockwise_feature_correlation <- function(
  matrices,
  max_features_per_view = NULL,
  use_abs = FALSE,
  method = c("pearson", "spearman", "kendall"),
  feature_order = c("input", "cluster")
) {
  method <- match.arg(method)
  feature_order <- match.arg(feature_order)
  matrices <- as_mimesys_view_list(matrices)

  matrices <- lapply(matrices, function(x) {
    x <- as.matrix(x)
    storage.mode(x) <- "double"
    x
  })

  if (!is.null(max_features_per_view)) {
    max_features_per_view <- unlist(max_features_per_view)
    if (is.null(names(max_features_per_view)) || any(!nzchar(names(max_features_per_view)))) {
      if (length(max_features_per_view) != 1L) {
        stop("max_features_per_view must be named or a single shared value.", call. = FALSE)
      }
      max_features_per_view <- stats::setNames(
        rep(max_features_per_view[[1]], length(matrices)),
        names(matrices)
      )
    }

    view_names <- names(matrices)
    matrices <- lapply(view_names, function(view) {
      cap <- if (view %in% names(max_features_per_view)) max_features_per_view[[view]] else NA_integer_
      .mimesys_keep_top_variable(matrices[[view]], cap)
    })
    names(matrices) <- view_names
  }

  if (any(vapply(matrices, ncol, integer(1)) == 0L)) {
    stop("Each matrix must contain at least one feature.", call. = FALSE)
  }

  rowname_ok <- all(vapply(matrices, function(x) !is.null(rownames(x)), logical(1)))
  if (rowname_ok) {
    common_samples <- Reduce(intersect, lapply(matrices, rownames))
    if (length(common_samples) < 2L) {
      stop("Matrices must share at least two sample rownames.", call. = FALSE)
    }
    matrices <- lapply(matrices, function(x) x[common_samples, , drop = FALSE])
  } else {
    n_rows <- vapply(matrices, nrow, integer(1))
    if (length(unique(n_rows)) != 1L) {
      stop("Matrices without rownames must have the same number of rows.", call. = FALSE)
    }
  }

  if (identical(feature_order, "cluster")) {
    matrices <- lapply(matrices, .mimesys_cluster_feature_order)
  }

  view_names <- names(matrices)
  feature_view <- unlist(lapply(view_names, function(view) {
    rep(view, ncol(matrices[[view]]))
  }), use.names = FALSE)
  feature_names <- unlist(lapply(view_names, function(view) {
    cols <- colnames(matrices[[view]])
    if (is.null(cols)) {
      cols <- paste0("feature_", seq_len(ncol(matrices[[view]])))
    }
    paste(view, cols, sep = "::")
  }), use.names = FALSE)

  x_all <- do.call(cbind, matrices)
  colnames(x_all) <- feature_names
  cmat <- stats::cor(x_all, use = "pairwise.complete.obs", method = method)
  if (isTRUE(use_abs)) {
    cmat <- abs(cmat)
  }

  df <- as.data.frame(as.table(cmat), stringsAsFactors = FALSE)
  names(df) <- c("feature_x", "feature_y", "correlation")

  feature_lookup <- data.frame(
    feature = feature_names,
    view = feature_view,
    feature_order = seq_along(feature_names),
    stringsAsFactors = FALSE
  )

  match_x <- match(df$feature_x, feature_lookup$feature)
  match_y <- match(df$feature_y, feature_lookup$feature)
  df$view_x <- feature_lookup$view[match_x]
  df$view_y <- feature_lookup$view[match_y]
  df$order_x <- feature_lookup$feature_order[match_x]
  df$order_y <- feature_lookup$feature_order[match_y]
  df$block <- paste(df$view_x, df$view_y, sep = "-")

  attr(df, "use_abs") <- isTRUE(use_abs)
  attr(df, "method") <- method
  attr(df, "feature_order") <- feature_order
  df
}

.mimesys_cluster_feature_order <- function(x) {
  if (ncol(x) < 3L) {
    return(x)
  }
  cmat <- suppressWarnings(stats::cor(x, use = "pairwise.complete.obs"))
  cmat[!is.finite(cmat)] <- 0
  diag(cmat) <- 1
  hc <- stats::hclust(stats::as.dist(1 - abs(cmat)), method = "average")
  x[, hc$order, drop = FALSE]
}

#' Plot blockwise feature-feature correlations across omics views
#'
#' @param corr_df Data frame returned by `blockwise_feature_correlation()`.
#' @param title Plot title.
#' @param use_abs Optional logical controlling the color scale. If `NULL`, this
#'   is read from `corr_df`.
#' @param layout Either `"matrix"` for one full matrix with block boundaries or
#'   `"faceted"` for one facet per view-pair block.
#' @param show_feature_diagonal Logical; if `FALSE`, the feature self-correlation
#'   diagonal is hidden so block structure is easier to see.
#' @param color_limit Optional positive number used as the symmetric color limit
#'   for signed correlations or upper color limit for absolute correlations.
#'
#' @return A `ggplot` object.
#' @export
plot_blockwise_feature_correlation <- function(
  corr_df,
  title = "Blockwise feature correlation",
  use_abs = NULL,
  layout = c("matrix", "faceted"),
  show_feature_diagonal = FALSE,
  color_limit = 0.5
) {
  layout <- match.arg(layout)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot_blockwise_feature_correlation() requires ggplot2.", call. = FALSE)
  }
  needed <- c("order_x", "order_y", "view_x", "view_y", "correlation")
  if (!is.data.frame(corr_df) || !all(needed %in% names(corr_df))) {
    stop("corr_df must come from blockwise_feature_correlation().", call. = FALSE)
  }
  if (is.null(use_abs)) {
    use_abs <- isTRUE(attr(corr_df, "use_abs"))
  }
  plot_df <- corr_df
  if (!isTRUE(show_feature_diagonal)) {
    plot_df$correlation[plot_df$feature_x == plot_df$feature_y] <- NA_real_
  }

  if (is.null(color_limit) || is.na(color_limit)) {
    color_limit <- if (isTRUE(use_abs)) 1 else 1
  }
  color_limit <- abs(as.numeric(color_limit[[1]]))
  if (!is.finite(color_limit) || color_limit <= 0) {
    stop("color_limit must be a positive number or NULL.", call. = FALSE)
  }

  fill_scale <- if (isTRUE(use_abs)) {
    ggplot2::scale_fill_gradient(
      low = "grey95",
      high = "#2166AC",
      limits = c(0, color_limit),
      na.value = "white"
    )
  } else {
    ggplot2::scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-color_limit, color_limit),
      na.value = "white"
    )
  }

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = order_x,
      y = order_y,
      fill = correlation
    )
  ) +
    ggplot2::geom_tile() +
    fill_scale +
    ggplot2::labs(
      title = title,
      x = "Features",
      y = "Features",
      fill = ifelse(isTRUE(use_abs), "|cor|", "cor")
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "white", color = "grey70"),
      plot.title = ggplot2::element_text(face = "bold")
    )

  if (identical(layout, "faceted")) {
    return(
      p +
        ggplot2::facet_grid(view_y ~ view_x, scales = "free", space = "free")
    )
  }

  block_limits <- stats::aggregate(
    order_x ~ view_x,
    data = unique(plot_df[, c("view_x", "order_x")]),
    FUN = max
  )
  block_limits <- block_limits[order(block_limits$order_x), , drop = FALSE]
  boundary <- block_limits$order_x[-nrow(block_limits)] + 0.5
  label_pos <- stats::aggregate(
    order_x ~ view_x,
    data = unique(plot_df[, c("view_x", "order_x")]),
    FUN = function(x) mean(range(x))
  )

  p +
    ggplot2::geom_vline(xintercept = boundary, color = "black", linewidth = 0.25) +
    ggplot2::geom_hline(yintercept = boundary, color = "black", linewidth = 0.25) +
    ggplot2::scale_x_continuous(
      breaks = label_pos$order_x,
      labels = label_pos$view_x,
      expand = c(0, 0)
    ) +
    ggplot2::scale_y_continuous(
      breaks = label_pos$order_x,
      labels = label_pos$view_x,
      expand = c(0, 0)
    ) +
    ggplot2::theme(
      axis.text = ggplot2::element_text(color = "black"),
      axis.ticks = ggplot2::element_blank()
    )
}

#' Summarize blockwise feature correlations
#'
#' @param corr_df Data frame returned by `blockwise_feature_correlation()`.
#' @param unique_pairs Logical; if `TRUE`, symmetric off-diagonal blocks are
#'   collapsed so each view pair appears once.
#'
#' @return A data frame with one row per view-pair block.
#' @export
summarize_blockwise_feature_correlation <- function(corr_df, unique_pairs = TRUE) {
  needed <- c("view_x", "view_y", "feature_x", "feature_y", "correlation")
  if (!is.data.frame(corr_df) || !all(needed %in% names(corr_df))) {
    stop("corr_df must come from blockwise_feature_correlation().", call. = FALSE)
  }

  corr_df$diagonal_block <- corr_df$view_x == corr_df$view_y
  if (isTRUE(unique_pairs)) {
    pair_a <- pmin(corr_df$view_x, corr_df$view_y)
    pair_b <- pmax(corr_df$view_x, corr_df$view_y)
    corr_df$view_x <- pair_a
    corr_df$view_y <- pair_b
  }
  nonself <- corr_df$feature_x != corr_df$feature_y
  dat <- corr_df[nonself, , drop = FALSE]

  out <- stats::aggregate(
    cbind(
      mean_cor = dat$correlation,
      mean_abs_cor = abs(dat$correlation),
      median_abs_cor = abs(dat$correlation)
    ) ~ view_x + view_y + diagonal_block,
    data = dat,
    FUN = function(x) c(mean = mean(x, na.rm = TRUE), median = stats::median(x, na.rm = TRUE))
  )

  data.frame(
    view_x = out$view_x,
    view_y = out$view_y,
    diagonal_block = out$diagonal_block,
    mean_cor = out$mean_cor[, "mean"],
    median_cor = out$mean_cor[, "median"],
    mean_abs_cor = out$mean_abs_cor[, "mean"],
    median_abs_cor = out$median_abs_cor[, "median"],
    stringsAsFactors = FALSE
  )
}

#' Plot a compact heatmap of block-level correlations
#'
#' @param block_summary Data frame returned by
#'   `summarize_blockwise_feature_correlation()`.
#' @param value Column to plot.
#' @param title Plot title.
#'
#' @return A `ggplot` object.
#' @export
plot_blockwise_correlation_summary <- function(
  block_summary,
  value = c("mean_abs_cor", "median_abs_cor", "mean_cor", "median_cor"),
  title = "Block-level feature correlation summary"
) {
  value <- match.arg(value)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot_blockwise_correlation_summary() requires ggplot2.", call. = FALSE)
  }
  needed <- c("view_x", "view_y", value)
  if (!is.data.frame(block_summary) || !all(needed %in% names(block_summary))) {
    stop("block_summary must come from summarize_blockwise_feature_correlation().", call. = FALSE)
  }
  plot_df <- block_summary
  plot_df$plot_value <- plot_df[[value]]

  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = view_x,
      y = view_y,
      fill = plot_value
    )
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.3f", plot_value)),
      size = 3.5
    ) +
    ggplot2::scale_fill_gradient(low = "grey95", high = "#2166AC") +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = title,
      x = "View",
      y = "View",
      fill = value
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "white", color = "grey70"),
      plot.title = ggplot2::element_text(face = "bold")
    )
}

#' Compare observed and simulated block-level correlations
#'
#' @param observed_summary,simulated_summary Data frames from
#'   `summarize_blockwise_feature_correlation()`.
#' @param value Column to compare.
#'
#' @return A data frame with observed, simulated, and difference values.
#' @export
compare_blockwise_correlation_summary <- function(
  observed_summary,
  simulated_summary,
  value = c("mean_abs_cor", "median_abs_cor", "mean_cor", "median_cor")
) {
  value <- match.arg(value)
  needed <- c("view_x", "view_y", "diagonal_block", value)
  if (!is.data.frame(observed_summary) || !all(needed %in% names(observed_summary))) {
    stop("observed_summary must come from summarize_blockwise_feature_correlation().", call. = FALSE)
  }
  if (!is.data.frame(simulated_summary) || !all(needed %in% names(simulated_summary))) {
    stop("simulated_summary must come from summarize_blockwise_feature_correlation().", call. = FALSE)
  }

  obs <- observed_summary[, needed, drop = FALSE]
  sim <- simulated_summary[, needed, drop = FALSE]
  names(obs)[names(obs) == value] <- "observed"
  names(sim)[names(sim) == value] <- "simulated"
  out <- merge(obs, sim, by = c("view_x", "view_y", "diagonal_block"), all = TRUE)
  out$difference <- out$simulated - out$observed
  out$value <- value
  out
}

#' Compare observed and simulated feature-correlation structure
#'
#' Correlates the vector of observed feature-feature correlations with the
#' matching vector of simulated feature-feature correlations within each view
#' block. This asks whether the simulated data preserve the pattern of the
#' correlation matrix, not only its average magnitude.
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
#' @export
compare_feature_correlation_structure <- function(observed_corr, simulated_corr, return_pairs = FALSE) {
  needed <- c("feature_x", "feature_y", "view_x", "view_y", "correlation")
  if (!is.data.frame(observed_corr) || !all(needed %in% names(observed_corr))) {
    stop("observed_corr must come from blockwise_feature_correlation().", call. = FALSE)
  }
  if (!is.data.frame(simulated_corr) || !all(needed %in% names(simulated_corr))) {
    stop("simulated_corr must come from blockwise_feature_correlation().", call. = FALSE)
  }

  obs <- observed_corr[, needed, drop = FALSE]
  sim <- simulated_corr[, needed, drop = FALSE]
  names(obs)[names(obs) == "correlation"] <- "observed_cor"
  names(sim)[names(sim) == "correlation"] <- "simulated_cor"
  dat <- merge(obs, sim, by = c("feature_x", "feature_y", "view_x", "view_y"))
  dat <- dat[dat$feature_x != dat$feature_y, , drop = FALSE]
  dat <- dat[is.finite(dat$observed_cor) & is.finite(dat$simulated_cor), , drop = FALSE]
  if (nrow(dat) == 0L) {
    empty_summary <- data.frame(
      view_x = character(),
      view_y = character(),
      diagonal_block = logical(),
      cor_pearson = numeric(),
      cor_spearman = numeric(),
      mae_cor = numeric(),
      rmse_cor = numeric(),
      n_pairs = integer(),
      stringsAsFactors = FALSE
    )
    if (isTRUE(return_pairs)) {
      return(list(summary = empty_summary, pairs = dat))
    }
    return(empty_summary)
  }

  keys <- split(seq_len(nrow(dat)), paste(dat$view_x, dat$view_y, sep = "\r"))
  rows <- lapply(keys, function(idx) {
    d <- dat[idx, , drop = FALSE]
    data.frame(
      view_x = d$view_x[[1]],
      view_y = d$view_y[[1]],
      diagonal_block = d$view_x[[1]] == d$view_y[[1]],
      cor_pearson = if (nrow(d) >= 2L) {
        suppressWarnings(stats::cor(d$observed_cor, d$simulated_cor, method = "pearson"))
      } else {
        NA_real_
      },
      cor_spearman = if (nrow(d) >= 2L) {
        suppressWarnings(stats::cor(d$observed_cor, d$simulated_cor, method = "spearman"))
      } else {
        NA_real_
      },
      mae_cor = mean(abs(d$simulated_cor - d$observed_cor), na.rm = TRUE),
      rmse_cor = sqrt(mean((d$simulated_cor - d$observed_cor)^2, na.rm = TRUE)),
      n_pairs = nrow(d),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  if (isTRUE(return_pairs)) {
    rownames(dat) <- NULL
    return(list(summary = out, pairs = dat))
  }
  out
}

#' Plot observed-vs-simulated block correlation differences
#'
#' @param comparison_df Data frame from `compare_blockwise_correlation_summary()`.
#' @param title Plot title.
#' @param difference_limit Symmetric color limit for simulated-minus-observed
#'   differences.
#'
#' @return A `ggplot` object.
#' @export
plot_blockwise_correlation_difference <- function(
  comparison_df,
  title = "Simulated minus observed block-level correlation",
  difference_limit = 0.05
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot_blockwise_correlation_difference() requires ggplot2.", call. = FALSE)
  }
  needed <- c("view_x", "view_y", "difference")
  if (!is.data.frame(comparison_df) || !all(needed %in% names(comparison_df))) {
    stop("comparison_df must come from compare_blockwise_correlation_summary().", call. = FALSE)
  }

  ggplot2::ggplot(
    comparison_df,
    ggplot2::aes(x = view_x, y = view_y, fill = difference)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%+.3f", difference)),
      size = 3.5
    ) +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-difference_limit, difference_limit),
      na.value = "grey90"
    ) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = title,
      x = "View",
      y = "View",
      fill = "Sim - Obs"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold")
    )
}

#' Plot crosstalk-gradient diagnostics
#'
#' @param diagnostics Data frame returned by `validate_template_simulation()`.
#' @param metrics Diagnostic columns to display.
#' @param title Plot title.
#'
#' @return A `ggplot` object.
#' @export
plot_crosstalk_gradient_diagnostic <- function(
  diagnostics,
  metrics = c("mean_abs_feature_cor", "signal_agreement"),
  title = "Crosstalk-gradient validation"
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot_crosstalk_gradient_diagnostic() requires ggplot2.", call. = FALSE)
  }
  needed <- c("view_a", "view_b", "crosstalk", metrics)
  if (!is.data.frame(diagnostics) || !all(needed %in% names(diagnostics))) {
    missing <- paste(setdiff(needed, names(diagnostics)), collapse = ", ")
    stop(sprintf("diagnostics is missing required columns: %s", missing), call. = FALSE)
  }

  diag <- diagnostics
  diag$view_pair <- paste(diag$view_a, diag$view_b, sep = "-")
  diag$crosstalk_numeric <- suppressWarnings(as.numeric(diag$crosstalk))
  diag$effect_alignment <- .mimesys_effect_alignment_label(diag)

  grouping <- c("crosstalk_numeric", "effect_alignment", "view_pair")
  key <- interaction(diag[, grouping, drop = FALSE], drop = TRUE, sep = "\r")
  split_idx <- split(seq_len(nrow(diag)), key)
  summary <- do.call(rbind, lapply(split_idx, function(idx) {
    first <- diag[idx[1], grouping, drop = FALSE]
    vals <- lapply(metrics, function(metric) {
      x <- diag[[metric]][idx]
      if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
    })
    names(vals) <- metrics
    data.frame(first, vals, stringsAsFactors = FALSE)
  }))
  rownames(summary) <- NULL

  metric_labels <- c(
    mean_abs_feature_cor = "Mean absolute feature correlation",
    signal_agreement = "Signal agreement",
    y_signal_cor = "Outcome-signal correlation"
  )
  long <- do.call(rbind, lapply(metrics, function(metric) {
    data.frame(
      crosstalk = summary$crosstalk_numeric,
      effect_alignment = summary$effect_alignment,
      view_pair = summary$view_pair,
      metric = if (metric %in% names(metric_labels)) metric_labels[[metric]] else metric,
      value = summary[[metric]],
      stringsAsFactors = FALSE
    )
  }))

  ggplot2::ggplot(
    long,
    ggplot2::aes(x = crosstalk, y = value, color = view_pair)
  ) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 2.2) +
    ggplot2::facet_grid(metric ~ effect_alignment, scales = "free_y") +
    ggplot2::labs(
      title = title,
      x = "Crosstalk",
      y = "Diagnostic value",
      color = "View pair"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "white", color = "grey70"),
      plot.title = ggplot2::element_text(face = "bold")
    )
}

.mimesys_template_diagnostic_plot_list <- function(diagnostics) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Template diagnostic plotting requires ggplot2.", call. = FALSE)
  }
  needed <- c(
    "rep", "view_a", "view_b", "crosstalk", "mean_abs_feature_cor",
    "signal_agreement", "y_signal_cor"
  )
  if (!is.data.frame(diagnostics) || !all(needed %in% names(diagnostics))) {
    missing <- paste(setdiff(needed, names(diagnostics)), collapse = ", ")
    stop(sprintf("diagnostics is missing required columns: %s", missing), call. = FALSE)
  }

  diag <- diagnostics
  diag$view_pair <- paste(diag$view_a, diag$view_b, sep = "-")
  diag$crosstalk_numeric <- suppressWarnings(as.numeric(diag$crosstalk))
  diag$crosstalk_factor <- factor(
    diag$crosstalk_numeric,
    levels = sort(unique(diag$crosstalk_numeric))
  )
  diag$effect_alignment <- .mimesys_effect_alignment_label(diag)

  summary <- summarize_template_diagnostics(diag)
  summary$view_pair <- paste(summary$view_a, summary$view_b, sep = "-")
  summary$crosstalk_numeric <- suppressWarnings(as.numeric(summary$crosstalk))
  summary$crosstalk_factor <- factor(
    summary$crosstalk_numeric,
    levels = sort(unique(summary$crosstalk_numeric))
  )
  summary$effect_alignment <- .mimesys_effect_alignment_label(summary)

  metric_labels <- c(
    mean_abs_feature_cor = "Mean absolute\nfeature correlation",
    signal_agreement = "Signal agreement",
    y_signal_cor = "Outcome-signal\ncorrelation"
  )

  make_long <- function(dat) {
    do.call(rbind, lapply(names(metric_labels), function(metric) {
      data.frame(
        rep = if ("rep" %in% names(dat)) dat$rep else NA_integer_,
        view_pair = dat$view_pair,
        crosstalk_numeric = dat$crosstalk_numeric,
        crosstalk_factor = dat$crosstalk_factor,
        effect_alignment = dat$effect_alignment,
        metric = metric_labels[[metric]],
        value = dat[[metric]],
        stringsAsFactors = FALSE
      )
    }))
  }

  summary_long <- make_long(summary)
  diag_long <- make_long(diag)
  metric_order <- unname(metric_labels)
  summary_long$metric <- factor(summary_long$metric, levels = metric_order)
  diag_long$metric <- factor(diag_long$metric, levels = metric_order)

  base_theme <- ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "white", color = "grey70"),
      strip.text = ggplot2::element_text(color = "black"),
      legend.position = "right",
      plot.title.position = "plot",
      plot.title = ggplot2::element_text(face = "bold")
    )

  line_plot <- ggplot2::ggplot(
    summary_long,
    ggplot2::aes(
      x = crosstalk_numeric,
      y = value,
      color = view_pair,
      group = interaction(metric, view_pair, effect_alignment)
    )
  ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2.2) +
    ggplot2::facet_grid(metric ~ effect_alignment, scales = "free_y") +
    ggplot2::labs(
      title = "Template simulation diagnostics across crosstalk",
      x = "Crosstalk",
      y = "Diagnostic value",
      color = "View pair"
    ) +
    base_theme

  heatmap_plot <- ggplot2::ggplot(
    summary,
    ggplot2::aes(
      x = crosstalk_factor,
      y = view_pair,
      fill = mean_abs_feature_cor
    )
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.4) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.3f", mean_abs_feature_cor)),
      size = 3
    ) +
    ggplot2::facet_wrap(~ effect_alignment) +
    ggplot2::scale_fill_gradient(low = "grey95", high = "#2166AC") +
    ggplot2::labs(
      title = "Cross-view feature correlation",
      x = "Crosstalk",
      y = "View pair",
      fill = "Mean |cor|"
    ) +
    base_theme +
    ggplot2::theme(panel.grid = ggplot2::element_blank())

  variability_plot <- ggplot2::ggplot(
    diag,
    ggplot2::aes(
      x = crosstalk_factor,
      y = signal_agreement
    )
  ) +
    ggplot2::geom_boxplot(
      ggplot2::aes(fill = effect_alignment),
      position = ggplot2::position_dodge(width = 0.75),
      width = 0.55,
      outlier.shape = NA,
      alpha = 0.25
    ) +
    ggplot2::geom_point(
      ggplot2::aes(color = effect_alignment),
      position = ggplot2::position_jitterdodge(jitter.width = 0.12, dodge.width = 0.75),
      alpha = 0.45,
      size = 1.2
    ) +
    ggplot2::facet_wrap(~ view_pair, nrow = 1) +
    ggplot2::labs(
      title = "Replicate-level variability in signal agreement",
      x = "Crosstalk",
      y = "Signal agreement",
      fill = "Effect setting",
      color = "Effect setting"
    ) +
    base_theme

  scatter_plot <- ggplot2::ggplot(
    diag,
    ggplot2::aes(
      x = mean_abs_feature_cor,
      y = signal_agreement,
      color = crosstalk_numeric
    )
  ) +
    ggplot2::geom_point(alpha = 0.75, size = 2.2) +
    ggplot2::facet_grid(effect_alignment ~ view_pair) +
    ggplot2::scale_color_gradient(low = "grey35", high = "#2166AC") +
    ggplot2::labs(
      title = "Feature correlation versus signal agreement",
      x = "Mean absolute cross-view feature correlation",
      y = "Signal agreement",
      color = "Crosstalk"
    ) +
    base_theme

  list(
    original_crosstalk_line = line_plot,
    original_feature_correlation_heatmap = heatmap_plot,
    original_replicate_variability = variability_plot,
    original_featurecor_signal_scatter = scatter_plot
  )
}

#' Save the full MimeSys validation figure set
#'
#' Saves the original crosstalk diagnostic plots plus added real-vs-simulated
#' validation figures as individual PDFs and as one multi-page PDF.
#'
#' @param diagnostics Data frame returned by `validate_template_simulation()`
#'   across crosstalk settings and replicates.
#' @param observed_matrices Named list of observed sample-by-feature matrices.
#' @param simulated_matrices Named list of simulated sample-by-feature matrices.
#' @param plot_dir Directory where PDFs should be saved.
#' @param prefix File prefix for saved PDFs.
#' @param max_features_per_view Optional feature cap used for PCA, marginal,
#'   and block-level diagnostics.
#' @param width,height PDF dimensions in inches.
#'
#' @return A named character vector of saved PDF paths.
#' @export
save_full_validation_figure_set <- function(
  diagnostics,
  observed_matrices,
  simulated_matrices,
  plot_dir,
  prefix = "mimesys_validation",
  max_features_per_view = NULL,
  width = 10,
  height = 7
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("save_full_validation_figure_set() requires ggplot2.", call. = FALSE)
  }
  if (missing(plot_dir) || length(plot_dir) != 1L || is.na(plot_dir) || !nzchar(plot_dir)) {
    stop("plot_dir must be a non-empty directory path supplied by the user.", call. = FALSE)
  }
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(plot_dir)) {
    stop(sprintf("Could not create plot_dir: %s", plot_dir), call. = FALSE)
  }

  original_plots <- .mimesys_template_diagnostic_plot_list(diagnostics)

  pca_diag <- multiomics_pca_diagnostic(
    observed_matrices = observed_matrices,
    simulated_matrices = simulated_matrices,
    max_features_per_view = max_features_per_view,
    equal_view_weight = TRUE
  )
  pca_plot <- plot_multiomics_pca_diagnostic(
    pca_diag,
    title = "Observed versus simulated multi-omics PCA"
  )

  marginal_diag <- marginal_distribution_diagnostic(
    observed_matrices = observed_matrices,
    simulated_matrices = simulated_matrices,
    max_features_per_view = max_features_per_view
  )
  marginal_plot <- plot_marginal_distribution_diagnostic(
    marginal_diag,
    metric = "ks_stat",
    title = "Observed versus simulated marginal distributions"
  )

  obs_corr <- blockwise_feature_correlation(
    observed_matrices,
    max_features_per_view = max_features_per_view
  )
  sim_corr <- blockwise_feature_correlation(
    simulated_matrices,
    max_features_per_view = max_features_per_view
  )
  obs_block_summary <- summarize_blockwise_feature_correlation(obs_corr)
  sim_block_summary <- summarize_blockwise_feature_correlation(sim_corr)

  observed_block_plot <- plot_blockwise_correlation_summary(
    obs_block_summary,
    value = "mean_cor",
    title = "Observed block-level mean signed correlation"
  )
  simulated_block_plot <- plot_blockwise_correlation_summary(
    sim_block_summary,
    value = "mean_cor",
    title = "Simulated block-level mean signed correlation"
  )
  block_comparison <- compare_blockwise_correlation_summary(
    obs_block_summary,
    sim_block_summary,
    value = "mean_cor"
  )
  block_difference_plot <- plot_blockwise_correlation_difference(
    block_comparison,
    title = "Simulated minus observed block-level mean signed correlation"
  )

  observed_abs_block_plot <- plot_blockwise_correlation_summary(
    obs_block_summary,
    value = "mean_abs_cor",
    title = "Observed block-level mean absolute correlation"
  )
  simulated_abs_block_plot <- plot_blockwise_correlation_summary(
    sim_block_summary,
    value = "mean_abs_cor",
    title = "Simulated block-level mean absolute correlation"
  )
  block_abs_comparison <- compare_blockwise_correlation_summary(
    obs_block_summary,
    sim_block_summary,
    value = "mean_abs_cor"
  )
  block_abs_difference_plot <- plot_blockwise_correlation_difference(
    block_abs_comparison,
    title = "Simulated minus observed block-level mean absolute correlation"
  )

  crosstalk_gradient_plot <- plot_crosstalk_gradient_diagnostic(
    diagnostics = diagnostics,
    metrics = c("mean_abs_feature_cor", "signal_agreement"),
    title = "Controlled crosstalk-gradient validation"
  )

  plots <- c(
    original_plots,
    list(
      pca = pca_plot,
      marginal_ks = marginal_plot,
      observed_block_signed_summary = observed_block_plot,
      simulated_block_signed_summary = simulated_block_plot,
      block_signed_difference = block_difference_plot,
      observed_block_absolute_summary = observed_abs_block_plot,
      simulated_block_absolute_summary = simulated_abs_block_plot,
      block_absolute_difference = block_abs_difference_plot,
      crosstalk_gradient = crosstalk_gradient_plot
    )
  )

  paths <- file.path(plot_dir, paste0(prefix, "_", names(plots), ".pdf"))
  for (i in seq_along(plots)) {
    ggplot2::ggsave(
      filename = paths[[i]],
      plot = plots[[i]],
      width = width,
      height = height,
      units = "in",
      device = grDevices::pdf
    )
  }

  master_path <- file.path(plot_dir, paste0(prefix, "_full_validation_figure_set.pdf"))
  grDevices::pdf(master_path, width = width, height = height, onefile = TRUE)
  tryCatch(
    {
      for (plot_name in names(plots)) {
        print(plots[[plot_name]])
      }
    },
    finally = grDevices::dev.off()
  )

  out <- stats::setNames(paths, names(plots))
  c(out, full_validation_figure_set = master_path)
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
#' @export
save_blockwise_feature_correlation_plot <- function(
  matrices,
  plot_dir,
  prefix = "blockwise_feature_correlation",
  title = "Blockwise feature correlation",
  max_features_per_view = NULL,
  use_abs = FALSE,
  method = c("pearson", "spearman", "kendall"),
  layout = c("matrix", "faceted"),
  show_feature_diagonal = FALSE,
  color_limit = 0.5,
  width = 10,
  height = 8
) {
  layout <- match.arg(layout)
  if (missing(plot_dir) || length(plot_dir) != 1L || is.na(plot_dir) || !nzchar(plot_dir)) {
    stop("plot_dir must be a non-empty directory path supplied by the user.", call. = FALSE)
  }
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(plot_dir)) {
    stop(sprintf("Could not create plot_dir: %s", plot_dir), call. = FALSE)
  }

  corr_df <- blockwise_feature_correlation(
    matrices = matrices,
    max_features_per_view = max_features_per_view,
    use_abs = use_abs,
    method = method
  )
  plot <- plot_blockwise_feature_correlation(
    corr_df,
    title = title,
    use_abs = use_abs,
    layout = layout,
    show_feature_diagonal = show_feature_diagonal,
    color_limit = color_limit
  )

  path <- file.path(plot_dir, paste0(prefix, ".pdf"))
  ggplot2::ggsave(
    filename = path,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = grDevices::pdf
  )
  path
}

#' Plot template simulation diagnostics
#'
#' Visualizes whether cross-view feature correlation, aligned signal, and
#' outcome signal behave as expected across crosstalk settings.
#'
#' @param diagnostics A data frame returned by `validate_template_simulation()`
#'   or `summarize_template_diagnostics()`.
#' @param metrics Character vector of diagnostic columns to plot.
#'
#' @return A `ggplot` object.
#' @export
plot_template_diagnostics <- function(
  diagnostics,
  metrics = c("mean_abs_feature_cor", "signal_agreement", "y_signal_cor")
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot_template_diagnostics() requires ggplot2.", call. = FALSE)
  }
  if (!is.data.frame(diagnostics)) {
    stop("diagnostics must be a data frame.", call. = FALSE)
  }

  needed <- c("view_a", "view_b", metrics)
  if (!all(needed %in% names(diagnostics))) {
    missing <- paste(setdiff(needed, names(diagnostics)), collapse = ", ")
    stop(sprintf("diagnostics is missing required columns: %s", missing), call. = FALSE)
  }

  if (!"crosstalk" %in% names(diagnostics)) {
    diagnostics$crosstalk <- NA_real_
  }
  diagnostics$effect_alignment <- .mimesys_effect_alignment_label(diagnostics)

  long <- do.call(rbind, lapply(metrics, function(metric) {
    data.frame(
      view_pair = paste(diagnostics$view_a, diagnostics$view_b, sep = "-"),
      crosstalk = suppressWarnings(as.numeric(diagnostics$crosstalk)),
      de_align = diagnostics$effect_alignment,
      metric = metric,
      value = diagnostics[[metric]],
      stringsAsFactors = FALSE
    )
  }))

  metric_labels <- c(
    mean_abs_feature_cor = "Mean absolute\nfeature correlation",
    signal_agreement = "Signal agreement",
    y_signal_cor = "Outcome-signal\ncorrelation"
  )
  long$metric <- factor(
    long$metric,
    levels = metrics,
    labels = ifelse(metrics %in% names(metric_labels), metric_labels[metrics], metrics)
  )

  ggplot2::ggplot(
    long,
    ggplot2::aes(
      x = crosstalk,
      y = value,
      color = view_pair,
      group = interaction(metric, view_pair, de_align)
    )
  ) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::facet_grid(metric ~ de_align, scales = "free_y") +
    ggplot2::labs(
      x = "Crosstalk",
      y = "Diagnostic value",
      color = "View pair"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "white", color = "grey70")
    )
}

#' Save publication-style template diagnostic plots
#'
#' Creates a standard plot set for validating template-based simulations:
#' a crosstalk diagnostic line plot, a feature-correlation heatmap, a replicate
#' variability plot, and a feature-correlation versus signal-agreement scatter
#' plot.
#'
#' @param diagnostics A data frame returned by `validate_template_simulation()`.
#'   This should usually contain multiple crosstalk settings and simulation
#'   replicates.
#' @param plot_dir Directory where PDF files should be saved.
#' @param prefix File prefix used for the saved PDFs.
#' @param width,height PDF dimensions in inches.
#' @param combined_pdf Logical; if `TRUE`, also save one multi-page PDF
#'   containing all diagnostic plots.
#'
#' @return A named character vector containing the saved PDF paths.
#' @export
save_template_diagnostic_plots <- function(
  diagnostics,
  plot_dir,
  prefix = "template_diagnostics",
  width = 10,
  height = 7,
  combined_pdf = TRUE
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("save_template_diagnostic_plots() requires ggplot2.", call. = FALSE)
  }
  if (missing(plot_dir) || length(plot_dir) != 1L || is.na(plot_dir) || !nzchar(plot_dir)) {
    stop("plot_dir must be a non-empty directory path supplied by the user.", call. = FALSE)
  }
  if (!is.data.frame(diagnostics)) {
    stop("diagnostics must be a data frame.", call. = FALSE)
  }

  needed <- c(
    "rep", "view_a", "view_b", "crosstalk", "mean_abs_feature_cor",
    "signal_agreement", "y_signal_cor"
  )
  if (!all(needed %in% names(diagnostics))) {
    missing <- paste(setdiff(needed, names(diagnostics)), collapse = ", ")
    stop(sprintf("diagnostics is missing required columns: %s", missing), call. = FALSE)
  }

  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(plot_dir)) {
    stop(sprintf("Could not create plot_dir: %s", plot_dir), call. = FALSE)
  }

  diag <- diagnostics
  diag$view_pair <- paste(diag$view_a, diag$view_b, sep = "-")
  diag$crosstalk_numeric <- suppressWarnings(as.numeric(diag$crosstalk))
  diag$crosstalk_factor <- factor(
    diag$crosstalk_numeric,
    levels = sort(unique(diag$crosstalk_numeric))
  )
  diag$effect_alignment <- .mimesys_effect_alignment_label(diag)

  summary <- summarize_template_diagnostics(diag)
  summary$view_pair <- paste(summary$view_a, summary$view_b, sep = "-")
  summary$crosstalk_numeric <- suppressWarnings(as.numeric(summary$crosstalk))
  summary$crosstalk_factor <- factor(
    summary$crosstalk_numeric,
    levels = sort(unique(summary$crosstalk_numeric))
  )
  summary$effect_alignment <- .mimesys_effect_alignment_label(summary)

  metric_labels <- c(
    mean_abs_feature_cor = "Mean absolute\nfeature correlation",
    signal_agreement = "Signal agreement",
    y_signal_cor = "Outcome-signal\ncorrelation"
  )

  make_long <- function(dat) {
    do.call(rbind, lapply(names(metric_labels), function(metric) {
      data.frame(
        rep = if ("rep" %in% names(dat)) dat$rep else NA_integer_,
        view_pair = dat$view_pair,
        crosstalk_numeric = dat$crosstalk_numeric,
        crosstalk_factor = dat$crosstalk_factor,
        effect_alignment = dat$effect_alignment,
        metric = metric_labels[[metric]],
        value = dat[[metric]],
        stringsAsFactors = FALSE
      )
    }))
  }

  summary_long <- make_long(summary)
  diag_long <- make_long(diag)
  metric_order <- unname(metric_labels)
  summary_long$metric <- factor(summary_long$metric, levels = metric_order)
  diag_long$metric <- factor(diag_long$metric, levels = metric_order)

  base_theme <- ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "white", color = "grey70"),
      strip.text = ggplot2::element_text(color = "black"),
      legend.position = "right",
      plot.title.position = "plot",
      plot.title = ggplot2::element_text(face = "bold")
    )

  line_plot <- ggplot2::ggplot(
    summary_long,
    ggplot2::aes(
      x = crosstalk_numeric,
      y = value,
      color = view_pair,
      group = interaction(metric, view_pair, effect_alignment)
    )
  ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2.2) +
    ggplot2::facet_grid(metric ~ effect_alignment, scales = "free_y") +
    ggplot2::labs(
      title = "Template simulation diagnostics across crosstalk",
      x = "Crosstalk",
      y = "Diagnostic value",
      color = "View pair"
    ) +
    base_theme

  heatmap_plot <- ggplot2::ggplot(
    summary,
    ggplot2::aes(
      x = crosstalk_factor,
      y = view_pair,
      fill = mean_abs_feature_cor
    )
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.4) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.3f", mean_abs_feature_cor)),
      size = 3
    ) +
    ggplot2::facet_wrap(~ effect_alignment) +
    ggplot2::scale_fill_gradient(low = "grey95", high = "#2166AC") +
    ggplot2::labs(
      title = "Cross-view feature correlation",
      x = "Crosstalk",
      y = "View pair",
      fill = "Mean |cor|"
    ) +
    base_theme +
    ggplot2::theme(panel.grid = ggplot2::element_blank())

  variability_plot <- ggplot2::ggplot(
    diag,
    ggplot2::aes(
      x = crosstalk_factor,
      y = signal_agreement
    )
  ) +
    ggplot2::geom_boxplot(
      ggplot2::aes(fill = effect_alignment),
      position = ggplot2::position_dodge(width = 0.75),
      width = 0.55,
      outlier.shape = NA,
      alpha = 0.25
    ) +
    ggplot2::geom_point(
      ggplot2::aes(color = effect_alignment),
      position = ggplot2::position_jitterdodge(jitter.width = 0.12, dodge.width = 0.75),
      alpha = 0.45,
      size = 1.2
    ) +
    ggplot2::facet_wrap(~ view_pair, nrow = 1) +
    ggplot2::labs(
      title = "Replicate-level variability in signal agreement",
      x = "Crosstalk",
      y = "Signal agreement",
      fill = "Effect setting",
      color = "Effect setting"
    ) +
    base_theme

  scatter_plot <- ggplot2::ggplot(
    diag,
    ggplot2::aes(
      x = mean_abs_feature_cor,
      y = signal_agreement,
      color = crosstalk_numeric
    )
  ) +
    ggplot2::geom_point(
      alpha = 0.75,
      size = 2.2
    ) +
    ggplot2::facet_grid(effect_alignment ~ view_pair) +
    ggplot2::scale_color_gradient(low = "grey35", high = "#2166AC") +
    ggplot2::labs(
      title = "Feature correlation versus signal agreement",
      x = "Mean absolute cross-view feature correlation",
      y = "Signal agreement",
      color = "Crosstalk"
    ) +
    base_theme

  plots <- list(
    crosstalk_line = line_plot,
    feature_correlation_heatmap = heatmap_plot,
    replicate_variability = variability_plot,
    featurecor_signal_scatter = scatter_plot
  )
  paths <- file.path(plot_dir, paste0(prefix, "_", names(plots), ".pdf"))

  for (i in seq_along(plots)) {
    ggplot2::ggsave(
      filename = paths[[i]],
      plot = plots[[i]],
      width = width,
      height = height,
      units = "in",
      device = grDevices::pdf
    )
  }

  out <- stats::setNames(paths, names(plots))

  if (isTRUE(combined_pdf)) {
    combined_path <- file.path(plot_dir, paste0(prefix, "_all_plots.pdf"))
    grDevices::pdf(combined_path, width = width, height = height, onefile = TRUE)
    tryCatch(
      {
        for (plot_name in names(plots)) {
          print(plots[[plot_name]])
        }
      },
      finally = grDevices::dev.off()
    )
    out <- c(out, all_plots = combined_path)
  }

  out
}
