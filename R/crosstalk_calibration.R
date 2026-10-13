#' Calibrate crosstalk to empirical cross-view correlation structure
#'
#' Sweeps candidate crosstalk values, simulates null data from the supplied
#' template, and chooses the value that minimizes the rescaled Frobenius distance
#' between empirical and simulated cross-view feature-correlation blocks.
#'
#' @param template A fitted \code{mimesys_template} or any object accepted by
#'   \code{fit_template()}.
#' @param crosstalk_grid Numeric vector of crosstalk values to evaluate.
#' @param nsample Number of simulated samples per replicate.
#' @param nrep Number of replicates per crosstalk value.
#' @param p Optional named integer vector of simulated features per view.
#' @param max_features_per_view Optional named cap used when computing
#'   correlation distances.
#' @param dependence Dependence simulator passed to \code{gen_simmba()}.
#' @param cor.shrinkage Shrinkage passed to \code{gen_simmba()} when
#'   \code{dependence = "empirical_copula"} or `"module_copula"`.
#' @param module.count,module.between.shrinkage Module-copula parameters passed
#'   to \code{gen_simmba()}.
#' @param view.models Optional named vector of per-view observation models passed
#'   to \code{gen_simmba()}.
#' @param seed Random seed.
#' @param ... Passed to \code{fit_template()} when \code{template} is raw data.
#'
#' @return A list with replicate-level distances, summarized distances, and the
#'   best crosstalk value.
#' @export
calibrate_crosstalk <- function(
  template,
  crosstalk_grid = seq(0, 1, by = 0.25),
  nsample = NULL,
  nrep = 5,
  p = NULL,
  max_features_per_view = p,
  dependence = c("factor", "empirical_copula", "module_copula"),
  cor.shrinkage = 0.05,
  module.count = NULL,
  module.between.shrinkage = 0.75,
  view.models = NULL,
  seed = 1,
  ...
) {
  dependence <- match.arg(dependence)
  observed_matrices <- if (inherits(template, "mimesys_template")) {
    lapply(template$fits, function(f) {
      x <- f$qmat
      colnames(x) <- rownames(f$Lambda)
      x
    })
  } else {
    as_mimesys_view_list(template)
  }
  tmpl <- if (inherits(template, "mimesys_template")) template else fit_template(template, ...)
  if (is.null(nsample)) {
    nsample <- tmpl$n0
  }

  rows <- list()
  idx <- 1L
  for (ct in crosstalk_grid) {
    sim <- gen_simmba(
      nsample = nsample,
      template = tmpl,
      p = p,
      crosstalk = ct,
      de.align = NULL,
      nrep = nrep,
      seed = seed + as.integer(round(10000 * ct)),
      dependence = dependence,
      cor.shrinkage = cor.shrinkage,
      module.count = module.count,
      module.between.shrinkage = module.between.shrinkage,
      view.models = view.models
    )
    for (rep_idx in seq_len(nrep)) {
      sim_mats <- simulation_view_matrices(
        sim,
        rep = rep_idx,
        split = "train",
        feature_names = "template_feature"
      )
      distances <- feature_correlation_frobenius(
        observed_matrices = observed_matrices,
        simulated_matrices = sim_mats,
        max_features_per_view = max_features_per_view
      )
      distances$crosstalk <- ct
      distances$rep <- rep_idx
      rows[[idx]] <- distances
      idx <- idx + 1L
    }
  }

  raw <- do.call(rbind, rows)
  summary <- .mimesys_summarize_frobenius(raw)
  cross <- summary[!summary$diagonal_block, , drop = FALSE]
  best <- cross[which.min(cross$mean), , drop = FALSE]

  list(
    best_crosstalk = best$crosstalk[[1]],
    best = best,
    summary = summary,
    raw = raw,
    dependence = dependence,
    cor.shrinkage = cor.shrinkage,
    module.count = module.count,
    module.between.shrinkage = module.between.shrinkage,
    view.models = view.models
  )
}

#' Plot crosstalk calibration distances
#'
#' Visualizes the output of \code{calibrate_crosstalk()} by plotting the mean
#' rescaled Frobenius distance between observed and simulated correlation blocks
#' across candidate crosstalk values. Lower values indicate better reconstruction
#' of the observed correlation structure.
#'
#' @param calibration A list returned by \code{calibrate_crosstalk()}.
#' @param include_diagonal Logical; if \code{FALSE}, only cross-view blocks are
#'   shown because those are the blocks used to choose the best crosstalk.
#' @param show_sd Logical; if \code{TRUE}, add \code{mean +/- sd} error bars when
#'   replicate-level standard deviations are available.
#'
#' @return A \code{ggplot} object.
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
    stop("calibration must be a list returned by calibrate_crosstalk().", call. = FALSE)
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
#'   summaries.
#' @export
summarize_simulation <- function(sim) {
  if (!is.list(sim) || is.null(sim$trainDat) || is.null(sim$testDat)) {
    stop("sim must be a list returned by gen_simmba().", call. = FALSE)
  }
  first <- sim$trainDat[[1]]
  if (is.null(first$feature_metadata) || is.null(first$sample_metadata)) {
    stop("sim does not look like a MimeSys simulation object.", call. = FALSE)
  }
  feature_counts <- table(first$feature_metadata$featureType)
  feature_counts <- stats::setNames(as.integer(feature_counts), names(feature_counts))
  train_n <- vapply(sim$trainDat, function(x) nrow(x$sample_metadata), integer(1))
  test_n <- vapply(sim$testDat, function(x) nrow(x$sample_metadata), integer(1))
  outcome_values <- unlist(lapply(sim$trainDat, function(x) x$sample_metadata$Y), use.names = FALSE)
  beta <- first$feature_metadata$beta
  active_features <- if (is.null(beta)) NA_integer_ else sum(beta != 0, na.rm = TRUE)

  signal_cor <- vapply(sim$trainDat, function(x) {
    sm <- x$sample_metadata
    if (!all(c("Y", "Xbeta") %in% names(sm))) {
      return(NA_real_)
    }
    suppressWarnings(stats::cor(sm$Y, sm$Xbeta, use = "complete.obs"))
  }, numeric(1))

  out <- list(
    nrep = length(sim$trainDat),
    train_n = stats::setNames(train_n, names(sim$trainDat)),
    test_n = stats::setNames(test_n, names(sim$testDat)),
    feature_counts = feature_counts,
    crosstalk = sim$crosstalk,
    dependence = sim$dependence,
    ygen.mode = sim$ygen.mode,
    outcome.type = sim$outcome.type,
    shared.latent.fun = sim$shared.latent.fun,
    shared.latent.n_factors = sim$shared.latent.n_factors,
    signal.share = sim$signal.share,
    friedman.share = sim$friedman.share,
    de.align = sim$de.align,
    active_features_first_rep = active_features,
    outcome_summary = summary(outcome_values),
    y_signal_cor_mean = mean(signal_cor, na.rm = TRUE),
    y_signal_cor_sd = if (sum(is.finite(signal_cor)) > 1L) stats::sd(signal_cor, na.rm = TRUE) else NA_real_
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
  cat("  Features: ", paste(names(x$feature_counts), x$feature_counts, sep = "=", collapse = ", "), "\n", sep = "")
  cat("  Crosstalk: ", paste(x$crosstalk, collapse = ", "), "\n", sep = "")
  cat("  Outcome: ", x$outcome.type, " | ", x$ygen.mode, "\n", sep = "")
  if (identical(x$ygen.mode, "SharedLatent")) {
    shared_fun <- if (is.null(x$shared.latent.fun)) "linear" else x$shared.latent.fun
    shared_n <- if (is.null(x$shared.latent.n_factors)) "all" else x$shared.latent.n_factors
    cat("  Shared latent function: ", shared_fun, "\n", sep = "")
    cat("  Shared latent factors: ", shared_n, "\n", sep = "")
  }
  cat("  Dependence: ", x$dependence, "\n", sep = "")
  if (!is.null(x$signal.share)) {
    cat("  Signal share: ", paste(names(x$signal.share), round(x$signal.share, 3), sep = "=", collapse = ", "), "\n", sep = "")
  }
  cat("  Active features in first replicate: ", x$active_features_first_rep, "\n", sep = "")
  cat("  Mean cor(Y, signal): ", round(x$y_signal_cor_mean, 3), "\n", sep = "")
  invisible(x)
}

#' Compute rescaled Frobenius correlation-block distances
#'
#' @param observed_matrices,simulated_matrices Named sample-by-feature matrix
#'   lists, MAE/SE objects, or compatible containers.
#' @param max_features_per_view Optional feature cap.
#'
#' @return Data frame with one row per upper-triangle view block.
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
  obs <- aligned$observed
  sim <- aligned$simulated
  views <- names(obs)
  obs_all <- do.call(cbind, obs)
  sim_all <- do.call(cbind, sim)
  obs_cor <- stats::cor(obs_all, use = "pairwise.complete.obs")
  sim_cor <- stats::cor(sim_all, use = "pairwise.complete.obs")

  p <- vapply(obs, ncol, integer(1))
  starts <- cumsum(c(1L, head(p, -1L)))
  ends <- cumsum(p)
  idx <- stats::setNames(Map(seq, starts, ends), views)

  rows <- list()
  k <- 1L
  for (i in seq_along(views)) {
    for (j in seq.int(i, length(views))) {
      view_x <- views[[i]]
      view_y <- views[[j]]
      block_obs <- obs_cor[idx[[view_x]], idx[[view_y]], drop = FALSE]
      block_sim <- sim_cor[idx[[view_x]], idx[[view_y]], drop = FALSE]
      if (identical(view_x, view_y)) {
        diag(block_obs) <- NA_real_
        diag(block_sim) <- NA_real_
      }
      diff <- block_sim - block_obs
      rows[[k]] <- data.frame(
        block = paste0("cor(", view_x, if (!identical(view_x, view_y)) paste0(",", view_y) else "", ")"),
        view_x = view_x,
        view_y = view_y,
        diagonal_block = identical(view_x, view_y),
        rescaled_frobenius = sqrt(sum(diff^2, na.rm = TRUE) / sum(is.finite(diff))),
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
    FUN = function(v) c(mean = mean(v), sd = if (length(v) > 1L) stats::sd(v) else NA_real_)
  )
  out <- do.call(data.frame, out)
  names(out) <- sub("rescaled_frobenius\\.", "", names(out))
  out
}
