#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) == 0) {
    return(default)
  }
  sub(paste0("^--", name, "="), "", hit[[1]])
}

base_dir <- path.expand(get_arg("base_dir", "~/Mimesys_validation"))
out_dir <- path.expand(get_arg("out_dir", file.path(base_dir, "tcga_validation_grid")))
plot_dir <- path.expand(get_arg("plot_dir", file.path(out_dir, "plots")))
task_dir <- file.path(out_dir, "tasks")
output_prefix <- get_arg("prefix", NULL)
pair_scatter_reps <- as.integer(get_arg("pair_scatter_reps", "1"))
max_pair_scatter_points <- as.integer(get_arg("max_pair_scatter_points", "500"))
save_pair_scatter <- as.logical(get_arg("save_pair_scatter", "TRUE"))
if (is.na(pair_scatter_reps) || pair_scatter_reps < 1L) {
  pair_scatter_reps <- 1L
}
if (is.na(max_pair_scatter_points) || max_pair_scatter_points < 1L) {
  max_pair_scatter_points <- 500L
}

suppressPackageStartupMessages({
  library(MimeSys)
})

save_plot_pdf <- function(plot, path, width = 11, height = 8) {
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

block_summary_long <- function(block_df, value_col) {
  out <- block_df[, c("view_x", "view_y", "diagonal_block", "crosstalk", value_col), drop = FALSE]
  names(out)[names(out) == value_col] <- "plot_value"
  out$crosstalk <- factor(out$crosstalk, levels = sort(unique(as.numeric(out$crosstalk))))
  out
}

plot_block_summary_by_crosstalk <- function(
  block_df,
  value_col,
  title,
  fill_title = "Block correlation",
  signed = FALSE
) {
  plot_df <- block_summary_long(block_df, value_col)
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = view_x, y = view_y, fill = plot_value)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f", plot_value)), size = 3) +
    ggplot2::facet_wrap(~ crosstalk, nrow = 1) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = title,
      x = "View",
      y = "View",
      fill = fill_title
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "white", color = "grey70"),
      plot.title = ggplot2::element_text(face = "bold")
    )
  if (isTRUE(signed)) {
    lim <- max(abs(plot_df$plot_value), na.rm = TRUE)
    if (!is.finite(lim) || lim <= 0) {
      lim <- 0.05
    }
    p + ggplot2::scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-lim, lim),
      na.value = "grey90"
    )
  } else {
    p + ggplot2::scale_fill_gradient(low = "grey95", high = "#2166AC", na.value = "grey90")
  }
}

plot_block_difference_by_crosstalk <- function(block_df, title) {
  plot_df <- block_summary_long(block_df, "difference")
  lim <- max(abs(plot_df$plot_value), na.rm = TRUE)
  if (!is.finite(lim) || lim <= 0) {
    lim <- 0.05
  }
  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = view_x, y = view_y, fill = plot_value)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%+.3f", plot_value)), size = 3) +
    ggplot2::facet_wrap(~ crosstalk, nrow = 1) +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-lim, lim),
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
      strip.background = ggplot2::element_rect(fill = "white", color = "grey70"),
      plot.title = ggplot2::element_text(face = "bold")
    )
}

plot_marginal_ks_by_crosstalk <- function(marginal_df, title) {
  plot_df <- marginal_df
  plot_df$crosstalk <- factor(plot_df$crosstalk, levels = sort(unique(as.numeric(plot_df$crosstalk))))
  ggplot2::ggplot(plot_df, ggplot2::aes(x = view, y = ks_stat, fill = view)) +
    ggplot2::geom_boxplot(width = 0.55, outlier.alpha = 0.18) +
    ggplot2::geom_jitter(width = 0.12, alpha = 0.18, size = 0.45) +
    ggplot2::facet_wrap(~ crosstalk, nrow = 1) +
    ggplot2::guides(fill = "none") +
    ggplot2::labs(
      title = title,
      x = "View",
      y = "KS statistic"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "white", color = "grey70"),
      plot.title = ggplot2::element_text(face = "bold")
    )
}

plot_pca_representatives_by_crosstalk <- function(pca_df, title) {
  plot_df <- pca_df
  reps <- stats::aggregate(grid_rep ~ crosstalk, data = plot_df, FUN = min)
  plot_df <- merge(plot_df, reps, by = "crosstalk", suffixes = c("", "_representative"))
  plot_df <- plot_df[plot_df$grid_rep == plot_df$grid_rep_representative, , drop = FALSE]
  plot_df$crosstalk <- factor(plot_df$crosstalk, levels = sort(unique(as.numeric(plot_df$crosstalk))))
  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = PC1, y = PC2, color = data_type)
  ) +
    ggplot2::geom_point(alpha = 0.6, size = 1.3) +
    ggplot2::stat_ellipse(linewidth = 0.6, alpha = 0.8, show.legend = FALSE) +
    ggplot2::facet_wrap(~ crosstalk, scales = "free") +
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

plot_modality_pca_representatives_by_crosstalk <- function(pca_df, title) {
  plot_df <- pca_df
  reps <- stats::aggregate(grid_rep ~ crosstalk + view, data = plot_df, FUN = min)
  plot_df <- merge(plot_df, reps, by = c("crosstalk", "view"), suffixes = c("", "_representative"))
  plot_df <- plot_df[plot_df$grid_rep == plot_df$grid_rep_representative, , drop = FALSE]
  plot_df$crosstalk <- factor(plot_df$crosstalk, levels = sort(unique(as.numeric(plot_df$crosstalk))))
  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = PC1, y = PC2, color = data_type)
  ) +
    ggplot2::geom_point(alpha = 0.6, size = 1.15) +
    ggplot2::stat_ellipse(linewidth = 0.55, alpha = 0.8, show.legend = FALSE) +
    ggplot2::facet_grid(view ~ crosstalk, scales = "free") +
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

plot_structure_metric_by_crosstalk <- function(structure_df, metric, title) {
  plot_df <- structure_df[, c("view_x", "view_y", "diagonal_block", "crosstalk", metric), drop = FALSE]
  names(plot_df)[names(plot_df) == metric] <- "plot_value"
  plot_df$crosstalk <- factor(plot_df$crosstalk, levels = sort(unique(as.numeric(plot_df$crosstalk))))
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = view_x, y = view_y, fill = plot_value)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f", plot_value)), size = 3) +
    ggplot2::facet_wrap(~ crosstalk, nrow = 1) +
    ggplot2::coord_equal() +
    ggplot2::labs(title = title, x = "View", y = "View", fill = metric) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "white", color = "grey70"),
      plot.title = ggplot2::element_text(face = "bold")
    )
  if (metric %in% c("cor_pearson", "cor_spearman")) {
    p + ggplot2::scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-1, 1),
      na.value = "grey90"
    )
  } else {
    p + ggplot2::scale_fill_gradient(low = "grey95", high = "#B2182B", na.value = "grey90")
  }
}

plot_structure_scatter_representatives <- function(pair_df, title, max_points_per_panel = 4000) {
  plot_df <- pair_df
  reps <- stats::aggregate(grid_rep ~ crosstalk + view_x + view_y, data = plot_df, FUN = min)
  plot_df <- merge(plot_df, reps, by = c("crosstalk", "view_x", "view_y"), suffixes = c("", "_representative"))
  plot_df <- plot_df[plot_df$grid_rep == plot_df$grid_rep_representative, , drop = FALSE]
  plot_df$view_pair <- paste(plot_df$view_x, plot_df$view_y, sep = "-")
  split_key <- interaction(plot_df$crosstalk, plot_df$view_pair, drop = TRUE)
  set.seed(1)
  keep <- unlist(lapply(split(seq_len(nrow(plot_df)), split_key), function(idx) {
    if (length(idx) <= max_points_per_panel) idx else sample(idx, max_points_per_panel)
  }), use.names = FALSE)
  plot_df <- plot_df[keep, , drop = FALSE]
  plot_df$crosstalk <- factor(plot_df$crosstalk, levels = sort(unique(as.numeric(plot_df$crosstalk))))

  ggplot2::ggplot(plot_df, ggplot2::aes(x = observed_cor, y = simulated_cor)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, color = "grey35", linewidth = 0.4) +
    ggplot2::geom_point(alpha = 0.18, size = 0.35, color = "#2166AC") +
    ggplot2::facet_grid(view_pair ~ crosstalk) +
    ggplot2::coord_equal(xlim = c(-1, 1), ylim = c(-1, 1)) +
    ggplot2::labs(
      title = title,
      x = "Observed feature-pair correlation",
      y = "Simulated feature-pair correlation"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "white", color = "grey70"),
      plot.title = ggplot2::element_text(face = "bold")
    )
}

files <- list.files(task_dir, pattern = "^task_[0-9]+\\.rds$", full.names = TRUE)
if (length(files) == 0) {
  stop("No task RDS files found in: ", task_dir)
}

extract_piece <- function(x, name) {
  if (is.data.frame(x) && identical(name, "diagnostics")) {
    return(x)
  }
  if (is.list(x) && !is.null(x[[name]])) {
    return(x[[name]])
  }
  NULL
}

standardize_task_df <- function(x) {
  if (is.null(x) || !is.data.frame(x)) {
    return(x)
  }
  if (!"dataset" %in% names(x)) {
    if ("cancer" %in% names(x)) {
      x$dataset <- x$cancer
    } else {
      x$dataset <- "dataset"
    }
  }
  if (!"cancer" %in% names(x)) {
    x$cancer <- x$dataset
  }
  if (!"feature_counts" %in% names(x)) {
    feature_cols <- grep("_features$", names(x), value = TRUE)
    if (length(feature_cols) > 0L) {
      x$feature_counts <- apply(x[, feature_cols, drop = FALSE], 1, function(z) {
        vals <- as.character(z)
        names(vals) <- sub("_features$", "", feature_cols)
        paste(paste(names(vals), vals, sep = "="), collapse = ",")
      })
    } else {
      x$feature_counts <- NA_character_
    }
  }
  x
}

safe_bind <- function(pieces) {
  pieces <- Filter(Negate(is.null), pieces)
  if (length(pieces) == 0L) {
    return(NULL)
  }
  do.call(rbind, lapply(pieces, standardize_task_df))
}

sample_pair_scatter <- function(pair_df, max_points_per_panel) {
  if (is.null(pair_df) || !is.data.frame(pair_df) || nrow(pair_df) == 0L) {
    return(NULL)
  }
  if (!all(c("observed_cor", "simulated_cor", "view_x", "view_y") %in% names(pair_df))) {
    return(NULL)
  }
  if ("grid_rep" %in% names(pair_df)) {
    pair_df <- pair_df[as.integer(pair_df$grid_rep) <= pair_scatter_reps, , drop = FALSE]
  }
  if (nrow(pair_df) == 0L) {
    return(NULL)
  }
  group_cols <- intersect(
    c("dataset", "feature_setting", "n_subjects", "de_align_label", "crosstalk", "grid_rep", "view_x", "view_y"),
    names(pair_df)
  )
  split_key <- if (length(group_cols) > 0L) {
    interaction(pair_df[, group_cols, drop = FALSE], drop = TRUE, sep = "\r")
  } else {
    rep("all", nrow(pair_df))
  }
  set.seed(1)
  keep <- unlist(lapply(split(seq_len(nrow(pair_df)), split_key), function(idx) {
    if (length(idx) <= max_points_per_panel) {
      idx
    } else {
      sample(idx, max_points_per_panel)
    }
  }), use.names = FALSE)
  pair_df[keep, , drop = FALSE]
}

sanitize_tag <- function(x) {
  x <- paste(x, collapse = "_")
  gsub("[^A-Za-z0-9_.-]+", "_", x)
}

diagnostics_pieces <- vector("list", length(files))
marginal_pieces <- vector("list", length(files))
pca_pieces <- vector("list", length(files))
modality_pca_pieces <- vector("list", length(files))
block_comparison_pieces <- vector("list", length(files))
block_structure_pieces <- vector("list", length(files))
block_pair_pieces <- vector("list", length(files))

for (file_idx in seq_along(files)) {
  x <- readRDS(files[[file_idx]])
  diagnostics_pieces[[file_idx]] <- extract_piece(x, "diagnostics")
  marginal_pieces[[file_idx]] <- extract_piece(x, "marginal_diagnostics")
  pca_pieces[[file_idx]] <- extract_piece(x, "pca_diagnostics")
  modality_pca_pieces[[file_idx]] <- extract_piece(x, "modality_pca_diagnostics")
  block_comparison_pieces[[file_idx]] <- extract_piece(x, "block_comparison")

  structure_out <- extract_piece(x, "block_structure_agreement")
  if (is.null(structure_out)) {
    structure_out <- extract_piece(x, "correlation_structure_diagnostics")
    if (!is.null(structure_out)) {
      if ("structure_cor" %in% names(structure_out) && !"cor_pearson" %in% names(structure_out)) {
        structure_out$cor_pearson <- structure_out$structure_cor
      }
      if ("mean_abs_error" %in% names(structure_out) && !"mae_cor" %in% names(structure_out)) {
        structure_out$mae_cor <- structure_out$mean_abs_error
      }
    }
  }
  block_structure_pieces[[file_idx]] <- structure_out

  if (isTRUE(save_pair_scatter)) {
    pair_out <- extract_piece(x, "block_structure_pairs")
    pair_out <- standardize_task_df(pair_out)
    block_pair_pieces[[file_idx]] <- sample_pair_scatter(pair_out, max_pair_scatter_points)
  }

  rm(x)
  if (file_idx %% 100L == 0L) {
    gc()
  }
}

diagnostics <- safe_bind(diagnostics_pieces)
marginal_diagnostics <- safe_bind(marginal_pieces)
pca_diagnostics <- safe_bind(pca_pieces)
modality_pca_diagnostics <- safe_bind(modality_pca_pieces)
block_comparisons <- safe_bind(block_comparison_pieces)
block_structure_agreement <- safe_bind(block_structure_pieces)
block_structure_pairs <- safe_bind(block_pair_pieces)

if (is.null(diagnostics) || nrow(diagnostics) == 0L) {
  stop("No diagnostics found in task RDS files.", call. = FALSE)
}

if (is.null(output_prefix) || !nzchar(output_prefix)) {
  output_prefix <- sanitize_tag(sort(unique(diagnostics$dataset)))
}

summary <- stats::aggregate(
  cbind(mean_abs_feature_cor, signal_agreement, y_signal_cor) ~
    dataset + feature_setting + n_subjects + de_align_label + crosstalk + view_a + view_b,
  data = diagnostics,
  FUN = function(x) mean(x, na.rm = TRUE)
)

robustness <- stats::aggregate(
  cbind(mean_abs_feature_cor, signal_agreement, y_signal_cor) ~
    dataset + feature_setting + n_subjects + de_align_label + crosstalk,
  data = diagnostics,
  FUN = function(x) mean(x, na.rm = TRUE)
)

saveRDS(diagnostics, file.path(out_dir, paste0(output_prefix, "_grid_diagnostics_raw.rds")))
write.csv(diagnostics, file.path(out_dir, paste0(output_prefix, "_grid_diagnostics_raw.csv")), row.names = FALSE)

saveRDS(summary, file.path(out_dir, paste0(output_prefix, "_grid_diagnostics_summary.rds")))
write.csv(summary, file.path(out_dir, paste0(output_prefix, "_grid_diagnostics_summary.csv")), row.names = FALSE)

saveRDS(robustness, file.path(out_dir, paste0(output_prefix, "_grid_robustness_table.rds")))
write.csv(robustness, file.path(out_dir, paste0(output_prefix, "_grid_robustness_table.csv")), row.names = FALSE)

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

if (!is.null(marginal_diagnostics)) {
  marginal_summary <- stats::aggregate(
    cbind(ks_stat, mean_diff, sd_ratio, median_diff) ~
      dataset + feature_setting + n_subjects + de_align_label + crosstalk + view,
    data = marginal_diagnostics,
    FUN = function(x) mean(x, na.rm = TRUE)
  )
  saveRDS(marginal_diagnostics, file.path(out_dir, paste0(output_prefix, "_marginal_diagnostics_raw.rds")))
  write.csv(marginal_diagnostics, file.path(out_dir, paste0(output_prefix, "_marginal_diagnostics_raw.csv")), row.names = FALSE)
  saveRDS(marginal_summary, file.path(out_dir, paste0(output_prefix, "_marginal_diagnostics_summary.rds")))
  write.csv(marginal_summary, file.path(out_dir, paste0(output_prefix, "_marginal_diagnostics_summary.csv")), row.names = FALSE)
}

if (!is.null(block_comparisons)) {
  if (!"value" %in% names(block_comparisons)) {
    block_comparisons$value <- "mean_abs_cor"
  }
  block_summary <- stats::aggregate(
    cbind(observed, simulated, difference) ~
      dataset + feature_setting + n_subjects + de_align_label + crosstalk +
      value + view_x + view_y + diagonal_block,
    data = block_comparisons,
    FUN = function(x) mean(x, na.rm = TRUE)
  )
  saveRDS(block_comparisons, file.path(out_dir, paste0(output_prefix, "_block_comparisons_raw.rds")))
  write.csv(block_comparisons, file.path(out_dir, paste0(output_prefix, "_block_comparisons_raw.csv")), row.names = FALSE)
  saveRDS(block_summary, file.path(out_dir, paste0(output_prefix, "_block_comparisons_summary.rds")))
  write.csv(block_summary, file.path(out_dir, paste0(output_prefix, "_block_comparisons_summary.csv")), row.names = FALSE)
}

if (!is.null(block_structure_agreement)) {
  needed <- c("cor_pearson", "cor_spearman", "mae_cor", "rmse_cor")
  for (nm in needed) {
    if (!nm %in% names(block_structure_agreement)) {
      block_structure_agreement[[nm]] <- NA_real_
    }
  }
  block_structure_summary <- stats::aggregate(
    cbind(cor_pearson, cor_spearman, mae_cor, rmse_cor, n_pairs) ~
      dataset + feature_setting + n_subjects + de_align_label + crosstalk +
      view_x + view_y + diagonal_block,
    data = block_structure_agreement,
    FUN = function(x) mean(x, na.rm = TRUE)
  )
  saveRDS(block_structure_agreement, file.path(out_dir, paste0(output_prefix, "_block_structure_agreement_raw.rds")))
  write.csv(block_structure_agreement, file.path(out_dir, paste0(output_prefix, "_block_structure_agreement_raw.csv")), row.names = FALSE)
  saveRDS(block_structure_summary, file.path(out_dir, paste0(output_prefix, "_block_structure_agreement_summary.rds")))
  write.csv(block_structure_summary, file.path(out_dir, paste0(output_prefix, "_block_structure_agreement_summary.csv")), row.names = FALSE)
}

if (!is.null(block_structure_pairs)) {
  saveRDS(block_structure_pairs, file.path(out_dir, paste0(output_prefix, "_block_structure_pairs_sampled.rds")))
}

if (!is.null(pca_diagnostics)) {
  saveRDS(pca_diagnostics, file.path(out_dir, paste0(output_prefix, "_pca_diagnostics_raw.rds")))
  write.csv(pca_diagnostics, file.path(out_dir, paste0(output_prefix, "_pca_diagnostics_raw.csv")), row.names = FALSE)
}

if (!is.null(modality_pca_diagnostics)) {
  saveRDS(modality_pca_diagnostics, file.path(out_dir, paste0(output_prefix, "_modality_pca_diagnostics_raw.rds")))
  write.csv(modality_pca_diagnostics, file.path(out_dir, paste0(output_prefix, "_modality_pca_diagnostics_raw.csv")), row.names = FALSE)
}

plot_groups <- unique(diagnostics[, c("dataset", "feature_setting", "n_subjects")])

for (i in seq_len(nrow(plot_groups))) {
  group <- plot_groups[i, , drop = FALSE]
  group_keep <- diagnostics$dataset == group$dataset[[1]] &
    diagnostics$feature_setting == group$feature_setting[[1]] &
    diagnostics$n_subjects == group$n_subjects[[1]]

  prefix <- paste(
    group$dataset[[1]],
    group$feature_setting[[1]],
    paste0("n", group$n_subjects[[1]]),
    "full_validation",
    sep = "_"
  )

  group_diag <- diagnostics[group_keep, , drop = FALSE]
  all_plots <- MimeSys:::.mimesys_template_diagnostic_plot_list(group_diag)
  all_plots$original_crosstalk_line <- NULL
  all_plots$original_feature_correlation_heatmap <- NULL

  if (!is.null(marginal_diagnostics)) {
    keep <- marginal_diagnostics$dataset == group$dataset[[1]] &
      marginal_diagnostics$feature_setting == group$feature_setting[[1]] &
      marginal_diagnostics$n_subjects == group$n_subjects[[1]]
    marg_group <- marginal_diagnostics[keep, , drop = FALSE]
    if (nrow(marg_group) > 0) {
      p_marg <- plot_marginal_ks_by_crosstalk(
        marg_group,
        title = paste(prefix, "marginal KS by view across crosstalk")
      )
      save_plot_pdf(
        p_marg,
        file.path(plot_dir, paste0(prefix, "_marginal_ks_by_crosstalk.pdf")),
        width = 12,
        height = 7
      )
      all_plots$marginal_ks_by_crosstalk <- p_marg
    }
  }

  if (!is.null(block_comparisons)) {
    keep <- block_comparisons$dataset == group$dataset[[1]] &
      block_comparisons$feature_setting == group$feature_setting[[1]] &
      block_comparisons$n_subjects == group$n_subjects[[1]]
    block_group <- block_comparisons[keep, , drop = FALSE]
    if (nrow(block_group) > 0) {
      if (!"value" %in% names(block_group)) {
        block_group$value <- "mean_abs_cor"
      }
      block_mean <- stats::aggregate(
        cbind(observed, simulated, difference) ~ value + crosstalk + view_x + view_y + diagonal_block,
        data = block_group,
        FUN = function(x) mean(x, na.rm = TRUE)
      )

      block_path <- file.path(plot_dir, paste0(prefix, "_block_validation_by_crosstalk.pdf"))
      grDevices::pdf(block_path, width = 12, height = 7, onefile = TRUE)
      tryCatch(
        {
          for (metric in intersect(c("mean_cor", "mean_abs_cor", "median_cor", "median_abs_cor"), unique(block_mean$value))) {
            metric_group <- block_mean[block_mean$value == metric, , drop = FALSE]
            signed_metric <- metric %in% c("mean_cor", "median_cor")
            observed_once <- stats::aggregate(
              observed ~ view_x + view_y + diagonal_block,
              data = metric_group,
              FUN = function(x) mean(x, na.rm = TRUE)
            )
            names(observed_once)[names(observed_once) == "observed"] <- metric

            p_obs <- plot_blockwise_correlation_summary(
              observed_once,
              value = metric,
              title = paste(prefix, "observed block-level", metric)
            )
            if (isTRUE(signed_metric)) {
              p_obs <- p_obs +
                ggplot2::scale_fill_gradient2(
                  low = "#2166AC",
                  mid = "white",
                  high = "#B2182B",
                  midpoint = 0,
                  na.value = "grey90"
                )
            }
            p_sim <- plot_block_summary_by_crosstalk(
              metric_group,
              value_col = "simulated",
              title = paste(prefix, "simulated block-level", metric),
              fill_title = metric,
              signed = signed_metric
            )
            p_diff <- plot_block_difference_by_crosstalk(
              metric_group,
              title = paste(prefix, "simulated minus observed block-level", metric)
            )
            print(p_obs)
            print(p_sim)
            print(p_diff)

            all_plots[[paste0("observed_block_", metric)]] <- p_obs
            all_plots[[paste0("simulated_block_", metric, "_by_crosstalk")]] <- p_sim
            all_plots[[paste0("block_difference_", metric, "_by_crosstalk")]] <- p_diff
          }
        },
        finally = grDevices::dev.off()
      )
    }
  }

  if (!is.null(block_structure_agreement) && nrow(block_structure_agreement) > 0) {
    keep <- block_structure_agreement$dataset == group$dataset[[1]] &
      block_structure_agreement$feature_setting == group$feature_setting[[1]] &
      block_structure_agreement$n_subjects == group$n_subjects[[1]]
    structure_group <- block_structure_agreement[keep, , drop = FALSE]
    if (nrow(structure_group) > 0) {
      for (nm in c("cor_pearson", "cor_spearman", "mae_cor", "rmse_cor")) {
        if (!nm %in% names(structure_group)) {
          structure_group[[nm]] <- NA_real_
        }
      }
      structure_mean <- stats::aggregate(
        cbind(cor_pearson, cor_spearman, mae_cor, rmse_cor, n_pairs) ~
          crosstalk + view_x + view_y + diagonal_block,
        data = structure_group,
        FUN = function(x) mean(x, na.rm = TRUE)
      )
      p_corr_corr <- plot_structure_metric_by_crosstalk(
        structure_mean,
        metric = "cor_pearson",
        title = paste(prefix, "corr-of-corr structure agreement")
      )
      p_corr_error <- plot_structure_metric_by_crosstalk(
        structure_mean,
        metric = "mae_cor",
        title = paste(prefix, "feature-correlation MAE")
      )
      save_plot_pdf(
        p_corr_corr,
        file.path(plot_dir, paste0(prefix, "_corr_of_corr_by_crosstalk.pdf")),
        width = 12,
        height = 7
      )
      save_plot_pdf(
        p_corr_error,
        file.path(plot_dir, paste0(prefix, "_correlation_error_by_crosstalk.pdf")),
        width = 12,
        height = 7
      )
      all_plots$corr_of_corr_by_crosstalk <- p_corr_corr
      all_plots$correlation_error_by_crosstalk <- p_corr_error
    }
  }

  if (!is.null(block_structure_pairs) && nrow(block_structure_pairs) > 0) {
    keep <- block_structure_pairs$dataset == group$dataset[[1]] &
      block_structure_pairs$feature_setting == group$feature_setting[[1]] &
      block_structure_pairs$n_subjects == group$n_subjects[[1]]
    pair_group <- block_structure_pairs[keep, , drop = FALSE]
    if (nrow(pair_group) > 0) {
      p_corr_scatter <- plot_structure_scatter_representatives(
        pair_group,
        title = paste(prefix, "observed versus simulated feature-pair correlations")
      )
      save_plot_pdf(
        p_corr_scatter,
        file.path(plot_dir, paste0(prefix, "_observed_vs_simulated_correlation_scatter.pdf")),
        width = 14,
        height = 9
      )
      all_plots$observed_vs_simulated_correlation_scatter <- p_corr_scatter
    }
  }

  if (!is.null(pca_diagnostics) && nrow(pca_diagnostics) > 0) {
    keep <- pca_diagnostics$dataset == group$dataset[[1]] &
      pca_diagnostics$feature_setting == group$feature_setting[[1]] &
      pca_diagnostics$n_subjects == group$n_subjects[[1]]
    pca_group <- pca_diagnostics[keep, , drop = FALSE]
    if (nrow(pca_group) > 0) {
      p_pca <- plot_pca_representatives_by_crosstalk(
        pca_group,
        title = paste(prefix, "observed versus simulated PCA representatives")
      )
      save_plot_pdf(
        p_pca,
        file.path(plot_dir, paste0(prefix, "_pca_representatives_by_crosstalk.pdf")),
        width = 12,
        height = 8
      )
      all_plots$pca_representatives_by_crosstalk <- p_pca
    }
  }

  if (!is.null(modality_pca_diagnostics) && nrow(modality_pca_diagnostics) > 0) {
    keep <- modality_pca_diagnostics$dataset == group$dataset[[1]] &
      modality_pca_diagnostics$feature_setting == group$feature_setting[[1]] &
      modality_pca_diagnostics$n_subjects == group$n_subjects[[1]]
    modality_pca_group <- modality_pca_diagnostics[keep, , drop = FALSE]
    if (nrow(modality_pca_group) > 0) {
      p_modality_pca <- plot_modality_pca_representatives_by_crosstalk(
        modality_pca_group,
        title = paste(prefix, "modality-wise observed versus simulated PCA representatives")
      )
      save_plot_pdf(
        p_modality_pca,
        file.path(plot_dir, paste0(prefix, "_modality_pca_representatives_by_crosstalk.pdf")),
        width = 14,
        height = 9
      )
      all_plots$modality_pca_representatives_by_crosstalk <- p_modality_pca
    }
  }

  full_path <- file.path(plot_dir, paste0(prefix, "_all_diagnostics.pdf"))
  grDevices::pdf(full_path, width = 12, height = 8, onefile = TRUE)
  tryCatch(
    {
      for (plot_name in names(all_plots)) {
        print(all_plots[[plot_name]])
      }
    },
    finally = grDevices::dev.off()
  )
}

cat("Combined task files:", length(files), "\n")
cat("Raw rows:", nrow(diagnostics), "\n")
cat("Saved combined outputs to:\n")
cat(out_dir, "\n")
cat("Saved plots to:\n")
cat(plot_dir, "\n")
