#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(MimeSys)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) == 0L) {
    return(default)
  }
  sub(paste0("^--", name, "="), "", hit[[1]])
}

split_arg <- function(x, sep = ",") {
  if (is.null(x) || !nzchar(x)) {
    return(character())
  }
  out <- trimws(strsplit(x, sep, fixed = TRUE)[[1]])
  out[nzchar(out)]
}

parse_counts <- function(x) {
  pieces <- split_arg(x, ",")
  out <- vector("list", length(pieces))
  names(out) <- character(length(pieces))
  for (i in seq_along(pieces)) {
    kv <- strsplit(pieces[[i]], "=", fixed = TRUE)[[1]]
    if (length(kv) != 2L) {
      stop("--feature_counts entries must look like view=count.", call. = FALSE)
    }
    names(out)[[i]] <- trimws(kv[[1]])
    out[[i]] <- as.integer(trimws(kv[[2]]))
  }
  unlist(out, use.names = TRUE)
}

parse_key_values <- function(x) {
  pieces <- split_arg(x, ",")
  if (length(pieces) == 0L) {
    return(NULL)
  }
  out <- vector("list", length(pieces))
  names(out) <- character(length(pieces))
  for (i in seq_along(pieces)) {
    kv <- strsplit(pieces[[i]], "=", fixed = TRUE)[[1]]
    if (length(kv) != 2L) {
      stop("Key-value entries must look like view=value.", call. = FALSE)
    }
    names(out)[[i]] <- trimws(kv[[1]])
    out[[i]] <- trimws(kv[[2]])
  }
  unlist(out, use.names = TRUE)
}

parse_numeric_grid <- function(x) {
  as.numeric(split_arg(x, ","))
}

keep_top_variable <- function(x, n_features) {
  vars <- apply(x, 2, stats::var, na.rm = TRUE)
  vars <- sort(vars[is.finite(vars)], decreasing = TRUE)
  keep <- names(vars)[seq_len(min(as.integer(n_features), length(vars)))]
  x[, keep, drop = FALSE]
}

apply_scenario_transform <- function(mats, scenario) {
  if (identical(scenario, "raw")) {
    return(mats)
  }
  if (identical(scenario, "clrlog")) {
    return(MimeSys:::.mimesys_apply_view_transforms(
      mats,
      view.transform = c(microbiome = "clr", metabolite = "log1p")
    ))
  }
  stop("Unknown scenario: ", scenario, call. = FALSE)
}

aligned_cor_matrix <- function(mats, caps) {
  views <- names(caps)
  mats <- lapply(views, function(view) keep_top_variable(mats[[view]], caps[[view]]))
  names(mats) <- views
  common_samples <- Reduce(intersect, lapply(mats, rownames))
  mats <- lapply(mats, function(x) x[common_samples, , drop = FALSE])
  x <- do.call(cbind, lapply(views, function(view) {
    mat <- mats[[view]]
    colnames(mat) <- paste(view, colnames(mat), sep = "::")
    mat
  }))
  corr <- stats::cor(x, use = "pairwise.complete.obs")
  list(corr = corr, mats = mats, views = views)
}

block_indices <- function(mats) {
  starts <- cumsum(c(1L, head(vapply(mats, ncol, integer(1)), -1L)))
  ends <- cumsum(vapply(mats, ncol, integer(1)))
  stats::setNames(Map(seq, starts, ends), names(mats))
}

rescaled_frobenius <- function(observed, simulated, idx, block_x, block_y) {
  block_obs <- observed[idx[[block_x]], idx[[block_y]], drop = FALSE]
  block_sim <- simulated[idx[[block_x]], idx[[block_y]], drop = FALSE]
  if (identical(block_x, block_y)) {
    diag(block_obs) <- NA_real_
    diag(block_sim) <- NA_real_
  }
  diff <- block_sim - block_obs
  sqrt(sum(diff^2, na.rm = TRUE) / sum(is.finite(diff)))
}

corr_heatmap_df <- function(corr, idx, view, source, scenario, crosstalk) {
  j <- idx[[view]]
  block <- corr[j, j, drop = FALSE]
  diag(block) <- NA_real_
  out <- as.data.frame(as.table(block), stringsAsFactors = FALSE)
  names(out) <- c("feature_x", "feature_y", "correlation")
  out$order_x <- match(out$feature_x, colnames(block))
  out$order_y <- match(out$feature_y, rownames(block))
  out$view <- view
  out$source <- source
  out$scenario <- scenario
  out$crosstalk <- crosstalk
  out
}

full_corr_heatmap_df <- function(corr, idx, source, scenario, crosstalk) {
  out <- as.data.frame(as.table(corr), stringsAsFactors = FALSE)
  names(out) <- c("feature_x", "feature_y", "correlation")
  out$order_x <- match(out$feature_x, colnames(corr))
  out$order_y <- match(out$feature_y, rownames(corr))
  out$source <- source
  out$scenario <- scenario
  out$crosstalk <- crosstalk
  out
}

block_boundary_df <- function(idx) {
  ends <- vapply(idx, max, integer(1))
  data.frame(boundary = head(ends, -1L) + 0.5)
}

block_label_df <- function(idx) {
  data.frame(
    view = names(idx),
    pos = vapply(idx, function(x) mean(range(x)), numeric(1)),
    stringsAsFactors = FALSE
  )
}

run_scenario <- function(template, scenario, caps, nsample, nrep, crosstalk_grid, seed,
                         dependence, cor_shrinkage, module_count, module_between_shrinkage,
                         view_models) {
  observed_mats <- apply_scenario_transform(template, scenario)
  observed <- aligned_cor_matrix(observed_mats, caps)
  tmpl <- fit_template(observed$mats)

  rows <- list()
  heatmaps <- list()
  full_heatmaps <- list()
  idx <- block_indices(observed$mats)
  for (view in observed$views) {
    heatmaps[[length(heatmaps) + 1L]] <- corr_heatmap_df(observed$corr, idx, view, "Empirical", scenario, NA_real_)
  }
  full_heatmaps[[length(full_heatmaps) + 1L]] <- full_corr_heatmap_df(observed$corr, idx, "Empirical", scenario, NA_real_)

  for (ct_idx in seq_along(crosstalk_grid)) {
    ct <- crosstalk_grid[[ct_idx]]
    for (rep_idx in seq_len(nrep)) {
      sim <- gen_simmba(
        nsample = nsample,
        template = tmpl,
        crosstalk = ct,
        de.align = NULL,
        nrep = 1,
        snr = 2,
        seed = seed + 1000L * ct_idx + rep_idx,
        dependence = dependence,
        cor.shrinkage = cor_shrinkage,
        module.count = module_count,
        module.between.shrinkage = module_between_shrinkage,
        view.models = view_models
      )
      sim_mats <- simulation_view_matrices(sim, rep = 1, split = "train", feature_names = "template_feature")
      sim_aligned <- aligned_cor_matrix(sim_mats, caps)
      sim_corr <- sim_aligned$corr[colnames(observed$corr), colnames(observed$corr), drop = FALSE]

      for (view_x in observed$views) {
        for (view_y in observed$views) {
          if (match(view_x, observed$views) <= match(view_y, observed$views)) {
            rows[[length(rows) + 1L]] <- data.frame(
              scenario = scenario,
              crosstalk = ct,
              rep = rep_idx,
              block = paste0("cor(", view_x, if (!identical(view_x, view_y)) paste0(",", view_y) else "", ")"),
              view_x = view_x,
              view_y = view_y,
              diagonal_block = identical(view_x, view_y),
              rescaled_frobenius = rescaled_frobenius(observed$corr, sim_corr, idx, view_x, view_y),
              stringsAsFactors = FALSE
            )
          }
        }
      }

      if (rep_idx == 1L) {
        for (view in observed$views) {
          heatmaps[[length(heatmaps) + 1L]] <- corr_heatmap_df(sim_corr, idx, view, "MimeSys", scenario, ct)
        }
        full_heatmaps[[length(full_heatmaps) + 1L]] <- full_corr_heatmap_df(sim_corr, idx, "MimeSys", scenario, ct)
      }
    }
  }

  frobenius <- do.call(rbind, rows)
  cross_summary <- stats::aggregate(
    rescaled_frobenius ~ crosstalk,
    data = frobenius[!frobenius$diagonal_block, , drop = FALSE],
    FUN = mean
  )
  best_crosstalk <- cross_summary$crosstalk[which.min(cross_summary$rescaled_frobenius)]

  heatmaps <- do.call(rbind, heatmaps)
  heatmaps <- heatmaps[
    heatmaps$source == "Empirical" | heatmaps$crosstalk == best_crosstalk,
    ,
    drop = FALSE
  ]
  heatmaps$crosstalk[heatmaps$source == "Empirical"] <- best_crosstalk

  full_heatmaps <- do.call(rbind, full_heatmaps)
  full_heatmaps <- full_heatmaps[
    full_heatmaps$source == "Empirical" | full_heatmaps$crosstalk == best_crosstalk,
    ,
    drop = FALSE
  ]
  full_heatmaps$crosstalk[full_heatmaps$source == "Empirical"] <- best_crosstalk

  list(
    frobenius = frobenius,
    heatmaps = heatmaps,
    full_heatmaps = full_heatmaps,
    idx = idx,
    best_crosstalk = data.frame(
      scenario = scenario,
      best_crosstalk = best_crosstalk,
      mean_cross_view_rescaled_frobenius = min(cross_summary$rescaled_frobenius),
      stringsAsFactors = FALSE
    )
  )
}

summarize_frobenius <- function(frobenius_raw) {
  out <- stats::aggregate(
    rescaled_frobenius ~ scenario + crosstalk + block + view_x + view_y + diagonal_block,
    data = frobenius_raw,
    FUN = function(x) c(mean = mean(x), sd = if (length(x) > 1L) stats::sd(x) else NA_real_)
  )
  out <- do.call(data.frame, out)
  names(out) <- sub("rescaled_frobenius\\.", "", names(out))
  out
}

summarize_crosstalk_calibration <- function(frobenius_raw) {
  cross <- frobenius_raw[!frobenius_raw$diagonal_block, , drop = FALSE]
  out <- stats::aggregate(
    rescaled_frobenius ~ scenario + crosstalk,
    data = cross,
    FUN = function(x) c(mean = mean(x), sd = if (length(x) > 1L) stats::sd(x) else NA_real_)
  )
  out <- do.call(data.frame, out)
  names(out) <- sub("rescaled_frobenius\\.", "", names(out))
  out
}

best_crosstalk_table <- function(calibration_summary) {
  do.call(rbind, lapply(split(calibration_summary, calibration_summary$scenario), function(x) {
    x[which.min(x$mean), , drop = FALSE]
  }))
}

make_pi_frobenius_report <- function(frobenius_summary, best_crosstalk) {
  best <- best_crosstalk[, c("scenario", "crosstalk"), drop = FALSE]
  out <- merge(frobenius_summary, best, by = c("scenario", "crosstalk"))
  out$block_type <- ifelse(out$diagonal_block, "within-view", "cross-view")
  out <- out[order(out$scenario, out$diagonal_block, out$block), ]
  out[, c(
    "scenario",
    "crosstalk",
    "block_type",
    "block",
    "view_x",
    "view_y",
    "mean",
    "sd"
  )]
}

plot_crosstalk_calibration <- function(calibration_summary, best_crosstalk, out_dir, feature_counts, nsample) {
  scenario_labels <- c(raw = "No transform", clrlog = "CLR microbiome + log metabolite")
  calibration_summary$scenario <- factor(calibration_summary$scenario,
    levels = names(scenario_labels),
    labels = unname(scenario_labels)
  )
  best_crosstalk$scenario <- factor(best_crosstalk$scenario,
    levels = names(scenario_labels),
    labels = unname(scenario_labels)
  )

  p <- ggplot(calibration_summary, aes(x = crosstalk, y = mean)) +
    geom_line(color = "#2B6CB0", linewidth = 0.6) +
    geom_point(color = "#2B6CB0", size = 2) +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.025, color = "grey35") +
    geom_vline(data = best_crosstalk, aes(xintercept = crosstalk), linetype = "dashed", color = "#B2182B") +
    facet_wrap(~scenario, nrow = 1) +
    labs(
      title = "FranzosaE 2019 empirical crosstalk calibration",
      subtitle = paste0("nsample=", nsample, ", features=", paste(paste(names(feature_counts), feature_counts, sep = "="), collapse = ", ")),
      x = "Crosstalk",
      y = "Cross-view rescaled Frobenius distance"
    ) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))

  ggsave(
    filename = file.path(out_dir, "franzosa_crosstalk_calibration.pdf"),
    plot = p,
    width = 8,
    height = 4.5,
    device = grDevices::pdf
  )
}

write_readme <- function(out_dir, best_crosstalk, feature_counts, nsample, nrep,
                         crosstalk_grid, dependence, cor_shrinkage,
                         module_count, module_between_shrinkage, view_models) {
  lines <- c(
    "FranzosaE 2019 crosstalk calibration",
    "",
    paste0("nsample: ", nsample),
    paste0("nrep per crosstalk: ", nrep),
    paste0("dependence: ", dependence),
    paste0("cor_shrinkage: ", cor_shrinkage),
    paste0("module_count: ", if (is.null(module_count)) "auto" else module_count),
    paste0("module_between_shrinkage: ", module_between_shrinkage),
    paste0("view_models: ", if (is.null(view_models)) "continuous" else paste(paste(names(view_models), view_models, sep = "="), collapse = ",")),
    paste0("features: ", paste(paste(names(feature_counts), feature_counts, sep = "="), collapse = ", ")),
    paste0("crosstalk grid: ", paste(crosstalk_grid, collapse = ", ")),
    "",
    "Best crosstalk is selected by minimizing the mean cross-view rescaled Frobenius distance between empirical and simulated feature-correlation blocks.",
    "",
    paste(capture.output(print(best_crosstalk, row.names = FALSE)), collapse = "\n")
  )
  writeLines(lines, file.path(out_dir, "README_crosstalk_calibration.txt"))
}

select_best_heatmap_data <- function(heatmap_df, best_crosstalk) {
  merge(heatmap_df, best_crosstalk[, c("scenario", "crosstalk")], by = c("scenario", "crosstalk"))
}

plot_pi_frobenius_table <- function(pi_frobenius_report) {
  report <- pi_frobenius_report
  report$transform <- ifelse(report$scenario == "clrlog", "CLR/log", "Raw")
  report$block_label <- ifelse(
    report$block_type == "cross-view",
    paste(report$view_x, report$view_y, sep = " vs "),
    report$view_x
  )
  report$mean <- sprintf("%.3f", report$mean)
  report$sd <- sprintf("%.3f", report$sd)
  report$crosstalk <- sprintf("%.2f", report$crosstalk)
  report <- report[order(report$transform, report$block_type, report$block_label), ]

  body <- data.frame(
    Transform = report$transform,
    Crosstalk = report$crosstalk,
    Block = report$block_label,
    Mean = report$mean,
    SD = report$sd,
    stringsAsFactors = FALSE
  )
  header <- as.data.frame(as.list(names(body)), stringsAsFactors = FALSE)
  names(header) <- names(body)
  table_data <- rbind(header, body)

  long <- do.call(rbind, lapply(seq_len(nrow(table_data)), function(i) {
    data.frame(
      row = i,
      col = seq_along(table_data),
      label = unname(unlist(table_data[i, ], use.names = FALSE)),
      header = i == 1L,
      stringsAsFactors = FALSE
    )
  }))
  long$row_plot <- max(long$row) - long$row + 1L

  ggplot(long, aes(x = col, y = row_plot)) +
    geom_tile(aes(fill = header), color = "grey75", linewidth = 0.25) +
    geom_text(aes(label = label, fontface = ifelse(header, "bold", "plain")), size = 3.1) +
    scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "grey90"), guide = "none") +
    scale_x_continuous(expand = c(0, 0), limits = c(0.5, ncol(table_data) + 0.5)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0.5, nrow(table_data) + 0.5)) +
    labs(
      title = "Rescaled Frobenius reconstruction error at calibrated crosstalk",
      subtitle = "Lower values indicate closer empirical-vs-simulated correlation blocks",
      x = NULL,
      y = NULL
    ) +
    theme_void(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0, size = 11),
      plot.subtitle = element_text(hjust = 0, size = 9),
      plot.margin = margin(4, 8, 4, 8)
    )
}

save_full_block_figure_with_table <- function(plot, table_plot, filename, width = 10, height = 12) {
  grDevices::pdf(filename, width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grid.newpage()
  print(plot, vp = grid::viewport(x = 0.5, y = 0.64, width = 1, height = 0.72))
  print(table_plot, vp = grid::viewport(x = 0.5, y = 0.13, width = 0.92, height = 0.22))
}

template_path <- path.expand(get_arg(
  "template_path",
  "inst/extdata/HuMMANet_FranzosaE_2019_template.rds"
))
out_dir <- path.expand(get_arg("out_dir", "~/Mimesys_validation/franzosa_local_block_reconstruction"))
feature_counts <- parse_counts(get_arg("feature_counts", "microbiome=50,metabolite=50"))
nsample <- as.integer(get_arg("nsample", "100"))
nrep <- as.integer(get_arg("nrep", "1"))
crosstalk_grid <- parse_numeric_grid(get_arg("crosstalk_grid", get_arg("crosstalk", "1")))
seed <- as.integer(get_arg("seed", "1"))
dependence <- get_arg("dependence", "factor")
cor_shrinkage <- as.numeric(get_arg("cor_shrinkage", "0.05"))
module_count_arg <- get_arg("module_count", "")
module_count <- if (nzchar(module_count_arg)) as.integer(module_count_arg) else NULL
module_between_shrinkage <- as.numeric(get_arg("module_between_shrinkage", "0.75"))
view_models <- parse_key_values(get_arg("view_models", ""))
scenarios <- split_arg(get_arg("scenarios", "raw,clrlog"), ",")
if (!all(scenarios %in% c("raw", "clrlog"))) {
  stop("--scenarios must contain only raw and/or clrlog.", call. = FALSE)
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
obj <- readRDS(template_path)
template <- if (is.list(obj) && !is.null(obj$template)) obj$template else obj
template <- as_mimesys_view_list(template)

res <- lapply(scenarios, run_scenario,
  template = template,
  caps = feature_counts,
  nsample = nsample,
  nrep = nrep,
  crosstalk_grid = crosstalk_grid,
  seed = seed,
  dependence = dependence,
  cor_shrinkage = cor_shrinkage,
  module_count = module_count,
  module_between_shrinkage = module_between_shrinkage,
  view_models = view_models
)
names(res) <- scenarios

frobenius_raw <- do.call(rbind, lapply(res, `[[`, "frobenius"))
frobenius_summary <- summarize_frobenius(frobenius_raw)
crosstalk_calibration <- summarize_crosstalk_calibration(frobenius_raw)
best_crosstalk <- best_crosstalk_table(crosstalk_calibration)
pi_frobenius_report <- make_pi_frobenius_report(frobenius_summary, best_crosstalk)

write.csv(frobenius_raw, file.path(out_dir, "franzosa_rescaled_frobenius_raw.csv"), row.names = FALSE)
write.csv(frobenius_summary, file.path(out_dir, "franzosa_rescaled_frobenius_summary.csv"), row.names = FALSE)
write.csv(crosstalk_calibration, file.path(out_dir, "franzosa_crosstalk_calibration_summary.csv"), row.names = FALSE)
write.csv(best_crosstalk, file.path(out_dir, "franzosa_best_crosstalk.csv"), row.names = FALSE)
write.csv(pi_frobenius_report, file.path(out_dir, "franzosa_pi_frobenius_report.csv"), row.names = FALSE)
plot_crosstalk_calibration(crosstalk_calibration, best_crosstalk, out_dir, feature_counts, nsample)
write_readme(
  out_dir, best_crosstalk, feature_counts, nsample, nrep,
  crosstalk_grid, dependence, cor_shrinkage,
  module_count, module_between_shrinkage, view_models
)

heatmap_df <- do.call(rbind, lapply(res, `[[`, "heatmaps"))
heatmap_df$source <- factor(heatmap_df$source, levels = c("Empirical", "MimeSys"))
scenario_labels <- c(raw = "No transform", clrlog = "CLR microbiome + log metabolite")
heatmap_df$scenario <- factor(heatmap_df$scenario, levels = scenarios, labels = unname(scenario_labels[scenarios]))

p <- ggplot(heatmap_df, aes(x = order_x, y = order_y, fill = correlation)) +
  geom_tile() +
  facet_grid(scenario + view ~ source, scales = "free", space = "free") +
  coord_equal() +
  scale_fill_gradient2(
    low = "#B2182B",
    mid = "white",
    high = "#1B7837",
    midpoint = 0,
    limits = c(-1, 1),
    na.value = "white"
  ) +
  labs(
    title = "FranzosaE 2019 feature-level correlation reconstruction",
    subtitle = paste0("nsample=", nsample, ", features=", paste(paste(names(feature_counts), feature_counts, sep = "="), collapse = ", ")),
    x = "Features",
    y = "Features",
    fill = "cor"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    strip.background = element_rect(fill = "white", color = "grey70"),
    plot.title = element_text(face = "bold")
  )

ggsave(
  filename = file.path(out_dir, "franzosa_feature_correlation_reconstruction.pdf"),
  plot = p,
  width = 10,
  height = 9,
  device = grDevices::pdf
)

full_heatmap_df <- do.call(rbind, lapply(res, `[[`, "full_heatmaps"))
full_heatmap_df$source <- factor(full_heatmap_df$source, levels = c("Empirical", "MimeSys"))
full_heatmap_df$scenario <- factor(full_heatmap_df$scenario, levels = scenarios, labels = unname(scenario_labels[scenarios]))
idx <- res[[1]]$idx
boundary_df <- block_boundary_df(idx)
label_df <- block_label_df(idx)

p_full <- ggplot(full_heatmap_df, aes(x = order_x, y = order_y, fill = correlation)) +
  geom_tile() +
  geom_vline(data = boundary_df, aes(xintercept = boundary), inherit.aes = FALSE, linewidth = 0.25) +
  geom_hline(data = boundary_df, aes(yintercept = boundary), inherit.aes = FALSE, linewidth = 0.25) +
  facet_grid(scenario ~ source, scales = "free", space = "free") +
  coord_equal() +
  scale_x_continuous(breaks = label_df$pos, labels = label_df$view, expand = c(0, 0)) +
  scale_y_continuous(breaks = label_df$pos, labels = label_df$view, expand = c(0, 0)) +
  scale_fill_gradient2(
    low = "#B2182B",
    mid = "white",
    high = "#1B7837",
    midpoint = 0,
    limits = c(-1, 1),
    na.value = "white"
  ) +
  labs(
    title = "FranzosaE 2019 full block feature-correlation reconstruction",
    subtitle = paste0("nsample=", nsample, ", features=", paste(paste(names(feature_counts), feature_counts, sep = "="), collapse = ", ")),
    x = "Features",
    y = "Features",
    fill = "cor"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    strip.background = element_rect(fill = "white", color = "grey70"),
    plot.title = element_text(face = "bold")
  )

p_frobenius_table <- plot_pi_frobenius_table(pi_frobenius_report)
save_full_block_figure_with_table(
  plot = p_full,
  table_plot = p_frobenius_table,
  filename = file.path(out_dir, "franzosa_full_block_feature_correlation_reconstruction.pdf"),
  width = 10,
  height = 12
)

saveRDS(
  list(
    frobenius_raw = frobenius_raw,
    frobenius_summary = frobenius_summary,
    crosstalk_calibration = crosstalk_calibration,
    best_crosstalk = best_crosstalk,
    pi_frobenius_report = pi_frobenius_report,
    heatmap_data = heatmap_df,
    full_heatmap_data = full_heatmap_df
  ),
  file.path(out_dir, "franzosa_block_correlation_reconstruction.rds")
)

cat("Saved outputs to:\n", out_dir, "\n", sep = "")
cat("\nBest crosstalk by transform:\n")
print(best_crosstalk)
cat("\nPI-ready Frobenius report:\n")
print(pi_frobenius_report)
cat("\nFrobenius summary:\n")
print(frobenius_summary)
