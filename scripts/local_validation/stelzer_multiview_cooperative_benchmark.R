#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) == 0L) return(default)
  sub(paste0("^--", name, "="), "", hit[[1]])
}

split_arg <- function(x, sep = ",") {
  if (is.null(x) || !nzchar(x)) return(character())
  out <- trimws(strsplit(x, sep, fixed = TRUE)[[1]])
  out[nzchar(out)]
}

load_mimesys <- function() {
  project_dir <- getwd()
  if (requireNamespace("devtools", quietly = TRUE) && file.exists(file.path(project_dir, "DESCRIPTION"))) {
    devtools::load_all(project_dir, quiet = TRUE)
  } else {
    suppressPackageStartupMessages(library(MimeSys))
  }
}

keep_top_variable <- function(x, n_features) {
  vars <- apply(x, 2, stats::var, na.rm = TRUE)
  vars <- vars[is.finite(vars)]
  vars <- sort(vars, decreasing = TRUE)
  n_keep <- min(as.integer(n_features), length(vars))
  if (n_keep <= 0L) stop("No variable features available.", call. = FALSE)
  x[, names(vars)[seq_len(n_keep)], drop = FALSE]
}

load_stelzer_views <- function(data_path) {
  env <- new.env(parent = emptyenv())
  load(data_path, envir = env)
  if (!exists("pcl", envir = env, inherits = FALSE)) {
    stop("Expected object `pcl` in ", data_path, call. = FALSE)
  }
  pcl <- get("pcl", envir = env, inherits = FALSE)
  feature_table <- pcl$feature_table
  sample_metadata <- pcl$sample_metadata
  feature_metadata <- pcl$feature_metadata
  subject_id <- sample_metadata$subjectID

  subject_feature <- t(apply(feature_table, 1, function(x) {
    tapply(x, subject_id, mean, na.rm = TRUE)
  }))

  feature_type <- feature_metadata$featureType
  names(feature_type) <- rownames(feature_metadata)

  list(
    Immunome = t(subject_feature[feature_type == "CYTOF", , drop = FALSE]),
    Metabolome = t(subject_feature[feature_type == "Metabolomics", , drop = FALSE]),
    Proteome = t(subject_feature[feature_type == "Proteomics", , drop = FALSE])
  )
}

prep_view_matrices <- function(train_x, test_x) {
  means <- colMeans(train_x, na.rm = TRUE)
  sds <- apply(train_x, 2, stats::sd, na.rm = TRUE)
  sds[!is.finite(sds) | sds == 0] <- 1

  train_scaled <- sweep(train_x, 2, means, "-")
  train_scaled <- sweep(train_scaled, 2, sds, "/")
  test_scaled <- sweep(test_x, 2, means, "-")
  test_scaled <- sweep(test_scaled, 2, sds, "/")

  train_scaled[!is.finite(train_scaled)] <- 0
  test_scaled[!is.finite(test_scaled)] <- 0
  list(train = train_scaled, test = test_scaled)
}

prep_views <- function(train_views, test_views) {
  out <- lapply(names(train_views), function(view) {
    prep_view_matrices(train_views[[view]], test_views[[view]])
  })
  names(out) <- names(train_views)
  list(
    train = lapply(out, `[[`, "train"),
    test = lapply(out, `[[`, "test")
  )
}

r2_score <- function(y, pred) {
  ok <- is.finite(y) & is.finite(pred)
  y <- y[ok]
  pred <- pred[ok]
  if (length(y) < 2L) return(NA_real_)
  denom <- sum((y - mean(y))^2)
  if (!is.finite(denom) || denom <= 0) return(NA_real_)
  1 - sum((y - pred)^2) / denom
}

fit_cv_multiview <- function(x_train, y_train, rho, foldid, type_measure = "mse", family = stats::gaussian()) {
  multiview::cv.multiview(
    x_list = x_train,
    y = y_train,
    rho = rho,
    family = family,
    type.measure = type_measure,
    nfolds = length(unique(foldid)),
    foldid = foldid
  )
}

cv_error_min <- function(cvfit) {
  suppressWarnings(min(cvfit$cvm, na.rm = TRUE))
}

fit_eval_fixed_rho <- function(x_train, y_train, x_test, y_test, rho, foldid, method_label) {
  fit <- fit_cv_multiview(x_train, y_train, rho = rho, foldid = foldid)
  pred <- as.numeric(stats::predict(fit, newx = x_test, s = "lambda.min"))
  data.frame(
    method = method_label,
    rho = rho,
    selected_rho = rho,
    lambda_min = fit$lambda.min,
    cv_mse = cv_error_min(fit),
    test_r2 = r2_score(y_test, pred),
    stringsAsFactors = FALSE
  )
}

fit_eval_cv_rho <- function(x_train, y_train, x_test, y_test, rho_grid, foldid) {
  fits <- lapply(rho_grid, function(rho) {
    tryCatch(
      fit_cv_multiview(x_train, y_train, rho = rho, foldid = foldid),
      error = function(e) e
    )
  })

  cv_mse <- vapply(fits, function(fit) {
    if (inherits(fit, "error")) return(Inf)
    cv_error_min(fit)
  }, numeric(1))

  if (all(!is.finite(cv_mse))) {
    return(data.frame(
      method = "rho_cv",
      rho = NA_real_,
      selected_rho = NA_real_,
      lambda_min = NA_real_,
      cv_mse = NA_real_,
      test_r2 = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  best_idx <- which.min(cv_mse)
  best_fit <- fits[[best_idx]]
  pred <- as.numeric(stats::predict(best_fit, newx = x_test, s = "lambda.min"))
  data.frame(
    method = "rho_cv",
    rho = NA_real_,
    selected_rho = rho_grid[[best_idx]],
    lambda_min = best_fit$lambda.min,
    cv_mse = cv_mse[[best_idx]],
    test_r2 = r2_score(y_test, pred),
    stringsAsFactors = FALSE
  )
}

make_foldid <- function(n, folds) {
  sample(rep(seq_len(folds), length.out = n))
}

summarize_results <- function(raw) {
  suppressPackageStartupMessages(library(dplyr))
  raw |>
    dplyr::group_by(
      dataset, scenario, feature_n, nsample, crosstalk, method
    ) |>
    dplyr::summarise(
      performance_mean = mean(test_r2, na.rm = TRUE),
      performance_sd = stats::sd(test_r2, na.rm = TRUE),
      selected_rho_mean = mean(selected_rho, na.rm = TRUE),
      selected_rho_sd = stats::sd(selected_rho, na.rm = TRUE),
      n_success = sum(is.finite(test_r2)),
      n_total = dplyr::n(),
      .groups = "drop"
    )
}

plot_performance_pages <- function(summary, out_file) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
  })

  scenario_order <- intersect(c("balanced_aligned", "balanced_random"), unique(as.character(summary$scenario)))
  scenario_labels <- c(
    balanced_aligned = "Balanced aligned signal",
    balanced_random = "Balanced random signal"
  )

  summary <- summary |>
    dplyr::mutate(
      scenario = factor(scenario, levels = scenario_order),
      feature_n = factor(feature_n, levels = c(100, 200, 500)),
      nsample = factor(nsample, levels = c(53, 100, 200, 500)),
      crosstalk = factor(crosstalk, levels = c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1)),
      method = factor(method, levels = c("rho_0", "rho_1", "rho_cv"),
                      labels = c("rho = 0", "rho = 1", "rho selected by CV"))
    )

  clean_theme <- ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 18),
      axis.title = ggplot2::element_text(face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, color = "black"),
      axis.text.y = ggplot2::element_text(color = "black"),
      strip.background = ggplot2::element_rect(fill = "white", color = "black", linewidth = 0.4),
      strip.text = ggplot2::element_text(face = "bold", color = "black"),
      panel.grid = ggplot2::element_blank(),
      legend.position = "bottom"
    )

  grDevices::pdf(out_file, width = 13, height = 9, bg = "white")
  on.exit(grDevices::dev.off(), add = TRUE)

  for (scen in scenario_order) {
    p <- summary |>
      dplyr::filter(as.character(scenario) == scen) |>
      ggplot2::ggplot(ggplot2::aes(x = crosstalk, y = performance_mean, fill = method)) +
      ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.75), width = 0.68, color = "white", linewidth = 0.15) +
      ggplot2::facet_grid(nsample ~ feature_n, scales = "free_y") +
      ggplot2::scale_fill_manual(values = c("rho = 0" = "#7A7A7A", "rho = 1" = "#00A9B7", "rho selected by CV" = "#3B6FB6")) +
      ggplot2::labs(
        title = scenario_labels[[scen]],
        x = "Crosstalk",
        y = expression("Mean test " * R^2),
        fill = NULL
      ) +
      clean_theme
    print(p)
  }
}

plot_rho_pages <- function(summary, out_file) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
  })

  rho_summary <- summary |>
    dplyr::filter(method == "rho_cv") |>
    dplyr::mutate(
      feature_n = factor(feature_n, levels = c(100, 200, 500)),
      nsample = factor(nsample, levels = c(53, 100, 200, 500)),
      crosstalk = factor(crosstalk, levels = c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1))
    )

  scenario_order <- intersect(c("balanced_aligned", "balanced_random"), unique(as.character(rho_summary$scenario)))
  scenario_labels <- c(
    balanced_aligned = "Balanced aligned signal",
    balanced_random = "Balanced random signal"
  )

  clean_theme <- ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 18),
      axis.title = ggplot2::element_text(face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, color = "black"),
      axis.text.y = ggplot2::element_text(color = "black"),
      strip.background = ggplot2::element_rect(fill = "white", color = "black", linewidth = 0.4),
      strip.text = ggplot2::element_text(face = "bold", color = "black"),
      panel.grid = ggplot2::element_blank(),
      legend.position = "none"
    )

  grDevices::pdf(out_file, width = 13, height = 9, bg = "white")
  on.exit(grDevices::dev.off(), add = TRUE)

  for (scen in scenario_order) {
    p <- rho_summary |>
      dplyr::filter(scenario == scen) |>
      ggplot2::ggplot(ggplot2::aes(x = crosstalk, y = selected_rho_mean, group = 1)) +
      ggplot2::geom_hline(yintercept = c(0, 1), linetype = "dashed", linewidth = 0.4, color = "grey50") +
      ggplot2::geom_line(color = "#3B6FB6", linewidth = 0.7) +
      ggplot2::geom_point(color = "#3B6FB6", size = 2) +
      ggplot2::facet_grid(nsample ~ feature_n) +
      ggplot2::scale_y_continuous(limits = c(0, max(1, max(rho_summary$selected_rho_mean, na.rm = TRUE)))) +
      ggplot2::labs(
        title = scenario_labels[[scen]],
        x = "Crosstalk",
        y = "Mean selected rho"
      ) +
      clean_theme
    print(p)
  }
}

main <- function() {
  load_mimesys()
  if (!requireNamespace("multiview", quietly = TRUE)) {
    stop("This benchmark requires the multiview package.", call. = FALSE)
  }

  default_data_path <- file.path(getwd(), "inst/extdata/StelzerDOS.RData")
  if (!file.exists(default_data_path)) {
    default_data_path <- "~/Mimesys_2.0/inst/extdata/StelzerDOS.RData"
  }
  data_path <- path.expand(get_arg("data_path", default_data_path))
  out_dir <- path.expand(get_arg("out_dir", "~/Mimesys_validation/stelzer_multiview_cooperative_benchmark"))
  feature_sizes <- as.integer(split_arg(get_arg("feature_sizes", "100,200,500")))
  sample_sizes <- as.integer(split_arg(get_arg("sample_sizes", "53,100,200,500")))
  crosstalk_grid <- as.numeric(split_arg(get_arg("crosstalk", "0,0.1,0.25,0.5,0.75,0.9,1")))
  rho_grid <- as.numeric(split_arg(get_arg("rho_grid", "0,0.1,0.25,0.5,0.75,1")))
  scenarios <- split_arg(get_arg("scenarios", "balanced_aligned,balanced_random"))
  nrep <- as.integer(get_arg("nrep", "20"))
  folds <- as.integer(get_arg("folds", "5"))
  seed <- as.integer(get_arg("seed", "30000"))
  snr <- as.numeric(get_arg("snr", "0.5"))
  de_prob <- as.numeric(get_arg("de_prob", "0.02"))
  de_fac_loc <- as.numeric(get_arg("de_fac_loc", "0.8"))
  de_fac_scale <- as.numeric(get_arg("de_fac_scale", "0.2"))
  dependence <- get_arg("dependence", "factor")

  if (any(is.na(feature_sizes)) || any(feature_sizes <= 0)) stop("--feature_sizes must be positive integers.", call. = FALSE)
  if (any(is.na(sample_sizes)) || any(sample_sizes <= 0)) stop("--sample_sizes must be positive integers.", call. = FALSE)
  if (any(is.na(crosstalk_grid))) stop("--crosstalk must be numeric.", call. = FALSE)
  if (any(is.na(rho_grid))) stop("--rho_grid must be numeric.", call. = FALSE)
  if (is.na(nrep) || nrep <= 0) stop("--nrep must be positive.", call. = FALSE)
  if (folds < 3) stop("--folds must be at least 3 for multiview::cv.multiview().", call. = FALSE)

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  views_full <- load_stelzer_views(data_path)
  scenario_defs <- list(
    balanced_aligned = list(signal_share = c(Immunome = 1, Metabolome = 1, Proteome = 1), de_align = 1),
    balanced_random = list(signal_share = c(Immunome = 1, Metabolome = 1, Proteome = 1), de_align = NULL)
  )
  unknown <- setdiff(scenarios, names(scenario_defs))
  if (length(unknown) > 0L) stop("Unknown scenario(s): ", paste(unknown, collapse = ", "), call. = FALSE)

  rows <- list()
  counter <- 1L
  for (feature_n in feature_sizes) {
    views <- lapply(views_full, keep_top_variable, n_features = feature_n)
    tmpl <- fit_template(views)
    for (nsample in sample_sizes) {
      for (crosstalk in crosstalk_grid) {
        for (scenario in scenarios) {
          scen <- scenario_defs[[scenario]]
          for (rep_idx in seq_len(nrep)) {
            this_seed <- seed + counter
            set.seed(this_seed)
            sim <- gen_simmba(
              nsample = nsample,
              template = tmpl,
              crosstalk = crosstalk,
              outcome.type = "continuous",
              de.align = scen$de_align,
              signal.share = scen$signal_share,
              dependence = dependence,
              de.prob = de_prob,
              de.facLoc = rep(de_fac_loc, length(tmpl$views)),
              de.facScale = rep(de_fac_scale, length(tmpl$views)),
              snr = snr,
              nrep = 1,
              seed = this_seed
            )

            x_train <- simulation_view_matrices(sim, rep = 1, split = "train", feature_names = "template_feature")
            x_test <- simulation_view_matrices(sim, rep = 1, split = "test", feature_names = "template_feature")
            y_train <- sim$trainDat[[1]]$sample_metadata$Y
            y_test <- sim$testDat[[1]]$sample_metadata$Y
            x_prepped <- prep_views(x_train, x_test)
            foldid <- make_foldid(length(y_train), folds)

            out <- tryCatch(
              rbind(
                fit_eval_fixed_rho(x_prepped$train, y_train, x_prepped$test, y_test, rho = 0, foldid = foldid, method_label = "rho_0"),
                fit_eval_fixed_rho(x_prepped$train, y_train, x_prepped$test, y_test, rho = 1, foldid = foldid, method_label = "rho_1"),
                fit_eval_cv_rho(x_prepped$train, y_train, x_prepped$test, y_test, rho_grid = rho_grid, foldid = foldid)
              ),
              error = function(e) {
                data.frame(
                  method = c("rho_0", "rho_1", "rho_cv"),
                  rho = c(0, 1, NA_real_),
                  selected_rho = c(0, 1, NA_real_),
                  lambda_min = NA_real_,
                  cv_mse = NA_real_,
                  test_r2 = NA_real_,
                  error = conditionMessage(e),
                  stringsAsFactors = FALSE
                )
              }
            )
            if (!"error" %in% names(out)) out$error <- NA_character_
            meta <- data.frame(
              dataset = "StelzerDOS",
              scenario = scenario,
              feature_n = feature_n,
              nsample = nsample,
              crosstalk = crosstalk,
              rep = rep_idx,
              seed = this_seed,
              folds = folds,
              snr = snr,
              de_prob = de_prob,
              dependence = dependence,
              stringsAsFactors = FALSE
            )
            rows[[counter]] <- cbind(meta[rep(1, nrow(out)), , drop = FALSE], out)
            message("Done: feature_n=", feature_n, ", nsample=", nsample,
                    ", crosstalk=", crosstalk, ", scenario=", scenario,
                    ", rep=", rep_idx)
            counter <- counter + 1L
          }
        }
      }
    }
  }

  raw <- do.call(rbind, rows)
  raw_file <- file.path(out_dir, "stelzer_multiview_cooperative_benchmark_raw.csv")
  summary_file <- file.path(out_dir, "stelzer_multiview_cooperative_benchmark_summary.csv")
  perf_pdf <- file.path(out_dir, "stelzer_multiview_performance_by_scenario.pdf")
  rho_pdf <- file.path(out_dir, "stelzer_multiview_selected_rho_by_scenario.pdf")

  utils::write.csv(raw, raw_file, row.names = FALSE)
  summary <- summarize_results(raw)
  utils::write.csv(summary, summary_file, row.names = FALSE)
  plot_performance_pages(summary, perf_pdf)
  plot_rho_pages(summary, rho_pdf)

  cat("Raw:", raw_file, "\n")
  cat("Summary:", summary_file, "\n")
  cat("Performance PDF:", perf_pdf, "\n")
  cat("Selected rho PDF:", rho_pdf, "\n")
}

main()
