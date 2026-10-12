test_that("template simulator preserves requested structure", {
  set.seed(1)
  n <- 60
  latent <- matrix(stats::rnorm(n * 2), n, 2)
  template <- list(
    mrna = latent %*% matrix(stats::rnorm(2 * 30), 2, 30) + matrix(stats::rnorm(n * 30), n, 30),
    protein = latent %*% matrix(stats::rnorm(2 * 20), 2, 20) + matrix(stats::rnorm(n * 20), n, 20)
  )
  rownames(template$mrna) <- rownames(template$protein) <- paste0("S", seq_len(n))

  tm <- fit_template(template, rank = 3, shared.rank = 1)
  sim <- gen_simmba(
    nsample = 80,
    template = tm,
    p = c(mrna = 25, protein = 15),
    crosstalk = 0.5,
    de.align = 1,
    nrep = 2,
    seed = 11
  )

  expect_s3_class(sim$template, "mimesys_template")
  expect_length(sim$trainDat, 2)
  expect_equal(nrow(sim$trainDat[[1]]$feature_table), 40)
  expect_true("beta" %in% names(sim$trainDat[[1]]$feature_metadata))
  expect_true("Y" %in% names(sim$trainDat[[1]]$sample_metadata))
})

test_that("feature counts can be scalar, positional, or named by view", {
  set.seed(10)
  n <- 50
  template <- list(
    mrna = matrix(stats::rnorm(n * 12), n, 12),
    protein = matrix(stats::rnorm(n * 10), n, 10)
  )
  rownames(template$mrna) <- rownames(template$protein) <- paste0("S", seq_len(n))
  tm <- fit_template(template, rank = 2, shared.rank = 1)

  scalar <- sim_from_template(tm, n = 30, p = 6)
  positional <- sim_from_template(tm, n = 30, p = c(5, 7))
  named <- sim_from_template(tm, n = 30, p = c(protein = 4, mrna = 8))

  scalar_counts <- table(scalar$feature_metadata$featureType)
  positional_counts <- table(positional$feature_metadata$featureType)
  named_counts <- table(named$feature_metadata$featureType)

  expect_equal(stats::setNames(as.integer(scalar_counts), names(scalar_counts)), c(mrna = 6L, protein = 6L))
  expect_equal(stats::setNames(as.integer(positional_counts), names(positional_counts)), c(mrna = 5L, protein = 7L))
  expect_equal(stats::setNames(as.integer(named_counts), names(named_counts)), c(mrna = 8L, protein = 4L))
  expect_error(sim_from_template(tm, n = 30, p = c(mrna = 5)), "missing view")
})

test_that("Friedman outcomes can control nonlinear signal allocation by view", {
  set.seed(12)
  n <- 50
  template <- list(
    mrna = matrix(stats::rnorm(n * 8), n, 8),
    protein = matrix(stats::rnorm(n * 8), n, 8)
  )
  rownames(template$mrna) <- rownames(template$protein) <- paste0("S", seq_len(n))

  tm <- fit_template(template, rank = 2, shared.rank = 1)
  sim <- gen_simmba(
    nsample = 60,
    template = tm,
    p = c(mrna = 8, protein = 8),
    ygen.mode = "Friedman",
    de.prob = 1,
    friedman.share = c(mrna = 1, protein = 0),
    nrep = 1,
    seed = 121
  )

  fm <- sim$trainDat[[1]]$feature_metadata
  expect_equal(sum(fm$friedman_term), 5)
  expect_true(all(fm$featureType[fm$friedman_term] == "mrna"))
  expect_equal(sim$friedman.share, c(mrna = 1, protein = 0))
})

test_that("Friedman inputs are rank transformed before nonlinear outcome generation", {
  x <- cbind(
    a = c(10, 20, 30, 40),
    b = c(1000, 100, 10, 1),
    c = c(5, 5, 5, 9),
    d = c(NA, 2, 4, 8),
    e = c(3, 3, 3, 3)
  )

  z <- .rank_to_unit_interval(x)

  expect_true(all(is.finite(z)))
  expect_true(all(z > 0 & z < 1))
  expect_equal(z[, "a"], c(0.2, 0.4, 0.6, 0.8))
  expect_equal(z[, "b"], c(0.8, 0.6, 0.4, 0.2))
  expect_equal(z[, "e"], rep(0.5, 4))
})

test_that("shared latent outcome mode generates finite shared signals", {
  set.seed(13)
  n <- 70
  shared <- stats::rnorm(n)
  template <- list(
    mrna = outer(shared, stats::rnorm(16)) + matrix(stats::rnorm(n * 16), n, 16),
    protein = outer(shared, stats::rnorm(16)) + matrix(stats::rnorm(n * 16), n, 16)
  )
  rownames(template$mrna) <- rownames(template$protein) <- paste0("S", seq_len(n))

  tm <- fit_template(template, rank = 3, shared.rank = 1)
  late_sim <- gen_simmba(
    nsample = 80,
    template = tm,
    p = c(mrna = 12, protein = 12),
    ygen.mode = "SharedLatent",
    crosstalk = 1,
    nrep = 1,
    seed = 131
  )

  expect_true(all(is.finite(late_sim$trainDat[[1]]$sample_metadata$Y)))
  expect_true("Xbeta_shared" %in% names(late_sim$trainDat[[1]]$sample_metadata))
  expect_equal(nrow(late_sim$trainDat[[1]]$truth$shared), nrow(late_sim$trainDat[[1]]$sample_metadata))
})

test_that("shared latent outcome uses all available shared factors", {
  set.seed(14)
  n <- 90
  shared <- matrix(stats::rnorm(n * 2), n, 2)
  template <- list(
    mrna = shared %*% matrix(stats::rnorm(2 * 18), 2, 18) + matrix(stats::rnorm(n * 18), n, 18),
    protein = shared %*% matrix(stats::rnorm(2 * 18), 2, 18) + matrix(stats::rnorm(n * 18), n, 18)
  )
  rownames(template$mrna) <- rownames(template$protein) <- paste0("S", seq_len(n))

  tm <- fit_template(template, rank = 4, shared.rank = 2)
  sim <- gen_simmba(
    nsample = 80,
    template = tm,
    p = c(mrna = 12, protein = 12),
    ygen.mode = "SharedLatent",
    crosstalk = 1,
    nrep = 1,
    seed = 141
  )

  shared_train <- sim$trainDat[[1]]$truth$shared
  shared_test <- sim$testDat[[1]]$truth$shared
  shared_all <- rbind(shared_train, shared_test)
  id_all <- c(
    sim$trainDat[[1]]$sample_metadata$subjectID,
    sim$testDat[[1]]$sample_metadata$subjectID
  )
  expected_all <- as.numeric(scale(rowSums(scale(shared_all)) / sqrt(ncol(shared_all))))
  names(expected_all) <- id_all
  expected <- unname(expected_all[sim$trainDat[[1]]$sample_metadata$subjectID])

  expect_equal(ncol(shared_train), 2)
  expect_equal(sim$trainDat[[1]]$sample_metadata$Xbeta_shared, expected, tolerance = 1e-8)
  expect_false(isTRUE(all.equal(expected, as.numeric(scale(shared_train[, 1])), tolerance = 1e-8)))
})

test_that("template diagnostics summarize signal and crosstalk", {
  set.seed(2)
  n <- 70
  shared <- stats::rnorm(n)
  template <- list(
    mrna = outer(shared, stats::rnorm(24)) + matrix(stats::rnorm(n * 24), n, 24),
    protein = outer(shared, stats::rnorm(18)) + matrix(stats::rnorm(n * 18), n, 18)
  )
  rownames(template$mrna) <- rownames(template$protein) <- paste0("S", seq_len(n))

  tm <- fit_template(template, rank = 3, shared.rank = 1)
  sim <- gen_simmba(
    nsample = 100,
    template = tm,
    p = c(mrna = 20, protein = 16),
    crosstalk = 1,
    de.align = 1,
    de.prob = 0.2,
    nrep = 3,
    seed = 22
  )
  diagnostics <- validate_template_simulation(sim)
  summary <- summarize_template_diagnostics(diagnostics)

  expect_equal(nrow(diagnostics), 3)
  expect_equal(nrow(summary), 1)
  expect_true(all(c("mean_abs_feature_cor", "signal_agreement", "y_signal_cor") %in% names(summary)))
  expect_true(is.finite(summary$mean_abs_feature_cor))
})

test_that("simulation summaries report key settings", {
  set.seed(23)
  n <- 50
  template <- list(
    mrna = matrix(stats::rnorm(n * 10), n, 10),
    protein = matrix(stats::rnorm(n * 10), n, 10)
  )
  rownames(template$mrna) <- rownames(template$protein) <- paste0("S", seq_len(n))

  sim <- gen_simmba(
    nsample = 60,
    template = fit_template(template, rank = 2, shared.rank = 1),
    p = 6,
    crosstalk = 0.25,
    nrep = 2,
    seed = 230
  )
  sm <- summarize_simulation(sim)

  expect_s3_class(sm, "mimesys_simulation_summary")
  expect_equal(sm$nrep, 2)
  expect_equal(sm$feature_counts, c(mrna = 6L, protein = 6L))
  expect_equal(sm$crosstalk, 0.25)
  expect_output(print(sm), "MimeSys simulation")
})

test_that("crosstalk increases observed cross-view dependence in template simulations", {
  set.seed(3)
  n <- 80
  shared <- matrix(stats::rnorm(n * 2), n, 2)
  template <- list(
    mrna = shared %*% matrix(stats::rnorm(2 * 35), 2, 35) + matrix(stats::rnorm(n * 35), n, 35),
    protein = shared %*% matrix(stats::rnorm(2 * 30), 2, 30) + matrix(stats::rnorm(n * 30), n, 30)
  )
  rownames(template$mrna) <- rownames(template$protein) <- paste0("S", seq_len(n))
  tm <- fit_template(template, rank = 4, shared.rank = 1)

  low <- sim_from_template(tm, n = 120, p = c(mrna = 25, protein = 25), crosstalk = 0)
  high <- sim_from_template(tm, n = 120, p = c(mrna = 25, protein = 25), crosstalk = 1)

  x_low <- t(as.matrix(low$feature_table))
  x_high <- t(as.matrix(high$feature_table))
  low_type <- low$feature_metadata$featureType
  high_type <- high$feature_metadata$featureType

  low_cor <- mean(abs(stats::cor(x_low[, low_type == "mrna"], x_low[, low_type == "protein"])))
  high_cor <- mean(abs(stats::cor(x_high[, high_type == "mrna"], x_high[, high_type == "protein"])))

  expect_gt(high_cor, low_cor)
})

test_that("crosstalk calibration plot returns a ggplot", {
  testthat::skip_if_not_installed("ggplot2")
  calibration <- list(
    best_crosstalk = 0.5,
    summary = data.frame(
      crosstalk = c(0, 0.5, 1, 0, 0.5, 1),
      block = rep(c("cor(mrna,protein)", "cor(mrna,mirna)"), each = 3),
      view_x = rep(c("mrna", "mrna"), each = 3),
      view_y = rep(c("protein", "mirna"), each = 3),
      diagonal_block = FALSE,
      mean = c(0.2, 0.1, 0.15, 0.3, 0.12, 0.2),
      sd = c(0.02, 0.01, 0.03, 0.02, 0.02, 0.03),
      stringsAsFactors = FALSE
    )
  )

  expect_s3_class(plot_crosstalk_calibration(calibration), "ggplot")
})

test_that("two-part copula view model preserves sparse compositional microbiome output", {
  set.seed(31)
  n <- 70
  shared <- stats::rnorm(n)
  microbe_latent <- outer(shared, stats::rnorm(18)) + matrix(stats::rnorm(n * 18), n, 18)
  microbiome <- exp(microbe_latent)
  microbiome[microbe_latent < 0.25] <- 0
  microbiome <- sweep(microbiome, 1, pmax(rowSums(microbiome), 1), "/")
  metabolite <- outer(shared, stats::rnorm(12)) + matrix(stats::rnorm(n * 12), n, 12)
  template <- list(microbiome = microbiome, metabolite = metabolite)
  rownames(template$microbiome) <- rownames(template$metabolite) <- paste0("S", seq_len(n))

  tm <- fit_template(template, rank = 3, shared.rank = 1)
  sim <- sim_from_template(
    tm,
    n = 50,
    p = c(microbiome = 12, metabolite = 10),
    crosstalk = 0.5,
    dependence = "empirical_copula",
    view.models = c(microbiome = "two_part_copula")
  )
  x <- t(as.matrix(sim$feature_table))
  microbe <- x[, sim$feature_metadata$featureType == "microbiome", drop = FALSE]
  metabolite <- x[, sim$feature_metadata$featureType == "metabolite", drop = FALSE]

  expect_true(all(sim$feature_metadata$view_model[sim$feature_metadata$featureType == "microbiome"] == "two_part_copula"))
  expect_true(mean(microbe == 0) > 0)
  expect_equal(as.numeric(rowSums(microbe)), rep(1, nrow(microbe)), tolerance = 1e-8)
  expect_true(any(metabolite < 0))
})

test_that("diagnostic plotting helper saves expected PDFs", {
  testthat::skip_if_not_installed("ggplot2")

  diagnostics <- expand.grid(
    rep = seq_len(3),
    crosstalk = c("0", "1"),
    view_pair = c("mrna-mirna", "mrna-protein"),
    stringsAsFactors = FALSE
  )
  pair_split <- strsplit(diagnostics$view_pair, "-", fixed = TRUE)
  diagnostics$view_a <- vapply(pair_split, `[`, character(1), 1)
  diagnostics$view_b <- vapply(pair_split, `[`, character(1), 2)
  diagnostics$de_align <- 1
  diagnostics$mean_abs_feature_cor <- ifelse(diagnostics$crosstalk == "1", 0.12, 0.05)
  diagnostics$signal_agreement <- ifelse(diagnostics$crosstalk == "1", 0.7, 0.02)
  diagnostics$y_signal_cor <- 0.8
  diagnostics$event_rate <- NA_real_
  diagnostics$observed_event_rate <- NA_real_

  paths <- save_template_diagnostic_plots(
    diagnostics = diagnostics,
    plot_dir = tempdir(),
    prefix = "test_template"
  )

  expect_named(paths, c(
    "crosstalk_line",
    "feature_correlation_heatmap",
    "replicate_variability",
    "featurecor_signal_scatter",
    "all_plots"
  ))
  expect_true(all(file.exists(paths)))
})

test_that("diagnostic summaries retain random-effect rows", {
  diagnostics <- data.frame(
    rep = c(1, 2),
    view_a = c("mrna", "mrna"),
    view_b = c("protein", "protein"),
    crosstalk = c("0", "0"),
    de_align = c(NA_real_, NA_real_),
    mean_abs_feature_cor = c(0.04, 0.06),
    signal_agreement = c(-0.01, 0.03),
    y_signal_cor = c(0.8, 0.82),
    stringsAsFactors = FALSE
  )

  summary <- summarize_template_diagnostics(diagnostics)

  expect_equal(nrow(summary), 1)
  expect_equal(summary$effect_alignment, "random effects")
  expect_equal(summary$signal_agreement, 0.01)
})

test_that("blockwise feature-correlation diagnostics support matrices and simulations", {
  testthat::skip_if_not_installed("ggplot2")

  set.seed(4)
  n <- 30
  shared <- stats::rnorm(n)
  matrices <- list(
    mrna = cbind(g1 = shared + stats::rnorm(n), g2 = stats::rnorm(n), g3 = stats::rnorm(n)),
    protein = cbind(p1 = shared + stats::rnorm(n), p2 = stats::rnorm(n))
  )
  rownames(matrices$mrna) <- rownames(matrices$protein) <- paste0("S", seq_len(n))

  corr_df <- blockwise_feature_correlation(
    matrices,
    max_features_per_view = c(mrna = 2, protein = 2)
  )

  expect_true(all(c("view_x", "view_y", "order_x", "order_y", "correlation") %in% names(corr_df)))
  expect_equal(nrow(corr_df), 16)
  expect_s3_class(
    plot_blockwise_feature_correlation(
      corr_df,
      show_feature_diagonal = FALSE,
      color_limit = 0.5
    ),
    "ggplot"
  )
  block_summary <- summarize_blockwise_feature_correlation(corr_df)
  expect_equal(nrow(block_summary), 3)
  expect_true(any(block_summary$view_x == "mrna" & block_summary$view_y == "protein"))
  expect_true(any(block_summary$diagonal_block))
  expect_s3_class(plot_blockwise_correlation_summary(block_summary), "ggplot")

  tm <- fit_template(matrices, rank = 2, shared.rank = 1)
  sim <- gen_simmba(
    nsample = 40,
    template = tm,
    p = c(mrna = 3, protein = 2),
    crosstalk = 0.5,
    de.align = 1,
    nrep = 1,
    seed = 44
  )

  sim_matrices <- simulation_view_matrices(sim, rep = 1)
  expect_named(sim_matrices, c("mrna", "protein"))

  path <- save_blockwise_feature_correlation_plot(
    matrices = sim_matrices,
    plot_dir = tempdir(),
    prefix = "test_blockwise",
    max_features_per_view = 2,
    show_feature_diagonal = FALSE,
    color_limit = 0.5
  )
  expect_true(file.exists(path))
})

test_that("high-level validation figures support PCA, marginals, block differences, and crosstalk gradients", {
  testthat::skip_if_not_installed("ggplot2")

  set.seed(5)
  n_obs <- 35
  n_sim <- 40
  observed <- list(
    mrna = matrix(stats::rnorm(n_obs * 8), n_obs, 8),
    protein = matrix(stats::rnorm(n_obs * 6), n_obs, 6)
  )
  simulated <- list(
    mrna = matrix(stats::rnorm(n_sim * 8), n_sim, 8),
    protein = matrix(stats::rnorm(n_sim * 6), n_sim, 6)
  )
  colnames(observed$mrna) <- colnames(simulated$mrna) <- paste0("g", seq_len(8))
  colnames(observed$protein) <- colnames(simulated$protein) <- paste0("p", seq_len(6))

  pca_diag <- multiomics_pca_diagnostic(
    observed,
    simulated,
    max_features_per_view = c(mrna = 5, protein = 4)
  )
  expect_true(all(c("PC1", "PC2", "data_type") %in% names(pca_diag$scores)))
  expect_s3_class(plot_multiomics_pca_diagnostic(pca_diag), "ggplot")

  marginal_diag <- marginal_distribution_diagnostic(
    observed,
    simulated,
    max_features_per_view = c(mrna = 5, protein = 4)
  )
  expect_true(all(c("view", "feature", "ks_stat") %in% names(marginal_diag)))
  expect_s3_class(plot_marginal_distribution_diagnostic(marginal_diag), "ggplot")

  obs_corr <- blockwise_feature_correlation(observed, max_features_per_view = 4)
  sim_corr <- blockwise_feature_correlation(simulated, max_features_per_view = 4)
  obs_summary <- summarize_blockwise_feature_correlation(obs_corr)
  sim_summary <- summarize_blockwise_feature_correlation(sim_corr)
  comparison <- compare_blockwise_correlation_summary(obs_summary, sim_summary)
  expect_true("difference" %in% names(comparison))
  expect_s3_class(plot_blockwise_correlation_difference(comparison), "ggplot")

  diagnostics <- expand.grid(
    rep = seq_len(2),
    crosstalk = c("0", "0.5", "1"),
    view_a = "mrna",
    view_b = "protein",
    stringsAsFactors = FALSE
  )
  diagnostics$de_align <- 1
  diagnostics$mean_abs_feature_cor <- as.numeric(diagnostics$crosstalk) * 0.1 + 0.05
  diagnostics$signal_agreement <- as.numeric(diagnostics$crosstalk) * 0.6

  expect_s3_class(plot_crosstalk_gradient_diagnostic(diagnostics), "ggplot")
})

test_that("full validation figure set saves individual and combined PDFs", {
  testthat::skip_if_not_installed("ggplot2")

  set.seed(6)
  n_obs <- 30
  n_sim <- 30
  observed <- list(
    mrna = matrix(stats::rnorm(n_obs * 6), n_obs, 6),
    protein = matrix(stats::rnorm(n_obs * 5), n_obs, 5)
  )
  simulated <- list(
    mrna = matrix(stats::rnorm(n_sim * 6), n_sim, 6),
    protein = matrix(stats::rnorm(n_sim * 5), n_sim, 5)
  )
  colnames(observed$mrna) <- colnames(simulated$mrna) <- paste0("g", seq_len(6))
  colnames(observed$protein) <- colnames(simulated$protein) <- paste0("p", seq_len(5))

  diagnostics <- expand.grid(
    rep = seq_len(2),
    crosstalk = c("0", "1"),
    view_a = "mrna",
    view_b = "protein",
    stringsAsFactors = FALSE
  )
  diagnostics$de_align <- 1
  diagnostics$mean_abs_feature_cor <- ifelse(diagnostics$crosstalk == "1", 0.12, 0.05)
  diagnostics$signal_agreement <- ifelse(diagnostics$crosstalk == "1", 0.5, 0.02)
  diagnostics$y_signal_cor <- 0.8

  paths <- save_full_validation_figure_set(
    diagnostics = diagnostics,
    observed_matrices = observed,
    simulated_matrices = simulated,
    plot_dir = tempdir(),
    prefix = "test_full_validation",
    max_features_per_view = 4,
    width = 6,
    height = 5
  )

  expect_true("full_validation_figure_set" %in% names(paths))
  expect_true(all(file.exists(paths)))
})
