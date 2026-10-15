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
  expect_equal(dim(tm$observed$mrna), dim(template$mrna))
  expect_equal(dim(tm$observed$protein), dim(template$protein))
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

test_that("MAE-first simulation stores data and truth in Bioconductor containers", {
  skip_if_not_installed("MultiAssayExperiment")
  skip_if_not_installed("SummarizedExperiment")

  set.seed(101)
  n <- 50
  sample_ids <- paste0("S", seq_len(n))
  mrna <- matrix(stats::rnorm(n * 12), n, 12,
                 dimnames = list(sample_ids, paste0("m", seq_len(12))))
  protein <- matrix(stats::rnorm(n * 10), n, 10,
                    dimnames = list(sample_ids, paste0("p", seq_len(10))))
  input_mae <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = MultiAssayExperiment::ExperimentList(
      mrna = SummarizedExperiment::SummarizedExperiment(assays = list(counts = t(mrna))),
      protein = SummarizedExperiment::SummarizedExperiment(assays = list(counts = t(protein)))
    )
  )

  sim_mae <- simulate_mimesys(
    template = input_mae,
    n = 60,
    p = c(mrna = 8, protein = 8),
    crosstalk = 0.5,
    outcome = list(model = "shared_latent", shared_latent_function = "linear"),
    template_args = list(rank = 3, shared.rank = 1),
    nrep = 1,
    seed = 102
  )

  expect_s4_class(sim_mae, "MultiAssayExperiment")
  expect_setequal(names(MultiAssayExperiment::experiments(sim_mae)), c("mrna", "protein"))
  expect_true(all(c("Y", "split") %in% colnames(MultiAssayExperiment::colData(sim_mae))))
  expect_true("beta" %in% colnames(SummarizedExperiment::rowData(
    MultiAssayExperiment::experiments(sim_mae)[["mrna"]]
  )))
  expect_true("simulation" %in% names(S4Vectors::metadata(sim_mae)))
})

test_that("MAE-first simulation returns a replicate collection when requested", {
  skip_if_not_installed("MultiAssayExperiment")
  set.seed(103)
  template <- list(
    mrna = matrix(stats::rnorm(50 * 10), 50, 10),
    protein = matrix(stats::rnorm(50 * 10), 50, 10)
  )
  rownames(template$mrna) <- rownames(template$protein) <- paste0("S", seq_len(50))

  sim <- simulate_mimesys(
    template, n = 60, p = 8, crosstalk = 0.5, nrep = 2,
    template_args = list(rank = 2, shared.rank = 1), seed = 104
  )
  expect_s3_class(sim, "mimesys_mae_simulation")
  expect_length(sim$replicates, 2)
  expect_s4_class(sim$replicates[[1]], "MultiAssayExperiment")
})

test_that("MAE-first defaults calibrate crosstalk and null output contains no outcome", {
  skip_if_not_installed("MultiAssayExperiment")
  set.seed(105)
  n <- 45
  latent <- stats::rnorm(n)
  template <- list(
    mrna = outer(latent, stats::rnorm(10)) + matrix(stats::rnorm(n * 10), n, 10),
    protein = outer(latent, stats::rnorm(10)) + matrix(stats::rnorm(n * 10), n, 10)
  )
  rownames(template$mrna) <- rownames(template$protein) <- paste0("S", seq_len(n))
  calibration <- list(grid = c(0, 0.5, 1), nrep = 1, max_features_per_view = 6)

  null_mae <- simulate_mimesys_null(
    template, n = 50, p = 8,
    template_args = list(rank = 2, shared.rank = 1),
    calibration = calibration, seed = 106
  )
  settings <- S4Vectors::metadata(null_mae)$simulation
  expect_s4_class(null_mae, "MultiAssayExperiment")
  expect_false("Y" %in% colnames(MultiAssayExperiment::colData(null_mae)))
  expect_true(is.numeric(settings$crosstalk))
  expect_true(is.list(settings$crosstalk_calibration))

  outcome_mae <- simulate_mimesys(
    template, n = 50, p = 8,
    template_args = list(rank = 2, shared.rank = 1),
    calibration = calibration, seed = 107
  )
  expect_true("Y" %in% colnames(MultiAssayExperiment::colData(outcome_mae)))
  expect_true(is.list(S4Vectors::metadata(outcome_mae)$simulation$crosstalk_calibration))
})

test_that("named presets are discoverable and can be overridden", {
  presets <- mimesys_presets()
  expect_setequal(
    presets$preset,
    c("default", "linear", "friedman", "friedman2", "shared_latent",
      "shared_latent_friedman", "shared_latent_friedman2")
  )

  set.seed(108)
  template <- list(
    mrna = matrix(stats::rnorm(50 * 10), 50, 10),
    protein = matrix(stats::rnorm(50 * 10), 50, 10)
  )
  rownames(template$mrna) <- rownames(template$protein) <- paste0("S", seq_len(50))
  sim_mae <- simulate_mimesys(
    template, n = 60, p = 8, crosstalk = 0.5,
    preset = "friedman",
    outcome = list(model = "linear"),
    template_args = list(rank = 2, shared.rank = 1), seed = 109
  )
  expect_equal(S4Vectors::metadata(sim_mae)$simulation$ygen.mode, "LM")
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

test_that("cross-view reconstruction error is the matched correlation-block RMSE", {
  observed <- list(
    view_a = cbind(a1 = 1:8, a2 = c(1, 3, 2, 5, 4, 7, 6, 8)),
    view_b = cbind(b1 = c(2, 1, 4, 3, 6, 5, 8, 7), b2 = c(8, 6, 7, 5, 4, 2, 3, 1))
  )
  simulated <- observed
  simulated$view_b[, "b1"] <- rev(simulated$view_b[, "b1"])

  result <- feature_correlation_frobenius(observed, simulated)
  cross <- result[!result$diagonal_block, , drop = FALSE]
  expected <- sqrt(mean((
    stats::cor(observed$view_a, observed$view_b) -
      stats::cor(simulated$view_a, simulated$view_b)
  )^2))

  expect_equal(nrow(cross), 1L)
  expect_equal(cross$rescaled_frobenius, expected)
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

test_that("shared latent outcome can use selected factors and nonlinear functions", {
  set.seed(15)
  n <- 110
  shared <- matrix(stats::rnorm(n * 5), n, 5)
  template <- list(
    mrna = shared %*% matrix(stats::rnorm(5 * 24), 5, 24) + matrix(stats::rnorm(n * 24), n, 24),
    protein = shared %*% matrix(stats::rnorm(5 * 24), 5, 24) + matrix(stats::rnorm(n * 24), n, 24)
  )
  rownames(template$mrna) <- rownames(template$protein) <- paste0("S", seq_len(n))

  tm <- fit_template(template, rank = 6, shared.rank = 5)
  sim_linear_3 <- gen_simmba(
    nsample = 80,
    template = tm,
    p = c(mrna = 14, protein = 14),
    ygen.mode = "SharedLatent",
    shared.latent.fun = "linear",
    shared.latent.n_factors = 3,
    crosstalk = 1,
    nrep = 1,
    seed = 151
  )
  sim_friedman <- gen_simmba(
    nsample = 80,
    template = tm,
    p = c(mrna = 14, protein = 14),
    ygen.mode = "SharedLatent",
    shared.latent.fun = "friedman",
    shared.latent.n_factors = 5,
    crosstalk = 1,
    nrep = 1,
    seed = 151
  )
  sim_friedman2 <- gen_simmba(
    nsample = 80,
    template = tm,
    p = c(mrna = 14, protein = 14),
    ygen.mode = "SharedLatent",
    shared.latent.fun = "friedman2",
    shared.latent.n_factors = 5,
    crosstalk = 1,
    nrep = 1,
    seed = 151
  )

  shared_train <- sim_linear_3$trainDat[[1]]$truth$shared
  shared_test <- sim_linear_3$testDat[[1]]$truth$shared
  shared_all <- rbind(shared_train, shared_test)
  id_all <- c(
    sim_linear_3$trainDat[[1]]$sample_metadata$subjectID,
    sim_linear_3$testDat[[1]]$sample_metadata$subjectID
  )
  expected_all <- as.numeric(scale(rowSums(scale(shared_all[, 1:3, drop = FALSE])) / sqrt(3)))
  names(expected_all) <- id_all
  expected <- unname(expected_all[sim_linear_3$trainDat[[1]]$sample_metadata$subjectID])

  expect_equal(sim_linear_3$shared.latent.fun, "linear")
  expect_equal(sim_linear_3$shared.latent.n_factors, 3)
  expect_equal(sim_linear_3$trainDat[[1]]$sample_metadata$Xbeta_shared, expected, tolerance = 1e-8)
  expect_true(all(is.finite(sim_friedman$trainDat[[1]]$sample_metadata$Xbeta_shared)))
  expect_true(all(is.finite(sim_friedman2$trainDat[[1]]$sample_metadata$Xbeta_shared)))
  expect_false(isTRUE(all.equal(
    sim_friedman$trainDat[[1]]$sample_metadata$Xbeta_shared,
    sim_friedman2$trainDat[[1]]$sample_metadata$Xbeta_shared,
    tolerance = 1e-8
  )))
})

test_that("nonlinear shared latent outcomes require five shared factors", {
  set.seed(16)
  n <- 70
  shared <- matrix(stats::rnorm(n * 2), n, 2)
  template <- list(
    mrna = shared %*% matrix(stats::rnorm(2 * 16), 2, 16) + matrix(stats::rnorm(n * 16), n, 16),
    protein = shared %*% matrix(stats::rnorm(2 * 16), 2, 16) + matrix(stats::rnorm(n * 16), n, 16)
  )
  rownames(template$mrna) <- rownames(template$protein) <- paste0("S", seq_len(n))

  tm <- fit_template(template, rank = 4, shared.rank = 2)
  expect_error(
    gen_simmba(
      nsample = 60,
      template = tm,
      p = c(mrna = 12, protein = 12),
      ygen.mode = "SharedLatent",
      shared.latent.fun = "friedman",
      crosstalk = 1,
      nrep = 1,
      seed = 161
    ),
    "requires at least five"
  )
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

test_that("2logp rank, input validation, and template print are available", {
  set.seed(71)
  n <- 50
  template <- list(
    mrna = matrix(stats::rnorm(n * 12), n, 12),
    protein = matrix(stats::rnorm(n * 10), n, 10)
  )
  rownames(template$mrna) <- rownames(template$protein) <- paste0("S", seq_len(n))

  input_report <- validate_mimesys_input(template)
  expect_equal(attr(input_report, "common_samples"), n)
  expect_true(all(c("missing_values", "constant_features") %in% names(input_report)))

  tm <- fit_template(template, rank = "2logp", shared.rank = 1)
  expect_equal(tm$fits$mrna$rank, ceiling(2 * log(12)))
  expect_equal(tm$provenance$rank_input, "2logp")
  expect_output(print(tm), "MimeSys template")
})

test_that("calibration and validation reports retain uncertainty diagnostics", {
  skip_if_not_installed("ggplot2")
  set.seed(72)
  n <- 45
  shared <- stats::rnorm(n)
  template <- list(
    mrna = outer(shared, stats::rnorm(9)) + matrix(stats::rnorm(n * 9), n, 9),
    protein = outer(shared, stats::rnorm(8)) + matrix(stats::rnorm(n * 8), n, 8)
  )
  rownames(template$mrna) <- rownames(template$protein) <- paste0("S", seq_len(n))
  tm <- fit_template(template, rank = 2, shared.rank = 1)

  calibration <- calibrate_crosstalk(
    tm, crosstalk_grid = c(0, 0.5, 1), nsample = 40, nrep = 2,
    p = 6, max_features_per_view = 5, seed = 73
  )
  expect_equal(sum(calibration$selection_frequency$selection_frequency), 1)
  expect_gte(nrow(calibration$near_optimal), 1)
  expect_true(is.list(calibration$provenance))

  report <- validate_template(
    tm, p = 6, nsample = 40, nrep = 1,
    crosstalk_grid = c(0, 0.5), max_features_per_view = 5, seed = 74
  )
  expect_s3_class(report, "mimesys_template_validation")
  expect_true(all(c("pearson", "rmse") %in% names(report$correlation_structure)))
  expect_s3_class(plot_template_validation(report, "correlation_structure"), "ggplot")
  expect_s3_class(plot_template_validation(report, "marginal"), "ggplot")
  expect_output(print(report), "MimeSys template validation")
})

test_that("MAE simulations retain generation provenance", {
  skip_if_not_installed("MultiAssayExperiment")
  set.seed(75)
  n <- 45
  template <- list(
    mrna = matrix(stats::rnorm(n * 8), n, 8),
    protein = matrix(stats::rnorm(n * 8), n, 8)
  )
  rownames(template$mrna) <- rownames(template$protein) <- paste0("S", seq_len(n))

  sim_mae <- simulate_mimesys(
    template, n = 40, p = 6, crosstalk = 0.5,
    template_args = list(rank = 2, shared.rank = 1), seed = 76
  )
  provenance <- S4Vectors::metadata(sim_mae)$simulation$provenance
  expect_equal(provenance$selected_crosstalk, 0.5)
  expect_true(all(c("package_version", "created_at", "template") %in% names(provenance)))
})
