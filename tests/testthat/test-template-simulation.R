test_that("template simulator preserves requested structure", {
    set.seed(1)
    n <- 60
    latent <- matrix(stats::rnorm(n * 2), n, 2)
    template <- list(
        mrna = latent %*% matrix(stats::rnorm(2 * 30), 2,
            30) + matrix(stats::rnorm(n * 30), n, 30),
        protein = latent %*% matrix(stats::rnorm(2 * 20), 2,
            20) + matrix(stats::rnorm(n * 20), n, 20)
    )
    rownames(template$mrna) <- rownames(template$protein) <- paste0("S",
        seq_len(n))

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

    expect_null(sim$template)
    expect_true(is.list(sim$template_provenance))
    kept <- gen_simmba(nsample = 30, template = tm, nrep = 1, seed = 1,
        return.template = TRUE)
    expect_s3_class(kept$template, "mimesys_template")
    expect_length(sim$trainDat, 2)
    expect_equal(nrow(sim$trainDat[[1]]$feature_table), 40)
    expect_true("beta" %in% names(sim$trainDat[[1]]$feature_metadata))
    expect_true("Y" %in% names(sim$trainDat[[1]]$sample_metadata))
})

test_that("MAE-first simulation stores data and truth in Bioconductor containers",
    {
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
        outcome = list(model = "shared_latent",
            shared_latent_function = "linear"),
        template_args = list(rank = 3, shared.rank = 1),
        nrep = 1,
        seed = 102
    )

    expect_s4_class(sim_mae, "MultiAssayExperiment")
    expect_setequal(names(MultiAssayExperiment::experiments(sim_mae)), c("mrna",
        "protein"))
    expect_true(all(c("Y",
        "split") %in% colnames(MultiAssayExperiment::colData(sim_mae))))
    expect_true("beta" %in% colnames(SummarizedExperiment::rowData(
        MultiAssayExperiment::experiments(sim_mae)[["mrna"]]
    )))
    expect_true("simulation" %in% names(S4Vectors::metadata(sim_mae)))
})

test_that("MAE-first simulation returns a replicate collection when requested",
    {
    skip_if_not_installed("MultiAssayExperiment")
    set.seed(103)
    template <- list(
        mrna = matrix(stats::rnorm(50 * 10), 50, 10),
        protein = matrix(stats::rnorm(50 * 10), 50, 10)
    )
    rownames(template$mrna) <- rownames(template$protein) <- paste0("S",
        seq_len(50))

    sim <- simulate_mimesys(
        template, n = 60, p = 8, crosstalk = 0.5, nrep = 2,
        template_args = list(rank = 2, shared.rank = 1), seed = 104
    )
    expect_s3_class(sim, "mimesys_mae_simulation")
    expect_length(sim$replicates, 2)
    expect_s4_class(sim$replicates[[1]], "MultiAssayExperiment")
})

test_that("MAE-first defaults calibrate crosstalk and null output contains no outcome",
    {
    skip_if_not_installed("MultiAssayExperiment")
    set.seed(105)
    n <- 45
    latent <- stats::rnorm(n)
    template <- list(
        mrna = outer(latent, stats::rnorm(10)) + matrix(stats::rnorm(n * 10), n,
            10),
        protein = outer(latent, stats::rnorm(10)) + matrix(stats::rnorm(n * 10),
            n, 10)
    )
    rownames(template$mrna) <- rownames(template$protein) <- paste0("S",
        seq_len(n))
    calibration <- list(grid = c(0, 0.5, 1), nrep = 1,
        max_features_per_view = 6)

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
            "shared_latent_friedman", "shared_latent_friedman2", "interaction")
    )

    set.seed(108)
    template <- list(
        mrna = matrix(stats::rnorm(50 * 10), 50, 10),
        protein = matrix(stats::rnorm(50 * 10), 50, 10)
    )
    rownames(template$mrna) <- rownames(template$protein) <- paste0("S",
        seq_len(50))
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
    rownames(template$mrna) <- rownames(template$protein) <- paste0("S",
        seq_len(n))
    tm <- fit_template(template, rank = 2, shared.rank = 1)

    scalar <- sim_from_template(tm, n = 30, p = 6)
    positional <- sim_from_template(tm, n = 30, p = c(5, 7))
    named <- sim_from_template(tm, n = 30, p = c(protein = 4, mrna = 8))

    scalar_counts <- table(scalar$feature_metadata$featureType)
    positional_counts <- table(positional$feature_metadata$featureType)
    named_counts <- table(named$feature_metadata$featureType)

    expect_equal(stats::setNames(as.integer(scalar_counts),
        names(scalar_counts)), c(mrna = 6L, protein = 6L))
    expect_equal(stats::setNames(as.integer(positional_counts),
        names(positional_counts)), c(mrna = 5L, protein = 7L))
    expect_equal(stats::setNames(as.integer(named_counts), names(named_counts)),
        c(mrna = 8L, protein = 4L))
    expect_error(sim_from_template(tm, n = 30, p = c(mrna = 5)), "missing view")
})

test_that("cross-view reconstruction error is the matched correlation-block RMSE",
    {
    observed <- list(
        view_a = cbind(a1 = 1:8, a2 = c(1, 3, 2, 5, 4, 7, 6, 8)),
        view_b = cbind(b1 = c(2, 1, 4, 3, 6, 5, 8, 7), b2 = c(8, 6, 7, 5, 4, 2,
            3, 1))
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
    rownames(template$mrna) <- rownames(template$protein) <- paste0("S",
        seq_len(n))

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

test_that("Friedman inputs are rank transformed before nonlinear outcome generation",
    {
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
        mrna = outer(shared, stats::rnorm(16)) + matrix(stats::rnorm(n * 16), n,
            16),
        protein = outer(shared, stats::rnorm(16)) + matrix(stats::rnorm(n * 16),
            n, 16)
    )
    rownames(template$mrna) <- rownames(template$protein) <- paste0("S",
        seq_len(n))

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
    expect_equal(nrow(late_sim$trainDat[[1]]$truth$shared),
        nrow(late_sim$trainDat[[1]]$sample_metadata))
})

test_that("shared latent outcome uses all available shared factors", {
    set.seed(14)
    n <- 90
    shared <- matrix(stats::rnorm(n * 2), n, 2)
    template <- list(
        mrna = shared %*% matrix(stats::rnorm(2 * 18), 2,
            18) + matrix(stats::rnorm(n * 18), n, 18),
        protein = shared %*% matrix(stats::rnorm(2 * 18), 2,
            18) + matrix(stats::rnorm(n * 18), n, 18)
    )
    rownames(template$mrna) <- rownames(template$protein) <- paste0("S",
        seq_len(n))

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
    expect_equal(sim$trainDat[[1]]$sample_metadata$Xbeta_shared, expected,
        tolerance = 1e-8)
    expect_false(isTRUE(all.equal(expected, as.numeric(scale(shared_train[,
        1])), tolerance = 1e-8)))
})

test_that("shared latent outcome can use selected factors and nonlinear functions",
    {
    set.seed(15)
    n <- 110
    shared <- matrix(stats::rnorm(n * 5), n, 5)
    template <- list(
        mrna = shared %*% matrix(stats::rnorm(5 * 24), 5,
            24) + matrix(stats::rnorm(n * 24), n, 24),
        protein = shared %*% matrix(stats::rnorm(5 * 24), 5,
            24) + matrix(stats::rnorm(n * 24), n, 24)
    )
    rownames(template$mrna) <- rownames(template$protein) <- paste0("S",
        seq_len(n))

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
    expected_all <- as.numeric(scale(rowSums(scale(shared_all[, 1:3,
        drop = FALSE])) / sqrt(3)))
    names(expected_all) <- id_all
    expected <- unname(expected_all[sim_linear_3$trainDat[[1]]$sample_metadata$subjectID])

    expect_equal(sim_linear_3$shared.latent.fun, "linear")
    expect_equal(sim_linear_3$shared.latent.n_factors, 3)
    expect_equal(sim_linear_3$trainDat[[1]]$sample_metadata$Xbeta_shared,
        expected, tolerance = 1e-8)
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
        mrna = shared %*% matrix(stats::rnorm(2 * 16), 2,
            16) + matrix(stats::rnorm(n * 16), n, 16),
        protein = shared %*% matrix(stats::rnorm(2 * 16), 2,
            16) + matrix(stats::rnorm(n * 16), n, 16)
    )
    rownames(template$mrna) <- rownames(template$protein) <- paste0("S",
        seq_len(n))

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
    rownames(template$mrna) <- rownames(template$protein) <- paste0("S",
        seq_len(n))

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

test_that("crosstalk increases observed cross-view dependence in template simulations",
    {
    set.seed(3)
    n <- 80
    shared <- matrix(stats::rnorm(n * 2), n, 2)
    template <- list(
        mrna = shared %*% matrix(stats::rnorm(2 * 35), 2,
            35) + matrix(stats::rnorm(n * 35), n, 35),
        protein = shared %*% matrix(stats::rnorm(2 * 30), 2,
            30) + matrix(stats::rnorm(n * 30), n, 30)
    )
    rownames(template$mrna) <- rownames(template$protein) <- paste0("S",
        seq_len(n))
    tm <- fit_template(template, rank = 4, shared.rank = 1)

    low <- sim_from_template(tm, n = 120, p = c(mrna = 25, protein = 25),
        crosstalk = 0)
    high <- sim_from_template(tm, n = 120, p = c(mrna = 25, protein = 25),
        crosstalk = 1)

    x_low <- t(as.matrix(low$feature_table))
    x_high <- t(as.matrix(high$feature_table))
    low_type <- low$feature_metadata$featureType
    high_type <- high$feature_metadata$featureType

    low_cor <- mean(abs(stats::cor(x_low[, low_type == "mrna"], x_low[,
        low_type == "protein"])))
    high_cor <- mean(abs(stats::cor(x_high[, high_type == "mrna"], x_high[,
        high_type == "protein"])))

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

test_that("two-part view model preserves sparse compositional microbiome output",
    {
    set.seed(31)
    n <- 120
    shared <- stats::rnorm(n)
    microbe_latent <- outer(shared,
        stats::rnorm(18)) + matrix(stats::rnorm(n * 18), n, 18)
    microbiome <- exp(microbe_latent)
    microbiome[microbe_latent < 0.25] <- 0
    microbiome[rowSums(microbiome) == 0, 1] <- 1
    microbiome <- sweep(microbiome, 1, rowSums(microbiome), "/")
    metabolite <- outer(shared, stats::rnorm(12)) + matrix(stats::rnorm(n * 12),
        n, 12)
    template <- list(microbiome = microbiome, metabolite = metabolite)
    rownames(template$microbiome) <- rownames(template$metabolite) <- paste0("S",
        seq_len(n))

    tm <- fit_template(template, rank = 3, shared.rank = 1)
    expect_equal(tm$fits$microbiome$composition, list(type = "closed", total = 1))

    ## Full feature set: samples are re-closed to the template total
    sim <- sim_from_template(
        tm, n = 400, crosstalk = 0.5,
        view.models = c(microbiome = "two_part_copula")
    )
    x <- t(as.matrix(sim$feature_table))
    microbe <- x[, sim$feature_metadata$featureType == "microbiome",
        drop = FALSE]
    metabolite <- x[, sim$feature_metadata$featureType == "metabolite",
        drop = FALSE]
    expect_true(all(sim$feature_metadata$view_model[sim$feature_metadata$featureType == "microbiome"] == "two_part_copula"))
    expect_true(mean(microbe == 0) > 0)
    expect_equal(as.numeric(rowSums(microbe)), rep(1, nrow(microbe)),
        tolerance = 1e-8)
    expect_true(any(metabolite < 0))

    ## Feature subset: no re-closing, and positive abundances are drawn from the
    ## whole positive distribution (not only its upper quantiles)
    sub <- sim_from_template(
        tm, n = 2000, p = c(microbiome = 12, metabolite = 10), crosstalk = 0.5,
        view.models = c(microbiome = "two_part_copula")
    )
    fm <- sub$feature_metadata
    xs <- t(as.matrix(sub$feature_table))[, fm$featureType == "microbiome",
        drop = FALSE]
    expect_true(any(abs(rowSums(xs) - 1) > 1e-6))
    src <- fm$template_feature[fm$featureType == "microbiome"]
    shift <- vapply(seq_along(src), function(j) {
        sim_pos <- xs[xs[, j] > 0, j]
        obs_pos <- tm$observed$microbiome[, src[j]]
        obs_pos <- obs_pos[obs_pos > 0]
        log10(stats::median(sim_pos)) - log10(stats::median(obs_pos))
        }, numeric(1))
    expect_lt(abs(stats::median(shift)), 0.1)
    zero_err <- abs(colMeans(xs == 0) - colMeans(tm$observed$microbiome[,
        src] == 0))
    expect_lt(mean(zero_err), 0.05)
})

test_that("binary outcomes hit the target prevalence regardless of feature scale",
    {
    set.seed(41)
    n <- 80
    latent <- stats::rnorm(n)
    template <- list(
        small = outer(latent, stats::rnorm(12)) + matrix(stats::rnorm(n * 12),
            n, 12),
        huge = 1e7 * exp(outer(latent,
            stats::rnorm(12)) + matrix(stats::rnorm(n * 12), n, 12))
    )
    rownames(template$small) <- rownames(template$huge) <- paste0("S",
        seq_len(n))
    tm <- fit_template(template, rank = 2, shared.rank = 1)

    for (prev in c(0.2, 0.5)) {
        sim <- gen_simmba(
            nsample = 400, template = tm, outcome.type = "binary",
            binary.prevalence = prev, de.prob = 0.5, nrep = 4, seed = 2
        )
        rates <- vapply(seq_len(4), function(k) {
            y <- c(sim$trainDat[[k]]$sample_metadata$Y,
                sim$testDat[[k]]$sample_metadata$Y)
            mean(y)
            }, numeric(1))
        expect_true(all(abs(rates - prev) < 0.08))
        eta <- c(sim$trainDat[[1]]$sample_metadata$Xbeta,
            sim$testDat[[1]]$sample_metadata$Xbeta)
        expect_true(all(is.finite(eta)))
        expect_equal(stats::var(eta), pi^2 / 3, tolerance = 1e-6)
    }

    ## Both views contribute on standardized features: neither predictor is
    ## orders of magnitude larger than the other
    sim <- gen_simmba(nsample = 300, template = tm, de.prob = 0.5, nrep = 1,
        seed = 3)
    d <- rbind(t(as.matrix(sim$trainDat[[1]]$feature_table)))
    fm <- sim$trainDat[[1]]$feature_metadata
    v <- vapply(c("small", "huge"), function(view) {
        sel <- fm$featureType == view
        stats::var(as.numeric(MimeSys:::.standardize_columns(d[,
            sel]) %*% fm$beta[sel]))
        }, numeric(1))
    expect_lt(max(v) / min(v), 50)

    fr <- gen_simmba(
        nsample = 400, template = tm, outcome.type = "binary", ygen.mode = "Friedman",
        de.prob = 0.5, nrep = 1, seed = 4
    )
    y <- c(fr$trainDat[[1]]$sample_metadata$Y,
        fr$testDat[[1]]$sample_metadata$Y)
    expect_true(abs(mean(y) - 0.5) < 0.08)
})

test_that("survival outcomes are finite and censoring can be targeted", {
    set.seed(42)
    n <- 80
    latent <- stats::rnorm(n)
    template <- list(
        a = outer(latent, stats::rnorm(10)) + matrix(stats::rnorm(n * 10), n,
            10),
        b = 5e6 * exp(outer(latent,
            stats::rnorm(10)) + matrix(stats::rnorm(n * 10), n, 10))
    )
    rownames(template$a) <- rownames(template$b) <- paste0("S", seq_len(n))
    tm <- fit_template(template, rank = 2, shared.rank = 1)

    sim <- gen_simmba(nsample = 300, template = tm, outcome.type = "survival",
        de.prob = 0.5, nrep = 2, seed = 5)
    sm <- sim$trainDat[[1]]$sample_metadata
    expect_true(all(is.finite(sm$time)) && all(sm$time > 0))
    expect_true(all(sm$status %in% c(0L, 1L)))
    expect_true(mean(sm$status) > 0 && mean(sm$status) < 1)

    for (rate in c(0, 0.3, 0.6)) {
        sim <- gen_simmba(nsample = 500, template = tm,
            outcome.type = "survival",
            cens.rate = rate, de.prob = 0.5, nrep = 1, seed = 6)
        status <- c(sim$trainDat[[1]]$sample_metadata$status,
            sim$testDat[[1]]$sample_metadata$status)
        expect_true(abs(mean(status == 0) - rate) < 0.02)
    }
    expect_error(gen_simmba(nsample = 50, template = tm,
        outcome.type = "survival",
        cens.rate = 1, nrep = 1), "cens.rate")
    expect_error(gen_simmba(nsample = 50, template = tm,
        outcome.type = "binary",
        binary.prevalence = 0, nrep = 1), "binary.prevalence")
})

test_that("MAE interface passes prevalence and censoring options through", {
    testthat::skip_if_not_installed("MultiAssayExperiment")
    set.seed(43)
    n <- 60
    latent <- stats::rnorm(n)
    template <- list(
        a = outer(latent, stats::rnorm(8)) + matrix(stats::rnorm(n * 8), n, 8),
        b = outer(latent, stats::rnorm(8)) + matrix(stats::rnorm(n * 8), n, 8)
    )
    rownames(template$a) <- rownames(template$b) <- paste0("S", seq_len(n))
    tm <- fit_template(template, rank = 2, shared.rank = 1)
    mae <- simulate_mimesys(tm, n = 400, crosstalk = 1,
        outcome = list(type = "binary", prevalence = 0.25),
        signal = list(de_prob = 0.5))
    y <- SummarizedExperiment::colData(mae)$Y
    expect_true(abs(mean(y) - 0.25) < 0.08)
    sim <- simulate_mimesys(tm, n = 50, crosstalk = 1,
        return_type = "simulation")
    expect_null(sim$template)
})

test_that("fitted templates do not keep per-sample factor scores", {
    set.seed(44)
    x <- list(a = matrix(stats::rnorm(40 * 6), 40, 6),
        b = matrix(stats::rnorm(40 * 6), 40, 6))
    tm <- fit_template(x, rank = 2, shared.rank = 1)
    expect_null(tm$fits$a$scores)
    expect_null(tm$fits$a$Z)
    expect_s3_class(tm, "mimesys_template")
})

test_that("blockwise feature-correlation diagnostics support matrices and simulations",
    {
    testthat::skip_if_not_installed("ggplot2")

    set.seed(4)
    n <- 30
    shared <- stats::rnorm(n)
    matrices <- list(
        mrna = cbind(g1 = shared + stats::rnorm(n), g2 = stats::rnorm(n),
            g3 = stats::rnorm(n)),
        protein = cbind(p1 = shared + stats::rnorm(n), p2 = stats::rnorm(n))
    )
    rownames(matrices$mrna) <- rownames(matrices$protein) <- paste0("S",
        seq_len(n))

    corr_df <- blockwise_feature_correlation(
        matrices,
        max_features_per_view = c(mrna = 2, protein = 2)
    )

    expect_true(all(c("view_x", "view_y", "order_x", "order_y",
        "correlation") %in% names(corr_df)))
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
    rownames(template$mrna) <- rownames(template$protein) <- paste0("S",
        seq_len(n))

    input_report <- validate_mimesys_input(template)
    expect_equal(attr(input_report, "common_samples"), n)
    expect_true(all(c("missing_values",
        "constant_features") %in% names(input_report)))

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
        mrna = outer(shared, stats::rnorm(9)) + matrix(stats::rnorm(n * 9), n,
            9),
        protein = outer(shared, stats::rnorm(8)) + matrix(stats::rnorm(n * 8),
            n, 8)
    )
    rownames(template$mrna) <- rownames(template$protein) <- paste0("S",
        seq_len(n))
    tm <- fit_template(template, rank = 2, shared.rank = 1)

    calibration <- calibrate_crosstalk(
        tm, crosstalk_grid = c(0, 0.5, 1), nsample = 40, nrep = 2,
        p = 6, max_features_per_view = 5, seed = 73
    )
    expect_equal(sum(calibration$selection_frequency$selection_frequency), 1)
    expect_gte(nrow(calibration$near_optimal), 1)
    expect_true(is.list(calibration$provenance))

    report <- suppressWarnings(validate_template(
        tm, p = 6, nsample = 40, nrep = 1,
        crosstalk_grid = c(0, 0.5), max_features_per_view = 5, seed = 74
    ))
    expect_true(all(report$marginal_diagnostics$mean_ks >= 0 & report$marginal_diagnostics$mean_ks <= 1))
    expect_true(all(c("mean_ks", "median_ks",
        "mean_std_quantile_rmse") %in% names(report$marginal_diagnostics)))
    expect_s3_class(report, "mimesys_template_validation")
    expect_true(all(c("pearson",
        "rmse") %in% names(report$correlation_structure)))
    expect_s3_class(plot_template_validation(report, "correlation_structure"),
        "ggplot")
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
    rownames(template$mrna) <- rownames(template$protein) <- paste0("S",
        seq_len(n))

    sim_mae <- simulate_mimesys(
        template, n = 40, p = 6, crosstalk = 0.5,
        template_args = list(rank = 2, shared.rank = 1), seed = 76
    )
    provenance <- S4Vectors::metadata(sim_mae)$simulation$provenance
    expect_equal(provenance$selected_crosstalk, 0.5)
    expect_true(all(c("package_version", "created_at",
        "template") %in% names(provenance)))
})


.mk_template <- function(seed, n = 80, p = c(a = 12, b = 12), shared_rank = 1) {
    set.seed(seed)
    latent <- stats::rnorm(n)
    out <- lapply(p, function(k) outer(latent,
        stats::rnorm(k)) + matrix(stats::rnorm(n * k), n, k))
    for (v in names(out)) rownames(out[[v]]) <- paste0("S", seq_len(n))
    out
}

test_that("calibration uses the fitted scale and respects the valid crosstalk range",
    {
    x <- .mk_template(51)
    x$b <- exp(x$b)
    tm <- fit_template(x, rank = 2, shared.rank = 1,
        view.transform = c(b = "log1p"))
    ## Raw input plus a transform must be compared on the transformed scale
    seen <- NULL
    real_align <- MimeSys:::.mimesys_align_observed_simulated
    testthat::local_mocked_bindings(
        .mimesys_align_observed_simulated = function(observed_matrices,
            simulated_matrices, ...) {
            seen <<- observed_matrices
            real_align(observed_matrices, simulated_matrices, ...)
        }
    )
    calibrate_crosstalk(x, crosstalk_grid = c(0, 1), nrep = 1, seed = 3,
        rank = 2, shared.rank = 1, view.transform = c(b = "log1p"))
    expect_equal(unname(seen$b), unname(log1p(x$b)))

    lim <- crosstalk_max(tm)
    expect_warning(
        cal <- calibrate_crosstalk(tm, crosstalk_grid = c(0, 1, lim + 1),
            nrep = 1, seed = 1),
        "Dropping crosstalk_grid"
    )
    expect_false(any(cal$raw$crosstalk > lim))
    expect_error(
        suppressWarnings(calibrate_crosstalk(tm, crosstalk_grid = lim + 1,
            nrep = 1)),
        "No crosstalk_grid values"
    )
})

test_that("Friedman feature sampling never leaves the candidate set", {
    set.seed(52)
    ft <- rep(c("a", "b"), c(50, 50))
    idx <- c(3, 7, 9, 11, 77)
    picks <- replicate(200, MimeSys:::.sample_friedman_features(idx, ft, c("a",
        "b"), c(a = 1, b = 0)))
    expect_true(all(picks %in% idx))
    expect_equal(MimeSys:::.sample_values(42, 1), 42)
})

test_that("outcome_feature and beta describe the features that drive the outcome",
    {
    tm <- fit_template(.mk_template(53), rank = 2, shared.rank = 1)
    run <- function(mode, ...) {
        s <- gen_simmba(nsample = 120, template = tm, ygen.mode = mode,
            de.prob = 0.4,
            nrep = 1, seed = 2, ...)
        s$trainDat[[1]]$feature_metadata
    }
    lm_fm <- run("LM")
    expect_true(any(lm_fm$beta != 0))
    expect_identical(lm_fm$outcome_feature, lm_fm$beta != 0)

    fr <- run("Friedman")
    expect_true(all(fr$beta == 0))
    expect_identical(fr$outcome_feature, fr$friedman_term)
    expect_equal(sum(fr$friedman_term), 5)

    fr2 <- run("Friedman2")
    expect_identical(fr2$outcome_feature, fr2$beta != 0 | fr2$friedman_term)

    sl <- run("SharedLatent")
    expect_true(all(sl$beta == 0))
    expect_false(any(sl$outcome_feature))

    it <- run("Interaction", interaction.weight = 1)
    expect_true(all(it$beta == 0))
    expect_identical(it$outcome_feature, it$interaction_term)
})

test_that("nonlinear.weight mixes standardized linear and Friedman parts", {
    tm <- fit_template(.mk_template(54), rank = 2, shared.rank = 1)
    ## w = 0: Xbeta is the standardized linear part; w = 1: the standardized Friedman part
    full0 <- gen_simmba(nsample = 150, template = tm, ygen.mode = "Friedman2",
        de.prob = 0.4,
        nonlinear.weight = 0, nrep = 1, seed = 3)
    all0 <- rbind(full0$trainDat[[1]]$sample_metadata,
        full0$testDat[[1]]$sample_metadata)
    expect_equal(all0$Xbeta, as.numeric(scale(all0$Xbeta_linear)),
        tolerance = 1e-8)
    full1 <- gen_simmba(nsample = 150, template = tm, ygen.mode = "Friedman2",
        de.prob = 0.4,
        nonlinear.weight = 1, nrep = 1, seed = 3)
    all1 <- rbind(full1$trainDat[[1]]$sample_metadata,
        full1$testDat[[1]]$sample_metadata)
    expect_equal(all1$Xbeta, as.numeric(scale(all1$Xbeta_friedman)),
        tolerance = 1e-8)
    expect_equal(stats::var(all1$Xbeta), 1, tolerance = 1e-8)
    expect_error(gen_simmba(nsample = 50, template = tm, nonlinear.weight = 2,
        nrep = 1),
        "nonlinear.weight")
})

test_that("cross-view interaction outcomes use feature pairs from different views",
    {
    tm <- fit_template(.mk_template(55, p = c(a = 10, b = 10, c = 10)),
        rank = 2, shared.rank = 1)
    s <- gen_simmba(nsample = 200, template = tm, ygen.mode = "Interaction",
        interaction.weight = 1, de.prob = 0.3, nrep = 1, seed = 4)
    tr <- s$trainDat[[1]]
    pairs <- tr$truth$interaction_pairs
    expect_equal(nrow(pairs), 3)
    expect_true(all(pairs$view_x != pairs$view_y))
    expect_setequal(paste(pairs$view_x, pairs$view_y), c("a b", "a c", "b c"))
    fm <- tr$feature_metadata
    expect_true(all(fm$interaction_term[match(c(pairs$feature_x,
        pairs$feature_y), fm$featureID)]))
    all_sm <- rbind(tr$sample_metadata, s$testDat[[1]]$sample_metadata)
    expect_equal(all_sm$Xbeta, as.numeric(scale(all_sm$Xbeta_interaction)),
        tolerance = 1e-8)

    s2 <- gen_simmba(nsample = 200, template = tm, ygen.mode = "Interaction",
        n.interactions = 5, nrep = 1, seed = 4)
    expect_equal(nrow(s2$trainDat[[1]]$truth$interaction_pairs), 5)

    one <- fit_template(list(a = .mk_template(56)$a), rank = 2)
    expect_error(gen_simmba(nsample = 50, template = one,
        ygen.mode = "Interaction", nrep = 1),
        "at least two views")
})

test_that("input guards for shared-latent and train/test settings", {
    tm <- fit_template(.mk_template(57), rank = 2, shared.rank = 1)
    expect_warning(
        gen_simmba(nsample = 60, template = tm, ygen.mode = "SharedLatent",
            crosstalk = 0, nrep = 1),
        "unrelated to the features"
    )
    expect_error(gen_simmba(nsample = 10, template = tm, p.train = 0.99,
        nrep = 1),
        "empty training or test set")
    expect_error(gen_simmba(nsample = 10, template = tm, p.train = 1, nrep = 1),
        "p.train")
})

test_that("back.transform returns views on their original measurement scale", {
    set.seed(58)
    n <- 120
    latent <- stats::rnorm(n)
    counts <- matrix(stats::rpois(n * 10, lambda = exp(2 + outer(latent,
        stats::runif(10)))), n, 10)
    comp_latent <- outer(latent,
        stats::rnorm(15)) + matrix(stats::rnorm(n * 15), n, 15)
    comp <- exp(comp_latent)
    comp[comp_latent < 0] <- 0
    comp[rowSums(comp) == 0, 1] <- 1
    comp <- sweep(comp, 1, rowSums(comp), "/")
    x <- list(counts = counts, comp = comp)
    for (v in names(x)) rownames(x[[v]]) <- paste0("S", seq_len(n))
    tm <- fit_template(x, rank = 2, shared.rank = 1,
        view.transform = c(counts = "log1p", comp = "clr"))

    fitted <- sim_from_template(tm, n = 500)
    orig <- sim_from_template(tm, n = 500, back.transform = TRUE)
    fm <- orig$feature_metadata
    expect_true(all(fm$scale == "original"))
    xo <- t(as.matrix(orig$feature_table))
    xc <- xo[, fm$featureType == "counts"]
    xp <- xo[, fm$featureType == "comp"]
    expect_true(all(xc >= 0) && all(abs(xc - round(xc)) < 1e-8))
    expect_lt(abs(stats::median(xc) - stats::median(counts)),
        0.25 * stats::median(counts))
    expect_lt(abs(mean(xp == 0) - mean(comp == 0)), 0.05)
    expect_equal(as.numeric(rowSums(xp)), rep(1, nrow(xp)), tolerance = 1e-8)
    expect_true(any(t(as.matrix(fitted$feature_table))[,
        fm$featureType == "comp"] < 0))

    sim <- gen_simmba(nsample = 100, template = tm, nrep = 1, seed = 1,
        back.transform = TRUE)
    expect_true(all(sim$trainDat[[1]]$feature_metadata$scale == "original"))
    expect_true(all(sim$trainDat[[1]]$feature_table >= 0))
})

test_that("simulate_mimesys returns the original scale by default", {
    testthat::skip_if_not_installed("MultiAssayExperiment")
    set.seed(59)
    n <- 80
    latent <- stats::rnorm(n)
    x <- list(
        a = exp(outer(latent, stats::rnorm(8)) + matrix(stats::rnorm(n * 8), n,
            8)),
        b = outer(latent, stats::rnorm(8)) + matrix(stats::rnorm(n * 8), n, 8)
    )
    for (v in names(x)) rownames(x[[v]]) <- paste0("S", seq_len(n))
    tm <- fit_template(x, rank = 2, shared.rank = 1,
        view.transform = c(a = "log1p"))
    mae <- simulate_mimesys(tm, n = 100, crosstalk = 1,
        null_model = list(back_transform = TRUE))
    a <- SummarizedExperiment::assay(MultiAssayExperiment::experiments(mae)[["a"]])
    expect_true(all(a > 0))
    expect_gt(stats::median(a), 0.5 * stats::median(x$a))
    mae_int <- simulate_mimesys(tm, n = 100, crosstalk = 1,
        preset = "interaction", null_model = list(back_transform = TRUE))
    pairs <- S4Vectors::metadata(mae_int)$truth$interaction_pairs
    expect_equal(nrow(pairs), 1)
    expect_true(pairs$view_x != pairs$view_y)
    ## The default returns the fitted (log1p) scale
    mae_fit <- simulate_mimesys(tm, n = 100, crosstalk = 1)
    a_fit <- SummarizedExperiment::assay(MultiAssayExperiment::experiments(mae_fit)[["a"]])
    expect_lt(max(a_fit), max(a))
})

test_that("calibration picks the value with the smallest mean over cross-view blocks", {
    tm <- fit_template(.mk_template(61, p = c(a = 10, b = 10, c = 10)), rank = 2,
        shared.rank = 1, seed = 1)
    cal <- calibrate_crosstalk(tm, crosstalk_grid = c(0, 0.5, 1), nrep = 3, seed = 2)
    by_rep <- stats::aggregate(rescaled_frobenius ~ crosstalk + rep,
        data = cal$raw[!cal$raw$diagonal_block, ], FUN = mean)
    overall <- stats::aggregate(rescaled_frobenius ~ crosstalk, data = by_rep, FUN = mean)
    expect_equal(cal$best_crosstalk, overall$crosstalk[which.min(overall$rescaled_frobenius)])
    expect_equal(cal$cross_view_summary$mean, overall$rescaled_frobenius)
    expect_true(cal$best_crosstalk %in% cal$near_optimal$crosstalk)
})

test_that("seeds make fitting and MAE simulation reproducible without touching global RNG", {
    x <- .mk_template(62)
    set.seed(5)
    state <- .Random.seed
    a <- fit_template(x, seed = 3)
    b <- fit_template(x, seed = 3)
    expect_identical(state, .Random.seed)
    a$provenance$fitted_at <- b$provenance$fitted_at <- NULL
    expect_equal(a, b)
    testthat::skip_if_not_installed("MultiAssayExperiment")
    n1 <- simulate_mimesys_null(a, n = 20, p = 5, crosstalk = 0.5, seed = 1)
    n2 <- simulate_mimesys_null(a, n = 20, p = 5, crosstalk = 0.5, seed = 1)
    expect_identical(SummarizedExperiment::assay(MultiAssayExperiment::experiments(n1)[[1]]),
        SummarizedExperiment::assay(MultiAssayExperiment::experiments(n2)[[1]]))
    s1 <- simulate_mimesys(x, n = 20, p = 5, crosstalk = 0.5, seed = 1)
    s2 <- simulate_mimesys(x, n = 20, p = 5, crosstalk = 0.5, seed = 1)
    expect_identical(SummarizedExperiment::assay(MultiAssayExperiment::experiments(s1)[[1]]),
        SummarizedExperiment::assay(MultiAssayExperiment::experiments(s2)[[1]]))
    expect_identical(state, .Random.seed)
})

test_that("per-view arguments are matched by name and checked", {
    x <- .mk_template(63, p = c(a = 8, b = 8, c = 8))
    tm <- fit_template(x, rank = c(c = 1, a = 3, b = 2), shared.rank = 1, seed = 1)
    expect_equal(unname(vapply(tm$fits, `[[`, integer(1), "rank")), c(3L, 2L, 1L))
    s <- sim_from_template(tm, n = 10, crosstalk = c(c = 0, a = 1, b = 0.5))
    expect_equal(unname(s$truth$crosstalk), c(1, 0.5, 0))
    expect_error(sim_from_template(tm, n = 10, crosstalk = c(1, 0.5)), "length 1 or one value per view")
    expect_error(sim_from_template(tm, n = 10, crosstalk = c(a = 1, z = 1, b = 1)), "unknown view")
    expect_error(fit_template(x, view.transform = c(aa = "log1p")), "unknown view")
    g <- gen_simmba(nsample = 60, template = tm, de.prob = c(c = 0, a = 0.5, b = 0), nrep = 1, seed = 1)
    fm <- g$trainDat[[1]]$feature_metadata
    expect_true(all(fm$beta[fm$featureType %in% c("b", "c")] == 0))
    expect_error(gen_simmba(nsample = 50, template = tm, de.prob = 1.5, nrep = 1), "de.prob")
    expect_error(gen_simmba(nsample = 50, template = tm, de.align = 2, nrep = 1), "de.align")
    expect_error(gen_simmba(nsample = 50, template = tm, nrep = 0), "nrep")
})

test_that("input problems give clear messages", {
    x <- .mk_template(64)
    dup <- x
    rownames(dup$a)[2] <- rownames(dup$a)[1]
    expect_error(fit_template(dup), "Duplicate sample IDs")
    na <- x
    na$a[, 1] <- NA
    na$a[1:50, 2] <- NA
    expect_message(expect_warning(fit_template(na, rank = 2, shared.rank = 1),
        "more than 50% missing"), "all-missing")
    expect_error(fit_template(x, shared.rank = 100), "exceeds")
    const <- x
    const$b[] <- 1
    expect_error(fit_template(const), "no non-constant features")
})

test_that("zero-weight components are not marked as outcome features", {
    tm <- fit_template(.mk_template(65), rank = 2, shared.rank = 1, seed = 1)
    it <- gen_simmba(nsample = 80, template = tm, ygen.mode = "Interaction",
        interaction.weight = 0, de.prob = 0.4, nrep = 1, seed = 1)
    fm <- it$trainDat[[1]]$feature_metadata
    expect_false(any(fm$interaction_term))
    expect_null(it$trainDat[[1]]$truth$interaction_pairs)
    fr <- gen_simmba(nsample = 80, template = tm, ygen.mode = "Friedman2",
        nonlinear.weight = 0, de.prob = 0.4, nrep = 1, seed = 1)
    expect_false(any(fr$trainDat[[1]]$feature_metadata$friedman_term))
})

test_that("single-view templates skip calibration", {
    one <- fit_template(list(a = .mk_template(66)$a), rank = 2, seed = 1)
    expect_error(calibrate_crosstalk(one, nrep = 1), "at least two views")
    rep <- expect_message(validate_template(one, nrep = 1), "Skipping crosstalk calibration")
    expect_null(rep$calibration)
    testthat::skip_if_not_installed("MultiAssayExperiment")
    expect_message(mae <- simulate_mimesys(one, n = 30, seed = 1), "Skipping")
    expect_s4_class(mae, "MultiAssayExperiment")
})

test_that("summaries describe binary, survival and Friedman outcomes", {
    tm <- fit_template(.mk_template(67), rank = 2, shared.rank = 1, seed = 1)
    sv <- gen_simmba(nsample = 100, template = tm, outcome.type = "survival",
        cens.rate = 0.3, de.prob = 0.4, nrep = 2, seed = 1)
    ss <- summarize_simulation(sv)
    expect_equal(mean(ss$event_rate), 0.7, tolerance = 0.02)
    expect_true(is.finite(ss$y_signal_cor_mean))
    bn <- summarize_simulation(gen_simmba(nsample = 100, template = tm, outcome.type = "binary",
        binary.prevalence = 0.3, de.prob = 0.4, nrep = 2, seed = 1))
    expect_true(abs(mean(bn$prevalence) - 0.3) < 0.1)
    fr <- summarize_simulation(gen_simmba(nsample = 100, template = tm, ygen.mode = "Friedman",
        de.prob = 0.4, nrep = 1, seed = 1))
    expect_equal(fr$active_features_first_rep, 5)
})

test_that("views keep template order and large-mean features are standardized", {
    x <- .mk_template(68, p = c(zeta = 6, alpha = 6))
    tm <- fit_template(x, rank = 2, shared.rank = 1, seed = 1)
    sim <- gen_simmba(nsample = 40, template = tm, nrep = 1, seed = 1)
    expect_identical(names(simulation_view_matrices(sim)), c("zeta", "alpha"))
    testthat::skip_if_not_installed("MultiAssayExperiment")
    mae <- as_mimesys_mae(sim)
    expect_identical(names(MultiAssayExperiment::experiments(mae)), c("zeta", "alpha"))
    big <- cbind(1e9 + stats::rnorm(50), stats::rnorm(50))
    expect_equal(apply(MimeSys:::.standardize_columns(big), 2, stats::sd), c(1, 1))
})

test_that("friedman.share works at the view level", {
    set.seed(69)
    ft <- rep(c("a", "b", "c"), c(30, 5, 30))
    picks <- replicate(400, MimeSys:::.sample_friedman_features(seq_along(ft), ft,
        c("a", "b", "c"), c(a = 0.25, b = 0.5, c = 0.25)))
    expect_gt(mean(ft[picks] == "b"), 0.4)
})

test_that("shared-rank selection stores p-values and flags borderline choices", {
    tm <- fit_template(.mk_template(70), rank = 3, nperm = 50, seed = 1)
    expect_true(length(tm$gcca.pvalue) %in% c(tm$shared.rank, tm$shared.rank + 1L))
    expect_length(tm$gcca.step.eigen, length(tm$gcca.pvalue))
    expect_true(all(tm$gcca.pvalue[seq_len(tm$shared.rank)] < 0.05))
    expect_true(all(tm$gcca.pvalue > 0 & tm$gcca.pvalue <= 1))
    expect_warning(MimeSys:::.warn_borderline_shared_rank(1L, c(0.04, 0.5), 200), "borderline")
    expect_silent(MimeSys:::.warn_borderline_shared_rank(1L, c(0.005, 0.6), 200))
})

test_that("shared.rank cannot exceed the smallest within-view rank", {
    x <- .mk_template(71)
    expect_error(fit_template(x, rank = 3, shared.rank = 4), "smallest within-view rank")
    tm <- fit_template(x, rank = 3, shared.rank = 3, seed = 1)
    expect_equal(tm$shared.rank, 3L)
    expect_error(gen_simmba(nsample = 30, template = tm, p = c(a = 2, b = 2),
        ygen.mode = "Friedman", nrep = 1), "at least five")
})

test_that("simulated compositions never exceed the template total", {
    set.seed(72)
    n <- 120
    latent <- stats::rnorm(n)
    lat <- outer(latent, stats::rnorm(20)) + matrix(stats::rnorm(n * 20), n, 20)
    comp <- exp(lat)
    comp[lat < 0] <- 0
    comp[rowSums(comp) == 0, 1] <- 1
    comp <- sweep(comp, 1, rowSums(comp), "/")
    other <- outer(latent, stats::rnorm(8)) + matrix(stats::rnorm(n * 8), n, 8)
    x <- list(comp = comp, other = other)
    for (v in names(x)) rownames(x[[v]]) <- paste0("S", seq_len(n))
    totals <- function(sim) {
        z <- t(as.matrix(sim$feature_table))
        rowSums(z[, sim$feature_metadata$featureType == "comp", drop = FALSE])
    }
    for (vm in list(NULL, c(comp = "two_part_copula"))) {
        tm <- fit_template(x, rank = 2, shared.rank = 1, seed = 1)
        expect_identical(tm$fits$comp$composition$type, "closed")
        more <- sim_from_template(tm, n = 300, p = c(comp = 40, other = 8), view.models = vm)
        expect_equal(unname(totals(more)), rep(1, 300), tolerance = 1e-8)
        fewer <- sim_from_template(tm, n = 300, p = c(comp = 8, other = 8), view.models = vm)
        expect_true(all(totals(fewer) <= 1 + 1e-8))
    }
    ## Template that is itself a subset of a composition
    sub <- list(comp = comp[, 1:10], other = other)
    tm2 <- fit_template(sub, rank = 2, shared.rank = 1, seed = 1)
    expect_identical(tm2$fits$comp$composition$type, "proportion")
    lim <- max(rowSums(comp[, 1:10]))
    s2 <- sim_from_template(tm2, n = 500, p = c(comp = 20, other = 8))
    expect_true(all(totals(s2) <= lim + 1e-8))
    ## CLR-fitted view, back-transformed
    tm3 <- fit_template(x, rank = 2, shared.rank = 1, seed = 1, view.transform = c(comp = "clr"))
    s3 <- sim_from_template(tm3, n = 300, p = c(comp = 30, other = 8), back.transform = TRUE)
    expect_equal(unname(totals(s3)), rep(1, 300), tolerance = 1e-8)
})

test_that("every view with de.prob > 0 carries at least one signal feature", {
    tm <- fit_template(.mk_template(73, p = c(a = 3, b = 3)), rank = 1, shared.rank = 1, seed = 1)
    zero <- vapply(seq_len(30), function(s) {
        g <- gen_simmba(nsample = 40, template = tm, de.prob = 0.05, nrep = 1, seed = s)
        fm <- g$trainDat[[1]]$feature_metadata
        any(tapply(fm$beta != 0, fm$featureType, sum) == 0)
    }, logical(1))
    expect_false(any(zero))
    g <- gen_simmba(nsample = 40, template = tm, de.prob = c(a = 0.5, b = 0), nrep = 1, seed = 1)
    fm <- g$trainDat[[1]]$feature_metadata
    expect_true(all(fm$beta[fm$featureType == "b"] == 0))
})

test_that("stability selection repeats the shared-rank choice", {
    x <- .mk_template(74)
    tm <- fit_template(x, rank = 3, nperm = 30, shared.rank.runs = 5, seed = 1)
    sel <- tm$shared.rank.selection
    expect_length(sel$runs, 5)
    expect_equal(tm$shared.rank, as.integer(round(mean(sel$runs))))
    expect_equal(sum(sel$distribution), 5)
    one <- suppressWarnings(fit_template(x, rank = 3, nperm = 30, shared.rank.runs = 1, seed = 1))
    expect_null(one$shared.rank.selection)
    fixed <- fit_template(x, rank = 3, shared.rank = 1, seed = 1)
    expect_null(fixed$shared.rank.selection)
    expect_error(fit_template(x, shared.rank.runs = 0), "shared.rank.runs")
})

test_that("gen_mediation plants mediators and returns calibrated truth", {
    tm <- fit_template(.mk_template(80, p = c(a = 12, b = 12, c = 12)), rank = 2,
        shared.rank = 1, seed = 1)
    g <- gen_mediation(nsample = 300, template = tm, n.mediators = c(a = 2, b = 2, c = 0),
        TIE = 0.6, DE = 0.3, n.mc = 1000, nrep = 2, seed = 1)
    expect_length(g$replicates, 2)
    r <- g$replicates[[1]]
    tr <- r$truth$mediation
    expect_equal(tr$NIE, 0.6, tolerance = 1e-8)
    expect_equal(tr$NDE, 0.3)
    expect_equal(tr$TE, 0.9, tolerance = 1e-8)
    expect_equal(unname(tr$NIE_by_view[["c"]]), 0)
    expect_equal(sum(tr$NIE_by_view) + tr$NIE_interaction, tr$NIE, tolerance = 1e-8)
    fm <- r$feature_metadata
    expect_equal(as.vector(table(fm$featureType[fm$mediator])[c("a", "b")]), c(2L, 2L))
    expect_false(any(fm$mediator[fm$featureType == "c"]))
    expect_equal(sum(fm$indirect), tr$NIE, tolerance = 1e-8)
    expect_true(all(r$sample_metadata$exposure %in% c(0, 1)))
    expect_true(all(is.finite(r$sample_metadata$Y)))
    ## Negative direction flips the sign
    neg <- gen_mediation(nsample = 100, template = tm, n.mediators = 2, TIE = 0.5,
        dir.b = -1, n.mc = 500, seed = 2)
    expect_equal(neg$replicates[[1]]$truth$mediation$NIE, -0.5, tolerance = 1e-8)
})

test_that("mediation effects are recovered by product of coefficients", {
    tm <- fit_template(.mk_template(81, n = 150, p = c(a = 10, b = 10)), rank = 2,
        shared.rank = 1, seed = 1)
    g <- gen_mediation(nsample = 2000, template = tm, n.mediators = 2, TIE = 1, DE = 0.5,
        n.mc = 2000, nrep = 3, seed = 3)
    est <- vapply(g$replicates, function(r) {
        sm <- r$sample_metadata
        x <- t(as.matrix(r$feature_table))
        M <- scale(x[, r$feature_metadata$mediator])
        a <- apply(M, 2, function(m) stats::coef(stats::lm(m ~ sm$exposure))[2])
        fit <- stats::lm(sm$Y ~ sm$exposure + M)
        c(nie = sum(a * stats::coef(fit)[-(1:2)]), de = unname(stats::coef(fit)[2]))
    }, numeric(2))
    expect_lt(abs(mean(est["nie", ]) - 1), 0.15)
    expect_lt(abs(mean(est["de", ]) - 0.5), 0.15)
})

test_that("null mediation scenarios have zero true indirect effect", {
    tm <- fit_template(.mk_template(82), rank = 2, shared.rank = 1, seed = 1)
    for (et in c("null_a", "null_b", "null_ab")) {
        g <- gen_mediation(nsample = 100, template = tm, n.mediators = 2, TIE = 1,
            effect.type = et, n.mc = 500, seed = 4)
        fm <- g$replicates[[1]]$feature_metadata
        expect_equal(g$replicates[[1]]$truth$mediation$NIE, 0)
        if (et == "null_a") expect_true(all(fm$a_path[fm$mediator] == 0))
        if (et == "null_b") expect_true(all(fm$b_path[fm$mediator] == 0))
    }
})

test_that("mediation supports continuous exposure, random directions and interactions", {
    tm <- fit_template(.mk_template(83), rank = 2, shared.rank = 1, seed = 1)
    g <- gen_mediation(nsample = 150, template = tm, n.mediators = 3, TIE = 0.6,
        exposure.type = "continuous", dir.a = "random", n.interactions = 2,
        n.mc = 500, seed = 5)
    tr <- g$replicates[[1]]$truth$mediation
    expect_equal(tr$pathway_effect_total, 0.6, tolerance = 1e-8)
    expect_equal(nrow(tr$interaction_pairs), 2)
    expect_false(all(g$replicates[[1]]$sample_metadata$exposure %in% c(0, 1)))
    expect_error(gen_mediation(nsample = 50, template = tm, n.mediators = 1, TIE = 5,
        n.mc = 500), "too large")
    expect_error(gen_mediation(nsample = 50, template = tm, n.mediators = c(a = 0, b = 0),
        n.mc = 500), "At least one mediator")
    expect_error(gen_mediation(nsample = 50, template = tm, dir.a = 2, n.mc = 500), "dir.a")
})

test_that("simulate_mediation returns an MAE with mediation truth", {
    testthat::skip_if_not_installed("MultiAssayExperiment")
    tm <- fit_template(.mk_template(84), rank = 2, shared.rank = 1, seed = 1)
    mae <- simulate_mediation(tm, n = 80, crosstalk = 1,
        mediation = list(n_mediators = c(a = 2, b = 0), tie = 0.4, n_mc = 500))
    cd <- MultiAssayExperiment::colData(mae)
    expect_true(all(c("exposure", "Y") %in% names(cd)))
    expect_equal(S4Vectors::metadata(mae)$truth$mediation$NIE, 0.4, tolerance = 1e-8)
    rd <- SummarizedExperiment::rowData(MultiAssayExperiment::experiments(mae)[["a"]])
    expect_equal(sum(rd$mediator), 2)
    expect_error(simulate_mediation(tm, n = 50, crosstalk = 1, mediation = list(bad = 1)),
        "unsupported")
})

test_that("mediation with a binary outcome hits prevalence and reports risk-scale truth", {
    tm <- fit_template(.mk_template(85), rank = 2, shared.rank = 1, seed = 1)
    g <- gen_mediation(nsample = 2000, template = tm, n.mediators = 2, TIE = 0.8, DE = 0.4,
        outcome.type = "binary", outcome.prevalence = 0.3, n.mc = 2000, seed = 6)
    r <- g$replicates[[1]]
    tr <- r$truth$mediation
    expect_true(all(r$sample_metadata$Y %in% c(0, 1)))
    expect_lt(abs(mean(r$sample_metadata$Y) - 0.3), 0.04)
    expect_equal(tr$NIE, 0.8, tolerance = 1e-8)
    expect_identical(tr$scale, "conditional log-odds (linear predictor)")
    expect_equal(tr$binary$log_odds_ratio, log(tr$binary$odds_ratio))
    rd <- tr$binary$risk_difference
    expect_equal(unname(rd["TE"]), unname(rd["NIE"] + rd["NDE"]), tolerance = 1e-10)
    expect_true(all(tr$binary$odds_ratio > 1))
    expect_error(gen_mediation(nsample = 50, template = tm, outcome.type = "binary",
        outcome.prevalence = 1, n.mc = 500), "outcome.prevalence")
})

test_that("continuous exposure reports the linear-projection indirect effect", {
    tm <- fit_template(.mk_template(86), rank = 2, shared.rank = 1, seed = 1)
    g <- gen_mediation(nsample = 100, template = tm, n.mediators = 2, TIE = 0.5,
        exposure.type = "continuous", n.mc = 1000, seed = 7)
    tr <- g$replicates[[1]]$truth$mediation
    fm <- g$replicates[[1]]$feature_metadata
    expect_true(is.finite(tr$NIE_linear))
    expect_equal(tr$NIE_linear, sum(fm$b_path * fm$a_slope), tolerance = 1e-8)
    gb <- gen_mediation(nsample = 100, template = tm, n.mediators = 2, TIE = 0.5,
        n.mc = 1000, seed = 7)
    expect_equal(gb$replicates[[1]]$truth$mediation$NIE_linear, 0.5, tolerance = 1e-8)
})

test_that("binary exposure keeps mediator latent scores exactly standard normal", {
    set.seed(90)
    z <- matrix(stats::rnorm(20000), 20000, 1)
    t_raw <- stats::rbinom(20000, 1, 0.1)
    mu <- 0.1
    sdt <- sqrt(0.1 * 0.9)
    support <- list(values = (c(0, 1) - mu) / sdt, probs = c(0.9, 0.1))
    tstar <- (t_raw - mu) / sdt
    des <- MimeSys:::.mediation_design(0.85, 1L, 1L, diag(1), list(), 0.3)
    hook <- MimeSys:::.mediation_stage_hook(des, 1L, 1L, function(r, s) tstar, support)
    out <- hook(list(z))[[1]][, 1]
    expect_gt(stats::ks.test(out, "pnorm")$p.value, 0.01)
    expect_gt(stats::cor(out, t_raw), 0.3)
    ## Two stages: a linked stage-2 mediator is standard normal too, also with
    ## correlated template scores
    z2 <- 0.4 * z + sqrt(1 - 0.4^2) * matrix(stats::rnorm(20000), 20000, 1)
    Sigma <- matrix(c(1, 0.4, 0.4, 1), 2, 2)
    des2 <- MimeSys:::.mediation_design(c(0.5, 0.4), c(1L, 2L), 2L, Sigma,
        list(matrix(1, 1, 1)), 0.3)
    hook2 <- MimeSys:::.mediation_stage_hook(des2, c(1L, 2L), c(1L, 1L),
        function(r, s) tstar, support)
    out2 <- hook2(list(z, z2))
    expect_gt(stats::ks.test(out2[[2]][, 1], "pnorm")$p.value, 0.01)
    expect_gt(stats::cor(out2[[1]][, 1], out2[[2]][, 1]), 0.4)
})

test_that("mediation truth covers closure effects, matched nulls and random interactions", {
    set.seed(91)
    n <- 120
    latent <- stats::rnorm(n)
    lat <- outer(latent, stats::rnorm(15)) + matrix(stats::rnorm(n * 15), n, 15)
    comp <- exp(lat)
    comp <- sweep(comp, 1, rowSums(comp), "/")
    other <- outer(latent, stats::rnorm(8)) + matrix(stats::rnorm(n * 8), n, 8)
    x <- list(comp = comp, other = other)
    for (v in names(x)) rownames(x[[v]]) <- paste0("S", seq_len(n))
    tm <- fit_template(x, rank = 2, shared.rank = 1, seed = 1)
    g <- gen_mediation(nsample = 50, template = tm, n.mediators = c(comp = 2, other = 2),
        TIE = 1, n.mc = 1000, seed = 1)
    fm <- g$replicates[[1]]$feature_metadata
    non_med_comp <- fm$featureType == "comp" & !fm$mediator
    non_med_other <- fm$featureType == "other" & !fm$mediator
    expect_true(any(abs(fm$exposure_effect[non_med_comp]) > 1e-6))
    expect_true(all(abs(fm$exposure_effect[non_med_other]) < 1e-10))
    expect_equal(fm$exposure_effect[fm$mediator], fm$a_effect[fm$mediator])

    sp <- gen_mediation(nsample = 50, template = tm, n.mediators = c(comp = 0, other = 3),
        TIE = 0.8, n.mc = 1000, seed = 2)
    na <- gen_mediation(nsample = 50, template = tm, n.mediators = c(comp = 0, other = 3),
        TIE = 0.8, effect.type = "null_a", n.mc = 1000, seed = 2)
    fs <- sp$replicates[[1]]$feature_metadata
    fa <- na$replicates[[1]]$feature_metadata
    expect_equal(fa$b_path[fa$mediator], fs$b_path[fs$mediator])
    expect_equal(na$replicates[[1]]$truth$mediation$NIE, 0)

    ri <- gen_mediation(nsample = 50, template = tm, n.mediators = c(comp = 0, other = 2),
        TIE = 0.5, dir.b = "random", n.interactions = 1, n.mc = 1000, seed = 3)
    expect_true(all(ri$replicates[[1]]$truth$mediation$interaction_pairs$coefficient != 0))
})

test_that("simulate_mediation returns data on the scale its truth refers to", {
    testthat::skip_if_not_installed("MultiAssayExperiment")
    tm <- fit_template(.mk_template(92), rank = 2, shared.rank = 1, seed = 1)
    mae <- simulate_mediation(tm, n = 60, crosstalk = 1,
        mediation = list(n_mediators = 2, tie = 0.5, n_mc = 500))
    md <- S4Vectors::metadata(mae)$simulation
    expect_equal(sum(names(md) == "crosstalk"), 1L)
    expect_identical(S4Vectors::metadata(mae)$truth$mediation$feature_scale, "fitted")
    rd <- SummarizedExperiment::rowData(MultiAssayExperiment::experiments(mae)[["a"]])
    expect_true(all(rd$scale == "fitted"))
    expect_error(simulate_mediation(tm, n = "a", crosstalk = 1), "n must be")
})

test_that("MAE experiments are matched through sampleMap primary IDs", {
    data(stelzer_mae, package = "MimeSys")
    e <- MultiAssayExperiment::experiments(stelzer_mae)
    a <- e[[1]][, rev(seq_len(ncol(e[[1]])))]
    b <- e[[2]]
    colnames(a) <- paste0("rna_", colnames(a))
    colnames(b) <- paste0("met_", colnames(b))
    sm <- S4Vectors::DataFrame(
        assay = factor(rep(c("A", "B"), c(ncol(a), ncol(b)))),
        primary = c(sub("rna_", "", colnames(a)), sub("met_", "", colnames(b))),
        colname = c(colnames(a), colnames(b)))
    mae <- MultiAssayExperiment::MultiAssayExperiment(
        MultiAssayExperiment::ExperimentList(list(A = a, B = b)),
        colData = MultiAssayExperiment::colData(stelzer_mae), sampleMap = sm)
    v <- as_mimesys_view_list(mae)
    expect_identical(rownames(v$A), rownames(v$B)[match(rownames(v$A), rownames(v$B))])
    expect_setequal(rownames(v$A), rownames(MultiAssayExperiment::colData(stelzer_mae)))
    tm <- fit_template(mae, shared.rank = 1)
    expect_equal(tm$n0, ncol(b))
})

test_that("sim_from_template is reproducible and internal callers keep drawing", {
    tm <- fit_template(.mk_template(81), rank = 2, shared.rank = 1)
    s1 <- sim_from_template(tm, n = 20)
    s2 <- sim_from_template(tm, n = 20)
    expect_identical(s1$feature_table, s2$feature_table)
    s3 <- sim_from_template(tm, n = 20, seed = 2)
    expect_false(identical(s1$feature_table, s3$feature_table))
    nul <- simulate_mimesys_null(tm, n = 20, crosstalk = 0.5, nrep = 2)
    expect_s3_class(nul, "mimesys_mae_simulation")
    x1 <- SummarizedExperiment::assay(MultiAssayExperiment::experiments(nul$replicates[[1]])[[1]])
    x2 <- SummarizedExperiment::assay(MultiAssayExperiment::experiments(nul$replicates[[2]])[[1]])
    expect_false(identical(x1, x2))
    raw <- simulate_mimesys_null(tm, n = 20, crosstalk = 0.5, return_type = "simulation")
    expect_length(raw$replicates, 1)
    cal <- calibrate_crosstalk(tm, crosstalk_grid = c(0, 1), nrep = 3)
    d <- cal$raw$rescaled_frobenius[cal$raw$crosstalk == 1 & !cal$raw$diagonal_block]
    expect_gt(stats::sd(d), 0)
})

test_that("inputs are validated", {
    tm <- fit_template(.mk_template(82), rank = 2, shared.rank = 1)
    expect_error(sim_from_template(tm, n = 2.5), "integer")
    expect_error(sim_from_template(tm, n = 10, feature.jitter = 2), "feature.jitter")
    expect_error(sim_from_template(tm, n = 10, p = c(a = 5, b = 5, zz = 3)), "unknown view")
    expect_error(simulate_mimesys(tm, n = 20, nrep = 1.5, crosstalk = 1), "nrep")
    expect_error(fit_template(.mk_template(82), nperm = 0), "nperm")
    expect_error(fit_template(.mk_template(82), var.explained = 0), "var.explained")
    expect_error(gen_simmba(nsample = 30, template = tm, ygen.mode = "Interaction",
        n.interactions = 1.5, nrep = 1), "n.interactions")
    expect_error(calibrate_crosstalk(tm, p = 1, nrep = 1), "at least two features")
    expect_message(simulate_mimesys_null(tm, n = 20, p = 1), "Skipping crosstalk calibration")
})

test_that("two-part model is restricted to non-negative data", {
    tm <- fit_template(.mk_template(83), rank = 2, shared.rank = 1)
    expect_error(sim_from_template(tm, n = 10, view.models = c(a = "two_part_copula")),
        "non-negative")
})

test_that("TIE = 0 gives a zero target", {
    data(stelzer_mae, package = "MimeSys")
    tm <- fit_template(stelzer_mae, shared.rank = 1)
    med <- gen_mediation(nsample = 40, template = tm, p = 10, n.mediators = 2,
        TIE = 0, DE = 0.2, n.mc = 500, nrep = 1)
    expect_identical(med$replicates[[1]]$truth$mediation$TIE_target, 0)
})

test_that("MAE experiments stored as matrices are transposed", {
    data(stelzer_mae, package = "MimeSys")
    e <- MultiAssayExperiment::experiments(stelzer_mae)
    mats <- lapply(as.list(e), SummarizedExperiment::assay)
    mae <- MultiAssayExperiment::MultiAssayExperiment(
        MultiAssayExperiment::ExperimentList(mats),
        colData = MultiAssayExperiment::colData(stelzer_mae))
    v <- as_mimesys_view_list(mae)
    expect_equal(dim(v[[1]]), rev(dim(mats[[1]])))
    expect_error(gen_simmba(nsample = 50.5, template = fit_template(mae,
        shared.rank = 1), nrep = 1), "nsample")
    expect_error(fit_template(mae, shared.rank.runs = 2.5), "shared.rank.runs")
})

test_that("a single stage reproduces parallel mediation", {
    data(stelzer_mae, package = "MimeSys")
    tm <- fit_template(stelzer_mae, shared.rank = 2)
    args <- list(nsample = 50, template = tm, p = 12,
        n.mediators = c(CYTOF = 2, Metabolomics = 2, Proteomics = 0), TIE = 0.5, n.mc = 800)
    a <- do.call(gen_mediation, args)$replicates[[1]]
    b <- do.call(gen_mediation, c(args, list(stages = list(c("CYTOF", "Metabolomics")))))$replicates[[1]]
    expect_identical(a$feature_table, b$feature_table)
    expect_identical(a$sample_metadata, b$sample_metadata)
    expect_equal(a$truth$mediation$NIE, b$truth$mediation$NIE)
})

test_that("sequential mediation decomposes the indirect effect by path", {
    data(franzosa_mae, package = "MimeSys")
    tf <- fit_template(franzosa_mae, view.transform = c(Metabolome = "log1p"), shared.rank = 2)
    st <- list(Microbe = "Microbiome", Metab = "Metabolome")
    m <- gen_mediation(nsample = 6000, template = tf,
        n.mediators = c(Microbiome = 3, Metabolome = 3), stages = st,
        stage.links = list(density = 0.5, strength = 0.3), TIE = 1, DE = 0.5, n.mc = 6000)
    r <- m$replicates[[1]]
    tr <- r$truth$mediation
    expect_equal(nrow(tr$path_effects), 3)
    expect_equal(sum(tr$path_effects$effect), tr$NIE, tolerance = 1e-8)
    expect_equal(tr$NIE, 1, tolerance = 1e-8)
    expect_equal(sum(tr$NIE_by_stage), tr$NIE, tolerance = 1e-8)
    expect_true(all(tr$path_effects$effect > 0))
    fm <- r$feature_metadata
    expect_true(all(fm$n_parents[fm$mediator & fm$stage == 2] >= 1))
    expect_true(all(fm$n_children[fm$mediator & fm$stage == 1] >= 1))
    expect_true(all(tr$links$parent_stage == 1 & tr$links$child_stage == 2))
    ## Stage-2 mediators keep their template distribution
    X <- t(as.matrix(r$feature_table))
    for (id in fm$featureID[fm$mediator & fm$stage == 2]) {
        tfeat <- fm$template_feature[fm$featureID == id]
        ks <- suppressWarnings(stats::ks.test(X[, id],
            tf$fits$Metabolome$qmat[, tfeat])$statistic)
        expect_lt(ks, 0.05)
    }
    ## Product-of-coefficients estimates of the paths
    sm <- r$sample_metadata
    xs <- scale(X[, fm$mediator])
    s1 <- fm$stage[fm$mediator] == 1
    M1 <- xs[, s1, drop = FALSE]
    M2 <- xs[, !s1, drop = FALSE]
    expo <- sm$exposure
    b <- coef(lm(sm$Y ~ expo + M1 + M2))
    a1 <- apply(M1, 2, function(x) coef(lm(x ~ expo))[2])
    f2 <- lapply(seq_len(ncol(M2)), function(j) coef(lm(M2[, j] ~ expo + M1)))
    a2 <- vapply(f2, `[`, numeric(1), 2)
    D <- t(vapply(f2, function(cc) cc[-(1:2)], numeric(ncol(M1))))
    b1 <- b[1 + 1 + seq_len(ncol(M1))]
    b2 <- b[1 + 1 + ncol(M1) + seq_len(ncol(M2))]
    pe <- stats::setNames(tr$path_effects$effect, tr$path_effects$path)
    expect_equal(sum(a1 * b1), pe[["T -> Microbe -> Y"]], tolerance = 0.25)
    expect_equal(sum(a2 * b2), pe[["T -> Metab -> Y"]], tolerance = 0.25)
    expect_equal(sum((D %*% a1) * b2), pe[["T -> Microbe -> Metab -> Y"]], tolerance = 0.3)
})

test_that("sequential scenarios: shares, pure chain, no links, nulls, patterns", {
    data(franzosa_mae, package = "MimeSys")
    tf <- fit_template(franzosa_mae, view.transform = c(Metabolome = "log1p"), shared.rank = 2)
    st <- list("Microbiome", "Metabolome")
    g <- function(...) gen_mediation(nsample = 60, template = tf,
        n.mediators = c(Microbiome = 3, Metabolome = 3), stages = st, TIE = 0.6,
        n.mc = 1000, ...)$replicates[[1]]
    sh <- g(stage.share = c(1, 3))$truth$mediation
    expect_equal(unname(sh$NIE_by_stage), c(0.15, 0.45), tolerance = 1e-8)
    ch <- g(stage.direct = c(1, 0), stage.share = c(0, 1))$truth$mediation
    pe <- stats::setNames(ch$path_effects$effect, ch$path_effects$path)
    expect_equal(pe[["T -> Microbiome -> Metabolome -> Y"]], 0.6, tolerance = 1e-8)
    expect_equal(unname(pe[c(1, 3)]), c(0, 0), tolerance = 1e-10)
    nl <- g(stage.links = list(strength = 0))
    expect_null(nl$truth$mediation$links)
    expect_equal(nl$truth$mediation$path_effects$effect[2], 0, tolerance = 1e-10)
    for (et in c("null_a", "null_b", "null_ab")) {
        z <- g(effect.type = et)$truth$mediation
        expect_equal(z$NIE, 0)
        expect_equal(z$path_effects$effect, c(0, 0, 0), tolerance = 1e-10)
    }
    for (sd in 1:5) {
        hub <- gen_mediation(nsample = 40, template = tf,
            n.mediators = c(Microbiome = 3, Metabolome = 3), stages = st, TIE = 0.6,
            n.mc = 200, seed = sd, stage.links = list(type = "hub", density = 0.34))$replicates[[1]]
        fm <- hub$feature_metadata
        expect_true(all(fm$n_parents[fm$mediator & fm$stage == 2] == 1))
        expect_equal(sum(fm$n_children > 0), 1)
    }
    expect_error(gen_mediation(nsample = 40, template = tf, n.mediators = 2,
        stages = c("Microbiome", "Metabolome")), "must be a list")
    expect_silent(g(stage.links = list(direction = "random")))
    cv <- g(stage.links = list(type = "convergent", density = 0.7))
    fm <- cv$feature_metadata
    expect_true(all(fm$n_parents[fm$mediator & fm$stage == 2] == 2))
    ex <- g(stage.links = list(type = "explicit",
        matrix = list(rbind(c(1, 0, 0), c(0, 1, 0), c(0, 0, -1)))))
    expect_equal(sign(ex$truth$mediation$links$coefficient), c(1, 1, -1))
    expect_error(g(stage.links = list(type = "explicit", matrix = list(diag(2)))), "3 x 3")
    bin <- g(outcome.type = "binary")$truth$mediation$binary
    expect_equal(sum(bin$path_risk_difference), unname(bin$risk_difference["NIE"]),
        tolerance = 1e-10)
    expect_error(gen_mediation(nsample = 60, template = tf, n.mediators = 2,
        stages = list("Microbiome")), "not in any stage")
    expect_error(gen_mediation(nsample = 60, template = tf, n.mediators = 2,
        stages = list("Microbiome", "Nope")), "unknown view")
    expect_error(g(stage.links = list(strength = 1)), "strength")
})

test_that("simulate_mediation passes the stage options", {
    data(franzosa_mae, package = "MimeSys")
    tf <- fit_template(franzosa_mae, shared.rank = 2)
    mae <- simulate_mediation(tf, n = 60, crosstalk = 1, mediation = list(
        n_mediators = c(Microbiome = 2, Metabolome = 2),
        stages = list("Microbiome", "Metabolome"), n_mc = 800, tie = 0.5))
    tr <- S4Vectors::metadata(mae)$truth$mediation
    expect_equal(nrow(tr$path_effects), 3)
    rd <- SummarizedExperiment::rowData(MultiAssayExperiment::experiments(mae)[["Metabolome"]])
    expect_true(all(c("stage", "n_parents") %in% colnames(rd)))
})
