## Multi-study simulation (gen_multistudy(), simulate_multistudy())

ms_template <- function() {
    data(stelzer_mae, package = "MimeSys", envir = environment())
    fit_template(stelzer_mae, shared.rank = 2, seed = 1)
}

pooled <- function(sim, rep = 1) {
    rbind(sim$trainDat[[rep]]$sample_metadata, sim$testDat[[rep]]$sample_metadata)
}

test_that("one homogeneous study reproduces gen_simmba()", {
    tm <- ms_template()
    settings <- list(
        list(outcome.type = "continuous"),
        list(outcome.type = "binary", binary.prevalence = 0.3,
            outcome.model = "friedman2"),
        list(outcome.type = "survival", cens.rate = 0.3,
            outcome.model = "interaction", signal.share = c(1, 1, 1))
    )
    for (a in settings) {
        s1 <- do.call(gen_simmba, c(list(n = 50, template = tm, p = 15, seed = 4), a))
        s2 <- do.call(gen_multistudy, c(list(n = 50, template = tm, p = 15, seed = 4), a))
        for (split in c("trainDat", "testDat")) {
            d1 <- s1[[split]]$Rep_1
            d2 <- s2[[split]]$Rep_1
            expect_identical(unname(as.matrix(d1$feature_table)),
                unname(as.matrix(d2$feature_table)))
            cols <- intersect(c("Xbeta", "Y", "time", "status"),
                names(d1$sample_metadata))
            expect_identical(d1$sample_metadata[cols],
                `rownames<-`(d2$sample_metadata[cols], rownames(d1$sample_metadata)))
            expect_identical(d1$feature_metadata$beta, d2$feature_metadata$beta_Study_1)
        }
    }
})

test_that("studies are labelled, share features and split as requested", {
    tm <- ms_template()
    sim <- gen_multistudy(n = c(40, 50, 30), template = tm, p = 10,
        holdout.studies = 3, seed = 1)
    tr <- sim$trainDat$Rep_1
    te <- sim$testDat$Rep_1
    expect_equal(as.vector(table(tr$sample_metadata$study)), c(28, 35))
    expect_equal(as.vector(table(te$sample_metadata$study)), c(12, 15, 30))
    ids <- c(tr$sample_metadata$subjectID, te$sample_metadata$subjectID)
    expect_false(anyDuplicated(ids) > 0)
    expect_true(all(startsWith(ids, paste0(c(tr$sample_metadata$study,
        te$sample_metadata$study), "_"))))
    expect_identical(colnames(tr$feature_table), tr$sample_metadata$subjectID)
    expect_equal(nrow(tr$truth$shared), nrow(tr$sample_metadata))
    expect_equal(sim$trainDat$Rep_1$truth$studies$holdout, c(FALSE, FALSE, TRUE))
    expect_equal(gen_multistudy(n = 20, nstudy = 3, template = tm, p = 5,
        seed = 1)$n, c(20L, 20L, 20L))

    ## train.fraction = 1 with a held-out study: all non-holdout samples train
    s1 <- gen_multistudy(n = c(20, 25), template = tm, p = 5,
        holdout.studies = 2, train.fraction = 1, seed = 2)
    expect_equal(nrow(s1$trainDat$Rep_1$sample_metadata), 20)
    expect_true(all(s1$testDat$Rep_1$sample_metadata$study == "Study_2"))

    mae <- as_mimesys_mae(sim)
    cd <- MultiAssayExperiment::colData(mae)
    expect_equal(as.vector(table(cd$study, cd$split)["Study_3", "test"]), 30)
    expect_true(is.matrix(S4Vectors::metadata(mae)$truth$beta))
})

test_that("batch effects change the features but not the biology or the outcome", {
    tm <- ms_template()
    b0 <- gen_multistudy(n = c(60, 60), template = tm, p = 20, seed = 4)
    b1 <- gen_multistudy(n = c(60, 60), template = tm, p = 20, seed = 4,
        batch.loc.sd = 1, batch.scale.sd = 0.3)
    expect_identical(pooled(b0)$Xbeta, pooled(b1)$Xbeta)
    expect_identical(pooled(b0)$Y, pooled(b1)$Y)
    expect_false(identical(b0$trainDat$Rep_1$feature_table,
        b1$trainDat$Rep_1$feature_table))
    bt <- b1$trainDat$Rep_1$truth$batch
    expect_equal(dim(bt$loc), c(60, 2))
    expect_true(all(bt$scale > 0))
    ## Values stay within the template support
    x0 <- as.matrix(b0$trainDat$Rep_1$feature_table)
    x1 <- as.matrix(b1$trainDat$Rep_1$feature_table)
    for (j in seq_len(nrow(x0))[1:5]) {
        tf <- b1$trainDat$Rep_1$feature_metadata$template_feature[j]
        v <- b1$trainDat$Rep_1$feature_metadata$featureType[j]
        q <- range(tm$fits[[v]]$qmat[, tf])
        expect_true(all(x1[j, ] >= q[1] - 1e-8 & x1[j, ] <= q[2] + 1e-8))
    }
})

test_that("case-mix shifts move study means and enter the outcome", {
    tm <- ms_template()
    sim <- gen_multistudy(n = 400, nstudy = 3, template = tm, p = 20,
        shift.sd = 1, outcome.model = "shared_latent", seed = 6)
    d <- pooled(sim)
    sh <- sim$trainDat$Rep_1$truth$shift
    expect_named(sh, c("Study_1", "Study_2", "Study_3"))
    expect_length(sh$Study_1$shared, 2)
    ## Study means of the outcome predictor follow the shared-factor shifts
    m <- tapply(d$Xbeta, d$study, mean)
    target <- vapply(sh, function(x) mean(x$shared), numeric(1))
    expect_gt(stats::cor(m, target), 0.9)
    ## The shared latent truth includes the shift
    shared <- rbind(sim$trainDat$Rep_1$truth$shared, sim$testDat$Rep_1$truth$shared)
    st <- c(sim$trainDat$Rep_1$sample_metadata$study, sim$testDat$Rep_1$sample_metadata$study)
    expect_equal(unname(colMeans(shared[st == "Study_3", ])), sh$Study_3$shared,
        tolerance = 0.2)
})

test_that("effect correlation and support overlap control the study coefficients", {
    tm <- ms_template()
    same <- gen_multistudy(n = 30, nstudy = 3, template = tm, p = 40, seed = 2)
    B <- same$trainDat$Rep_1$truth$beta
    expect_true(all(B[, 1] == B[, 2] & B[, 2] == B[, 3]))

    sim <- gen_multistudy(n = 30, nstudy = 3, template = tm, p = 600,
        signal.prob = 0.3, effect.cor = 0.5, support.overlap = 0.6, seed = 2)
    fm <- sim$trainDat$Rep_1$feature_metadata
    B <- sim$trainDat$Rep_1$truth$beta
    ## Same number of signal features in every study as in the common set
    expect_equal(unname(colSums(B != 0)), rep(sum(fm$common_signal), 3))
    for (s in 1:3) {
        expect_equal(mean(fm$common_signal[B[, s] != 0]), 0.6, tolerance = 0.02)
    }
    both <- fm$common_signal & B[, 1] != 0 & B[, 2] != 0
    expect_equal(stats::cor(B[both, 1], B[both, 2]), 0.5, tolerance = 0.15)
    expect_equal(fm$beta, unname(rowMeans(B)))

    ind <- gen_multistudy(n = 30, nstudy = 2, template = tm, p = 100,
        support.overlap = 0, seed = 3)
    Bi <- ind$trainDat$Rep_1$truth$beta
    expect_false(any(Bi[, 1] != 0 & ind$trainDat$Rep_1$feature_metadata$common_signal))
})

test_that("study coefficients keep alignment, means and signal.share balance", {
    tm <- ms_template()
    ## signal.share uses one factor per view, so effect.cor = 1 stays identical
    sh <- gen_multistudy(n = c(30, 30), template = tm, p = 20,
        signal.share = c(1, 1, 1), seed = 1)
    B <- sh$trainDat$Rep_1$truth$beta
    expect_identical(B[, 1], B[, 2])
    ## New study-specific signal features follow signal.align = 1
    al <- gen_multistudy(n = 30, nstudy = 3, template = tm, p = 200,
        signal.align = 1, support.overlap = 0.3, seed = 2)
    fm <- al$trainDat$Rep_1$feature_metadata
    Ba <- al$trainDat$Rep_1$truth$beta
    new <- Ba != 0 & !fm$common_signal
    expect_true(any(new))
    expect_true(all(sign(Ba[new]) == sign(matrix(fm$shared_loading, nrow(Ba), 3)[new])))
    ## Mixing keeps the mean effect when effects are all positive
    m <- vapply(c(1, 0.3), function(r) {
        x <- gen_multistudy(n = 20, nstudy = 4, template = tm, p = 300,
            signal.prob = 0.3, signal.neg.prob = 0, effect.cor = r, seed = 3)
        Bx <- x$trainDat$Rep_1$truth$beta
        mean(Bx[x$trainDat$Rep_1$feature_metadata$common_signal, ])
    }, numeric(1))
    expect_equal(m[2], m[1], tolerance = 0.1)
})

test_that("binary prevalence is calibrated on the pooled sample with study intercepts", {
    tm <- ms_template()
    sim <- gen_multistudy(n = 1500, nstudy = 4, template = tm, p = 15,
        outcome.type = "binary", binary.prevalence = 0.3, intercept.sd = 1,
        seed = 3)
    d <- pooled(sim)
    expect_equal(mean(stats::plogis(d$Xbeta)), 0.3, tolerance = 1e-6)
    expect_true(is.numeric(sim$trainDat$Rep_1$truth$outcome_intercept))
    prev <- tapply(d$Y, d$study, mean)
    a <- sim$trainDat$Rep_1$truth$studies$intercept
    expect_gt(stats::cor(as.vector(prev), a), 0.9)
})

test_that("study SNRs vary with snr.sd and keep their within-study meaning", {
    tm <- ms_template()
    sim <- gen_multistudy(n = 3000, nstudy = 3, template = tm, p = 15,
        snr = 1, snr.sd = 0.7, seed = 5)
    d <- pooled(sim)
    realized <- vapply(split(d, d$study), function(x) {
        stats::var(x$Xbeta) / stats::var(x$Y - x$Xbeta)
    }, numeric(1))
    expect_equal(unname(realized), sim$trainDat$Rep_1$truth$studies$snr,
        tolerance = 0.1)
})

test_that("survival censoring is solved within each study", {
    tm <- ms_template()
    sim <- gen_multistudy(n = 300, nstudy = 3, template = tm, p = 15,
        outcome.type = "survival", cens.rate = 0.25, intercept.sd = 0.5, seed = 7)
    d <- pooled(sim)
    expect_equal(as.vector(tapply(1 - d$status, d$study, mean)), rep(0.25, 3),
        tolerance = 0.01)
})

test_that("nonlinear terms are shared across studies", {
    tm <- ms_template()
    fr <- gen_multistudy(n = 40, nstudy = 3, template = tm, p = 15,
        outcome.model = "friedman", effect.cor = 0, support.overlap = 0, seed = 9)
    expect_equal(sum(fr$trainDat$Rep_1$feature_metadata$friedman_term), 5)
    expect_true(all(fr$trainDat$Rep_1$truth$beta == 0))
    it <- gen_multistudy(n = 40, nstudy = 2, template = tm, p = 15,
        outcome.model = "interaction", signal.share = c(1, 0, 1), seed = 9)
    fm <- it$trainDat$Rep_1$feature_metadata
    expect_true(all(it$trainDat$Rep_1$truth$beta[fm$featureType == tm$views[2], ] == 0))
    expect_equal(sum(fm$interaction_term), 6)
})

test_that("two-part views keep zeros with batch effects and back-transform", {
    data(franzosa_mae, package = "MimeSys", envir = environment())
    tf <- fit_template(franzosa_mae, view.transform = c(Metabolome = "log1p"),
        shared.rank = 2, seed = 3)
    mic <- setdiff(tf$views, "Metabolome")[1]
    sim <- suppressWarnings(gen_multistudy(n = 40, nstudy = 2, template = tf,
        p = 30, view.models = stats::setNames("two_part_copula", mic),
        batch.loc.sd = 1, shift.sd = 0.5, back.transform = TRUE, seed = 2))
    d <- sim$trainDat$Rep_1
    x <- as.matrix(d$feature_table[d$feature_metadata$featureType == mic, ])
    expect_true(any(x == 0))
    expect_true(all(x >= 0))
    expect_true(all(d$feature_metadata$scale == "original"))
})

test_that("simulate_multistudy() returns one MAE with study labels", {
    tm <- ms_template()
    mae <- simulate_multistudy(tm, n = c(30, 40), p = 8, crosstalk = 1,
        heterogeneity = list(effect_cor = 0.5, batch_loc_sd = 0.3,
            holdout_studies = 2),
        outcome = list(type = "binary", prevalence = 0.4), seed = 1)
    cd <- MultiAssayExperiment::colData(mae)
    expect_equal(ncol(mae[[1]]), 70)
    expect_true(all(cd$split[cd$study == "Study_2"] == "test"))
    expect_equal(S4Vectors::metadata(mae)$simulation$nstudy, 2L)
    expect_s3_class(simulate_multistudy(tm, n = 20, nstudy = 2, p = 5,
        crosstalk = 1, nrep = 2, seed = 1), "mimesys_mae_simulation")
    expect_error(simulate_multistudy(tm, n = 20, nstudy = 2, p = 5, crosstalk = 1,
        heterogeneity = list(effect_corr = 0.5)), "unsupported option")
})

test_that("gen_multistudy() validates its arguments", {
    tm <- ms_template()
    expect_error(gen_multistudy(n = c(20, 30), nstudy = 3, template = tm, p = 5),
        "one value per study")
    expect_error(gen_multistudy(n = 1, template = tm, p = 5), "at least 2")
    expect_error(gen_multistudy(n = 20, nstudy = 2, template = tm, p = 5,
        effect.cor = 1.5), "effect.cor")
    expect_error(gen_multistudy(n = 20, nstudy = 2, template = tm, p = 5,
        shift.sd = -1), "shift.sd")
    expect_error(gen_multistudy(n = 20, nstudy = 2, template = tm, p = 5,
        holdout.studies = 3), "holdout.studies")
    expect_error(gen_multistudy(n = 20, nstudy = 2, template = tm, p = 5,
        holdout.studies = 1:2), "every study")
    expect_error(gen_multistudy(n = 20, nstudy = 2, template = tm, p = 5,
        train.fraction = 1), "holdout.studies")
})
