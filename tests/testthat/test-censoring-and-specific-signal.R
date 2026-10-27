data(stelzer_mae)
tm_cs <- fit_template(stelzer_mae, shared.rank = 1, seed = 1)

surv_meta <- function(sim) {
    rbind(sim$trainDat$Rep_1$sample_metadata, sim$testDat$Rep_1$sample_metadata)
}

test_that("every censoring type reaches the target censoring proportion", {
    for (ty in c("uniform", "exponential", "administrative")) {
        for (dep in c(0, 1)) {
            sim <- gen_simmba(300, tm_cs, p = 10, outcome.type = "survival",
                cens.rate = 0.35, cens.type = ty, cens.dependence = dep,
                cens.dropout = if (ty == "administrative") 0.2 else 0, seed = 2)
            m <- surv_meta(sim)
            expect_equal(1 - mean(m$status), 0.35, tolerance = 1 / 300 + 1e-8)
        }
    }
})

test_that("default censoring settings are not recorded and give the old draws", {
    sim <- gen_simmba(60, tm_cs, p = 10, outcome.type = "survival", cens.rate = 0.3, seed = 3)
    expect_null(sim$cens.type)
    expect_null(sim$shared.latent.weight)
    ## explicit defaults are the same as leaving them out
    sim2 <- gen_simmba(60, tm_cs, p = 10, outcome.type = "survival", cens.rate = 0.3,
        cens.type = "uniform", cens.dependence = 0, shared.latent.weight = 1, seed = 3)
    sim$template_provenance <- sim2$template_provenance <- NULL
    expect_identical(sim, sim2)
})

test_that("administrative censoring uses a fixed study end", {
    sim <- gen_simmba(200, tm_cs, p = 10, outcome.type = "survival", cens.rate = 0.4,
        cens.type = "administrative", cens.accrual = 0, seed = 4)
    m <- surv_meta(sim)
    ## no accrual and no dropout: every censored subject is censored at the end
    expect_lt(diff(range(m$time[m$status == 0])), 1e-8)
    expect_true(all(m$time[m$status == 1] <= max(m$time[m$status == 0]) + 1e-8))
    expect_identical(sim$cens.type, "administrative")
    expect_identical(sim$cens.accrual, 0)
})

test_that("informative censoring censors higher-risk subjects earlier", {
    cens_cor <- function(dep) {
        m <- surv_meta(gen_simmba(800, tm_cs, p = 10, outcome.type = "survival",
            cens.rate = 0.4, cens.dependence = dep, seed = 5))
        cor(1 - m$status, m$Xbeta)
    }
    expect_gt(cens_cor(1.5), cens_cor(0) + 0.2)
})

test_that("user censoring times are used as given or rescaled to the target", {
    sim <- gen_simmba(100, tm_cs, p = 10, outcome.type = "survival", cens.times = 0.5, seed = 6)
    m <- surv_meta(sim)
    expect_true(all(m$time <= 0.5))
    expect_true(all(m$time[m$status == 0] == 0.5))
    ## one time per sample, Inf = never censored
    tt <- rep(c(Inf, 1), 50)
    sim <- gen_simmba(100, tm_cs, p = 10, outcome.type = "survival", cens.times = tt, seed = 6)
    m <- surv_meta(sim)
    idx <- as.integer(sub("\\D+", "", m$subjectID))
    expect_true(all(m$status[is.infinite(tt[idx])] == 1L))
    expect_true(all(m$time[!is.infinite(tt[idx])] <= 1))
    sim <- gen_simmba(200, tm_cs, p = 10, outcome.type = "survival",
        cens.times = function(n) stats::rgamma(n, 2), cens.rate = 0.25, seed = 6)
    expect_equal(1 - mean(surv_meta(sim)$status), 0.25, tolerance = 1 / 200 + 1e-8)
    expect_error(gen_simmba(100, tm_cs, p = 10, outcome.type = "survival",
        cens.times = c(1, 2, 3), seed = 6), "cens.times has 3 values")
    expect_error(gen_simmba(100, tm_cs, p = 10, outcome.type = "survival",
        cens.times = -1, seed = 6), "positive")
    expect_error(gen_simmba(100, tm_cs, p = 10, outcome.type = "survival",
        cens.times = rep(Inf, 100), cens.rate = 0.3, seed = 6), "cannot be reached")
})

test_that("censoring types other than uniform need cens.rate", {
    expect_error(gen_simmba(60, tm_cs, p = 10, outcome.type = "survival",
        cens.type = "exponential", seed = 1), "needs cens.rate")
    expect_error(gen_simmba(60, tm_cs, p = 10, outcome.type = "survival",
        cens.rate = c(0.1, 0.2), seed = 1), "single number")
    expect_error(gen_simmba(60, tm_cs, p = 10, outcome.type = "survival",
        cens.rate = 0.2, cens.dropout = 1, seed = 1), "cens.dropout")
})

test_that("gen_multistudy takes one censoring proportion per study", {
    sim <- gen_multistudy(c(150, 150, 150), tm_cs, p = 10, outcome.type = "survival",
        cens.rate = c(Study_2 = 0.5, Study_1 = 0.1, Study_3 = 0.3),
        cens.type = "exponential", seed = 7)
    m <- surv_meta(sim)
    rate <- tapply(1 - m$status, m$study, mean)
    expect_equal(as.numeric(rate), c(0.1, 0.5, 0.3), tolerance = 1 / 150 + 1e-8)
    expect_error(gen_multistudy(c(50, 50), tm_cs, p = 10, outcome.type = "survival",
        cens.rate = c(0.1, 0.2, 0.3), seed = 7), "one per study")
    expect_error(gen_multistudy(c(50, 50), tm_cs, p = 10, outcome.type = "survival",
        cens.times = list(1), seed = 7), "one entry per study")
    sim <- gen_multistudy(c(50, 50), tm_cs, p = 10, outcome.type = "survival",
        cens.times = list(0.3, 1), seed = 7)
    m <- surv_meta(sim)
    expect_true(all(m$time[m$study == "Study_1"] <= 0.3))
})

test_that("simulate wrappers pass the new outcome options", {
    sim <- simulate_mimesys(tm_cs, n = 120, p = 10, crosstalk = 1,
        outcome = list(model = "shared_latent", type = "survival", cens_rate = 0.3,
            cens_type = "administrative", cens_dropout = 0.1,
            shared_latent_weight = 0.5),
        return_type = "simulation", seed = 8)
    expect_identical(sim$cens.type, "administrative")
    expect_identical(sim$shared.latent.weight, 0.5)
    expect_false(is.null(sim$trainDat$Rep_1$sample_metadata$Xbeta_specific))
    ms <- simulate_mimesys(design = "multistudy", tm_cs, n = 60, nstudy = 2, p = 10, crosstalk = 1,
        outcome = list(type = "survival", cens_rate = c(0.2, 0.4)),
        return_type = "simulation", seed = 8)
    m <- surv_meta(ms)
    expect_equal(as.numeric(tapply(1 - m$status, m$study, mean)), c(0.2, 0.4),
        tolerance = 1 / 60 + 1e-8)
})

test_that("shared.latent.weight mixes shared and view-specific latent signal", {
    cm <- min(crosstalk_max(tm_cs, per.view = TRUE))
    fit_r2 <- function(w) {
        sim <- gen_simmba(1200, tm_cs, p = 40, crosstalk = cm, outcome.model = "shared_latent",
            shared.latent.weight = w, train.fraction = 0.5, seed = 9)
        d <- sim$trainDat$Rep_1
        md <- d$sample_metadata
        X <- t(as.matrix(d$feature_table))
        ft <- d$feature_metadata$featureType
        list(sim = sim, md = md,
            shared_r2 = stats::cor(md$Xbeta, md$Xbeta_shared)^2,
            view_r2 = vapply(tm_cs$views, function(v)
                summary(stats::lm(md$Xbeta ~ X[, ft == v]))$r.squared, numeric(1)))
    }
    one <- fit_r2(1)
    expect_null(one$md$Xbeta_specific)
    half <- fit_r2(0.5)
    expect_equal(half$shared_r2, 0.5, tolerance = 0.08)
    spec <- half$sim$trainDat$Rep_1$truth$specific_scores
    expect_identical(dim(spec), c(600L, 3L))
    ## view-specific scores are independent across views and of the shared factor
    cc <- stats::cor(spec)
    expect_lt(max(abs(cc[upper.tri(cc)])), 0.1)
    expect_lt(max(abs(stats::cor(spec, half$md$Xbeta_shared))), 0.1)
    none <- fit_r2(0)
    expect_lt(none$shared_r2, 0.05)
    ## each view explains much less of a view-specific outcome than of a shared one
    expect_true(all(none$view_r2 < 0.6))
    expect_true(all(one$view_r2 > 0.7))
    ## signal.share sets each view's part of the view-specific signal
    sim <- gen_simmba(1200, tm_cs, p = 40, crosstalk = cm, outcome.model = "shared_latent",
        shared.latent.weight = 0, signal.share = c(1, 0, 0), seed = 9)
    md <- sim$trainDat$Rep_1$sample_metadata
    spec <- sim$trainDat$Rep_1$truth$specific_scores
    expect_gt(stats::cor(md$Xbeta, spec[, 1]), 0.99)
    expect_error(gen_simmba(60, tm_cs, p = 10, outcome.model = "shared_latent",
        shared.latent.weight = 2, seed = 1), "shared.latent.weight")
})

test_that("view-specific signal works with several studies and per-study templates", {
    sim <- gen_multistudy(c(80, 80), tm_cs, p = 10, outcome.model = "shared_latent",
        shared.latent.weight = 0.3, shift.sd = 0.3, outcome.type = "binary", seed = 10)
    spec <- sim$trainDat$Rep_1$truth$specific_scores
    expect_identical(ncol(spec), 3L)
    expect_identical(rownames(spec), sim$trainDat$Rep_1$sample_metadata$subjectID)
    tl <- fit_study_templates(list(A = stelzer_mae, B = stelzer_mae), shared.rank = 1)
    sim <- gen_multistudy(c(60, 60), tl, p = 10, outcome.model = "shared_latent",
        shared.latent.weight = 0.5, outcome.type = "survival", cens.rate = c(0.2, 0.5),
        seed = 11)
    m <- surv_meta(sim)
    expect_equal(as.numeric(tapply(1 - m$status, m$study, mean)), c(0.2, 0.5),
        tolerance = 1 / 60 + 1e-8)
})

test_that("baseline.shape gives Weibull proportional-hazards event times", {
    for (k in c(0.5, 2)) {
        sim <- gen_simmba(1500, tm_cs, p = 10, outcome.type = "survival",
            baseline.shape = k, cens.rate = 0.3, keep.true.times = TRUE, seed = 12)
        m <- surv_meta(sim)
        fit <- survival::survreg(survival::Surv(m$event_time) ~ m$Xbeta, dist = "weibull")
        ## log T = -(log h0 + eta) / k + W / k: survreg scale 1/k, slope -1/k
        expect_equal(fit$scale, 1 / k, tolerance = 0.06)
        expect_equal(unname(stats::coef(fit)[2]), -1 / k, tolerance = 0.06)
        expect_equal(1 - mean(m$status), 0.3, tolerance = 1 / 1500 + 1e-8)
        expect_identical(sim$baseline.shape, k)
    }
    expect_error(gen_simmba(60, tm_cs, p = 10, outcome.type = "survival",
        baseline.shape = 0, seed = 1), "baseline.shape")
})

test_that("keep.true.times stores the latent event and censoring times", {
    sim <- gen_simmba(200, tm_cs, p = 10, outcome.type = "survival", cens.rate = 0.4,
        cens.type = "administrative", keep.true.times = TRUE, seed = 13)
    m <- surv_meta(sim)
    expect_equal(m$time, pmin(m$event_time, m$cens_time))
    expect_identical(m$status, as.integer(m$cens_time >= m$event_time))
    ## off by default, and switching it on changes nothing else
    sim0 <- gen_simmba(200, tm_cs, p = 10, outcome.type = "survival", cens.rate = 0.4,
        cens.type = "administrative", seed = 13)
    m0 <- surv_meta(sim0)
    expect_null(m0$event_time)
    expect_identical(m0[, c("time", "status", "Xbeta")], m[, c("time", "status", "Xbeta")])
    ms <- simulate_mimesys(design = "multistudy", tm_cs, n = 50, nstudy = 2, p = 10, crosstalk = 1,
        outcome = list(type = "survival", cens_rate = 0.3, baseline_shape = 1.5,
            keep_true_times = TRUE), return_type = "simulation", seed = 13)
    expect_false(is.null(ms$trainDat$Rep_1$sample_metadata$event_time))
    expect_identical(ms$baseline.shape, 1.5)
})
