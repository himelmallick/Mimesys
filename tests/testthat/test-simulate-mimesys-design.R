data(stelzer_mae)
tm_d <- fit_template(stelzer_mae, shared.rank = 1, seed = 1)

strip_time <- function(x) {
    x$template_provenance <- NULL
    x$provenance <- NULL
    x
}

test_that("each design gives the result of its engine", {
    sim <- simulate_mimesys(tm_d, n = 40, p = 6, crosstalk = 0.5,
        outcome = list(type = "binary"), return_type = "simulation", seed = 2)
    eng <- gen_simmba(40, tm_d, p = 6, crosstalk = 0.5, outcome.type = "binary",
        signal.align = 1, seed = 2)
    expect_identical(sim$trainDat, eng$trainDat)

    ms <- simulate_mimesys(tm_d, design = "multistudy", n = c(30, 20), p = 6,
        crosstalk = 0.5, heterogeneity = list(shift_sd = 0.2),
        return_type = "simulation", seed = 3)
    eng <- gen_multistudy(c(30, 20), tm_d, p = 6, crosstalk = 0.5, shift.sd = 0.2,
        signal.align = 1, seed = 3)
    expect_identical(ms$trainDat, eng$trainDat)

    nul <- simulate_mimesys(tm_d, design = "null", n = 30, p = 6, crosstalk = 0.5,
        return_type = "simulation", seed = 4)
    expect_named(nul, c("replicates", "settings", "template"))
    expect_null(nul$replicates$Rep_1$sample_metadata$Y)

    med <- simulate_mimesys(tm_d, design = "mediation", n = 60, p = 6, crosstalk = 1,
        mediation = list(n_mediators = 2, n_mc = 300), return_type = "simulation", seed = 5)
    eng <- gen_mediation(60, tm_d, p = 6, crosstalk = 1, n.mediators = 2, n.mc = 300,
        seed = 5)
    expect_identical(med$replicates, eng$replicates)
})

test_that("options a design does not use are errors", {
    expect_error(simulate_mimesys(tm_d, n = 30, p = 5, crosstalk = 1,
        mediation = list(n_mediators = 2)), "does not use `mediation`")
    expect_error(simulate_mimesys(tm_d, n = 30, p = 5, crosstalk = 1, design = "null",
        outcome = "linear"), "does not use `outcome`")
    expect_error(simulate_mimesys(tm_d, n = 30, p = 5, crosstalk = 1, design = "mediation",
        heterogeneity = list(shift_sd = 1)), "heterogeneity")
    expect_error(simulate_mimesys(tm_d, n = 30, p = 5, crosstalk = 1,
        heterogeneity = list(shift_sd = 1)), "design \"multistudy\"")
    expect_error(simulate_mimesys(tm_d, n = c(30, 40), p = 5, crosstalk = 1),
        "design = \"multistudy\"")
    expect_error(simulate_mimesys(tm_d, n = 30, p = 5, crosstalk = 1, design = "bad"),
        "should be one of")
})

test_that("the default design is outcome", {
    a <- simulate_mimesys(tm_d, n = 40, p = 6, crosstalk = 0.5, seed = 6)
    b <- simulate_mimesys(tm_d, n = 40, p = 6, crosstalk = 0.5, design = "outcome", seed = 6)
    expect_identical(MultiAssayExperiment::colData(a), MultiAssayExperiment::colData(b))
    expect_true("Y" %in% names(MultiAssayExperiment::colData(a)))
})
