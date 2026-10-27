#' Simulate multi-study multi-omics data with between-study heterogeneity
#'
#' Simulates several studies from one fitted template. All studies measure the
#' same simulated features. Studies can differ in their samples (case mix), in
#' technical batch effects, in which features carry the outcome signal and how
#' strongly, and in their outcome intercept and signal-to-noise ratio. With
#' one study and no heterogeneity, the result is identical to
#' [gen_simmba()] with the same seed.
#'
#' @param n Samples per study: one value for every study (with `nstudy`) or one
#'   value per study.
#' @param template The template(s) to simulate from:
#'   * one fitted `mimesys_template` (see [fit_template()]) or raw data for
#'     one (a `MultiAssayExperiment` or a named list of sample x feature
#'     matrices, fitted with `fit_template(template, ...)`), used for every
#'     study; or
#'   * one template per study: a list of fitted templates with the same views
#'     and features (from [fit_study_templates()]), or a list of raw data sets,
#'     one per study, fitted here with `fit_study_templates(template, ...)`.
#'     Each study then has its own cohort's marginals, factor structure and
#'     cross-view dependence on the same simulated features.
#' @param nstudy Number of studies. Defaults to `length(n)`.
#' @param effect.cor Correlation between studies of the effect sizes of the
#'   common signal features, in \[0, 1\]. A study coefficient is
#'   `mu + sqrt(effect.cor) * (b - mu) + sqrt(1 - effect.cor) * (b_s - mu)`,
#'   with `b` the common draw, `b_s` an independent draw from the same
#'   effect-size distribution and `mu` the mean of that distribution, so the
#'   mean and variance of the coefficients are kept (the two-peaked +/- shape
#'   of the effect sizes is smoothed). The fresh draws have their own signs,
#'   so lowering `effect.cor` also loosens the sign agreement with the shared
#'   direction set by `signal.align`. `1` (default) gives the same effects in
#'   every study.
#' @param support.overlap Fraction of each study's signal features (per view)
#'   that come from the common signal set, in \[0, 1\]. The remaining signal
#'   features are drawn from the other features of that view (following
#'   `signal.align`, like the common draw), with their own effect sizes, so
#'   every study has as many signal features as the common set (fewer only
#'   when a view has fewer non-signal features than needed). `1` (default)
#'   gives the same signal features in every study; `0` gives
#'   study-specific signal features.
#' @param intercept.sd Standard deviation of the study intercepts. Continuous
#'   outcomes: in units of the within-study SD of the signal. Binary outcomes:
#'   on the logit scale; a common intercept is then solved so the pooled
#'   prevalence equals `binary.prevalence`. Survival outcomes: on the
#'   log-hazard scale (a study-specific baseline hazard).
#' @param snr.sd Standard deviation of the study log signal-to-noise ratios:
#'   study `s` uses `snr * exp(e_s)` with `e_s ~ N(0, snr.sd^2)`. Applies to
#'   every outcome type.
#' @param shift.sd Case-mix shift: standard deviation of a study-specific mean
#'   shift of the latent factor scores (shared and view-specific), in units of
#'   the factor SD. The shift changes the features through the template
#'   loadings, keeps their correlation structure, and changes the outcome
#'   through the same biology (it is part of the signal).
#' @param batch.loc.sd,batch.scale.sd Technical batch effects: per study and
#'   feature, the latent normal score `z` becomes `a + b * z` with
#'   `a ~ N(0, batch.loc.sd^2)` and `log(b) ~ N(0, batch.scale.sd^2)` before
#'   it is mapped through the template marginals, so values stay on the
#'   template's support (including zeros of two-part views). Batch effects are
#'   added after the outcome is generated, so they do not carry signal.
#' @param holdout.studies Optional study numbers whose samples all go to the
#'   test split (for example, leave-one-study-out or new-study evaluation).
#'   The other studies are split by `train.fraction`.
#' @param train.fraction Proportion of samples of each non-holdout study
#'   assigned to the training split. May be `1` when `holdout.studies` is
#'   supplied.
#' @param crosstalk Multiplier on template cross-view dependence (see
#'   [gen_simmba()]). With one template per study, a list of one value per
#'   study is also accepted.
#' @param seed Random seed. The global random-number state is restored when the
#'   function returns. Use `NULL` to use the current random state.
#' @param cens.rate Target censoring proportion for survival outcomes (see
#'   [gen_simmba()]), solved within each study: one value for every study, or
#'   one value per study (in study order, or named by study) to give the
#'   studies different follow-up.
#' @param cens.times Optional censoring times (see [gen_simmba()]): one time,
#'   one per sample of all studies (in study order), a function of `n` (called
#'   once per study), or a list with one such entry per study.
#' @inheritParams gen_simmba
#'
#' @details
#' Samples are simulated study by study; study 1 chooses the simulated features
#' and the other studies reuse them. With one template per study, the same
#' named features (and the same feature pairs for features added beyond the
#' template's own, `p` above the template size) are simulated through each
#' study's own template, so real between-cohort differences in distributions
#' and correlations carry over; the heterogeneity arguments add to them.
#' Signal positions under `signal.align` follow study 1's shared loadings.
#' When the templates have different ranks, the latent truth
#' (`truth$scores`, `truth$shared`) is not comparable across studies and is
#' not returned. Each study's case-mix shift and batch effects are drawn right
#' after its latent scores.
#'
#' The outcome is generated on the pooled samples. Features are standardized
#' over all studies (before batch effects), so study differences in the
#' features carry through to the predictor. Friedman features, interaction
#' pairs and shared-latent factors are the same in every study; the
#' heterogeneity in coefficients (`effect.cor`, `support.overlap`) applies to
#' the linear part. With `signal.share`, each view's coefficients are rescaled
#' by one factor for all studies, computed from the pooled within-study
#' variance of that view's predictor (so `effect.cor = 1` still gives the same
#' coefficients in every study).
#'
#' Case-mix shifts move the shared factors even when `crosstalk = 0`; with
#' `outcome.model = "shared_latent"` the outcome then differs between studies
#' without any feature reflecting it.
#'
#' `snr` keeps its [gen_simmba()] meaning within a study. For continuous
#' outcomes the noise variance of study `s` is `V / snr_s`, with `V` the
#' pooled within-study variance of the predictor. For binary and survival
#' outcomes the predictor is centred on its pooled mean, divided by the pooled
#' within-study SD and multiplied by `sqrt(snr_s * pi^2 / 3)` or
#' `sqrt(snr_s * pi^2 / 6)`. Binary outcomes then get the study intercepts
#' and a common intercept solved so the pooled prevalence equals
#' `binary.prevalence`.
#' Survival censoring (`cens.rate`) is solved within each study, so studies
#' can be given different censoring proportions.
#'
#' @return A list in the [gen_simmba()] format: `trainDat` and `testDat` hold
#'   one pooled data set per replicate, with `sample_metadata$study` and
#'   sample IDs prefixed by study (`Study_1_Subject1`). `feature_metadata`
#'   adds one coefficient column per study (`beta_Study_1`, ...),
#'   `common_signal` (in the common signal set) and `beta` (the average over
#'   studies). `truth$studies` lists each study's size, holdout flag,
#'   intercept and SNR (for binary outcomes `truth$outcome_intercept` is the
#'   common intercept the study intercepts are added to); `truth$beta` is the
#'   feature x study coefficient matrix;
#'   `truth$shift` and `truth$batch` hold the drawn case-mix shifts and batch
#'   effects. Use [as_mimesys_mae()] to convert a replicate.
#' @examples
#' data(stelzer_mae)
#' tmpl <- fit_template(stelzer_mae, shared.rank = 2)
#' sim <- gen_multistudy(
#'     n = c(40, 50, 30), template = tmpl, p = 10,
#'     effect.cor = 0.5, support.overlap = 0.7, shift.sd = 0.3,
#'     batch.loc.sd = 0.2, holdout.studies = 3, seed = 1
#' )
#' table(sim$trainDat$Rep_1$sample_metadata$study)
#' table(sim$testDat$Rep_1$sample_metadata$study)
#' sim$trainDat$Rep_1$truth$studies
#' @export
gen_multistudy <- function(n,
    template,
    nstudy = NULL,
    effect.cor = 1,
    support.overlap = 1,
    intercept.sd = 0,
    snr.sd = 0,
    shift.sd = 0,
    batch.loc.sd = 0,
    batch.scale.sd = 0,
    holdout.studies = NULL,
    p = NULL,
    crosstalk = 1,
    feature.select = c("random", "variance"),
    feature.jitter = 0.2,
    view.models = NULL,
    snr = 1,
    train.fraction = 0.7,
    signal.prob = 0.1,
    signal.neg.prob = 0.5,
    effect.meanlog = 1,
    effect.sdlog = 0.4,
    signal.align = NULL,
    signal.share = NULL,
    friedman.share = NULL,
    outcome.model = c("linear", "friedman", "friedman2", "shared_latent", "interaction"),
    shared.latent.function = c("linear", "friedman", "friedman2"),
    shared.latent.n.factors = NULL,
    shared.latent.weight = 1,
    nonlinear.weight = 0.5,
    interaction.weight = 0.5,
    n.interactions = NULL,
    outcome.type = c("continuous", "binary", "survival"),
    binary.prevalence = 0.5,
    baseline.hazard = 1,
    baseline.shape = 1,
    cens.rate = NULL,
    cens.lower = 1,
    cens.upper = 3,
    cens.type = c("uniform", "exponential", "administrative"),
    cens.times = NULL,
    cens.dependence = 0,
    cens.accrual = 0.5,
    cens.dropout = 0,
    keep.true.times = FALSE,
    nrep = 1,
    seed = 1234,
    return.template = FALSE,
    back.transform = FALSE,
    ...) {

    if (!is.null(seed)) withr::local_seed(seed)
    ## One template for all studies, or one per study (fitted, or raw data
    ## fitted here with fit_study_templates())
    tmpl_list <- .study_template_list(template, list(...))
    if (!is.null(tmpl_list) && is.null(nstudy) && length(n) == 1L) {
        nstudy <- length(tmpl_list)
    }
    n_by_study <- .study_sizes(n, nstudy)
    S <- length(n_by_study)
    if (!is.null(tmpl_list) && length(tmpl_list) != S) {
        stop("template has ", length(tmpl_list), " study templates but there are ", S,
            " studies; give one template per study.", call. = FALSE)
    }
    if (is.list(crosstalk)) {
        if (is.null(tmpl_list) || length(crosstalk) != S) {
            stop("A list of crosstalk values needs one template per study and one ",
                "value per study.", call. = FALSE)
        }
    }
    study_names <- if (!is.null(tmpl_list)) names(tmpl_list) else paste0("Study_", seq_len(S))
    het <- .check_heterogeneity(list(effect.cor = effect.cor,
        support.overlap = support.overlap, intercept.sd = intercept.sd,
        snr.sd = snr.sd, shift.sd = shift.sd, batch.loc.sd = batch.loc.sd,
        batch.scale.sd = batch.scale.sd))
    holdout <- .check_holdout_studies(holdout.studies, S)
    .check_fraction(train.fraction, "train.fraction", include.lower = FALSE)
    if (train.fraction >= 1 && length(holdout) == 0L) {
        stop("train.fraction = 1 needs holdout.studies; otherwise no sample is left for testing.",
            call. = FALSE)
    }
    for (s in setdiff(seq_len(S), holdout)) {
        n_train <- round(n_by_study[s] * train.fraction)
        if (n_train < 1L || (train.fraction < 1 && n_train > n_by_study[s] - 1L)) {
            stop("train.fraction = ", train.fraction, " with n = ", n_by_study[s],
                " in ", study_names[s], " leaves an empty training or test set.",
                call. = FALSE)
        }
    }
    setup <- .simmba_setup(if (is.null(tmpl_list)) template else tmpl_list[[1]],
        list(...), p = p, crosstalk = if (is.list(crosstalk)) crosstalk[[1]] else crosstalk,
        feature.select = feature.select, feature.jitter = feature.jitter,
        view.models = view.models, snr = snr, signal.prob = signal.prob,
        signal.neg.prob = signal.neg.prob, effect.meanlog = effect.meanlog,
        effect.sdlog = effect.sdlog, signal.align = signal.align,
        signal.share = signal.share, friedman.share = friedman.share,
        outcome.model = outcome.model,
        shared.latent.function = shared.latent.function,
        shared.latent.n.factors = shared.latent.n.factors,
        shared.latent.weight = shared.latent.weight,
        nonlinear.weight = nonlinear.weight,
        interaction.weight = interaction.weight, n.interactions = n.interactions,
        outcome.type = outcome.type, binary.prevalence = binary.prevalence,
        baseline.hazard = baseline.hazard, cens.rate = cens.rate,
        cens.lower = cens.lower, cens.upper = cens.upper,
        cens.type = cens.type, cens.times = cens.times,
        cens.dependence = cens.dependence, cens.accrual = cens.accrual,
        cens.dropout = cens.dropout, baseline.shape = baseline.shape,
        keep.true.times = keep.true.times,
        back.transform = back.transform)
    tmpl <- setup$template
    cfg <- setup$cfg
    cfg$crosstalk <- crosstalk
    cfg$cens.rate <- .per_study_values(cfg$cens.rate, study_names, "cens.rate")
    if (is.list(cfg$cens.times) && length(cfg$cens.times) != S) {
        stop("A list of cens.times needs one entry per study (", S, ").", call. = FALSE)
    }
    if (!is.null(tmpl_list)) {
        if (identical(cfg$ygen.mode, "SharedLatent") &&
            length(unique(vapply(tmpl_list, function(t) t$shared.rank, integer(1)))) > 1L) {
            stop("outcome.model = \"shared_latent\" with one template per study needs ",
                "the same shared.rank in every template; fit them with a fixed ",
                "shared.rank (e.g. fit_study_templates(x, shared.rank = 2)).", call. = FALSE)
        }
        tmpl <- tmpl_list
    }
    nrep <- .check_count(nrep, "nrep", min = 1L)

    trainDat <- testDat <- agreement.emp <- vector("list", nrep)
    names(trainDat) <- names(testDat) <- names(agreement.emp) <- paste("Rep",
        seq_len(nrep), sep = "_")
    for (k in seq_len(nrep)) {
        res <- .simulate_outcome_replicate(tmpl, n_by_study, cfg, het = het,
            study_names = study_names)
        pcl <- res$pcl
        pcl$truth$studies <- data.frame(study = study_names, n = n_by_study,
            holdout = seq_len(S) %in% holdout, intercept = res$intercept,
            snr = res$snr, stringsAsFactors = FALSE)
        if (!is.null(res$outcome_intercept)) {
            pcl$truth$outcome_intercept <- res$outcome_intercept
        }
        pcl$truth$beta <- matrix(res$beta, ncol = S,
            dimnames = list(pcl$feature_metadata$featureID, study_names))
        if (het$shift.sd > 0) {
            pcl$truth$shift <- stats::setNames(res$shift, study_names)
        }
        if (het$batch.loc.sd > 0 || het$batch.scale.sd > 0) {
            pcl$truth$batch <- list(
                loc = .batch_matrix(res$batch, "loc", pcl$feature_metadata$featureID, study_names),
                scale = .batch_matrix(res$batch, "scale", pcl$feature_metadata$featureID, study_names))
        }
        agreement.emp[[k]] <- res$agreement

        ## Train / test split within each study; holdout studies are test only
        tr.row <- te.row <- integer(0)
        for (s in seq_len(S)) {
            rows <- which(res$study == s)
            if (s %in% holdout) {
                te.row <- c(te.row, rows)
                next
            }
            tr <- sample.int(length(rows), size = round(length(rows) * train.fraction),
                replace = FALSE)
            tr.row <- c(tr.row, rows[tr])
            te.row <- c(te.row, rows[-tr])
        }
        trainDat[[k]] <- .subset_pcl_samples(pcl, tr.row)
        testDat[[k]] <- .subset_pcl_samples(pcl, te.row)
    }

    list(trainDat = trainDat, testDat = testDat,
        template = if (isTRUE(return.template)) tmpl else NULL,
        template_provenance = if (is.null(tmpl_list)) tmpl$provenance else
            lapply(tmpl_list, `[[`, "provenance"),
        p = p, crosstalk = crosstalk,
        view.models = view.models,
        cancor = if (is.null(tmpl_list)) .realized_cancor(tmpl, crosstalk) else
            stats::setNames(lapply(seq_len(S), function(s) .realized_cancor(tmpl_list[[s]],
                if (is.list(crosstalk)) crosstalk[[s]] else crosstalk)), study_names),
        agreement.emp = agreement.emp,
        nstudy = S, n = n_by_study, holdout.studies = holdout,
        effect.cor = effect.cor, support.overlap = support.overlap,
        intercept.sd = intercept.sd, snr.sd = snr.sd, shift.sd = shift.sd,
        batch.loc.sd = batch.loc.sd, batch.scale.sd = batch.scale.sd,
        snr = snr, train.fraction = train.fraction,
        signal.prob = cfg$signal.prob, signal.neg.prob = cfg$signal.neg.prob,
        effect.meanlog = cfg$effect.meanlog, effect.sdlog = cfg$effect.sdlog,
        signal.align = signal.align, signal.share = cfg$signal.share,
        friedman.share = cfg$friedman.share,
        shared.latent.function = cfg$shared.latent.function,
        shared.latent.n.factors = cfg$shared.latent.n.factors,
        nonlinear.weight = nonlinear.weight, back.transform = back.transform,
        interaction.weight = interaction.weight,
        n.interactions = n.interactions,
        nrep = nrep, seed = seed, outcome.model = cfg$outcome.model,
        outcome.type = cfg$outcome.type, binary.prevalence = binary.prevalence,
        baseline.hazard = baseline.hazard, cens.rate = cens.rate,
        cens.lower = cens.lower, cens.upper = cens.upper) |>
        .add_nondefault_settings(cfg)
}

## ---- multi-study helpers -------------------------------------------------

## A setting given once or once per study (in study order or named by study)
.per_study_values <- function(x, study_names, name) {
    if (is.null(x) || length(x) == 1L) return(x)
    S <- length(study_names)
    if (length(x) != S) {
        stop(name, " has ", length(x), " values for ", S,
            " studies; give one value or one per study.", call. = FALSE)
    }
    if (!is.null(names(x))) {
        if (!setequal(names(x), study_names)) {
            stop("The names of ", name, " must be the study names: ",
                paste(study_names, collapse = ", "), ".", call. = FALSE)
        }
        x <- x[study_names]
    }
    unname(x)
}

## Heterogeneity draws use a side stream: a seed derived from the current
## random state (without advancing it) plus a tag. The main stream, and so the
## base data, coefficients and outcome noise, is then the same at every
## heterogeneity level for a given seed.
.with_side_stream <- function(tag, expr) {
    if (!exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
        stats::runif(1)
    }
    state <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    side <- (sum(as.numeric(state) %% 2147483647) + tag) %% 2147483647
    withr::with_seed(as.integer(side), expr)
}

.study_sizes <- function(n, nstudy) {
    if (!is.numeric(n) || length(n) < 1L || any(!is.finite(n)) ||
        any(abs(n - round(n)) > 1e-8) || any(n < 2)) {
        stop("n must contain whole numbers of at least 2.", call. = FALSE)
    }
    n <- as.integer(round(n))
    if (is.null(nstudy)) return(n)
    nstudy <- .check_count(nstudy, "nstudy", min = 1L)
    if (length(n) == 1L) return(rep(n, nstudy))
    if (length(n) != nstudy) {
        stop("n must be a single value or one value per study (length nstudy = ",
            nstudy, ").", call. = FALSE)
    }
    n
}

.check_heterogeneity <- function(het) {
    .check_fraction(het$effect.cor, "effect.cor")
    .check_fraction(het$support.overlap, "support.overlap")
    for (name in c("intercept.sd", "snr.sd", "shift.sd", "batch.loc.sd",
        "batch.scale.sd")) {
        x <- het[[name]]
        if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 0) {
            stop(name, " must be a single non-negative number.", call. = FALSE)
        }
    }
    het
}

.check_holdout_studies <- function(holdout, S) {
    if (is.null(holdout) || length(holdout) == 0L) return(integer(0))
    if (!is.numeric(holdout) || any(!is.finite(holdout)) ||
        any(abs(holdout - round(holdout)) > 1e-8) || any(holdout < 1) ||
        any(holdout > S) || anyDuplicated(holdout)) {
        stop("holdout.studies must be distinct study numbers between 1 and ", S, ".",
            call. = FALSE)
    }
    holdout <- sort(as.integer(round(holdout)))
    if (length(holdout) == S) {
        stop("holdout.studies cannot contain every study; at least one study must ",
            "supply training samples.", call. = FALSE)
    }
    holdout
}

## Coefficients of every study (feature x study). With effect.cor = 1 and
## support.overlap = 1 every column is the common draw and no random numbers
## are used. Study-specific signal features follow signal.align like the
## common draw (picked by squared shared loading, signed by the loading).
## Effect sizes are mixed around the mean of the effect-size distribution, so
## the mean and variance of the coefficients are kept and their correlation
## between studies is effect.cor.
.study_coefficients <- function(beta0, feature_type, views, S, het, cfg,
    loading = rep(0, length(beta0))) {
    B <- matrix(beta0, nrow = length(beta0), ncol = S)
    if (het$effect.cor >= 1 && het$support.overlap >= 1) return(B)
    r <- het$effect.cor
    for (s in seq_len(S)) {
        for (m in seq_along(views)) {
            idx <- which(feature_type == views[m])
            common <- idx[beta0[idx] != 0]
            k <- length(common)
            if (k == 0L) next
            draw_effects <- function(size) {
                log2(.lnorm_factors(size, cfg$signal.neg.prob[m],
                    cfg$effect.meanlog[m], cfg$effect.sdlog[m]))
            }
            b <- numeric(length(beta0))
            keep <- common
            if (het$support.overlap < 1) {
                n_keep <- as.integer(round(het$support.overlap * k))
                keep <- .sample_values(common, n_keep)
                if (n_keep == 0L) keep <- integer(0)
                others <- setdiff(idx, common)
                n_new <- min(k - n_keep, length(others))
                if (n_new > 0L) {
                    n_al <- if (is.null(cfg$signal.align)) 0L else
                        as.integer(round(cfg$signal.align * n_new))
                    pos_al <- if (n_al > 0L) .sample_values(others, n_al,
                        prob = loading[others]^2 + 1e-12) else integer(0)
                    rest <- setdiff(others, pos_al)
                    new_pos <- c(pos_al, .sample_values(rest, n_new - n_al))
                    b[new_pos] <- draw_effects(n_new)
                    if (n_al > 0L) {
                        sg <- sign(loading[pos_al])
                        sg[sg == 0] <- sign(b[pos_al])[sg == 0]
                        b[pos_al] <- abs(b[pos_al]) * sg
                    }
                }
            }
            if (length(keep) > 0L) {
                mu <- .effect_mean(cfg$signal.neg.prob[m], cfg$effect.meanlog[m],
                    cfg$effect.sdlog[m])
                b[keep] <- if (r >= 1) beta0[keep] else
                    mu + sqrt(r) * (beta0[keep] - mu) +
                        sqrt(1 - r) * (draw_effects(length(keep)) - mu)
            }
            B[idx, s] <- b[idx]
        }
    }
    B
}

## Case-mix shift of study latent factor scores. The shared factors move by
## d_sh and the view-specific innovations by d_E, so each view's scores move by
## d_sh B_m' sqrt(gamma_m) + d_E S_m (S_m = sqrt(I - gamma_m B_m B_m')) and the
## latent feature scores by that shift times the loadings.
.draw_case_mix_shift <- function(Z_list, ctx, sd) {
    M <- length(Z_list)
    d_shared <- stats::rnorm(ctx$r, 0, sd)
    d_scores <- vector("list", M)
    for (m in seq_len(M)) {
        d_E <- stats::rnorm(ctx$ks[m], 0, sd)
        if (ctx$r == 0 || ctx$gamma[m] == 0) {
            d_m <- d_E
        } else {
            Bm <- sqrt(ctx$gamma[m]) * ctx$B[[m]]
            d_m <- as.numeric(Bm %*% d_shared +
                .sqrt_psd(diag(ctx$ks[m]) - Bm %*% t(Bm)) %*% d_E)
        }
        d_scores[[m]] <- d_m
        Z_list[[m]] <- sweep(Z_list[[m]], 2,
            as.numeric(ctx$lam[[m]] %*% d_m), "+")
    }
    list(Z = Z_list, shift = list(shared = d_shared, scores = unlist(d_scores)))
}

## Batch effects on the latent normal scores: z -> a + b z per feature
.draw_batch_effects <- function(Z_list, loc_sd, scale_sd) {
    loc <- scale <- vector("list", length(Z_list))
    for (m in seq_along(Z_list)) {
        p_m <- ncol(Z_list[[m]])
        a <- if (loc_sd > 0) stats::rnorm(p_m, 0, loc_sd) else rep(0, p_m)
        b <- if (scale_sd > 0) exp(stats::rnorm(p_m, 0, scale_sd)) else rep(1, p_m)
        Z_list[[m]] <- sweep(sweep(Z_list[[m]], 2, b, "*"), 2, a, "+")
        loc[[m]] <- a
        scale[[m]] <- b
    }
    list(Z = Z_list, batch = list(loc = unlist(loc), scale = unlist(scale)))
}

## Mean of a spike-in coefficient log2(.lnorm_factors()). The coefficient is
## d * |L| with L = log(f) / log(2) ~ N(meanlog, sdlog^2) / log(2) and an
## independent sign d (-1 with probability neg.prob), so its mean is
## (1 - 2 neg.prob) E|L| (folded normal mean).
.effect_mean <- function(neg.prob, meanlog, sdlog) {
    m <- meanlog / log(2)
    s <- sdlog / log(2)
    folded <- if (s > 0) s * sqrt(2 / pi) * exp(-m^2 / (2 * s^2)) +
        m * (1 - 2 * stats::pnorm(-m / s)) else abs(m)
    (1 - 2 * neg.prob) * folded
}

.batch_matrix <- function(batch, which, feature_ids, study_names) {
    matrix(unlist(lapply(batch, `[[`, which)), ncol = length(study_names),
        dimnames = list(feature_ids, study_names))
}

## Samples `rows` of a pooled pcl (feature table, metadata and latent truth)
.subset_pcl_samples <- function(pcl, rows) {
    out <- pcl
    out$sample_metadata <- pcl$sample_metadata[rows, , drop = FALSE]
    out$feature_table <- pcl$feature_table[, rows, drop = FALSE]
    if (!is.null(pcl$truth$scores)) out$truth$scores <- pcl$truth$scores[rows, , drop = FALSE]
    if (!is.null(pcl$truth$shared)) out$truth$shared <- pcl$truth$shared[rows, , drop = FALSE]
    if (!is.null(pcl$truth$specific_scores)) {
        out$truth$specific_scores <- pcl$truth$specific_scores[rows, , drop = FALSE]
    }
    out
}

#' Simulate multi-study multi-omics data as a MultiAssayExperiment
#'
#' MAE-first interface to [gen_multistudy()], with the same presets, outcome
#' and signal options as [simulate_mimesys()]. All studies of a replicate are
#' returned in one `MultiAssayExperiment`; `colData()` holds the `study` and
#' `split` of every sample.
#'
#' @inheritParams simulate_mimesys
#' @param template One template (or data set) for every study, or one per
#'   study: a list of templates from [fit_study_templates()] or a list of data
#'   sets (fitted with `fit_study_templates(template, ...)`, using
#'   `template_args`). See [gen_multistudy()].
#' @param n Samples per study: one value for every study (with `nstudy`) or one
#'   value per study.
#' @param nstudy Number of studies. Defaults to `length(n)`.
#' @param heterogeneity Named list of between-study differences. Options are
#'   `effect_cor`, `support_overlap`, `intercept_sd`, `snr_sd`, `shift_sd`,
#'   `batch_loc_sd`, `batch_scale_sd` and `holdout_studies` (see
#'   `effect.cor`, `support.overlap`, `intercept.sd`, `snr.sd`, `shift.sd`,
#'   `batch.loc.sd`, `batch.scale.sd` and `holdout.studies` in
#'   [gen_multistudy()]). The defaults give homogeneous studies.
#' @param calibration Named list controlling automatic crosstalk calibration
#'   (see [simulate_mimesys()]). With one template, calibration uses the
#'   average study size; with one template per study, each study's crosstalk
#'   is calibrated on its own template and size.
#' @param train_fraction Fraction of each non-holdout study's samples labelled
#'   `"train"`.
#' @param return_type `"mae"` (default) returns MAE-oriented output. Use
#'   `"simulation"` to return the underlying [gen_multistudy()] result.
#'
#' @return A `MultiAssayExperiment`, a `mimesys_mae_simulation` list (several
#'   replicates), or the [gen_multistudy()] result when
#'   `return_type = "simulation"`. The study-level truth is in
#'   `metadata(mae)$truth`.
#' @examples
#' data(stelzer_mae)
#' tmpl <- fit_template(stelzer_mae, shared.rank = 2)
#' mae <- simulate_multistudy(tmpl, n = c(40, 50, 30), p = 10, crosstalk = 1,
#'     heterogeneity = list(effect_cor = 0.5, shift_sd = 0.3,
#'         holdout_studies = 3), seed = 1)
#' table(MultiAssayExperiment::colData(mae)$study,
#'     MultiAssayExperiment::colData(mae)$split)
#' @export
simulate_multistudy <- function(template,
    n,
    nstudy = NULL,
    p = NULL,
    crosstalk = "calibrate",
    preset = "default",
    outcome = NULL,
    signal = list(),
    heterogeneity = list(),
    null_model = list(),
    template_args = list(),
    calibration = list(),
    train_fraction = 0.7,
    nrep = 1,
    seed = 1234,
    return_type = c("mae", "simulation")) {
    return_type <- match.arg(return_type)
    n_by_study <- .study_sizes(n, nstudy)
    .check_count(nrep, "nrep", min = 1L)
    if (!is.null(seed)) withr::local_seed(seed)
    opts <- .mimesys_outcome_signal_options(preset, outcome, signal)
    outcome <- opts$outcome
    signal <- opts$signal
    het <- .mimesys_merge_options(heterogeneity, list(
        effect_cor = 1, support_overlap = 1, intercept_sd = 0, snr_sd = 0,
        shift_sd = 0, batch_loc_sd = 0, batch_scale_sd = 0,
        holdout_studies = NULL
    ), "heterogeneity")

    tmpl_list <- .study_template_list(template, template_args)
    if (!is.null(tmpl_list) && is.null(nstudy) && length(n) == 1L) {
        n_by_study <- .study_sizes(n, length(tmpl_list))
    }
    setup <- .mimesys_wrapper_setup(
        if (is.null(tmpl_list)) template else tmpl_list[[1]],
        n = if (is.null(tmpl_list)) max(2L, as.integer(round(mean(n_by_study)))) else n_by_study[1],
        p = p, nrep = nrep,
        crosstalk = crosstalk, null_model = null_model,
        template_args = template_args, calibration = calibration, seed = seed,
        null_preset = opts$preset$null_model)
    fitted_template <- setup$template
    null_model <- setup$null_model
    resolved <- setup$crosstalk
    if (!is.null(tmpl_list)) {
        if (length(tmpl_list) != length(n_by_study)) {
            stop("template has ", length(tmpl_list), " study templates but n has ",
                length(n_by_study), " studies.", call. = FALSE)
        }
        ## Crosstalk resolved (calibrated) for each study's template
        res_list <- c(list(resolved), lapply(seq_along(tmpl_list)[-1], function(s) {
            .mimesys_resolve_crosstalk(tmpl_list[[s]], n = n_by_study[s], p = p,
                crosstalk = crosstalk, calibration = calibration,
                null_model = null_model, seed = seed)
        }))
        names(res_list) <- names(tmpl_list)
        fitted_template <- tmpl_list
        resolved <- list(value = lapply(res_list, `[[`, "value"),
            calibration = lapply(res_list, `[[`, "calibration"))
    }
    sim <- gen_multistudy(
        n = n_by_study, template = fitted_template,
        effect.cor = het$effect_cor, support.overlap = het$support_overlap,
        intercept.sd = het$intercept_sd, snr.sd = het$snr_sd,
        shift.sd = het$shift_sd, batch.loc.sd = het$batch_loc_sd,
        batch.scale.sd = het$batch_scale_sd,
        holdout.studies = het$holdout_studies,
        p = p, crosstalk = resolved$value,
        feature.select = null_model$feature_select,
        feature.jitter = null_model$feature_jitter,
        view.models = null_model$view_models,
        back.transform = isTRUE(null_model$back_transform),
        snr = outcome$snr, train.fraction = train_fraction,
        signal.prob = signal$prob, signal.neg.prob = signal$neg_prob,
        effect.meanlog = signal$effect_meanlog, effect.sdlog = signal$effect_sdlog,
        signal.align = signal$align, signal.share = signal$share,
        friedman.share = outcome$friedman_share, outcome.model = outcome$model,
        shared.latent.function = outcome$shared_latent_function,
        shared.latent.n.factors = outcome$shared_latent_n_factors,
        shared.latent.weight = outcome$shared_latent_weight,
        nonlinear.weight = outcome$nonlinear_weight,
        interaction.weight = outcome$interaction_weight,
        n.interactions = outcome$n_interactions,
        outcome.type = outcome$type, binary.prevalence = outcome$prevalence,
        baseline.hazard = outcome$baseline_hazard, cens.rate = outcome$cens_rate,
        cens.lower = outcome$cens_lower, cens.upper = outcome$cens_upper,
        cens.type = outcome$cens_type, cens.times = outcome$cens_times,
        cens.dependence = outcome$cens_dependence,
        cens.accrual = outcome$cens_accrual, cens.dropout = outcome$cens_dropout,
        baseline.shape = outcome$baseline_shape,
        keep.true.times = isTRUE(outcome$keep_true_times),
        nrep = as.integer(nrep), seed = seed
    )
    sim$crosstalk_calibration <- resolved$calibration
    sim$provenance <- .mimesys_wrapper_provenance(
        if (is.null(tmpl_list)) fitted_template else list(provenance = lapply(tmpl_list, `[[`, "provenance")),
        resolved, purpose = "multistudy", preset = preset, seed = seed)
    if (identical(return_type, "simulation")) return(sim)

    maes <- lapply(seq_len(sim$nrep), function(i) as_mimesys_mae(sim, rep = i))
    .mimesys_collect_maes(maes, .mimesys_simulation_settings(sim))
}

#' Fit one template per study on common views and features
#'
#' Fits a [fit_template()] template to each study's data, restricted to the
#' views and features that every study has, so the templates can be used
#' together in [gen_multistudy()] and [simulate_multistudy()] (one template
#' per study). Features that a fit drops in any study (all missing, or
#' constant) are removed from every study and the templates are refitted, so
#' all templates end up with the same features in the same order.
#'
#' @param x A list with one data set per study: `MultiAssayExperiment`s,
#'   `SummarizedExperiment`s or named lists of sample x feature matrices (see
#'   [as_mimesys_view_list()]). Names, if given, label the studies.
#' @param ... Passed to [fit_template()] for every study (for example
#'   `view.transform`, `rank` or `shared.rank`).
#' @param min.features Smallest number of common features a view must keep.
#'
#' @return A `mimesys_template_list`: a named list of `mimesys_template`s with
#'   identical views and features.
#' @examples
#' data(stelzer_mae)
#' ## Two "cohorts" from one data set, for illustration
#' half <- split(colnames(stelzer_mae[[1]]), rep(1:2, length.out = ncol(stelzer_mae[[1]])))
#' cohorts <- lapply(half, function(ids) stelzer_mae[, ids])
#' tmpls <- fit_study_templates(cohorts, shared.rank = 1)
#' tmpls
#' sim <- gen_multistudy(n = c(30, 30), template = tmpls, p = 10, seed = 1)
#' table(sim$trainDat$Rep_1$sample_metadata$study)
#' @export
fit_study_templates <- function(x, ..., min.features = 2L) {
    if (!is.list(x) || inherits(x, "mimesys_template") || length(x) < 2L) {
        stop("x must be a list with one data set per study (at least two).",
            call. = FALSE)
    }
    min.features <- .check_count(min.features, "min.features", min = 1L)
    study_names <- names(x)
    if (is.null(study_names) || any(!nzchar(study_names)) || anyDuplicated(study_names)) {
        study_names <- paste0("Study_", seq_along(x))
    }
    view_lists <- lapply(x, as_mimesys_view_list)
    views <- Reduce(intersect, lapply(view_lists, names))
    if (length(views) == 0L) {
        stop("The studies have no view in common (views are matched by name).",
            call. = FALSE)
    }
    all_views <- unique(unlist(lapply(view_lists, names)))
    if (length(all_views) > length(views)) {
        message("Using the views common to every study: ", paste(views, collapse = ", "),
            "; dropped: ", paste(setdiff(all_views, views), collapse = ", "), ".")
    }
    usable <- function(X) {
        X <- as.matrix(X)
        colnames(X)[apply(X, 2, function(col) {
            col <- col[is.finite(col)]
            length(col) > 0L && length(unique(col)) > 1L
        })]
    }
    features <- stats::setNames(lapply(views, function(v) {
        Reduce(intersect, lapply(view_lists, function(l) usable(l[[v]])))
    }), views)
    for (pass in seq_len(5L)) {
        short <- views[lengths(features) < min.features]
        if (length(short) > 0L) {
            stop("Too few features common to every study in view(s) ",
                paste(short, collapse = ", "), " (need ", min.features,
                "). Check that feature names match across studies.", call. = FALSE)
        }
        fits <- lapply(view_lists, function(l) {
            fit_template(lapply(stats::setNames(views, views), function(v) {
                l[[v]][, features[[v]], drop = FALSE]
            }), ...)
        })
        kept <- stats::setNames(lapply(views, function(v) {
            Reduce(intersect, lapply(fits, function(f) colnames(f$fits[[v]]$qmat)))
        }), views)
        if (identical(kept, features)) break
        features <- kept
        if (pass == 5L) {
            stop("Could not find features kept by every study's template.", call. = FALSE)
        }
    }
    dropped <- sum(vapply(views, function(v) {
        length(unique(unlist(lapply(view_lists, function(l) colnames(l[[v]]))))) -
            length(features[[v]])
    }, numeric(1)))
    if (dropped > 0) {
        message("Kept ", sum(lengths(features)), " features common to every study (",
            paste(sprintf("%s: %d", views, lengths(features)), collapse = ", "),
            "); dropped ", dropped, ".")
    }
    names(fits) <- study_names
    structure(fits, class = c("mimesys_template_list", "list"))
}

#' @export
print.mimesys_template_list <- function(x, ...) {
    first <- x[[1]]
    cat("MimeSys templates for", length(x), "studies\n")
    cat("  Views (features): ", paste(sprintf("%s (%d)", first$views,
        vapply(first$fits, function(f) ncol(f$qmat), integer(1))), collapse = ", "),
        "\n", sep = "")
    for (s in names(x)) {
        cat("  ", s, ": ", x[[s]]$n0, " samples, shared rank ", x[[s]]$shared.rank,
            ", view ranks ", paste(vapply(x[[s]]$fits, function(f) f$rank, integer(1)),
                collapse = "/"), "\n", sep = "")
    }
    invisible(x)
}

## The per-study templates in `template`, or NULL for a single template.
## Raw data sets (one per study) are fitted with fit_study_templates().
.study_template_list <- function(template, fit_args = list()) {
    if (inherits(template, "mimesys_template") ||
        inherits(template, "MultiAssayExperiment") ||
        inherits(template, "SummarizedExperiment") || !is.list(template) ||
        length(template) == 0L) {
        return(NULL)
    }
    is_tmpl <- vapply(template, inherits, logical(1), "mimesys_template")
    is_data <- vapply(template, function(e) {
        inherits(e, "MultiAssayExperiment") || inherits(e, "SummarizedExperiment") ||
            (is.list(e) && !is.data.frame(e))
    }, logical(1))
    if (all(is_tmpl)) {
        tl <- template
    } else if (all(is_data)) {
        tl <- do.call(fit_study_templates, c(list(template), fit_args))
    } else if (any(is_tmpl) || any(is_data)) {
        stop("template must be one template (or data set) for every study, or a ",
            "list with one template (or data set) per study.", call. = FALSE)
    } else {
        return(NULL)   # a named list of matrices: one data set
    }
    if (is.null(names(tl)) || any(!nzchar(names(tl))) || anyDuplicated(names(tl))) {
        names(tl) <- paste0("Study_", seq_along(tl))
    }
    .check_study_templates(tl)
    tl
}

.check_study_templates <- function(tl) {
    ref <- tl[[1]]
    for (s in seq_along(tl)[-1]) {
        t <- tl[[s]]
        if (!identical(t$views, ref$views)) {
            stop("Study templates must have the same views in the same order; ",
                names(tl)[s], " differs. Use fit_study_templates().", call. = FALSE)
        }
        same_features <- vapply(seq_along(ref$views), function(m) {
            identical(colnames(t$fits[[m]]$qmat), colnames(ref$fits[[m]]$qmat))
        }, logical(1))
        if (!all(same_features)) {
            stop("Study templates must have the same features in every view; ",
                names(tl)[s], " differs in ", paste(ref$views[!same_features],
                    collapse = ", "), ". Use fit_study_templates().", call. = FALSE)
        }
        vt <- function(x) if (is.null(x$view.transform)) "none" else x$view.transform
        if (!identical(vt(t), vt(ref))) {
            stop("Study templates must use the same view.transform; ", names(tl)[s],
                " differs.", call. = FALSE)
        }
    }
    invisible(TRUE)
}
