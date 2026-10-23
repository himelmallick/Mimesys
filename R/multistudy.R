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
#' @param template A fitted `mimesys_template` (see [fit_template()]) or a named
#'   list of sample x feature matrices, which is then fitted with
#'   `fit_template(template, ...)`.
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
#' @param seed Random seed. The global random-number state is restored when the
#'   function returns. Use `NULL` to use the current random state.
#' @inheritParams gen_simmba
#'
#' @details
#' Samples are simulated study by study from the template; study 1 chooses the
#' simulated features and the other studies reuse them. Each study's case-mix
#' shift and batch effects are drawn right after its latent scores.
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
#' Survival censoring (`cens.rate`) is solved within each study.
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
    nonlinear.weight = 0.5,
    interaction.weight = 0.5,
    n.interactions = NULL,
    outcome.type = c("continuous", "binary", "survival"),
    binary.prevalence = 0.5,
    baseline.hazard = 1,
    cens.rate = NULL,
    cens.lower = 1,
    cens.upper = 3,
    nrep = 1,
    seed = 1234,
    return.template = FALSE,
    back.transform = FALSE,
    ...) {

    if (!is.null(seed)) withr::local_seed(seed)
    n_by_study <- .study_sizes(n, nstudy)
    S <- length(n_by_study)
    study_names <- paste0("Study_", seq_len(S))
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
    setup <- .simmba_setup(template, list(...), p = p, crosstalk = crosstalk,
        feature.select = feature.select, feature.jitter = feature.jitter,
        view.models = view.models, snr = snr, signal.prob = signal.prob,
        signal.neg.prob = signal.neg.prob, effect.meanlog = effect.meanlog,
        effect.sdlog = effect.sdlog, signal.align = signal.align,
        signal.share = signal.share, friedman.share = friedman.share,
        outcome.model = outcome.model,
        shared.latent.function = shared.latent.function,
        shared.latent.n.factors = shared.latent.n.factors,
        nonlinear.weight = nonlinear.weight,
        interaction.weight = interaction.weight, n.interactions = n.interactions,
        outcome.type = outcome.type, binary.prevalence = binary.prevalence,
        baseline.hazard = baseline.hazard, cens.rate = cens.rate,
        cens.lower = cens.lower, cens.upper = cens.upper,
        back.transform = back.transform)
    tmpl <- setup$template
    cfg <- setup$cfg
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
        template_provenance = tmpl$provenance, p = p, crosstalk = crosstalk,
        view.models = view.models,
        cancor = .realized_cancor(tmpl, crosstalk),
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
        cens.lower = cens.lower, cens.upper = cens.upper)
}

## ---- multi-study helpers -------------------------------------------------

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
    out$truth$scores <- pcl$truth$scores[rows, , drop = FALSE]
    out$truth$shared <- pcl$truth$shared[rows, , drop = FALSE]
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
#'   (see [simulate_mimesys()]); calibration uses the average study size.
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

    setup <- .mimesys_wrapper_setup(template,
        n = max(2L, as.integer(round(mean(n_by_study)))), p = p, nrep = nrep,
        crosstalk = crosstalk, null_model = null_model,
        template_args = template_args, calibration = calibration, seed = seed,
        null_preset = opts$preset$null_model)
    fitted_template <- setup$template
    null_model <- setup$null_model
    resolved <- setup$crosstalk
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
        nonlinear.weight = outcome$nonlinear_weight,
        interaction.weight = outcome$interaction_weight,
        n.interactions = outcome$n_interactions,
        outcome.type = outcome$type, binary.prevalence = outcome$prevalence,
        baseline.hazard = outcome$baseline_hazard, cens.rate = outcome$cens_rate,
        cens.lower = outcome$cens_lower, cens.upper = outcome$cens_upper,
        nrep = as.integer(nrep), seed = seed
    )
    sim$crosstalk_calibration <- resolved$calibration
    sim$provenance <- .mimesys_wrapper_provenance(fitted_template, resolved,
        purpose = "multistudy", preset = preset, seed = seed)
    if (identical(return_type, "simulation")) return(sim)

    maes <- lapply(seq_len(sim$nrep), function(i) as_mimesys_mae(sim, rep = i))
    .mimesys_collect_maes(maes, .mimesys_simulation_settings(sim))
}
