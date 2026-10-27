#' Generate simulated multi-omics data
#'
#' Generates null multiview data from a real template (see \code{fit_template()}
#' and \code{sim_from_template()}), spikes in signal, generates the outcome, and
#' splits the data into training and testing sets.
#'
#' @param n Number of samples to simulate in each replicate before the
#'   train/test split.
#' @param template Either a \code{mimesys_template} from \code{fit_template()}
#'   (recommended: fit once, reuse) or a named list of sample x feature
#'   matrices,
#'   which is then fitted with \code{fit_template(template, ...)}.
#' @param p Feature counts per view. Use \code{NULL} to keep template sizes, a
#'   single integer to use the same count for every view, an unnamed vector in
#'   template view order, or a named vector with one value per view.
#' @param crosstalk Multiplier on template cross-view dependence. Use a scalar
#'   for the same value in every view or one value per view. A value of 0
#'   removes
#'   shared latent structure, 1 uses the fitted template level, and values above
#'   1 amplify shared structure up to \code{crosstalk_max(template)}.
#' @param feature.select Feature choice rule passed to
#'   \code{sim_from_template()}
#'   when \code{p} is below the template feature count. Use \code{"random"} or
#'   \code{"variance"}.
#' @param feature.jitter Offspring-feature perturbation passed to
#'   \code{sim_from_template()} when \code{p} exceeds the template feature
#'   count.
#'   Values near 0 make offspring closer to their parent features; larger values
#'   add more loading-level perturbation.
#' @param view.models Optional named vector of per-view observation models
#'   passed
#'   to \code{sim_from_template()}. Use `"two_part_copula"` for sparse
#'   compositional microbiome views; unspecified views use `"continuous"`.
#' @param snr Signal-to-noise ratio. For continuous outcomes it is
#'   \code{Var(eta) / Var(noise)}. For binary and survival outcomes the linear
#'   predictor is standardized and rescaled so that \code{snr} has the same
#'   latent-scale meaning: \code{Var(eta) / (pi^2 / 3)} for the logistic latent
#'   variable and \code{Var(eta) / (pi^2 / 6)} for the log event time. Larger
#'   values create easier prediction tasks.
#' @param train.fraction Proportion of samples assigned to the training split.
#' @param signal.prob Probability that a feature receives a nonzero spike-in
#'   coefficient; scalar, one value per view, or named by view. Every view with
#'   \code{signal.prob > 0} gets at least one nonzero coefficient; use 0 for a view
#'   that should carry no linear signal.
#' @param signal.neg.prob Probability that a nonzero spike-in coefficient is
#'   negative; scalar or one value per view.
#' @param effect.meanlog,effect.sdlog Mean and SD on the log scale of the
#'   log-normal effect-size factor \eqn{f}; a signal coefficient is
#'   \eqn{\pm |\log_2 f|}. Scalar or one value per view.
#' @param signal.align NULL draws signal positions and signs independently per
#'   view. A value between 0 and 1 places that fraction of the signal features
#'   on features loading on the shared cross-view direction (sampled with
#'   probability proportional to the squared loading), with signs following
#'   that direction, so view-specific signals agree to the extent the views
#'   cross-talk.
#' @param signal.share Optional numeric vector controlling the relative marginal
#'   variance contribution of each view's linear predictor. If NULL, spike-in
#'   coefficients are left on their original scale. If supplied, values are
#'   recycled or matched by view name, normalized to sum to one, and each view's
#'   coefficients are rescaled so its view-specific linear predictor has
#'   variance
#'   equal to the requested share. For example, \code{c(1, 1, 1)} creates a
#'   balanced complementary-signal setting where no view is intended to
#'   dominate.
#' @param friedman.share Optional numeric vector controlling which views supply
#'   the five nonlinear Friedman features when \code{outcome.model} is
#'   \code{"friedman"} or \code{"friedman2"}: a single value, one value per view
#'   or a vector named by view, normalized to sum to one. Each of the five
#'   features comes from a view chosen with probability equal to its share
#'   (regardless of how many candidate features the view has), then uniformly
#'   among that view's spiked features. If NULL, \code{signal.share} is used
#'   when available; otherwise the five features are sampled uniformly from all
#'   spiked features.
#' @param outcome.model Outcome model: "linear", "friedman", "friedman2",
#'   "shared_latent", or "interaction". Friedman models apply the Friedman formula
#'   after rank-transforming the selected simulated features to the unit
#'   interval; "friedman2" mixes it with the linear predictor using
#'   \code{nonlinear.weight}. "shared_latent" uses the simulated shared latent
#'   factors. "interaction" mixes the linear predictor with products of feature
#'   pairs taken from different views, using \code{interaction.weight}.
#' @param shared.latent.function Response function used when
#'   \code{outcome.model = "shared_latent"}. \code{"linear"} combines selected shared
#'   latent factors as a standardized equal-weight score, preserving the
#'   previous
#'   behavior. \code{"friedman"} applies the Friedman nonlinear function to the
#'   first five selected shared latent factors after rank-transforming them to
#'   the unit interval. \code{"friedman2"} mixes the linear shared-latent score
#'   with the Friedman nonlinear component using \code{nonlinear.weight}.
#' @param shared.latent.n.factors Number of shared latent factors to use for
#'   \code{outcome.model = "shared_latent"}. If \code{NULL}, all available shared
#'   factors are used. Friedman shared-latent functions require at least five
#'   selected shared factors.
#' @param shared.latent.weight Share of signal variance from the shared latent
#'   factors, between 0 and 1, for \code{outcome.model = "shared_latent"}. The
#'   rest comes from view-specific latent signal: in each view, the part of the
#'   view's latent factor scores that is not shared with the other views
#'   (independent across views), summed with weights proportional to each
#'   factor's loading strength. Views contribute equally, or by
#'   \code{signal.share} when it is given. Both parts are standardized and
#'   combined as \code{sqrt(w) * shared + sqrt(1 - w) * view_specific}. The
#'   default \code{1} is a purely shared outcome (every view carries the same
#'   signal); smaller values make the views complementary, since each view
#'   then carries outcome information the others do not. The view-specific
#'   part is stored in \code{sample_metadata$Xbeta_specific} and the per-view
#'   scores in \code{truth$specific_scores}.
#' @param nonlinear.weight Share of signal variance from the nonlinear Friedman
#'   part, between 0 and 1, for \code{outcome.model = "friedman2"} and
#'   \code{shared.latent.function = "friedman2"}. Both parts are standardized and
#'   combined as \code{sqrt(1 - w) * linear + sqrt(w) * nonlinear}.
#' @param interaction.weight Share of signal variance from cross-view
#'   interactions, between 0 and 1, for \code{outcome.model = "interaction"}.
#'   \code{1} gives a pure interaction outcome with no main effects.
#' @param n.interactions Number of cross-view feature pairs for
#'   \code{outcome.model = "interaction"}. Defaults to the number of view pairs.
#'   Pairs cycle over view pairs; each feature is drawn from that view's
#'   nonzero-coefficient features when available.
#' @param outcome.type Outcome type: \code{"continuous"}, \code{"binary"}, or
#'   \code{"survival"}. Continuous and binary outcomes are stored as
#'   \code{Y}; survival outcomes are stored as \code{time} and \code{status}.
#' @param binary.prevalence Target proportion of cases for binary outcomes.
#'   The logistic intercept is solved in each replicate so that the expected
#'   prevalence equals this value.
#' @param baseline.hazard Multiplicative scale on the baseline hazard for survival
#'   outcomes.
#' @param baseline.shape Weibull shape \eqn{k} of the baseline hazard for
#'   survival outcomes (default 1, a constant hazard: exponential event times).
#'   Event times then have hazard
#'   \code{baseline.hazard * k * t^(k - 1) * exp(eta)}, still proportional
#'   hazards: \code{k > 1} gives a hazard that rises over follow-up,
#'   \code{k < 1} one that falls. \code{snr} keeps its meaning on the log-time
#'   scale, since signal and error variance of \code{log T} both scale by
#'   \eqn{1/k^2}.
#' @param keep.true.times Logical; if \code{TRUE}, survival outcomes also store
#'   the latent event time and censoring time of every sample in
#'   \code{sample_metadata$event_time} and \code{sample_metadata$cens_time}
#'   (\code{time} is their minimum). Default \code{FALSE}.
#' @param cens.rate Optional target censoring proportion for survival outcomes.
#'   If supplied, the censoring times of the chosen \code{cens.type} (or the
#'   supplied \code{cens.times}) are multiplied by a common factor \code{c},
#'   solved so that this proportion of subjects is censored; \code{0} means no
#'   censoring. If \code{NULL}, uniform censoring times are uniform on
#'   \code{(cens.lower, cens.upper)} and supplied \code{cens.times} are used
#'   as they are; the other censoring types then need \code{cens.rate}.
#' @param cens.lower Lower bound for uniform censoring time when
#'   \code{cens.rate} is \code{NULL}.
#' @param cens.upper Upper bound for uniform censoring time when
#'   \code{cens.rate} is \code{NULL}.
#' @param cens.type Censoring mechanism for survival outcomes.
#'   \code{"uniform"} (default): censoring times uniform on \code{(0, c)}.
#'   \code{"exponential"}: exponential censoring times (constant loss to
#'   follow-up), mean \code{c}. \code{"administrative"}: staggered entry and a
#'   fixed study end; subjects enter uniformly over the first
#'   \code{cens.accrual} fraction of the study, the study ends at time
#'   \code{c}, and subjects still at risk are censored then, so the censoring
#'   time is \code{c * (1 - cens.accrual * U)}; with \code{cens.dropout > 0}
#'   subjects can also drop out earlier. In every case \code{c} is solved for
#'   \code{cens.rate}.
#' @param cens.times Optional user-supplied censoring times, used instead of
#'   \code{cens.type}: a single time for every subject (e.g. a fixed end of
#'   follow-up), one positive time per simulated sample (in sample order, before
#'   the train/test split; \code{Inf} for never censored), or a function of
#'   \code{n} returning \code{n} times. With \code{cens.rate} the times are
#'   rescaled by a common factor to reach the target, so they give only the
#'   shape of the censoring distribution.
#' @param cens.dependence Informative censoring strength (default 0,
#'   independent censoring). Each censoring time is multiplied by
#'   \code{exp(-cens.dependence * z)}, where \code{z} is the standardized
#'   linear predictor, so with a positive value higher-risk subjects are
#'   censored earlier (and with a negative value, later). This violates the
#'   independent-censoring assumption of most survival methods. Censoring is
#'   still solved for \code{cens.rate}.
#' @param cens.accrual Fraction of the study period over which subjects enter
#'   when \code{cens.type = "administrative"}, between 0 and 1 (default 0.5).
#'   \code{0} gives a single fixed censoring time.
#' @param cens.dropout Probability that a subject drops out (is lost to
#'   follow-up) before the study end when \code{cens.type = "administrative"},
#'   between 0 and 1 (default 0). Dropout times are exponential.
#' @param nrep Number of independent simulation replicates (default 1).
#' @param seed Random seed. The global random-number state is restored when the
#'   function returns. Use \code{NULL} to use the current random state.
#' @param return.template Logical; if \code{TRUE}, the fitted template is
#'   returned in \code{$template}. The template contains the source data (see
#'   \code{fit_template()}), so the default \code{FALSE} returns only its
#'   provenance in \code{$template_provenance}, making the result safe to share.
#' @param back.transform Logical; if \code{TRUE}, returned features of views
#'   fitted with a \code{view.transform} are mapped back to the original
#'   measurement scale (see \code{sim_from_template()}). Outcomes are always
#'   generated from the fitted (transformed) scale.
#' @param ... Passed to \code{fit_template()} when \code{template} is a raw list
#'
#' @details The ground truth for which features drive the outcome is in
#'   \code{feature_metadata$outcome_feature}. \code{beta} is nonzero only for
#'   features that enter the outcome through the linear predictor; in
#'   \code{"friedman"}, \code{"shared_latent"} and pure-interaction outcomes it
#'   is
#'   zero for every feature. \code{friedman_term} and \code{interaction_term}
#'   mark the nonlinear and interaction features, and
#'   \code{truth$interaction_pairs} lists the interacting pairs.
#'
#'   Linear predictor. Within each replicate, every simulated feature is
#'   standardized to mean 0 and SD 1 across the simulated samples before it is
#'   multiplied by its coefficient, so \code{beta} is an effect per standard
#'   deviation and views on very different measurement scales contribute
#'   comparably. The linear predictor is \code{Xs \%*\% beta}, with \code{Xs}
#'   the standardized feature matrix (with \code{signal.share}, each view's
#'   part is rescaled first). For continuous outcomes the noise variance is
#'   then set to \code{Var(eta) / snr}; for binary and survival outcomes the
#'   predictor is standardized and rescaled to match \code{snr}, and the
#'   final predictor is stored in \code{sample_metadata$Xbeta}. To reproduce the predictor from the returned data,
#'   standardize the features of the full replicate (training and testing
#'   samples together) on the scale they were simulated on.
#'
#' @return List with trainDat, testDat, simulation parameters, the template
#'   provenance (or the template itself if \code{return.template = TRUE}),
#'   realized canonical correlations (\code{cancor}) and per-replicate
#'   empirical agreement between view-specific linear predictors
#'   (\code{agreement.emp}).
#' @examples
#' data(stelzer_mae)
#' tmpl <- fit_template(stelzer_mae, shared.rank = 2)
#' sim <- gen_simmba(
#'     n = 60, template = tmpl, p = 10, nrep = 2,
#'     outcome.model = "interaction", outcome.type = "binary", seed = 1
#' )
#' head(sim$trainDat$Rep_1$sample_metadata)
#' sim$trainDat$Rep_1$truth$interaction_pairs
#' @export
gen_simmba <- function(n,
    template,
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

    n <- .check_count(n, "n", min = 2L)
    if (!is.null(seed)) withr::local_seed(seed)
    nrep <- .check_count(nrep, "nrep", min = 1L)
    .check_prevalence(train.fraction, "train.fraction")
    n_train <- round(n * train.fraction)
    if (n_train < 1L || n_train > n - 1L) {
        stop("train.fraction = ", train.fraction, " with n = ", n,
            " leaves an empty training or test set.", call. = FALSE)
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
    if (length(cfg$cens.rate) > 1L) {
        stop("cens.rate must be NULL or a single number in [0, 1).", call. = FALSE)
    }

    trainDat <- testDat <- agreement.emp <- vector("list", nrep)
    names(trainDat) <- names(testDat) <- names(agreement.emp) <- paste("Rep",
        seq_len(nrep), sep = "_")

    for (k in seq_len(nrep)) {

        ## 1-5. Null data, spike-in coefficients, predictor and outcome (one
        ## study, see .simulate_outcome_replicate())
        res <- .simulate_outcome_replicate(tmpl, n, cfg)
        pcl <- res$pcl
        agreement.emp[[k]] <- res$agreement

        ## 6. Train / test split
        tr.row <- sample.int(n, size = round(n * train.fraction),
            replace = FALSE)
        trainDat[[k]] <- .subset_pcl_samples(pcl, tr.row)
        testDat[[k]]  <- .subset_pcl_samples(pcl, -tr.row)
    }

    out <- list(trainDat = trainDat, testDat = testDat,
        template = if (isTRUE(return.template)) tmpl else NULL,
        template_provenance = tmpl$provenance, p = p, crosstalk = crosstalk,
        view.models = view.models,
        cancor = .realized_cancor(tmpl, crosstalk),
        agreement.emp = agreement.emp,
        snr = snr, train.fraction = train.fraction,
        signal.prob = cfg$signal.prob, signal.neg.prob = cfg$signal.neg.prob,
        effect.meanlog = cfg$effect.meanlog,
        effect.sdlog = cfg$effect.sdlog, signal.align = signal.align,
        signal.share = cfg$signal.share, friedman.share = cfg$friedman.share,
        shared.latent.function = cfg$shared.latent.function,
        shared.latent.n.factors = cfg$shared.latent.n.factors,
        nonlinear.weight = nonlinear.weight, back.transform = back.transform,
        interaction.weight = interaction.weight,
        n.interactions = n.interactions,
        nrep = nrep, seed = seed, outcome.model = cfg$outcome.model,
        outcome.type = cfg$outcome.type, binary.prevalence = binary.prevalence,
        baseline.hazard = baseline.hazard, cens.rate = cens.rate,
        cens.lower = cens.lower, cens.upper = cens.upper)
    ## New settings are listed only when they differ from their defaults, so
    ## results made with the defaults keep their earlier structure
    .add_nondefault_settings(out, cfg)
}

## Settings added after the first release are recorded only when they are
## not at their defaults (they are also in sample_metadata and truth).
.add_nondefault_settings <- function(out, cfg) {
    if (cfg$shared.latent.weight < 1) out$shared.latent.weight <- cfg$shared.latent.weight
    if (!identical(cfg$cens.type, "uniform")) out$cens.type <- cfg$cens.type
    if (!is.null(cfg$cens.times)) out$cens.times <- cfg$cens.times
    if (cfg$cens.dependence != 0) out$cens.dependence <- cfg$cens.dependence
    if (identical(cfg$cens.type, "administrative")) {
        out$cens.accrual <- cfg$cens.accrual
        out$cens.dropout <- cfg$cens.dropout
    }
    if (cfg$baseline.shape != 1) out$baseline.shape <- cfg$baseline.shape
    if (isTRUE(cfg$keep.true.times)) out$keep.true.times <- TRUE
    out
}

## Argument checks and template set-up shared by gen_simmba() and
## gen_multistudy(). Returns the fitted template and the validated settings
## used by .simulate_outcome_replicate().
.simmba_setup <- function(template, fit_args, p, crosstalk, feature.select,
    feature.jitter, view.models, snr, signal.prob, signal.neg.prob,
    effect.meanlog, effect.sdlog, signal.align, signal.share, friedman.share,
    outcome.model, shared.latent.function, shared.latent.n.factors,
    nonlinear.weight, interaction.weight, n.interactions, outcome.type,
    binary.prevalence, baseline.hazard, cens.rate, cens.lower, cens.upper,
    back.transform, shared.latent.weight = 1,
    cens.type = "uniform", cens.times = NULL, cens.dependence = 0,
    cens.accrual = 0.5, cens.dropout = 0, baseline.shape = 1,
    keep.true.times = FALSE) {
    outcome.model  <- match.arg(outcome.model, c("linear", "friedman",
        "friedman2", "shared_latent", "interaction"))
    ## Internal codes for the outcome models
    ygen.mode <- c(linear = "LM", friedman = "Friedman", friedman2 = "Friedman2",
        shared_latent = "SharedLatent", interaction = "Interaction")[[outcome.model]]
    shared.latent.function <- match.arg(shared.latent.function,
        c("linear", "friedman", "friedman2"))
    outcome.type   <- match.arg(outcome.type, c("continuous", "binary", "survival"))
    feature.select <- match.arg(feature.select, c("random", "variance"))
    if (!is.null(shared.latent.n.factors)) {
        shared.latent.n.factors <- .check_count(shared.latent.n.factors,
            "shared.latent.n.factors", min = 1L)
    }
    .check_fraction(nonlinear.weight, "nonlinear.weight")
    .check_fraction(interaction.weight, "interaction.weight")
    if (!is.null(n.interactions)) .check_count(n.interactions, "n.interactions", min = 1L)
    if (!is.numeric(snr) || length(snr) != 1L || !is.finite(snr) || snr <= 0) {
        stop("snr must be a single positive number.", call. = FALSE)
    }
    .check_prevalence(binary.prevalence, "binary.prevalence")
    if (!is.null(cens.rate) && (!is.numeric(cens.rate) || length(cens.rate) < 1L ||
        any(!is.finite(cens.rate)) || any(cens.rate < 0) || any(cens.rate >= 1))) {
        stop("cens.rate must be NULL or a number in [0, 1).", call. = FALSE)
    }
    .check_fraction(shared.latent.weight, "shared.latent.weight")
    if (!is.numeric(baseline.shape) || length(baseline.shape) != 1L ||
        !is.finite(baseline.shape) || baseline.shape <= 0) {
        stop("baseline.shape must be a single positive number.", call. = FALSE)
    }
    if (!is.logical(keep.true.times) || length(keep.true.times) != 1L ||
        is.na(keep.true.times)) {
        stop("keep.true.times must be TRUE or FALSE.", call. = FALSE)
    }
    cens.type <- match.arg(cens.type, c("uniform", "exponential", "administrative"))
    .check_censoring_times(cens.times)
    if (!is.numeric(cens.dependence) || length(cens.dependence) != 1L ||
        !is.finite(cens.dependence)) {
        stop("cens.dependence must be a single finite number.", call. = FALSE)
    }
    .check_fraction(cens.accrual, "cens.accrual")
    if (!is.numeric(cens.dropout) || length(cens.dropout) != 1L ||
        !is.finite(cens.dropout) || cens.dropout < 0 || cens.dropout >= 1) {
        stop("cens.dropout must be a number in [0, 1).", call. = FALSE)
    }
    if (identical(outcome.type, "survival") && is.null(cens.rate) &&
        is.null(cens.times) && !identical(cens.type, "uniform")) {
        stop("cens.type = \"", cens.type, "\" needs cens.rate (the target censoring ",
            "proportion), which sets the time scale of the censoring times.", call. = FALSE)
    }

    ## 0. Fit the template once (skipped if a fitted template is supplied)
    tmpl  <- .resolve_template(template, fit_args)
    views <- tmpl$views
    M     <- length(views)
    signal.prob     <- .per_view_values(signal.prob, views, "signal.prob")
    signal.neg.prob <- .per_view_values(signal.neg.prob, views, "signal.neg.prob")
    effect.meanlog   <- .per_view_values(effect.meanlog, views, "effect.meanlog")
    effect.sdlog <- .per_view_values(effect.sdlog, views, "effect.sdlog")
    .check_unit_interval(signal.prob, "signal.prob")
    .check_unit_interval(signal.neg.prob, "signal.neg.prob")
    if (any(!is.finite(effect.meanlog))) stop("effect.meanlog must be finite.", call. = FALSE)
    if (any(!is.finite(effect.sdlog)) || any(effect.sdlog < 0)) {
        stop("effect.sdlog must be finite and non-negative.", call. = FALSE)
    }
    if (!is.null(signal.align)) .check_unit_interval(signal.align, "signal.align", single = TRUE)
    gamma_named <- .per_view_values(crosstalk, views, "crosstalk")
    signal.share <- .normalize_signal_share(signal.share, views, "signal.share")
    friedman.share <- .normalize_signal_share(friedman.share, views,
        "friedman.share")
    if (is.null(friedman.share)) {
        friedman.share <- signal.share
    }
    if (identical(ygen.mode, "Interaction") && M < 2L) {
        stop("outcome.model = \"interaction\" needs at least two views.",
            call. = FALSE)
    }
    if (identical(ygen.mode, "SharedLatent")) {
        if (tmpl$shared.rank == 0L) {
            warning("The template has no shared latent factors; SharedLatent outcomes fall back ",
                "to the average of the view-specific linear predictors.", call. = FALSE)
        } else if (all(gamma_named == 0) && shared.latent.weight > 0) {
            warning("outcome.model = \"shared_latent\" with crosstalk = 0: shared factors do not enter ",
                "any simulated feature, so the outcome is unrelated to the features.", call. = FALSE)
        }
    }
    list(template = tmpl, cfg = list(
        p = p, crosstalk = crosstalk, feature.select = feature.select,
        feature.jitter = feature.jitter, view.models = view.models,
        back.transform = back.transform, snr = snr,
        signal.prob = signal.prob, signal.neg.prob = signal.neg.prob,
        effect.meanlog = effect.meanlog, effect.sdlog = effect.sdlog,
        signal.align = signal.align, signal.share = signal.share,
        friedman.share = friedman.share, outcome.model = outcome.model,
        ygen.mode = ygen.mode, shared.latent.function = shared.latent.function,
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
        keep.true.times = keep.true.times))
}


## ---- outcome helpers ----------------------------------------------------

.standardize_columns <- function(x) {
    x <- as.matrix(x)
    ## Two-pass (centre, then scale) to avoid cancellation for features with
    ## large means
    out <- sweep(x, 2, colMeans(x), "-")
    s <- sqrt(colSums(out^2) / max(nrow(x) - 1, 1))
    ok <- is.finite(s) & s > 0
    out[, ok] <- sweep(out[, ok, drop = FALSE], 2, s[ok], "/")
    out[, !ok] <- 0
    out
}

## Censoring times. The original rule (cens.type "uniform", no supplied times,
## independent censoring) is kept as it was: uniform(cens.lower, cens.upper),
## or uniform(0, c) with c solved so that the requested proportion of subjects
## is censored. Otherwise a censoring time is c * V: V is drawn from the
## chosen mechanism (or supplied), optionally multiplied by
## exp(-dependence * z) for informative censoring, and c is solved for
## cens.rate (c = 1 without a target).
.draw_censoring <- function(T0, cens.rate = NULL, cens.lower = 1,
    cens.upper = 3, cens.type = "uniform", cens.times = NULL,
    cens.dependence = 0, eta = NULL, cens.accrual = 0.5, cens.dropout = 0) {
    n <- length(T0)
    if (identical(cens.type, "uniform") && is.null(cens.times) &&
        cens.dependence == 0) {
        if (is.null(cens.rate)) {
            return(stats::runif(n, cens.lower, cens.upper))
        }
        if (cens.rate == 0) {
            return(rep(Inf, n))
        }
        u <- stats::runif(n)
        f <- function(log_c) mean(u * exp(log_c) < T0) - cens.rate
        lo <- log(min(T0)) - 10
        hi <- log(max(T0)) + 10
        log_c <- stats::uniroot(f, interval = c(lo, hi), tol = 1e-10)$root
        return(u * exp(log_c))
    }
    if (!is.null(cens.rate) && cens.rate == 0) {
        return(rep(Inf, n))
    }
    V <- if (!is.null(cens.times)) {
        as.numeric(cens.times)
    } else switch(cens.type,
        uniform = if (is.null(cens.rate)) stats::runif(n, cens.lower, cens.upper) else
            stats::runif(n),
        exponential = stats::rexp(n),
        administrative = {
            V <- 1 - cens.accrual * stats::runif(n)
            if (cens.dropout > 0) {
                ## dropout before the study end (time 1) with probability cens.dropout
                V <- pmin(V, stats::rexp(n) / -log(1 - cens.dropout))
            }
            V
        })
    if (cens.dependence != 0) {
        if (is.null(eta)) stop("Informative censoring needs the linear predictor.", call. = FALSE)
        V <- V * exp(-cens.dependence * .standardize_predictor(eta))
    }
    if (is.null(cens.rate)) return(V)
    .solve_censoring_scale(T0, V, cens.rate)
}

## Smallest common factor c with mean(c * V < T0) equal to the target (to
## within one subject), searched on the log scale.
.solve_censoring_scale <- function(T0, V, cens.rate) {
    fin <- is.finite(V)
    max_rate <- mean(fin)
    if (cens.rate > max_rate) {
        stop("cens.rate = ", cens.rate, " cannot be reached: only ",
            round(100 * max_rate, 1), "% of the supplied censoring times are finite.",
            call. = FALSE)
    }
    f <- function(log_c) mean(V * exp(log_c) < T0) - cens.rate
    lo <- log(min(T0)) - log(max(V[fin])) - 10
    hi <- log(max(T0)) - log(min(V[fin])) + 10
    log_c <- stats::uniroot(f, interval = c(lo, hi), tol = 1e-10)$root
    V * exp(log_c)
}

## User censoring times: NULL, a function, or positive numbers (Inf allowed).
## The length is checked against the number of samples when they are drawn.
.check_censoring_times <- function(cens.times) {
    if (is.null(cens.times) || is.function(cens.times)) return(invisible(TRUE))
    if (is.list(cens.times)) {
        for (x in cens.times) .check_censoring_times(x)
        return(invisible(TRUE))
    }
    if (!is.numeric(cens.times) || length(cens.times) == 0L || anyNA(cens.times) ||
        any(cens.times <= 0)) {
        stop("cens.times must be NULL, a function of n, or positive numbers ",
            "(Inf for never censored).", call. = FALSE)
    }
    invisible(TRUE)
}

## Censoring times supplied for the rows of one study: a function of n, one
## time for all, one per study sample, or one per sample of all studies
.resolve_censoring_times <- function(cens.times, rows, n_total, study) {
    if (is.null(cens.times)) return(NULL)
    if (is.list(cens.times) && !is.function(cens.times)) cens.times <- cens.times[[study]]
    n <- length(rows)
    x <- if (is.function(cens.times)) cens.times(n) else cens.times
    .check_censoring_times(x)
    if (is.function(x) || !is.numeric(x)) {
        stop("A cens.times function must return numbers.", call. = FALSE)
    }
    if (length(x) == 1L) return(rep(x, n))
    if (length(x) == n) return(as.numeric(x))
    if (length(x) == n_total) return(as.numeric(x[rows]))
    stop("cens.times has ", length(x), " values; give one time, one per sample (",
        n_total, ") or a function of n.", call. = FALSE)
}

## ---- spike-in helpers ---------------------------------------------------

## Same distribution as splatter:::getLNormFactors() restricted to selected
## features, without a ::: call
.lnorm_factors <- function(n, neg.prob, fac.loc, fac.scale) {
    if (n == 0) return(numeric(0))
    dir  <- (-1)^stats::rbinom(n, 1, neg.prob)
    facs <- stats::rlnorm(n, fac.loc, fac.scale)
    dir[facs < 1] <- -dir[facs < 1]
    facs^dir
}

## Log fold changes for one view. de.align = NULL reproduces the original
## splatter spike-in. Otherwise a fraction de.align of DE features is drawn with
## probability proportional to shared_loading^2 and signed by shared_loading.
.draw_beta_view <- function(p_m, shared_loading, de.prob, neg.prob, fac.loc,
    fac.scale, de.align = NULL) {
    beta <- numeric(p_m)
    if (p_m == 0L || de.prob <= 0) return(beta)
    ## At least one signal feature per view with de.prob > 0, so a replicate
    ## never ends up with no outcome signal by chance
    n_de <- max(1L, stats::rbinom(1, p_m, de.prob))
    n_al <- if (is.null(de.align)) 0L else as.integer(round(de.align * n_de))
    pos_al <- if (n_al > 0) sample.int(p_m, n_al,
        prob = shared_loading^2 + 1e-12) else integer(0)
    rest   <- setdiff(seq_len(p_m), pos_al)
    pos_un <- rest[sample.int(length(rest), n_de - n_al)]
    pos    <- c(pos_al, pos_un)
    beta[pos] <- log2(.lnorm_factors(length(pos), neg.prob, fac.loc, fac.scale))
    if (n_al > 0) {
        aligned_sign <- sign(shared_loading[pos_al])
        fallback_sign <- sign(beta[pos_al])
        fallback_sign[fallback_sign == 0 | !is.finite(fallback_sign)] <- 1
        aligned_sign[aligned_sign == 0 | !is.finite(aligned_sign)] <- fallback_sign[aligned_sign == 0 | !is.finite(aligned_sign)]
        beta[pos_al] <- abs(beta[pos_al]) * aligned_sign
    }
    beta
}

.normalize_signal_share <- function(signal.share, views,
    arg_name = "signal.share") {
    if (is.null(signal.share)) {
        return(NULL)
    }
    if (!is.numeric(signal.share) || length(signal.share) == 0L) {
        stop(arg_name, " must be a numeric vector.", call. = FALSE)
    }
    if (!is.null(names(signal.share))) {
        missing <- setdiff(views, names(signal.share))
        if (length(missing) > 0L) {
            stop(arg_name, " is missing view(s): ", paste(missing,
                collapse = ", "), call. = FALSE)
        }
        signal.share <- signal.share[views]
    } else {
        signal.share <- .per_view_values(signal.share, views, arg_name)
    }
    if (any(!is.finite(signal.share)) || any(signal.share < 0)) {
        stop(arg_name, " values must be finite and non-negative.",
            call. = FALSE)
    }
    if (sum(signal.share) <= 0) {
        stop("At least one ", arg_name, " value must be positive.",
            call. = FALSE)
    }
    signal.share / sum(signal.share)
}

## Sample `size` elements of `x` (safe when length(x) == 1, unlike sample()).
.sample_values <- function(x, size, prob = NULL) {
    x[sample.int(length(x), size, prob = prob)]
}

.sample_friedman_features <- function(nonzero_index, feature_type, views,
    friedman.share) {
    if (is.null(friedman.share)) {
        return(.sample_values(nonzero_index, 5))
    }
    ## View-level shares: each of the five Friedman slots picks a view with
    ## probability proportional to its share (among views that still have
    ## unused candidates), then a feature uniformly within that view.
    pool <- split(nonzero_index, factor(as.character(feature_type[nonzero_index]),
        levels = views))
    share <- friedman.share[views]
    share[!is.finite(share)] <- 0
    out <- integer(0)
    for (slot in seq_len(5L)) {
        avail <- lengths(pool) > 0L & share > 0
        if (!any(avail)) break
        v <- if (sum(avail) == 1L) which(avail) else
            sample(which(avail), 1L, prob = share[avail])
        pick <- .sample_values(pool[[v]], 1L)
        out <- c(out, pick)
        pool[[v]] <- setdiff(pool[[v]], pick)
    }
    if (length(out) < 5L) {
        ## Not enough candidates in the views with positive share
        fallback <- setdiff(nonzero_index, out)
        out <- c(out, .sample_values(fallback, 5L - length(out)))
    }
    out
}

## Standardize two signals and mix them so that w is the share of signal
## variance from `b` (exact when the two parts are uncorrelated).
.mix_signals <- function(a, b, w) {
    if (w <= 0) return(.standardize_predictor(a))
    if (w >= 1) return(.standardize_predictor(b))
    sqrt(1 - w) * .standardize_predictor(a) + sqrt(w) * .standardize_predictor(b)
}

## Cross-view feature pairs for interaction outcomes. Pairs cycle over the
## view pairs; within a view, features with nonzero beta are preferred.
.sample_interaction_pairs <- function(beta, feature_type, views,
    n.interactions = NULL) {
    view_pairs <- utils::combn(views, 2)
    if (is.null(n.interactions)) n.interactions <- ncol(view_pairs)
    n.interactions <- as.integer(n.interactions)
    pool <- lapply(views, function(v) {
        idx <- which(feature_type == v)
        nz <- idx[beta[idx] != 0]
        if (length(nz) > 0L) nz else idx
    })
    names(pool) <- views
    used <- integer(0)
    i <- j <- integer(n.interactions)
    for (t in seq_len(n.interactions)) {
        vp <- view_pairs[, (t - 1L) %% ncol(view_pairs) + 1L]
        pick <- function(v) {
            cand <- setdiff(pool[[v]], used)
            if (length(cand) == 0L) cand <- pool[[v]]
            .sample_values(cand, 1)
        }
        i[t] <- pick(vp[[1]])
        j[t] <- pick(vp[[2]])
        used <- c(used, i[t], j[t])
    }
    list(i = i, j = j)
}

.standardize_predictor <- function(x) {
    x <- as.numeric(x)
    s <- stats::sd(x)
    if (!is.finite(s) || s <= 0) {
        return(x - mean(x, na.rm = TRUE))
    }
    as.numeric(scale(x))
}

.rank_to_unit_interval <- function(x) {
    x <- as.matrix(x)
    out <- apply(x, 2, function(col) {
        ok <- is.finite(col)
        z <- rep(0.5, length(col))
        if (sum(ok) > 1L) {
            z[ok] <- rank(col[ok], ties.method = "average") / (sum(ok) + 1)
        }
        z
    })
    out <- as.matrix(out)
    dimnames(out) <- dimnames(x)
    out
}

.friedman_response <- function(x) {
    x <- as.matrix(x)
    if (ncol(x) < 5L) {
        stop("Friedman response requires at least five predictors.",
            call. = FALSE)
    }
    10 * sin(pi * x[, 1] * x[, 2]) + 20 * (x[, 3] - 0.5)^2 + 10 * x[,
        4] + 5 * x[, 5]
}

.shared_latent_predictor <- function(pcl, view_pred, fun = c("linear",
    "friedman", "friedman2"),
    n_factors = NULL, nonlinear.weight = 0.5) {
    fun <- match.arg(fun)
    shared <- pcl$truth$shared
    if (is.matrix(shared) && ncol(shared) > 0L) {
        if (!is.null(n_factors)) {
            if (n_factors > ncol(shared)) {
                stop(
                    "shared.latent.n.factors = ", n_factors,
                    " but only ", ncol(shared), " shared latent factor(s) are available.",
                    call. = FALSE
                )
            }
            shared <- shared[, seq_len(n_factors), drop = FALSE]
        }
        shared_scaled <- scale(shared)
        shared_scaled[!is.finite(shared_scaled)] <- 0
        eta_linear <- .standardize_predictor(rowSums(shared_scaled) / sqrt(ncol(shared_scaled)))

        if (identical(fun, "linear")) {
            return(eta_linear)
        }
        if (ncol(shared_scaled) < 5L) {
            stop(
                "shared.latent.function = '", fun,
                "' requires at least five selected shared latent factors; only ",
                ncol(shared_scaled), " are available. Increase shared.rank in fit_template() ",
                "or use shared.latent.function = 'linear'.",
                call. = FALSE
            )
        }
        eta_friedman <- .standardize_predictor(.friedman_response(.rank_to_unit_interval(shared_scaled[,
            seq_len(5), drop = FALSE])))
        if (identical(fun, "friedman")) {
            return(eta_friedman)
        }
        return(.mix_signals(eta_linear, eta_friedman, nonlinear.weight))
    }

    if (!identical(fun, "linear")) {
        stop("Nonlinear shared latent outcomes require simulated shared latent factors.",
            call. = FALSE)
    }
    view_pred <- as.matrix(view_pred)
    if (ncol(view_pred) == 0L) {
        stop("Cannot construct shared latent predictor without shared scores or view predictors.",
            call. = FALSE)
    }
    out <- .standardize_predictor(rowMeans(scale(view_pred), na.rm = TRUE))
    attr(out, "fallback") <- TRUE
    out
}

.rescale_signal_by_view <- function(beta, view_pred, feature_type, views,
    signal.share) {
    beta <- as.numeric(beta)
    view_pred <- as.matrix(view_pred)
    colnames(view_pred) <- views

    for (view in views) {
        target <- signal.share[[view]]
        sel <- feature_type == view
        if (target == 0) {
            beta[sel] <- 0
            view_pred[, view] <- 0
            next
        }
        current_var <- stats::var(view_pred[, view])
        if (!is.finite(current_var) || current_var <= 0) {
            warning("Cannot rescale signal for view '", view,
                "' because its linear predictor variance is zero.", call. = FALSE)
            next
        }
        scale_factor <- sqrt(target / current_var)
        beta[sel] <- beta[sel] * scale_factor
        view_pred[, view] <- view_pred[, view] * scale_factor
    }

    list(beta = beta, view_pred = view_pred)
}
