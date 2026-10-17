#' Generate simulated multi-omics data
#'
#' Generates null multiview data from a real template (see \code{fit_template()}
#' and \code{sim_from_template()}), spikes in signal, generates the outcome, and
#' splits the data into training and testing sets.
#'
#' @param nsample Number of samples to simulate in each replicate before the
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
#' @param p.train Proportion of samples assigned to the training split.
#' @param de.prob Probability that a feature receives a nonzero spike-in
#'   coefficient; scalar, one value per view, or named by view. Every view with
#'   \code{de.prob > 0} gets at least one nonzero coefficient; use 0 for a view
#'   that should carry no linear signal.
#' @param de.downProb Probability that a nonzero spike-in coefficient is
#'   negative; scalar or one value per view.
#' @param de.facLoc Location parameter for the log-normal effect-size factor;
#'   scalar or one value per view.
#' @param de.facScale Scale parameter for the log-normal effect-size factor;
#'   scalar or one value per view.
#' @param de.align NULL keeps the original spike-in (DE positions and signs
#'   drawn
#'   independently per view). A value between 0 and 1 places that fraction of DE
#'   features on features loading on the shared cross-view direction (sampled
#'   with
#'   probability proportional to the squared loading) with signs following that
#'   direction, so view-specific signals agree to the extent the views
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
#'   the five nonlinear Friedman features when \code{ygen.mode} is
#'   \code{"Friedman"} or \code{"Friedman2"}: a single value, one value per view
#'   or a vector named by view, normalized to sum to one. Each of the five
#'   features comes from a view chosen with probability equal to its share
#'   (regardless of how many candidate features the view has), then uniformly
#'   among that view's spiked features. If NULL, \code{signal.share} is used
#'   when available; otherwise the five features are sampled uniformly from all
#'   spiked features.
#' @param ygen.mode Y generation mode: "LM", "Friedman", "Friedman2",
#'   "SharedLatent", or "Interaction". Friedman modes apply the Friedman formula
#'   after rank-transforming the selected simulated features to the unit
#'   interval; "Friedman2" mixes it with the linear predictor using
#'   \code{nonlinear.weight}. "SharedLatent" uses the simulated shared latent
#'   factors. "Interaction" mixes the linear predictor with products of feature
#'   pairs taken from different views, using \code{interaction.weight}.
#' @param shared.latent.fun Response function used when
#'   \code{ygen.mode = "SharedLatent"}. \code{"linear"} combines selected shared
#'   latent factors as a standardized equal-weight score, preserving the
#'   previous
#'   behavior. \code{"friedman"} applies the Friedman nonlinear function to the
#'   first five selected shared latent factors after rank-transforming them to
#'   the unit interval. \code{"friedman2"} mixes the linear shared-latent score
#'   with the Friedman nonlinear component using \code{nonlinear.weight}.
#' @param shared.latent.n_factors Number of shared latent factors to use for
#'   \code{ygen.mode = "SharedLatent"}. If \code{NULL}, all available shared
#'   factors are used. Friedman shared-latent functions require at least five
#'   selected shared factors.
#' @param nonlinear.weight Share of signal variance from the nonlinear Friedman
#'   part, between 0 and 1, for \code{ygen.mode = "Friedman2"} and
#'   \code{shared.latent.fun = "friedman2"}. Both parts are standardized and
#'   combined as \code{sqrt(1 - w) * linear + sqrt(w) * nonlinear}.
#' @param interaction.weight Share of signal variance from cross-view
#'   interactions, between 0 and 1, for \code{ygen.mode = "Interaction"}.
#'   \code{1} gives a pure interaction outcome with no main effects.
#' @param n.interactions Number of cross-view feature pairs for
#'   \code{ygen.mode = "Interaction"}. Defaults to the number of view pairs.
#'   Pairs cycle over view pairs; each feature is drawn from that view's
#'   nonzero-coefficient features when available.
#' @param outcome.type Outcome type: \code{"continuous"}, \code{"binary"}, or
#'   \code{"survival"}. Continuous and binary outcomes are stored as
#'   \code{Y}; survival outcomes are stored as \code{time} and \code{status}.
#' @param binary.prevalence Target proportion of cases for binary outcomes.
#'   The logistic intercept is solved in each replicate so that the expected
#'   prevalence equals this value.
#' @param surv.hscale Multiplicative scale on the baseline hazard for survival
#'   outcomes.
#' @param cens.rate Optional target censoring proportion for survival outcomes.
#'   If supplied, censoring times are uniform on \code{(0, c)} with \code{c}
#'   solved so that this proportion of subjects is censored; \code{0} means no
#'   censoring. If \code{NULL}, censoring times are uniform on
#'   \code{(cens.lower, cens.upper)}.
#' @param cens.lower Lower bound for uniform censoring time when
#'   \code{cens.rate} is \code{NULL}.
#' @param cens.upper Upper bound for uniform censoring time when
#'   \code{cens.rate} is \code{NULL}.
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
#'   \code{"Friedman"}, \code{"SharedLatent"} and pure-interaction outcomes it
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
#'     nsample = 60, template = tmpl, p = 10, nrep = 2,
#'     ygen.mode = "Interaction", outcome.type = "binary", seed = 1
#' )
#' head(sim$trainDat$Rep_1$sample_metadata)
#' sim$trainDat$Rep_1$truth$interaction_pairs
#' @export
gen_simmba <- function(nsample,
    template,
    p = NULL,
    crosstalk = 1,
    feature.select = c("random", "variance"),
    feature.jitter = 0.2,
    view.models = NULL,
    snr = 1,
    p.train = 0.7,
    de.prob = 0.1,
    de.downProb = 0.5,
    de.facLoc = 1,
    de.facScale = 0.4,
    de.align = NULL,
    signal.share = NULL,
    friedman.share = NULL,
    ygen.mode = c("LM", "Friedman", "Friedman2", "SharedLatent", "Interaction"),
    shared.latent.fun = c("linear", "friedman", "friedman2"),
    shared.latent.n_factors = NULL,
    nonlinear.weight = 0.5,
    interaction.weight = 0.5,
    n.interactions = NULL,
    outcome.type = c("continuous", "binary", "survival"),
    binary.prevalence = 0.5,
    surv.hscale = 1,
    cens.rate = NULL,
    cens.lower = 1,
    cens.upper = 3,
    nrep = 1,
    seed = 1234,
    return.template = FALSE,
    back.transform = FALSE,
    ...) {

    .check_count(nsample, "nsample", min = 2L)
    if (!is.null(seed)) withr::local_seed(seed)
    ygen.mode      <- match.arg(ygen.mode)
    shared.latent.fun <- match.arg(shared.latent.fun)
    outcome.type   <- match.arg(outcome.type)
    feature.select <- match.arg(feature.select)
    if (!is.null(shared.latent.n_factors)) {
        if (!is.numeric(shared.latent.n_factors) || length(shared.latent.n_factors) != 1L ||
            !is.finite(shared.latent.n_factors) || shared.latent.n_factors < 1 ||
            abs(shared.latent.n_factors - round(shared.latent.n_factors)) > 1e-8) {
            stop("shared.latent.n_factors must be NULL or a positive integer.",
                call. = FALSE)
        }
        shared.latent.n_factors <- as.integer(shared.latent.n_factors)
    }
    .check_weight(nonlinear.weight, "nonlinear.weight")
    .check_weight(interaction.weight, "interaction.weight")
    if (!is.null(n.interactions) && (!is.numeric(n.interactions) || length(n.interactions) != 1L ||
        !is.finite(n.interactions) || n.interactions < 1 ||
        abs(n.interactions - round(n.interactions)) > 1e-8)) {
        stop("n.interactions must be NULL or a positive integer.",
            call. = FALSE)
    }
    if (!is.numeric(nrep) || length(nrep) != 1L || !is.finite(nrep) || nrep < 1 ||
        abs(nrep - round(nrep)) > 1e-8) {
        stop("nrep must be a single positive integer.", call. = FALSE)
    }
    nrep <- as.integer(round(nrep))
    if (!is.numeric(nsample) || length(nsample) != 1L || !is.finite(nsample) || nsample < 2) {
        stop("nsample must be a single integer of at least 2.", call. = FALSE)
    }
    nsample <- as.integer(nsample)
    if (!is.numeric(p.train) || length(p.train) != 1L || !is.finite(p.train) ||
        p.train <= 0 || p.train >= 1) {
        stop("p.train must be strictly between 0 and 1.", call. = FALSE)
    }
    n_train <- round(nsample * p.train)
    if (n_train < 1L || n_train > nsample - 1L) {
        stop("p.train = ", p.train, " with nsample = ", nsample,
            " leaves an empty training or test set.", call. = FALSE)
    }
    if (!is.numeric(snr) || length(snr) != 1L || !is.finite(snr) || snr <= 0) {
        stop("snr must be a single positive number.", call. = FALSE)
    }
    if (!is.numeric(binary.prevalence) || length(binary.prevalence) != 1L ||
        !is.finite(binary.prevalence) || binary.prevalence <= 0 || binary.prevalence >= 1) {
        stop("binary.prevalence must be strictly between 0 and 1.",
            call. = FALSE)
    }
    if (!is.null(cens.rate) && (!is.numeric(cens.rate) || length(cens.rate) != 1L ||
        !is.finite(cens.rate) || cens.rate < 0 || cens.rate >= 1)) {
        stop("cens.rate must be NULL or a number in [0, 1).", call. = FALSE)
    }

    ## 0. Fit the template once (skipped if a fitted template is supplied)
    tmpl  <- if (inherits(template,
        "mimesys_template")) template else fit_template(template, ...)
    views <- tmpl$views
    M     <- length(views)
    de.prob     <- .per_view_values(de.prob, views, "de.prob")
    de.downProb <- .per_view_values(de.downProb, views, "de.downProb")
    de.facLoc   <- .per_view_values(de.facLoc, views, "de.facLoc")
    de.facScale <- .per_view_values(de.facScale, views, "de.facScale")
    .check_unit_interval(de.prob, "de.prob")
    .check_unit_interval(de.downProb, "de.downProb")
    if (any(!is.finite(de.facLoc))) stop("de.facLoc must be finite.", call. = FALSE)
    if (any(!is.finite(de.facScale)) || any(de.facScale < 0)) {
        stop("de.facScale must be finite and non-negative.", call. = FALSE)
    }
    if (!is.null(de.align)) .check_unit_interval(de.align, "de.align", single = TRUE)
    gamma_named <- .per_view_values(crosstalk, views, "crosstalk")
    signal.share <- .normalize_signal_share(signal.share, views, "signal.share")
    friedman.share <- .normalize_signal_share(friedman.share, views,
        "friedman.share")
    if (is.null(friedman.share)) {
        friedman.share <- signal.share
    }
    if (identical(ygen.mode, "Interaction") && M < 2L) {
        stop("ygen.mode = 'Interaction' needs at least two views.",
            call. = FALSE)
    }
    if (identical(ygen.mode, "SharedLatent")) {
        gamma_check <- gamma_named
        if (tmpl$shared.rank == 0L) {
            warning("The template has no shared latent factors; SharedLatent outcomes fall back ",
                "to the average of the view-specific linear predictors.", call. = FALSE)
        } else if (all(gamma_check == 0)) {
            warning("ygen.mode = 'SharedLatent' with crosstalk = 0: shared factors do not enter ",
                "any simulated feature, so the outcome is unrelated to the features.", call. = FALSE)
        }
    }

    trainDat <- testDat <- agreement.emp <- vector("list", nrep)
    names(trainDat) <- names(testDat) <- names(agreement.emp) <- paste("Rep",
        seq_len(nrep), sep = "_")

    for (k in seq_len(nrep)) {

        ## 1. Null multiview data from the template
        pcl <- .sim_from_template_impl(tmpl, n = nsample, p = p,
            crosstalk = crosstalk,
            feature.select = feature.select,
            feature.jitter = feature.jitter,
            view.models = view.models,
            back.transform = if (isTRUE(back.transform)) "both" else FALSE)
        X  <- t(as.matrix(pcl$feature_table))
        fm <- pcl$feature_metadata

        ## 2. Spike-in coefficients, view by view in template order
        ## (table(featureType) would sort views alphabetically and misassign de.*)
        beta0 <- unlist(lapply(seq_len(M), function(m) {
            sel <- fm$featureType == views[m]
            .draw_beta_view(sum(sel), fm$shared_loading[sel], de.prob[m],
                de.downProb[m],
                de.facLoc[m], de.facScale[m], de.align)
        }))
        ## 3. Linear predictor on standardized features, so coefficients are
        ## per-SD effects and views on very different measurement scales (e.g.
        ## counts vs. intensities in the millions) contribute comparably. Then
        ## optional view-level signal balancing.
        Xs <- .standardize_columns(X)
        view_pred <- vapply(views, function(v) {
            sel <- fm$featureType == v
            as.numeric(Xs[, sel, drop = FALSE] %*% beta0[sel])
            }, numeric(nsample))
        view_pred <- matrix(view_pred, nrow = nsample, dimnames = list(NULL,
            views))
        if (!is.null(signal.share)) {
            scaled <- .rescale_signal_by_view(beta0, view_pred, fm$featureType,
                views, signal.share)
            beta0 <- scaled$beta
            view_pred <- scaled$view_pred
        }
        eta_lin <- rowSums(view_pred)
        shared_res <- .shared_latent_predictor(
            pcl,
            view_pred,
            fun = shared.latent.fun,
            n_factors = shared.latent.n_factors,
            nonlinear.weight = nonlinear.weight
        )
        eta_shared <- as.numeric(shared_res)
        shared_fallback <- isTRUE(attr(shared_res, "fallback"))
        A <- suppressWarnings(stats::cor(view_pred))
        A[!is.finite(A)] <- 0
        agreement.emp[[k]] <- A

        nfeat <- nrow(fm)
        friedman_term <- interaction_term <- rep(FALSE, nfeat)
        interaction_pairs <- NULL

        ## Optional non-linear Friedman component
        if (ygen.mode %in% c("Friedman", "Friedman2")) {
            if (length(beta0) < 5L) {
                stop("Friedman outcomes need at least five simulated features in total; ",
                    "got ", length(beta0), ". Increase p.", call. = FALSE)
            }
            nonzero_index <- which(beta0 != 0)
            candidate_index <- if (length(nonzero_index) >= 5L) nonzero_index else seq_along(beta0)
            friedman_index <- .sample_friedman_features(candidate_index,
                fm$featureType, views, friedman.share)
            friedman_term[friedman_index] <- TRUE
            X.friedman <- .rank_to_unit_interval(X[, friedman_index,
                drop = FALSE])
            Xbeta.friedman <- .friedman_response(X.friedman)
            pcl$sample_metadata$Xbeta_friedman <- Xbeta.friedman
        }

        ## Optional cross-view interaction component
        if (identical(ygen.mode, "Interaction")) {
            interaction_pairs <- .sample_interaction_pairs(beta0,
                fm$featureType, views, n.interactions)
            eta_int <- rowSums(
                Xs[, interaction_pairs$i, drop = FALSE] * Xs[,
                    interaction_pairs$j, drop = FALSE]
            )
            interaction_term[c(interaction_pairs$i,
                interaction_pairs$j)] <- TRUE
            interaction_pairs <- data.frame(
                feature_x = fm$featureID[interaction_pairs$i],
                view_x = fm$featureType[interaction_pairs$i],
                feature_y = fm$featureID[interaction_pairs$j],
                view_y = fm$featureType[interaction_pairs$j],
                stringsAsFactors = FALSE
            )
            pcl$sample_metadata$Xbeta_interaction <- eta_int
        }

        pcl$sample_metadata$Xbeta_linear <- eta_lin
        pcl$sample_metadata$Xbeta_shared <- eta_shared

        ## 4. Final linear predictor
        eta <- switch(ygen.mode,
            LM           = eta_lin,
            Friedman     = Xbeta.friedman,
            Friedman2    = .mix_signals(eta_lin, Xbeta.friedman,
                nonlinear.weight),
            SharedLatent = eta_shared,
            Interaction  = .mix_signals(eta_lin, eta_int, interaction.weight))
        pcl$sample_metadata$Xbeta <- eta

        ## Ground truth: beta is kept only where the linear predictor enters Y
        uses_linear <- switch(ygen.mode,
            LM = TRUE,
            Friedman = FALSE,
            Friedman2 = nonlinear.weight < 1,
            SharedLatent = shared_fallback,
            Interaction = interaction.weight < 1)
        if (!uses_linear) beta0[] <- 0
        pcl$feature_metadata$beta <- beta0
        pcl$feature_metadata$friedman_term <- friedman_term
        pcl$feature_metadata$interaction_term <- interaction_term
        ## Components with zero weight do not enter the outcome
        if (identical(ygen.mode, "Friedman2") && nonlinear.weight <= 0) {
            friedman_term[] <- FALSE
        }
        if (identical(ygen.mode, "Interaction") && interaction.weight <= 0) {
            interaction_term[] <- FALSE
            interaction_pairs <- NULL
        }
        pcl$feature_metadata$friedman_term <- friedman_term
        pcl$feature_metadata$interaction_term <- interaction_term
        pcl$feature_metadata$outcome_feature <- beta0 != 0 | friedman_term | interaction_term
        pcl$truth$interaction_pairs <- interaction_pairs

        ## 5. Outcome
        if (outcome.type == "continuous") {
            sigma2 <- as.vector(stats::var(eta) / snr)
            pcl$sample_metadata$Y <- as.vector(eta + stats::rnorm(nsample) * sqrt(sigma2))
        } else if (outcome.type == "binary") {
            eta <- .standardize_predictor(eta) * sqrt(snr * pi^2 / 3)
            eta <- eta + .logistic_intercept(eta, binary.prevalence)
            pcl$sample_metadata$Xbeta <- eta
            pcl$sample_metadata$Y <- stats::rbinom(nsample, size = 1,
                prob = stats::plogis(eta))
        } else if (outcome.type == "survival") {
            eta <- .standardize_predictor(eta) * sqrt(snr * pi^2 / 6)
            pcl$sample_metadata$Xbeta <- eta
            h  <- as.vector(surv.hscale * exp(eta))
            T0 <- stats::rexp(nsample, rate = h)
            C  <- .draw_censoring(T0, cens.rate, cens.lower, cens.upper)
            pcl$sample_metadata$time   <- pmin(T0, C)
            pcl$sample_metadata$status <- as.integer(C >= T0)
        }

        ## Optionally return features on the original measurement scale. Outcomes
        ## above were generated on the fitted scale.
        if (!is.null(pcl$original_table)) {
            pcl$feature_table <- pcl$original_table
            pcl$original_table <- NULL
            pcl$feature_metadata$scale <- "original"
        }

        ## 6. Train / test split
        train <- test <- pcl
        tr.row <- sample.int(nsample, size = round(nsample * p.train),
            replace = FALSE)
        train$sample_metadata <- pcl$sample_metadata[tr.row, , drop = FALSE]
        test$sample_metadata  <- pcl$sample_metadata[-tr.row, , drop = FALSE]
        train$feature_table   <- pcl$feature_table[, tr.row, drop = FALSE]
        test$feature_table    <- pcl$feature_table[, -tr.row, drop = FALSE]
        train$truth$scores    <- pcl$truth$scores[tr.row, , drop = FALSE]
        test$truth$scores     <- pcl$truth$scores[-tr.row, , drop = FALSE]
        train$truth$shared    <- pcl$truth$shared[tr.row, , drop = FALSE]
        test$truth$shared     <- pcl$truth$shared[-tr.row, , drop = FALSE]

        trainDat[[k]] <- train
        testDat[[k]]  <- test
    }

    out <- list(trainDat = trainDat, testDat = testDat,
        template = if (isTRUE(return.template)) tmpl else NULL,
        template_provenance = tmpl$provenance, p = p, crosstalk = crosstalk,
        view.models = view.models,
        cancor = .realized_cancor(tmpl, crosstalk),
        agreement.emp = agreement.emp,
        snr = snr, p.train = p.train,
        de.prob = de.prob, de.downProb = de.downProb, de.facLoc = de.facLoc,
        de.facScale = de.facScale, de.align = de.align,
        signal.share = signal.share, friedman.share = friedman.share,
        shared.latent.fun = shared.latent.fun,
        shared.latent.n_factors = shared.latent.n_factors,
        nonlinear.weight = nonlinear.weight, back.transform = back.transform,
        interaction.weight = interaction.weight,
        n.interactions = n.interactions,
        nrep = nrep, seed = seed, ygen.mode = ygen.mode,
        outcome.type = outcome.type, binary.prevalence = binary.prevalence,
        surv.hscale = surv.hscale, cens.rate = cens.rate,
        cens.lower = cens.lower, cens.upper = cens.upper)
    out
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

## Intercept b such that mean(plogis(b + eta)) equals the target prevalence
.logistic_intercept <- function(eta, prevalence) {
    f <- function(b) mean(stats::plogis(b + eta)) - prevalence
    stats::uniroot(f, interval = c(-100, 100), tol = 1e-10)$root
}

## Censoring times: uniform(cens.lower, cens.upper), or uniform(0, c) with c
## solved so that the requested proportion of subjects is censored.
.draw_censoring <- function(T0, cens.rate = NULL, cens.lower = 1,
    cens.upper = 3) {
    n <- length(T0)
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
    u * exp(log_c)
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

.check_unit_interval <- function(x, name, single = FALSE) {
    if (!is.numeric(x) || (single && length(x) != 1L) || any(!is.finite(x)) ||
        any(x < 0) || any(x > 1)) {
        stop(name, " must be ", if (single) "a single number" else "number(s)",
            " between 0 and 1.", call. = FALSE)
    }
    invisible(x)
}

.check_weight <- function(w, name) {
    if (!is.numeric(w) || length(w) != 1L || !is.finite(w) || w < 0 || w > 1) {
        stop(name, " must be a single number between 0 and 1.", call. = FALSE)
    }
    invisible(w)
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
                    "shared.latent.n_factors = ", n_factors,
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
                "shared.latent.fun = '", fun,
                "' requires at least five selected shared latent factors; only ",
                ncol(shared_scaled), " are available. Increase shared.rank in fit_template() ",
                "or use shared.latent.fun = 'linear'.",
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
