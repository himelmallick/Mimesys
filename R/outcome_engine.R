## Outcome engine shared by gen_simmba() (one study) and gen_multistudy().
##
## One replicate: simulate every study's features from the template (all
## studies share the same simulated features), draw the coefficients, build the
## linear / nonlinear / shared-latent predictor on the pooled samples and draw
## the outcome. With a single study and no heterogeneity, the random draws and
## the arithmetic are exactly those of the original single-study gen_simmba().
## Heterogeneity is drawn from side streams (see .with_side_stream()), so it
## never changes the main random stream.
##
## Returns the pooled pcl (samples in study order, before the train/test
## split), the study index of every sample, the empirical agreement between
## view predictors and the per-study parameters that were drawn.

.no_heterogeneity <- function() {
    list(effect.cor = 1, support.overlap = 1, intercept.sd = 0, snr.sd = 0,
        shift.sd = 0, batch.loc.sd = 0, batch.scale.sd = 0)
}

.simulate_outcome_replicate <- function(tmpl, n_by_study, cfg,
    het = .no_heterogeneity(), study_names = NULL) {
    S <- length(n_by_study)
    multi <- !is.null(study_names)
    ## One template for every study, or one template per study (same views and
    ## features, e.g. from fit_study_templates()); crosstalk per study
    tmpls <- if (inherits(tmpl, "mimesys_template")) rep(list(tmpl), S) else tmpl
    per_template <- !inherits(tmpl, "mimesys_template")
    crosstalk_s <- if (is.list(cfg$crosstalk)) cfg$crosstalk else rep(list(cfg$crosstalk), S)
    tmpl <- tmpls[[1]]
    views <- tmpl$views
    M <- length(views)

    ## 1. Features of every study. Study 1 chooses the simulated features;
    ## later studies reuse them. Study-level case-mix shifts and batch effects
    ## are drawn inside the hooks, right after that study's latent draws.
    pcls <- shifts <- batches <- vector("list", S)
    structure <- NULL
    for (s in seq_len(S)) {
        draws <- new.env(parent = emptyenv())
        shift_hook <- if (het$shift.sd > 0) function(Z_list, ctx) {
            res <- .with_side_stream(1000L + s,
                .draw_case_mix_shift(Z_list, ctx, het$shift.sd))
            draws$shift <- res$shift
            res$Z
        } else NULL
        batch_hook <- if (het$batch.loc.sd > 0 || het$batch.scale.sd > 0) {
            function(Z_list, ctx) {
                res <- .with_side_stream(2000L + s,
                    .draw_batch_effects(Z_list, het$batch.loc.sd,
                        het$batch.scale.sd))
                draws$batch <- res$batch
                res$Z
            }
        } else NULL
        structure_s <- if (per_template && s > 1L) {
            .map_structure(structure, tmpls[[1]], tmpls[[s]], cfg$feature.jitter)
        } else structure
        pcl_s <- .sim_from_template_impl(tmpls[[s]], n = n_by_study[s], p = cfg$p,
            crosstalk = crosstalk_s[[s]],
            feature.select = cfg$feature.select,
            feature.jitter = cfg$feature.jitter,
            view.models = cfg$view.models,
            back.transform = if (isTRUE(cfg$back.transform)) "both" else FALSE,
            .structure = structure_s, .latent_hook = shift_hook,
            .observed_hook = batch_hook)
        if (is.null(structure)) structure <- attr(pcl_s, "structure")
        if (!is.null(draws$shift)) {
            ## Shifted factor scores are part of the truth
            pcl_s$truth$shared <- sweep(pcl_s$truth$shared, 2,
                draws$shift$shared, "+")
            pcl_s$truth$scores <- sweep(pcl_s$truth$scores, 2,
                draws$shift$scores, "+")
        }
        pcls[[s]] <- pcl_s
        shifts[s] <- list(draws$shift)
        batches[s] <- list(draws$batch)
    }
    sid <- rep(seq_len(S), n_by_study)
    pcl <- if (S == 1L && !multi) pcls[[1]] else
        .pool_study_pcls(pcls, study_names)
    n <- length(sid)
    X  <- t(as.matrix(pcl$feature_table))
    fm <- pcl$feature_metadata
    ygen.mode <- cfg$ygen.mode

    ## 2. Spike-in coefficients, view by view in template order
    ## (table(featureType) would sort views alphabetically and misassign de.*)
    beta0 <- unlist(lapply(seq_len(M), function(m) {
        sel <- fm$featureType == views[m]
        .draw_beta_view(sum(sel), fm$shared_loading[sel], cfg$signal.prob[m],
            cfg$signal.neg.prob[m],
            cfg$effect.meanlog[m], cfg$effect.sdlog[m], cfg$signal.align)
    }))
    ## Study coefficients: the common draw, or study-specific deviations
    B <- if (het$effect.cor >= 1 && het$support.overlap >= 1) {
        matrix(beta0, nrow = length(beta0), ncol = S)
    } else {
        .with_side_stream(3000L, .study_coefficients(beta0, fm$featureType,
            views, S, het, cfg, loading = fm$shared_loading))
    }

    ## 3. Linear predictor on standardized features (pooled over studies), so
    ## coefficients are per-SD effects and views on very different measurement
    ## scales contribute comparably. Then optional view-level signal balancing
    ## (one factor per view for all studies).
    Xs <- .standardize_columns(X)
    view_pred <- matrix(0, nrow = n, ncol = M, dimnames = list(NULL, views))
    for (s in seq_len(S)) {
        rows <- which(sid == s)
        Xs_s <- if (S == 1L) Xs else Xs[rows, , drop = FALSE]
        vp <- vapply(views, function(v) {
            sel <- fm$featureType == v
            as.numeric(Xs_s[, sel, drop = FALSE] %*% B[sel, s])
        }, numeric(length(rows)))
        view_pred[rows, ] <- matrix(vp, nrow = length(rows))
    }
    if (!is.null(cfg$signal.share)) {
        scaled <- .rescale_signal_pooled(B, view_pred, fm$featureType, views,
            cfg$signal.share, sid)
        B <- scaled$beta
        view_pred <- scaled$view_pred
    }
    eta_lin <- rowSums(view_pred)
    ## Common signal features (a view with signal.share 0 carries none); the
    ## Friedman features and interaction pairs are drawn from these.
    beta_common <- beta0
    if (!is.null(cfg$signal.share)) {
        beta_common[fm$featureType %in% views[cfg$signal.share[views] == 0]] <- 0
    }
    shared_res <- .shared_latent_predictor(
        pcl,
        view_pred,
        fun = cfg$shared.latent.function,
        n_factors = cfg$shared.latent.n.factors,
        nonlinear.weight = cfg$nonlinear.weight
    )
    eta_shared <- as.numeric(shared_res)
    shared_fallback <- isTRUE(attr(shared_res, "fallback"))
    ## Optional view-specific latent signal mixed into the shared-latent
    ## outcome (nothing is computed at the default weight of 1)
    w_shared <- if (is.null(cfg$shared.latent.weight)) 1 else cfg$shared.latent.weight
    eta_specific <- specific_scores <- NULL
    if (identical(ygen.mode, "SharedLatent") && w_shared < 1) {
        specific_scores <- do.call(rbind, lapply(seq_len(S), function(s) {
            .view_specific_scores(pcls[[s]]$truth, tmpls[[s]])
        }))
        specific_scores <- apply(specific_scores, 2, .pooled_standardize, sid = sid)
        specific_scores <- matrix(specific_scores, nrow = n, dimnames = list(NULL, views))
        share <- if (is.null(cfg$signal.share)) rep(1 / M, M) else cfg$signal.share[views]
        eta_specific <- .pooled_standardize(as.numeric(specific_scores %*% sqrt(share)), sid)
    }
    A <- suppressWarnings(stats::cor(view_pred))
    A[!is.finite(A)] <- 0

    nfeat <- nrow(fm)
    friedman_term <- interaction_term <- rep(FALSE, nfeat)
    interaction_pairs <- NULL

    ## Optional non-linear Friedman component (same features in every study)
    if (ygen.mode %in% c("Friedman", "Friedman2")) {
        if (length(beta0) < 5L) {
            stop("Friedman outcomes need at least five simulated features in total; ",
                "got ", length(beta0), ". Increase p.", call. = FALSE)
        }
        nonzero_index <- which(beta_common != 0)
        candidate_index <- if (length(nonzero_index) >= 5L) nonzero_index else seq_along(beta0)
        friedman_index <- .sample_friedman_features(candidate_index,
            fm$featureType, views, cfg$friedman.share)
        friedman_term[friedman_index] <- TRUE
        X.friedman <- .rank_to_unit_interval(X[, friedman_index,
            drop = FALSE])
        Xbeta.friedman <- .friedman_response(X.friedman)
        pcl$sample_metadata$Xbeta_friedman <- Xbeta.friedman
    }

    ## Optional cross-view interaction component (same pairs in every study)
    if (identical(ygen.mode, "Interaction")) {
        interaction_pairs <- .sample_interaction_pairs(beta_common,
            fm$featureType, views, cfg$n.interactions)
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
    if (!is.null(eta_specific)) {
        pcl$sample_metadata$Xbeta_specific <- eta_specific
        rownames(specific_scores) <- pcl$sample_metadata$subjectID
        pcl$truth$specific_scores <- specific_scores
    }

    ## 4. Final predictor
    eta <- switch(ygen.mode,
        LM           = eta_lin,
        Friedman     = Xbeta.friedman,
        Friedman2    = .mix_signals(eta_lin, Xbeta.friedman,
            cfg$nonlinear.weight),
        SharedLatent = if (is.null(eta_specific)) eta_shared else
            .mix_signals(eta_specific, eta_shared, w_shared),
        Interaction  = .mix_signals(eta_lin, eta_int, cfg$interaction.weight))
    pcl$sample_metadata$Xbeta <- eta

    ## Ground truth: beta is kept only where the linear predictor enters Y
    uses_linear <- switch(ygen.mode,
        LM = TRUE,
        Friedman = FALSE,
        Friedman2 = cfg$nonlinear.weight < 1,
        SharedLatent = shared_fallback,
        Interaction = cfg$interaction.weight < 1)
    if (!uses_linear) B[] <- 0
    ## Components with zero weight do not enter the outcome
    if (identical(ygen.mode, "Friedman2") && cfg$nonlinear.weight <= 0) {
        friedman_term[] <- FALSE
    }
    if (identical(ygen.mode, "Interaction") && cfg$interaction.weight <= 0) {
        interaction_term[] <- FALSE
        interaction_pairs <- NULL
    }
    if (multi) {
        pcl$feature_metadata$beta <- rowMeans(B)
        for (s in seq_len(S)) {
            pcl$feature_metadata[[paste0("beta_", study_names[s])]] <- B[, s]
        }
        pcl$feature_metadata$common_signal <- beta_common != 0 & uses_linear
    } else {
        pcl$feature_metadata$beta <- B[, 1]
    }
    pcl$feature_metadata$friedman_term <- friedman_term
    pcl$feature_metadata$interaction_term <- interaction_term
    pcl$feature_metadata$outcome_feature <- rowSums(B != 0) > 0 |
        friedman_term | interaction_term
    pcl$truth$interaction_pairs <- interaction_pairs

    ## 5. Study intercepts and signal-to-noise ratios
    alpha <- if (het$intercept.sd > 0) .with_side_stream(4000L,
        stats::rnorm(S, 0, het$intercept.sd)) else rep(0, S)
    snr_s <- if (het$snr.sd > 0) exp(log(cfg$snr) + .with_side_stream(5000L,
        stats::rnorm(S, 0, het$snr.sd))) else rep(cfg$snr, S)
    has_alpha <- any(alpha != 0)

    ## 6. Outcome, study by study
    outcome_intercept <- NULL
    if (cfg$outcome.type == "continuous") {
        v_ref <- .within_study_variance(eta, sid)
        Y <- numeric(n)
        if (has_alpha) {
            eta <- eta + alpha[sid] * sqrt(v_ref)
            pcl$sample_metadata$Xbeta <- eta
        }
        for (s in seq_len(S)) {
            rows <- which(sid == s)
            sigma2 <- as.vector(v_ref / snr_s[s])
            Y[rows] <- as.vector(eta[rows] + stats::rnorm(length(rows)) * sqrt(sigma2))
        }
        pcl$sample_metadata$Y <- Y
    } else if (cfg$outcome.type == "binary") {
        eta <- .pooled_standardize(eta, sid) * sqrt(snr_s[sid] * pi^2 / 3)
        ## Study intercepts first, then the common intercept that gives the
        ## target prevalence in the pooled sample
        if (has_alpha) eta <- eta + alpha[sid]
        outcome_intercept <- .logistic_intercept(eta, cfg$binary.prevalence)
        eta <- eta + outcome_intercept
        pcl$sample_metadata$Xbeta <- eta
        Y <- integer(n)
        for (s in seq_len(S)) {
            rows <- which(sid == s)
            Y[rows] <- stats::rbinom(length(rows), size = 1,
                prob = stats::plogis(eta[rows]))
        }
        pcl$sample_metadata$Y <- Y
    } else if (cfg$outcome.type == "survival") {
        eta <- .pooled_standardize(eta, sid) * sqrt(snr_s[sid] * pi^2 / 6)
        if (has_alpha) eta <- eta + alpha[sid]
        pcl$sample_metadata$Xbeta <- eta
        h <- as.vector(cfg$baseline.hazard * exp(eta))
        time <- numeric(n)
        status <- integer(n)
        ## Censoring target: one value, or one per study
        cens_rate <- cfg$cens.rate
        if (length(cens_rate) > 1L && length(cens_rate) != S) {
            stop("cens.rate has ", length(cens_rate), " values for ", S,
                " studies; give one value or one per study.", call. = FALSE)
        }
        cens_type <- if (is.null(cfg$cens.type)) "uniform" else cfg$cens.type
        cens_dep <- if (is.null(cfg$cens.dependence)) 0 else cfg$cens.dependence
        shape <- if (is.null(cfg$baseline.shape)) 1 else cfg$baseline.shape
        keep_times <- isTRUE(cfg$keep.true.times)
        if (keep_times) event_time <- cens_time <- numeric(n)
        for (s in seq_len(S)) {
            rows <- which(sid == s)
            T0 <- stats::rexp(length(rows), rate = h[rows])
            ## Weibull baseline: E^(1/k) for E ~ Exp(h) has hazard h k t^(k-1)
            if (shape != 1) T0 <- T0^(1 / shape)
            C  <- .draw_censoring(T0,
                if (length(cens_rate) > 1L) cens_rate[[s]] else cens_rate,
                cfg$cens.lower, cfg$cens.upper, cens.type = cens_type,
                cens.times = .resolve_censoring_times(cfg$cens.times, rows, n, s),
                cens.dependence = cens_dep, eta = eta[rows],
                cens.accrual = cfg$cens.accrual, cens.dropout = cfg$cens.dropout)
            time[rows] <- pmin(T0, C)
            status[rows] <- as.integer(C >= T0)
            if (keep_times) {
                event_time[rows] <- T0
                cens_time[rows] <- C
            }
        }
        pcl$sample_metadata$time   <- time
        pcl$sample_metadata$status <- status
        if (keep_times) {
            pcl$sample_metadata$event_time <- event_time
            pcl$sample_metadata$cens_time <- cens_time
        }
    }

    ## 7. Returned features: batch effects (if any) are added after the outcome
    ## was generated, so they are pure technical noise. Optionally return the
    ## original measurement scale; outcomes always come from the fitted scale.
    if (!is.null(pcl$observed_table)) {
        pcl$feature_table <- pcl$observed_table
        if (!is.null(pcl$observed_original_table)) {
            pcl$original_table <- pcl$observed_original_table
        }
    }
    pcl$observed_table <- pcl$observed_original_table <- NULL
    if (!is.null(pcl$original_table)) {
        pcl$feature_table <- pcl$original_table
        pcl$original_table <- NULL
        pcl$feature_metadata$scale <- "original"
    }

    list(pcl = pcl, study = sid, agreement = A, beta = B, common_beta = beta_common,
        intercept = alpha, snr = snr_s, shift = shifts, batch = batches,
        outcome_intercept = outcome_intercept)
}

## Pool the per-study simulations (same features) into one pcl with study
## labels and unique sample IDs.
.pool_study_pcls <- function(pcls, study_names) {
    first <- pcls[[1]]
    ids <- unlist(lapply(seq_along(pcls), function(s) {
        paste0(study_names[s], "_", pcls[[s]]$sample_metadata$subjectID)
    }))
    bind_tables <- function(name) {
        tabs <- lapply(pcls, `[[`, name)
        if (is.null(tabs[[1]])) return(NULL)
        out <- do.call(cbind, tabs)
        colnames(out) <- ids
        out
    }
    sample_metadata <- data.frame(subjectID = ids,
        study = rep(study_names, vapply(pcls, function(x) nrow(x$sample_metadata), integer(1))),
        row.names = ids, stringsAsFactors = FALSE)
    out <- list(feature_table = bind_tables("feature_table"),
        sample_metadata = sample_metadata,
        feature_metadata = first$feature_metadata,
        truth = first$truth)
    ## Latent truth: pooled when every study has the same latent dimensions
    ## (always with one template); with per-study templates of different
    ## ranks it is not comparable across studies and is left out
    for (name in c("scores", "shared")) {
        parts <- lapply(pcls, function(x) x$truth[[name]])
        same <- length(unique(vapply(parts, ncol, integer(1)))) == 1L
        out$truth[name] <- list(if (same) do.call(rbind, parts) else NULL)
    }
    out$truth$crosstalk <- lapply(pcls, function(x) x$truth$crosstalk)
    out$truth$cancor <- lapply(pcls, function(x) x$truth$cancor)
    names(out$truth$crosstalk) <- names(out$truth$cancor) <- study_names
    for (name in c("original_table", "observed_table", "observed_original_table")) {
        out[name] <- list(bind_tables(name))
    }
    out
}

## Pooled within-study variance (the variance of a single study when S = 1)
.within_study_variance <- function(x, sid) {
    if (length(unique(sid)) == 1L) return(stats::var(x))
    dev <- x - stats::ave(x, sid)
    sum(dev^2) / (length(x) - length(unique(sid)))
}

## Centre on the pooled mean and scale by the pooled within-study SD, so study
## differences in the predictor are kept. One study: ordinary standardization.
.pooled_standardize <- function(x, sid) {
    if (length(unique(sid)) == 1L) return(.standardize_predictor(x))
    x <- as.numeric(x)
    v <- .within_study_variance(x, sid)
    if (!is.finite(v) || v <= 0) return(x - mean(x))
    (x - mean(x)) / sqrt(v)
}

## View-level signal balancing with one factor per view for all studies: the
## factor makes the pooled within-study variance of the view predictor equal
## to its share. With one study this is .rescale_signal_by_view().
.rescale_signal_pooled <- function(B, view_pred, feature_type, views,
    signal.share, sid) {
    B <- as.matrix(B)
    view_pred <- as.matrix(view_pred)
    colnames(view_pred) <- views
    for (view in views) {
        target <- signal.share[[view]]
        sel <- feature_type == view
        if (target == 0) {
            B[sel, ] <- 0
            view_pred[, view] <- 0
            next
        }
        current_var <- .within_study_variance(view_pred[, view], sid)
        if (!is.finite(current_var) || current_var <= 0) {
            warning("Cannot rescale signal for view '", view,
                "' because its linear predictor variance is zero.", call. = FALSE)
            next
        }
        scale_factor <- sqrt(target / current_var)
        B[sel, ] <- B[sel, ] * scale_factor
        view_pred[, view] <- view_pred[, view] * scale_factor
    }
    list(beta = B, view_pred = view_pred)
}

## View-specific latent signal of one study: in each view, the part of the
## latent factor scores that is not shared with the other views,
## F_m - sqrt(gamma_m) S B_m' = E_m R_m (independent of the shared variates and
## across views), combined over the view's factors with weights proportional
## to the factor's loading strength (column sums of squared loadings). Returns
## an n x M matrix.
.view_specific_scores <- function(truth, tmpl) {
    ks <- vapply(tmpl$fits, function(f) f$rank, integer(1))
    block <- rep(seq_along(ks), ks)
    gamma <- truth$crosstalk
    r <- if (is.matrix(truth$shared)) ncol(truth$shared) else 0L
    out <- vapply(seq_along(ks), function(m) {
        Fm <- truth$scores[, block == m, drop = FALSE]
        if (r > 0L && gamma[[m]] > 0) {
            Fm <- Fm - truth$shared %*% t(sqrt(gamma[[m]]) * tmpl$B[[m]])
        }
        w <- colSums(tmpl$fits[[m]]$Lambda^2)
        as.numeric(Fm %*% (w / sqrt(sum(w^2))))
    }, numeric(nrow(truth$scores)))
    matrix(out, ncol = length(ks), dimnames = list(NULL, tmpl$views))
}
