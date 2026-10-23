#' Fit a template model for null multiview simulation
#'
#' Learns, from any matched multiview template, the three ingredients needed to
#' generate realistic null data: (i) each feature's marginal distribution
#' (empirical quantile function, so skewness, bounds, ties and zero inflation
#' carry over), (ii) within-view dependence (low-rank factor model on normal
#' scores), and (iii) cross-view dependence (a low-dimensional shared latent
#' space estimated by generalized CCA across the view-level factor scores).
#'
#' The model is a Gaussian copula with factor structure. For view m,
#' \deqn{f_m = \sqrt{\gamma_m} B_m s + (I - \gamma_m B_m B_m^T)^{1/2}
#' \epsilon_m,}
#' \deqn{z_m = \Lambda_m f_m + e_m, \quad x_{mj} = Q_{mj}(\Phi(z_{mj})),}
#' with shared variates \eqn{s \sim N(0, I_r)}, \eqn{Var(f_m) = I} and
#' \eqn{e_{mj} \sim N(0, 1 - ||\Lambda_{mj}||^2)}. The cross-view factor
#' covariance
#' is \eqn{\sqrt{\gamma_m \gamma_{m'}} B_m B_{m'}^T}, so the crosstalk weights
#' \eqn{\gamma} change how strongly views talk to each other without touching
#' marginals or within-view structure. \eqn{\gamma = 1} reproduces the
#' template's
#' shared structure and \eqn{\gamma = 0} gives independent views.
#'
#' @param template Named list of numeric matrices, samples in rows and features
#'   in columns. Rows are matched by rownames when present, otherwise by
#'   position.
#' @param rank Optional factors per view (scalar or one per view). Use
#'   \code{"2logp"} to retain \code{ceiling(2 * log(p_m))} factors in view
#'   \code{m}, subject to \code{max.rank}. If NULL, the smallest rank explaining
#'   \code{var.explained} of normal-score variance is selected, capped at
#'   \code{max.rank}.
#' @param var.explained Target variance fraction when \code{rank} is NULL.
#' @param max.rank Cap on the per-view rank (default floor(n0 / 3)). The cap
#'   also applies to an explicit \code{rank}.
#' @param shared.rank Number of shared variates r. If NULL, chosen by
#'   sequential (deflated) parallel analysis: at each step the largest
#'   generalized CCA eigenvalue of the current view scores is compared with
#'   the permutation distribution of the largest eigenvalue obtained after
#'   permuting samples independently within each view. If it is significant
#'   (p < 0.05), that shared direction is removed from every view and the test
#'   is repeated; the rank is the number of significant steps. Comparing each
#'   eigenvalue with the null of the eigenvalue of the same index instead
#'   over-selects, because a real shared factor pushes all later observed
#'   eigenvalues up.
#' @param seed Random seed (default \code{1234}, the same default as every
#'   other MimeSys function). Tie-breaking in the normal scores and the
#'   shared-rank permutations are random, so the seed makes the fit
#'   reproducible. The caller's random-number state is restored afterwards;
#'   use \code{NULL} to draw from the current random-number stream.
#' @param shared.rank.runs Number of stability-selection runs used when
#'   \code{shared.rank = NULL} (default 1). The parallel-analysis choice is
#'   repeated this many times, each with fresh tie-breaking and permutations,
#'   and the rounded mean of the selected ranks is used; the runs, their
#'   distribution and the share of runs equal to the chosen rank are stored in
#'   \code{shared.rank.selection}, with a message when fewer than half agree.
#'   The default single run warns when the choice is borderline; use, e.g., 20
#'   runs to stabilize a borderline choice.
#' @param nperm Permutations per parallel-analysis step (default 200). The
#'   permutation p-value of each step is stored in \code{gcca.pvalue} (one
#'   per step tested, so the selected rank plus the first non-significant
#'   step), the 95th percentile of each step's null in \code{gcca.threshold},
#'   the tested eigenvalues in \code{gcca.step.eigen}, and a warning is given
#'   when the selected rank rests on a borderline step.
#' @param view.transform Optional transform per view before fitting. Supported
#'   values are `"none"`, `"clr"`, and `"log1p"`. For microbiome data, use
#'   `"clr"` to fit and simulate on centered log-ratio scale.
#' @return An object of class \code{mimesys_template}. The template keeps the
#'   empirical matrices used for fitting (\code{observed}) and each feature's
#'   sorted values (the empirical quantile functions), because crosstalk
#'   calibration and marginal reconstruction need them. Treat a fitted template
#'   with the same care as the source data: do not share it, or simulation
#'   output that embeds it, when the source data are confidential.
#' @examples
#' data(stelzer_mae)
#' tmpl <- fit_template(
#'     stelzer_mae,
#'     view.transform = c(Metabolomics = "log1p", Proteomics = "log1p")
#' )
#' tmpl
#' @export
fit_template <- function(template, rank = NULL, var.explained = 0.7,
    max.rank = NULL, shared.rank = NULL, nperm = 200,
    view.transform = NULL, seed = 1234, shared.rank.runs = 1) {

    ## Random tie-breaking and the shared-rank permutations use random
    ## numbers; a seed makes the fit reproducible without touching the
    ## caller's random-number state.
    if (!is.null(seed)) withr::local_seed(seed)
    .check_count(shared.rank.runs, "shared.rank.runs", min = 1L)
    shared.rank.runs <- as.integer(round(shared.rank.runs))
    .check_count(nperm, "nperm", min = 1L)
    .check_fraction(var.explained, "var.explained", include.lower = FALSE)
    template <- as_mimesys_view_list(template)
    args <- list(template = template, rank = rank, var.explained = var.explained,
        max.rank = max.rank, nperm = nperm, view.transform = view.transform)
    if (!is.null(shared.rank) || length(template) < 2L || shared.rank.runs == 1L) {
        return(do.call(.fit_template_once, c(args, list(shared.rank = shared.rank))))
    }
    ## Stability selection: repeat the parallel-analysis choice of the shared
    ## rank (each run redraws tie-breaking and permutations) and use the
    ## rounded mean. The mean is far less sensitive than the mode when the
    ## runs spread over several neighbouring values.
    runs <- vapply(seq_len(shared.rank.runs), function(b) {
        do.call(.fit_template_once, c(args, list(shared.rank = NULL,
            .warn = FALSE, .quiet = TRUE)))$shared.rank
        }, integer(1))
    counts <- table(factor(runs, levels = sort(unique(runs))))
    r <- as.integer(round(mean(runs)))
    agreement <- mean(runs == r)
    if (agreement < 0.5) {
        message("Shared rank ", r, " (mean of ", shared.rank.runs,
            " stability-selection runs; distribution ",
            paste(names(counts), counts, sep = ": ", collapse = ", "),
            "). The data do not separate neighbouring ranks clearly; for a ",
            "simulation study, consider fixing shared.rank.")
    }
    out <- do.call(.fit_template_once, c(args, list(shared.rank = NULL,
        .forced_rank = r, .warn = FALSE)))
    out$shared.rank.selection <- list(
        runs = runs,
        distribution = stats::setNames(as.integer(counts), names(counts)),
        agreement = agreement
    )
    out$provenance$shared_rank_runs <- shared.rank.runs
    out
}

.fit_template_once <- function(template, rank = NULL, var.explained = 0.7,
    max.rank = NULL, shared.rank = NULL, nperm = 200,
    view.transform = NULL, .forced_rank = NULL, .warn = TRUE, .quiet = FALSE) {

    rank_input <- rank
    template <- as_mimesys_view_list(template)
    raw_template <- template
    template <- .mimesys_apply_view_transforms(template,
        view.transform = view.transform)
    view_transform <- attr(template, "view.transform")
    transform_info <- attr(template, "transform.info")
    stopifnot(is.list(template), length(template) >= 1)
    if (is.null(names(template))) names(template) <- paste0("view",
        seq_along(template))
    views <- names(template)
    M <- length(views)
    template <- lapply(template, as.matrix)

    ## Match samples across views
    rn <- lapply(template, rownames)
    dup <- views[vapply(rn, function(r) !is.null(r) && anyDuplicated(r) > 0L, logical(1))]
    if (length(dup) > 0L) {
        stop("Duplicate sample IDs in view(s): ", paste(dup, collapse = ", "),
            ". Each sample must appear once per view.", call. = FALSE)
    }
    if (all(!vapply(rn, is.null, logical(1)))) {
        common <- Reduce(intersect, rn)
        if (length(common) < 10) {
            stop("Only ", length(common), " samples are shared across views (at least 10 ",
                "are needed). Samples are matched by sample ID: the column names of each ",
                "experiment or the row names of each matrix.", call. = FALSE)
        }
        template <- lapply(template, function(x) x[common, , drop = FALSE])
        raw_template <- lapply(raw_template, function(x) x[common, ,
            drop = FALSE])
    } else if (length(unique(vapply(template, nrow, integer(1)))) != 1) {
        stop("Views have different numbers of rows and no rownames to match on.")
    }
    n0 <- nrow(template[[1]])
    if (is.null(max.rank)) max.rank <- max(1, floor(n0 / 3))
    rank_2logp <- is.character(rank) && length(rank) == 1L && identical(rank,
        "2logp")
    if (is.character(rank) && !rank_2logp) {
        stop("rank must be NULL, numeric, or \"2logp\".", call. = FALSE)
    }
    if (!is.null(rank) && !rank_2logp) {
        rank <- .per_view_values(rank, views, "rank")
        if (any(!is.finite(rank)) || any(rank < 1)) {
            stop("rank values must be positive integers.", call. = FALSE)
        }
    }

    fits <- vector("list", M)
    observed <- vector("list", M)
    names(fits) <- views
    names(observed) <- views
    for (m in seq_len(M)) {
        X <- template[[m]]
        if (is.null(colnames(X))) colnames(X) <- paste0(views[m], "_",
            seq_len(ncol(X)))
        ## Drop all-missing features, warn about heavily imputed ones,
        ## median-impute the rest, and drop constant features
        miss <- colMeans(is.na(X))
        if (any(miss == 1)) {
            if (!.quiet) message("Dropping ", sum(miss == 1), " all-missing feature(s) from view ",
                views[m], ".")
            X <- X[, miss < 1, drop = FALSE]
            miss <- miss[miss < 1]
        }
        if (any(miss > 0.5) && !.quiet) {
            warning(sum(miss > 0.5), " feature(s) in view ", views[m], " are more than 50% ",
                "missing; median imputation puts most of their mass at one value. ",
                "Consider removing them before fitting.", call. = FALSE)
        }
        if (anyNA(X)) {
            X <- .impute_column_medians(X)
        }
        keep <- colSums(X != X[rep(1L, nrow(X)), , drop = FALSE]) > 0
        X <- X[, keep, drop = FALSE]
        if (ncol(X) < 1L) {
            stop("View ", views[m], " has no non-constant features after removing ",
                "missing and constant features.", call. = FALSE)
        }
        ## Original-scale quantile functions for transformed views, so simulated
        ## values can be returned on the measurement scale (back.transform).
        qmat_raw <- discrete_raw <- NULL
        if (!identical(view_transform[[views[m]]], "none")) {
            R0 <- as.matrix(raw_template[[m]])
            storage.mode(R0) <- "double"
            if (is.null(colnames(R0))) colnames(R0) <- colnames(template[[m]])
            R0 <- R0[, colnames(X), drop = FALSE]
            if (anyNA(R0)) {
                R0 <- .impute_column_medians(R0)
            }
            qmat_raw <- .sort_columns(R0)
            discrete_raw <- .integer_columns(R0)
        }
        ## Composition constraint (closed totals or bounded relative
        ## abundances), from the view on its original measurement scale
        composition <- .composition_info(if (is.null(qmat_raw)) X else R0)
        ## Retain the exact empirical matrix used to fit this view. Crosstalk
        ## calibration compares null simulations against this observed structure.
        observed[[m]] <- X

        ## (i) Marginals: sorted values define the empirical quantile function
        qmat <- .sort_columns(X)
        discrete <- .integer_columns(X)

        ## Normal scores with random tie breaking (randomized PIT, handles zeros)
        Z <- apply(X, 2, function(v) stats::qnorm(rank(v,
            ties.method = "random") / (n0 + 1)))
        Z <- scale(Z)
        two_part <- .fit_two_part_marginals(X)

        ## (ii) Within-view factor model via truncated SVD
        d  <- svd(Z, nu = 0, nv = 0)$d
        ve <- cumsum(d^2) / sum(d^2)
        k  <- if (rank_2logp) {
            ceiling(2 * log(ncol(X)))
        } else if (!is.null(rank)) {
            rank[m]
        } else {
            which(ve >= var.explained)[1]
        }
        k  <- as.integer(max(1, min(k, max.rank, ncol(X) - 1)))
        sv <- svd(Z, nu = k, nv = k)
        Lambda <- sv$v[, seq_len(k),
            drop = FALSE] %*% diag(sv$d[seq_len(k)] / sqrt(n0 - 1), k)
        h <- rowSums(Lambda^2)
        big <- h > 0.99
        if (any(big)) Lambda[big, ] <- Lambda[big, ,
            drop = FALSE] * sqrt(0.99 / h[big])
        rownames(Lambda) <- colnames(X)

        fits[[m]] <- list(qmat = qmat, discrete = discrete, Lambda = Lambda,
            qmat_raw = qmat_raw, discrete_raw = discrete_raw,
            composition = composition,
            two_part = two_part,
            scores = sv$u[, seq_len(k), drop = FALSE] * sqrt(n0 - 1),
            rank = k, var.explained = ve[k],
            raw.var = .col_vars(X))
    }

    ## (iii) Shared cross-view space by MAXVAR generalized CCA. Scores are
    ## orthonormal within view, so the top eigenvectors of the correlation matrix
    ## of the stacked scores define the shared variates.
    B <- NULL; eig_obs <- NULL; thr <- NULL; pval <- NULL; step_eig <- NULL
    r <- 0L
    if (M > 1) {
        Fall <- do.call(cbind, lapply(fits, `[[`, "scores"))
        eg <- eigen(stats::cor(Fall), symmetric = TRUE)
        eig_obs <- eg$values
        if (is.null(shared.rank)) {
            ## Sequential (deflated) parallel analysis: at each step, test the
            ## largest GCCA eigenvalue of the current scores against the
            ## permutation distribution of the largest eigenvalue, then
            ## remove that shared direction from every view and repeat
            pa <- .sequential_parallel_analysis(lapply(fits, `[[`, "scores"), nperm,
                kmax = min(vapply(fits, `[[`, integer(1), "rank")))
            r <- pa$rank
            thr <- pa$threshold
            pval <- pa$pvalue
            step_eig <- pa$step_eigen
            if (!is.null(.forced_rank)) {
                r <- as.integer(.forced_rank)
            } else if (.warn) {
                .warn_borderline_shared_rank(r, pval, nperm)
            }
        } else {
            if (!is.numeric(shared.rank) || length(shared.rank) != 1L ||
                !is.finite(shared.rank) || shared.rank < 0 ||
                    abs(shared.rank - round(shared.rank)) > 1e-8) {
                stop("shared.rank must be NULL or a single non-negative integer.",
                    call. = FALSE)
            }
            r <- as.integer(round(shared.rank))
            r_max <- min(vapply(fits, `[[`, integer(1), "rank"))
            if (r > r_max) {
                stop("shared.rank = ", r, " exceeds the smallest within-view rank (",
                    r_max, "); the shared space cannot have more dimensions than any ",
                    "view's factor model. Use shared.rank <= ", r_max,
                    " or increase rank.", call. = FALSE)
            }
        }
        r_max <- min(vapply(fits, `[[`, integer(1), "rank"))
        if (r > r_max) {
            if (!.quiet) message("Parallel analysis selected ", r, " shared factors; capping at the ",
                "smallest within-view rank (", r_max, ").")
            r <- r_max
        }
        if (r > 0) {
            S <- scale(Fall %*% eg$vectors[, seq_len(r), drop = FALSE])
            S <- qr.Q(qr(S)) * sqrt(n0 - 1)                      # sample-orthonormal
            B <- lapply(fits, function(f) crossprod(f$scores,
                S) / (n0 - 1))   # k_m x r
            names(B) <- views
        }
    }

    ## Per-sample factor scores are only needed to estimate the shared space; drop
    ## them so the fitted template stores as little sample-level data as possible.
    fits <- lapply(fits, function(f) {
        f$scores <- NULL
        f
    })

    provenance <- list(
        package_version = as.character(utils::packageVersion("MimeSys")),
        fitted_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
        rank_input = rank_input,
        var_explained = var.explained,
        max_rank = max.rank,
        shared_rank_input = shared.rank,
        nperm = nperm,
        view_transform = view_transform,
        input_dimensions = data.frame(
            view = views,
            samples = vapply(observed, nrow, integer(1)),
            features = vapply(observed, ncol, integer(1)),
            row.names = NULL
        )
    )
    structure(list(views = views, n0 = n0, fits = fits, observed = observed,
        B = B, shared.rank = r,
        gcca.eigen = eig_obs, gcca.threshold = thr, gcca.pvalue = pval,
        gcca.step.eigen = step_eig,
        view.transform = view_transform, transform.info = transform_info,
        provenance = provenance),
        class = "mimesys_template")
}

#' @export
print.mimesys_template <- function(x, ...) {
    cat("MimeSys template\n")
    cat("  Views: ", paste(x$views, collapse = ", "), "\n", sep = "")
    cat("  Matched samples: ", x$n0, "\n", sep = "")
    ranks <- vapply(x$fits, function(f) f$rank, integer(1))
    features <- vapply(x$fits, function(f) nrow(f$Lambda), integer(1))
    cat("  Features: ", paste(names(features), features, sep = "=",
        collapse = ", "), "\n", sep = "")
    cat("  Within-view ranks: ", paste(names(ranks), ranks, sep = "=",
        collapse = ", "), "\n", sep = "")
    cat("  Shared rank: ", x$shared.rank, "\n", sep = "")
    cat("  Calibration-ready: ", !is.null(x$observed), "\n", sep = "")
    if (!is.null(x$provenance$package_version)) {
        cat("  MimeSys version: ", x$provenance$package_version, "\n", sep = "")
    }
    invisible(x)
}


#' Simulate null multiview data from a fitted template
#'
#' @param template A \code{mimesys_template} from \code{fit_template()}.
#' @param n Number of samples to generate (need not match the template).
#' @param p Feature counts per view. Use \code{NULL} to keep template sizes, a
#'   single integer to use the same count for every view, an unnamed vector in
#'   template view order, or a named vector with one value per view. Fewer
#'   features than the template subsamples features; more adds offspring
#'   features
#'   (see \code{feature.jitter}).
#' @param crosstalk Crosstalk weights gamma, scalar or one per view. 0 gives
#'   independent views, 1 the template level, values above 1 amplify it up to
#'   \code{crosstalk_max()}. Unequal weights let some views carry more of the
#'   shared structure than others (the analogue of t_x != t_z in Ding et al.
#'   2022).
#' @param feature.select Feature choice when p is below template size:
#'   "random" or "variance" (highest variance first, on the fitted scale, i.e.
#'   after any \code{view.transform}).
#' @param feature.jitter Between 0 and 1. An offspring feature copies a parent's
#'   marginal and gets loadings sqrt(1 - jitter) * parent + sqrt(jitter) *
#'   another
#'   feature, so it resembles but does not duplicate its parent.
#' @param view.models Optional named vector of per-view observation models. Use
#'   `"two_part_copula"` for sparse microbiome views: presence/absence is drawn
#'   from each feature's zero probability and, given presence, abundance from
#'   that feature's positive values. Other views default to `"continuous"`.
#'   The two-part model is for sparse non-negative abundances only: views with
#'   negative values, or fitted with `view.transform = "clr"`, give an error.
#'   See Details for the composition constraint.
#' @param back.transform Logical; if \code{TRUE}, views fitted with a
#'   \code{view.transform} are returned on the original measurement scale. Each
#'   simulated feature's latent uniform is mapped through that feature's
#'   original-scale empirical quantile function, which is the exact inverse for
#'   \code{log1p} and, for \code{clr}, reproduces each feature's original
#'   marginal (including zeros) with dependence from the CLR-scale model.
#'   The default \code{FALSE} returns the fitted (transformed) scale, which is what
#'   calibration and validation compare against.
#' @param seed Random seed (default `1234`, the same default as every other
#'   MimeSys simulator). The caller's random-number state is restored on exit.
#'   Use `NULL` to draw from the current random-number stream.
#' @details Composition constraint. For a view whose original-scale data are a
#'   composition, simulated samples never exceed the template total. If every
#'   template sample has the same total (a closed composition, e.g. relative
#'   abundances summing to 1) and at least the template's features are
#'   simulated, including offspring features, every simulated sample is
#'   renormalized to that total; offspring features therefore dilute the
#'   others. If fewer features are simulated, samples above the total are
#'   scaled down to it. If the template is itself a subset of a composition
#'   (values in \[0, 1\] with sample totals at most 1), samples above the
#'   largest observed total are scaled down to it. The constraint applies to
#'   original-scale output (untransformed views, or transformed views with
#'   \code{back.transform = TRUE}).
#' @return A pcl-style list (feature_table features x samples, sample_metadata,
#'   feature_metadata) plus \code{truth}: factor scores, shared variates,
#'   crosstalk
#'   weights and implied latent canonical correlations. feature_metadata has
#'   \code{shared_loading}, the correlation of each feature's normal score with
#'   shared variate 1 at crosstalk = 1 (used for spike-in alignment).
#' @examples
#' data(stelzer_mae)
#' tmpl <- fit_template(stelzer_mae, shared.rank = 2)
#' null <- sim_from_template(tmpl, n = 30, p = 10, crosstalk = 0.5)
#' dim(null$feature_table)
#' head(null$feature_metadata)
#' @export
sim_from_template <- function(template, n, p = NULL, crosstalk = 1,
    feature.select = c("random", "variance"),
    feature.jitter = 0.2,
    view.models = NULL,
    back.transform = FALSE,
    seed = 1234) {
    .check_count(n, "n", min = 1L)
    if (!is.null(seed)) withr::local_seed(seed)
    .sim_from_template_impl(template, n = n, p = p, crosstalk = crosstalk,
        feature.select = feature.select, feature.jitter = feature.jitter,
        view.models = view.models, back.transform = back.transform)
}

.sim_from_template_impl <- function(tmpl, n, p = NULL, crosstalk = 1,
    feature.select = c("random", "variance"),
    feature.jitter = 0.2,
    view.models = NULL,
    back.transform = FALSE,
    .structure = NULL,
    .latent_hook = NULL,
    .observed_hook = NULL) {

    ## .structure (internal): reuse the parent features and loadings from an
    ## earlier call (attr(x, "structure")), so counterfactual simulations share
    ## the same features. .latent_hook (internal): function(Z_list) applied to
    ## the list of every view's latent normal scores before they are mapped to
    ## the data scale (used to inject exposure effects for mediation
    ## simulation).
    stopifnot(inherits(tmpl, "mimesys_template"))
    feature.select <- match.arg(feature.select)
    .check_count(n, "n", min = 1L)
    .check_fraction(feature.jitter, "feature.jitter")
    views <- tmpl$views
    M <- length(views)
    p0 <- vapply(tmpl$fits, function(f) nrow(f$Lambda), integer(1))
    p <- .normalize_feature_counts(p, views, p0, arg_name = "p")
    view.models <- .normalize_view_models(view.models, views)
    vt <- tmpl$view.transform
    tp_views <- views[view.models == "two_part_copula"]
    if (length(tp_views) > 0L) {
        ## The two-part model is for sparse non-negative (e.g. microbiome)
        ## abundances: zeros mark absence, positive values are abundances.
        bad <- if (is.null(vt)) character(0) else tp_views[vt[tp_views] %in% "clr"]
        if (length(bad) > 0L) {
            stop("two_part_copula is meant for sparse non-negative data, but view(s) ",
                paste(bad, collapse = ", "), " were fitted on the CLR scale, which has no zeros. ",
                "Fit those views without view.transform (or with log1p) to use the two-part model.",
                call. = FALSE)
        }
        neg <- tp_views[vapply(tp_views, function(v) {
            any(tmpl$fits[[v]]$qmat < 0, na.rm = TRUE)
        }, logical(1))]
        if (length(neg) > 0L) {
            stop("two_part_copula is meant for sparse non-negative (e.g. microbiome) data, ",
                "but view(s) ", paste(neg, collapse = ", "), " contain negative values. ",
                "Use the default \"continuous\" model for those views.", call. = FALSE)
        }
        nozero <- tp_views[vapply(tp_views, function(v) {
            all(tmpl$fits[[v]]$two_part$zero_prob == 0)
        }, logical(1))]
        if (length(nozero) > 0L) {
            warning("two_part_copula is meant for sparse data, but view(s) ",
                paste(nozero, collapse = ", "), " have no zeros; the two-part model ",
                "then reduces to the continuous model.", call. = FALSE)
        }
    }

    ks <- vapply(tmpl$fits, function(f) f$rank, integer(1))
    block <- rep(seq_len(M), ks)
    r <- tmpl$shared.rank
    gamma <- .per_view_values(crosstalk, views, "crosstalk")
    if (!is.numeric(gamma) || any(!is.finite(gamma))) stop("crosstalk must be finite numbers.", call. = FALSE)
    if (any(gamma < 0)) stop("crosstalk must be non-negative.")
    lim <- crosstalk_max(tmpl, per.view = TRUE)
    if (any(gamma > lim + 1e-10)) {
        stop("crosstalk exceeds the valid limit for: ",
            paste(sprintf("%s (max %.2f)", views[gamma > lim],
                lim[gamma > lim]), collapse = ", "))
    }
    sample_ids <- paste0("Subject", seq_len(n))

    X_list <- meta_list <- Xorig_list <- vector("list", M)
    parent_list <- lam_list <- offspring_list <- vector("list", M)
    if (!is.null(.structure)) {
        parent_list <- .structure$parent
        lam_list <- .structure$lam
        offspring_list <- .structure$offspring
    }
    for (m in seq_len(M)) {
        if (!is.null(.structure)) next
        fit <- tmpl$fits[[m]]
        ## Template parents for the requested number of features
        if (p[m] <= p0[m]) {
            parent <- if (feature.select == "random") sort(sample.int(p0[m],
                p[m])) else
            sort(order(fit$raw.var, decreasing = TRUE)[seq_len(p[m])])
            Lam <- fit$Lambda[parent, , drop = FALSE]
            offspring <- rep(FALSE, p[m])
        } else {
            extra  <- sample.int(p0[m], p[m] - p0[m], replace = TRUE)
            other  <- sample.int(p0[m], p[m] - p0[m], replace = TRUE)
            parent <- c(seq_len(p0[m]), extra)
            Lam <- rbind(fit$Lambda,
                sqrt(1 - feature.jitter) * fit$Lambda[extra, , drop = FALSE] +
                    sqrt(feature.jitter) * fit$Lambda[other, , drop = FALSE])
            h <- rowSums(Lam^2)
            Lam[h > 0.99, ] <- Lam[h > 0.99, ,
                drop = FALSE] * sqrt(0.99 / h[h > 0.99])
            offspring <- c(rep(FALSE, p0[m]), rep(TRUE, length(extra)))
        }
        parent_list[[m]] <- parent
        lam_list[[m]] <- Lam
        offspring_list[[m]] <- offspring
    }

    ## Latent factor scores: shared part plus view-specific remainder
    shared_out <- matrix(stats::rnorm(n * max(r, 1)), n, max(r, 1))[,
        seq_len(r), drop = FALSE]
    scores_out <- do.call(cbind, lapply(seq_len(M), function(m) {
        E <- matrix(stats::rnorm(n * ks[m]), n, ks[m])
        if (r == 0 || gamma[m] == 0) return(E)
        Bm <- sqrt(gamma[m]) * tmpl$B[[m]]
        shared_out %*% t(Bm) + E %*% .sqrt_psd(diag(ks[m]) - Bm %*% t(Bm))
    }))

    ## Latent normal scores of every view first, so that a latent hook can act
    ## on all views jointly (sequential mediation links features across
    ## views). The random draws keep the order of the per-view loop.
    ## (The composition step's random picks for empty samples are drawn here
    ## too, in their original order.)
    Z_list <- abundance_noise <- picks_fitted <- picks_orig <- vector("list", M)
    for (m in seq_len(M)) {
        Lam <- lam_list[[m]]
        psi <- 1 - rowSums(Lam^2)
        Z_list[[m]] <- scores_out[, block == m, drop = FALSE] %*% t(Lam) +
            sweep(matrix(stats::rnorm(n * p[m]), n, p[m]), 2, sqrt(psi), "*")
        if (identical(view.models[[m]], "two_part_copula")) {
            abundance_noise[[m]] <- matrix(stats::rnorm(n * p[m], sd = 0.15), n, p[m])
        }
        comp <- tmpl$fits[[m]]$composition
        closed <- !is.null(comp) && identical(comp$type, "closed") &&
            is.finite(comp$total) && p[m] >= p0[m]
        transformed_m <- !is.null(tmpl$view.transform) &&
            !identical(tmpl$view.transform[[views[m]]], "none")
        if (closed && !transformed_m) {
            picks_fitted[[m]] <- sample.int(p[m], n, replace = TRUE)
        }
        if (closed && transformed_m &&
            (isTRUE(back.transform) || identical(back.transform, "both"))) {
            picks_orig[[m]] <- sample.int(p[m], n, replace = TRUE)
        }
    }
    ## Context for the latent hooks: loadings, cross-view coupling and the
    ## latent scores drawn above.
    hook_context <- list(lam = lam_list, B = tmpl$B, gamma = gamma, ks = ks,
        r = r, p = p, views = views)
    .call_hook <- function(hook, Z) {
        if (length(formals(hook)) >= 2L) hook(Z, hook_context) else hook(Z)
    }
    if (is.function(.latent_hook)) Z_list <- .call_hook(.latent_hook, Z_list)
    ## .observed_hook (internal): a second set of latent scores (e.g. with
    ## batch effects) mapped through the same marginals. The base scores still
    ## give feature_table; the second set gives observed_table.
    Z_observed <- if (is.function(.observed_hook)) .call_hook(.observed_hook, Z_list) else NULL

    ## Map latent scores to the data scale. No random numbers are drawn here,
    ## so the mapping can be applied to more than one set of latent scores.
    .map_views <- function(Z_all) {
        X_out <- Xorig_out <- scale_out <- vector("list", M)
        for (m in seq_len(M)) {
            fit <- tmpl$fits[[m]]
            parent <- parent_list[[m]]

            U <- NULL
            ## Map the normal scores through each parent's empirical quantiles
            Z <- Z_all[[m]]
            if (identical(view.models[[m]], "two_part_copula")) {
                ## A little extra noise decouples abundance from presence;
                ## rescaling keeps unit variance so pnorm() stays uniform.
                Z_abundance <- (Z + abundance_noise[[m]]) / sqrt(1 + 0.15^2)
                Xm <- .simulate_two_part_view(
                    fit = fit,
                    parent = parent,
                    U_presence = stats::pnorm(Z),
                    U_abundance = stats::pnorm(Z_abundance)
                )
            } else {
                U <- stats::pnorm(Z)
                Xm <- .emp_quantile_matrix(fit$qmat, parent, U, fit$discrete)
            }

            transformed <- !is.null(tmpl$view.transform) &&
                !identical(tmpl$view.transform[[views[m]]], "none")
            if (!transformed) {
                ## Fitted scale is the original scale: enforce the composition
                ## constraint (no sample above the template total)
                Xm <- .apply_composition(Xm, fit$composition, p[m], p0[m],
                    pick = picks_fitted[[m]])
            }
            view_scale <- "fitted"
            if (isTRUE(back.transform) || identical(back.transform, "both")) {
                Xo <- .mimesys_original_scale_view(fit, parent, Xm, U,
                    info = tmpl$transform.info[[views[m]]],
                    two_part = identical(view.models[[m]], "two_part_copula"),
                    full_features = p[m] >= p0[m], view = views[m])
                if (transformed) {
                    Xo <- .apply_composition(Xo, fit$composition, p[m], p0[m],
                        pick = picks_orig[[m]])
                }
                if (identical(back.transform, "both")) {
                    dimnames(Xo) <- list(sample_ids, paste0(views[m], "_",
                        seq_len(p[m])))
                    Xorig_out[[m]] <- Xo
                } else {
                    Xm <- Xo
                    view_scale <- "original"
                }
            }
            ids <- paste0(views[m], "_", seq_len(p[m]))
            dimnames(Xm) <- list(sample_ids, ids)
            X_out[[m]] <- Xm
            scale_out[[m]] <- view_scale
        }
        list(X = X_out, Xorig = Xorig_out, scale = scale_out)
    }
    mapped <- .map_views(Z_list)
    X_list <- mapped$X
    Xorig_list <- mapped$Xorig
    for (m in seq_len(M)) {
        fit <- tmpl$fits[[m]]
        ids <- paste0(views[m], "_", seq_len(p[m]))
        meta_list[[m]] <- data.frame(
            featureID = ids, featureType = views[m],
            template_feature = colnames(fit$qmat)[parent_list[[m]]],
            offspring = offspring_list[[m]],
            view_model = view.models[[m]],
            scale = mapped$scale[[m]],
            transform = if (is.null(tmpl$view.transform)) "none" else tmpl$view.transform[[views[m]]],
            shared_loading = if (r > 0) as.numeric(lam_list[[m]] %*% tmpl$B[[m]][,
                1]) else 0,
            stringsAsFactors = FALSE
        )
    }

    feature_metadata <- do.call(rbind, meta_list)
    rownames(feature_metadata) <- feature_metadata$featureID

    out <- list(feature_table = as.data.frame(t(do.call(cbind, X_list))),
        sample_metadata = data.frame(subjectID = sample_ids,
            row.names = sample_ids),
        feature_metadata = feature_metadata,
        truth = list(scores = scores_out, shared = shared_out,
            crosstalk = gamma,
            cancor = .realized_cancor(tmpl, gamma),
            view.models = stats::setNames(view.models, views)))
    if (identical(back.transform, "both")) {
        out$original_table <- as.data.frame(t(do.call(cbind, Xorig_list)))
    }
    if (!is.null(Z_observed)) {
        observed <- .map_views(Z_observed)
        out$observed_table <- as.data.frame(t(do.call(cbind, observed$X)))
        if (identical(back.transform, "both")) {
            out$observed_original_table <- as.data.frame(t(do.call(cbind,
                observed$Xorig)))
        }
    }
    attr(out, "structure") <- list(parent = parent_list, lam = lam_list,
        offspring = offspring_list)
    out
}

## Composition type of a view on its original scale:
##   "closed":     non-negative, every sample has the same total (e.g. 1 or 100)
##   "proportion": values in [0, 1] with sample totals at most 1 (relative
##                 abundances of a subset of taxa); `total` is the largest
##                 observed sample total
##   "none":       anything else (no constraint)
.composition_info <- function(x) {
    x <- as.matrix(x)
    fin <- x[is.finite(x)]
    if (length(fin) == 0L || any(fin < 0)) return(list(type = "none", total = NA_real_))
    rs <- rowSums(x, na.rm = TRUE)
    if (mean(rs) > 0 && stats::sd(rs) <= 1e-6 * mean(rs)) {
        return(list(type = "closed", total = mean(rs)))
    }
    if (all(fin <= 1) && max(rs) <= 1 + 1e-8 && ncol(x) > 1L) {
        return(list(type = "proportion", total = max(rs)))
    }
    list(type = "none", total = NA_real_)
}

## Keep simulated samples within the template's composition total.
## Closed compositions: when at least the template's features are simulated
## (including offspring), every sample is renormalized to the total; when a
## subset is simulated, samples above the total are scaled down to it.
## Proportions (subset of a composition): samples above the largest observed
## total are scaled down to it.
.apply_composition <- function(x, comp, p, p0, pick = NULL) {
    if (is.null(comp) || identical(comp$type, "none") || !is.finite(comp$total)) {
        return(x)
    }
    total <- comp$total
    rs <- rowSums(x)
    if (identical(comp$type, "closed") && p >= p0) {
        ok <- rs > 0
        x[ok, ] <- x[ok, , drop = FALSE] * (total / rs[ok])
        ## An all-zero sample puts its total on one random feature. The draws
        ## are made for every sample so the random stream does not depend on
        ## which samples are empty (keeps paired simulations in step).
        if (is.null(pick)) pick <- sample.int(ncol(x), nrow(x), replace = TRUE)
        for (i in which(!ok)) x[i, pick[i]] <- total
        return(x)
    }
    over <- rs > total
    if (any(over)) {
        x[over, ] <- x[over, , drop = FALSE] * (total / rs[over])
    }
    x
}

## Map one simulated view to the original measurement scale. Continuous views
## reuse the simulated uniforms with the original-scale quantile functions
## (exact for monotone transforms such as log1p; for CLR this keeps each
## feature's original marginal, including its zeros, with dependence from the
## CLR-scale model). Two-part views use the exact log1p inverse.
.mimesys_original_scale_view <- function(fit, parent, Xm, U, info, two_part,
    full_features, view) {
    if (is.null(info) || identical(info$type, "none")) {
        return(Xm)
    }
    if (two_part) {
        if (identical(info$type, "log1p")) {
            return(.mimesys_inverse_log1p(Xm, info))
        }
        warning("Cannot back-transform two-part view ", view, " fitted on the ",
            info$type,
            " scale; returning the fitted scale.", call. = FALSE)
        return(Xm)
    }
    if (is.null(fit$qmat_raw) || is.null(U)) {
        return(Xm)
    }
    .emp_quantile_matrix(fit$qmat_raw, parent, U, fit$discrete_raw)
}


#' Largest valid crosstalk weight
#'
#' The construction needs \eqn{\gamma_m B_m B_m^T \le I}, so the limit for view
#' m
#' is \eqn{1 / \sigma_{max}(B_m)^2}, which is always at least 1.
#' @param template A \code{mimesys_template}.
#' @param per.view Return one limit per view instead of the overall scalar
#'   limit.
#' @return A single number (the smallest per-view limit) or, if
#'   \code{per.view = TRUE}, a named numeric vector with one limit per view.
#' @examples
#' data(stelzer_mae)
#' template <- fit_template(stelzer_mae, shared.rank = 2)
#' crosstalk_max(template)
#' crosstalk_max(template, per.view = TRUE)
#' @export
crosstalk_max <- function(template, per.view = FALSE) {
    if (is.null(template$B)) {
        lim <- setNames(rep(Inf, length(template$views)), template$views)
    } else {
        lim <- vapply(template$B, function(b) 1 / max(svd(b, nu = 0, nv = 0)$d)^2,
            numeric(1))
    }
    if (per.view) lim else min(lim)
}


## Sequential parallel analysis for the number of shared variates. Comparing
## the i-th observed eigenvalue with the i-th null eigenvalue over-selects once
## real shared structure exists (the remaining eigenvalues behave like the
## *largest* null eigenvalues), so each step compares the current largest
## eigenvalue with the null distribution of the largest eigenvalue, after the
## shared directions found so far have been projected out of every view.
.sequential_parallel_analysis <- function(scores, nperm, kmax, alpha = 0.05) {
    n <- nrow(scores[[1]])
    top_eigen <- function(sc) {
        eigen(stats::cor(do.call(cbind, sc)), symmetric = TRUE,
            only.values = TRUE)$values[1]
    }
    pvalue <- threshold <- step_eigen <- numeric(0)
    r <- 0L
    current <- scores
    while (r < kmax && all(vapply(current, ncol, integer(1)) > 0L)) {
        obs <- top_eigen(current)
        null <- replicate(nperm, top_eigen(lapply(current,
            function(f) f[sample.int(n), , drop = FALSE])))
        p <- (1 + sum(null >= obs)) / (nperm + 1)
        pvalue <- c(pvalue, p)
        threshold <- c(threshold, stats::quantile(null, 1 - alpha, names = FALSE))
        step_eigen <- c(step_eigen, obs)
        if (p >= alpha) break
        r <- r + 1L
        ## Deflate: remove the leading shared variate from each view and
        ## re-orthonormalize the remaining scores
        all_sc <- do.call(cbind, current)
        v <- eigen(stats::cor(all_sc), symmetric = TRUE)$vectors[, 1]
        s <- scale(all_sc %*% v)
        current <- lapply(current, function(f) {
            f <- f - s %*% (crossprod(s, f) / sum(s^2))
            k <- ncol(f) - 1L
            if (k < 1L) return(f[, 0, drop = FALSE])
            sv <- svd(f, nu = k, nv = 0)
            sv$u[, seq_len(k), drop = FALSE] * sqrt(n - 1)
        })
    }
    list(rank = r, pvalue = pvalue, threshold = threshold, step_eigen = step_eigen)
}

## Warn when the parallel-analysis choice of shared rank rests on a borderline
## factor: the last retained factor is only just significant, or the first
## rejected factor is only just non-significant.
.warn_borderline_shared_rank <- function(r, pval, nperm) {
    k <- length(pval)
    last_in <- if (r > 0L) pval[[r]] else NA_real_
    first_out <- if (r < k) pval[[r + 1L]] else NA_real_
    borderline <- (is.finite(last_in) && last_in > 0.025) ||
        (is.finite(first_out) && first_out < 0.08)
    if (borderline) {
        warning(
            "Shared rank ", r, " was chosen by parallel analysis, but the choice is ",
            "borderline (permutation p-value ",
            if (is.finite(last_in)) sprintf("%.3f for factor %d", last_in, r) else "",
            if (is.finite(last_in) && is.finite(first_out)) " and " else "",
            if (is.finite(first_out)) sprintf("%.3f for factor %d", first_out, r + 1L) else "",
            "). The rank may change with the seed; consider setting shared.rank ",
            "explicitly, a seed, or a larger nperm (currently ", nperm, ").",
            call. = FALSE
        )
    }
    invisible(NULL)
}

## ---- internal helpers ---------------------------------------------------

.normalize_view_models <- function(view.models, views) {
    valid <- c("continuous", "two_part_copula")
    if (is.null(view.models)) {
        return(rep("continuous", length(views)))
    }
    if (is.null(names(view.models)) || any(!nzchar(names(view.models)))) {
        if (length(view.models) != 1L) {
            stop("view.models must be named by view, or be a single shared model.",
                call. = FALSE)
        }
        view.models <- stats::setNames(rep(view.models[[1]], length(views)),
            views)
    }
    out <- vapply(views, function(view) {
        model <- if (view %in% names(view.models)) view.models[[view]] else "continuous"
        model <- tolower(gsub("-", "_", model))
        if (model %in% c("none", "default", "standard")) {
            model <- "continuous"
        }
        if (!model %in% valid) {
            stop("Unknown view model '", model, "' for view ", view, ".",
                call. = FALSE)
        }
        model
        }, character(1))
    unname(out)
}

.normalize_feature_counts <- function(p, views, default, arg_name = "p") {
    if (is.null(p)) {
        p <- default
        names(p) <- views
        return(as.integer(p))
    }
    if (is.list(p) && !is.data.frame(p)) {
        p <- unlist(p, use.names = TRUE)
    }
    if (!is.numeric(p) || length(p) == 0L) {
        stop(arg_name,
            " must be NULL, a numeric vector, or a list coercible to a numeric vector.", call. = FALSE)
    }
    if (length(p) == 1L && (is.null(names(p)) || !nzchar(names(p)))) {
        p <- stats::setNames(rep(p[[1]], length(views)), views)
    } else if (is.null(names(p)) || all(!nzchar(names(p)))) {
        if (length(p) != length(views)) {
            stop(
                arg_name, " must be length 1, length ", length(views),
                ", or named with the template view names: ", paste(views,
                    collapse = ", "),
                call. = FALSE
            )
        }
        names(p) <- views
    }
    missing <- setdiff(views, names(p))
    if (length(missing) > 0L) {
        stop(arg_name, " is missing view(s): ", paste(missing, collapse = ", "),
            call. = FALSE)
    }
    unknown <- setdiff(names(p), views)
    if (length(unknown) > 0L) {
        stop(arg_name, " has unknown view name(s): ", paste(unknown, collapse = ", "),
            ". Views are: ", paste(views, collapse = ", "), ".", call. = FALSE)
    }
    p <- p[views]
    if (any(!is.finite(p)) || any(p <= 0) || any(abs(p - round(p)) > 1e-8)) {
        stop(arg_name, " values must be positive finite integers.",
            call. = FALSE)
    }
    as.integer(p)
}

.fit_two_part_marginals <- function(X) {
    positive <- X > 0
    zero_prob <- colMeans(!positive)
    prevalence <- colMeans(positive)
    positive_values <- lapply(seq_len(ncol(X)), function(j) {
        vals <- sort(X[positive[, j], j])
        if (length(vals) == 0L) {
            vals <- 0
        }
        vals
    })
    positive_discrete <- vapply(positive_values, function(v) {
        all(abs(v - round(v)) < 1e-8)
        }, logical(1))
    list(
        zero_prob = zero_prob,
        prevalence = prevalence,
        positive_values = positive_values,
        positive_discrete = positive_discrete
    )
}

.simulate_two_part_view <- function(fit, parent, U_presence, U_abundance) {
    tp <- fit$two_part
    p <- length(parent)
    n <- nrow(U_presence)
    zp <- matrix(tp$zero_prob[parent], n, p, byrow = TRUE)
    present <- U_presence > zp
    ## Condition on presence: rescale the abundance uniform to the part of
    ## (0, 1) above the zero probability, so positive values are drawn from the
    ## whole positive distribution rather than only its upper quantiles.
    u <- ifelse(zp < 1, (U_abundance - zp) / pmax(1 - zp, 1e-12), 0.5)
    u <- pmin(pmax(u, 1e-8), 1 - 1e-8)
    Xm <- matrix(vapply(seq_len(p), function(j) {
        pj <- parent[[j]]
        .emp_quantile(tp$positive_values[[pj]], u[, j],
            tp$positive_discrete[[pj]])
        }, numeric(n)), n, p)
    Xm[!present] <- 0
    Xm
}

.sqrt_psd <- function(A) {
    e <- eigen(A, symmetric = TRUE)
    e$vectors %*% diag(sqrt(pmax(e$values, 0)),
        length(e$values)) %*% t(e$vectors)
}

## Latent canonical correlations between each pair of views implied by gamma
.realized_cancor <- function(tmpl, gamma) {
    views <- tmpl$views
    M <- length(views)
    gamma <- .per_view_values(gamma, views, "crosstalk")
    out <- list()
    if (is.null(tmpl$B) || M < 2) return(out)
    for (a in seq_len(M - 1L)) for (b in seq.int(a + 1L, M)) {
        C <- sqrt(gamma[a] * gamma[b]) * tmpl$B[[a]] %*% t(tmpl$B[[b]])
        out[[paste(views[a], views[b], sep = ":")]] <- svd(C, nu = 0,
            nv = 0)$d[seq_len(tmpl$shared.rank)]
    }
    out
}

## Vectorized column helpers
.col_vars <- function(x) {
    n <- nrow(x)
    centered <- x - rep(colMeans(x), each = n)
    colSums(centered^2) / (n - 1)
}

.sort_columns <- function(x) {
    out <- matrix(x[order(col(x), x)], nrow(x), ncol(x))
    dimnames(out) <- list(NULL, colnames(x))
    out
}

.integer_columns <- function(x) {
    colSums(abs(x - round(x)) >= 1e-8) == 0
}

.impute_column_medians <- function(x) {
    for (j in which(colSums(is.na(x)) > 0)) {
        x[is.na(x[, j]), j] <- stats::median(x[, j], na.rm = TRUE)
    }
    x
}

## Empirical quantiles for many features at once. `qmat` holds each template
## feature's sorted values in a column; column j of `U` is mapped through
## template feature `cols[j]`. Same rule as .emp_quantile().
.emp_quantile_matrix <- function(qmat, cols, U, discrete) {
    ## A column loop is faster here than whole-matrix indexing: each step is
    ## already vectorized over samples and avoids large index matrices.
    n <- nrow(U)
    matrix(vapply(seq_along(cols), function(j) {
        .emp_quantile(qmat[, cols[j]], U[, j], discrete[cols[j]])
        }, numeric(n)), n, length(cols))
}


## Empirical quantile function: step (type 1) for integer-valued features so
## counts stay integers, linear interpolation (type 7) otherwise
.emp_quantile <- function(sorted, u, discrete) {
    n0 <- length(sorted)
    if (discrete) return(sorted[pmin(pmax(ceiling(u * n0), 1), n0)])
    h  <- (n0 - 1) * u + 1
    lo <- floor(h)
    hi <- pmin(lo + 1, n0)
    sorted[lo] + (h - lo) * (sorted[hi] - sorted[lo])
}
