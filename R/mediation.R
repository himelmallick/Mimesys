#' Simulate multi-view mediation data from a template
#'
#' Simulates an exposure, a set of mediating features spread over one or more
#' views, and an outcome, on top of realistic null multi-view data from a
#' fitted template. The design is domain-agnostic: any view and any feature
#' can act as a mediator.
#'
#' The model has three parts.
#'
#' 1. **Exposure.** `T` is binary (`exposure.type = "binary"`, with
#'    `P(T = 1) = exposure.prevalence`) or continuous (standard normal). Its
#'    standardized version is `T* = (T - mu_T) / sd_T`.
#' 2. **Exposure to mediators (a-paths).** For each mediator `j`, the latent
#'    normal score from the template is replaced by
#'    `Z_j <- a_j T* + sqrt(1 - a_j^2) Z_j`, before it is mapped through the
#'    feature's empirical quantile function. For a binary exposure this sum is
#'    a two-component normal mixture, so it is passed through the mixture CDF
#'    and back to the normal scale. Either way the latent score stays exactly
#'    standard normal, so every mediator keeps its template distribution
#'    (including zeros and bounds); the exposure shifts where each sample falls
#'    in it. `a_j` is the mixing weight of the exposure: for a continuous
#'    exposure it equals the latent correlation between exposure and
#'    mediator; for a binary exposure the mixture correction lowers that
#'    correlation somewhat (the realized effect is reported as `a_effect`).
#' 3. **Mediators to outcome (b-paths).** With mediators standardized by their
#'    template mean and SD (`M~_j`),
#'    `Y = DE * T + sum_j b_j M~_j + sum_k g_k M~_k1 M~_k2 + sigma * e`, where the
#'    optional interaction terms use pairs of mediators (possibly from
#'    different views).
#'
#' Path sizes are set from the requested total indirect effect: `a_j = b_j`
#' in absolute value, with signs from `dir.a` and `dir.b`. The b-paths are then
#' rescaled so that the population natural indirect effect equals
#' `TIE * dir.a * dir.b` exactly, computed by Monte Carlo with `n.mc` paired
#' counterfactual draws per replicate. With a `"random"` direction (or, for
#' sequential mediation, negative links) the
#' pathway-specific effects can cancel, so their summed absolute size is set
#' to `TIE` instead and the net NIE is reported as truth. The true natural indirect effect (NIE), natural direct effect
#' (NDE = `DE`) and total effect (TE), the indirect effect carried by each view,
#' and each mediator's a-path, b-path and indirect contribution are returned as
#' truth. Effects are for the contrast `T = 1` versus `T = 0` (one SD for a
#' continuous exposure), in outcome units. For a continuous exposure the
#' mediators respond non-linearly to the exposure, so the truth also includes
#' `NIE_linear`, the indirect effect implied by each mediator's population
#' linear regression on the exposure (`a_slope`), which is what
#' product-of-coefficients estimators target; for a binary exposure the two
#' coincide.
#'
#' For a binary outcome, `NIE`, `NDE` and `TE` (and `TIE`, `DE`) are
#' conditional effects on the linear-predictor (log-odds) scale. Because odds
#' ratios are non-collapsible, the natural (marginal counterfactual) effects
#' differ; they are given in `truth$mediation$binary` as risk differences,
#' odds ratios and log odds ratios.
#'
#' In views that are closed compositions, renormalizing each sample after the
#' mediators are shifted changes the other features of the view too. The
#' column `exposure_effect` gives the true exposure effect on every feature
#' (standardized units, fitted scale); it is 0 for non-mediators except
#' through such closure.
#'
#' All truth refers to the fitted (transformed) scale; see `back.transform`.
#'
#' **Sequential mediation.** With `stages`, the mediators are grouped into
#' ordered stages of views, for example `list("Microbiome", "Metabolome")`
#' for exposure -> microbes -> metabolites -> outcome. Each mediator of stage
#' `s > 1` is linked to mediators of stage `s - 1` (its parents): its latent
#' score becomes `W_j = a_j T* + sum_k d_jk W_k + c_j Z_j`, where `W_k` are
#' the parents' latent scores (after their own exposure effect), `Z_j` is the
#' mediator's template latent score and `c_j` keeps the variance exactly 1
#' (allowing for the template correlation between the views), so every
#' mediator still keeps its template distribution. Links can be many-to-many:
#' one parent may feed several children and one child may have several
#' parents. `stage.links` sets the pattern, the share of each child's latent
#' variance explained by its parents, and the link signs; `stage.direct` scales
#' the direct exposure paths of each stage (0 for a stage reached only
#' through the previous one). The links are defined on the latent scale, so
#' in a closed composition the child depends on its parent before the
#' renormalization of the parent's view. For skewed or sparse parents (for
#' example microbiome relative abundances) the observed values are a
#' non-linear function of the latent score, so linear estimators on the
#' observed scale attribute part of the chain to the direct path.
#'
#' The natural indirect effect is split into natural path-specific effects,
#' one for every path `T -> stage r -> ... -> stage s -> Y`: for two stages,
#' `T -> M1 -> Y`, `T -> M1 -> M2 -> Y` and `T -> M2 -> Y`. They are computed
#' by the paired Monte Carlo with cross-world counterfactuals, switching the
#' paths on one at a time (later-starting paths first: `T -> M2 -> Y`, then
#' `T -> M1 -> M2 -> Y`, then `T -> M1 -> Y`), so they add up to the NIE. By
#' default all b-paths share one calibration scale, as without stages;
#' `stage.share` instead fixes the share of the indirect effect carried by
#' the mediators of each stage (exact without interactions). For example,
#' `stage.direct = c(1, 0)` with `stage.share = c(0, 1)` gives a pure chain
#' `T -> M1 -> M2 -> Y`, and `stage.links = list(strength = 0)` gives two
#' parallel mediator sets whose only association is the template crosstalk,
#' a null for the link between them.
#'
#' `effect.type` gives null scenarios for type I error studies: `"null_a"`
#' (mediators affect the outcome but not the exposure), `"null_b"` (the exposure
#' affects the mediators but they do not affect the outcome) and `"null_ab"`
#' (neither). The true NIE is 0 in all three. In `"null_a"` the b-paths are
#' calibrated as in the corresponding spiked scenario, so type I error and
#' power runs are matched on effect sizes.
#'
#' @param n Number of samples per replicate.
#' @param template A fitted `mimesys_template` or input accepted by
#'   [fit_template()].
#' @param p Feature counts per view; see [sim_from_template()].
#' @param crosstalk Cross-view dependence; see [sim_from_template()].
#' @param n.mediators Number of mediators per view: one value for every view,
#'   one per view, or a vector named by view (use 0 for a view without
#'   mediators).
#' @param TIE Target total (natural) indirect effect, a non-negative number.
#'   Its sign is the product of the a- and b-path directions.
#' @param DE Direct effect of the exposure on the outcome.
#' @param dir.a,dir.b Directions of the a- and b-paths: `1`, `-1`, or
#'   `"random"` (an independent random sign for each mediator).
#' @param effect.type `"spiked"`, `"null_a"`, `"null_b"` or `"null_ab"`.
#' @param exposure.type `"binary"` or `"continuous"`.
#' @param exposure.prevalence `P(T = 1)` for a binary exposure.
#' @param n.interactions Number of pairwise interactions between mediators in
#'   the outcome model (default 0). Each coefficient has the mean absolute
#'   b-path size, signed by `dir.b` (a random sign per interaction when
#'   `dir.b = "random"`), before calibration. The sign of an interaction's
#'   contribution to the indirect effect also depends on how the two
#'   mediators respond to the exposure, so it need not match the direction of
#'   the linear pathways; the net NIE is calibrated either way and the
#'   interaction part is reported as `NIE_interaction`.
#' @param stages `NULL` (default: every mediator view in one stage, i.e. parallel
#'   mediators) or a list of character vectors of view names in causal order,
#'   for sequential mediation; names of the list, if given, label the stages.
#'   Every view with mediators must be in a stage, and each stage needs at
#'   least one mediator. With a single stage the mediators are parallel and
#'   the stage-specific outputs are omitted.
#' @param stage.links Named list for the links between consecutive stages:
#'   `type` (`"random"`: each parent-child pair linked with probability
#'   `density`, with every child having at least one parent and every parent
#'   at least one child; `"hub"`: a `density` fraction of the parents, at
#'   least one, each feed many children, every child having one hub parent;
#'   `"convergent"`: every child has `max(2, density * n_parents)` parents;
#'   `"explicit"`: give the links yourself in `matrix`), `density` (default
#'   0.3), `strength` (share of each child's latent variance explained by its
#'   parents, in \[0, 1); one value or one per link between stages; default
#'   0.3), `direction` (sign of the links: `1`, `-1` or `"random"`) and
#'   `matrix` (for `"explicit"`: a list with one matrix per link between
#'   stages, rows the mediators of the later stage and columns those of the
#'   earlier stage, in feature order; non-zero entries are links, their
#'   values relative weights and signs; `direction` is not used).
#' @param stage.direct Multiplier on the direct exposure paths of each stage,
#'   one value or one per stage (default 1). Use 0 for a stage the exposure
#'   reaches only through the previous stage.
#' @param stage.share Optional share of the indirect effect carried by the
#'   mediators of each stage (one non-negative value per stage, normalized to
#'   sum to 1), counting every path that ends at that stage. `NULL` (default)
#'   uses one calibration scale for all b-paths.
#' @param sigma Standard deviation of the outcome noise (continuous outcome).
#' @param outcome.type `"continuous"` or `"binary"`. For a binary outcome,
#'   `Y ~ Bernoulli(logit^-1(b0 + DE * T + sum_j b_j M~_j + ...))`: `TIE` and
#'   `DE` are on the log-odds scale (`NIE`, `NDE`, `TE` in the truth), the
#'   intercept `b0` is solved so that the population prevalence equals
#'   `outcome.prevalence`, and counterfactual risk-difference and odds-ratio
#'   effects are reported in `truth$mediation$binary`.
#' @param outcome.prevalence Target `P(Y = 1)` for a binary outcome.
#' @param min.variability Mediators are drawn from features whose most common
#'   value covers at most `1 - min.variability` of the template samples, so
#'   near-constant features (for example, almost always zero) are not used.
#' @param n.mc Monte Carlo sample size used to compute the true effects.
#' @param feature.select,feature.jitter,view.models,back.transform Passed to
#'   [sim_from_template()]. The outcome is always computed from the fitted
#'   (transformed) scale.
#' @param nrep Number of replicates.
#' @param seed Random seed (the caller's random state is restored), or `NULL`.
#' @param return.template Include the fitted template in the result.
#' @param ... Passed to [fit_template()] when `template` is not fitted.
#'
#' @return A list with `replicates` (one pcl-style list per replicate),
#'   `settings`, `template_provenance` and `template` (when requested). Each
#'   replicate has:
#'   \describe{
#'     \item{feature_table, sample_metadata}{Features by samples; sample
#'       metadata with `exposure` and `Y`.}
#'     \item{feature_metadata}{Per feature: `mediator`, `a_path` (latent
#'       correlation with the exposure), `a_effect` (true exposure effect on
#'       the standardized mediator for the contrast), `a_slope` (population
#'       linear slope on the exposure), `b_path`, `indirect`
#'       (`b_path * a_effect`), `exposure_effect` (true exposure effect on
#'       every feature) and `interaction_term`; with `stages`, also
#'       `stage`, `n_parents` and `n_children`. With stages, `a_path` is the
#'       direct exposure path and `a_effect` the total exposure effect
#'       (direct and through earlier stages).}
#'     \item{truth$mediation}{`NIE`, `NDE`, `TE`; `TIE_target` (`NA` with negative links or a
#'       random direction); `pathway_effect_total`; `NIE_linear`;
#'       `NIE_by_view`; `NIE_interaction`; with `stages`, `path_effects`
#'       (natural path-specific effects, which sum to `NIE`), `NIE_by_stage`,
#'       `stages` and `links` (one row per parent -> child link with its
#'       latent coefficient and `approx_chain_effect`, a latent-scale linear
#'       approximation of the effect carried along that link; the exact
#'       totals are in `path_effects`); `interaction_pairs`; `contrast`;
#'       `effect.type`, `exposure.type`, `outcome.type`; `scale` and
#'       `feature_scale`; and, for a binary outcome, `binary` with the
#'       intercept, counterfactual risks, natural effects as risk
#'       differences, odds ratios and log odds ratios, and the path-specific
#'       risk differences (`path_risk_difference`).}
#'   }
#' @examples
#' data(stelzer_mae)
#' tmpl <- fit_template(stelzer_mae, shared.rank = 2, seed = 1)
#' med <- gen_mediation(
#'     n = 200, template = tmpl, p = 20,
#'     n.mediators = c(CYTOF = 2, Metabolomics = 2, Proteomics = 0),
#'     TIE = 0.5, DE = 0.2, n.mc = 2000, seed = 1
#' )
#' rep1 <- med$replicates[[1]]
#' rep1$truth$mediation[c("NIE", "NDE", "TE")]
#' rep1$truth$mediation$NIE_by_view
#' fm <- rep1$feature_metadata
#' fm[fm$mediator, c("featureType", "a_path", "a_effect", "b_path", "indirect")]
#'
#' ## Sequential mediation: exposure -> CyTOF -> metabolome -> outcome
#' seq_med <- gen_mediation(
#'     n = 200, template = tmpl, p = 20,
#'     n.mediators = c(CYTOF = 3, Metabolomics = 3, Proteomics = 0),
#'     stages = list(Immune = "CYTOF", Metabolic = "Metabolomics"),
#'     stage.links = list(type = "random", density = 0.5, strength = 0.3),
#'     TIE = 0.5, n.mc = 2000, seed = 1
#' )
#' seq_med$replicates[[1]]$truth$mediation$path_effects
#' seq_med$replicates[[1]]$truth$mediation$links
#' @export
gen_mediation <- function(n,
    template,
    p = NULL,
    crosstalk = 1,
    n.mediators = 5,
    TIE = 1,
    DE = 0,
    dir.a = 1,
    dir.b = 1,
    effect.type = c("spiked", "null_a", "null_b", "null_ab"),
    exposure.type = c("binary", "continuous"),
    exposure.prevalence = 0.5,
    n.interactions = 0,
    stages = NULL,
    stage.links = list(),
    stage.direct = 1,
    stage.share = NULL,
    sigma = 1,
    outcome.type = c("continuous", "binary"),
    outcome.prevalence = 0.5,
    min.variability = 0.2,
    n.mc = 5000,
    feature.select = c("random", "variance"),
    feature.jitter = 0.2,
    view.models = NULL,
    back.transform = FALSE,
    nrep = 1,
    seed = 1234,
    return.template = FALSE,
    ...) {
    if (!is.null(seed)) withr::local_seed(seed)
    effect.type <- match.arg(effect.type)
    exposure.type <- match.arg(exposure.type)
    outcome.type <- match.arg(outcome.type)
    if (identical(outcome.type, "binary")) {
        .check_prevalence(outcome.prevalence, "outcome.prevalence")
    }
    feature.select <- match.arg(feature.select)
    .check_count(n, "n", min = 4L)
    .check_count(nrep, "nrep", min = 1L)
    .check_count(n.mc, "n.mc", min = 100L)
    if (!is.numeric(TIE) || length(TIE) != 1L || !is.finite(TIE) || TIE < 0) {
        stop("TIE must be a single non-negative number; its sign comes from dir.a and dir.b.",
            call. = FALSE)
    }
    if (!is.numeric(DE) || length(DE) != 1L || !is.finite(DE)) {
        stop("DE must be a single finite number.", call. = FALSE)
    }
    if (!is.numeric(sigma) || length(sigma) != 1L || !is.finite(sigma) || sigma < 0) {
        stop("sigma must be a single non-negative number.", call. = FALSE)
    }
    .check_unit_interval(min.variability, "min.variability", single = TRUE)
    if (identical(exposure.type, "binary")) {
        .check_prevalence(exposure.prevalence, "exposure.prevalence")
    }
    .check_count(n.interactions, "n.interactions", min = 0L)
    .check_direction(dir.a, "dir.a")
    .check_direction(dir.b, "dir.b")

    tmpl <- .resolve_template(template, list(...))
    views <- tmpl$views
    n_med <- .per_view_values(n.mediators, views, "n.mediators")
    if (!is.numeric(n_med) || any(!is.finite(n_med)) || any(n_med < 0) ||
        any(abs(n_med - round(n_med)) > 1e-8)) {
        stop("n.mediators must be non-negative integers.", call. = FALSE)
    }
    n_med <- stats::setNames(as.integer(round(n_med)), views)
    K <- sum(n_med)
    if (K < 1L) stop("At least one mediator is needed (n.mediators).", call. = FALSE)
    if (n.interactions > choose(K, 2)) {
        stop("n.interactions = ", n.interactions, " exceeds the ", choose(K, 2),
            " possible pairs of ", K, " mediators.", call. = FALSE)
    }

    ## Exposure scale and contrast
    if (identical(exposure.type, "binary")) {
        mu_t <- exposure.prevalence
        sd_t <- sqrt(exposure.prevalence * (1 - exposure.prevalence))
    } else {
        mu_t <- 0
        sd_t <- 1
    }
    contrast <- c(t0 = 0, t1 = 1)
    ## Support of the standardized exposure, used to keep the mediators'
    ## latent scores exactly standard normal for a binary exposure
    exposure_support <- if (identical(exposure.type, "binary")) {
        list(values = (c(0, 1) - mu_t) / sd_t,
            probs = c(1 - exposure.prevalence, exposure.prevalence))
    } else {
        NULL
    }

    ## Stages (sequential mediation) and their links
    stage_info <- .mediation_stage_of_view(stages, views, n_med)
    n_stages <- stage_info$n_stages
    stage_names <- stage_info$names
    link_opt <- .mediation_link_options(stage.links, n_stages)
    sd_mult <- stage.direct
    if (!is.numeric(sd_mult) || !(length(sd_mult) %in% c(1L, n_stages)) ||
        any(!is.finite(sd_mult)) || any(sd_mult < 0)) {
        stop("stage.direct must be non-negative: one value, or one per stage.", call. = FALSE)
    }
    sd_mult <- rep_len(sd_mult, n_stages)
    if (!is.null(stage.share)) {
        if (!is.numeric(stage.share) || length(stage.share) != n_stages ||
            any(!is.finite(stage.share)) || any(stage.share < 0) || sum(stage.share) <= 0) {
            stop("stage.share must give one non-negative value per stage, not all zero.",
                call. = FALSE)
        }
        stage.share <- stage.share / sum(stage.share)
    }
    med_view <- rep(views, n_med)
    med_view_idx <- match(med_view, views)
    med_stage <- unname(stage_info$stage_of_view[med_view])
    paths <- .mediation_paths(n_stages)
    n_paths <- length(paths)
    path_s <- vapply(paths, `[[`, numeric(1), "s")
    path_r <- vapply(paths, `[[`, numeric(1), "r")
    ## Counterfactual configurations: configuration i has paths 1..i switched
    ## to the exposed level; 0 is all unexposed, n_paths all exposed
    cfgs <- lapply(0:n_paths, function(i) {
        m <- matrix(0, n_stages, n_stages)
        for (h in seq_len(i)) m[paths[[h]][["r"]], paths[[h]][["s"]]] <- 1
        m
    })

    ## Path size: with a = b = s, NIE is about K s^2 / sd_t before calibration
    s <- sqrt(TIE * sd_t / K)
    if (s >= 0.95) {
        stop("TIE = ", TIE, " is too large for ", K, " mediators: the a-paths would ",
            "need a mixing weight of ", sprintf("%.2f", s), " with the exposure. ",
            "Use more mediators or a smaller TIE.", call. = FALSE)
    }
    signed_tie <- 0
    ## Net calibration needs every path to have the sign dir.a * dir.b, i.e.
    ## positive links (explicit links: all weights positive)
    links_positive <- n_stages == 1L || if (identical(link_opt$type, "explicit")) {
        all(unlist(link_opt$matrix) >= 0)
    } else {
        is.numeric(link_opt$direction) && link_opt$direction == 1
    }
    fixed_dirs <- !identical(dir.a, "random") && !identical(dir.b, "random") &&
        links_positive

    replicates <- vector("list", nrep)
    names(replicates) <- paste("Rep", seq_len(nrep), sep = "_")
    for (k in seq_len(nrep)) {
        ## 1. Features, mediators and paths
        structure <- attr(.sim_from_template_impl(tmpl, n = 2L, p = p, crosstalk = crosstalk,
            feature.select = feature.select, feature.jitter = feature.jitter,
            view.models = view.models), "structure")
        med_pos <- .select_mediators(tmpl, structure$parent, n_med, min.variability)
        med_col <- unlist(med_pos[n_med > 0L], use.names = FALSE)
        p_counts <- lengths(structure$parent)
        offsets <- c(0L, cumsum(p_counts))[seq_along(views)]
        rows <- unlist(lapply(seq_along(views), function(m) offsets[m] + med_pos[[m]]),
            use.names = FALSE)
        alpha_sp <- .path_signs(dir.a, K) * s * sd_mult[med_stage]
        beta_sp <- .path_signs(dir.b, K) * s
        pairs <- if (n.interactions > 0) {
            cmb <- utils::combn(K, 2)
            cmb[, sample.int(ncol(cmb), n.interactions), drop = FALSE]
        } else {
            matrix(integer(0), 2, 0)
        }
        ## Interaction coefficients: the mean b-path size, signed by dir.b
        ## (random sign per interaction when dir.b is "random")
        gamma_sp <- mean(abs(beta_sp)) * .path_signs(dir.b, ncol(pairs))
        ## Links between consecutive stages
        patterns <- if (n_stages > 1L) {
            c(list(NULL), lapply(2:n_stages, function(st) {
                .mediation_link_pattern(sum(med_stage == st), sum(med_stage == st - 1L),
                    link_opt, st - 1L)
            }))[-1L]
        } else {
            list()
        }
        gam_view <- .per_view_values(crosstalk, views, "crosstalk")
        Sigma <- if (n_stages > 1L) {
            .mediation_latent_cov(tmpl, structure, gam_view, med_view_idx, med_col)
        } else {
            diag(K)
        }
        design_sp <- .mediation_design(alpha_sp, med_stage, n_stages, Sigma, patterns,
            link_opt$strength)

        ## Standardization constants (template mean and SD) for every feature
        std <- .feature_standardization(tmpl, structure$parent)

        ## 2. Paired counterfactual Monte Carlo (same random stream for every
        ## exposure value, so differences isolate the exposure effect)
        ## The Monte Carlo runs in chunks (each with its own seed, shared across
        ## exposure values) so that only per-feature means and the mediator
        ## columns are kept, not the full n.mc x features matrix.
        chunk_size <- 1000L
        chunks <- split(seq_len(n.mc), ceiling(seq_len(n.mc) / chunk_size))
        mc_seeds <- sample.int(.Machine$integer.max, length(chunks))
        ## `cfg`: a stage x stage matrix of 0/1 contrast levels (see cfgs), or
        ## a vector of exposure values used for every path
        draw_all <- function(cfg, design) {
            t_std <- function(t_value) (t_value - mu_t) / sd_t
            sums <- 0
            med <- vector("list", length(chunks))
            for (ci in seq_along(chunks)) {
                ii <- chunks[[ci]]
                t_fun <- if (is.matrix(cfg)) {
                    function(r, s) t_std(contrast[[cfg[r, s] + 1L]])
                } else {
                    tv <- t_std(rep_len(cfg, n.mc)[ii])
                    function(r, s) tv
                }
                hook <- .mediation_stage_hook(design, med_view_idx, med_col, t_fun,
                    exposure_support)
                sim <- withr::with_seed(mc_seeds[ci], .sim_from_template_impl(tmpl,
                    n = length(ii), p = p, crosstalk = crosstalk,
                    view.models = view.models, .structure = structure,
                    .latent_hook = hook))
                x <- .standardize_features(sim, std)
                sums <- sums + colSums(x)
                med[[ci]] <- x[, rows, drop = FALSE]
            }
            list(mean = sums / n.mc, med = do.call(rbind, med))
        }
        draw_cfgs <- function(design) {
            ## The all-exposed configuration first, as in earlier versions
            out <- vector("list", n_paths + 1L)
            out[[n_paths + 1L]] <- draw_all(cfgs[[n_paths + 1L]], design)
            for (i in rev(seq_len(n_paths))) out[[i]] <- draw_all(cfgs[[i]], design)
            out
        }
        int_mean <- function(x) {
            if (ncol(pairs) == 0) return(numeric(0))
            colMeans(x$med[, pairs[1, ], drop = FALSE] * x$med[, pairs[2, ], drop = FALSE])
        }
        ## Per-path change of each mediator mean and of each interaction mean
        path_deltas <- function(draws) {
            lapply(seq_len(n_paths), function(i) {
                list(med = colMeans(draws[[i + 1L]]$med) - colMeans(draws[[i]]$med),
                    int = int_mean(draws[[i + 1L]]) - int_mean(draws[[i]]))
            })
        }
        ## Contribution of each path (summed over mediators, or summed
        ## absolute contributions when `absolute`)
        path_effects <- function(deltas, beta, gamma, absolute = FALSE) {
            f <- if (absolute) abs else identity
            vapply(deltas, function(d) sum(f(beta * d$med)) + sum(f(gamma * d$int)),
                numeric(1))
        }

        ## Calibrate the b-paths as in the spiked scenario (also for null_a, so
        ## null and spiked runs have matched b-path sizes)
        beta <- beta_sp
        gamma <- gamma_sp
        if (effect.type %in% c("spiked", "null_a") && TIE > 0) {
            cal_draws <- draw_cfgs(design_sp)
            dl <- path_deltas(cal_draws)
            ## Fixed directions: calibrate the net indirect effect to +/- TIE.
            ## Random directions: pathway effects can cancel, so calibrate the
            ## total magnitude of the pathway-specific effects instead
            target <- if (fixed_dirs) TIE * dir.a * dir.b else TIE
            absolute <- !fixed_dirs
            if (is.null(stage.share)) {
                nie_raw <- sum(path_effects(dl, beta_sp, gamma_sp, absolute))
                if (!is.finite(nie_raw) || abs(nie_raw) < 1e-8) {
                    stop("The simulated indirect effect is essentially zero; check the mediators.",
                        call. = FALSE)
                }
                beta <- beta_sp * target / nie_raw
                gamma <- gamma_sp * target / nie_raw
            } else {
                ## Per-stage scale: the stage-s b-paths carry stage.share[s]
                ## of the indirect effect (exact without interactions)
                lambda <- numeric(n_stages)
                for (st in seq_len(n_stages)) {
                    if (stage.share[st] == 0) next
                    raw_s <- sum(path_effects(dl[path_s == st], beta_sp * (med_stage == st),
                        rep(0, length(gamma_sp)), absolute))
                    if (!is.finite(raw_s) || abs(raw_s) < 1e-8) {
                        stop("The exposure does not reach the stage-", st, " mediators ",
                            "(check stage.direct and stage.links), so they cannot carry ",
                            "a share of the indirect effect.", call. = FALSE)
                    }
                    lambda[st] <- stage.share[st] * target / raw_s
                }
                beta <- beta_sp * lambda[med_stage]
                gamma <- if (ncol(pairs) > 0) {
                    gamma_sp * (lambda[med_stage[pairs[1, ]]] + lambda[med_stage[pairs[2, ]]]) / 2
                } else {
                    gamma_sp
                }
                if (ncol(pairs) > 0) {
                    total <- sum(path_effects(dl, beta, gamma, absolute))
                    if (abs(total) > 1e-8) {
                        beta <- beta * target / total
                        gamma <- gamma * target / total
                    }
                }
            }
            signed_tie <- if (fixed_dirs) target else NA_real_
        }
        alpha <- if (effect.type %in% c("spiked", "null_b")) alpha_sp else rep(0, K)
        if (effect.type %in% c("null_b", "null_ab")) {
            beta <- rep(0, K)
            gamma <- rep(0, ncol(pairs))
        }
        ## Same links in every scenario; only the exposure paths change
        design <- if (identical(alpha, alpha_sp)) design_sp else {
            .mediation_design(alpha, med_stage, n_stages, Sigma, patterns,
                link_opt$strength, D = design_sp$D)
        }

        ## True effects for the scenario as simulated (the spiked scenario
        ## reuses the calibration draws, which used the same a-paths)
        draws <- if (identical(effect.type, "spiked") && TIE > 0) cal_draws else draw_cfgs(design)
        x1 <- draws[[n_paths + 1L]]
        x0 <- draws[[1L]]
        dl <- path_deltas(draws)
        exposure_effect <- x1$mean - x0$mean
        a_effect <- exposure_effect[rows]
        int_effect <- int_mean(x1) - int_mean(x0)
        nie_view <- vapply(views, function(v) {
            sum((beta * a_effect)[med_view == v])
            }, numeric(1))
        nie_int <- sum(gamma * int_effect)
        path_eff <- path_effects(dl, beta, gamma)
        nie <- sum(nie_view) + nie_int
        de_effect <- DE * (contrast[["t1"]] - contrast[["t0"]])
        lin <- function(m, t_value) {
            out <- DE * t_value + as.numeric(m %*% beta)
            if (ncol(pairs) > 0) {
                out <- out + as.numeric((m[, pairs[1, ], drop = FALSE] *
                    m[, pairs[2, ], drop = FALSE]) %*% gamma)
            }
            out
        }
        ## Continuous exposure: population linear projection of each mediator on
        ## the exposure (what product-of-coefficients methods estimate); the
        ## 0 -> 1 contrast above can differ because mediators respond
        ## non-linearly to the exposure
        a_slope <- a_effect
        t_pop <- m_pop <- NULL
        if (identical(exposure.type, "continuous")) {
            t_pop <- stats::rnorm(n.mc)
            m_pop <- draw_all(t_pop, design)$med
            a_slope <- as.numeric(stats::cov(t_pop, m_pop)) / stats::var(t_pop)
        }
        nie_linear <- if (ncol(pairs) == 0) sum(beta * a_slope) else NA_real_
        binary_truth <- NULL
        intercept <- 0
        if (identical(outcome.type, "binary")) {
            m1 <- x1$med
            m0 <- x0$med
            ## Intercept so that the population outcome prevalence is as
            ## requested (exposure drawn from its distribution)
            intercept <- if (identical(exposure.type, "binary")) {
                eta1 <- lin(m1, contrast[["t1"]])
                eta0 <- lin(m0, contrast[["t0"]])
                .logistic_intercept(c(eta1, eta0), outcome.prevalence,
                    weights = c(rep(exposure.prevalence / length(eta1), length(eta1)),
                        rep((1 - exposure.prevalence) / length(eta0), length(eta0))))
            } else {
                .logistic_intercept(lin(m_pop, t_pop), outcome.prevalence)
            }
            ## Counterfactual (natural) effects: p(t, M(t'))
            p_cfg <- vapply(draws, function(d) {
                mean(stats::plogis(intercept + lin(d$med, contrast[["t1"]])))
            }, numeric(1))
            p11 <- p_cfg[[n_paths + 1L]]
            p10 <- p_cfg[[1L]]
            p00 <- mean(stats::plogis(intercept + lin(m0, contrast[["t0"]])))
            odds <- function(q) q / (1 - q)
            or <- c(NIE = odds(p11) / odds(p10), NDE = odds(p10) / odds(p00),
                TE = odds(p11) / odds(p00))
            binary_truth <- list(
                intercept = intercept,
                risk = c(p11 = p11, p10 = p10, p00 = p00),
                risk_difference = c(NIE = p11 - p10, NDE = p10 - p00, TE = p11 - p00),
                odds_ratio = or,
                log_odds_ratio = log(or),
                path_risk_difference = stats::setNames(diff(p_cfg),
                    vapply(paths, function(h) .mediation_path_label(h[["r"]], h[["s"]],
                        stage_names), character(1)))
            )
        }

        ## 3. Observed data
        exposure <- if (identical(exposure.type, "binary")) {
            stats::rbinom(n, 1L, exposure.prevalence)
        } else {
            stats::rnorm(n)
        }
        t_star <- (exposure - mu_t) / sd_t
        pcl <- .sim_from_template_impl(tmpl, n = n, p = p, crosstalk = crosstalk,
            view.models = view.models,
            back.transform = if (isTRUE(back.transform)) "both" else FALSE,
            .structure = structure,
            .latent_hook = .mediation_stage_hook(design, med_view_idx, med_col,
                function(r, s) t_star, exposure_support))
        mt <- .standardize_features(pcl, std)[, rows, drop = FALSE]
        eta <- lin(mt, exposure)
        y <- if (identical(outcome.type, "binary")) {
            stats::rbinom(n, 1L, stats::plogis(intercept + eta))
        } else {
            eta + sigma * stats::rnorm(n)
        }
        if (!is.null(pcl$original_table)) {
            pcl$feature_table <- pcl$original_table
            pcl$original_table <- NULL
            pcl$feature_metadata$scale <- "original"
        }
        pcl$sample_metadata$exposure <- exposure
        pcl$sample_metadata$Y <- y

        ## Feature-level truth
        fm <- pcl$feature_metadata
        active_pairs <- if (ncol(pairs) > 0) gamma != 0 else logical(0)
        fm$mediator <- FALSE
        fm$a_path <- fm$a_effect <- fm$a_slope <- fm$b_path <- fm$indirect <- 0
        fm$interaction_term <- FALSE
        fm$exposure_effect <- exposure_effect
        fm$mediator[rows] <- TRUE
        fm$a_path[rows] <- alpha
        fm$a_effect[rows] <- a_effect
        fm$a_slope[rows] <- a_slope
        fm$b_path[rows] <- beta
        fm$indirect[rows] <- beta * a_effect
        if (any(active_pairs)) {
            fm$interaction_term[rows[unique(as.vector(pairs[, active_pairs, drop = FALSE]))]] <- TRUE
        }
        links_df <- NULL
        if (n_stages > 1L) {
            fm$stage <- NA_integer_
            fm$stage[rows] <- med_stage
            fm$n_parents <- fm$n_children <- 0L
            links_df <- do.call(rbind, lapply(2:n_stages, function(st) {
                Ds <- design$D[[st]]
                ch <- rows[med_stage == st]
                pa <- rows[med_stage == st - 1L]
                nz <- which(Ds != 0, arr.ind = TRUE)
                if (nrow(nz) == 0L) return(NULL)
                data.frame(
                    parent = fm$featureID[pa[nz[, 2]]],
                    child = fm$featureID[ch[nz[, 1]]],
                    parent_view = fm$featureType[pa[nz[, 2]]],
                    child_view = fm$featureType[ch[nz[, 1]]],
                    parent_stage = st - 1L, child_stage = st,
                    coefficient = Ds[nz],
                    ## Latent-scale linear approximation: the exposure moves the
                    ## parent's latent score by mu / sd_t for the contrast,
                    ## the link passes d of it to the child, whose b-path
                    ## carries it to the outcome
                    approx_chain_effect = design$mu[match(pa[nz[, 2]], rows)] /
                        sd_t * (contrast[["t1"]] - contrast[["t0"]]) * Ds[nz] *
                        fm$b_path[ch[nz[, 1]]],
                    stringsAsFactors = FALSE
                )
            }))
            if (!is.null(links_df)) {
                rownames(links_df) <- NULL
                np <- table(factor(links_df$child, levels = fm$featureID))
                nc <- table(factor(links_df$parent, levels = fm$featureID))
                fm$n_parents <- as.integer(np[fm$featureID])
                fm$n_children <- as.integer(nc[fm$featureID])
            }
        }
        pcl$feature_metadata <- fm

        path_table <- data.frame(
            path = vapply(paths, function(h) .mediation_path_label(h[["r"]], h[["s"]],
                stage_names), character(1)),
            from_stage = as.integer(path_r), to_stage = as.integer(path_s),
            effect = path_eff,
            stringsAsFactors = FALSE
        )
        pcl$truth$mediation <- list(
            NIE = nie,
            NDE = de_effect,
            TE = nie + de_effect,
            TIE_target = if (identical(effect.type, "spiked")) signed_tie else 0,
            pathway_effect_total = sum(path_effects(dl, beta, gamma, absolute = TRUE)),
            NIE_linear = nie_linear,
            NIE_by_view = nie_view,
            NIE_interaction = nie_int,
            path_effects = if (n_stages > 1L) path_table else NULL,
            NIE_by_stage = if (n_stages > 1L) {
                stats::setNames(vapply(seq_len(n_stages), function(st) {
                    sum(path_eff[path_s == st])
                }, numeric(1)), stage_names)
            } else {
                NULL
            },
            links = links_df,
            stages = if (n_stages > 1L) stats::setNames(lapply(seq_len(n_stages), function(st) {
                views[stage_info$stage_of_view == st]
            }), stage_names) else NULL,
            interaction_pairs = if (ncol(pairs) > 0) data.frame(
                feature_x = fm$featureID[rows[pairs[1, ]]],
                feature_y = fm$featureID[rows[pairs[2, ]]],
                coefficient = gamma,
                stringsAsFactors = FALSE
            ) else NULL,
            contrast = contrast,
            effect.type = effect.type,
            exposure.type = exposure.type,
            outcome.type = outcome.type,
            scale = if (identical(outcome.type, "binary")) {
                "conditional log-odds (linear predictor)"
            } else {
                "outcome"
            },
            feature_scale = "fitted",
            binary = binary_truth
        )
        replicates[[k]] <- pcl
    }

    list(
        replicates = replicates,
        template = if (isTRUE(return.template)) tmpl else NULL,
        template_provenance = tmpl$provenance,
        settings = list(
            n = n, p = p, crosstalk = crosstalk,
            n.mediators = n_med, TIE = TIE, DE = DE, dir.a = dir.a, dir.b = dir.b,
            effect.type = effect.type, exposure.type = exposure.type,
            exposure.prevalence = if (identical(exposure.type, "binary")) exposure.prevalence else NULL,
            n.interactions = n.interactions, stages = stages,
            stage.links = if (n_stages > 1L) link_opt else NULL,
            stage.direct = stage.direct, stage.share = stage.share, sigma = sigma,
            outcome.type = outcome.type,
            outcome.prevalence = if (identical(outcome.type, "binary")) outcome.prevalence else NULL,
            min.variability = min.variability, n.mc = n.mc,
            view.models = view.models, back.transform = back.transform,
            nrep = nrep, seed = seed
        )
    )
}

## ---- mediation helpers --------------------------------------------------

.check_direction <- function(x, name) {
    ok <- length(x) == 1L && (identical(x, "random") ||
        (is.numeric(x) && x %in% c(-1, 1)))
    if (!ok) stop(name, " must be 1, -1 or \"random\".", call. = FALSE)
    invisible(x)
}

.path_signs <- function(direction, K) {
    if (identical(direction, "random")) {
        sample(c(-1, 1), K, replace = TRUE)
    } else {
        rep(as.numeric(direction), K)
    }
}

## Share of template samples not at a feature's most common value
.feature_variability <- function(sorted) {
    r <- rle(sorted)
    1 - max(r$lengths) / length(sorted)
}

.select_mediators <- function(tmpl, parent_list, n_med, min.variability) {
    views <- tmpl$views
    out <- lapply(seq_along(views), function(m) {
        k <- n_med[[m]]
        if (k == 0L) return(integer(0))
        qmat <- tmpl$fits[[m]]$qmat
        parent <- parent_list[[m]]
        var_ok <- vapply(parent, function(j) .feature_variability(qmat[, j]),
            numeric(1)) >= min.variability
        eligible <- which(var_ok)
        if (length(eligible) < k) {
            stop("View ", views[m], " has only ", length(eligible), " eligible feature(s) for ",
                k, " mediators (features must vary in at least ",
                round(100 * min.variability), "% of template samples). ",
                "Lower n.mediators or min.variability, or increase p.", call. = FALSE)
        }
        sort(.sample_values(eligible, k))
    })
    stats::setNames(out, views)
}

## Template mean and SD of every simulated feature (fitted scale), in
## feature order
.feature_standardization <- function(tmpl, parent_list) {
    mu <- sds <- list()
    for (m in seq_along(tmpl$views)) {
        q <- tmpl$fits[[m]]$qmat[, parent_list[[m]], drop = FALSE]
        mu[[m]] <- colMeans(q)
        sdm <- apply(q, 2, stats::sd)
        sdm[!is.finite(sdm) | sdm == 0] <- 1
        sds[[m]] <- sdm
    }
    list(mean = unlist(mu, use.names = FALSE), sd = unlist(sds, use.names = FALSE))
}

## Samples x features matrix of standardized features (fitted scale)
.standardize_features <- function(sim, std) {
    x <- t(as.matrix(sim$feature_table))
    sweep(sweep(x, 2, std$mean, "-"), 2, std$sd, "/")
}

## design = "mediation" of simulate_mimesys(): MAE interface to gen_mediation()
.simulate_mediation_mae <- function(template,
    n,
    p = NULL,
    crosstalk = "calibrate",
    mediation = list(),
    null_model = list(),
    template_args = list(),
    calibration = list(),
    nrep = 1,
    seed = 1234,
    return_type = c("mae", "simulation")) {
    return_type <- match.arg(return_type)
    if (!is.null(seed)) withr::local_seed(seed)
    mediation_defaults <- list(
        n_mediators = 5, tie = 1, de = 0, dir_a = 1, dir_b = 1,
        effect_type = "spiked", exposure_type = "binary",
        exposure_prevalence = 0.5, n_interactions = 0, sigma = 1,
        outcome_type = "continuous", outcome_prevalence = 0.5,
        min_variability = 0.2, n_mc = 5000,
        stages = NULL, stage_links = list(), stage_direct = 1, stage_share = NULL
    )
    mediation <- .mimesys_merge_options(mediation, mediation_defaults, "mediation")
    setup <- .mimesys_wrapper_setup(template, n = n, p = p, nrep = nrep,
        crosstalk = crosstalk, null_model = null_model,
        template_args = template_args, calibration = calibration, seed = seed,
        min_n = 4L)
    tmpl <- setup$template
    null_model <- setup$null_model
    resolved <- setup$crosstalk
    sim <- gen_mediation(
        n = n, template = tmpl, p = p, crosstalk = resolved$value,
        n.mediators = mediation$n_mediators, TIE = mediation$tie, DE = mediation$de,
        dir.a = mediation$dir_a, dir.b = mediation$dir_b,
        effect.type = mediation$effect_type, exposure.type = mediation$exposure_type,
        exposure.prevalence = mediation$exposure_prevalence,
        n.interactions = mediation$n_interactions,
        stages = mediation$stages, stage.links = mediation$stage_links,
        stage.direct = mediation$stage_direct, stage.share = mediation$stage_share,
        sigma = mediation$sigma,
        outcome.type = mediation$outcome_type,
        outcome.prevalence = mediation$outcome_prevalence,
        min.variability = mediation$min_variability, n.mc = mediation$n_mc,
        feature.select = null_model$feature_select,
        feature.jitter = null_model$feature_jitter,
        view.models = null_model$view_models,
        back.transform = isTRUE(null_model$back_transform),
        nrep = nrep, seed = NULL
    )
    sim$crosstalk_calibration <- resolved$calibration
    if (identical(return_type, "simulation")) return(sim)
    settings <- sim$settings
    settings$crosstalk <- resolved$value
    settings$purpose <- "mediation"
    settings$seed <- seed
    settings$crosstalk_calibration <- resolved$calibration
    settings$provenance <- .mimesys_wrapper_provenance(tmpl, resolved)
    maes <- lapply(sim$replicates, function(rep_k) {
        .mimesys_replicate_to_mae(rep_k, simulation = settings, split = "all")
    })
    .mimesys_collect_maes(maes, settings)
}
