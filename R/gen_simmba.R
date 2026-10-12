#' Generate simulated multi-omics data for benchmarking
#'
#' Generates null multiview data from a real template (see \code{fit_template()}
#' and \code{sim_from_template()}), spikes in signal, generates the outcome, and
#' splits the data into training and testing sets.
#'
#' @param nsample Number of samples to simulate in each replicate before the
#'   train/test split.
#' @param template Either a \code{mimesys_template} from \code{fit_template()}
#'   (recommended: fit once, reuse) or a named list of sample x feature matrices,
#'   which is then fitted with \code{fit_template(template, ...)}.
#' @param p Feature counts per view. Use \code{NULL} to keep template sizes, a
#'   single integer to use the same count for every view, an unnamed vector in
#'   template view order, or a named vector with one value per view.
#' @param crosstalk Multiplier on template cross-view dependence. Use a scalar
#'   for the same value in every view or one value per view. A value of 0 removes
#'   shared latent structure, 1 uses the fitted template level, and values above
#'   1 amplify shared structure up to \code{crosstalk_max(template)}.
#' @param feature.select Feature choice rule passed to \code{sim_from_template()}
#'   when \code{p} is below the template feature count. Use \code{"random"} or
#'   \code{"variance"}.
#' @param feature.jitter Offspring-feature perturbation passed to
#'   \code{sim_from_template()} when \code{p} exceeds the template feature count.
#'   Values near 0 make offspring closer to their parent features; larger values
#'   add more loading-level perturbation.
#' @param dependence Dependence simulator passed to \code{sim_from_template()}.
#'   `"factor"` keeps the current low-rank template simulator; `"empirical_copula"`
#'   uses the empirical latent feature correlation more directly; and
#'   `"module_copula"` smooths the empirical latent correlation by feature
#'   modules.
#' @param cor.shrinkage Shrinkage toward independence for
#'   \code{dependence = "empirical_copula"} or `"module_copula"`.
#' @param module.count Number of within-view feature modules for
#'   \code{dependence = "module_copula"}. If `NULL`, a data-adaptive default is
#'   used per view.
#' @param module.between.shrinkage How strongly to smooth correlations between
#'   modules in \code{dependence = "module_copula"}.
#' @param view.models Optional named vector of per-view observation models passed
#'   to \code{sim_from_template()}. Use `"two_part_copula"` for sparse
#'   compositional microbiome views; unspecified views use `"continuous"`.
#' @param snr Signal-to-noise ratio for continuous outcomes. Larger values
#'   create easier prediction tasks.
#' @param p.train Proportion of samples assigned to the training split.
#' @param de.prob Probability that a feature receives a nonzero spike-in
#'   coefficient; scalar or one value per view.
#' @param de.downProb Probability that a nonzero spike-in coefficient is
#'   negative; scalar or one value per view.
#' @param de.facLoc Location parameter for the log-normal effect-size factor;
#'   scalar or one value per view.
#' @param de.facScale Scale parameter for the log-normal effect-size factor;
#'   scalar or one value per view.
#' @param de.align NULL keeps the original spike-in (DE positions and signs drawn
#'   independently per view). A value between 0 and 1 places that fraction of DE
#'   features on features loading on the shared cross-view direction (sampled with
#'   probability proportional to the squared loading) with signs following that
#'   direction, so view-specific signals agree to the extent the views cross-talk.
#' @param signal.share Optional numeric vector controlling the relative marginal
#'   variance contribution of each view's linear predictor. If NULL, spike-in
#'   coefficients are left on their original scale. If supplied, values are
#'   recycled or matched by view name, normalized to sum to one, and each view's
#'   coefficients are rescaled so its view-specific linear predictor has variance
#'   equal to the requested share. For example, \code{c(1, 1, 1)} creates a
#'   balanced complementary-signal setting where no view is intended to dominate.
#' @param friedman.share Optional numeric vector controlling which views supply
#'   the five nonlinear Friedman features when \code{ygen.mode} is
#'   \code{"Friedman"} or \code{"Friedman2"}. Values are recycled or matched by
#'   view name and normalized to sum to one. If NULL, \code{signal.share} is used
#'   when available; otherwise Friedman features are sampled uniformly from all
#'   nonzero spiked features.
#' @param ygen.mode Y generation mode: "LM", "Friedman", "Friedman2",
#'   or "SharedLatent". Friedman modes apply the Friedman formula after
#'   rank-transforming the selected simulated features to the unit interval.
#'   "SharedLatent" uses all simulated shared latent factors, combined as a
#'   standardized equal-weight score.
#' @param outcome.type Outcome type: \code{"continuous"}, \code{"binary"}, or
#'   \code{"survival"}. Continuous and binary outcomes are stored as
#'   \code{Y}; survival outcomes are stored as \code{time} and \code{status}.
#' @param surv.hscale Multiplicative scale on hazard for survival outcomes.
#' @param cens.lower Lower bound for uniform censoring time.
#' @param cens.upper Upper bound for uniform censoring time.
#' @param nrep Number of independent simulation replicates.
#' @param seed Random seed.
#' @param ... Passed to \code{fit_template()} when \code{template} is a raw list
#'
#' @return List with trainDat, testDat, simulation parameters, the fitted
#'   template, realized canonical correlations (\code{cancor}) and per-replicate
#'   empirical agreement between view-specific linear predictors (\code{agreement.emp}).
#' @export
gen_simmba <- function(nsample,
                       template,
                       p = NULL,
                       crosstalk = 1,
                       feature.select = c("random", "variance"),
                       feature.jitter = 0.2,
                       dependence = c("factor", "empirical_copula", "module_copula"),
                       cor.shrinkage = 0.05,
                       module.count = NULL,
                       module.between.shrinkage = 0.75,
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
                       ygen.mode = c("LM", "Friedman", "Friedman2", "SharedLatent"),
                       outcome.type = c("continuous", "binary", "survival"),
                       surv.hscale = 1,
                       cens.lower = 1,
                       cens.upper = 3,
                       nrep = 100,
                       seed = 1234,
                       ...) {

  set.seed(seed)
  ygen.mode      <- match.arg(ygen.mode)
  outcome.type   <- match.arg(outcome.type)
  feature.select <- match.arg(feature.select)
  dependence     <- match.arg(dependence)

  ## 0. Fit the template once (skipped if a fitted template is supplied)
  tmpl  <- if (inherits(template, "mimesys_template")) template else fit_template(template, ...)
  views <- tmpl$views
  M     <- length(views)
  de.prob     <- rep_len(de.prob, M)
  de.downProb <- rep_len(de.downProb, M)
  de.facLoc   <- rep_len(de.facLoc, M)
  de.facScale <- rep_len(de.facScale, M)
  signal.share <- .normalize_signal_share(signal.share, views, "signal.share")
  friedman.share <- .normalize_signal_share(friedman.share, views, "friedman.share")
  if (is.null(friedman.share)) {
    friedman.share <- signal.share
  }

  trainDat <- testDat <- agreement.emp <- vector("list", nrep)
  names(trainDat) <- names(testDat) <- names(agreement.emp) <- paste("Rep", 1:nrep, sep = "_")

  for (k in seq_len(nrep)) {

    ## 1. Null multiview data from the template
    pcl <- sim_from_template(tmpl, n = nsample, p = p, crosstalk = crosstalk,
                             feature.select = feature.select,
                             feature.jitter = feature.jitter,
                             dependence = dependence,
                             cor.shrinkage = cor.shrinkage,
                             module.count = module.count,
                             module.between.shrinkage = module.between.shrinkage,
                             view.models = view.models)
    X  <- t(as.matrix(pcl$feature_table))
    fm <- pcl$feature_metadata

    ## 2. Spike-in coefficients, view by view in template order
    ## (table(featureType) would sort views alphabetically and misassign de.*)
    beta0 <- unlist(lapply(seq_len(M), function(m) {
      sel <- fm$featureType == views[m]
      .draw_beta_view(sum(sel), fm$shared_loading[sel], de.prob[m], de.downProb[m],
                      de.facLoc[m], de.facScale[m], de.align)
    }))
    ## 3. Optional view-level signal balancing, then linear predictor.
    view_pred <- sapply(views, function(v) {
      sel <- fm$featureType == v
      as.numeric(X[, sel, drop = FALSE] %*% beta0[sel])
    })
    if (!is.null(signal.share)) {
      scaled <- .rescale_signal_by_view(beta0, view_pred, fm$featureType, views, signal.share)
      beta0 <- scaled$beta
      view_pred <- scaled$view_pred
    }
    pcl$feature_metadata$beta <- beta0

    eta_lin <- rowSums(view_pred)
    eta_shared <- .shared_latent_predictor(pcl, view_pred)
    A <- suppressWarnings(stats::cor(view_pred))
    A[!is.finite(A)] <- 0
    agreement.emp[[k]] <- A

    ## Optional non-linear Friedman component (beta0 itself is unchanged)
    if (ygen.mode %in% c("Friedman", "Friedman2")) {
      nonzero_index <- which(beta0 != 0)
      candidate_index <- if (length(nonzero_index) >= 5L) nonzero_index else seq_along(beta0)
      friedman_index <- .sample_friedman_features(candidate_index, fm$featureType, views, friedman.share)
      pcl$feature_metadata$friedman_term <- seq_len(nrow(pcl$feature_metadata)) %in% friedman_index
      X.friedman <- .rank_to_unit_interval(X[, friedman_index, drop = FALSE])
      f <- function(x) {
        10 * sin(pi * x[, 1] * x[, 2]) + 20 * (x[, 3] - 0.5)^2 + 10 * x[, 4] + 5 * x[, 5]
      }
      Xbeta.friedman <- f(X.friedman)
    } else {
      pcl$feature_metadata$friedman_term <- FALSE
    }

    pcl$sample_metadata$Xbeta_linear <- eta_lin
    pcl$sample_metadata$Xbeta_shared <- eta_shared

    ## 4. Final linear predictor
    eta <- switch(ygen.mode,
                  LM           = eta_lin,
                  Friedman     = Xbeta.friedman,
                  Friedman2    = eta_lin + Xbeta.friedman,
                  SharedLatent = eta_shared)
    pcl$sample_metadata$Xbeta <- eta

    ## 5. Outcome
    if (outcome.type == "continuous") {
      sigma2 <- as.vector(stats::var(eta) / snr)
      pcl$sample_metadata$Y <- as.vector(eta + stats::rnorm(nsample) * sqrt(sigma2))
    } else if (outcome.type == "binary") {
      pcl$sample_metadata$Y <- stats::rbinom(nsample, size = 1, prob = stats::plogis(eta))
    } else if (outcome.type == "survival") {
      h  <- as.vector(surv.hscale * exp(eta))
      T0 <- stats::rexp(nsample, rate = h)
      C  <- stats::runif(nsample, cens.lower, cens.upper)
      pcl$sample_metadata$time   <- pmin(T0, C)
      pcl$sample_metadata$status <- as.integer(C >= T0)
    }

    ## 6. Train / test split
    train <- test <- pcl
    tr.row <- sample.int(nsample, size = round(nsample * p.train), replace = FALSE)
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

  list(trainDat = trainDat, testDat = testDat,
       template = tmpl, p = p, crosstalk = crosstalk,
       dependence = dependence, cor.shrinkage = cor.shrinkage,
       module.count = module.count, module.between.shrinkage = module.between.shrinkage,
       view.models = view.models,
       cancor = .realized_cancor(tmpl, crosstalk),
       agreement.emp = agreement.emp,
       snr = snr, p.train = p.train,
       de.prob = de.prob, de.downProb = de.downProb, de.facLoc = de.facLoc,
       de.facScale = de.facScale, de.align = de.align,
       signal.share = signal.share, friedman.share = friedman.share,
       nrep = nrep, seed = seed, ygen.mode = ygen.mode,
       outcome.type = outcome.type, surv.hscale = surv.hscale,
       cens.lower = cens.lower, cens.upper = cens.upper)
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
  n_de <- stats::rbinom(1, p_m, de.prob)
  if (n_de == 0) return(beta)
  n_al <- if (is.null(de.align)) 0L else as.integer(round(de.align * n_de))
  pos_al <- if (n_al > 0) sample.int(p_m, n_al, prob = shared_loading^2 + 1e-12) else integer(0)
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

.normalize_signal_share <- function(signal.share, views, arg_name = "signal.share") {
  if (is.null(signal.share)) {
    return(NULL)
  }
  if (!is.numeric(signal.share) || length(signal.share) == 0L) {
    stop(arg_name, " must be a numeric vector.", call. = FALSE)
  }
  if (!is.null(names(signal.share))) {
    missing <- setdiff(views, names(signal.share))
    if (length(missing) > 0L) {
      stop(arg_name, " is missing view(s): ", paste(missing, collapse = ", "), call. = FALSE)
    }
    signal.share <- signal.share[views]
  } else {
    signal.share <- rep_len(signal.share, length(views))
    names(signal.share) <- views
  }
  if (any(!is.finite(signal.share)) || any(signal.share < 0)) {
    stop(arg_name, " values must be finite and non-negative.", call. = FALSE)
  }
  if (sum(signal.share) <= 0) {
    stop("At least one ", arg_name, " value must be positive.", call. = FALSE)
  }
  signal.share / sum(signal.share)
}

.sample_friedman_features <- function(nonzero_index, feature_type, views, friedman.share) {
  if (is.null(friedman.share)) {
    return(sample(nonzero_index, 5))
  }

  view_weight <- friedman.share[as.character(feature_type[nonzero_index])]
  view_weight[!is.finite(view_weight)] <- 0
  positive <- nonzero_index[view_weight > 0]

  if (length(positive) >= 5L) {
    return(sample(positive, 5, prob = view_weight[view_weight > 0]))
  }

  fallback <- setdiff(nonzero_index, positive)
  c(positive, sample(fallback, 5L - length(positive)))
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

.shared_latent_predictor <- function(pcl, view_pred) {
  shared <- pcl$truth$shared
  if (is.matrix(shared) && ncol(shared) > 0L) {
    shared_scaled <- scale(shared)
    shared_scaled[!is.finite(shared_scaled)] <- 0
    return(.standardize_predictor(rowSums(shared_scaled) / sqrt(ncol(shared_scaled))))
  }

  view_pred <- as.matrix(view_pred)
  if (ncol(view_pred) == 0L) {
    stop("Cannot construct shared latent predictor without shared scores or view predictors.", call. = FALSE)
  }
  .standardize_predictor(rowMeans(scale(view_pred), na.rm = TRUE))
}

.rescale_signal_by_view <- function(beta, view_pred, feature_type, views, signal.share) {
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
      warning("Cannot rescale signal for view '", view, "' because its linear predictor variance is zero.", call. = FALSE)
      next
    }
    scale_factor <- sqrt(target / current_var)
    beta[sel] <- beta[sel] * scale_factor
    view_pred[, view] <- view_pred[, view] * scale_factor
  }

  list(beta = beta, view_pred = view_pred)
}
