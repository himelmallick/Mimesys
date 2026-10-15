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
#' \deqn{f_m = \sqrt{\gamma_m} B_m s + (I - \gamma_m B_m B_m^T)^{1/2} \epsilon_m,}
#' \deqn{z_m = \Lambda_m f_m + e_m, \quad x_{mj} = Q_{mj}(\Phi(z_{mj})),}
#' with shared variates \eqn{s \sim N(0, I_r)}, \eqn{Var(f_m) = I} and
#' \eqn{e_{mj} \sim N(0, 1 - ||\Lambda_{mj}||^2)}. The cross-view factor covariance
#' is \eqn{\sqrt{\gamma_m \gamma_{m'}} B_m B_{m'}^T}, so the crosstalk weights
#' \eqn{\gamma} change how strongly views talk to each other without touching
#' marginals or within-view structure. \eqn{\gamma = 1} reproduces the template's
#' shared structure and \eqn{\gamma = 0} gives independent views.
#'
#' @param template Named list of numeric matrices, samples in rows and features
#'   in columns. Rows are matched by rownames when present, otherwise by position.
#' @param rank Optional factors per view (scalar or one per view). Use
#'   \code{"2logp"} to retain \code{ceiling(2 * log(p_m))} factors in view
#'   \code{m}, subject to \code{max.rank}. If NULL, the smallest rank explaining
#'   \code{var.explained} of normal-score variance is selected, capped at
#'   \code{max.rank}.
#' @param var.explained Target variance fraction when \code{rank} is NULL.
#' @param max.rank Cap on the per-view rank (default floor(n0 / 3)).
#' @param shared.rank Number of shared variates r. If NULL, chosen by parallel
#'   analysis: generalized CCA eigenvalues above the 95th percentile of those
#'   obtained after permuting samples independently within each view. This
#'   filters out the spurious cross-correlation a small template always shows.
#' @param nperm Permutations for the parallel analysis.
#' @param view.transform Optional transform per view before fitting. Supported
#'   values are `"none"`, `"clr"`, and `"log1p"`. For microbiome data, use
#'   `"clr"` to fit and simulate on centered log-ratio scale.
#' @return An object of class \code{mimesys_template}.
#' @export
fit_template <- function(template, rank = NULL, var.explained = 0.7,
                         max.rank = NULL, shared.rank = NULL, nperm = 50,
                         view.transform = NULL) {

  rank_input <- rank
  template <- as_mimesys_view_list(template)
  template <- .mimesys_apply_view_transforms(template, view.transform = view.transform)
  view_transform <- attr(template, "view.transform")
  stopifnot(is.list(template), length(template) >= 1)
  if (is.null(names(template))) names(template) <- paste0("view", seq_along(template))
  views <- names(template)
  M <- length(views)
  template <- lapply(template, as.matrix)

  ## Match samples across views
  rn <- lapply(template, rownames)
  if (all(!vapply(rn, is.null, logical(1)))) {
    common <- Reduce(intersect, rn)
    if (length(common) < 10) stop("Fewer than 10 samples shared across views.")
    template <- lapply(template, function(x) x[common, , drop = FALSE])
  } else if (length(unique(vapply(template, nrow, integer(1)))) != 1) {
    stop("Views have different numbers of rows and no rownames to match on.")
  }
  n0 <- nrow(template[[1]])
  if (is.null(max.rank)) max.rank <- max(1, floor(n0 / 3))
  rank_2logp <- is.character(rank) && length(rank) == 1L && identical(rank, "2logp")
  if (is.character(rank) && !rank_2logp) {
    stop("rank must be NULL, numeric, or \"2logp\".", call. = FALSE)
  }
  if (!is.null(rank) && !rank_2logp) rank <- rep_len(rank, M)

  fits <- vector("list", M)
  observed <- vector("list", M)
  names(fits) <- views
  names(observed) <- views
  for (m in seq_len(M)) {
    X <- template[[m]]
    if (is.null(colnames(X))) colnames(X) <- paste0(views[m], "_", seq_len(ncol(X)))
    ## Median-impute missing values, drop constant features
    if (anyNA(X)) {
      for (j in which(colSums(is.na(X)) > 0)) X[is.na(X[, j]), j] <- stats::median(X[, j], na.rm = TRUE)
    }
    keep <- apply(X, 2, function(v) length(unique(v)) > 1)
    X <- X[, keep, drop = FALSE]
    ## Retain the exact empirical matrix used to fit this view. Crosstalk
    ## calibration compares null simulations against this observed structure.
    observed[[m]] <- X

    ## (i) Marginals: sorted values define the empirical quantile function
    qmat <- apply(X, 2, sort)
    discrete <- apply(X, 2, function(v) all(abs(v - round(v)) < 1e-8))

    ## Normal scores with random tie breaking (randomized PIT, handles zeros)
    Z <- apply(X, 2, function(v) stats::qnorm(rank(v, ties.method = "random") / (n0 + 1)))
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
    Lambda <- sv$v[, 1:k, drop = FALSE] %*% diag(sv$d[1:k] / sqrt(n0 - 1), k)
    h <- rowSums(Lambda^2)
    big <- h > 0.99
    if (any(big)) Lambda[big, ] <- Lambda[big, , drop = FALSE] * sqrt(0.99 / h[big])
    rownames(Lambda) <- colnames(X)

    fits[[m]] <- list(qmat = qmat, discrete = discrete, Lambda = Lambda,
                      Z = Z,
                      two_part = two_part,
                      scores = sv$u[, 1:k, drop = FALSE] * sqrt(n0 - 1),
                      rank = k, var.explained = ve[k],
                      raw.var = apply(X, 2, stats::var))
  }

  ## (iii) Shared cross-view space by MAXVAR generalized CCA. Scores are
  ## orthonormal within view, so the top eigenvectors of the correlation matrix
  ## of the stacked scores define the shared variates.
  B <- NULL; eig_obs <- NULL; thr <- NULL; r <- 0L
  if (M > 1) {
    Fall <- do.call(cbind, lapply(fits, `[[`, "scores"))
    eg <- eigen(stats::cor(Fall), symmetric = TRUE)
    eig_obs <- eg$values
    if (is.null(shared.rank)) {
      eig_null <- replicate(nperm, {
        Fp <- do.call(cbind, lapply(fits, function(f) f$scores[sample.int(n0), , drop = FALSE]))
        eigen(stats::cor(Fp), symmetric = TRUE, only.values = TRUE)$values
      })
      thr <- apply(eig_null, 1, stats::quantile, probs = 0.95)
      above <- eig_obs > thr
      r <- if (!above[1]) 0L else if (all(above)) length(above) else which(!above)[1] - 1L
    } else {
      r <- as.integer(shared.rank)
    }
    if (r > 0) {
      S <- scale(Fall %*% eg$vectors[, seq_len(r), drop = FALSE])
      S <- qr.Q(qr(S)) * sqrt(n0 - 1)                      # sample-orthonormal
      B <- lapply(fits, function(f) crossprod(f$scores, S) / (n0 - 1))   # k_m x r
      names(B) <- views
    }
  }

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
                 gcca.eigen = eig_obs, gcca.threshold = thr,
                 view.transform = view_transform, provenance = provenance),
            class = "mimesys_template")
}

#' @export
print.mimesys_template <- function(x, ...) {
  cat("MimeSys template\n")
  cat("  Views: ", paste(x$views, collapse = ", "), "\n", sep = "")
  cat("  Matched samples: ", x$n0, "\n", sep = "")
  ranks <- vapply(x$fits, function(f) f$rank, integer(1))
  features <- vapply(x$fits, function(f) nrow(f$Lambda), integer(1))
  cat("  Features: ", paste(names(features), features, sep = "=", collapse = ", "), "\n", sep = "")
  cat("  Within-view ranks: ", paste(names(ranks), ranks, sep = "=", collapse = ", "), "\n", sep = "")
  cat("  Shared rank: ", x$shared.rank, "\n", sep = "")
  cat("  Calibration-ready: ", !is.null(x$observed), "\n", sep = "")
  if (!is.null(x$provenance$package_version)) {
    cat("  MimeSys version: ", x$provenance$package_version, "\n", sep = "")
  }
  invisible(x)
}


#' Simulate null multiview data from a fitted template
#'
#' @param tmpl A \code{mimesys_template} from \code{fit_template()}.
#' @param n Number of samples to generate (need not match the template).
#' @param p Feature counts per view. Use \code{NULL} to keep template sizes, a
#'   single integer to use the same count for every view, an unnamed vector in
#'   template view order, or a named vector with one value per view. Fewer
#'   features than the template subsamples features; more adds offspring features
#'   (see \code{feature.jitter}).
#' @param crosstalk Crosstalk weights gamma, scalar or one per view. 0 gives
#'   independent views, 1 the template level, values above 1 amplify it up to
#'   \code{crosstalk_max()}. Unequal weights let some views carry more of the
#'   shared structure than others (the analogue of t_x != t_z in Ding et al. 2022).
#' @param feature.select Feature choice when p is below template size:
#'   "random" or "variance" (highest raw variance first).
#' @param feature.jitter Between 0 and 1. An offspring feature copies a parent's
#'   marginal and gets loadings sqrt(1 - jitter) * parent + sqrt(jitter) * another
#'   feature, so it resembles but does not duplicate its parent.
#' @param dependence Latent dependence simulator. `"factor"` uses the low-rank
#'   template factor model; `"empirical_copula"` uses the empirical latent feature
#'   correlation; `"module_copula"` smooths empirical latent correlations by
#'   within-view modules.
#' @param cor.shrinkage Shrinkage toward independence for copula dependence modes.
#' @param module.count,module.between.shrinkage Module-copula tuning parameters.
#' @param view.models Optional named vector of per-view observation models. Use
#'   `"two_part_copula"` for sparse compositional microbiome views; other views
#'   default to `"continuous"`.
#' @return A pcl-style list (feature_table features x samples, sample_metadata,
#'   feature_metadata) plus \code{truth}: factor scores, shared variates, crosstalk
#'   weights and implied latent canonical correlations. feature_metadata has
#'   \code{shared_loading}, the correlation of each feature's normal score with
#'   shared variate 1 at crosstalk = 1 (used for spike-in alignment).
#' @export
sim_from_template <- function(tmpl, n, p = NULL, crosstalk = 1,
                              feature.select = c("random", "variance"),
                              feature.jitter = 0.2,
                              dependence = c("factor", "empirical_copula", "module_copula"),
                              cor.shrinkage = 0.05,
                              module.count = NULL,
                              module.between.shrinkage = 0.75,
                              view.models = NULL) {

  stopifnot(inherits(tmpl, "mimesys_template"))
  feature.select <- match.arg(feature.select)
  dependence <- match.arg(dependence)
  views <- tmpl$views
  M <- length(views)
  p0 <- vapply(tmpl$fits, function(f) nrow(f$Lambda), integer(1))
  p <- .normalize_feature_counts(p, views, p0, arg_name = "p")
  view.models <- .normalize_view_models(view.models, views)

  ks <- vapply(tmpl$fits, function(f) f$rank, integer(1))
  block <- rep(seq_len(M), ks)
  r <- tmpl$shared.rank
  gamma <- setNames(rep_len(crosstalk, M), views)
  if (any(gamma < 0)) stop("crosstalk must be non-negative.")
  lim <- crosstalk_max(tmpl, per.view = TRUE)
  if (any(gamma > lim + 1e-10)) {
    stop("crosstalk exceeds the valid limit for: ",
         paste(sprintf("%s (max %.2f)", views[gamma > lim], lim[gamma > lim]), collapse = ", "))
  }
  sample_ids <- paste0("Subject", seq_len(n))

  X_list <- meta_list <- vector("list", M)
  parent_list <- lam_list <- offspring_list <- vector("list", M)
  for (m in seq_len(M)) {
    fit <- tmpl$fits[[m]]
    ## Template parents for the requested number of features
    if (p[m] <= p0[m]) {
      parent <- if (feature.select == "random") sort(sample.int(p0[m], p[m])) else
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
      Lam[h > 0.99, ] <- Lam[h > 0.99, , drop = FALSE] * sqrt(0.99 / h[h > 0.99])
      offspring <- c(rep(FALSE, p0[m]), rep(TRUE, length(extra)))
    }
    parent_list[[m]] <- parent
    lam_list[[m]] <- Lam
    offspring_list[[m]] <- offspring
  }

  latent_layout <- .build_latent_layout(tmpl, parent_list, p, view.models)

  if (dependence %in% c("empirical_copula", "module_copula")) {
    z_template <- latent_layout$template
    latent_cor <- stats::cor(z_template, use = "pairwise.complete.obs")
    latent_cor[!is.finite(latent_cor)] <- 0
    latent_cor <- .apply_crosstalk_to_latent_cor(latent_cor, latent_layout$latent_p, gamma)
    if (identical(dependence, "module_copula")) {
      latent_cor <- .module_denoise_correlation(
        latent_cor,
        p = latent_layout$latent_p,
        module.count = module.count,
        between.shrinkage = module.between.shrinkage
      )
    }
    latent_cor <- .shrink_correlation(latent_cor, cor.shrinkage)
    Z_all <- matrix(stats::rnorm(n * sum(latent_layout$latent_p)), n, sum(latent_layout$latent_p)) %*%
      .sqrt_psd(latent_cor)
    shared_out <- matrix(numeric(0), n, 0)
    scores_out <- Z_all
  } else {
    ## Latent factor scores: shared part plus view-specific remainder
    shared_out <- matrix(stats::rnorm(n * max(r, 1)), n, max(r, 1))[, seq_len(r), drop = FALSE]
    scores_out <- do.call(cbind, lapply(seq_len(M), function(m) {
      E <- matrix(stats::rnorm(n * ks[m]), n, ks[m])
      if (r == 0 || gamma[m] == 0) return(E)
      Bm <- sqrt(gamma[m]) * tmpl$B[[m]]
      shared_out %*% t(Bm) + E %*% .sqrt_psd(diag(ks[m]) - Bm %*% t(Bm))
    }))
  }

  latent_starts <- cumsum(c(1L, head(latent_layout$latent_p, -1L)))
  latent_ends <- cumsum(latent_layout$latent_p)
  for (m in seq_len(M)) {
    fit <- tmpl$fits[[m]]
    parent <- parent_list[[m]]
    Lam <- lam_list[[m]]
    offspring <- offspring_list[[m]]

    if (dependence %in% c("empirical_copula", "module_copula")) {
      Z <- Z_all[, latent_starts[m]:latent_ends[m], drop = FALSE]
    } else {
      psi <- 1 - rowSums(Lam^2)
      ## Normal-score data, then map through each parent's empirical quantiles
      Z <- scores_out[, block == m, drop = FALSE] %*% t(Lam) +
        sweep(matrix(stats::rnorm(n * p[m]), n, p[m]), 2, sqrt(psi), "*")
    }
    if (identical(view.models[[m]], "two_part_copula")) {
      if (dependence %in% c("empirical_copula", "module_copula")) {
        Z_presence <- Z[, seq_len(p[m]), drop = FALSE]
        Z_abundance <- Z[, p[m] + seq_len(p[m]), drop = FALSE]
      } else {
        Z_presence <- Z
        Z_abundance <- Z + matrix(stats::rnorm(n * p[m], sd = 0.15), n, p[m])
      }
      Xm <- .simulate_two_part_view(
        fit = fit,
        parent = parent,
        U_presence = stats::pnorm(Z_presence),
        U_abundance = stats::pnorm(Z_abundance)
      )
    } else {
      U <- stats::pnorm(Z)
      Xm <- matrix(vapply(seq_len(p[m]), function(j) {
        .emp_quantile(fit$qmat[, parent[j]], U[, j], fit$discrete[parent[j]])
      }, numeric(n)), n, p[m])
    }

    ids <- paste0(views[m], "_", seq_len(p[m]))
    dimnames(Xm) <- list(sample_ids, ids)
    X_list[[m]] <- Xm
    meta_list[[m]] <- data.frame(
      featureID = ids, featureType = views[m],
      template_feature = colnames(fit$qmat)[parent], offspring = offspring,
      view_model = view.models[[m]],
      shared_loading = if (r > 0) as.numeric(Lam %*% tmpl$B[[m]][, 1]) else 0,
      stringsAsFactors = FALSE
    )
  }

  feature_metadata <- do.call(rbind, meta_list)
  rownames(feature_metadata) <- feature_metadata$featureID

  list(feature_table = as.data.frame(t(do.call(cbind, X_list))),
       sample_metadata = data.frame(subjectID = sample_ids, row.names = sample_ids),
       feature_metadata = feature_metadata,
       truth = list(scores = scores_out, shared = shared_out, crosstalk = gamma,
                    cancor = .realized_cancor(tmpl, gamma),
                    dependence = dependence,
                    view.models = stats::setNames(view.models, views)))
}


#' Largest valid crosstalk weight
#'
#' The construction needs \eqn{\gamma_m B_m B_m^T \le I}, so the limit for view m
#' is \eqn{1 / \sigma_{max}(B_m)^2}, which is always at least 1.
#' @param tmpl A \code{mimesys_template}.
#' @param per.view Return one limit per view instead of the overall scalar limit.
#' @export
crosstalk_max <- function(tmpl, per.view = FALSE) {
  if (is.null(tmpl$B)) {
    lim <- setNames(rep(Inf, length(tmpl$views)), tmpl$views)
  } else {
    lim <- vapply(tmpl$B, function(b) 1 / max(svd(b, nu = 0, nv = 0)$d)^2, numeric(1))
  }
  if (per.view) lim else min(lim)
}


## ---- internal helpers ---------------------------------------------------

.normalize_view_models <- function(view.models, views) {
  valid <- c("continuous", "two_part_copula")
  if (is.null(view.models)) {
    return(rep("continuous", length(views)))
  }
  if (is.null(names(view.models)) || any(!nzchar(names(view.models)))) {
    if (length(view.models) != 1L) {
      stop("view.models must be named by view, or be a single shared model.", call. = FALSE)
    }
    view.models <- stats::setNames(rep(view.models[[1]], length(views)), views)
  }
  out <- vapply(views, function(view) {
    model <- if (view %in% names(view.models)) view.models[[view]] else "continuous"
    model <- tolower(gsub("-", "_", model))
    if (model %in% c("none", "default", "standard")) {
      model <- "continuous"
    }
    if (!model %in% valid) {
      stop("Unknown view model '", model, "' for view ", view, ".", call. = FALSE)
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
    stop(arg_name, " must be NULL, a numeric vector, or a list coercible to a numeric vector.", call. = FALSE)
  }
  if (length(p) == 1L && (is.null(names(p)) || !nzchar(names(p)))) {
    p <- stats::setNames(rep(p[[1]], length(views)), views)
  } else if (is.null(names(p)) || all(!nzchar(names(p)))) {
    if (length(p) != length(views)) {
      stop(
        arg_name, " must be length 1, length ", length(views),
        ", or named with the template view names: ", paste(views, collapse = ", "),
        call. = FALSE
      )
    }
    names(p) <- views
  }
  missing <- setdiff(views, names(p))
  if (length(missing) > 0L) {
    stop(arg_name, " is missing view(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  p <- p[views]
  if (any(!is.finite(p)) || any(p <= 0) || any(abs(p - round(p)) > 1e-8)) {
    stop(arg_name, " values must be positive finite integers.", call. = FALSE)
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
  presence_Z <- apply(positive, 2, function(v) {
    stats::qnorm(rank(as.numeric(v), ties.method = "random") / (length(v) + 1))
  })
  presence_Z <- scale(presence_Z)
  abundance_Z <- matrix(NA_real_, nrow(X), ncol(X))
  for (j in seq_len(ncol(X))) {
    z <- stats::rnorm(nrow(X))
    pos <- positive[, j]
    if (sum(pos) > 1L) {
      z[pos] <- stats::qnorm(rank(X[pos, j], ties.method = "random") / (sum(pos) + 1))
    }
    abundance_Z[, j] <- z
  }
  abundance_Z <- scale(abundance_Z)
  presence_Z[!is.finite(presence_Z)] <- 0
  abundance_Z[!is.finite(abundance_Z)] <- 0
  colnames(presence_Z) <- colnames(abundance_Z) <- colnames(X)
  list(
    zero_prob = zero_prob,
    prevalence = prevalence,
    positive_values = positive_values,
    positive_discrete = positive_discrete,
    presence_Z = presence_Z,
    abundance_Z = abundance_Z
  )
}

.build_latent_layout <- function(tmpl, parent_list, p, view.models) {
  pieces <- vector("list", length(tmpl$views))
  latent_p <- integer(length(tmpl$views))
  for (m in seq_along(tmpl$views)) {
    fit <- tmpl$fits[[m]]
    parent <- parent_list[[m]]
    if (identical(view.models[[m]], "two_part_copula")) {
      pieces[[m]] <- cbind(
        fit$two_part$presence_Z[, parent, drop = FALSE],
        fit$two_part$abundance_Z[, parent, drop = FALSE]
      )
      latent_p[[m]] <- 2L * p[[m]]
    } else {
      pieces[[m]] <- fit$Z[, parent, drop = FALSE]
      latent_p[[m]] <- p[[m]]
    }
  }
  list(template = do.call(cbind, pieces), latent_p = stats::setNames(latent_p, tmpl$views))
}

.simulate_two_part_view <- function(fit, parent, U_presence, U_abundance) {
  p <- length(parent)
  n <- nrow(U_presence)
  Xm <- matrix(0, n, p)
  for (j in seq_len(p)) {
    pj <- parent[[j]]
    present <- U_presence[, j] > fit$two_part$zero_prob[[pj]]
    if (any(present)) {
      Xm[present, j] <- .emp_quantile(
        fit$two_part$positive_values[[pj]],
        pmin(pmax(U_abundance[present, j], 1e-8), 1 - 1e-8),
        fit$two_part$positive_discrete[[pj]]
      )
    }
  }
  row_total <- rowSums(Xm)
  nonzero <- row_total > 0
  if (any(nonzero)) {
    Xm[nonzero, ] <- sweep(Xm[nonzero, , drop = FALSE], 1, row_total[nonzero], "/")
  }
  empty <- which(!nonzero)
  if (length(empty) > 0L) {
    for (i in empty) {
      j <- sample.int(p, 1)
      Xm[i, j] <- 1
    }
  }
  Xm
}

.sqrt_psd <- function(A) {
  e <- eigen(A, symmetric = TRUE)
  e$vectors %*% diag(sqrt(pmax(e$values, 0)), length(e$values)) %*% t(e$vectors)
}

.shrink_correlation <- function(R, shrinkage = 0.05) {
  shrinkage <- max(0, min(1, shrinkage))
  R[!is.finite(R)] <- 0
  R <- (R + t(R)) / 2
  diag(R) <- 1
  R <- (1 - shrinkage) * R + shrinkage * diag(nrow(R))
  e <- eigen(R, symmetric = TRUE)
  vals <- pmax(e$values, 1e-8)
  R <- e$vectors %*% diag(vals, length(vals)) %*% t(e$vectors)
  d <- sqrt(diag(R))
  R <- sweep(sweep(R, 1, d, "/"), 2, d, "/")
  R[!is.finite(R)] <- 0
  diag(R) <- 1
  R
}

.apply_crosstalk_to_latent_cor <- function(R, p, gamma) {
  p <- as.integer(p)
  starts <- cumsum(c(1L, head(p, -1L)))
  ends <- cumsum(p)
  for (a in seq_along(p)) {
    for (b in seq_along(p)) {
      if (a != b) {
        R[starts[a]:ends[a], starts[b]:ends[b]] <-
          sqrt(gamma[[a]] * gamma[[b]]) * R[starts[a]:ends[a], starts[b]:ends[b]]
      }
    }
  }
  diag(R) <- 1
  R
}

.module_denoise_correlation <- function(R, p, module.count = NULL, between.shrinkage = 0.75) {
  between.shrinkage <- max(0, min(1, between.shrinkage))
  p <- as.integer(p)
  starts <- cumsum(c(1L, head(p, -1L)))
  ends <- cumsum(p)

  for (m in seq_along(p)) {
    j <- starts[m]:ends[m]
    if (length(j) < 6L) {
      next
    }
    R_view <- R[j, j, drop = FALSE]
    R_view[!is.finite(R_view)] <- 0
    diag(R_view) <- 1
    k <- if (is.null(module.count)) {
      max(2L, min(12L, floor(sqrt(length(j)))))
    } else if (length(module.count) == 1L || is.null(names(module.count))) {
      as.integer(module.count[[1]])
    } else {
      as.integer(module.count[[m]])
    }
    k <- max(2L, min(k, floor(length(j) / 2)))
    hc <- stats::hclust(stats::as.dist(1 - abs(R_view)), method = "average")
    module <- stats::cutree(hc, k = k)
    for (a in unique(module)) {
      for (b in unique(module)) {
        if (a >= b) {
          next
        }
        ia <- which(module == a)
        ib <- which(module == b)
        block <- R_view[ia, ib, drop = FALSE]
        target <- mean(block[is.finite(block)], na.rm = TRUE)
        if (!is.finite(target)) {
          target <- 0
        }
        R_view[ia, ib] <- (1 - between.shrinkage) * block + between.shrinkage * target
        R_view[ib, ia] <- t(R_view[ia, ib, drop = FALSE])
      }
    }
    diag(R_view) <- 1
    R[j, j] <- R_view
  }
  R
}

## Latent canonical correlations between each pair of views implied by gamma
.realized_cancor <- function(tmpl, gamma) {
  views <- tmpl$views
  M <- length(views)
  gamma <- rep_len(gamma, M)
  out <- list()
  if (is.null(tmpl$B) || M < 2) return(out)
  for (a in 1:(M - 1)) for (b in (a + 1):M) {
    C <- sqrt(gamma[a] * gamma[b]) * tmpl$B[[a]] %*% t(tmpl$B[[b]])
    out[[paste(views[a], views[b], sep = ":")]] <- svd(C, nu = 0, nv = 0)$d[seq_len(tmpl$shared.rank)]
  }
  out
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
