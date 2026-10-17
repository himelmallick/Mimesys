## Sequential (staged) mediation: stages of views, links between the
## mediators of consecutive stages, and the latent construction that keeps
## every mediator's latent score standard normal.
##
## Notation (all on the latent normal scale). For mediator j with direct
## exposure weight a_j, stage-1 mediators are
##     W_j = a_j T* + sqrt(1 - a_j^2) Z_j,
## and stage-s mediators (s > 1) are
##     W_j = a_j T* + sum_k d_jk W_k + c_j Z_j,
## where k runs over the mediators of stage s - 1, Z_j is the mediator's own
## template latent score (correlated with the others through the template),
## and c_j > 0 makes Var(W_j) = 1. Each W_j is linear in T* and the Z's:
## W_j = mu_j T* + L_j Z, so given T it is normal with mean mu_j T* and
## variance 1 - mu_j^2.

## Stage index for every view (0 = no mediators)
.mediation_stage_of_view <- function(stages, views, n_med) {
    if (is.null(stages)) {
        out <- stats::setNames(ifelse(n_med > 0L, 1L, 0L), views)
        return(list(stage_of_view = out, n_stages = 1L, names = "stage1"))
    }
    if (is.character(stages)) {
        stop("stages must be a list (one element per stage), e.g. ",
            "list(\"Microbiome\", \"Metabolome\"); a character vector is ambiguous.",
            call. = FALSE)
    }
    if (!is.list(stages) || length(stages) < 1L ||
        !all(vapply(stages, is.character, logical(1)))) {
        stop("stages must be NULL or a list of character vectors of view names, ",
            "in causal order (e.g. list(\"Microbiome\", \"Metabolome\")).", call. = FALSE)
    }
    all_v <- unlist(stages, use.names = FALSE)
    unknown <- setdiff(all_v, views)
    if (length(unknown) > 0L) {
        stop("stages has unknown view name(s): ", paste(unknown, collapse = ", "),
            ". Views are: ", paste(views, collapse = ", "), ".", call. = FALSE)
    }
    if (anyDuplicated(all_v)) {
        stop("Each view can belong to only one stage; repeated: ",
            paste(unique(all_v[duplicated(all_v)]), collapse = ", "), ".", call. = FALSE)
    }
    out <- stats::setNames(integer(length(views)), views)
    for (s in seq_along(stages)) out[stages[[s]]] <- s
    outside <- views[n_med > 0L & out == 0L]
    if (length(outside) > 0L) {
        stop("View(s) ", paste(outside, collapse = ", "), " have mediators but are not ",
            "in any stage. Add them to stages or set their n.mediators to 0.", call. = FALSE)
    }
    out[n_med == 0L] <- 0L
    per_stage <- vapply(seq_along(stages), function(s) sum(n_med[out == s]), numeric(1))
    if (any(per_stage == 0)) {
        stop("Every stage needs at least one mediator; stage(s) ",
            paste(which(per_stage == 0), collapse = ", "), " have none.", call. = FALSE)
    }
    nm <- names(stages)
    if (is.null(nm) || any(!nzchar(nm))) {
        nm <- vapply(stages, paste, character(1), collapse = "+")
    }
    list(stage_of_view = out, n_stages = length(stages), names = nm)
}

## Defaults and checks for the link settings
.mediation_link_options <- function(links, n_stages) {
    defaults <- list(type = "random", density = 0.3, strength = 0.3,
        direction = 1, matrix = NULL)
    if (is.null(links)) links <- list()
    if (!is.list(links)) stop("stage.links must be a named list.", call. = FALSE)
    unknown <- setdiff(names(links), names(defaults))
    if (length(unknown) > 0L) {
        stop("Unknown stage.links option(s): ", paste(unknown, collapse = ", "),
            ". Options are: ", paste(names(defaults), collapse = ", "), ".", call. = FALSE)
    }
    links <- utils::modifyList(defaults, links)
    links$type <- match.arg(links$type, c("random", "hub", "convergent", "explicit"))
    .check_fraction(links$density, "stage.links$density", include.lower = FALSE)
    n_tr <- max(n_stages - 1L, 1L)
    st <- links$strength
    if (!is.numeric(st) || !(length(st) %in% c(1L, n_tr)) || any(!is.finite(st)) ||
        any(st < 0) || any(st >= 1)) {
        stop("stage.links$strength must be in [0, 1): one value, or one per link ",
            "between consecutive stages.", call. = FALSE)
    }
    links$strength <- rep_len(st, n_tr)
    .check_direction(links$direction, "stage.links$direction")
    if (identical(links$type, "explicit") && n_stages > 1L) {
        if (!is.list(links$matrix) || length(links$matrix) != n_stages - 1L) {
            stop("stage.links$type = \"explicit\" needs stage.links$matrix: a list with ",
                "one matrix per link between consecutive stages.", call. = FALSE)
        }
    }
    links
}

## Link pattern (rows: children in stage s, columns: parents in stage s - 1),
## signed, before scaling to the requested strength
.mediation_link_pattern <- function(n_child, n_parent, links, transition) {
    if (identical(links$type, "explicit")) {
        mat <- links$matrix[[transition]]
        if (!is.matrix(mat) || !is.numeric(mat) || !identical(dim(mat),
            c(as.integer(n_child), as.integer(n_parent))) || any(!is.finite(mat))) {
            stop("stage.links$matrix[[", transition, "]] must be a finite numeric ",
                n_child, " x ", n_parent, " matrix (rows: mediators of stage ",
                transition + 1L, ", columns: mediators of stage ", transition, ").",
                call. = FALSE)
        }
        return(mat)
    }
    pat <- matrix(0, n_child, n_parent)
    if (identical(links$type, "random")) {
        pat[] <- stats::runif(length(pat)) < links$density
        for (j in which(rowSums(pat) == 0)) pat[j, sample.int(n_parent, 1L)] <- 1
        for (k in which(colSums(pat) == 0)) pat[sample.int(n_child, 1L), k] <- 1
    } else if (identical(links$type, "hub")) {
        n_hub <- max(1L, min(n_parent, round(links$density * n_parent)))
        hubs <- .sample_values(seq_len(n_parent), n_hub)
        child_hub <- hubs[sample.int(n_hub, n_child, replace = TRUE)]
        first <- seq_len(min(n_hub, n_child))
        child_hub[sample.int(n_child)[first]] <- hubs[first]
        pat[cbind(seq_len(n_child), match(child_hub, seq_len(n_parent)))] <- 1
    } else {
        n_par <- min(n_parent, max(2L, round(links$density * n_parent)))
        for (j in seq_len(n_child)) pat[j, .sample_values(seq_len(n_parent), n_par)] <- 1
    }
    pat * matrix(.path_signs(links$direction, length(pat)), n_child, n_parent)
}

## Latent correlation between the template scores of the mediators (before
## any exposure effect), from the fitted loadings and crosstalk
.mediation_latent_cov <- function(tmpl, structure, gamma, med_view_idx, med_col) {
    K <- length(med_view_idx)
    lam <- structure$lam
    r <- tmpl$shared.rank
    load <- lapply(seq_len(K), function(j) lam[[med_view_idx[j]]][med_col[j], ])
    out <- diag(K)
    for (i in seq_len(K)) {
        for (j in seq_len(K)) {
            if (i >= j) next
            mi <- med_view_idx[i]
            mj <- med_view_idx[j]
            cf <- if (mi == mj) {
                diag(length(load[[i]]))
            } else if (r > 0) {
                sqrt(gamma[mi] * gamma[mj]) * tmpl$B[[mi]] %*% t(tmpl$B[[mj]])
            } else {
                matrix(0, length(load[[i]]), length(load[[j]]))
            }
            out[i, j] <- out[j, i] <- as.numeric(load[[i]] %*% cf %*% load[[j]])
        }
    }
    out
}

## Latent design for one scenario. D (scaled link weights, list by stage) is
## computed when `D` is NULL, from the link patterns, so that the parents
## explain `strength` of each child's latent variance; pass the D of the
## spiked scenario to keep the links identical across scenarios.
.mediation_design <- function(alpha, med_stage, n_stages, Sigma, patterns,
    strength, D = NULL) {
    K <- length(alpha)
    mu <- numeric(K)
    L <- matrix(0, K, K)
    cvec <- numeric(K)
    s1 <- which(med_stage == 1L)
    if (any(abs(alpha[s1]) >= 1)) stop("Exposure paths must be below 1 in size.", call. = FALSE)
    mu[s1] <- alpha[s1]
    cvec[s1] <- sqrt(1 - alpha[s1]^2)
    L[cbind(s1, s1)] <- cvec[s1]
    build_D <- is.null(D)
    if (build_D) D <- vector("list", n_stages)
    if (n_stages > 1L) {
        for (s in 2:n_stages) {
            ch <- which(med_stage == s)
            pa <- which(med_stage == s - 1L)
            if (build_D) {
                pat <- patterns[[s - 1L]]
                ## Covariance of the parents' latent scores
                R <- outer(mu[pa], mu[pa]) + L[pa, , drop = FALSE] %*% Sigma %*%
                    t(L[pa, , drop = FALSE])
                Ds <- pat
                for (j in seq_along(ch)) {
                    w <- pat[j, ]
                    v <- as.numeric(w %*% R %*% w)
                    Ds[j, ] <- if (v > 0 && strength[s - 1L] > 0) {
                        w * sqrt(strength[s - 1L] / v)
                    } else {
                        0
                    }
                }
                D[[s]] <- Ds
            }
            Ds <- D[[s]]
            for (j in seq_along(ch)) {
                cj <- ch[j]
                mu0 <- alpha[cj] + sum(Ds[j, ] * mu[pa])
                L0 <- as.numeric(Ds[j, , drop = FALSE] %*% L[pa, , drop = FALSE])
                v0 <- as.numeric(L0 %*% Sigma %*% L0)
                cross <- as.numeric(L0 %*% Sigma[, cj])
                disc <- cross^2 + 1 - mu0^2 - v0
                cc <- -cross + sqrt(max(disc, 0))
                if (disc <= 0 || cc <= 1e-3) {
                    stop("The exposure path and links into a stage-", s, " mediator explain ",
                        "all of its latent variance. Lower stage.links$strength, ",
                        "stage.direct or TIE.", call. = FALSE)
                }
                mu[cj] <- mu0
                L[cj, ] <- L0
                L[cj, cj] <- L[cj, cj] + cc
                cvec[cj] <- cc
            }
        }
    }
    list(alpha = alpha, mu = mu, c = cvec, D = D, med_stage = med_stage,
        n_stages = n_stages)
}

## Latent hook for staged mediation. `t_fun(r, s)` gives the standardized
## exposure (length n, or 1) that enters the stage-r mediators within the
## chain whose stage-s mediators are observed: each stage's mediators are
## built from their own copy of the chain, which allows the cross-world
## exposures needed for natural path-specific effects. With every
## t_fun(r, s) equal to the observed exposure this is the ordinary model.
.mediation_stage_hook <- function(design, med_view_idx, med_col, t_fun,
    support = NULL) {
    force(design)
    force(t_fun)
    force(support)
    function(Z_list) {
        n <- nrow(Z_list[[1]])
        K <- length(med_view_idx)
        Zm <- vapply(seq_len(K), function(j) Z_list[[med_view_idx[j]]][, med_col[j]],
            numeric(n))
        Zm <- matrix(Zm, n, K)
        st <- design$med_stage
        for (s in seq_len(design$n_stages)) {
            W <- matrix(0, n, K)
            for (q in seq_len(s)) {
                idx <- which(st == q)
                tq <- rep_len(t_fun(q, s), n)
                Wq <- outer(tq, design$alpha[idx]) +
                    sweep(Zm[, idx, drop = FALSE], 2, design$c[idx], "*")
                if (q > 1L) {
                    pa <- which(st == q - 1L)
                    Wq <- Wq + W[, pa, drop = FALSE] %*% t(design$D[[q]])
                }
                W[, idx] <- Wq
            }
            idx <- which(st == s)
            out <- W[, idx, drop = FALSE]
            if (!is.null(support)) {
                ## Binary exposure: W is a two-component normal mixture (mean
                ## mu * T*, variance 1 - mu^2 given T); map it through that
                ## mixture CDF and back so the factual latent score is exactly
                ## N(0, 1). The same fixed map is used in every counterfactual.
                for (j in seq_along(idx)) {
                    mj <- design$mu[idx[j]]
                    if (mj == 0) next
                    sdj <- sqrt(1 - mj^2)
                    u <- support$probs[1] * stats::pnorm((out[, j] - mj * support$values[1]) / sdj) +
                        support$probs[2] * stats::pnorm((out[, j] - mj * support$values[2]) / sdj)
                    out[, j] <- stats::qnorm(pmin(pmax(u, 1e-12), 1 - 1e-12))
                }
            }
            for (j in seq_along(idx)) {
                Z_list[[med_view_idx[idx[j]]]][, med_col[idx[j]]] <- out[, j]
            }
        }
        Z_list
    }
}

## Ordered paths (r, s): exposure enters at stage r and reaches the outcome
## through the mediators of stage s. They are switched on one at a time in
## this order to decompose the natural indirect effect: later-starting paths
## first, as in the usual sequential decomposition (e.g. for two stages:
## T -> M2 -> Y, then T -> M1 -> M2 -> Y, then T -> M1 -> Y).
.mediation_paths <- function(n_stages) {
    out <- list()
    for (s in rev(seq_len(n_stages))) {
        for (r in rev(seq_len(s))) out[[length(out) + 1L]] <- c(r = r, s = s)
    }
    out
}

.mediation_path_label <- function(r, s, stage_names) {
    paste(c("T", stage_names[r:s], "Y"), collapse = " -> ")
}
