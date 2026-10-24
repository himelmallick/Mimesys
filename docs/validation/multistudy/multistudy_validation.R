## Validation of multi-study simulation in MimeSys (gen_multistudy()).
##
## Each heterogeneity setting is swept over a grid. Every grid level uses the
## same seeds (common random numbers), and MimeSys draws heterogeneity from a
## separate random stream, so across levels the base samples, common
## coefficients and noise are identical; only the heterogeneity changes.
## Expected: every target metric moves monotonically with its setting, and
## settings that should not affect a quantity leave it unchanged.
##
## Run:  Rscript multistudy_validation.R
## Options (environment variables): TEMPLATE ("stelzer", default;
## "franzosa": microbiome and metabolome, both with the two-part model; or
## "ibd": three studies, each from its own real cohort's template - PRISM,
## Validation and Lloyd-Price baseline from HuMMANet, read from DATA_DIR),
## NREP (replicates per level, default 20), CORES (default 1), OUT_DIR
## (default ./multistudy_validation_<template>), PLOT_ONLY=1 (redraw figures
## from saved results).

suppressPackageStartupMessages({
    library(MimeSys)
    library(glmnet)
    library(ggplot2)
    library(patchwork)
})

nrep <- as.integer(Sys.getenv("NREP", "20"))
cores <- as.integer(Sys.getenv("CORES", "1"))
template_name <- match.arg(Sys.getenv("TEMPLATE", "stelzer"), c("stelzer", "franzosa", "ibd"))
out_dir <- Sys.getenv("OUT_DIR", paste0("multistudy_validation_", template_name))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
papply <- function(X, FUN) {
    if (cores > 1L) parallel::mclapply(X, FUN, mc.cores = cores) else lapply(X, FUN)
}

nstudy <- 2L
if (template_name == "ibd") {
    ## One template per cohort: PRISM, Validation, Lloyd-Price baseline (earliest
    ## visit with paired data per subject). Taxa present in >= 10% of samples of
    ## every cohort; the 100 most abundant metabolites detected in >= 80% of
    ## samples of every cohort (harmonized metabolite names).
    data_dir <- Sys.getenv("DATA_DIR", "ibd_trial")
    rd <- function(f) suppressWarnings(readRDS(file.path(data_dir, f)))
    raw <- list(PRISM = rd("franzosae_2019.rds"), Validation = rd("franzosae_2019__22d0f995.rds"),
        LloydPrice = rd("lloydprice_2019.rds"))
    tomat <- function(df) {
        m <- as.matrix(df[, -1]); rownames(m) <- df$sample_id; storage.mode(m) <- "double"; m
    }
    ibd <- lapply(raw, function(x) list(Microbiome = tomat(x$taxaAbundanceProfile),
        Metabolome = tomat(x$harmonizedMetaboliteProfile)))
    lm_meta <- raw$LloydPrice$MetadataProfile
    subj <- sub("_[0-9]+$", "", lm_meta$Subject_ID)
    o <- order(subj, as.integer(sub("^.*_", "", lm_meta$Subject_ID)))
    keep <- lm_meta$sample_id[o][!duplicated(subj[o])]
    ibd$LloydPrice <- lapply(ibd$LloydPrice, function(x) x[keep, ])
    ct <- Reduce(intersect, lapply(ibd, function(x) colnames(x$Microbiome)))
    taxa <- ct[apply(sapply(ibd, function(x) colMeans(x$Microbiome[, ct] > 0)), 1, min) >= 0.10]
    cm <- Reduce(intersect, lapply(ibd, function(x) colnames(x$Metabolome)))
    det <- sapply(ibd, function(x) colMeans(!is.na(x$Metabolome[, cm]) & x$Metabolome[, cm] > 0))
    cand <- cm[apply(det, 1, min) >= 0.80]
    abund <- rowMeans(sapply(ibd, function(x) apply(log1p(x$Metabolome[, cand]), 2,
        stats::median, na.rm = TRUE)))
    metab <- cand[order(abund, decreasing = TRUE)][seq_len(min(100, length(cand)))]
    tm <- fit_study_templates(lapply(ibd, function(x) list(Microbiome = x$Microbiome[, taxa],
        Metabolome = x$Metabolome[, metab])), view.transform = c(Metabolome = "log1p"))
    print(tm)
    view_models <- c(Microbiome = "two_part_copula", Metabolome = "continuous")
    p_view <- NULL   # the templates' own features
    nstudy <- length(tm)
} else if (template_name == "stelzer") {
    data(stelzer_mae)
    tm <- fit_template(stelzer_mae, shared.rank = 2, seed = 1)
    view_models <- NULL
    p_view <- 60
} else {
    ## Sparse views (about 29% zeros each): two-part copula for both
    data(franzosa_mae)
    tm <- fit_template(franzosa_mae, view.transform = c(Metabolome = "log1p"),
        shared.rank = 2, seed = 1)
    view_models <- stats::setNames(rep("two_part_copula", length(tm$views)), tm$views)
    p_view <- 40   # the template's own features (no offspring)
}

## Studies of 150; study 1 is split 70/30, the others are held out (all test).
base_args <- list(template = tm, n = rep(150, nstudy), p = p_view, crosstalk = 1,
    signal.prob = 0.1, snr = 2, holdout.studies = seq_len(nstudy)[-1],
    train.fraction = 0.7, view.models = view_models)
study_pairs <- utils::combn(nstudy, 2, simplify = FALSE)
## Mean of f(i, j) over all pairs of studies
pair_mean <- function(f) mean(vapply(study_pairs, function(pr) f(pr[1], pr[2]), numeric(1)))
sim <- function(seed, ...) {
    do.call(gen_multistudy, utils::modifyList(base_args,
        list(seed = seed, ...)))
}

## ---- helpers ---------------------------------------------------------------
pooled_data <- function(s) {
    tr <- s$trainDat$Rep_1
    te <- s$testDat$Rep_1
    list(X = t(cbind(as.matrix(tr$feature_table), as.matrix(te$feature_table))),
        meta = rbind(tr$sample_metadata, te$sample_metadata),
        truth = tr$truth, fm = tr$feature_metadata,
        train = tr, test = te)
}
## Share of each feature's variance explained by study, averaged over features
study_r2 <- function(X, study) {
    X <- scale(X, scale = FALSE)
    g <- apply(X, 2, function(x) stats::ave(x, study))
    mean(colSums(g^2) / pmax(colSums(X^2), 1e-12))
}
## Within-study feature correlation, averaged over studies
within_cor <- function(X, study) {
    Reduce(`+`, lapply(split(seq_len(nrow(X)), study), function(i) {
        suppressWarnings(stats::cor(X[i, ]))
    })) / length(unique(study))
}
## Lasso trained on study 1's training samples; R^2 (squared correlation, so
## a study intercept does not count against transfer) on study 1's test
## samples (within) and on the held-out studies (transfer, averaged)
lasso_r2 <- function(d) {
    xtr <- t(as.matrix(d$train$feature_table))
    fit <- glmnet::cv.glmnet(xtr, d$train$sample_metadata$Y, nfolds = 5)
    xte <- t(as.matrix(d$test$feature_table))
    pred <- as.numeric(stats::predict(fit, xte, s = "lambda.min"))
    y <- d$test$sample_metadata$Y
    st <- d$test$sample_metadata$study
    r2 <- function(k) {
        p <- pred[st == k]
        if (stats::sd(p) == 0) return(0)
        stats::cor(p, y[st == k])^2
    }
    studies <- d$truth$studies$study
    c(within_R2 = r2(studies[1]),
        transfer_R2 = mean(vapply(studies[-1], r2, numeric(1))))
}
row <- function(experiment, level, rep, values) {
    data.frame(experiment = experiment, level = level, rep = rep,
        metric = names(values), value = unname(values))
}

plot_only <- identical(Sys.getenv("PLOT_ONLY"), "1")   # redraw from saved results
if (!plot_only) {
## ---- A. effect.cor (the old rho.beta) --------------------------------------
grid_rho <- seq(0, 1, by = 0.1)
res_A <- do.call(rbind, papply(seq_len(nrep), function(r) {
    do.call(rbind, lapply(grid_rho, function(rho) {
        s <- sim(1000 + r, effect.cor = rho)
        d <- pooled_data(s)
        B <- d$truth$beta
        cs <- d$fm$common_signal
        Xs <- scale(d$X)
        lp <- Xs %*% B
        row("A. effect.cor", rho, r, c(
            coef_cor = pair_mean(function(i, j) stats::cor(B[cs, i], B[cs, j])),
            predictor_cor = pair_mean(function(i, j) stats::cor(lp[, i], lp[, j])),
            lasso_r2(d)))
    }))
}))

## ---- B. support.overlap ----------------------------------------------------
grid_omega <- seq(0, 1, by = 0.1)
res_B <- do.call(rbind, papply(seq_len(nrep), function(r) {
    do.call(rbind, lapply(grid_omega, function(om) {
        s <- sim(2000 + r, support.overlap = om)
        d <- pooled_data(s)
        B <- d$truth$beta
        Xs <- scale(d$X)
        lp <- Xs %*% B
        row("B. support.overlap", om, r, c(
            support_jaccard = pair_mean(function(i, j) {
                sum(B[, i] != 0 & B[, j] != 0) / sum(B[, i] != 0 | B[, j] != 0)
            }),
            predictor_cor = pair_mean(function(i, j) stats::cor(lp[, i], lp[, j])),
            lasso_r2(d)))
    }))
}))

## ---- C-E. null features: case-mix shift and batch effects ------------------
grid_sd <- c(0, 0.25, 0.5, 0.75, 1, 1.5)
null_metrics <- function(experiment, arg, seed_base) {
    do.call(rbind, papply(seq_len(nrep), function(r) {
        ref <- NULL
        do.call(rbind, lapply(grid_sd, function(v) {
            s <- do.call(sim, c(list(seed = seed_base + r), stats::setNames(list(v), arg)))
            d <- pooled_data(s)
            st <- d$meta$study
            wc <- within_cor(d$X, st)
            sds <- vapply(split(seq_len(nrow(d$X)), st), function(i) {
                apply(d$X[i, ], 2, stats::sd)
            }, numeric(ncol(d$X)))
            ok <- apply(sds > 0, 1, all)
            studies <- colnames(sds)
            if (is.null(ref)) {
                ## Reference: an independent replicate at level 0 (sampling noise)
                s2 <- do.call(sim, c(list(seed = seed_base + r + 1e5),
                    stats::setNames(list(0), arg)))
                d2 <- pooled_data(s2)
                wc2 <- within_cor(d2$X, d2$meta$study)
                ref <<- list(wc = wc, xbeta = d$meta$Xbeta[order(d$meta$subjectID)],
                    noise = mean(abs((wc2 - wc)[upper.tri(wc)]), na.rm = TRUE))
            }
            row(experiment, v, r, c(
                study_R2 = study_r2(d$X, st),
                sd_log_ratio = pair_mean(function(i, j) mean(abs(log(sds[ok, i] / sds[ok, j])))),
                zero_frac_diff = pair_mean(function(i, j) {
                    mean(abs(colMeans(d$X[st == studies[i], ] == 0) -
                        colMeans(d$X[st == studies[j], ] == 0)))
                }),
                constant_features = sum(!ok),
                within_cor_change = mean(abs((wc - ref$wc)[upper.tri(wc)]), na.rm = TRUE),
                within_cor_noise = ref$noise,
                outcome_change = max(abs(d$meta$Xbeta[order(d$meta$subjectID)] - ref$xbeta))))
        }))
    }))
}
res_C <- null_metrics("C. shift.sd (case mix)", "shift.sd", 3000)
res_D <- null_metrics("D. batch.loc.sd", "batch.loc.sd", 4000)
res_E <- null_metrics("E. batch.scale.sd", "batch.scale.sd", 5000)

res <- rbind(res_A, res_B, res_C, res_D, res_E)
utils::write.csv(res, file.path(out_dir, "multistudy_validation_raw.csv"), row.names = FALSE)
} else {
    res <- utils::read.csv(file.path(out_dir, "multistudy_validation_raw.csv"))
}

## ---- summary and monotonicity --------------------------------------------
summ <- do.call(rbind, lapply(split(res, list(res$experiment, res$metric, res$level),
    drop = TRUE), function(x) {
    data.frame(experiment = x$experiment[1], metric = x$metric[1], level = x$level[1],
        mean = mean(x$value), se = stats::sd(x$value) / sqrt(nrow(x)))
}))
summ <- summ[order(summ$experiment, summ$metric, summ$level), ]
mono <- do.call(rbind, lapply(split(res, list(res$experiment, res$metric), drop = TRUE),
    function(x) {
        m <- tapply(x$value, x$level, mean)
        st <- suppressWarnings(stats::cor.test(x$level, x$value, method = "spearman"))
        data.frame(experiment = x$experiment[1], metric = x$metric[1],
            mean_at_min = m[[1]], mean_at_max = m[[length(m)]],
            means_monotone_increasing = all(diff(m) > 0),
            n_decreasing_steps = sum(diff(m) < 0),
            spearman_rho = unname(st$estimate), spearman_p = st$p.value)
    }))
rownames(mono) <- NULL
utils::write.csv(summ, file.path(out_dir, "multistudy_validation_summary.csv"), row.names = FALSE)
utils::write.csv(mono, file.path(out_dir, "multistudy_validation_monotonicity.csv"), row.names = FALSE)
print(mono, digits = 3)

## ---- figure ----------------------------------------------------------------
labels <- c(coef_cor = "Coefficient correlation\nbetween studies",
    predictor_cor = "Predictor correlation\n(study 1 vs study 2 coefficients)",
    support_jaccard = "Signal-feature overlap\n(Jaccard)",
    within_R2 = "Lasso R2, same study\n(control: should stay flat)",
    transfer_R2 = "Lasso R2, train study 1,\ntest held-out studies",
    study_R2 = "Feature variance\nexplained by study",
    sd_log_ratio = "Mean |log SD ratio|\nbetween studies",
    within_cor_change = "Change in within-study\ncorrelations",
    outcome_change = "Max change in outcome\npredictor vs setting 0",
    zero_frac_diff = "Mean |difference in zero\nfraction| between studies",
    constant_features = "Features constant\nwithin a study")
panel_row <- function(exp, metrics) {
    s <- summ[summ$experiment == exp & summ$metric %in% metrics, ]
    r <- res[res$experiment == exp & res$metric %in% metrics, ]
    s$metric_label <- factor(labels[s$metric], levels = labels[metrics])
    r$metric_label <- factor(labels[r$metric], levels = labels[metrics])
    noise <- summ[summ$experiment == exp & summ$metric == "within_cor_noise", ]
    g <- ggplot(s, aes(level, mean)) +
        geom_point(data = r, aes(level, value), colour = "grey75", size = 0.5,
            position = position_jitter(width = 0.01 * diff(range(s$level)), height = 0)) +
        geom_ribbon(aes(ymin = mean - 1.96 * se, ymax = mean + 1.96 * se),
            fill = "#2a6fdb", alpha = 0.25) +
        geom_line(colour = "#2a6fdb", linewidth = 0.8) +
        geom_point(colour = "#2a6fdb", size = 1.4) +
        facet_wrap(~metric_label, nrow = 1, scales = "free_y") +
        labs(title = exp, x = NULL, y = NULL) +
        theme_bw(base_size = 9) +
        theme(strip.background = element_blank(), panel.grid.minor = element_blank(),
            plot.title = element_text(face = "bold", size = 10))
    if (nrow(noise) > 0 && "within_cor_change" %in% metrics) {
        g <- g + geom_hline(data = data.frame(
            metric_label = factor(labels[["within_cor_change"]], levels = labels[metrics]),
            y = mean(noise$mean)), aes(yintercept = y), linetype = "dashed",
            colour = "#c0392b")
    }
    if ("coef_cor" %in% metrics) {
        g <- g + geom_abline(data = data.frame(metric_label = factor(labels[["coef_cor"]],
            levels = labels[metrics])), aes(intercept = 0, slope = 1),
            linetype = "dotted", colour = "grey30")
    }
    g
}
caption <- sprintf(paste0("Template: %s. Mean and 95%% CI over %d replicates (grey: replicates); same seeds at every level.\n",
    "Dotted: identity line. Red dashed: difference between two independent replicates."),
    template_name, nrep)
sig <- patchwork::wrap_plots(
    panel_row("A. effect.cor", c("coef_cor", "predictor_cor", "transfer_R2", "within_R2")),
    panel_row("B. support.overlap", c("support_jaccard", "predictor_cor", "transfer_R2", "within_R2")),
    ncol = 1) + patchwork::plot_annotation(caption = caption)
ggsave(file.path(out_dir, "signal_heterogeneity.png"), sig, width = 10, height = 5.6, dpi = 150)
nm <- c("study_R2", "sd_log_ratio", "within_cor_change", "outcome_change")
if (template_name %in% c("franzosa", "ibd")) nm <- c(nm, "zero_frac_diff", "constant_features")
nul <- patchwork::wrap_plots(
    panel_row("C. shift.sd (case mix)", nm),
    panel_row("D. batch.loc.sd", nm),
    panel_row("E. batch.scale.sd", nm),
    ncol = 1) + patchwork::plot_annotation(caption = caption)
ggsave(file.path(out_dir, "null_feature_heterogeneity.png"), nul,
    width = 2.5 * length(nm), height = 8, dpi = 150)
