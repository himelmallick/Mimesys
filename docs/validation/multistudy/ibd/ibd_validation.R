## Validation of per-study templates in gen_multistudy() on three IBD cohorts
## from HuMMANet: FranzosaE_2019 (PRISM), FranzosaE_2019_Validation and
## LloydPrice_2019 (HMP2; baseline = earliest visit with paired data per
## subject).
##
## Questions:
##   1. Fidelity: is each simulated study closer to its own cohort than to the
##      other cohorts (marginals, zero fractions, correlation structure)?
##   2. Between-cohort differences: do the real differences in feature means
##      between cohorts carry over to the simulated studies?
##   3. Cohort identity: does a classifier trained on the real cohorts assign
##      simulated samples to the right cohort (and vice versa)?
## A single-template simulation (every study from the PRISM template) is the
## comparator for 2 and 3.
##
## Run: Rscript ibd_validation.R   (DATA_DIR = folder with the three HuMMANet
## .rds files; NREP replicates, default 20; OUT_DIR output folder)

suppressPackageStartupMessages({
    library(MimeSys)
    library(glmnet)
    library(ggplot2)
    library(patchwork)
})
data_dir <- Sys.getenv("DATA_DIR", "ibd_trial")
nrep <- as.integer(Sys.getenv("NREP", "20"))
out_dir <- Sys.getenv("OUT_DIR", "ibd_validation")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---- data ----------------------------------------------------------------
rd <- function(f) suppressWarnings(readRDS(file.path(data_dir, f)))
raw <- list(PRISM = rd("franzosae_2019.rds"),
    Validation = rd("franzosae_2019__22d0f995.rds"),
    LloydPrice = rd("lloydprice_2019.rds"))
tomat <- function(df) {
    m <- as.matrix(df[, -1])
    rownames(m) <- df$sample_id
    storage.mode(m) <- "double"
    m
}
ibd <- lapply(raw, function(x) list(Microbiome = tomat(x$taxaAbundanceProfile),
    Metabolome = tomat(x$harmonizedMetaboliteProfile), meta = x$MetadataProfile))
## Lloyd-Price baseline: earliest available visit per subject (Subject_ID is
## <subject>_<visit>)
m <- ibd$LloydPrice$meta
subj <- sub("_[0-9]+$", "", m$Subject_ID)
vis <- as.integer(sub("^.*_", "", m$Subject_ID))
o <- order(subj, vis)
keep <- m$sample_id[o][!duplicated(subj[o])]
for (v in c("Microbiome", "Metabolome")) ibd$LloydPrice[[v]] <- ibd$LloydPrice[[v]][keep, ]

## Common features: taxa present in >= 10% of samples of every cohort; the 100
## most abundant metabolites detected in >= 80% of samples of every cohort
ct <- Reduce(intersect, lapply(ibd, function(x) colnames(x$Microbiome)))
taxa <- ct[apply(sapply(ibd, function(x) colMeans(x$Microbiome[, ct] > 0)), 1, min) >= 0.10]
cm <- Reduce(intersect, lapply(ibd, function(x) colnames(x$Metabolome)))
det <- sapply(ibd, function(x) colMeans(!is.na(x$Metabolome[, cm]) & x$Metabolome[, cm] > 0))
cand <- cm[apply(det, 1, min) >= 0.80]
abund <- rowMeans(sapply(ibd, function(x) apply(log1p(x$Metabolome[, cand]), 2,
    stats::median, na.rm = TRUE)))
metab <- cand[order(abund, decreasing = TRUE)][seq_len(min(100, length(cand)))]
cohorts <- lapply(ibd, function(x) list(Microbiome = x$Microbiome[, taxa],
    Metabolome = x$Metabolome[, metab]))
study_n <- vapply(cohorts, function(x) nrow(x$Microbiome), integer(1))
cat("Cohorts:", paste(names(study_n), study_n, collapse = ", "), "| taxa", length(taxa),
    "| metabolites", length(metab), "\n")

## ---- templates -------------------------------------------------------------
tl <- fit_study_templates(cohorts, view.transform = c(Metabolome = "log1p"))
print(tl)
vm <- c(Microbiome = "two_part_copula", Metabolome = "continuous")

## Real data on the fitted scale (taxa proportions, log1p metabolites); the
## template's feature order
real <- lapply(cohorts, function(x) {
    cbind(x$Microbiome[, colnames(tl[[1]]$fits$Microbiome$qmat)],
        log1p(x$Metabolome[, colnames(tl[[1]]$fits$Metabolome$qmat)]))
})
nt <- ncol(tl[[1]]$fits$Microbiome$qmat)
is_taxa <- seq_len(ncol(real[[1]])) <= nt

sim_once <- function(seed, template) {
    s <- gen_multistudy(n = study_n, template = template, view.models = vm,
        crosstalk = 1, seed = seed)
    X <- rbind(t(as.matrix(s$trainDat$Rep_1$feature_table)),
        t(as.matrix(s$testDat$Rep_1$feature_table)))
    st <- c(s$trainDat$Rep_1$sample_metadata$study, s$testDat$Rep_1$sample_metadata$study)
    stats::setNames(lapply(levels(factor(st, levels = unique(st))), function(k) X[st == k, ]),
        unique(st))
}

## ---- distances between two data sets ---------------------------------------
ks_mean <- function(A, B, cols) {
    mean(vapply(which(cols), function(j) {
        a <- A[, j][is.finite(A[, j])]
        b <- B[, j][is.finite(B[, j])]
        suppressWarnings(stats::ks.test(a, b)$statistic)
    }, numeric(1)))
}
scor <- function(X) suppressWarnings(stats::cor(X, method = "spearman", use = "pairwise.complete.obs"))
cor_dist <- function(CA, CB, block) {
    d <- abs(CA - CB)[block]
    mean(d, na.rm = TRUE)
}
blocks <- local({
    p <- length(is_taxa)
    ut <- upper.tri(matrix(0, p, p))
    list(within_taxa = ut & outer(is_taxa, is_taxa, "&"),
        within_metab = ut & outer(!is_taxa, !is_taxa, "&"),
        cross_view = outer(is_taxa, !is_taxa, "&"))
})
dists <- function(A, B, CA = scor(A), CB = scor(B)) {
    c(KS_taxa = ks_mean(A, B, is_taxa), KS_metab = ks_mean(A, B, !is_taxa),
        zero_taxa = mean(abs(colMeans(A[, is_taxa] == 0, na.rm = TRUE) -
            colMeans(B[, is_taxa] == 0, na.rm = TRUE))),
        cor_within_taxa = cor_dist(CA, CB, blocks$within_taxa),
        cor_within_metab = cor_dist(CA, CB, blocks$within_metab),
        cor_cross_view = cor_dist(CA, CB, blocks$cross_view))
}
real_cor <- lapply(real, scor)

## ---- 1. fidelity: real cohort x simulated study ----------------------------
fid <- do.call(rbind, lapply(seq_len(nrep), function(r) {
    sims <- sim_once(1000 + r, tl)
    do.call(rbind, lapply(names(real), function(a) do.call(rbind, lapply(names(sims), function(b) {
        v <- dists(real[[a]], sims[[b]], real_cor[[a]], scor(sims[[b]]))
        data.frame(rep = r, real = a, sim = b, metric = names(v), value = unname(v))
    }))))
}))
## Sampling-noise floor: two independent simulations of the same cohort (same
## n) differ by this much by chance alone
floor_res <- do.call(rbind, lapply(seq_len(nrep), function(r) {
    s1 <- sim_once(4000 + r, tl)
    s2 <- sim_once(5000 + r, tl)
    do.call(rbind, lapply(names(s1), function(a) {
        v <- dists(s1[[a]], s2[[a]])
        data.frame(rep = r, real = a, metric = names(v), value = unname(v))
    }))
}))
real_real <- do.call(rbind, lapply(names(real), function(a) do.call(rbind,
    lapply(setdiff(names(real), a), function(b) {
        v <- dists(real[[a]], real[[b]], real_cor[[a]], real_cor[[b]])
        data.frame(real = a, other = b, metric = names(v), value = unname(v))
    }))))

## ---- 2. between-cohort differences ----------------------------------------
smd <- function(A, B) {
    ma <- colMeans(A, na.rm = TRUE); mb <- colMeans(B, na.rm = TRUE)
    sp <- sqrt((apply(A, 2, stats::var, na.rm = TRUE) + apply(B, 2, stats::var, na.rm = TRUE)) / 2)
    ifelse(sp > 0, (ma - mb) / sp, 0)
}
pairs <- utils::combn(names(real), 2, simplify = FALSE)
real_smd <- lapply(pairs, function(pr) smd(real[[pr[1]]], real[[pr[2]]]))
single <- stats::setNames(rep(list(tl$PRISM), 3), names(tl))   # comparator
diff_res <- do.call(rbind, lapply(seq_len(nrep), function(r) {
    out <- list()
    for (kind in c("Per-study templates", "One template")) {
        sims <- sim_once(2000 + r, if (kind == "One template") single else tl)
        for (i in seq_along(pairs)) {
            pr <- pairs[[i]]
            s <- smd(sims[[pr[1]]], sims[[pr[2]]])
            out[[length(out) + 1]] <- data.frame(rep = r, kind = kind,
                pair = paste(pr, collapse = " vs "),
                r_taxa = stats::cor(real_smd[[i]][is_taxa], s[is_taxa]),
                r_metab = stats::cor(real_smd[[i]][!is_taxa], s[!is_taxa]))
        }
    }
    do.call(rbind, out)
}))
smd_example <- local({
    sims <- sim_once(2001, tl)
    do.call(rbind, lapply(seq_along(pairs), function(i) {
        pr <- pairs[[i]]
        data.frame(pair = paste(pr, collapse = " vs "), real = real_smd[[i]],
            sim = smd(sims[[pr[1]]], sims[[pr[2]]]),
            view = ifelse(is_taxa, "Microbiome", "Metabolome"))
    }))
})

## ---- 3. cohort classifier ----------------------------------------------------
prep <- function(X, med) {
    X[, is_taxa] <- asin(sqrt(pmin(pmax(X[, is_taxa], 0), 1)))
    for (j in which(colSums(is.na(X)) > 0)) X[is.na(X[, j]), j] <- med[j]
    X
}
lab_real <- rep(names(real), vapply(real, nrow, integer(1)))
Xr <- do.call(rbind, real)
med <- apply(Xr, 2, stats::median, na.rm = TRUE)
Xr <- prep(Xr, med)
set.seed(1)
fit_real <- glmnet::cv.glmnet(Xr, factor(lab_real), family = "multinomial", nfolds = 5)
cls <- do.call(rbind, lapply(seq_len(min(nrep, 10)), function(r) {
    out <- list()
    for (kind in c("Per-study templates", "One template")) {
        sims <- sim_once(3000 + r, if (kind == "One template") single else tl)
        Xs <- prep(do.call(rbind, sims), med)
        lab_sim <- rep(names(sims), vapply(sims, nrow, integer(1)))
        acc_rs <- mean(predict(fit_real, Xs, s = "lambda.min", type = "class") == lab_sim)
        fit_sim <- glmnet::cv.glmnet(Xs, factor(lab_sim), family = "multinomial", nfolds = 5)
        acc_sr <- mean(predict(fit_sim, Xr, s = "lambda.min", type = "class") == lab_real)
        out[[length(out) + 1]] <- data.frame(rep = r, kind = kind,
            train_real_test_sim = acc_rs, train_sim_test_real = acc_sr)
    }
    do.call(rbind, out)
}))
chance <- max(table(lab_real)) / length(lab_real)

saveRDS(list(fid = fid, floor = floor_res, real_real = real_real, diff = diff_res, smd_example = smd_example,
    cls = cls, chance = chance, study_n = study_n, ranks = lapply(tl, `[[`, "shared.rank")),
    file.path(out_dir, "ibd_validation_results.rds"))

## ---- summaries ---------------------------------------------------------------
own <- aggregate(value ~ real + metric, fid[fid$real == fid$sim, ], mean)
other_sim <- aggregate(value ~ real + metric, fid[fid$real != fid$sim, ], mean)
other_real <- aggregate(value ~ real + metric, real_real, mean)
summ <- merge(merge(own, other_sim, by = c("real", "metric"), suffixes = c("_own_template", "_other_templates")),
    stats::setNames(other_real, c("real", "metric", "value_other_real_cohorts")), by = c("real", "metric"))
print(summ, digits = 3)
print(aggregate(cbind(r_taxa, r_metab) ~ kind + pair, diff_res, mean), digits = 3)
print(aggregate(cbind(train_real_test_sim, train_sim_test_real) ~ kind, cls, mean), digits = 3)
cat("chance (majority class):", round(chance, 3), "\n")
utils::write.csv(summ, file.path(out_dir, "fidelity_summary.csv"), row.names = FALSE)
