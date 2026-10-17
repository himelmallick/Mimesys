## Build the small example datasets shipped in data/.
##
## Sources (kept in the GitHub repository under inst/extdata, excluded from the
## package build because of their size):
##   StelzerDOS.RData        Stelzer et al. (2021) multi-omics data on the onset
##                           of labor (CyTOF immunome, plasma metabolome and
##                           proteome), 3 timepoints for each of 53 subjects.
##   FranzosaE_2019_raw.rds  Franzosa et al. (2019) gut microbiome (taxonomic
##                           relative abundances) and stool metabolome profiles.
##
## Run from the package root.

library(MultiAssayExperiment)
library(SummarizedExperiment)

to_se <- function(x) {
    ## x: samples in rows, features in columns
    SummarizedExperiment(assays = list(counts = t(x)))
}

## ---- Stelzer: first timepoint per subject, 40 most variable features/view
e <- new.env()
load("inst/extdata/StelzerDOS.RData", envir = e)
pcl <- e$pcl
x <- t(as.matrix(pcl$feature_table))
view <- pcl$feature_metadata[[2]]
subject <- sub("_.*", "", rownames(x))
first <- !duplicated(subject)
x <- x[first, , drop = FALSE]
rownames(x) <- subject[first]
views <- lapply(split(seq_len(ncol(x)), view), function(j) {
    m <- x[, j, drop = FALSE]
    v <- apply(m, 2, stats::var)
    m[, order(v, decreasing = TRUE)[seq_len(40)], drop = FALSE]
})
stelzer_mae <- MultiAssayExperiment(
    experiments = ExperimentList(lapply(views, to_se)),
    colData = S4Vectors::DataFrame(
        subject = rownames(x),
        row.names = rownames(x)
    )
)

## ---- Franzosa: 40 most prevalent taxa and 40 most variable metabolites
f <- readRDS("inst/extdata/FranzosaE_2019_raw.rds")
taxa <- as.matrix(f$taxaAbundanceProfile[, -1])
rownames(taxa) <- f$taxaAbundanceProfile$sample_id
metab <- as.matrix(f$harmonizedMetaboliteProfile[, -1])
rownames(metab) <- f$harmonizedMetaboliteProfile$sample_id
metab <- metab[rownames(taxa), , drop = FALSE]
taxa <- taxa[, order(colMeans(taxa > 0), decreasing = TRUE)[seq_len(40)]]
metab <- metab[, colSums(is.na(metab)) == 0, drop = FALSE]
mv <- apply(log1p(metab), 2, stats::var)
metab <- metab[, order(mv, decreasing = TRUE)[seq_len(40)]]
meta <- f$MetadataProfile
rownames(meta) <- meta$sample_id
franzosa_mae <- MultiAssayExperiment(
    experiments = ExperimentList(list(
        Microbiome = to_se(taxa),
        Metabolome = to_se(metab)
    )),
    colData = S4Vectors::DataFrame(
        disease = meta[rownames(taxa), "diseaseCat"],
        row.names = rownames(taxa)
    )
)

save(stelzer_mae, file = "data/stelzer_mae.rda", compress = "xz")
save(franzosa_mae, file = "data/franzosa_mae.rda", compress = "xz")
