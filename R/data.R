#' Stelzer et al. multi-omics example data
#'
#' A small `MultiAssayExperiment` derived from the multi-omics study of the
#' onset of labor by Stelzer et al. (2021): one timepoint (the first) for each
#' of 53 subjects and the 40 most variable features of each of three views.
#' Assays are features by samples.
#'
#' @format A `MultiAssayExperiment` with three experiments
#'   (`SummarizedExperiment`, assay `counts`):
#'   \describe{
#'     \item{CYTOF}{40 single-cell mass cytometry immune features.}
#'     \item{Metabolomics}{40 plasma metabolite intensities.}
#'     \item{Proteomics}{40 plasma protein measurements.}
#'   }
#'   `colData` holds the subject identifier.
#' @source Stelzer IA, et al. (2021). Integrated trajectories of the maternal
#'   metabolome, proteome, and immunome predict labor onset. *Science
#'   Translational Medicine* 13(592):eabd9898. The script that builds the
#'   object is in `inst/scripts/make-example-data.R`.
#' @return A `MultiAssayExperiment` (loaded with `data(stelzer_mae)`).
#' @examples
#' data(stelzer_mae)
#' stelzer_mae
"stelzer_mae"

#' Franzosa et al. microbiome and metabolome example data
#'
#' A small `MultiAssayExperiment` derived from the gut microbiome and stool
#' metabolome study of inflammatory bowel disease by Franzosa et al. (2019):
#' 155 samples, the 40 most prevalent microbial species (relative abundances)
#' and the 40 most variable metabolites. Assays are features by samples.
#'
#' @format A `MultiAssayExperiment` with two experiments
#'   (`SummarizedExperiment`, assay `counts`):
#'   \describe{
#'     \item{Microbiome}{40 species relative abundances (sparse,
#'   zero-inflated).}
#'     \item{Metabolome}{40 metabolite abundances.}
#'   }
#'   `colData` holds the disease group (`CD`, `UC`, `Control`).
#' @source Franzosa EA, et al. (2019). Gut microbiome structure and metabolic
#'   activity in inflammatory bowel disease. *Nature Microbiology*
#'   4(2):293-305. The script that builds the object is in
#'   `inst/scripts/make-example-data.R`.
#' @return A `MultiAssayExperiment` (loaded with `data(franzosa_mae)`).
#' @examples
#' data(franzosa_mae)
#' franzosa_mae
"franzosa_mae"
