#' MimeSys: template-based simulation of multi-omics data
#'
#' MimeSys fits a template to an empirical multi-omics dataset and simulates
#' new datasets that keep each feature's marginal distribution and the
#' within- and cross-view dependence of the template, with controllable sample
#' size, feature counts, cross-view dependence (crosstalk) and outcome signal.
#'
#' The main entry points are [fit_template()], [simulate_mimesys()] (returns a
#' `MultiAssayExperiment`), [gen_simmba()] (lower-level replicate lists),
#' [calibrate_crosstalk()], [validate_template()], and, for mediation,
#' [gen_mediation()] and [simulate_mediation()]. See
#' `vignette("mimesys-workflow", package = "MimeSys")`.
#'
#' @keywords internal
"_PACKAGE"
