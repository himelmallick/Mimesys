#' MimeSys: template-based simulation of multi-omics data
#'
#' MimeSys fits a template to an empirical multi-omics dataset and simulates
#' new datasets that keep each feature's marginal distribution and the
#' within- and cross-view dependence of the template, with controllable sample
#' size, feature counts, cross-view dependence (crosstalk) and outcome signal.
#'
#' The main entry points are [fit_template()] and [simulate_mimesys()], which
#' runs every simulation design (outcome, null, multi-study and mediation) and
#' returns a `MultiAssayExperiment`; [calibrate_crosstalk()] and
#' [validate_template()] check the template. The engines [gen_simmba()],
#' [gen_multistudy()], [gen_mediation()] and [sim_from_template()] return
#' lower-level replicate lists. See
#' `vignette("mimesys-workflow", package = "MimeSys")`.
#'
#' @keywords internal
"_PACKAGE"
