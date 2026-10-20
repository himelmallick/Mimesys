## Internal argument checks shared by the exported functions.

.check_count <- function(x, name, min = 1L) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
        abs(x - round(x)) > 1e-8 || x < min) {
        stop(name, " must be a single integer of at least ", min, ".",
            call. = FALSE)
    }
    invisible(as.integer(round(x)))
}

.check_fraction <- function(x, name, lower = 0, upper = 1,
    include.lower = TRUE) {
    ok <- is.numeric(x) && length(x) == 1L && is.finite(x) && x <= upper &&
        (if (include.lower) x >= lower else x > lower)
    if (!ok) {
        stop(name, " must be a single number in ",
            if (include.lower) "[" else "(", lower, ", ", upper, "].",
            call. = FALSE)
    }
    invisible(x)
}

## Number(s) in [0, 1] (single = TRUE: exactly one number)
.check_unit_interval <- function(x, name, single = FALSE) {
    if (!is.numeric(x) || (single && length(x) != 1L) || any(!is.finite(x)) ||
        any(x < 0) || any(x > 1)) {
        stop(name, " must be ", if (single) "a single number" else "number(s)",
            " between 0 and 1.", call. = FALSE)
    }
    invisible(x)
}

## A single proportion strictly between 0 and 1 (prevalences)
.check_prevalence <- function(x, name) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0 || x >= 1) {
        stop(name, " must be strictly between 0 and 1.", call. = FALSE)
    }
    invisible(x)
}

## ---- shared helpers for the simulators --------------------------------------

## A fitted template, or fit one from raw input (arguments in `fit_args`)
.resolve_template <- function(template, fit_args = list(), args_name = "...") {
    if (inherits(template, "mimesys_template")) return(template)
    if (!is.list(fit_args)) {
        stop(args_name, " must be a named list.", call. = FALSE)
    }
    do.call(fit_template, c(list(template = template), fit_args))
}

## Intercept b such that the (weighted) mean of plogis(b + eta) equals the
## target prevalence
.logistic_intercept <- function(eta, prevalence, weights = NULL) {
    if (is.null(weights)) weights <- rep(1 / length(eta), length(eta))
    f <- function(b) sum(weights * stats::plogis(b + eta)) - prevalence
    stats::uniroot(f, interval = c(-100, 100), tol = 1e-10)$root
}

## Common set-up of the MAE wrappers: argument checks, null-model options,
## template and crosstalk. The caller sets the seed first.
.mimesys_wrapper_setup <- function(template, n, p, nrep, crosstalk, null_model,
    template_args, calibration, seed, null_preset = list(), min_n = 2L) {
    .check_count(n, "n", min = min_n)
    .check_count(nrep, "nrep", min = 1L)
    null_defaults <- utils::modifyList(list(
        feature_select = "random", feature_jitter = 0.2,
        view_models = NULL, back_transform = FALSE
    ), null_preset)
    null_model <- .mimesys_merge_options(null_model, null_defaults, "null_model")
    tmpl <- .resolve_template(template, template_args, "template_args")
    resolved <- .mimesys_resolve_crosstalk(
        tmpl, n = as.integer(n), p = p, crosstalk = crosstalk,
        calibration = calibration, null_model = null_model, seed = seed
    )
    list(template = tmpl, null_model = null_model, crosstalk = resolved)
}

## Provenance recorded by the MAE wrappers
.mimesys_wrapper_provenance <- function(tmpl, resolved, ...) {
    c(
        list(
            package_version = as.character(utils::packageVersion("MimeSys")),
            created_at = format(Sys.time(), tz = "UTC", usetz = TRUE)
        ),
        list(...),
        list(
            selected_crosstalk = resolved$value,
            calibration_grid = if (is.null(resolved$calibration)) NULL else
                resolved$calibration$provenance$crosstalk.grid,
            template = tmpl$provenance
        )
    )
}

## One MAE, or a mimesys_mae_simulation collection for several replicates
.mimesys_collect_maes <- function(maes, settings) {
    if (length(maes) == 1L) return(maes[[1]])
    structure(list(replicates = unname(maes), settings = settings),
        class = "mimesys_mae_simulation")
}
