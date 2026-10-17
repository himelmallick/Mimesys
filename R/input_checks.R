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
