#' Coerce multi-omics data to MimeSys view matrices
#'
#' Converts common multi-omics containers to the internal MimeSys
#' representation:
#' a named list of numeric matrices with samples in rows and features in
#' columns.
#' `MultiAssayExperiment` and `SummarizedExperiment` assays are transposed from
#' the Bioconductor convention of features in rows and samples in columns.
#'
#' @param x A `MultiAssayExperiment`, a `SummarizedExperiment`, or a named list
#'   of matrices/data frames/`SummarizedExperiment` objects.
#' @param assay Assay name or index to extract from each `SummarizedExperiment`.
#'   If `NULL`, the first assay is used.
#' @param single_view_name Name to use when `x` is a single
#'   `SummarizedExperiment`.
#'
#' @return A named list of numeric sample-by-feature matrices.
#' @examples
#' data(stelzer_mae)
#' views <- as_mimesys_view_list(stelzer_mae)
#' vapply(views, dim, integer(2))
#' @export
as_mimesys_view_list <- function(x, assay = NULL, single_view_name = "view1") {
    if (inherits(x, "MultiAssayExperiment")) {
        if (!requireNamespace("MultiAssayExperiment", quietly = TRUE)) {
            stop("MultiAssayExperiment support requires the MultiAssayExperiment package.",
                call. = FALSE)
        }
        mae <- x
        x <- as.list(MultiAssayExperiment::experiments(mae))
        ## Match samples across experiments by primary ID: rename each
        ## experiment's columns to their primary IDs through the sampleMap
        smap <- MultiAssayExperiment::sampleMap(mae)
        for (a in names(x)) {
            rows <- as.character(smap$assay) == a
            primary <- stats::setNames(as.character(smap$primary[rows]),
                as.character(smap$colname[rows]))
            cn <- colnames(x[[a]])
            if (is.null(cn) || !all(cn %in% names(primary))) {
                warning("Experiment ", a, " has columns missing from the sampleMap; ",
                    "matching its samples by column name instead of primary ID.",
                    call. = FALSE)
                next
            }
            new <- unname(primary[cn])
            if (anyDuplicated(new)) {
                stop("Experiment ", a, " has several columns for the same primary sample (",
                    paste(utils::head(unique(new[duplicated(new)]), 3), collapse = ", "),
                    "). Keep one column per primary sample.", call. = FALSE)
            }
            colnames(x[[a]]) <- new
        }
        ## Experiments stored as plain matrices follow the Bioconductor
        ## convention too (features in rows): transpose them like assays
        x <- lapply(x, function(e) {
            if (inherits(e, "SummarizedExperiment")) e else t(as.matrix(e))
        })
    } else if (inherits(x, "SummarizedExperiment")) {
        x <- stats::setNames(list(x), single_view_name)
    }

    if (!is.list(x) || length(x) == 0L) {
        stop("x must be a MultiAssayExperiment, SummarizedExperiment, or named list.",
            call. = FALSE)
    }
    if (is.null(names(x)) || any(!nzchar(names(x)))) {
        names(x) <- paste0("view", seq_along(x))
    }

    out <- lapply(seq_along(x), function(i) {
        view <- names(x)[[i]]
        .mimesys_one_view_matrix(x[[i]], assay = assay, view = view)
    })
    names(out) <- make.unique(names(x), sep = "_")
    ## Views without sample IDs are matched by position, which is only
    ## possible when every view has the same number of samples
    generated <- vapply(out, function(m) isTRUE(attr(m, "mimesys_generated_ids")),
        logical(1))
    out <- lapply(out, function(m) {
        attr(m, "mimesys_generated_ids") <- NULL
        m
    })
    sizes <- vapply(out, nrow, integer(1))
    if (any(generated) && length(unique(sizes)) > 1L) {
        stop("View(s) ", paste(names(out)[generated], collapse = ", "),
            " have no sample IDs (row names), and the views have different numbers ",
            "of samples (", paste(sprintf("%s: %d", names(out), sizes), collapse = ", "),
            "), so samples cannot be matched across views. Add sample IDs as row ",
            "names (or use a MultiAssayExperiment).", call. = FALSE)
    }
    out
}

.mimesys_one_view_matrix <- function(x, assay = NULL, view = "view") {
    if (inherits(x, "SummarizedExperiment")) {
        if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
            stop("SummarizedExperiment support requires the SummarizedExperiment package.",
                call. = FALSE)
        }
        assay_names <- SummarizedExperiment::assayNames(x)
        assay_id <- assay
        if (is.null(assay_id)) {
            assay_id <- 1L
        } else if (is.character(assay_id)) {
            hit <- match(assay_id, assay_names)
            if (is.na(hit)) {
                stop("Assay '", assay_id, "' not found for view ", view, ".",
                    call. = FALSE)
            }
            assay_id <- hit
        }
        mat <- t(as.matrix(SummarizedExperiment::assay(x, assay_id)))
    } else {
        mat <- as.matrix(x)
    }

    storage.mode(mat) <- "double"
    if (!is.numeric(mat)) {
        stop("View ", view, " must be coercible to a numeric matrix.",
            call. = FALSE)
    }
    if (nrow(mat) < 2L || ncol(mat) < 1L) {
        stop("View ", view,
            " must contain at least two samples and one feature.", call. = FALSE)
    }
    if (is.null(rownames(mat))) {
        rownames(mat) <- paste0("sample_", seq_len(nrow(mat)))
        attr(mat, "mimesys_generated_ids") <- TRUE
    }
    if (is.null(colnames(mat))) {
        colnames(mat) <- paste0(view, "_feature_", seq_len(ncol(mat)))
    }
    colnames(mat) <- make.unique(as.character(colnames(mat)), sep = "_")
    mat
}

.mimesys_apply_view_transforms <- function(matrices, view.transform = NULL) {
    if (is.null(view.transform)) {
        attr(matrices, "view.transform") <- stats::setNames(rep("none",
            length(matrices)), names(matrices))
        attr(matrices, "transform.info") <- stats::setNames(
            lapply(matrices, function(x) list(type = "none")), names(matrices)
        )
        return(matrices)
    }
    unknown <- setdiff(names(view.transform), c("", names(matrices)))
    if (length(unknown) > 0L) {
        stop("view.transform has unknown view name(s): ", paste(unknown, collapse = ", "),
            ". Views are: ", paste(names(matrices), collapse = ", "), ".", call. = FALSE)
    }
    if (is.null(names(view.transform)) || any(!nzchar(names(view.transform)))) {
        if (length(view.transform) != 1L) {
            stop("view.transform must be named by view, or be a single shared transform.",
                call. = FALSE)
        }
        view.transform <- stats::setNames(rep(view.transform[[1]],
            length(matrices)), names(matrices))
    }

    applied <- stats::setNames(rep("none", length(matrices)), names(matrices))
    info <- stats::setNames(vector("list", length(matrices)), names(matrices))
    for (view in names(matrices)) {
        transform <- if (view %in% names(view.transform)) view.transform[[view]] else "none"
        transform <- tolower(as.character(transform))
        if (transform %in% c("raw", "none")) transform <- "none"
        if (identical(transform, "log")) transform <- "log1p"
        applied[[view]] <- transform
        raw <- matrices[[view]]
        matrices[[view]] <- switch(transform,
            none = raw,
            clr = .mimesys_clr_transform(raw),
            log1p = .mimesys_log1p_transform(raw),
            stop("Unknown view transform '", transform, "' for view ", view,
                ".", call. = FALSE)
        )
        info[[view]] <- .mimesys_transform_info(raw, matrices[[view]],
            transform)
        attr(matrices[[view]], "pseudocount") <- NULL
        attr(matrices[[view]], "shift") <- NULL
    }
    attr(matrices, "view.transform") <- applied
    attr(matrices, "transform.info") <- info
    matrices
}

.mimesys_clr_transform <- function(x, pseudocount = NULL) {
    x <- as.matrix(x)
    storage.mode(x) <- "double"
    if (any(x < 0, na.rm = TRUE)) {
        stop("CLR transform requires non-negative values.", call. = FALSE)
    }
    positive <- x[is.finite(x) & x > 0]
    if (is.null(pseudocount)) {
        pseudocount <- if (length(positive) > 0L) min(positive,
            na.rm = TRUE) / 2 else 1e-6
    }
    if (!is.finite(pseudocount) || pseudocount <= 0) {
        pseudocount <- 1e-6
    }
    log_x <- log(x + pseudocount)
    out <- sweep(log_x, 1, rowMeans(log_x, na.rm = TRUE), "-")
    out[!is.finite(out)] <- NA_real_
    dimnames(out) <- dimnames(x)
    attr(out, "pseudocount") <- pseudocount
    out
}

.mimesys_log1p_transform <- function(x) {
    x <- as.matrix(x)
    storage.mode(x) <- "double"
    min_x <- suppressWarnings(min(x, na.rm = TRUE))
    shift <- 0
    if (is.finite(min_x) && min_x < -1) {
        shift <- abs(min_x) + 1
        x <- x + shift
    }
    out <- log1p(x)
    out[!is.finite(out)] <- NA_real_
    dimnames(out) <- dimnames(x)
    attr(out, "shift") <- shift
    out
}

## Summary information needed to map simulated values back to the original
## measurement scale: per-view constants only.
.mimesys_transform_info <- function(raw, transformed, transform) {
    raw <- as.matrix(raw)
    finite <- raw[is.finite(raw)]
    integer_valued <- length(finite) > 0L && all(abs(finite - round(finite)) < 1e-8)
    list(
        type = transform,
        shift = if (identical(transform, "log1p")) attr(transformed,
            "shift") else 0,
        pseudocount = if (identical(transform, "clr")) attr(transformed,
            "pseudocount") else NULL,
        integer = integer_valued,
        min = if (length(finite)) min(finite) else NA_real_
    )
}

## Exact inverse of the log1p transform (used for two-part views, where the
## positive values are generated on the fitted scale).
.mimesys_inverse_log1p <- function(x, info) {
    out <- expm1(x) - info$shift
    if (isTRUE(info$integer)) out <- round(out)
    if (is.finite(info$min)) out <- pmax(out, info$min)
    dimnames(out) <- dimnames(x)
    out
}



## Expand a per-view argument to one value per view, in template view order.
## Accepts a single value (used for every view), an unnamed vector with one
## value per view (template order), or a vector/list named by view. Named
## values are matched by name; unknown names and wrong lengths are errors.
.per_view_values <- function(x, views, arg_name) {
    if (is.list(x) && !is.data.frame(x)) x <- unlist(x, use.names = TRUE)
    M <- length(views)
    nm <- names(x)
    if (is.null(nm) || all(!nzchar(nm))) {
        if (length(x) == 1L) return(stats::setNames(rep(x, M), views))
        if (length(x) == M) return(stats::setNames(x, views))
        stop(arg_name, " must have length 1 or one value per view (", M,
            "), or be named by view: ", paste(views, collapse = ", "), ".",
            call. = FALSE)
    }
    if (any(!nzchar(nm))) {
        stop(arg_name, " must be either fully named by view or unnamed.", call. = FALSE)
    }
    unknown <- setdiff(nm, views)
    if (length(unknown) > 0L) {
        stop(arg_name, " has unknown view name(s): ", paste(unknown, collapse = ", "),
            ". Views are: ", paste(views, collapse = ", "), ".", call. = FALSE)
    }
    missing <- setdiff(views, nm)
    if (length(missing) > 0L) {
        stop(arg_name, " is missing view(s): ", paste(missing, collapse = ", "), ".",
            call. = FALSE)
    }
    x[views]
}
