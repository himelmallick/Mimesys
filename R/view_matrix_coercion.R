#' Coerce multi-omics data to MimeSys view matrices
#'
#' Converts common multi-omics containers to the internal MimeSys representation:
#' a named list of numeric matrices with samples in rows and features in columns.
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
#' @export
as_mimesys_view_list <- function(x, assay = NULL, single_view_name = "view1") {
  if (inherits(x, "MultiAssayExperiment")) {
    if (!requireNamespace("MultiAssayExperiment", quietly = TRUE)) {
      stop("MultiAssayExperiment support requires the MultiAssayExperiment package.", call. = FALSE)
    }
    x <- as.list(MultiAssayExperiment::experiments(x))
  } else if (inherits(x, "SummarizedExperiment")) {
    x <- stats::setNames(list(x), single_view_name)
  }

  if (!is.list(x) || length(x) == 0L) {
    stop("x must be a MultiAssayExperiment, SummarizedExperiment, or named list.", call. = FALSE)
  }
  if (is.null(names(x)) || any(!nzchar(names(x)))) {
    names(x) <- paste0("view", seq_along(x))
  }

  out <- lapply(seq_along(x), function(i) {
    view <- names(x)[[i]]
    .mimesys_one_view_matrix(x[[i]], assay = assay, view = view)
  })
  names(out) <- make.unique(names(x), sep = "_")
  out
}

.mimesys_one_view_matrix <- function(x, assay = NULL, view = "view") {
  if (inherits(x, "SummarizedExperiment")) {
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
      stop("SummarizedExperiment support requires the SummarizedExperiment package.", call. = FALSE)
    }
    assay_names <- SummarizedExperiment::assayNames(x)
    assay_id <- assay
    if (is.null(assay_id)) {
      assay_id <- 1L
    } else if (is.character(assay_id)) {
      hit <- match(assay_id, assay_names)
      if (is.na(hit)) {
        stop("Assay '", assay_id, "' not found for view ", view, ".", call. = FALSE)
      }
      assay_id <- hit
    }
    mat <- t(as.matrix(SummarizedExperiment::assay(x, assay_id)))
  } else {
    mat <- as.matrix(x)
  }

  storage.mode(mat) <- "double"
  if (!is.numeric(mat)) {
    stop("View ", view, " must be coercible to a numeric matrix.", call. = FALSE)
  }
  if (nrow(mat) < 2L || ncol(mat) < 1L) {
    stop("View ", view, " must contain at least two samples and one feature.", call. = FALSE)
  }
  if (is.null(rownames(mat))) {
    rownames(mat) <- paste0("sample_", seq_len(nrow(mat)))
  }
  if (is.null(colnames(mat))) {
    colnames(mat) <- paste0(view, "_feature_", seq_len(ncol(mat)))
  }
  colnames(mat) <- make.unique(as.character(colnames(mat)), sep = "_")
  mat
}

.mimesys_apply_view_transforms <- function(matrices, view.transform = NULL) {
  if (is.null(view.transform)) {
    attr(matrices, "view.transform") <- stats::setNames(rep("none", length(matrices)), names(matrices))
    return(matrices)
  }
  if (is.null(names(view.transform)) || any(!nzchar(names(view.transform)))) {
    if (length(view.transform) != 1L) {
      stop("view.transform must be named by view, or be a single shared transform.", call. = FALSE)
    }
    view.transform <- stats::setNames(rep(view.transform[[1]], length(matrices)), names(matrices))
  }

  applied <- stats::setNames(rep("none", length(matrices)), names(matrices))
  for (view in names(matrices)) {
    transform <- if (view %in% names(view.transform)) view.transform[[view]] else "none"
    transform <- tolower(as.character(transform))
    applied[[view]] <- transform
    matrices[[view]] <- switch(transform,
      none = matrices[[view]],
      raw = matrices[[view]],
      clr = .mimesys_clr_transform(matrices[[view]]),
      log = .mimesys_log1p_transform(matrices[[view]]),
      log1p = .mimesys_log1p_transform(matrices[[view]]),
      stop("Unknown view transform '", transform, "' for view ", view, ".", call. = FALSE)
    )
  }
  attr(matrices, "view.transform") <- applied
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
    pseudocount <- if (length(positive) > 0L) min(positive, na.rm = TRUE) / 2 else 1e-6
  }
  if (!is.finite(pseudocount) || pseudocount <= 0) {
    pseudocount <- 1e-6
  }
  log_x <- log(x + pseudocount)
  out <- sweep(log_x, 1, rowMeans(log_x, na.rm = TRUE), "-")
  out[!is.finite(out)] <- NA_real_
  dimnames(out) <- dimnames(x)
  out
}

.mimesys_log1p_transform <- function(x) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  min_x <- suppressWarnings(min(x, na.rm = TRUE))
  if (is.finite(min_x) && min_x < -1) {
    x <- x + abs(min_x) + 1
  }
  out <- log1p(x)
  out[!is.finite(out)] <- NA_real_
  dimnames(out) <- dimnames(x)
  out
}
