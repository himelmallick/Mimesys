#' Run TCGA multi-omics MimeSys validation
#'
#' Downloads or reuses TCGA mRNA and miRNA data from GDC, downloads or reuses
#' RPPA protein data from UCSC Xena, builds a matched three-view template, runs
#' the MimeSys crosstalk validation grid, and saves diagnostics and plots in a
#' cancer-specific output directory.
#'
#' @param cancer TCGA cancer code, such as `"HNSC"`, `"BRCA"`, `"LUAD"`,
#'   `"KIRC"`, or `"UCEC"`.
#' @param base_dir Base output directory. A cancer-specific subdirectory is
#'   created inside this directory.
#' @param crosstalk_grid Numeric crosstalk values to evaluate.
#' @param feature_sets Named list of feature counts per view. Each element must
#'   be a named integer vector with names `mrna`, `mirna`, and `protein`.
#' @param sample_sizes Integer vector of simulated subject counts to evaluate.
#'   If `NULL`, uses the number of matched template samples.
#' @param nrep Number of simulation replicates per grid setting.
#' @param snr Signal-to-noise ratio passed to `gen_simmba()`.
#' @param seed Base seed used for the validation grid.
#' @param files.per.chunk Number of files per GDC download chunk.
#' @param plot_width,plot_height PDF dimensions in inches.
#'
#' @return A list containing output directories, diagnostics, summaries, and
#'   plot paths.
#' @export
run_tcga_mimesys_validation <- function(
  cancer,
  base_dir = path.expand("~/Mimesys_validation"),
  crosstalk_grid = c(0, 0.25, 0.5, 0.75, 1),
  feature_sets = list(
    small = c(mrna = 500, mirna = 100, protein = 50),
    medium = c(mrna = 1000, mirna = 250, protein = 100),
    large = c(mrna = 2000, mirna = 500, protein = 200)
  ),
  sample_sizes = NULL,
  nrep = 25,
  snr = 2,
  seed = 1000,
  files.per.chunk = 20,
  plot_width = 10,
  plot_height = 7
) {
  for (pkg in c("TCGAbiolinks", "SummarizedExperiment", "UCSCXenaTools", "stringr")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("run_tcga_mimesys_validation() requires package '%s'.", pkg), call. = FALSE)
    }
  }

  cancer <- toupper(cancer)
  project <- cancer
  gdc_project <- paste0("TCGA-", cancer)

  cancer_dir <- file.path(path.expand(base_dir), cancer)
  gdc_dir <- file.path(cancer_dir, paste0("GDCdata_", cancer))
  xena_dir <- file.path(cancer_dir, paste0("xena_tcga_", cancer))
  cache_dir <- file.path(cancer_dir, paste0("prepared_", cancer))
  plot_dir <- file.path(cancer_dir, "plots")

  for (dir_now in c(cancer_dir, gdc_dir, xena_dir, cache_dir, plot_dir)) {
    dir.create(dir_now, recursive = TRUE, showWarnings = FALSE)
  }

  tcga_sample_id <- function(x) {
    stringr::str_extract(
      as.character(x),
      "^TCGA-[A-Za-z0-9]{2}-[A-Za-z0-9]{4}-[0-9]{2}"
    )
  }

  collapse_duplicate_samples <- function(x) {
    if (!anyDuplicated(rownames(x))) {
      return(x)
    }
    ids <- rownames(x)
    summed <- rowsum(x, group = ids, reorder = FALSE, na.rm = TRUE)
    counts <- table(ids)
    summed / as.numeric(counts[rownames(summed)])
  }

  remove_bad_features <- function(x) {
    enough_values <- apply(x, 2, function(z) sum(is.finite(z)) >= 2)
    x <- x[, enough_values, drop = FALSE]
    vars <- apply(x, 2, stats::var, na.rm = TRUE)
    keep <- is.finite(vars) & vars > 0
    x[, keep, drop = FALSE]
  }

  keep_top_variable <- function(x, n_features) {
    vars <- apply(x, 2, stats::var, na.rm = TRUE)
    vars <- vars[is.finite(vars)]
    vars <- sort(vars, decreasing = TRUE)
    n_keep <- min(n_features, length(vars))
    if (n_keep == 0) {
      stop("No variable features available.", call. = FALSE)
    }
    x[, names(vars)[seq_len(n_keep)], drop = FALSE]
  }

  check_tcga_ids <- function(x, name) {
    ok <- grepl("^TCGA-[A-Za-z0-9]{2}-[A-Za-z0-9]{4}-[0-9]{2}$", x)
    message(sprintf("%s: %s/%s valid TCGA sample IDs", name, sum(ok), length(x)))
    invisible(ok)
  }

  extract_xena_rppa <- function(xena_data) {
    if (is.data.frame(xena_data) || is.matrix(xena_data)) {
      return(as.data.frame(xena_data, check.names = FALSE))
    }

    if (is.list(xena_data)) {
      nm <- names(xena_data)
      if (!is.null(nm)) {
        rppa_idx <- grep("RPPA|protein", nm, ignore.case = TRUE)
        if (length(rppa_idx) > 0) {
          return(as.data.frame(xena_data[[rppa_idx[1]]], check.names = FALSE))
        }
      }

      tabular <- vapply(
        xena_data,
        function(z) is.data.frame(z) || is.matrix(z),
        logical(1)
      )
      if (sum(tabular) == 1L) {
        return(as.data.frame(xena_data[[which(tabular)[1]]], check.names = FALSE))
      }
      if (length(xena_data) == 1L) {
        return(as.data.frame(xena_data[[1]], check.names = FALSE))
      }
    }

    stop(
      "Could not identify RPPA/protein object from Xena. ",
      "Inspect the object returned by UCSCXenaTools::XenaPrepare().",
      call. = FALSE
    )
  }

  message("Preparing ", cancer, " mRNA data...")
  rna_cache <- file.path(cache_dir, paste0(cancer, "_rna_se.rds"))
  if (file.exists(rna_cache)) {
    rna_se <- readRDS(rna_cache)
  } else {
    query_rna <- TCGAbiolinks::GDCquery(
      project = gdc_project,
      data.category = "Transcriptome Profiling",
      data.type = "Gene Expression Quantification",
      workflow.type = "STAR - Counts",
      sample.type = "Primary Tumor"
    )
    rna_se <- tryCatch(
      TCGAbiolinks::GDCprepare(query_rna, directory = gdc_dir),
      error = function(e) {
        TCGAbiolinks::GDCdownload(
          query_rna,
          directory = gdc_dir,
          method = "api",
          files.per.chunk = files.per.chunk
        )
        TCGAbiolinks::GDCprepare(query_rna, directory = gdc_dir)
      }
    )
    saveRDS(rna_se, rna_cache)
  }

  assay_names <- SummarizedExperiment::assayNames(rna_se)
  rna_counts <- if ("unstranded" %in% assay_names) {
    SummarizedExperiment::assay(rna_se, "unstranded")
  } else {
    SummarizedExperiment::assay(rna_se, 1)
  }
  rna_row_data <- SummarizedExperiment::rowData(rna_se)
  gene_names <- if ("gene_name" %in% colnames(rna_row_data)) {
    as.character(rna_row_data$gene_name)
  } else {
    rownames(rna_counts)
  }
  rownames(rna_counts) <- make.unique(gene_names)
  rna_mat <- t(log2(rna_counts + 1))
  rownames(rna_mat) <- tcga_sample_id(rownames(rna_mat))
  if (any(is.na(rownames(rna_mat)))) {
    stop("Some RNA sample IDs could not be parsed.", call. = FALSE)
  }
  rna_mat <- collapse_duplicate_samples(rna_mat)
  check_tcga_ids(rownames(rna_mat), "RNA")

  message("Preparing ", cancer, " miRNA data...")
  mirna_cache <- file.path(cache_dir, paste0(cancer, "_mirna_obj.rds"))
  if (file.exists(mirna_cache)) {
    mirna_obj <- readRDS(mirna_cache)
  } else {
    query_mirna <- TCGAbiolinks::GDCquery(
      project = gdc_project,
      data.category = "Transcriptome Profiling",
      data.type = "miRNA Expression Quantification",
      workflow.type = "BCGSC miRNA Profiling",
      sample.type = "Primary Tumor"
    )
    mirna_obj <- tryCatch(
      TCGAbiolinks::GDCprepare(query_mirna, directory = gdc_dir),
      error = function(e) {
        TCGAbiolinks::GDCdownload(
          query_mirna,
          directory = gdc_dir,
          method = "api",
          files.per.chunk = files.per.chunk
        )
        TCGAbiolinks::GDCprepare(query_mirna, directory = gdc_dir)
      }
    )
    saveRDS(mirna_obj, mirna_cache)
  }

  make_mirna_matrix <- function(obj, value_type = "rpm") {
    if (inherits(obj, "SummarizedExperiment")) {
      assay_names <- SummarizedExperiment::assayNames(obj)
      rpm_idx <- grep("rpm|million", assay_names, ignore.case = TRUE)
      mat <- if (value_type == "rpm" && length(rpm_idx) > 0) {
        SummarizedExperiment::assay(obj, rpm_idx[1])
      } else {
        SummarizedExperiment::assay(obj, 1)
      }
      mat <- as.matrix(mat)
      mode(mat) <- "numeric"
      mat <- t(mat)
      rownames(mat) <- tcga_sample_id(rownames(mat))
      if (any(is.na(rownames(mat)))) {
        stop("Some miRNA sample IDs could not be parsed.", call. = FALSE)
      }
      return(log2(collapse_duplicate_samples(mat) + 1))
    }

    df <- as.data.frame(obj, check.names = FALSE)
    if (!"miRNA_ID" %in% names(df)) {
      stop("Could not find miRNA_ID column.", call. = FALSE)
    }
    prefix <- if (value_type == "rpm") {
      "reads_per_million_miRNA_mapped_"
    } else {
      "read_count_"
    }
    value_cols <- grep(paste0("^", prefix), names(df), value = TRUE)
    if (length(value_cols) == 0) {
      stop("No miRNA columns found with prefix: ", prefix, call. = FALSE)
    }
    mat <- as.matrix(df[, value_cols, drop = FALSE])
    mode(mat) <- "numeric"
    rownames(mat) <- make.unique(as.character(df$miRNA_ID))
    sample_ids <- tcga_sample_id(sub(paste0("^", prefix), "", value_cols))
    if (any(is.na(sample_ids))) {
      stop("Some miRNA sample IDs could not be parsed.", call. = FALSE)
    }
    colnames(mat) <- sample_ids
    mat <- t(mat)
    log2(collapse_duplicate_samples(mat) + 1)
  }

  mirna_mat <- make_mirna_matrix(mirna_obj, value_type = "rpm")
  check_tcga_ids(rownames(mirna_mat), "miRNA")

  message("Preparing ", cancer, " RPPA data...")
  rppa_cache <- file.path(cache_dir, paste0(cancer, "_rppa.rds"))
  if (file.exists(rppa_cache)) {
    rppa <- readRDS(rppa_cache)
  } else {
    existing_rppa <- list.files(
      xena_dir,
      pattern = "RPPA.*\\.gz$",
      recursive = TRUE,
      full.names = TRUE
    )
    if (length(existing_rppa) > 0) {
      rppa <- utils::read.delim(
        gzfile(existing_rppa[1]),
        header = TRUE,
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
    } else {
      xena_query <- UCSCXenaTools::getTCGAdata(
        project = project,
        clinical = FALSE,
        mRNASeq = FALSE,
        miRNASeq = FALSE,
        RPPAArray = TRUE,
        download = FALSE
      )
      xq <- UCSCXenaTools::XenaQuery(xena_query$Xena)
      xd <- UCSCXenaTools::XenaDownload(xq, destdir = xena_dir)
      xena_data <- UCSCXenaTools::XenaPrepare(xd)
      rppa <- extract_xena_rppa(xena_data)
    }
    saveRDS(rppa, rppa_cache)
  }

  rppa <- as.data.frame(rppa, check.names = FALSE)
  protein_name_col <- if ("sample" %in% names(rppa)) {
    "sample"
  } else {
    names(rppa)[1]
  }
  protein_names <- make.unique(as.character(rppa[[protein_name_col]]))
  rppa_numeric <- rppa[, names(rppa) != protein_name_col, drop = FALSE]
  protein_mat <- as.matrix(rppa_numeric)
  mode(protein_mat) <- "numeric"
  rownames(protein_mat) <- protein_names
  protein_mat <- t(protein_mat)
  rownames(protein_mat) <- tcga_sample_id(rownames(protein_mat))
  if (any(is.na(rownames(protein_mat)))) {
    stop("Some RPPA sample IDs could not be parsed.", call. = FALSE)
  }
  protein_mat <- collapse_duplicate_samples(protein_mat)
  check_tcga_ids(rownames(protein_mat), "RPPA")

  common_samples <- Reduce(
    intersect,
    list(rownames(rna_mat), rownames(mirna_mat), rownames(protein_mat))
  )
  message("Three-way matched samples: ", length(common_samples))
  if (length(common_samples) == 0) {
    stop("No common samples across mRNA, miRNA, and RPPA.", call. = FALSE)
  }
  if (is.null(sample_sizes)) {
    sample_sizes <- length(common_samples)
  }
  sample_sizes <- as.integer(sample_sizes)
  if (any(is.na(sample_sizes)) || any(sample_sizes < 10)) {
    stop("sample_sizes must contain integers >= 10.", call. = FALSE)
  }

  template <- list(
    mrna = rna_mat[common_samples, , drop = FALSE],
    mirna = mirna_mat[common_samples, , drop = FALSE],
    protein = protein_mat[common_samples, , drop = FALSE]
  )
  stopifnot(identical(rownames(template$mrna), rownames(template$mirna)))
  stopifnot(identical(rownames(template$mrna), rownames(template$protein)))

  template_clean <- lapply(template, remove_bad_features)
  template_default <- list(
    mrna = keep_top_variable(template_clean$mrna, feature_sets$large[["mrna"]]),
    mirna = keep_top_variable(template_clean$mirna, feature_sets$large[["mirna"]]),
    protein = keep_top_variable(template_clean$protein, feature_sets$large[["protein"]])
  )

  saveRDS(template, file.path(cache_dir, paste0(cancer, "_multiomics_template_full.rds")))
  saveRDS(template_clean, file.path(cache_dir, paste0(cancer, "_multiomics_template_clean.rds")))
  saveRDS(template_default, file.path(cache_dir, paste0(cancer, "_multiomics_template_large.rds")))

  template_feature_sets <- lapply(feature_sets, function(counts) {
    list(
      mrna = keep_top_variable(template_clean$mrna, counts[["mrna"]]),
      mirna = keep_top_variable(template_clean$mirna, counts[["mirna"]]),
      protein = keep_top_variable(template_clean$protein, counts[["protein"]])
    )
  })
  for (nm in names(template_feature_sets)) {
    stopifnot(identical(rownames(template_feature_sets[[nm]]$mrna), rownames(template_feature_sets[[nm]]$mirna)))
    stopifnot(identical(rownames(template_feature_sets[[nm]]$mrna), rownames(template_feature_sets[[nm]]$protein)))
  }

  de_align_grid <- list(random_effects = NULL, aligned_effects = 1)
  validation_results <- list()
  counter <- 1L

  for (feature_setting in names(template_feature_sets)) {
    for (n_subjects in sample_sizes) {
      for (de_name in names(de_align_grid)) {
        for (ct in crosstalk_grid) {
          message(
            sprintf(
              "Running %s | %s | n=%s | %s | crosstalk=%s",
              cancer, feature_setting, n_subjects, de_name, ct
            )
          )
          template_input <- template_feature_sets[[feature_setting]]
          sim <- gen_simmba(
            nsample = n_subjects,
            template = fit_template(template_input),
            crosstalk = ct,
            de.align = de_align_grid[[de_name]],
            nrep = nrep,
            snr = snr,
            seed = seed + counter
          )
          diag_now <- validate_template_simulation(sim)
          diag_now$cancer <- cancer
          diag_now$feature_setting <- feature_setting
          diag_now$n_subjects <- n_subjects
          diag_now$de_align_label <- de_name
          diag_now$setting_id <- counter
          validation_results[[counter]] <- diag_now
          counter <- counter + 1L
        }
      }
    }
  }

  diag_all <- do.call(rbind, validation_results)
  raw_rds <- file.path(cancer_dir, paste0(cancer, "_validation_diagnostics_raw.rds"))
  raw_csv <- file.path(cancer_dir, paste0(cancer, "_validation_diagnostics_raw.csv"))
  saveRDS(diag_all, raw_rds)
  utils::write.csv(diag_all, raw_csv, row.names = FALSE)

  diag_summary <- stats::aggregate(
    cbind(mean_abs_feature_cor, signal_agreement, y_signal_cor) ~
      cancer + feature_setting + n_subjects + de_align_label + crosstalk + view_a + view_b,
    data = diag_all,
    FUN = function(x) mean(x, na.rm = TRUE)
  )
  summary_rds <- file.path(cancer_dir, paste0(cancer, "_validation_diagnostics_summary.rds"))
  summary_csv <- file.path(cancer_dir, paste0(cancer, "_validation_diagnostics_summary.csv"))
  saveRDS(diag_summary, summary_rds)
  utils::write.csv(diag_summary, summary_csv, row.names = FALSE)

  robustness_table <- stats::aggregate(
    cbind(mean_abs_feature_cor, signal_agreement, y_signal_cor) ~
      cancer + feature_setting + n_subjects + de_align_label + crosstalk,
    data = diag_all,
    FUN = function(x) mean(x, na.rm = TRUE)
  )
  robustness_rds <- file.path(cancer_dir, paste0(cancer, "_robustness_table.rds"))
  robustness_csv <- file.path(cancer_dir, paste0(cancer, "_robustness_table.csv"))
  saveRDS(robustness_table, robustness_rds)
  utils::write.csv(robustness_table, robustness_csv, row.names = FALSE)

  plot_grid <- expand.grid(
    feature_setting = names(template_feature_sets),
    n_subjects = sample_sizes,
    stringsAsFactors = FALSE
  )
  plot_paths <- lapply(seq_len(nrow(plot_grid)), function(i) {
    feature_setting <- plot_grid$feature_setting[[i]]
    n_subjects <- plot_grid$n_subjects[[i]]
    diag_feature <- diag_all[
      diag_all$feature_setting == feature_setting &
        diag_all$n_subjects == n_subjects,
      ,
      drop = FALSE
    ]
    save_template_diagnostic_plots(
      diagnostics = diag_feature,
      plot_dir = plot_dir,
      prefix = paste0(cancer, "_", feature_setting, "_n", n_subjects, "_validation"),
      width = plot_width,
      height = plot_height,
      combined_pdf = TRUE
    )
  })
  names(plot_paths) <- paste(plot_grid$feature_setting, paste0("n", plot_grid$n_subjects), sep = "_")

  session_file <- file.path(cancer_dir, paste0(cancer, "_validation_sessionInfo.txt"))
  capture.output(utils::sessionInfo(), file = session_file)

  out <- list(
    cancer = cancer,
    directories = list(
      cancer_dir = cancer_dir,
      gdc_dir = gdc_dir,
      xena_dir = xena_dir,
      cache_dir = cache_dir,
      plot_dir = plot_dir
    ),
    sample_counts = c(
      mrna = nrow(rna_mat),
      mirna = nrow(mirna_mat),
      protein = nrow(protein_mat),
      matched = length(common_samples)
    ),
    sample_sizes = sample_sizes,
    template_dimensions = lapply(template_feature_sets, function(x) vapply(x, dim, integer(2))),
    diagnostics = diag_all,
    summary = diag_summary,
    robustness_table = robustness_table,
    files = list(
      raw_rds = raw_rds,
      raw_csv = raw_csv,
      summary_rds = summary_rds,
      summary_csv = summary_csv,
      robustness_rds = robustness_rds,
      robustness_csv = robustness_csv,
      session_info = session_file
    ),
    plot_paths = plot_paths
  )

  message("Validation complete for ", cancer, ". Results saved to: ", cancer_dir)
  out
}
