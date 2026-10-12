\name{run_tcga_mimesys_validation}
\alias{run_tcga_mimesys_validation}
\title{Run TCGA multi-omics MimeSys validation}
\usage{
run_tcga_mimesys_validation(
  cancer,
  base_dir = path.expand("~/Mimesys_validation"),
  crosstalk_grid = c(0, 0.25, 0.5, 0.75, 1),
  feature_sets = list(small = c(mrna = 500, mirna = 100, protein = 50), medium =
    c(mrna = 1000, mirna = 250, protein = 100), large = c(mrna = 2000, mirna =
    500, protein = 200)),
  sample_sizes = NULL,
  nrep = 25,
  snr = 2,
  seed = 1000,
  files.per.chunk = 20,
  plot_width = 10,
  plot_height = 7
)
}
\arguments{
\item{cancer}{TCGA cancer code, such as \code{"HNSC"}, \code{"BRCA"},
\code{"LUAD"}, \code{"KIRC"}, or \code{"UCEC"}.}

\item{base_dir}{Base output directory. A cancer-specific subdirectory is
created inside this directory.}

\item{crosstalk_grid}{Numeric crosstalk values to evaluate.}

\item{feature_sets}{Named list of feature counts per view. Each element must be
a named integer vector with names \code{mrna}, \code{mirna}, and
\code{protein}.}

\item{sample_sizes}{Integer vector of simulated subject counts to evaluate. If
\code{NULL}, uses the number of matched template samples.}

\item{nrep}{Number of simulation replicates per grid setting.}

\item{snr}{Signal-to-noise ratio passed to \code{gen_simmba()}.}

\item{seed}{Base seed used for the validation grid.}

\item{files.per.chunk}{Number of files per GDC download chunk.}

\item{plot_width, plot_height}{PDF dimensions in inches.}
}
\value{
A list containing output directories, diagnostics, summaries, and plot paths.
}
\description{
Downloads or reuses TCGA mRNA and miRNA data from GDC, downloads or reuses RPPA
protein data from UCSC Xena, builds a matched three-view template, runs the
MimeSys crosstalk validation grid, and saves diagnostics and plots in a
cancer-specific output directory.
}
