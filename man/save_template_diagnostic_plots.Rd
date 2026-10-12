\name{save_template_diagnostic_plots}
\alias{save_template_diagnostic_plots}
\title{Save publication-style template diagnostic plots}
\usage{
save_template_diagnostic_plots(
  diagnostics,
  plot_dir,
  prefix = "template_diagnostics",
  width = 10,
  height = 7,
  combined_pdf = TRUE
)
}
\arguments{
\item{diagnostics}{A data frame returned by
\code{validate_template_simulation()}. This should usually contain multiple
crosstalk settings and simulation replicates.}

\item{plot_dir}{Directory where PDF files should be saved.}

\item{prefix}{File prefix used for the saved PDFs.}

\item{width, height}{PDF dimensions in inches.}

\item{combined_pdf}{Logical; if \code{TRUE}, also save one multi-page PDF
containing all diagnostic plots.}
}
\value{
A named character vector containing the saved PDF paths.
}
\description{
Creates a standard plot set for validating template-based simulations: a
crosstalk diagnostic line plot, a feature-correlation heatmap, a replicate
variability plot, and a feature-correlation versus signal-agreement scatter
plot.
}
