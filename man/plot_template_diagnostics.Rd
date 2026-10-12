\name{plot_template_diagnostics}
\alias{plot_template_diagnostics}
\title{Plot template simulation diagnostics}
\usage{
plot_template_diagnostics(
  diagnostics,
  metrics = c("mean_abs_feature_cor", "signal_agreement", "y_signal_cor")
)
}
\arguments{
\item{diagnostics}{A data frame returned by \code{validate_template_simulation()}
or \code{summarize_template_diagnostics()}.}

\item{metrics}{Character vector of diagnostic columns to plot.}
}
\value{
A \code{ggplot} object.
}
\description{
Visualizes whether cross-view feature correlation, aligned signal, and outcome
signal behave as expected across crosstalk settings.
}
