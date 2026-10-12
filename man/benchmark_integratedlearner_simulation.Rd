\name{benchmark_integratedlearner_simulation}
\alias{benchmark_integratedlearner_simulation}
\title{Benchmark IntegratedLearner on a MimeSys simulation}
\usage{
benchmark_integratedlearner_simulation(
  sim,
  outcome_type = c("continuous", "binary"),
  folds = 2,
  base_learner = "SL.randomForest",
  meta_learner = "SL.nnls.auc",
  run_concat = TRUE,
  run_stacked = TRUE,
  extra_metadata = NULL
)
}
\arguments{
\item{sim}{A list returned by template-based \code{gen_simmba()}.}

\item{outcome_type}{Either \code{"continuous"} or \code{"binary"}.}

\item{folds}{Number of cross-validation folds used inside
\code{IntegratedLearner()}.}

\item{base_learner}{Base learner passed to \code{IntegratedLearner()}.}

\item{meta_learner}{Meta learner passed to \code{IntegratedLearner()} for
stacked late fusion.}

\item{run_concat, run_stacked}{Logical flags passed to \code{IntegratedLearner()}.}

\item{extra_metadata}{Optional one-row data frame of simulation setting metadata
to attach to every output row.}
}
\value{
A data frame with one row per replicate and learner.
}
\description{
Converts each train/test replicate from \code{gen_simmba()} into
\code{MultiAssayExperiment} objects, runs \code{IntegratedLearner()}, and returns
predictive performance for single-view and fusion learners.
}
