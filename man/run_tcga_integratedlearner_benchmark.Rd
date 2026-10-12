\name{run_tcga_integratedlearner_benchmark}
\alias{run_tcga_integratedlearner_benchmark}
\title{Run IntegratedLearner benchmark using cached TCGA MimeSys templates}
\usage{
run_tcga_integratedlearner_benchmark(
  cancers = c("HNSC", "LUAD"),
  base_dir = path.expand("~/Mimesys_validation"),
  outcome_types = c("continuous", "binary"),
  feature_setting = "medium",
  feature_sets = list(small = c(mrna = 500, mirna = 100, protein = 50), medium =
    c(mrna = 1000, mirna = 250, protein = 100), large = c(mrna = 2000, mirna =
    500, protein = 200)),
  crosstalk_grid = c(0, 0.25, 0.5, 0.75, 1),
  nrep = 5,
  snr = 2,
  folds = 2,
  base_learner = "SL.randomForest",
  meta_learner = "SL.nnls.auc",
  seed = 5000
)
}
\arguments{
\item{cancers}{Character vector of TCGA cancer codes.}

\item{base_dir}{Base directory containing cancer-specific MimeSys validation
folders.}

\item{outcome_types}{Outcomes to benchmark: \code{"continuous"} and/or
\code{"binary"}.}

\item{feature_setting}{Which feature-size setting to use.}

\item{feature_sets}{Named list of feature counts per view.}

\item{crosstalk_grid}{Numeric crosstalk values.}

\item{nrep}{Number of simulated replicates per setting.}

\item{snr}{Signal-to-noise ratio for \code{gen_simmba()}.}

\item{folds}{Number of folds used inside \code{IntegratedLearner()}.}

\item{base_learner}{Base learner passed to \code{IntegratedLearner()}.}

\item{meta_learner}{Meta learner passed to \code{IntegratedLearner()} for
stacked late fusion.}

\item{seed}{Base seed.}
}
\value{
A list containing raw results, summary table, and output paths.
}
\description{
Loads TCGA templates created by \code{run_tcga_mimesys_validation()}, simulates
continuous and/or binary outcomes under a crosstalk/effect-alignment grid, runs
IntegratedLearner, and saves raw and summarized benchmark results.
}
