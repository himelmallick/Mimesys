# MimeSys

MimeSys is a template-based simulator for multi-omics data. It learns each
feature's distribution, the dependence within each omics view and the
dependence between views from a real data set (the *template*), then generates
synthetic multi-view data with any sample size and feature count, a single dial
for cross-view dependence (*crosstalk*), and planted signal with known truth:
outcomes (continuous, binary or survival), several heterogeneous studies, or
exposure-mediator-outcome structures (parallel or sequential mediation). It is
domain-agnostic and works with `MultiAssayExperiment` objects.

## Installation

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("himelmallick/Mimesys")
library(MimeSys)
```

## Quick start

`simulate_mimesys()` runs every simulation; `design` chooses what is
simulated.

```r
library(MimeSys)
library(MultiAssayExperiment)
data(stelzer_mae)

# Fit the template once and reuse it
template <- fit_template(
  stelzer_mae,
  view.transform = c(Metabolomics = "log1p", Proteomics = "log1p")
)

# Outcome study (default design); crosstalk is calibrated to the template
sim_mae <- simulate_mimesys(template, n = 200, p = 30,
  outcome = list(model = "linear", type = "survival", cens_rate = 0.3))
colData(sim_mae)[, c("time", "status", "split")]
rowData(experiments(sim_mae)[["Proteomics"]])       # feature truth (beta, ...)

# Null data, no outcome
null_mae <- simulate_mimesys(template, design = "null", n = 200, p = 200)

# Several studies with between-study heterogeneity; study 3 held out
ms_mae <- simulate_mimesys(template, design = "multistudy",
  n = c(150, 200, 100), p = 30,
  heterogeneity = list(effect_cor = 0.5, shift_sd = 0.3, holdout_studies = 3))
table(ms_mae$study, ms_mae$split)

# Mediation with exact natural effects
med_mae <- simulate_mimesys(template, design = "mediation", n = 300,
  mediation = list(n_mediators = c(CYTOF = 3, Metabolomics = 3, Proteomics = 0),
                   tie = 1, de = 0.5))
metadata(med_mae)$truth$mediation[c("NIE", "NDE", "TE")]
```

Every function that uses random numbers has `seed = 1234` by default and
leaves your own random-number state unchanged.

## Parameters

"Per view" settings accept one value for every view, one value per view in
template order, or a vector named by view.

### `simulate_mimesys()`

| Argument | Default | Description |
|---|---|---|
| `template` | required | `MultiAssayExperiment`, `SummarizedExperiment`, named list of matrices (samples in rows) or a fitted template. For `"multistudy"` also one template or data set per study (see `fit_study_templates()`). |
| `n` | required | Samples per replicate. For `"multistudy"`: one value per study, or one value with `nstudy`. |
| `p` | `NULL` | Features per view (per view). `NULL`: template sizes. Fewer subsamples features; more adds "offspring" features that copy a parent's marginal. |
| `crosstalk` | `"calibrate"` | Cross-view dependence. `"calibrate"` picks the value that best matches the template's cross-view correlations; a number skips calibration (`0` independent views, `1` template level, up to `crosstalk_max()`). |
| `design` | `"outcome"` | `"outcome"`, `"null"`, `"multistudy"` or `"mediation"` (table below). |
| `preset` | `"default"` | Outcome preset (`mimesys_presets()`): `default`/`linear`, `friedman`, `friedman2`, `shared_latent`, `shared_latent_friedman`, `shared_latent_friedman2`, `interaction`. Explicit options override it. |
| `outcome` | `NULL` | Outcome options (list below), or just the model name. |
| `signal` | `list()` | Signal options (list below). |
| `heterogeneity` | `list()` | Between-study differences (list below). |
| `nstudy` | `length(n)` | Number of studies. |
| `mediation` | `list()` | Mediation options (list below). |
| `null_model` | `list()` | Feature-generation options (list below). |
| `template_args` | `list()` | Passed to `fit_template()` when `template` is raw data. |
| `calibration` | `list()` | Crosstalk calibration options (list below). |
| `train_fraction` | `0.7` | Share of samples labelled `"train"` (per non-holdout study for `"multistudy"`). |
| `nrep` | `1` | Replicates. One returns a `MultiAssayExperiment`; more return a list of them (`$replicates`). |
| `seed` | `1234` | Random seed; your own random-number state is restored afterwards. |
| `return_type` | `"mae"` | `"mae"`, or `"simulation"` for the underlying engine result. |

Each design reads only its own option lists; supplying one it does not use is
an error.

| `design` | Simulates | Engine | Option lists |
|---|---|---|---|
| `"outcome"` | features with planted outcome signal and a train/test split | `gen_simmba()` | `preset`, `outcome`, `signal`, `train_fraction` |
| `"null"` | features only | `sim_from_template()` | none |
| `"multistudy"` | several heterogeneous studies, pooled in one MAE with a `study` column | `gen_multistudy()` | `preset`, `outcome`, `signal`, `heterogeneity`, `nstudy`, `train_fraction` |
| `"mediation"` | exposure, mediators, outcome and exact natural effects | `gen_mediation()` | `mediation` |

**`outcome`**

| Option | Default | Description |
|---|---|---|
| `model` | `"linear"` | `"linear"`, `"friedman"` (nonlinear Friedman function of five signal features), `"friedman2"` (linear plus Friedman), `"shared_latent"` (driven by the shared cross-view factors) or `"interaction"` (linear plus cross-view feature products). |
| `type` | `"continuous"` | `"continuous"`, `"binary"` or `"survival"`. |
| `snr` | `1` | Signal-to-noise ratio (on the latent logistic or log-time scale for binary and survival outcomes). |
| `prevalence` | `0.5` | Target case proportion (binary). |
| `friedman_share` | `NULL` | Per view: which views supply the Friedman features (defaults to `signal$share`). |
| `nonlinear_weight` | `0.5` | Share of signal variance from the nonlinear part (`"friedman2"`). |
| `interaction_weight` | `0.5` | Share of signal variance from interactions (`"interaction"`; 1 = no main effects). |
| `n_interactions` | `NULL` | Number of cross-view feature pairs (default: one per view pair). |
| `shared_latent_function` | `"linear"` | For `"shared_latent"`: `"linear"`, `"friedman"` or `"friedman2"` (the Friedman forms need 5 shared factors). |
| `shared_latent_n_factors` | `NULL` | For `"shared_latent"`: shared factors used (`NULL`: all). |
| `shared_latent_weight` | `1` | For `"shared_latent"`: share of signal from the shared factors; the rest is view-specific, so lower values make the views complementary. |
| `baseline_hazard` | `1` | Baseline hazard $h_0$ (survival). |
| `baseline_shape` | `1` | Weibull shape $k$ (1 = exponential; $k > 1$ rising, $k < 1$ falling hazard). |
| `cens_rate` | `NULL` | Target censoring proportion; `NULL` uses uniform censoring on (`cens_lower`, `cens_upper`) = (1, 3). For `"multistudy"`, one value or one per study. |
| `cens_type` | `"uniform"` | `"uniform"`, `"exponential"` (constant loss to follow-up) or `"administrative"` (staggered entry, fixed study end); scaled to reach `cens_rate`. |
| `cens_accrual`, `cens_dropout` | `0.5`, `0` | For `"administrative"`: share of the study over which subjects enter; probability of dropping out before the end. |
| `cens_times` | `NULL` | Your own censoring times: one time, one per sample, or a function of `n` (rescaled to `cens_rate` if given). |
| `cens_dependence` | `0` | Informative censoring: with a positive value higher-risk subjects are censored earlier. |
| `keep_true_times` | `FALSE` | Also return the latent `event_time` and `cens_time`. |

**`signal`**

| Option | Default | Description |
|---|---|---|
| `prob` | `0.1` | Per view: probability that a feature carries signal (at least one signal feature per view with `prob > 0`). |
| `neg_prob` | `0.5` | Per view: probability that a signal effect is negative. |
| `effect_meanlog`, `effect_sdlog` | `1`, `0.4` | Per view: log-normal effect-size factor $f$; the effect size is $\lvert\log_2 f\rvert$. |
| `align` | `1` | Share of signal features placed on, and signed by, the shared cross-view direction (`NULL`: random). |
| `share` | `NULL` | Per view: relative variance of each view's signal (`c(1, 1, 1)` balanced, `c(1, 0, 0)` one view only). |

**`heterogeneity`** (`"multistudy"`)

| Option | Default | Description |
|---|---|---|
| `effect_cor` | `1` | Correlation of the effect sizes between studies. |
| `support_overlap` | `1` | Share of each study's signal features taken from the common set (the rest are study-specific). |
| `intercept_sd` | `0` | SD of study intercepts (continuous: signal-SD units; binary: logit; survival: log hazard). |
| `snr_sd` | `0` | SD of the study log signal-to-noise ratios. |
| `shift_sd` | `0` | Case mix: SD of study shifts of the latent factors (biological; enters the outcome). |
| `batch_loc_sd`, `batch_scale_sd` | `0`, `0` | Technical batch effects per study and feature (added after the outcome). |
| `holdout_studies` | `NULL` | Studies whose samples all go to the test split. |

**`mediation`** (`"mediation"`)

| Option | Default | Description |
|---|---|---|
| `n_mediators` | `5` | Mediators per view (per view; 0 for a view without mediators). |
| `tie`, `de` | `1`, `0` | Target total natural indirect effect (calibrated exactly) and direct effect. |
| `dir_a`, `dir_b` | `1`, `1` | Directions of exposure-to-mediator and mediator-to-outcome paths: `1`, `-1` or `"random"`. |
| `effect_type` | `"spiked"` | `"spiked"`, or the nulls `"null_a"`, `"null_b"`, `"null_ab"`. |
| `exposure_type`, `exposure_prevalence` | `"binary"`, `0.5` | `"binary"` or `"continuous"` exposure; P(exposure = 1). |
| `outcome_type`, `outcome_prevalence` | `"continuous"`, `0.5` | `"continuous"` or `"binary"` outcome (effects on the log-odds scale); target P(Y = 1). |
| `sigma` | `1` | Outcome noise SD (continuous). |
| `n_interactions` | `0` | Pairwise mediator interactions in the outcome. |
| `min_variability` | `0.2` | Mediators vary in at least this share of template samples. |
| `n_mc` | `5000` | Monte Carlo draws for the true effects. |
| `stages` | `NULL` | Sequential mediation: views in causal order, e.g. `list("Microbiome", "Metabolome")`. |
| `stage_links` | `list()` | Links between stages: `type` (`"random"`, `"hub"`, `"convergent"`, `"explicit"`), `density` (0.3), `strength` (0.3), `direction` (1), `matrix`. |
| `stage_direct`, `stage_share` | `1`, `NULL` | Per stage: multiplier on direct exposure paths; share of the indirect effect. |

**`null_model`**, **`calibration`**

| Option | Default | Description |
|---|---|---|
| `feature_select` | `"random"` | When `p` is below the template size: `"random"` or `"variance"`. |
| `feature_jitter` | `0.2` | How far offspring features' loadings move from their parent's (0 to 1). |
| `view_models` | `NULL` | Per view: `"continuous"` or `"two_part_copula"` (sparse non-negative microbiome views). |
| `back_transform` | `FALSE` | Return transformed views on the original measurement scale. |
| `grid` | `seq(0, 1, 0.25)` | Candidate crosstalk values for calibration. |
| `nrep`, `nsample` | `5`, template size | Calibration replicates per value and their sample size. |
| `max_features_per_view` | `NULL` | Most variable features used in the calibration comparison. |

### `fit_template()`

| Argument | Default | Description |
|---|---|---|
| `template` | required | `MultiAssayExperiment`, `SummarizedExperiment` or named list of matrices; samples matched by ID (at least 10 shared). Fit repeated measures on one timepoint. |
| `rank` | `NULL` | Factors per view (per view). `NULL`: smallest rank explaining `var.explained`; `"2logp"`: `ceiling(2 * log(p))`. |
| `var.explained` | `0.7` | Share of normal-score variance the factors explain when `rank = NULL`. |
| `max.rank` | `floor(n0 / 3)` | Cap on the per-view rank. |
| `shared.rank` | `NULL` | Number of shared cross-view factors (`NULL`: sequential parallel analysis). |
| `nperm` | `200` | Permutations per parallel-analysis step. |
| `view.transform` | `NULL` | Per view: `"none"`, `"log1p"` or `"clr"`, applied before fitting. |
| `shared.rank.runs` | `1` | Repeat the shared-rank choice and use the rounded mean (stability selection). |
| `seed` | `1234` | Random seed. |

### Other functions

| Function | Arguments (default) |
|---|---|
| `fit_study_templates()` | `x` (one data set per study), `...` (passed to `fit_template()`), `min.features` (2): one template per study on the views and features all studies share |
| `calibrate_crosstalk()` | `template`, `crosstalk.grid` (`seq(0, 1, 0.25)`), `n` (template size), `nrep` (5), `p`, `feature.select`, `feature.jitter`, `max.features.per.view`, `view.models`, `seed` |
| `validate_template()` | `template`, `p`, `n`, `nrep`, `crosstalk.grid`, `max.features.per.view` (100), `view.models`, `template.args`, `seed` |
| `crosstalk_max()` | `template`, `per.view` (`FALSE`) |
| `validate_mimesys_input()` | `x`, `assay` (`NULL`), `min_common_samples` (10) |

### Engines

`gen_simmba()` (outcome), `gen_multistudy()` (several studies),
`gen_mediation()` (mediation) and `sim_from_template()` (null data) take the
settings above as separate arguments with dotted names:
`outcome$model` is `outcome.model`, `type` is `outcome.type`, `prevalence` is
`binary.prevalence`, `signal$prob` is `signal.prob` (likewise `neg_prob`,
`align`, `share` with the `signal.` prefix), `tie`/`de` are `TIE`/`DE`, and
every other option is the same name with `.` for `_` (`cens_rate` is
`cens.rate`, `shift_sd` is `shift.sd`, `n_mediators` is `n.mediators`,
`view_models` is `view.models`). The engines take a numeric `crosstalk`
(default 1; no calibration), `signal.align` defaults to `NULL` (random
positions) rather than 1, and they add `return.template` (default `FALSE`);
they
return lists instead of MAEs; `as_mimesys_mae()` converts a replicate.

## All functions

| Function | What it does |
|---|---|
| `simulate_mimesys()` | Main interface: every design, returns `MultiAssayExperiment`s |
| `validate_mimesys_input()` | Checks samples, features, missing values and constant features before fitting |
| `fit_template()`, `fit_study_templates()` | Fit a template (marginals, per-view factor models, shared cross-view space); one per study on common features |
| `calibrate_crosstalk()`, `validate_template()`, `crosstalk_max()` | Choose crosstalk, check simulated against real correlations and marginals, largest valid crosstalk |
| `as_mimesys_mae()`, `simulation_view_matrices()`, `as_mimesys_view_list()` | Conversions between simulation results, MAEs and matrices |
| `blockwise_feature_correlation()`, `compare_feature_correlation_structure()`, `summarize_simulation()`, `plot_*()` | Diagnostics and plots |

## Good to know

- **Scale.** Templates are fitted on the transformed scale (for example
  `log1p`), and data are returned on that scale unless
  `null_model = list(back_transform = TRUE)`; outcomes and true effects are
  defined on it.
- **Truth.** Feature-level truth is in `rowData()`, sample-level truth
  (outcome, predictor, study, split) in `colData()`, and settings, calibration
  and latent truth in `metadata()`.
- **Privacy.** A fitted template contains the source data; treat it like the
  source data. Simulation results contain only its provenance.

## Documentation

- Walkthrough: `vignette("mimesys-workflow", package = "MimeSys")`
- Algorithms and equations: [`docs/MimeSys_algorithms.md`](docs/MimeSys_algorithms.md)
  ([PDF](docs/MimeSys_algorithms.pdf))
- Example data: `?stelzer_mae` (CyTOF, metabolome, proteome; 53 subjects) and
  `?franzosa_mae` (gut microbiome and metabolome; 155 samples)
- Changes: `NEWS.md`

## License

MIT; see `LICENSE`.
