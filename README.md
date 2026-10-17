# MimeSys

MimeSys is a template-based simulator for multi-omics data generation and validation. It learns layer-specific distributions and cross-layer dependence from an empirical dataset, then generates synthetic multi-view data with controlled sample size, feature count, cross-omics correlation, differential-effect alignment, and outcome-generation structure.

The current package is designed around two linked goals:

1. Generate realistic null multi-omics data that preserves important empirical features of the source dataset.
2. Spike outcomes in controlled ways so downstream method evaluation can be performed outside the core package.

## Installation

```r
# Install the development version from GitHub
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("himelmallick/Mimesys")

# Load the package
library(MimeSys)
```

Two small example datasets ship with the package (see `?stelzer_mae` and
`?franzosa_mae`), and `vignette("mimesys-workflow", package = "MimeSys")` walks
through a complete, runnable workflow.

For local development from the package root:

```r
pak::local_install(".")
# or, during active development
devtools::load_all(".")
```

MimeSys imports `MultiAssayExperiment`, `SummarizedExperiment`, `S4Vectors` and `withr`. Plotting helpers use `ggplot2` when installed.

## Basic Workflow: MultiAssayExperiment

MimeSys is designed for Bioconductor containers. The recommended workflow uses a
`MultiAssayExperiment` (MAE), with one `SummarizedExperiment` per omics view.
Assays follow the usual Bioconductor orientation: features in rows and samples
in columns.

```r
library(MimeSys)
library(MultiAssayExperiment)

# Example data: CyTOF, metabolome and proteome for 53 subjects
data(stelzer_mae)

# Fit once and reuse the template across simulation settings.
template <- fit_template(
  stelzer_mae,
  view.transform = c(Metabolomics = "log1p", Proteomics = "log1p"),
  seed = 1
)

# Simulate one MAE. Defaults calibrate crosstalk, use aligned spike-ins and a
# continuous outcome, and return each view on the scale it was fitted on
# (here log1p for Metabolomics and Proteomics; see back_transform below).
# Pass signal = list(share = ...) to balance the signal across views.
sim_mae <- simulate_mimesys(
  template = template,
  n = 200,
  p = c(CYTOF = 30, Metabolomics = 30, Proteomics = 30),
  outcome = "linear",
  signal = list(share = c(1, 1, 1)),
  seed = 1
)

# Standard Bioconductor accessors
experiments(sim_mae)
colData(sim_mae)[, c("Y", "split")]
rowData(experiments(sim_mae)[["Proteomics"]])
metadata(sim_mae)$simulation
```

`simulate_mimesys()` returns a `MultiAssayExperiment` for one replicate. Its
simulated values are stored in the `"simulated"` assay of each experiment.
Outcome and train/test membership are in `colData()`, feature-level spike-in
truth is in each experiment's `rowData()`, and complete simulation settings and
latent truth are in `metadata()`.

Start by checking an empirical input before fitting a template:

```r
input_report <- validate_mimesys_input(stelzer_mae)
input_report
attr(input_report, "common_samples")
```

This reports samples, features, missing values, constant features, and the
number of sample IDs shared across views. `fit_template()` will stop on
structurally incompatible inputs, while this report makes the reason visible
before fitting.

By default, MimeSys first calibrates crosstalk on simulated **null** data. It
compares observed and simulated cross-view feature-correlation blocks over a
small grid, selects the value with the lowest reconstruction error averaged
over all cross-view blocks,
and then uses that value for the requested outcome simulation. The complete
calibration result is retained in
`metadata(sim_mae)$simulation$crosstalk_calibration`. Supply a numeric
`crosstalk` only when intentionally overriding that empirical calibration.
Calibration uses the template's original sample size and all matched empirical
features by default; set `calibration$nsample` or
`calibration$max_features_per_view` only when a faster approximation is needed.

Use `simulate_mimesys_null()` when only null multi-omics features are needed:

```r
null_mae <- simulate_mimesys_null(
  template,
  n = 200,
  p = 200,
  seed = 1
)
```

This returns an MAE with no `Y`, no spike-in coefficients, and the selected
crosstalk plus calibration diagnostics in `metadata(null_mae)$simulation`.

Every fitted template, calibration result, and MAE simulation also records the
MimeSys version, UTC creation time, random seed, selected crosstalk, and source
template provenance. For a simulated MAE, use
`metadata(sim_mae)$simulation$provenance`.

Use compact lists only when you need to depart from defaults:

Shared-latent Friedman outcomes need at least five shared factors. The example
template above has one (chosen by parallel analysis), so fix `shared.rank` for
these scenarios:

```r
template5 <- fit_template(
  stelzer_mae,
  view.transform = c(Metabolomics = "log1p", Proteomics = "log1p"),
  shared.rank = 5
)

sim_mae <- simulate_mimesys(
  template5,
  n = 200,
  p = 200,
  crosstalk = 0.75,
  outcome = list(
    model = "shared_latent",
    shared_latent_function = "friedman2",
    shared_latent_n_factors = 5
  ),
  signal = list(
    alignment = 1,
    share = c(CYTOF = 1, Metabolomics = 1, Proteomics = 1)
  ),
  seed = 1
)
```

Built-in presets make common outcome scenarios explicit:

```r
mimesys_presets()

sim_mae <- simulate_mimesys(
  template,
  n = 200,
  p = 200,
  preset = "interaction",
  seed = 1
)

sim_mae <- simulate_mimesys(
  template5,
  n = 200,
  p = 200,
  preset = "shared_latent_friedman2",
  seed = 1
)
```

Available presets are `"linear"`, `"friedman"`, `"friedman2"`,
`"shared_latent"`, `"shared_latent_friedman"`,
`"shared_latent_friedman2"`, and `"interaction"`. The default is `"default"`, an aligned linear
outcome. Explicit fields in `outcome`, `signal`, or `null_model` override their
preset values.

For `nrep > 1`, the function returns a `mimesys_mae_simulation` object whose
`$replicates` element is a list of MAEs, one per replicate. Set `return_type = "simulation"`
to access the lower-level legacy result from `gen_simmba()`.
`simulate_mimesys_null()` takes the same `nrep` and `return_type` arguments.

Every simulator, `fit_template()`, `calibrate_crosstalk()` and
`validate_template()` use the same default `seed = 1234` and restore the
caller's random-number state afterwards; pass `seed = NULL` to draw from the
current random-number stream instead.

## Input Data

`fit_template()` accepts a `MultiAssayExperiment` directly. It also accepts a
single `SummarizedExperiment` or a named list of matrices for interoperability.
Matrices use the internal MimeSys orientation, samples in rows and features in
columns; MAE and SummarizedExperiment assays are transposed automatically.

```r
views <- list(
  CYTOF        = matrix_1,
  Metabolomics = matrix_2,
  Proteomics = matrix_3
)
```

Samples are matched across views by sample ID: for an MAE, the primary sample
IDs from its `sampleMap()` (so experiments may use different column names for
the same sample); for a list of matrices, the row names of each matrix.
Samples not present in all views are dropped.

The template assumes **one sample per subject**. For repeated-measures or
longitudinal data, use a single timepoint (or one sample per subject) to fit the
template. Treating several samples from the same subject as independent inflates
the estimated cross-view structure; for example, fitting all 150 samples of the
Stelzer data (3 timepoints for each of 53 subjects, all features) finds 18
shared factors, while one timepoint per subject finds 2.

Names matter. View names are carried through the template, the simulation output, and the outcome-generation metadata. They are also used by arguments such as `p`, `signal.share`, `friedman.share`, `view.models`, and named `view.transform` vectors.

## Fitting A Template

`fit_template()` learns the empirical structure that will later be used to simulate new data.

```r
template <- fit_template(
  template = views,
  rank = NULL,
  var.explained = 0.7,
  max.rank = NULL,
  shared.rank = NULL,
  nperm = 200,
  shared.rank.runs = 1,
  view.transform = NULL,
  seed = 1234
)
```

### `template`

A named list of omics matrices. Rows are samples and columns are features. The argument is called `template` because these empirical matrices are used as the template for null-data generation.

### `rank`

The number of latent factors to estimate from the empirical data.

If `rank` is supplied, the template uses that number of factors (one value,
one per view, or a vector named by view), capped at `max.rank` and at one fewer
than the number of features in the view. If
`rank = NULL`, MimeSys chooses the number of factors using `var.explained`.
Alternatively, `rank = "2logp"` uses `ceiling(2 * log(p_m))` factors for view
`m`, subject to `max.rank`. This is a useful high-dimensional default when a
fixed, feature-count-aware rank is preferred over variance-threshold selection.

### `var.explained`

The target cumulative variance explained when `rank = NULL`.

For example, `var.explained = 0.7` asks the template to retain enough latent factors to explain approximately 70% of the variation of each view's normal scores (chosen
separately for each view).

This argument is ignored when `rank` is explicitly supplied.

### `view.transform`

Optional transformation applied before fitting the template.

Supported values are:

- `"none"`: use the matrix as supplied.
- `"log1p"`: apply `log1p(x)` to each value.
- `"clr"`: apply centered log-ratio transformation, usually for compositional microbiome-like data.
- `NULL`: equivalent to `"none"` for all views.

This can be a single value used for all views:

```r
template <- fit_template(
  views,
  view.transform = "log1p"
)
```

It can also be a named vector giving separate transformations by view:

```r
template <- fit_template(
  views,
  view.transform = c(
    Microbiome = "clr",
    Metabolome = "log1p",
    Proteome = "none"
  )
)
```

Views not listed in a named `view.transform` vector use `"none"`.

The template is fitted, calibrated and validated on the transformed scale, and
by default every simulator returns that fitted scale. Simulated data can be
returned on the original measurement scale with `back.transform = TRUE` in
`sim_from_template()`, `gen_simmba()` and `gen_mediation()`, or
`null_model = list(back_transform = TRUE)` in `simulate_mimesys()`,
`simulate_mimesys_null()` and `simulate_mediation()`. Each
simulated feature's latent value is mapped through that feature's
original-scale empirical quantiles: the exact inverse for `log1p` (integer data
stay integers) and, for `clr`, each feature's original distribution including
its zeros. Compositional views never exceed the template total (see the
microbiome section). Outcomes are always generated from the transformed scale.

### `max.rank`

Upper bound on the number of factors retained per view. It also caps an
explicitly supplied `rank`.

If `NULL`, the default cap is approximately one third of the empirical sample size.

### `shared.rank`

Number of shared latent variates to retain across views.

If `NULL`, MimeSys chooses the shared rank by sequential (deflated) parallel
analysis. At each step, the largest generalized-CCA eigenvalue of the current
view scores is compared with the permutation distribution of the largest
eigenvalue after permuting samples independently within each view. If it is
significant (p < 0.05), that shared direction is removed from every view and
the test is repeated; the shared rank is the number of significant steps. This
avoids treating finite-sample noise as cross-view structure. (Comparing each
eigenvalue with the null of the same-index eigenvalue over-selects, because one
real shared factor lifts all later observed eigenvalues.)

### `nperm`

Number of permutations used when `shared.rank = NULL`.

Larger values give a more stable null threshold for shared-rank selection but take longer (default 200). The p-value, null threshold and tested eigenvalue of each step are stored in `template$gcca.pvalue`, `template$gcca.threshold` and `template$gcca.step.eigen`.

### `shared.rank.runs`

Number of stability-selection runs when `shared.rank = NULL` (default 1). The parallel-analysis choice depends on random tie-breaking and permutations; with more than one run it is repeated this many times and the rounded mean is used. The runs and their distribution are stored in `template$shared.rank.selection`, and a message is shown when fewer than half the runs agree with the chosen rank. A single run warns when the choice is borderline; use, for example, 20 runs to stabilize such a choice. Fixing `shared.rank` skips the selection.

### `seed`

Random seed for a reproducible fit (default `1234`, as in every MimeSys function). The caller's random-number state is restored afterwards; `NULL` uses the current random-number stream.

## What The Template Stores

The fitted template stores:

- View names, feature names and transformation choices.
- Each feature's sorted values (its empirical quantile function), on the fitted
  scale and, for transformed views, on the original scale.
- Two-part (zero / positive) marginal summaries for sparse views.
- Within-view factor loadings and ranks.
- The shared cross-view loading matrices, the GCCA eigenvalues, and the
  per-step permutation thresholds and p-values used to choose the shared rank.
- The matrices the template was fitted on (`observed`), used by calibration
  and validation.
- Provenance (package version, fitting settings, input dimensions).

Because it keeps the fitted matrices and each feature's values, a fitted
template should be handled with the same care as the source data.

The template is intentionally reusable. Fit it once from the empirical dataset, then generate many synthetic datasets across a grid of sample sizes, feature counts, crosstalk values, and outcome-generation scenarios.

## Generating Simulation Datasets

`gen_simmba()` is the high-level simulation function. It calls the null simulator, spikes feature effects, generates outcomes, and splits each replicate into training and test data.

```r
sim <- gen_simmba(
  nsample = 200,
  template = template,
  p = 200,
  crosstalk = 0.5,
  feature.select = "random",
  feature.jitter = 0.2,
  view.models = NULL,
  snr = 1,
  p.train = 0.7,
  de.prob = 0.1,
  de.downProb = 0.5,
  de.facLoc = 1,
  de.facScale = 0.4,
  de.align = NULL,
  signal.share = NULL,
  friedman.share = NULL,
  ygen.mode = "LM",
  outcome.type = "continuous",
  binary.prevalence = 0.5,
  surv.hscale = 1,
  cens.rate = NULL,
  cens.lower = 1,
  cens.upper = 3,
  nrep = 1,
  seed = 1234,
  return.template = FALSE,
  back.transform = FALSE
)
```

`sim_from_template()` is the lower-level null feature generator. Use it directly when you only need synthetic omics matrices without outcome generation or train/test splitting.

```r
pcl <- sim_from_template(
  tmpl = template,
  n = 200,
  p = 200,
  crosstalk = 0.5,
  feature.select = "random",
  feature.jitter = 0.2,
  view.models = NULL
)
```

The sections below describe the high-level `gen_simmba()` arguments. Arguments that are passed directly to `sim_from_template()` are noted.

### `nsample`

Number of samples to simulate in each replicate before train/test splitting.

This can be smaller, equal to, or larger than the empirical sample size. In simulation studies, this is usually varied to test whether downstream results stabilize as sample size increases.

### `template`

Either a fitted `mimesys_template` returned by `fit_template()` or a raw named list of sample-by-feature matrices.

Using a fitted template is recommended when running grids:

```r
template <- fit_template(views, rank = 10)
sim <- gen_simmba(nsample = 200, template = template)
```

If a raw list is supplied, `gen_simmba()` fits the template internally and passes extra arguments through `...` to `fit_template()`.

### `p`

Number of features to simulate per omics layer.

`p` can be supplied in several forms:

```r
# Same number of features in every layer
p = 200

# Unnamed vector in template view order
p = c(100, 200, 500)

# Named vector
p = c(CYTOF = 100, Metabolomics = 200, Proteomics = 500)

# Named list
p = list(CYTOF = 100, Metabolomics = 200, Proteomics = 500)
```

If `p = NULL`, MimeSys uses the original number of features available in each layer of the template.

### `crosstalk`

Controls the strength of cross-omics dependence during null data generation.

Conceptually:

- Low `crosstalk` means the views are more independent.
- High `crosstalk` means the views share more latent structure and are more correlated across layers.

This is a simulation-control parameter, not an automatically estimated biological truth. It changes how strongly latent factors are shared across layers. Empirically, higher crosstalk should generally increase cross-view correlation, but high crosstalk can create redundancy rather than complementarity.

`crosstalk` can be a scalar shared by all views or a vector with one value per view. Values must be non-negative and cannot exceed `crosstalk_max(template)`.

### `feature.select`

Feature-selection rule used when the requested feature count is smaller than the number available in the template.

Supported values:

- `"random"`: randomly sample template features.
- `"variance"`: choose the highest-variance empirical features first (variance on the fitted scale, after any `view.transform`).

This argument is passed to `sim_from_template()`.

### `feature.jitter`

Amount of random perturbation applied when generating offspring features beyond the empirical template feature count.

Use `0` for closer copies of parent features and larger values for more feature-level variation. Typical values are between `0` and `1`.

This argument is passed to `sim_from_template()`.

### `view.models`

Optional named vector controlling the observation model for each view.

Supported values:

- `"continuous"`: standard empirical-quantile emission.
- `"two_part_copula"`: two-stage sparse microbiome-like emission. It is for
  sparse non-negative abundances only (zeros mean absence): views with negative
  values or fitted with `view.transform = "clr"` give an error, and a view with
  no zeros gives a warning. Dense views such as the Stelzer data should use
  `"continuous"`.

Example:

```r
view.models = c(
  Microbiome = "two_part_copula",
  Metabolome = "continuous"
)
```

Unspecified views use `"continuous"`.

### `snr`

Signal-to-noise ratio for outcome generation.

Higher `snr` creates easier prediction tasks. Lower `snr` creates harder tasks.

### `p.train`

Proportion of samples assigned to the training split.

For example, `p.train = 0.7` creates a 70% training and 30% test split within each replicate.

### `de.prob`

Probability that a feature is assigned a nonzero spike-in effect.

Can be a scalar or one value per view. Every view with `de.prob > 0` gets at
least one nonzero coefficient; use `0` for a view that should carry no linear
signal.

### `de.downProb`

Probability that a nonzero spiked feature has a negative effect.

Can be a scalar or one value per view.

### `de.facLoc`

Location parameter for the log-normal effect-size multiplier.

Can be a scalar or one value per view.

### `de.facScale`

Scale parameter for the log-normal effect-size multiplier.

Can be a scalar or one value per view.

### `de.align`

Controls how outcome-relevant feature effects are oriented relative to the shared cross-view structure.

`NULL` uses the original independent spike-in behavior. A numeric value between `0` and `1` places that fraction of differential features on features that load strongly on the shared cross-view direction.

More detail is provided in the differential-effect alignment section below.

### `signal.share`

Optional numeric vector controlling the relative marginal variance contribution of each view-specific linear predictor.

The linear predictor is computed on standardized features: within each
replicate every simulated feature is centred and scaled to SD 1 before it is
multiplied by its coefficient, so coefficients are per-SD effects and views on
different measurement scales are comparable.

If `NULL`, each view contributes according to its own coefficients (so a view
with more nonzero coefficients or larger effects contributes more).

If supplied, values are recycled or matched by view name, normalized to sum to one, and each view's coefficients are rescaled so that its view-specific linear predictor has variance equal to its share. For example, `c(1, 1, 1)` gives each of three views an equal share of the linear signal.

### `friedman.share`

Optional numeric vector controlling which views supply the five features that enter the nonlinear Friedman function when `ygen.mode` is `"Friedman"` or `"Friedman2"`. Each of the five is drawn from a view chosen with probability equal to its share, then uniformly among that view's nonzero spiked features.

If `NULL`, MimeSys uses `signal.share` when available; otherwise Friedman features are sampled uniformly from all nonzero spiked features.

### `ygen.mode`

Outcome-generation mode.

Supported values:

- `"LM"`
- `"Friedman"`
- `"Friedman2"`
- `"SharedLatent"`
- `"Interaction"`

These are described in detail below. `nonlinear.weight` (default `0.5`) sets the
nonlinear share for `"Friedman2"` and `shared.latent.fun = "friedman2"`;
`interaction.weight` (default `0.5`) and `n.interactions` control
`"Interaction"`.

### `outcome.type`

Observed outcome type.

Supported values:

- `"continuous"`
- `"binary"`
- `"survival"`

### `binary.prevalence`

Target proportion of cases for binary outcomes (default `0.5`). The logistic intercept is solved in every replicate so that the expected prevalence matches this value.

### `surv.hscale`

Multiplicative scale on the baseline hazard for survival outcomes.

### `cens.rate`

Optional target censoring proportion for survival outcomes. When supplied, censoring times are uniform on `(0, c)` with `c` solved so that this proportion of subjects is censored; `cens.rate = 0` gives no censoring. When `NULL` (default), `cens.lower` and `cens.upper` are used.

### `cens.lower` and `cens.upper`

Lower and upper bounds for uniform censoring times in survival simulations when `cens.rate` is `NULL`.

### `nrep`

Number of independent simulation replicates (default `1`).

Each replicate contains a training dataset and a test dataset generated under the same scenario.

### `seed`

Random seed for reproducible simulation (default `1234`). The caller's
random-number state is restored afterwards.

### `return.template`

Whether to include the fitted template in the result (default `FALSE`). A fitted template stores the source data it was fitted on (`template$observed` and each feature's sorted values), so it should be handled with the same care as the source data. By default `gen_simmba()` returns only `template_provenance`, so the simulation result can be shared without exposing the source data.

### `...`

Additional arguments passed to `fit_template()` when `template` is a raw list of matrices rather than an already fitted template.

## Lower-Level Null Simulation

Use `sim_from_template()` when you only need null feature matrices and do not want outcomes.

```r
pcl <- sim_from_template(
  tmpl = template,
  n = 200,
  p = 200,
  crosstalk = 0.5,
  view.models = c(Microbiome = "two_part_copula")
)
```

The returned object is a PCL-style list with:

- `feature_table`: features by samples.
- `sample_metadata`: sample-level metadata.
- `feature_metadata`: feature-level metadata including feature type, parent template feature, offspring status, view model, shared loading, and the output scale (`"fitted"` or `"original"`).
- `truth`: latent scores, shared variates, crosstalk weights, canonical correlations, and view models.

## Crosstalk Calibration

If you want crosstalk to be selected by the simulator instead of hand-picked, use calibration.

```r
cal <- calibrate_crosstalk(
  template,
  crosstalk_grid = c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1),
  nsample = 200,
  p = 200,
  nrep = 20,
  seed = 1
)

plot_crosstalk_calibration(cal)
```

`calibrate_crosstalk()` runs the null simulator over a grid of crosstalk values and summarizes the resulting cross-view correlation behavior. `plot_crosstalk_calibration()` visualizes whether increasing crosstalk produces the expected increase in empirical cross-view correlation.

Calibration includes uncertainty diagnostics. `cal$selection_frequency` gives
the fraction of null replicates selecting each candidate value, and
`cal$near_optimal` lists values within one standard error of the smallest mean
cross-block distance. These make it clear when the selected value is decisive
versus when several crosstalk values reconstruct the template equally well.

Calibration always compares simulations with the matrices stored in the fitted
template, so data fitted with a `view.transform` are calibrated on the
transformed scale. Grid values above `crosstalk_max(template)` are dropped with
a warning, and a warning is given when the best value sits at the top of a grid
that stops well short of `crosstalk_max(template)`.

## Template validation report

`validate_template()` creates one compact report before a simulation study. It
calibrates crosstalk, simulates a null dataset at the selected value, and checks
cross-block reconstruction error, correlation-pattern preservation, and
marginal distribution preservation (per-feature Kolmogorov-Smirnov distance and SD-scaled quantile error, which are comparable across views).

```r
report <- validate_template(
  template,
  p = c(CYTOF = 200, Metabolomics = 200, Proteomics = 200),
  nrep = 10,
  crosstalk_grid = seq(0, 1, by = 0.1),
  seed = 1
)

report
plot_template_validation(report, "calibration")
plot_template_validation(report, "correlation_structure")
plot_template_validation(report, "marginal")
```

The calibration objective uses the rescaled Frobenius RMSE between matching
empirical and simulated **cross-view** feature-correlation blocks. It therefore
requires shared feature names; simulated offspring features are not used in this
comparison.

Useful arguments:

- `crosstalk_grid`: values to evaluate.
- `nsample`: number of samples in each simulated replicate.
- `p`: number of features per layer.
- `nrep`: number of simulation replicates.
- `seed`: reproducibility.

## Differential-Effect Alignment

`de.align` controls how outcome-relevant feature effects are oriented relative to the shared cross-view structure.

This argument is most important for linear and additive outcome generation.

### `de.align = 1`

Outcome effects are aligned with the shared latent direction. This means signal is more consistent with the latent structure that links views together.

This can help integration methods when multiple views contain useful signal in a coordinated way.

### `de.align = 0`

Outcome effects are randomly oriented rather than aligned with the shared latent direction.

This creates signal that is less tied to the cross-view structure.

### Fractional values

The implementation supports numeric values between `0` and `1`. Fractional values interpolate between random and aligned signal directions.

In many simulation grids, the interpretable endpoints are:

- `0`: random effect orientation.
- `1`: aligned effect orientation.

## Signal Sharing

`signal.share` controls how much of the outcome signal is assigned to each omics layer.

```r
signal.share = c(CYTOF = 1, Metabolomics = 1, Proteomics = 1)
```

The values are normalized internally, so the above means that the outcome signal is balanced across the three views.

Single-layer scenarios can be created by assigning all signal to one view:

```r
signal.share = c(CYTOF = 1, Metabolomics = 0, Proteomics = 0)
```

This is useful for validating that single-layer methods win when only one layer carries the outcome signal.

Balanced scenarios are useful for generating outcomes where signal is distributed across layers:

```r
signal.share = c(CYTOF = 1, Metabolomics = 1, Proteomics = 1)
```

Important distinction:

- `crosstalk` controls feature-level cross-view dependence.
- `signal.share` controls outcome-signal allocation across views.
- `de.align` controls whether the spiked outcome signal is aligned with shared latent structure.

These are separate simulation knobs.

## Outcome Generation

MimeSys currently supports five outcome-generation modes through `ygen.mode`.

```r
ygen.mode = "LM"
ygen.mode = "Friedman"
ygen.mode = "Friedman2"
ygen.mode = "SharedLatent"
ygen.mode = "Interaction"
```

Every mode is computed from standardized features. The ground truth for which
features drive the outcome is `feature_metadata$outcome_feature`. `beta` is
nonzero only for features that enter through the linear predictor, so it is
zero for every feature in `"Friedman"`, `"SharedLatent"` and pure-interaction
(`interaction.weight = 1`) outcomes. `friedman_term` and `interaction_term`
mark the nonlinear and interaction features. In `"SharedLatent"` the outcome
depends on the shared latent factors (`truth$shared`) rather than on individual
features, so no feature is marked; components with weight 0 are not marked.

### `LM`

Linear additive outcome generation.

For each layer \(m\), MimeSys constructs a layer-specific linear predictor:

```text
eta_m = X_m beta_m
```

The final latent outcome score is an additive weighted sum across layers:

```text
eta = sum_m w_m eta_m
```

where `w_m` comes from the normalized `signal.share` vector.

Noise is added according to `snr`.

This mode is useful for scenarios where signal is additive across views.

### `Friedman`

Nonlinear Friedman-style outcome generation.

MimeSys first converts selected layer-level predictors to rank-based unit-scale variables and then applies a Friedman-type nonlinear function involving main effects and interactions.

This mode is useful when the outcome is not well represented by a simple linear sum of features.

### `Friedman2`

Hybrid additive plus nonlinear Friedman outcome generation.

This mode mixes the additive layer signal with a Friedman-style nonlinear component. Both parts are standardized and combined as

```text
eta = sqrt(1 - w) * linear + sqrt(w) * friedman,   w = nonlinear.weight (default 0.5)
```

so `w` is roughly the share of signal from the nonlinear part (exactly, when the two parts are uncorrelated; the Friedman features are drawn from the spiked features, so the parts can be correlated). It is useful for testing methods under partially linear and partially nonlinear outcome structure.

### `SharedLatent`

Outcome generation from the shared latent structure.

MimeSys uses the shared latent factors generated during null data simulation. Instead of using only one factor, the current implementation combines all available shared factors by standardizing them and summing across factors:

```text
eta_shared = rowSums(scale(shared_factors)) / sqrt(number_of_factors)
```

The outcome is then generated from this shared latent score with noise controlled by `snr`. By default, `SharedLatent` uses all available shared factors and a linear equal-weight score. You can control this more directly with:

```r
shared.latent.fun = "linear"     # default
shared.latent.fun = "friedman"   # nonlinear Friedman response on shared factors
shared.latent.fun = "friedman2"  # linear shared score mixed with Friedman response (nonlinear.weight)
shared.latent.n_factors = 5      # first K shared factors to use
```

The nonlinear shared-latent functions rank-transform the selected shared factors to the unit interval before applying the Friedman response. They require at least five selected shared factors, so fit the template with `shared.rank >= 5` when using `shared.latent.fun = "friedman"` or `"friedman2"`.

This mode is useful when the true outcome is driven by a cross-view latent process rather than view-specific feature effects. With `crosstalk = 0` the shared factors do not enter any simulated feature, so the outcome is unrelated to the features; MimeSys warns in that case.

### `Interaction`

Cross-view interaction outcome. MimeSys picks `n.interactions` feature pairs (default: one per pair of views), each pair taking one feature from each of two different views (preferring spiked features), and forms

```text
eta_int = sum_k x_{a_k} * x_{b_k}        (standardized features)
eta     = sqrt(1 - w) * linear + sqrt(w) * eta_int,   w = interaction.weight (default 0.5)
```

With `interaction.weight = 1` the outcome has no main effects, so only methods that model features from different views jointly can recover it. This is the most direct test of whether integrating views helps. The pairs are stored in `truth$interaction_pairs`.

## Outcome Types

`outcome.type` controls the final observed outcome.

```r
outcome.type = "continuous"
outcome.type = "binary"
outcome.type = "survival"
```

### Continuous

The latent score is returned as a continuous outcome with added noise.

The outcome is stored in sample metadata as `Y`.

### Binary

The linear predictor is standardized, scaled so that `snr` equals `Var(eta) / (pi^2 / 3)` (the latent-logistic signal-to-noise ratio), and shifted by an intercept solved to hit `binary.prevalence`. Outcomes are sampled as `Y ~ Bernoulli(logistic(eta))`.

This is useful for classification-style simulations.

### Survival

The linear predictor is standardized and scaled so that `snr` equals `Var(eta) / (pi^2 / 6)` (the signal-to-noise ratio on the log event-time scale). Event times are exponential with hazard `surv.hscale * exp(eta)`, censoring follows `cens.rate` or `cens.lower`/`cens.upper`, and sample metadata stores `time` and `status`.

This is useful for survival-style simulations.

## Signal-To-Noise Ratio

`snr` controls outcome signal strength relative to noise. For continuous outcomes it is `Var(eta) / Var(noise)`; for binary and survival outcomes it has the same meaning on the latent scale (see Outcome Types).

Higher `snr` creates easier prediction tasks. Lower `snr` creates harder tasks.

For simulation studies, it is often useful to keep `snr` fixed while varying sample size, feature count, crosstalk, and signal-sharing scenario.

## Microbiome Two-Stage Emission Model

MimeSys supports microbiome-like data generation through a two-stage emission model. This is used when a view is treated as sparse, compositional, or count-like.

The key idea is that microbiome data have two linked parts:

1. Whether a feature is present or absent.
2. How abundant it is when present.

The simulator separates these two processes.

### Stage 1: Presence Or Absence

For each sample \(i\) and feature \(j\), the simulator generates a presence indicator:

```text
Z_ij = 1 if feature j is present in sample i
Z_ij = 0 otherwise
```

The probability of presence can depend on the simulated latent structure and feature-specific parameters learned from the empirical template.

This step controls sparsity and zero inflation.

### Stage 2: Positive Abundance

Conditional on presence, MimeSys generates a positive abundance value:

```text
A_ij | Z_ij = 1
```

This positive component is generated from the layer-specific latent structure and empirical feature behavior.

The observed value is:

```text
X_ij = Z_ij * A_ij
```

So if the feature is absent, the observed value is zero. If the feature is present, the observed value is the generated positive abundance.

Presence is decided by comparing the feature's latent uniform with its empirical zero probability `pi_j`. Given presence, the abundance uniform is rescaled to `(u - pi_j) / (1 - pi_j)` before it is mapped through the feature's positive-value quantiles, so positive abundances reproduce the full positive distribution.

Composition views never exceed the template total. If every template sample
sums to the same total (a closed composition, such as relative abundances
summing to 1) and at least the template's features are simulated, including
offspring features, every simulated sample is renormalized to that total, so
offspring features dilute the others. If fewer features are simulated, samples
above the total are scaled down to it. If the template is itself a subset of a
composition (values in [0, 1] with sample totals at most 1), samples above the
largest observed total are scaled down to it. Other views are unconstrained.

### Why This Matters

This model better reflects microbiome-like data than a single Gaussian emission because it preserves:

- Sparsity.
- Zero inflation.
- Positive abundance structure.
- Dependence between presence and abundance.
- Cross-view dependence through the shared latent simulator.

The two-stage emission is part of null feature generation. It is separate from outcome generation, which is controlled by `ygen.mode`.

## Mediation Simulation

`gen_mediation()` and `simulate_mediation()` add an exposure, mediating
features and an outcome on top of the template's null data, for benchmarking
single- and multi-view mediation methods. The design is domain-agnostic: any
view and any feature can be a mediator.

```r
med_mae <- simulate_mediation(
  template,
  n = 300,
  mediation = list(
    n_mediators = c(CYTOF = 3, Metabolomics = 3, Proteomics = 0),
    tie = 1,                 # target total (natural) indirect effect
    de = 0.5,                # direct effect
    dir_a = 1, dir_b = 1,    # path directions (1, -1 or "random")
    effect_type = "spiked",  # or "null_a", "null_b", "null_ab"
    exposure_type = "binary",
    outcome_type = "continuous",  # or "binary" (TIE and DE on the log-odds scale)
    n_interactions = 0       # pairwise mediator interactions in the outcome
  )
)
colData(med_mae)[, c("exposure", "Y")]
metadata(med_mae)$truth$mediation[c("NIE", "NDE", "TE", "NIE_by_view")]
```

How it works:

- **Exposure:** `T` is binary (`exposure_prevalence`) or standard normal.
- **a-paths:** each mediator's latent normal score is mixed with the
  standardized exposure, `Z_j <- a_j T* + sqrt(1 - a_j^2) Z_j` (for a binary
  exposure this mixture is mapped back to an exact standard normal), before
  mapping through the feature's empirical quantiles. Mediators keep their
  template distribution (zeros, bounds, skew); the exposure shifts where
  samples fall in it, so for sparse features it changes both presence and
  abundance.
- **Closed compositions:** renormalizing a sample after its mediators shift
  changes the other features of that view too. `rowData()$exposure_effect`
  gives the true exposure effect on every feature, so non-mediators affected
  only through closure can be accounted for in mediator-selection benchmarks.
- **b-paths:** `Y = DE * T + sum_j b_j M~_j + sum_k g_k M~_k1 M~_k2 + sigma * e`,
  with mediators standardized by their template mean and SD.
- **Calibration and truth:** b-paths are rescaled so that the population
  natural indirect effect, computed by paired counterfactual Monte Carlo
  (`n_mc` draws), equals `tie * dir_a * dir_b`. The true NIE, NDE and TE, the
  indirect effect carried by each view, the interaction part, and each
  mediator's latent a-path (`a_path`), observed exposure effect (`a_effect`),
  b-path and indirect contribution are returned. Effects are for `T = 1`
  versus `T = 0` (one SD for a continuous exposure), in outcome units.
- **Binary outcome:** `Y ~ Bernoulli(logit^-1(b0 + DE * T + sum_j b_j M~_j))`.
  `tie` and `de` (and the reported `NIE`, `NDE`, `TE`) are conditional effects
  on the linear-predictor (log-odds) scale, and `b0` is solved for
  `outcome_prevalence`. Because odds ratios are non-collapsible, the natural
  counterfactual effects differ; they are reported in `truth$mediation$binary`
  as risk differences, odds ratios and log odds ratios.
- **Scale:** the outcome, paths and truth are defined on the fitted
  (transformed) scale, and `simulate_mediation()` returns features on that
  scale by default (`null_model = list(back_transform = TRUE)` returns the
  original scale, in which case analyse e.g. `log1p` of `log1p`-fitted views).
- **Null scenarios are matched:** in `"null_a"` the b-paths have the same size
  as in the corresponding spiked scenario.
- **Continuous exposure:** mediators respond non-linearly to the exposure, so
  the truth also gives `NIE_linear` (from each mediator's population linear
  regression on the exposure, `a_slope`), which is what product-of-coefficients
  estimators target. For a binary exposure the two coincide.
- **Nulls:** `effect_type = "null_a"` (mediators affect `Y` only),
  `"null_b"` (exposure affects mediators only) and `"null_ab"` give zero
  indirect effect for type I error studies.

### Sequential mediation

Give `stages`, a list of views in causal order, to simulate chains such as
exposure -> microbes -> metabolites -> outcome. Mediators of a later stage are
linked to mediators of the previous stage, many-to-many: one microbe can feed
several metabolites and one metabolite can have several microbe parents.

```r
data(franzosa_mae)
micro_template <- fit_template(franzosa_mae, view.transform = c(Metabolome = "log1p"))

seq_mae <- simulate_mediation(
  micro_template,
  n = 500,
  mediation = list(
    n_mediators = c(Microbiome = 5, Metabolome = 5),
    stages = list(Microbe = "Microbiome", Metabolite = "Metabolome"),
    stage_links = list(
      type = "random",   # or "hub", "convergent", "explicit"
      density = 0.3,     # link probability (or hub fraction / parents per child)
      strength = 0.3,    # share of a metabolite's latent variance from its parents
      direction = 1      # link signs: 1, -1 or "random"
    ),
    stage_direct = 1,    # multiplier on each stage's direct exposure paths
    stage_share = NULL,  # optional share of the indirect effect per stage
    tie = 1, de = 0.5
  )
)
truth <- metadata(seq_mae)$truth$mediation
truth$path_effects   # T -> Metabolite -> Y, T -> Microbe -> Metabolite -> Y, T -> Microbe -> Y
truth$links          # microbe -> metabolite links with their coefficients
```

- **Model:** a stage-2 mediator's latent score is
  `W_j = a_j T* + sum_k d_jk W_k + c_j Z_j`, with parents `k` in stage 1 and
  `c_j` chosen (allowing for the template crosstalk between the views) so that
  it stays exactly standard normal; mediators keep their template
  distributions.
- **Truth:** natural path-specific effects for every path, computed by
  cross-world paired Monte Carlo and adding up to the NIE (`path_effects`,
  `NIE_by_stage`), per-feature `stage`, `n_parents` and `n_children`, and the
  link list (`links`) for scoring link recovery.
- **Scenarios:** `stage_direct = c(1, 0)` with `stage_share = c(0, 1)` gives
  a pure chain `T -> M1 -> M2 -> Y`; `stage_links = list(strength = 0)` gives
  unlinked mediator sets that are still correlated through crosstalk (a null
  for the link); the `effect_type` nulls apply as before. Lowering
  `crosstalk` shows how much of a method's false links come from it.
- **Data scale:** links act on the latent scores. For skewed or sparse
  parents (e.g. microbiome relative abundances) the observed values are a
  non-linear function of them, so linear estimators attribute part of the
  chain `T -> M1 -> M2 -> Y` to the direct `T -> M2 -> Y` path; transformed
  features (ranks, log ratios) recover it better.
- With a single stage (or `stages = NULL`) the mediators are parallel and the
  stage-specific outputs (`path_effects`, `links`, `stage`) are omitted.

## Additive Versus Shared-Latent Outcome Generation

Two outcome-generation modes are especially important for current simulation studies.

### Additive Outcome Generation

Additive outcomes use the layer-specific feature matrices directly:

```text
eta = w_1 X_1 beta_1 + w_2 X_2 beta_2 + ... + w_M X_M beta_M
```

This supports scenarios where each layer contributes its own piece of the outcome.

Expected behavior:

- Single-layer signal: best single layer should perform well.
- High crosstalk may not always help, because shared information can become redundant.

### Shared-Latent Outcome Generation

Shared-latent outcomes use the common latent factors that couple the views:

```text
eta = f(shared latent factors)
```

In the current implementation:

```text
eta = rowSums(scale(shared_factors)) / sqrt(number_of_factors)
```

Expected behavior:

- The outcome is intrinsically multi-view because it is generated from the structure shared across layers.

## Inspecting Simulated Data

The simulation object contains replicate-level train and test data.

```r
sim$trainDat[[1]]
sim$testDat[[1]]
```

Within each replicate, `feature_table` is stored as features by samples. To recover sample-by-feature matrices split by omics view, use `simulation_view_matrices()`.

```r
train_views <- simulation_view_matrices(
  sim,
  rep = 1,
  split = "train",
  feature_names = "featureID"
)

names(train_views)
train_views$CYTOF
train_views$Metabolomics
train_views$Proteomics
```

Sample metadata contains outcomes.

```r
sim$trainDat[[1]]$sample_metadata
sim$testDat[[1]]$sample_metadata
```

For continuous and binary outcomes, the main outcome column is usually `Y`.

```r
y_train <- sim$trainDat[[1]]$sample_metadata$Y
y_test  <- sim$testDat[[1]]$sample_metadata$Y
```

## Summarizing A Simulation

Use `summarize_simulation()` to inspect the generated simulation.

```r
sim_summary <- summarize_simulation(sim)
sim_summary
```

The summary reports:

- Number of views.
- View names.
- Sample size.
- Feature count per view.
- Number of replicates.
- Crosstalk setting.
- Outcome mode.
- Outcome type.
- Signal sharing.
- Differential-effect alignment.

## Core Validation

MimeSys keeps validation focused on two checks: blockwise feature-correlation plots and reconstruction error for observed versus simulated correlation blocks.

For blockwise feature-correlation plots:

```r
observed_corr <- blockwise_feature_correlation(
  views,
  max_features_per_view = 100
)

sim_views <- simulation_view_matrices(
  sim,
  rep = 1,
  split = "train",
  feature_names = "template_feature"
)

sim_corr <- blockwise_feature_correlation(
  sim_views,
  max_features_per_view = 100
)

plot_blockwise_feature_correlation(observed_corr)
plot_blockwise_feature_correlation(sim_corr)
```

For reconstruction error:

```r
feature_correlation_frobenius(
  observed_matrices = views,
  simulated_matrices = sim_views,
  max_features_per_view = 100
)
```

For crosstalk-specific reconstruction-error sweeps, use:

```r
cal <- calibrate_crosstalk(
  template,
  crosstalk_grid = c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1),
  nsample = 200,
  p = 200,
  nrep = 20
)

plot_crosstalk_calibration(cal)
```


## Outcome-Generation Scenarios

A small set of outcome-generation scenarios is usually enough to check that the simulator can create distinct signal regimes.

### Single-Layer Signal

```r
signal.share = c(CYTOF = 1, Metabolomics = 0, Proteomics = 0)
```

Expected behavior:

- This validates that the simulator can create view-specific signal.

### Balanced Additive Signal

```r
signal.share = c(CYTOF = 1, Metabolomics = 1, Proteomics = 1)
ygen.mode = "LM"
```

Expected behavior:

- Fusion can outperform the best single layer when each layer carries useful information.
- The advantage should become clearer as sample size increases.

### Shared-Latent Signal

```r
ygen.mode = "SharedLatent"
```

Expected behavior:

- The outcome is driven by cross-view latent structure.
- The outcome is intrinsically tied to cross-view latent structure.

### Nonlinear Signal

```r
ygen.mode = "Friedman"
ygen.mode = "Friedman2"
```

Expected behavior:

- Linear methods may not fully capture the signal.
- Flexible learners may perform better.
- Fusion behavior depends on how the nonlinear signal is distributed across views.

## Example: Balanced Additive Scenario

```r
sim <- gen_simmba(
  template = template,
  nsample = 200,
  p = c(CYTOF = 200, Metabolomics = 200, Proteomics = 200),
  crosstalk = 0.5,
  nrep = 50,
  ygen.mode = "LM",
  outcome.type = "continuous",
  signal.share = c(CYTOF = 1, Metabolomics = 1, Proteomics = 1),
  de.align = 1,
  snr = 1,
  seed = 1
)
```

## Example: Shared-Latent Scenario

```r
sim <- gen_simmba(
  template = template,
  nsample = 200,
  p = 200,
  crosstalk = 0.75,
  nrep = 50,
  ygen.mode = "SharedLatent",
  outcome.type = "continuous",
  snr = 1,
  seed = 1
)
```

## Example: Nonlinear Friedman Scenario

```r
sim <- gen_simmba(
  template = template,
  nsample = 200,
  p = 200,
  crosstalk = 0.5,
  nrep = 50,
  ygen.mode = "Friedman",
  outcome.type = "continuous",
  signal.share = c(CYTOF = 1, Metabolomics = 1, Proteomics = 1),
  de.align = 1,
  snr = 1,
  seed = 1
)
```

