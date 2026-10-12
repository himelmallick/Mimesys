# MimeSys

MimeSys is a template-based simulator for multi-omics method benchmarking. It learns layer-specific distributions and cross-layer dependence from an empirical dataset, then generates synthetic multi-view data with controlled sample size, feature count, cross-omics correlation, differential-effect alignment, and outcome-generation structure.

The current package is designed around two linked goals:

1. Generate realistic null multi-omics data that preserves important empirical features of the source dataset.
2. Spike outcomes in controlled ways so benchmarking methods can be tested under scenarios where early, intermediate, late, or single-view methods are expected to behave differently.

## Installation

```r
# Install from GitHub
install.packages("pak")
pak::pak("himelmallick/MimeSys")

# Load the package
library(MimeSys)
```

For local development from the package root:

```r
pak::local_install(".")
# or, during active development
devtools::load_all(".")
```

Core packages used by MimeSys include `Matrix`, `MASS`, `mvtnorm`, `ggplot2`, and `dplyr`. Some benchmark scripts use additional modeling packages such as `SuperLearner`, `glmnet`, `xgboost`, or `multiview`, depending on the experiment being run.

## Basic Workflow

The typical workflow has three steps.

```r
library(MimeSys)

# 1. Split an empirical multi-omics dataset into named omics views
views <- list(
  Immunome   = immunome_matrix,
  Metabolome = metabolome_matrix,
  Proteome   = proteome_matrix
)

# 2. Fit a simulation template
template <- fit_template(
  views,
  rank = 10,
  view.transform = "none"
)

# 3. Generate benchmark-ready simulated data from the fitted template
sim <- gen_simmba(
  template = template,
  nsample = 200,
  p = c(Immunome = 200, Metabolome = 200, Proteome = 200),
  crosstalk = 0.5,
  nrep = 20,
  ygen.mode = "LM",
  outcome.type = "continuous",
  signal.share = c(Immunome = 1, Metabolome = 1, Proteome = 1),
  de.align = 1,
  snr = 1
)

summarize_simulation(sim)
```

The returned object contains simulated training and test data for each replicate, along with metadata describing the simulation settings.

## Input Data

Input data should be supplied as a named list of matrices or matrix-like objects.

```r
views <- list(
  Immunome   = matrix_1,
  Metabolome = matrix_2,
  Proteome   = matrix_3
)
```

Each matrix should have samples in rows and features in columns. All views should contain the same samples in the same order unless you have already handled matching externally.

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
  nperm = 50,
  view.transform = NULL
)
```

### `template`

A named list of omics matrices. Rows are samples and columns are features. The argument is called `template` because these empirical matrices are used as the template for null-data generation.

### `rank`

The number of latent factors to estimate from the empirical data.

If `rank` is supplied, the template uses exactly that number of factors. If `rank = NULL`, MimeSys chooses the number of factors using `var.explained`.

### `var.explained`

The target cumulative variance explained when `rank = NULL`.

For example, `var.explained = 0.7` asks the template to retain enough latent factors to explain approximately 70% of the variation in the concatenated transformed data.

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

### `max.rank`

Upper bound on the number of factors retained per view when `rank = NULL`.

If `NULL`, the default cap is approximately one third of the empirical sample size.

### `shared.rank`

Number of shared latent variates to retain across views.

If `NULL`, MimeSys chooses the shared rank using a permutation-based parallel-analysis procedure. This helps avoid treating finite-sample noise as real cross-view structure.

### `nperm`

Number of permutations used when `shared.rank = NULL`.

Larger values give a more stable null threshold for shared-rank selection but take longer.

## What The Template Stores

The fitted template stores:

- View names.
- Feature names.
- Transformation choices.
- Empirical means and scales.
- Latent factor structure.
- View-specific loading matrices.
- Residual variation.
- Empirical dependence mode information used for null simulation.

The template is intentionally reusable. Fit it once from the empirical dataset, then generate many synthetic datasets across a grid of sample sizes, feature counts, crosstalk values, and outcome-generation scenarios.

## Generating Benchmark Datasets

`gen_simmba()` is the high-level function used for benchmarking. It calls the null simulator, spikes feature effects, generates outcomes, and splits each replicate into training and test data.

```r
sim <- gen_simmba(
  nsample = 200,
  template = template,
  p = 200,
  crosstalk = 0.5,
  feature.select = "random",
  feature.jitter = 0.2,
  dependence = "factor",
  cor.shrinkage = 0.05,
  module.count = NULL,
  module.between.shrinkage = 0.75,
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
  surv.hscale = 1,
  cens.lower = 1,
  cens.upper = 3,
  nrep = 10,
  seed = 1234
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
  dependence = "factor",
  cor.shrinkage = 0.05,
  module.count = NULL,
  module.between.shrinkage = 0.75,
  view.models = NULL
)
```

The sections below describe the high-level `gen_simmba()` arguments. Arguments that are passed directly to `sim_from_template()` are noted.

### `nsample`

Number of samples to simulate in each replicate before train/test splitting.

This can be smaller, equal to, or larger than the empirical sample size. In benchmarking, this is usually varied to test whether method behavior stabilizes as sample size increases.

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
p = c(Immunome = 100, Metabolome = 200, Proteome = 500)

# Named list
p = list(Immunome = 100, Metabolome = 200, Proteome = 500)
```

If `p = NULL`, MimeSys uses the original number of features available in each layer of the template.

### `crosstalk`

Controls the strength of cross-omics dependence during null data generation.

Conceptually:

- Low `crosstalk` means the views are more independent.
- High `crosstalk` means the views share more latent structure and are more correlated across layers.

This is a simulation-control parameter, not an automatically estimated biological truth. It changes how strongly latent factors are shared across layers. Empirically, higher crosstalk should generally increase cross-view correlation, but it does not necessarily imply that fusion methods will perform better. High crosstalk can create redundancy rather than complementarity.

`crosstalk` can be a scalar shared by all views or a vector with one value per view. Values must be non-negative and cannot exceed `crosstalk_max(template)`.

### `feature.select`

Feature-selection rule used when the requested feature count is smaller than the number available in the template.

Supported values:

- `"random"`: randomly sample template features.
- `"variance"`: choose the highest-variance empirical features first.

This argument is passed to `sim_from_template()`.

### `feature.jitter`

Amount of random perturbation applied when generating offspring features beyond the empirical template feature count.

Use `0` for closer copies of parent features and larger values for more feature-level variation. Typical values are between `0` and `1`.

This argument is passed to `sim_from_template()`.

### `dependence`

Controls the null dependence simulator.

Supported values:

- `"factor"`: use the low-rank template factor model.
- `"empirical_copula"`: use the empirical latent feature correlation more directly.
- `"module_copula"`: smooth empirical latent correlations through within-view modules.

This argument is passed to `sim_from_template()`.

### `cor.shrinkage`

Shrinkage toward independence for `"empirical_copula"` and `"module_copula"`.

Higher values make the latent correlation matrix more conservative.

### `module.count`

Number of within-view modules used when `dependence = "module_copula"`.

If `NULL`, MimeSys uses a data-adaptive default per view.

### `module.between.shrinkage`

How strongly to smooth correlations between modules when `dependence = "module_copula"`.

Larger values shrink between-module correlations more strongly.

### `view.models`

Optional named vector controlling the observation model for each view.

Supported values:

- `"continuous"`: standard empirical-quantile emission.
- `"two_part_copula"`: two-stage sparse microbiome-like emission.

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

Can be a scalar or one value per view.

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

If `NULL`, spike-in coefficients are left on their original scale.

If supplied, values are recycled or matched by view name, normalized to sum to one, and used to rescale view-specific coefficients.

### `friedman.share`

Optional numeric vector controlling which views supply the five nonlinear Friedman features when `ygen.mode` is `"Friedman"` or `"Friedman2"`.

If `NULL`, MimeSys uses `signal.share` when available; otherwise Friedman features are sampled uniformly from all nonzero spiked features.

### `ygen.mode`

Outcome-generation mode.

Supported values:

- `"LM"`
- `"Friedman"`
- `"Friedman2"`
- `"SharedLatent"`

These are described in detail below.

### `outcome.type`

Observed outcome type.

Supported values:

- `"continuous"`
- `"binary"`
- `"survival"`

### `surv.hscale`

Multiplicative hazard scale for survival outcomes.

### `cens.lower` and `cens.upper`

Lower and upper bounds for uniform censoring times in survival simulations.

### `nrep`

Number of independent simulation replicates.

Each replicate contains a training dataset and a test dataset generated under the same scenario.

### `seed`

Random seed for reproducible simulation.

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
- `feature_metadata`: feature-level metadata including feature type, parent template feature, offspring status, view model, and shared loading.
- `truth`: latent scores, shared variates, crosstalk weights, canonical correlations, dependence mode, and view models.

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

In many benchmark scripts we use the interpretable endpoints:

- `0`: random effect orientation.
- `1`: aligned effect orientation.

## Signal Sharing

`signal.share` controls how much of the outcome signal is assigned to each omics layer.

```r
signal.share = c(Immunome = 1, Metabolome = 1, Proteome = 1)
```

The values are normalized internally, so the above means that the outcome signal is balanced across the three views.

Single-layer scenarios can be created by assigning all signal to one view:

```r
signal.share = c(Immunome = 1, Metabolome = 0, Proteome = 0)
```

This is useful for validating that single-layer methods win when only one layer carries the outcome signal.

Balanced scenarios are useful for testing whether fusion methods can benefit from distributed signal:

```r
signal.share = c(Immunome = 1, Metabolome = 1, Proteome = 1)
```

Important distinction:

- `crosstalk` controls feature-level cross-view dependence.
- `signal.share` controls outcome-signal allocation across views.
- `de.align` controls whether the spiked outcome signal is aligned with shared latent structure.

These are separate simulation knobs.

## Outcome Generation

MimeSys currently supports four outcome-generation modes through `ygen.mode`.

```r
ygen.mode = "LM"
ygen.mode = "Friedman"
ygen.mode = "Friedman2"
ygen.mode = "SharedLatent"
```

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

This mode retains the additive layer signal and adds a Friedman-style nonlinear component. It is useful for testing methods under partially linear and partially nonlinear outcome structure.

### `SharedLatent`

Outcome generation from the shared latent structure.

MimeSys uses the shared latent factors generated during null data simulation. Instead of using only one factor, the current implementation combines all available shared factors by standardizing them and summing across factors:

```text
eta_shared = rowSums(scale(shared_factors)) / sqrt(number_of_factors)
```

The outcome is then generated from this shared latent score with noise controlled by `snr`.

This mode is useful when the true outcome is driven by a cross-view latent process rather than view-specific feature effects.

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

The latent score is converted to event probabilities and sampled as a binary outcome.

This is useful for classification benchmarks.

### Survival

The latent score is converted into survival-like outcome metadata. The generated sample metadata contains survival fields such as time and event status.

This is useful for survival-model benchmarking.

## Signal-To-Noise Ratio

`snr` controls outcome signal strength relative to noise.

Higher `snr` creates easier prediction tasks. Lower `snr` creates harder tasks.

For benchmarking, it is often useful to keep `snr` fixed while varying sample size, feature count, crosstalk, and signal-sharing scenario.

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

### Why This Matters

This model better reflects microbiome-like data than a single Gaussian emission because it preserves:

- Sparsity.
- Zero inflation.
- Positive abundance structure.
- Dependence between presence and abundance.
- Cross-view dependence through the shared latent simulator.

The two-stage emission is part of null feature generation. It is separate from outcome generation, which is controlled by `ygen.mode`.

## Additive Versus Shared-Latent Outcome Generation

Two outcome-generation modes are especially important for current benchmarking.

### Additive Outcome Generation

Additive outcomes use the layer-specific feature matrices directly:

```text
eta = w_1 X_1 beta_1 + w_2 X_2 beta_2 + ... + w_M X_M beta_M
```

This supports scenarios where each layer contributes its own piece of the outcome.

Expected behavior:

- Single-layer signal: best single layer should perform well.
- Balanced signal across layers: fusion may outperform the best single layer.
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
- Intermediate or late fusion can become more competitive when the shared latent process is recoverable.
- Estimated cooperative `rho` does not have to increase monotonically with crosstalk. Crosstalk controls cross-view correlation in the simulated features, while `rho` is selected by predictive cross-validation.

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
train_views$Immunome
train_views$Metabolome
train_views$Proteome
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

## Correlation Diagnostics

MimeSys includes tools for validating whether the simulated data preserve or manipulate correlation structure as expected.

```r
validate_template_simulation(sim)
```

This computes diagnostics comparing correlation behavior across simulated views.

For crosstalk-specific diagnostics, use:

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

In validation experiments, useful expected patterns are:

- Reconstruction error should decrease as sample size increases.
- Reconstruction error should remain stable as the number of features increases while sample size is fixed.
- Cross-view correlation should generally increase with crosstalk.
- Predictive fusion advantage may be strongest when signal is shared across layers and not purely redundant.

## Benchmarking Scenarios

A small set of benchmark scenarios is usually enough to validate the simulator.

### Single-Layer Signal

```r
signal.share = c(Immunome = 1, Metabolome = 0, Proteome = 0)
```

Expected behavior:

- The informative single layer should perform best.
- Fusion should not be forced to win.
- This validates that the simulator can create view-specific signal.

### Balanced Additive Signal

```r
signal.share = c(Immunome = 1, Metabolome = 1, Proteome = 1)
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
- Intermediate and late fusion become natural competitors to early fusion.

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
  p = c(Immunome = 200, Metabolome = 200, Proteome = 200),
  crosstalk = 0.5,
  nrep = 50,
  ygen.mode = "LM",
  outcome.type = "continuous",
  signal.share = c(Immunome = 1, Metabolome = 1, Proteome = 1),
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
  signal.share = c(Immunome = 1, Metabolome = 1, Proteome = 1),
  de.align = 1,
  snr = 1,
  seed = 1
)
```

## Recommended Reporting Language

When describing a MimeSys simulation, report:

- Empirical dataset used to fit the template.
- Views included.
- Sample sizes simulated.
- Number of features per layer.
- Crosstalk grid.
- Outcome-generation mode.
- Outcome type.
- Signal-sharing vector.
- Differential-effect alignment.
- Signal-to-noise ratio.
- Number of replicates.
- Learners and fusion methods used for benchmarking.

Example:

```text
We fit a MimeSys template using the Stelzer multi-omics dataset with Immunome,
Metabolome, and Proteome views. Synthetic datasets were generated with 100, 200,
and 500 features per layer; sample sizes of 53, 100, 200, and 500; and crosstalk
values from 0 to 1. Outcomes were generated either additively across layers or
from the shared latent structure. For additive outcomes, signal sharing was
balanced across the three layers unless otherwise noted, and differential effects
were either aligned with the shared latent structure or randomly oriented.
```

## Local Checks

Run tests with:

```r
devtools::test()
```

Run a package check with:

```sh
R CMD build .
R CMD check --no-manual --no-build-vignettes MimeSys_*.tar.gz
```

## Notes On Interpretation

Crosstalk is not the same as fusion advantage. Crosstalk controls the cross-view correlation structure in the simulated features. Fusion advantage depends on how outcome-relevant signal is distributed across layers and whether the information is complementary or redundant.
