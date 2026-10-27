# MimeSys 0.99.0

## Survival censoring and view-specific signal (new)

* `cens.type` in `gen_simmba()` and `gen_multistudy()` (`outcome$cens_type`
  in the `simulate_*()` wrappers) chooses the censoring mechanism:
  `"uniform"` (default, unchanged), `"exponential"` (constant loss to
  follow-up) or `"administrative"` (staggered entry over the first
  `cens.accrual` fraction of the study, a fixed study end and optional
  exponential dropout, `cens.dropout`). Each type is scaled to reach
  `cens.rate`.
* `cens.times` takes user censoring times: one time for everyone, one per
  sample, or a function of `n`. With `cens.rate` they are rescaled to the
  target, so they set only the shape of the censoring distribution.
* `cens.dependence` adds informative censoring: censoring times are
  multiplied by `exp(-cens.dependence * z)`, `z` the standardized linear
  predictor, so higher-risk subjects are censored earlier (a scenario that
  breaks the independent-censoring assumption).
* `gen_multistudy()` and `simulate_multistudy()` accept one `cens.rate` per
  study (in study order or named by study), and a list of per-study
  `cens.times`.
* `shared.latent.weight` (`outcome$shared_latent_weight`) for
  `outcome.model = "shared_latent"`: the share of signal variance from the
  shared factors. The rest comes from view-specific latent signal (the part
  of each view's latent factor scores that is not shared with the other
  views, independent across views), so lower values make the views carry
  complementary outcome information. The view-specific part is returned in
  `sample_metadata$Xbeta_specific` and `truth$specific_scores`.
* `baseline.shape` (`outcome$baseline_shape`) gives a Weibull baseline
  hazard, `baseline.hazard * k * t^(k - 1)`: still proportional hazards, with
  a hazard that rises (`k > 1`) or falls (`k < 1`) over follow-up. `snr`
  keeps its log-time meaning.
* `keep.true.times = TRUE` (`outcome$keep_true_times`) stores each sample's
  latent `event_time` and `cens_time` in `sample_metadata`.
* All new arguments default to the previous behaviour; simulations with the
  defaults are bit-identical to before, and new settings are recorded in the
  result only when they differ from their defaults.

## Multi-study simulation (new)

* `gen_multistudy()` simulates several studies from one template, with
  case-mix shifts of the latent factors (`shift.sd`), technical batch effects
  on the latent scale (`batch.loc.sd`, `batch.scale.sd`), study-specific
  signal features (`support.overlap`) and effect sizes (`effect.cor`), and
  study intercepts and SNRs (`intercept.sd`, `snr.sd`). Studies can be held
  out for testing (`holdout.studies`). The predictor is standardized on the
  pooled samples so study differences are kept; binary prevalence is
  calibrated on the pooled sample and survival censoring within each study.
  The full study-level truth is returned. Heterogeneity is drawn from a
  separate random stream, so the base data stay the same across
  heterogeneity levels for a given seed.
* `simulate_multistudy()` is the MAE-first wrapper, with the presets and
  options of `simulate_mimesys()` plus a `heterogeneity` list.
* One template per study: `fit_study_templates()` fits a template to each
  cohort on the views and features every cohort has, and `gen_multistudy()` /
  `simulate_multistudy()` accept that list (or the raw data sets) so each
  study keeps its own cohort's marginals and correlation structure on the
  same simulated features. Crosstalk can be set, or calibrated, per study.
* Bug fixes: a list of matrices without sample IDs (row names) whose views have
  different numbers of samples is now an error; it used to be truncated to
  the smallest view without a warning. `compare_feature_correlation_structure()`
  now counts each pair of distinct features once: self-correlations are left
  out and each cross-view block appears once (it used to keep the 1s on the
  diagonal, which inflated the within-view scores, and list every pair in
  both orientations).
* Internal: `gen_simmba()` now runs on the same outcome engine with one study;
  its results are unchanged (verified bit-for-bit). `sim_from_template()`'s
  internal mapping can map a second set of latent scores (used for batch
  effects).

## Argument names made consistent (breaking)

* The core functions now use one name for one concept, in `dot.case`:
  `nsample` is now `n` in `gen_simmba()`, `gen_mediation()`,
  `calibrate_crosstalk()` and `validate_template()` (as in
  `sim_from_template()` and the wrappers); `tmpl` is now `template` in
  `sim_from_template()` and `crosstalk_max()`; `crosstalk_grid`,
  `max_features_per_view` and `template_args` are now `crosstalk.grid`,
  `max.features.per.view` and `template.args` in `calibrate_crosstalk()` and
  `validate_template()`; `shared.latent.fun` and `shared.latent.n_factors` are
  now `shared.latent.function` and `shared.latent.n.factors` in `gen_simmba()`
  (and in its result). The old names are no longer accepted. The option lists
  of the MAE wrappers (`snake_case`) and the correlation-diagnostic functions
  are unchanged.
* Clearer names for the signal and outcome arguments of `gen_simmba()`:
  `de.prob` -> `signal.prob`, `de.downProb` -> `signal.neg.prob`,
  `de.facLoc` -> `effect.meanlog`, `de.facScale` -> `effect.sdlog`,
  `de.align` -> `signal.align`, `ygen.mode` -> `outcome.model` (values now
  `"linear"`, `"friedman"`, `"friedman2"`, `"shared_latent"`,
  `"interaction"`, as in `simulate_mimesys()`), `p.train` -> `train.fraction`,
  `surv.hscale` -> `baseline.hazard`; the result fields are renamed to match.
  In `simulate_mimesys()`, the `signal` options are now `prob`, `neg_prob`,
  `effect_meanlog`, `effect_sdlog`, `align` and `share`, and the outcome
  option `surv_hscale` is now `baseline_hazard`.
* README rewritten around parameter tables, with the outcome and mediation
  models written out as equations.
* Internal: code shared by `gen_simmba()`, `gen_mediation()` and the MAE
  wrappers now lives in common helpers (argument checks, template
  resolution, the logistic-intercept solver, and the wrappers' set-up,
  provenance and output assembly). Results are unchanged.

## New: sequential mediation

* `gen_mediation()` and `simulate_mediation()` gain `stages` (ordered stages
  of views, e.g. exposure -> microbiome -> metabolome -> outcome),
  `stage.links` (many-to-many links between the mediators of consecutive
  stages: random, hub, convergent or explicit patterns, link strength as a
  share of latent variance, and link signs), `stage.direct` (direct exposure
  paths per stage) and `stage.share` (share of the indirect effect per
  stage). Mediators keep their template distributions. The truth adds
  natural path-specific effects for every path (`path_effects`,
  `NIE_by_stage`, path-specific risk differences for a binary outcome), the
  list of links, and per-feature `stage`, `n_parents` and `n_children`.
  Without `stages`, results are unchanged.
* Internally, the latent scores of all views are now built before any is
  mapped to the data scale, so a latent effect can link views; the random
  stream is unchanged.

## Review fixes

* Defaults are consistent: every simulator returns the fitted scale by default
  (`back_transform = FALSE` in `simulate_mimesys()` and
  `simulate_mimesys_null()`), and `fit_template()`, `sim_from_template()`
  (new `seed` argument), `gen_simmba()`, `gen_mediation()`,
  `calibrate_crosstalk()`, `validate_template()` and the MAE wrappers all use
  `seed = 1234`. `gen_simmba()` now defaults to `nrep = 1`.
* `simulate_mimesys_null()` gains `nrep` and `return_type`.
* MultiAssayExperiment input is matched by the primary sample IDs in
  `sampleMap()`, so experiments with different column names for the same
  sample are aligned correctly.
* Input checks: integer `n`, `nrep`, `nperm`, `n.interactions` and
  `shared.latent.n_factors`; `feature.jitter` in \[0, 1\]; `var.explained`
  in (0, 1\]; unknown view names in `p` are an error.
* Crosstalk calibration needs at least two features per view; with fewer,
  `calibrate_crosstalk()` stops with a clear message and the MAE wrappers skip
  calibration. Calibration is faster: observed correlation blocks are
  computed once, only the needed blocks are computed, and the complete-data
  correlation path is used when there are no missing values.
* Internal callers no longer reset the seed on every replicate (each
  calibration and null replicate is now an independent draw).
* In the MAE wrappers, crosstalk calibration uses `seed + 1`, so its
  replicates do not share random numbers with the simulation, and
  `simulate_mimesys(seed = s)` matches `gen_simmba(seed = s)` at the selected
  crosstalk. `simulate_mediation()` records its seed in the metadata.
* MAE experiments stored as plain matrices are transposed like assays.
* The mediation Monte Carlo runs in chunks of 1000 draws and keeps only
  per-feature means and mediator columns, cutting memory use; seeded mediation
  results differ from earlier development versions.
* `gen_mediation()` with `TIE = 0` reports `TIE_target = 0`.
* `two_part_copula` is for sparse non-negative (microbiome) data: views with
  negative values or fitted on the CLR scale are an error, and views without
  zeros give a warning. The abundance latent variable is rescaled to unit
  variance, so its uniform is exact. Feature metadata `scale` now reports the
  scale actually returned for each view.
* `validate_template()` chooses the top-variable features on the observed data
  and uses the same features for the simulated data.
* Documentation: linear-predictor standardization (`gen_simmba()` Details),
  `max.rank` capping an explicit `rank`, README examples use the
  `stelzer_mae` view names and runnable presets.

## New: mediation simulation

* Review fixes for mediation: a binary exposure now keeps mediator latent
  scores exactly standard normal (mixture-CDF correction), so mediators keep
  their template distributions; `simulate_mediation()` now returns the fitted
  scale by default, matching the scale of the truth (`feature_scale`); for a
  binary outcome `NIE`/`NDE`/`TE` are labelled conditional log-odds effects
  and natural log odds ratios are added to `truth$mediation$binary`; new
  `exposure_effect` column gives the true exposure effect on every feature
  (non-zero for non-mediators only through compositional closure); interaction
  coefficients no longer cancel with `dir.b = "random"`; `null_a` b-paths are
  calibrated as in the spiked scenario; duplicate `crosstalk` in MAE metadata
  removed; `n` validated; paired Monte Carlo stays in step for empty samples
  in closed compositions. This last change alters the random stream for
  closed compositional views, so seeded simulations of such views do not
  reproduce results from earlier development versions.

* `gen_mediation()` and `simulate_mediation()` simulate an exposure (binary or
  continuous), mediators in any views (`n.mediators` per view), and an outcome
  with direct effect, optional mediator interactions and null scenarios
  (`null_a`, `null_b`, `null_ab`). Exposure effects are injected on the latent
  normal-score scale, so mediators keep their template distributions. The
  b-paths are calibrated so the population natural indirect effect equals
  `TIE`, and the true NIE, NDE, TE, per-view indirect effects and
  per-mediator paths are returned. Continuous or binary outcome (binary: TIE
  and DE on the log-odds scale, intercept solved for `outcome.prevalence`,
  risk-difference and odds-ratio effects reported). For a continuous exposure
  the truth also includes `NIE_linear`, the linear-projection estimand.
  Domain-agnostic.

## Review fixes

* Compositions: simulated samples of a compositional view never exceed the
  template total. Closed compositions are renormalized to their total when at
  least the template's features are simulated (offspring features now dilute
  the others instead of pushing totals above 1), and capped at the total when
  fewer are simulated; for a template that is a subset of a composition,
  samples are capped at the largest observed total. Applies to untransformed
  views and to back-transformed `clr`/`log1p` views.
* Every view with `de.prob > 0` now gets at least one spiked feature, so a
  replicate can no longer have zero outcome signal by chance. Results for
  small `p` can change.

* `calibrate_crosstalk()` chose the grid value that best fit a single
  cross-view block. It now averages each replicate's distance over all
  cross-view blocks and picks the smallest mean; `cross_view_summary` holds
  the mean, SD and SE per grid value and `near_optimal` is based on it. Earlier
  calibrated crosstalk values (including the `crosstalk = "calibrate"` default)
  can change.
* Seeds: `fit_template()` gains `seed`; `simulate_mimesys()`,
  `simulate_mimesys_null()` and `validate_template()` now seed template
  fitting and simulation, so results are reproducible and the caller's random
  state is restored.
* Shared rank: when `shared.rank = NULL`, the rank is chosen by sequential
  (deflated) parallel analysis: the largest GCCA eigenvalue is tested against
  the permutation null of the largest eigenvalue, the accepted direction is
  removed from every view, and the test is repeated. The previous
  index-matched test over-selected (2, 3 and 5 true shared factors were
  estimated as 6-15); the sequential test recovers them in most runs. On the
  examples (default seed), the Franzosa rank went from 8 to 3 (2 for some
  seeds) and the full Stelzer data at one timepoint from 5 to 2. Per-step p-values, null thresholds and tested eigenvalues are in
  `gcca.pvalue`, `gcca.threshold` and `gcca.step.eigen`; a warning flags
  borderline choices. Optional stability selection (`shared.rank.runs > 1`,
  default 1) repeats the choice and uses the rounded mean (stored in
  `shared.rank.selection`). `nperm` default raised from 50 to 200.
  `shared.rank` cannot exceed the smallest within-view rank (an explicit
  value above it is an error; an automatic choice above it is capped with a
  message).
* Per-view arguments (`crosstalk`, `rank`, `de.prob`, `de.downProb`,
  `de.facLoc`, `de.facScale`, `signal.share`, `friedman.share`) are matched by
  view name; wrong lengths and unknown names are errors. Unknown names in
  `view.transform` are errors.
* `friedman.share` now works at the view level: each Friedman feature comes
  from a view chosen with probability equal to its share.
* Truth flags: interaction or Friedman features are not marked when their
  weight is 0.
* Single-view templates and templates without shared factors skip
  calibration in `simulate_mimesys()`, `simulate_mimesys_null()` and
  `validate_template()`.
* `summarize_simulation()` reports prevalence (binary), event rate and a
  time-signal correlation (survival), and counts `outcome_feature`.
* Input checks: duplicate sample IDs, all-missing features (dropped with a
  message), features more than 50% missing (warning), views with no
  non-constant features, and out-of-range `de.prob`, `de.align` and `nrep`.
* Views keep the template order in MAE output and `simulation_view_matrices()`.
* Feature standardization is numerically stable for features with large means.
* LICENSE copyright holder is now "MimeSys authors"; `mediation_spiking/` is
  excluded from the package build.

## Bioconductor preparation

* Version 0.99.0 with `Authors@R` (Nalin Arora maintainer; Himel Mallick and
  Ziyu Liu authors), `URL`, `BugReports`, valid `biocViews`, and `stats`,
  `utils` and `withr` declared in `Imports`.
* New example datasets `stelzer_mae` and `franzosa_mae` (small
  `MultiAssayExperiment` objects; built by `inst/scripts/make-example-data.R`).
  The large raw files in `inst/extdata` are excluded from the package build.
* Every exported function has a runnable example, and the vignette is
  rewritten as a fully evaluated `BiocStyle` workflow on the example data.
* Seeds are set with `withr::local_seed()`, so `gen_simmba()`,
  `calibrate_crosstalk()` and `validate_template()` no longer change the
  global random-number state; `seed = NULL` uses the current state.
* Code uses 4-space indentation, `seq_len()` and vectorized column helpers
  (sorting, variances, integer checks, constant-feature filtering), with
  output identical to the previous version.
* `calibrate_crosstalk()` warns when the template has no shared factors.

# MimeSys 0.0.0.9000

## Breaking changes

* Removed the `"empirical_copula"` and `"module_copula"` dependence modes and
  their arguments (`dependence`, `cor.shrinkage`, `module.count`,
  `module.between.shrinkage`; `null_model$dependence`, `cor_shrinkage`,
  `module_count`, `module_between_shrinkage`). On held-out data they did not
  reproduce correlation structure better than the factor model, and several
  features (offspring jitter, shared-latent outcomes, `de.align`) did not
  apply to them.
* `gen_simmba()` no longer returns the fitted template by default, because the
  template stores the source data. Use `return.template = TRUE` to keep it;
  `template_provenance` is always returned.

## Bug fixes

* Outcomes are now generated from standardized features. Previously the linear
  predictor used raw feature values, so views on large measurement scales
  dominated the signal, binary outcomes were near-deterministic with
  uncontrolled prevalence, and survival outcomes could be `NaN`.
* Binary outcomes: new `binary.prevalence` (default 0.5) sets the case
  proportion via a solved intercept; `snr` is now the latent-logistic
  signal-to-noise ratio.
* Survival outcomes: `snr` is now the log-time signal-to-noise ratio; new
  `cens.rate` targets the censoring proportion.
* `two_part_copula` microbiome view: positive abundances are now drawn
  conditional on presence (previously inflated about 3-fold).
* `simulate_mimesys()` accepts `outcome = list(prevalence = , cens_rate = )`.

## New features

* New `ygen.mode = "Interaction"` (preset `"interaction"`): the outcome mixes
  the linear predictor with products of feature pairs from different views,
  controlled by `interaction.weight` and `n.interactions`. Pairs are stored in
  `truth$interaction_pairs`.
* New `nonlinear.weight` (default 0.5) sets the nonlinear share for
  `"Friedman2"` and `shared.latent.fun = "friedman2"`. Previously `Friedman2`
  added the raw parts, so the nonlinear part's share depended on scale.
* `back.transform` in `sim_from_template()` and `gen_simmba()` returns views
  fitted with `view.transform` on the original measurement scale (exact for
  `log1p`; for `clr`, each feature's original distribution including zeros).
  Every simulator returns the fitted (transformed) scale by default; use
  `back.transform = TRUE` (or `null_model = list(back_transform = TRUE)`) for
  the original scale.
* `feature_metadata$outcome_feature` gives the features that drive the
  outcome in every mode; `beta` is now zero where the outcome does not use the
  linear predictor (`Friedman`, `SharedLatent`, pure `Interaction`).
* `validate_template()` marginal diagnostics add scale-free `mean_ks`,
  `median_ks` and `mean_std_quantile_rmse`; the marginal plot uses `mean_ks`.

## More bug fixes

* `calibrate_crosstalk()` with raw input and a `view.transform` compared
  transformed simulations with untransformed data; it now always uses the
  template's stored matrices.
* `calibrate_crosstalk()` drops grid values above `crosstalk_max()` with a
  warning instead of failing, and warns when the optimum is at the top of a
  grid that stops well short of the limit.
* Friedman feature sampling could pick features outside the allowed set when
  only one fallback feature was available (`sample()` on a length-one vector).
* `gen_simmba()` validates `p.train` (no empty train/test sets) and warns when
  a `SharedLatent` outcome cannot be related to the features (crosstalk 0 or
  no shared factors).

## Documentation

* The template assumes one sample per subject; README and vignette now say to
  fit repeated-measures data on a single timepoint.

## Performance

* `blockwise_feature_correlation()` is vectorized (about 100x faster for a few
  hundred features), which also speeds up `validate_template()`.
* Fitted templates no longer store per-sample factor scores or latent normal
  scores.
