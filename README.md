# MimeSys

MimeSys is a template-based simulator for multi-omics data. It learns each
feature's distribution, the dependence within each omics view and the
dependence between views from a real data set (the *template*), then generates
synthetic multi-view data with any sample size and feature count, a single dial
for cross-view dependence (*crosstalk*), and planted signal with known truth:
outcomes (linear, nonlinear, shared-latent or interaction; continuous, binary
or survival) or exposure-mediator-outcome structures (parallel or sequential
mediation). It is domain-agnostic and works with `MultiAssayExperiment`
objects.

- Full walkthrough: `vignette("mimesys-workflow", package = "MimeSys")`
- Algorithms in detail: [`docs/MimeSys_algorithms.md`](docs/MimeSys_algorithms.md)
  (PDF: [`docs/MimeSys_algorithms.pdf`](docs/MimeSys_algorithms.pdf))
- Example data: `?stelzer_mae` (CyTOF, metabolome, proteome; 53 subjects) and
  `?franzosa_mae` (gut microbiome and metabolome; 155 samples)

## Installation

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("himelmallick/Mimesys")
library(MimeSys)
```

For development from the package root: `devtools::load_all(".")`.
MimeSys imports `MultiAssayExperiment`, `SummarizedExperiment`, `S4Vectors`
and `withr`; plotting helpers use `ggplot2` when installed.

## Quick start

```r
library(MimeSys)
library(MultiAssayExperiment)
data(stelzer_mae)

# 1. Check the input, then fit the template once and reuse it
validate_mimesys_input(stelzer_mae)
template <- fit_template(
  stelzer_mae,
  view.transform = c(Metabolomics = "log1p", Proteomics = "log1p")
)

# 2. simulate_mimesys() runs every simulation; `design` chooses what:
#    "outcome" (default), "null", "multistudy" or "mediation".
#    Simulate an outcome study. Crosstalk is calibrated to the template by
#    default; views are returned on the fitted (here log1p) scale.
sim_mae <- simulate_mimesys(
  template,
  n = 200,
  p = c(CYTOF = 30, Metabolomics = 30, Proteomics = 30),
  outcome = "linear",
  signal = list(share = c(1, 1, 1))    # equal signal in each view
)
colData(sim_mae)[, c("Y", "split")]    # outcome and train/test split
rowData(experiments(sim_mae)[["Proteomics"]])   # per-feature truth (beta, ...)
metadata(sim_mae)$simulation           # settings, calibration, provenance

# 3. Null data only (no outcome)
null_mae <- simulate_mimesys(template, design = "null", n = 200, p = 200)

# 4. Mediation (see "Mediation" below)
med_mae <- simulate_mimesys(
  template, design = "mediation", n = 300,
  mediation = list(n_mediators = c(CYTOF = 3, Metabolomics = 3, Proteomics = 0),
                   tie = 1, de = 0.5)
)
metadata(med_mae)$truth$mediation[c("NIE", "NDE", "TE")]

# 5. Several studies with between-study heterogeneity (see "Multi-study
#    simulation" below); study 3 is held out for testing
ms_mae <- simulate_mimesys(
  template, design = "multistudy", n = c(150, 200, 100), p = 30,
  heterogeneity = list(effect_cor = 0.5, support_overlap = 0.7,
                       shift_sd = 0.3, batch_loc_sd = 0.3,
                       holdout_studies = 3)
)
table(colData(ms_mae)$study, colData(ms_mae)$split)
metadata(ms_mae)$truth$studies         # study sizes, intercepts, SNRs
```

Every function that uses random numbers has `seed = 1234` by default and leaves
your own random-number state unchanged; use `seed = NULL` to draw from the
current stream.

## Functions

| Function | What it does |
|---|---|
| `validate_mimesys_input()` | Samples, features, missing values and constant features per view, before fitting |
| `fit_template()` | Learns the template: marginals, per-view factor models, shared cross-view space |
| `sim_from_template()` | Null multi-view data (no outcome) |
| `gen_simmba()` | Null data plus planted outcome signal and a train/test split |
| `gen_multistudy()` | Several studies with case-mix, batch, coefficient, intercept and SNR heterogeneity, from one template or one template per study |
| `fit_study_templates()` | One template per study (real cohort) on the views and features every study has |
| `gen_mediation()` | Exposure, mediators and outcome with exact natural (path-specific) effects |
| `calibrate_crosstalk()` | Chooses crosstalk so simulated cross-view correlations match the template |
| `validate_template()` | Calibration plus correlation and marginal diagnostics |
| `crosstalk_max()` | Largest valid crosstalk per view |
| `simulate_mimesys()` | **Main interface.** Fits if needed, calibrates crosstalk, simulates and returns `MultiAssayExperiment` objects; `design = "outcome"` (default), `"null"`, `"multistudy"` or `"mediation"` |
| `mimesys_presets()` | Outcome presets for `simulate_mimesys()` |
| `as_mimesys_mae()`, `simulation_view_matrices()`, `as_mimesys_view_list()` | Conversions between simulation results, MAEs and lists of matrices |
| `blockwise_feature_correlation()`, `feature_correlation_frobenius()`, `compare_feature_correlation_structure()`, `summarize_simulation()`, `plot_*()` | Diagnostics and plots |

## Parameter reference

"Per view" arguments accept one value for every view, one value per view in
template order, or a vector named by view.

### `fit_template()`

| Argument | Default | Description |
|---|---|---|
| `template` | required | `MultiAssayExperiment`, `SummarizedExperiment`, or named list of matrices (samples in rows). Samples are matched across views by ID (for an MAE, the `sampleMap()` primary IDs); at least 10 shared samples. Fit repeated-measures data on one timepoint. |
| `rank` | `NULL` | Factors per view (per view). `NULL`: smallest rank explaining `var.explained`; `"2logp"`: `ceiling(2 * log(p))`. |
| `var.explained` | `0.7` | Share of normal-score variance the factors must explain when `rank = NULL`. |
| `max.rank` | `floor(n0 / 3)` | Cap on the per-view rank (also applies to an explicit `rank`). |
| `shared.rank` | `NULL` | Number of shared cross-view factors. `NULL`: chosen by sequential parallel analysis. |
| `nperm` | `200` | Permutations per parallel-analysis step. |
| `view.transform` | `NULL` | Per view: `"none"`, `"log1p"` or `"clr"`, applied before fitting. |
| `shared.rank.runs` | `1` | Stability selection: repeat the shared-rank choice this many times and use the rounded mean. |
| `seed` | `1234` | Random seed. |

### `sim_from_template()`

| Argument | Default | Description |
|---|---|---|
| `template` | required | Fitted template. |
| `n` | required | Number of samples to simulate. |
| `p` | `NULL` | Features per view (per view). `NULL`: template sizes. Fewer subsamples features; more adds "offspring" features. |
| `crosstalk` | `1` | Cross-view dependence (per view): `0` independent views, `1` template level, up to `crosstalk_max()`. |
| `feature.select` | `"random"` | When `p` is below the template size: `"random"` or `"variance"` (most variable first). |
| `feature.jitter` | `0.2` | In [0, 1]: how far offspring features' loadings move away from their parent's. |
| `view.models` | `NULL` | Per view: `"continuous"` (default) or `"two_part_copula"` (sparse non-negative microbiome views: presence and abundance as two linked parts). |
| `back.transform` | `FALSE` | Return transformed views on the original measurement scale. |
| `seed` | `1234` | Random seed. |

### `gen_simmba()`

Takes `n`, `template`, `p`, `crosstalk`, `feature.select`, `feature.jitter`,
`view.models`, `back.transform` and `seed` as above, plus:

| Argument | Default | Description |
|---|---|---|
| `signal.prob` | `0.1` | Per view: probability that a feature carries signal (nonzero coefficient); at least one signal feature per view with `signal.prob > 0`. |
| `signal.neg.prob` | `0.5` | Per view: probability that a signal coefficient is negative. |
| `effect.meanlog`, `effect.sdlog` | `1`, `0.4` | Per view: mean and SD on the log scale of the log-normal effect-size factor $f$; the coefficient size is $\lvert\log_2 f\rvert$. |
| `signal.align` | `NULL` | In [0, 1]: share of signal features placed on, and signed by, the shared cross-view direction. `NULL`: random positions and signs. |
| `signal.share` | `NULL` | Per view: relative variance of each view's linear predictor (e.g. `c(1, 1, 1)` balanced, `c(1, 0, 0)` one view only). |
| `friedman.share` | `NULL` | Per view: which views supply the five Friedman features (defaults to `signal.share`). |
| `outcome.model` | `"linear"` | `"linear"`, `"friedman"`, `"friedman2"`, `"shared_latent"` or `"interaction"` (see equations below). |
| `shared.latent.function` | `"linear"` | For `"shared_latent"`: `"linear"`, `"friedman"` or `"friedman2"` (the Friedman forms need at least 5 shared factors). |
| `shared.latent.n.factors` | `NULL` | For `"shared_latent"`: number of shared factors used (`NULL`: all). |
| `shared.latent.weight` | `1` | For `"shared_latent"`: share of signal variance from the shared factors; the rest is view-specific latent signal, so lower values make the views complementary. |
| `nonlinear.weight` | `0.5` | Share of signal variance from the nonlinear part (`"friedman2"`, `shared.latent.function = "friedman2"`). |
| `interaction.weight` | `0.5` | Share of signal variance from cross-view interactions (`"interaction"`; 1 = no main effects). |
| `n.interactions` | `NULL` | Number of cross-view feature pairs for `"interaction"` (default: one per view pair). |
| `outcome.type` | `"continuous"` | `"continuous"`, `"binary"` or `"survival"`. |
| `snr` | `1` | Signal-to-noise ratio (on the latent logistic or log-time scale for binary and survival outcomes). |
| `binary.prevalence` | `0.5` | Target case proportion for binary outcomes. |
| `baseline.hazard` | `1` | Baseline hazard $h_0$ for survival outcomes. |
| `baseline.shape` | `1` | Weibull shape $k$: hazard $h_0 k t^{k-1} e^{\eta}$ (1 = exponential; $k > 1$ rising, $k < 1$ falling hazard). |
| `keep.true.times` | `FALSE` | Also store the latent `event_time` and `cens_time` of every sample. |
| `cens.rate` | `NULL` | Target censoring proportion; `NULL` uses uniform censoring on (`cens.lower`, `cens.upper`) = (1, 3). |
| `cens.type` | `"uniform"` | Censoring mechanism: `"uniform"`, `"exponential"` (constant loss to follow-up) or `"administrative"` (staggered entry over the first `cens.accrual` = 0.5 of the study, fixed study end, dropout probability `cens.dropout` = 0); scaled to reach `cens.rate`. |
| `cens.times` | `NULL` | User censoring times: one time, one per sample, or a function of `n`; rescaled to `cens.rate` when it is given. |
| `cens.dependence` | `0` | Informative censoring: censoring times are multiplied by $e^{-\delta z}$ ($z$ the standardized predictor), so with $\delta > 0$ higher-risk subjects are censored earlier. |
| `train.fraction` | `0.7` | Share of samples in the training split. |
| `nrep` | `1` | Number of replicates. |
| `return.template` | `FALSE` | Include the fitted template (which contains the source data) in the result. |
| `...` | | Passed to `fit_template()` when `template` is raw data. |

### `gen_multistudy()`

Takes every `gen_simmba()` argument (with the same defaults) and returns the
same format, with all studies of a replicate pooled. In addition:

| Argument | Default | Description |
|---|---|---|
| `template` | required | One template (or data set) for every study, or one per study: a list from `fit_study_templates()`, or a list of data sets (fitted here). |
| `n` | required | Samples per study: one value (with `nstudy`) or one per study. |
| `nstudy` | `length(n)` | Number of studies. |
| `effect.cor` | `1` | Correlation of the common signal features' effect sizes between studies. |
| `support.overlap` | `1` | Fraction of each study's signal features taken from the common set (the rest are study-specific). |
| `intercept.sd` | `0` | SD of study intercepts (continuous: in signal-SD units; binary: logit scale; survival: log hazard). |
| `snr.sd` | `0` | SD of the study log-SNRs. |
| `shift.sd` | `0` | Case mix: SD of study mean shifts of the latent factor scores (biological; enters the outcome). |
| `batch.loc.sd`, `batch.scale.sd` | `0` | Technical batch effects per study and feature, on the latent scale (added after the outcome). |
| `holdout.studies` | `NULL` | Studies whose samples all go to the test split. |
| `cens.rate` | `NULL` | Survival censoring target, solved within each study: one value or one per study (in order or named by study). `cens.times` may also be a list with one entry per study. |
| `train.fraction` | `0.7` | Training share within each non-holdout study (`1` allowed with holdout studies). |

### `calibrate_crosstalk()` and `validate_template()`

| Argument | Default | Description |
|---|---|---|
| `template` | required | Fitted template or raw input (`validate_template()` fits raw input with `template.args`). |
| `crosstalk.grid` | `seq(0, 1, 0.25)` | Candidate crosstalk values (values above `crosstalk_max()` are dropped). |
| `n` | template size | Samples per simulated replicate. |
| `nrep` | `5` | Replicates per grid value. |
| `p`, `feature.select`, `feature.jitter`, `view.models` | as above | Simulation settings. |
| `max.features.per.view` | `NULL` / `100` | Most variable observed features used in the comparison (calibration / validation). |
| `seed` | `1234` | Random seed. |

### `gen_mediation()`

Takes `n`, `template`, `p`, `crosstalk`, `feature.select`, `feature.jitter`,
`view.models`, `back.transform`, `nrep`, `seed` and `return.template` as above,
plus:

| Argument | Default | Description |
|---|---|---|
| `n.mediators` | `5` | Mediators per view (per view; 0 for a view without mediators). |
| `TIE` | `1` | Target total natural indirect effect (calibrated exactly). |
| `DE` | `0` | Direct effect of the exposure on the outcome. |
| `dir.a`, `dir.b` | `1` | Directions of the exposure-to-mediator and mediator-to-outcome paths: `1`, `-1` or `"random"`. |
| `effect.type` | `"spiked"` | `"spiked"`, or the nulls `"null_a"`, `"null_b"`, `"null_ab"`. |
| `exposure.type` | `"binary"` | `"binary"` or `"continuous"`. |
| `exposure.prevalence` | `0.5` | P(exposure = 1) for a binary exposure. |
| `n.interactions` | `0` | Pairwise mediator interactions in the outcome model. |
| `outcome.type` | `"continuous"` | `"continuous"` or `"binary"` (effects on the log-odds scale). |
| `outcome.prevalence` | `0.5` | Target P(Y = 1) for a binary outcome. |
| `sigma` | `1` | Outcome noise SD (continuous outcome). |
| `min.variability` | `0.2` | Mediators are drawn from features that vary in at least this share of template samples. |
| `n.mc` | `5000` | Monte Carlo draws used to compute the true effects. |
| `stages` | `NULL` | Sequential mediation: list of views in causal order, e.g. `list("Microbiome", "Metabolome")`. |
| `stage.links` | `list()` | Links between consecutive stages: `type` (`"random"`, `"hub"`, `"convergent"`, `"explicit"`), `density` (0.3), `strength` (0.3), `direction` (1), `matrix`. |
| `stage.direct` | `1` | Per stage: multiplier on the direct exposure paths (0: reached only through the previous stage). |
| `stage.share` | `NULL` | Per stage: share of the indirect effect carried by that stage's mediators. |

### `simulate_mimesys()`

`simulate_mimesys(template, n, p, crosstalk = "calibrate", design = "outcome",
preset = "default", outcome, signal, heterogeneity, nstudy, mediation,
null_model, template_args, calibration, train_fraction = 0.7, nrep = 1,
seed = 1234, return_type = c("mae", "simulation"))`.

| `design` | Simulates | Engine | Option lists it reads |
|---|---|---|---|
| `"outcome"` (default) | features with planted outcome signal, train/test split | `gen_simmba()` | `preset`, `outcome`, `signal`, `train_fraction` |
| `"null"` | features only | `sim_from_template()` | none |
| `"multistudy"` | several studies with heterogeneity, pooled in one MAE | `gen_multistudy()` | `preset`, `outcome`, `signal`, `heterogeneity`, `nstudy`, `train_fraction` |
| `"mediation"` | exposure, mediators, outcome and exact effects | `gen_mediation()` | `mediation` |

`template`, `n`, `p`, `crosstalk`, `null_model`, `template_args`,
`calibration`, `nrep`, `seed` and `return_type` apply to every design. An
option list that the chosen design does not use is an error. With
`"multistudy"`, `n` is the size of each study, and `template` may hold one
template (or data set) per study, in which case crosstalk is calibrated for
each study separately. The option lists use `snake_case` and map to the
arguments above:

| List | Options (default) | Maps to |
|---|---|---|
| `outcome` | `model` (`"linear"`; also `"friedman"`, `"friedman2"`, `"shared_latent"`, `"interaction"`), `type` (`"continuous"`), `snr` (1), `friedman_share`, `shared_latent_function` (`"linear"`), `shared_latent_n_factors`, `nonlinear_weight` (0.5), `interaction_weight` (0.5), `n_interactions`, `prevalence` (0.5), `baseline_hazard` (1), `cens_rate`, `cens_lower` (1), `cens_upper` (3), `shared_latent_weight` (1), `cens_type` (`"uniform"`), `cens_times`, `cens_dependence` (0), `cens_accrual` (0.5), `cens_dropout` (0), `baseline_shape` (1), `keep_true_times` (`FALSE`) | `gen_simmba()` outcome arguments |
| `signal` | `prob` (0.1), `neg_prob` (0.5), `effect_meanlog` (1), `effect_sdlog` (0.4), `align` (1), `share` | `signal.prob`, `signal.neg.prob`, `effect.meanlog`, `effect.sdlog`, `signal.align`, `signal.share` |
| `null_model` | `feature_select` (`"random"`), `feature_jitter` (0.2), `view_models`, `back_transform` (`FALSE`) | null-data arguments |
| `calibration` | `grid` (`seq(0, 1, 0.25)`), `nrep` (5), `nsample` (template size), `max_features_per_view` | `calibrate_crosstalk()` |
| `mediation` | `n_mediators` (5), `tie` (1), `de` (0), `dir_a`, `dir_b` (1), `effect_type` (`"spiked"`), `exposure_type` (`"binary"`), `exposure_prevalence` (0.5), `n_interactions` (0), `sigma` (1), `outcome_type` (`"continuous"`), `outcome_prevalence` (0.5), `min_variability` (0.2), `n_mc` (5000), `stages`, `stage_links`, `stage_direct` (1), `stage_share` | `gen_mediation()` |
| `heterogeneity` | `effect_cor` (1), `support_overlap` (1), `intercept_sd` (0), `snr_sd` (0), `shift_sd` (0), `batch_loc_sd` (0), `batch_scale_sd` (0), `holdout_studies` | `gen_multistudy()` |
| `template_args` | | `fit_template()`, when `template` is raw data |

`crosstalk = "calibrate"` (default) runs `calibrate_crosstalk()` first; a number
skips it. Presets (`mimesys_presets()`): `default`/`linear`, `friedman`,
`friedman2`, `shared_latent`, `shared_latent_friedman`,
`shared_latent_friedman2`, `interaction`; explicit options override them. One
replicate returns a `MultiAssayExperiment`; `nrep > 1` returns a
`mimesys_mae_simulation` list (`$replicates`).

## Key concepts

**Template model.** Each feature keeps its empirical distribution (zeros,
bounds, skewness, ties). Dependence is a Gaussian copula: a low-rank factor
model within each view and a few shared factors across views, estimated by
generalized CCA. The number of shared factors is chosen by sequential parallel
analysis.

**Crosstalk.** Scales only the cross-view dependence; marginals and
within-view structure stay fixed. `0` gives independent views, `1` the
template's level. Calibration picks the value whose null simulations best match
the template's cross-view correlation blocks; check
`near_optimal` and `selection_frequency` for how decisive it is.

**Signal regimes.** `signal.share` sets how much outcome signal each view
carries (complementary vs single-view signal), `signal.align` whether the
signal follows the shared cross-view direction (redundant, aligned signal),
and `outcome.model` its form (linear, Friedman nonlinear, shared-latent,
cross-view interactions); see the equations below. The truth is in `rowData()` (`beta`, `outcome_feature`,
`friedman_term`, `interaction_term`) and `metadata(x)$truth`.

**Sparse and compositional views.** `view.models = c(Microbiome =
"two_part_copula")` generates presence and abundance as two linked parts
(non-negative views only). Compositional views (rows summing to a constant, or
relative abundances of a subset of taxa) never exceed the template total.

**Scales.** Templates are fitted, calibrated and validated on the transformed
scale, and simulators return that scale unless `back.transform = TRUE`
(`null_model = list(back_transform = TRUE)`).

**Privacy.** A fitted template contains the source data (`template$observed`
and each feature's sorted values); treat it like the source data. Results
include only the template's provenance unless `return.template = TRUE`.

## How outcomes are generated

`gen_simmba()` (and `simulate_mimesys()`) generate the outcome from the
simulated features in four steps. Notation: $x_{ij}$ is feature $j$ of sample
$i$ on the fitted scale, $v(j)$ its view, and $\mathrm{std}(\cdot)$ centres a
vector and scales it to SD 1.

**1. Standardize the features** (within each replicate, over all samples):

$$\tilde x_{ij} = \frac{x_{ij} - \bar x_j}{\mathrm{sd}(x_j)} .$$

**2. Draw the coefficients.** In view $m$ with $p_m$ features, the number of
signal features is $\max\{1,\ \mathrm{Binomial}(p_m,\ \texttt{signal.prob}_m)\}$.
For each signal feature,

$$f_j \sim \mathrm{LogNormal}(\texttt{effect.meanlog},\ \texttt{effect.sdlog}),
\qquad d_j = \begin{cases}-1 & \text{with prob. } \texttt{signal.neg.prob}\\ +1 & \text{otherwise}\end{cases},
\qquad \beta_j = d_j\,\lvert\log_2 f_j\rvert ,$$

and $\beta_j = 0$ for the other features. With `signal.align` $=\alpha$, a share
$\alpha$ of the signal features is drawn with probability proportional to
$\ell_j^2$, where $\ell_j$ is the feature's loading on the shared cross-view
direction (`shared_loading`), and gets the sign of $\ell_j$.

**3. Build the predictor.** Each view's linear predictor is
$\eta_{im} = \sum_{j \in m} \beta_j \tilde x_{ij}$. With `signal.share`
$=(w_1,\dots,w_M)$ (normalized to sum to 1), each view's coefficients are
rescaled so that $\mathrm{Var}(\eta_m) = w_m$. The linear predictor is
$\eta^{\text{lin}}_i = \sum_m \eta_{im}$, and the outcome model gives $\eta_i$:

| `outcome.model` | Predictor $\eta_i$ |
|---|---|
| `"linear"` | $\eta_i = \eta^{\text{lin}}_i$ |
| `"friedman"` | $\eta_i = F(u_{i1},\dots,u_{i5})$, with $F(u) = 10\sin(\pi u_1 u_2) + 20(u_3 - 0.5)^2 + 10u_4 + 5u_5$ and $u_{ik} = \mathrm{rank}(x_{i j_k})/(n+1)$ for five signal features $j_1,\dots,j_5$ |
| `"friedman2"` | $\eta = \sqrt{1-w}\ \mathrm{std}(\eta^{\text{lin}}) + \sqrt{w}\ \mathrm{std}(F)$, $w = $ `nonlinear.weight` |
| `"shared_latent"` | $\eta = \mathrm{std}\big(\sum_{k=1}^{K} \tilde s_{ik} / \sqrt{K}\big)$ for the first $K$ simulated shared factors $s$ (`shared.latent.function = "linear"`); `"friedman"` uses $F$ of the first five rank-transformed $s_k$, `"friedman2"` mixes the two with $w$. With `shared.latent.weight` $= w_s < 1$: $\eta = \sqrt{w_s}\,\mathrm{std}(\eta^{\text{sh}}) + \sqrt{1-w_s}\,\mathrm{std}\big(\sum_m \sqrt{\pi_m}\,\mathrm{std}(P_m \omega_m)\big)$, where $P_m = F_m - \sqrt{\gamma_m} S B_m^\top$ is view $m$'s non-shared latent part, $\omega_m \propto$ the factors' loading strengths and $\pi_m$ = `signal.share` (equal by default) |
| `"interaction"` | $\eta = \sqrt{1-w}\ \mathrm{std}(\eta^{\text{lin}}) + \sqrt{w}\ \mathrm{std}\big(\sum_t \tilde x_{i a_t}\tilde x_{i b_t}\big)$, pairs $(a_t, b_t)$ from different views, $w = $ `interaction.weight` |

**4. Generate the outcome.**

| `outcome.type` | Outcome | Meaning of `snr` |
|---|---|---|
| `"continuous"` | $Y_i = \eta_i + \varepsilon_i,\quad \varepsilon_i \sim N\big(0,\ \mathrm{Var}(\eta)/\mathrm{snr}\big)$ | $\mathrm{Var}(\eta)/\mathrm{Var}(\varepsilon)$ |
| `"binary"` | $\eta^\ast_i = b_0 + \mathrm{std}(\eta)_i\sqrt{\mathrm{snr}\cdot\pi^2/3},\quad Y_i \sim \mathrm{Bernoulli}\big(\mathrm{logit}^{-1}(\eta^\ast_i)\big)$, with $b_0$ solved so that the mean probability equals `binary.prevalence` | signal-to-noise on the latent logistic scale (error variance $\pi^2/3$) |
| `"survival"` | $\eta^\ast_i = \mathrm{std}(\eta)_i\sqrt{\mathrm{snr}\cdot\pi^2/6},\quad T_i = E_i^{1/k},\ E_i \sim \mathrm{Exponential}\big(h_0 e^{\eta^\ast_i}\big)$ (Weibull hazard $h_0 k t^{k-1} e^{\eta^\ast_i}$, $k$ = `baseline.shape`, default 1: exponential), $h_0 = $ `baseline.hazard`; censoring $C_i \sim U(\texttt{cens.lower}, \texttt{cens.upper})$, or $C_i = \kappa V_i e^{-\delta z_i}$ with $V_i$ from `cens.type` (uniform $U(0,1)$, exponential $\mathrm{Exp}(1)$, administrative $\min\{1 - a U_i,\ D_i\}$) or `cens.times`, $\delta$ = `cens.dependence` and $\kappa$ solved for `cens.rate`; $\texttt{time} = \min(T_i, C_i)$, $\texttt{status} = \mathbf 1(T_i \le C_i)$ | signal-to-noise on the log-time scale (error variance $\pi^2/6$) |

The truth is stored with the data: $\beta_j$ in `rowData()$beta` (0 where the
linear predictor does not enter the outcome), the features that drive the
outcome in `outcome_feature`, the predictor and its parts in `colData()`
(`Xbeta`, `Xbeta_linear`, ...), and the interaction pairs in
`metadata(x)$truth$interaction_pairs`.

## Multi-study simulation

`gen_multistudy()` simulates $S$ studies, either all from one template or
each from its own template (one real cohort per study). Every study measures
the same simulated features (study 1 chooses them).

**One template per study.** `fit_study_templates()` fits a template to each
cohort on the views and features that every cohort has (features constant or
missing in any cohort are dropped everywhere), so the templates line up:

```r
tmpls <- fit_study_templates(list(CohortA = mae_a, CohortB = mae_b),
                             view.transform = c(Metabolome = "log1p"))
sim <- gen_multistudy(n = c(200, 150), template = tmpls, p = 50)
```

Each study is then simulated from its own cohort's marginals, factor loadings
and cross-view structure, on the same named features (features added beyond
the template size are built from the same feature pairs in every study). The
real between-cohort differences in distributions and correlations carry over,
and the heterogeneity settings below add to them. Signal coefficients are
drawn once on the shared features; with `signal.align` they follow study 1's
shared direction. Crosstalk can be given per study as a list. If the cohort
templates have different ranks, the latent scores are not comparable and
`truth$scores`/`truth$shared` are left out (and `shared_latent` outcomes need
templates fitted with the same `shared.rank`).

In addition, study $s$ can differ from the others in five ways:

| Heterogeneity | Model | Argument |
|---|---|---|
| Case mix | Shared factors $F_i + \delta^{sh}_s$ and view-specific innovations $E_i + \delta^{E}_s$, $\delta \sim N(0, \sigma_{shift}^2 I)$; features follow through the loadings | `shift.sd` |
| Batch | Latent feature score $z_{ij} \to a_{sj} + b_{sj} z_{ij}$, $a_{sj} \sim N(0, \sigma_a^2)$, $\log b_{sj} \sim N(0, \sigma_b^2)$, then the template marginals | `batch.loc.sd`, `batch.scale.sd` |
| Signal features | A fraction $\omega$ of each study's signal features come from the common set; the rest are drawn from the other features of that view (following `signal.align`) | `support.overlap` |
| Effect sizes | $\beta_{sj} = \mu + \sqrt{\rho}\, (\beta_j - \mu) + \sqrt{1 - \rho}\, (\beta^{*}_{sj} - \mu)$, with $\beta^{*}_{sj}$ an independent draw from the effect-size distribution and $\mu$ the mean of that distribution, so the mean and variance are kept and $\mathrm{cor}(\beta_{sj}, \beta_{tj}) = \rho$ | `effect.cor` |
| Baseline and noise | Study intercept $\alpha_s \sim N(0, \sigma_\alpha^2)$ and $\mathrm{SNR}_s = \mathrm{SNR}\, e^{\epsilon_s}$, $\epsilon_s \sim N(0, \sigma_{snr}^2)$ | `intercept.sd`, `snr.sd` |

The predictor is built on the pooled samples, so study differences are kept.
Features are standardized over all studies before batch effects (with
`signal.share`, each view is rescaled by one factor for all studies), and
$\eta_i = \tilde x_i^\top \beta_{s(i)}$ (plus the Friedman, interaction or
shared-latent part, which is the same in every study). With $V$ the pooled
within-study variance of $\eta$:

- continuous: $Y_i = \eta_i + \alpha_{s}\sqrt{V} + \varepsilon_i$, $\varepsilon_i \sim N(0, V / \mathrm{SNR}_s)$;
- binary: $\mathrm{logit}\, P(Y_i = 1) = c + \alpha_s + \sqrt{\mathrm{SNR}_s\, \pi^2/3}\, (\eta_i - \bar\eta) / \sqrt{V}$, with $c$ solved so the pooled prevalence equals `binary.prevalence`;
- survival: hazard $h_0 \exp\{\alpha_s + \sqrt{\mathrm{SNR}_s\, \pi^2/6}\, (\eta_i - \bar\eta)/\sqrt{V}\}$, with censoring (`cens.rate`) solved within each study.

Batch effects are applied after the outcome is drawn, so they are pure
technical noise; case-mix shifts are biology and move the outcome too. All
heterogeneity is drawn from a separate random stream, so with the same seed
the base data, coefficients and outcome noise are identical at every
heterogeneity level (for example, a sweep over `batch.loc.sd` changes only
the batch effects). With one study and no heterogeneity the result equals
`gen_simmba()`.

The truth is in `truth$studies` (size, holdout flag, intercept, SNR),
`truth$beta` (feature x study coefficients, also `rowData()$beta_Study_1`, ...),
`truth$shift`, `truth$batch` and, for binary outcomes, `truth$outcome_intercept`
(the common intercept $c$); `rowData()$common_signal` marks the common
signal features. `colData()$study` labels the samples, and sample IDs start
with the study (`Study_1_Subject1`).

## Mediation

`gen_mediation()` / `simulate_mimesys(design = "mediation")` simulate an exposure $T$ (binary
or continuous; $T^\ast$ is its standardized version), mediators $M_j$ in any
views, and an outcome. The exposure acts on each mediator's latent normal score
$Z_j$, so mediators keep their template distributions ($Q_j$ is the feature's
empirical quantile function):

$$W_j = a_j T^\ast + \sqrt{1 - a_j^2}\, Z_j, \qquad M_j = Q_j\big(\Phi(W_j)\big)$$

(for a binary exposure $W_j$ is first mapped through its mixture CDF back to an
exact standard normal). With $\tilde M_j$ the mediator standardized by its
template mean and SD, the outcome is

$$Y = \mathrm{DE}\cdot T + \sum_j b_j \tilde M_j + \sum_k g_k \tilde M_{k_1}\tilde M_{k_2} + \sigma\varepsilon
\quad\text{or}\quad
\mathrm{logit}\,P(Y=1) = b_0 + \mathrm{DE}\cdot T + \sum_j b_j \tilde M_j + \dots$$

and the natural effects are

$$\mathrm{NIE} = \sum_j b_j \Delta_j + \sum_k g_k \Delta_k, \quad \Delta_j = E[\tilde M_j(1)] - E[\tilde M_j(0)], \qquad \mathrm{NDE} = \mathrm{DE}, \qquad \mathrm{TE} = \mathrm{NIE} + \mathrm{NDE}.$$

The b-paths are calibrated so that NIE equals `TIE` exactly, with the
$\Delta$'s computed by paired counterfactual Monte Carlo; the truth
(`NIE`, `NDE`, `TE`, per-view and per-mediator effects, and natural risk
differences and odds ratios for a binary outcome) is in
`metadata(x)$truth$mediation` and `rowData()`. The `null_a`, `null_b` and
`null_ab` scenarios have zero indirect effect and matched path sizes.

**Sequential mediation.** With `stages`, mediators of a later stage are linked,
many-to-many, to mediators $k$ of the previous stage (e.g. microbes ->
metabolites):

$$W_j = a_j T^\ast + \sum_k d_{jk} W_k + c_j Z_j ,$$

with $c_j$ chosen so that $\mathrm{Var}(W_j) = 1$. The indirect effect is split
into natural path-specific effects (for two stages $T \to M_1 \to Y$,
$T \to M_1 \to M_2 \to Y$ and $T \to M_2 \to Y$), which add up to the NIE.

```r
data(franzosa_mae)
micro_template <- fit_template(franzosa_mae, view.transform = c(Metabolome = "log1p"))
seq_mae <- simulate_mimesys(micro_template, design = "mediation", n = 500, mediation = list(
  n_mediators = c(Microbiome = 5, Metabolome = 5),
  stages = list(Microbe = "Microbiome", Metabolite = "Metabolome"),
  stage_links = list(type = "random", density = 0.3, strength = 0.3),
  tie = 1, de = 0.5
))
truth <- metadata(seq_mae)$truth$mediation
truth$path_effects   # natural path-specific effects; they sum to the NIE
truth$links          # microbe -> metabolite links and their coefficients
```

`stage_direct = c(1, 0)` with `stage_share = c(0, 1)` gives a pure chain;
`stage_links = list(strength = 0)` gives unlinked stages that are still
correlated through crosstalk. Note that crosstalk acts as mediator-mediator
confounding, and skewed parents make the observed data nonlinear in the links,
so regression estimates of the chain can be biased even though the truth is
exact (see `docs/MimeSys_algorithms.md`, section 12.6).


## License

MIT; see `LICENSE`.
