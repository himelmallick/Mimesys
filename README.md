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

## What you can control

| Design | Main options |
|---|---|
| all | `n`, `p` (features per view), `crosstalk` (`"calibrate"` or a number; 0 = independent views, 1 = template level), `null_model` (e.g. two-part model for sparse microbiome views), `nrep`, `seed` |
| `"outcome"` | `outcome`: model (`linear`, `friedman`, `friedman2`, `shared_latent`, `interaction`), type (`continuous`, `binary`, `survival`), `snr`, prevalence, censoring (rate, mechanism, informative, user times), Weibull shape, shared versus view-specific latent signal; `signal`: how much signal each view carries and whether it follows the shared direction; `preset` (see `mimesys_presets()`) |
| `"multistudy"` | the outcome options plus `heterogeneity`: effect-size correlation, signal-feature overlap, intercept and SNR spread, case-mix shift, batch effects, holdout studies; one template per study via `fit_study_templates()` |
| `"mediation"` | `mediation`: number of mediators per view, indirect and direct effects, exposure and outcome types, null scenarios, sequential stages and their links |

Every option is described in `?simulate_mimesys`, and the engines
`gen_simmba()`, `gen_multistudy()`, `gen_mediation()` and
`sim_from_template()` give full control.

## Other functions

| Function | What it does |
|---|---|
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
