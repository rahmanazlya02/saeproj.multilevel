# saeproj.multilevel

<!-- badges: start -->
<!-- badges: end -->

## Author

Nazlya Rahma Susanto

## Maintainer

Nazlya Rahma Susanto <susantonazlya@gmail.com>

## Description

The **saeproj.multilevel** package provides tools for *Small Area Estimation* (SAE) using a projection estimator with a linear mixed-effects working model.

The method is designed for two-survey settings:

- a smaller model survey, `data_model`, that contains the response variable and auxiliary predictors;
- a larger projection survey, `data_proj`, that contains auxiliary predictors but does not necessarily contain the response variable.

A linear mixed-effects model is fitted on the model survey using `lme4::lmer()`. Random effects are included at this fitting stage to account for the hierarchical structure of the model-survey data, such as observations nested within areas. After the model parameters are estimated, synthetic predictions are generated for all units in the projection survey using the fixed-effect component of the fitted model, excluding group-specific random-effect predictions or BLUPs. This prediction strategy allows estimates to be produced for areas or domains that are not observed in the model survey.

Domain-level mean estimates are obtained by aggregating these predictions using survey design-weighted means from the `survey` package.

The package supports both synthetic and bias-corrected projection estimators.

The main function is:

```r
sae_ml_linear()
```

## Installation

You can install the development version of `saeproj.multilevel` from GitHub with:

```r
# install.packages("devtools")
devtools::install_github("rahmanazlya02/saeproj.multilevel")
```
If you are developing the package locally, install it from the package directory with:

```r
devtools::install()
```

## Dependencies

The package requires:

- `lme4` — for fitting linear mixed-effects models;
- `survey` — for survey-weighted aggregation of predictions.

Both packages are listed under `Imports` and will be installed automatically.

## Example

### Simulated data

```r
library(saeproj.multilevel)

set.seed(42)

n_area  <- 6
n_model <- 120
n_proj  <- 500

area_model <- sample(paste0("A", 1:n_area), n_model, replace = TRUE)

data_model <- data.frame(
  kabkot = area_model,
  educ   = sample(1:3, n_model, replace = TRUE),
  age    = runif(n_model, 20, 60),
  weight = runif(n_model, 0.5, 2)
)

area_eff <- setNames(rnorm(n_area, 0, 1.5), paste0("A", 1:n_area))

data_model$income <- with(
  data_model,
  5 + 0.8 * educ + 0.05 * age + area_eff[kabkot] + rnorm(n_model, 0, 1)
)

data_proj <- data.frame(
  kabkot = sample(paste0("A", 1:n_area), n_proj, replace = TRUE),
  educ   = sample(1:3, n_proj, replace = TRUE),
  age    = runif(n_proj, 20, 60),
  weight = runif(n_proj, 0.5, 2)
)
```

### Bias-corrected projection estimator

```r
result <- sae_ml_linear(
  formula    = income ~ educ + age + (1 | kabkot),
  data_model = data_model,
  data_proj  = data_proj,
  domain     = "kabkot",
  weight     = "weight",
  estimator  = "bias_corrected"
)
```

## Primary output: domain-level estimates

The main output is stored in:

```r
result$estimates
```

This table contains the final domain-level estimates:

```r
result$estimates
```

The columns are:

| Column | Description |
|---|---|
| domain variable(s) | Domain identifier column(s), based on the `domain` argument |
| `estimate` | Final domain mean estimate |
| `variance` | Estimated variance of the final domain estimate |
| `se` | Standard error, computed as `sqrt(variance)` |
| `rse` | Relative standard error in percent |

`as.data.frame(result)` returns the same table and can be used for downstream workflows such as `write.csv()`, `dplyr`, or `ggplot2`.

```r
as.data.frame(result)
```

## Estimation components

For advanced users or debugging, the full estimation breakdown is available in:

```r
result$estimation_details
```

This table contains:

| Column | Description |
|---|---|
| domain variable(s) | Domain identifier column(s) |
| `estimate_synthetic` | Synthetic, projection-only domain estimate |
| `variance_synthetic` | Variance of the synthetic estimate |
| `correction` | Empirical residual correction |
| `variance_correction` | Variance of the residual correction |
| `estimate_final` | Final estimate, computed as `estimate_synthetic + correction` |
| `variance_final` | Final variance, computed as `variance_synthetic + variance_correction` |
| `se_final` | Standard error of the final estimate |
| `rse_final` | Relative standard error of the final estimate |
| `n_model` | Number of units in the domain in `data_model` |
| `n_proj` | Number of units in the domain in `data_proj` |

Example:

```r
result$estimation_details
```

## Print and summary

A concise output can be displayed with:

```r
print(result)
```

A detailed output can be displayed with:

```r
summary(result)
```

`summary()` displays the full `lmer` model output, key diagnostics, final estimates, conditional notes, and a pointer to `result$estimation_details`.

## Synthetic estimator

The synthetic estimator uses projection only and does not add the empirical residual correction.

```r
result_syn <- sae_ml_linear(
  formula    = income ~ educ + age + (1 | kabkot),
  data_model = data_model,
  data_proj  = data_proj,
  domain     = "kabkot",
  weight     = "weight",
  estimator  = "synthetic"
)

result_syn$estimates
```

## Retaining unit-level predictions

Set `keep_unit = TRUE` to store unit-level predictions from `data_proj` and model residuals from `data_model`.

```r
result_ku <- sae_ml_linear(
  formula    = income ~ educ + age + (1 | kabkot),
  data_model = data_model,
  data_proj  = data_proj,
  domain     = "kabkot",
  weight     = "weight",
  keep_unit  = TRUE
)

head(result_ku$unit_projection)
head(result_ku$unit_model_residual)
```

When `keep_unit = TRUE`:

- `result_ku$unit_projection` contains the domain variable(s) and `.y_hat`;
- `result_ku$unit_model_residual` contains the domain variable(s), response variable, `.y_hat_model`, and `.resid`.

## Model diagnostics

Model diagnostics are stored in:

```r
result$diagnostics
```

The diagnostics include:

```r
result$diagnostics$fixed_effects
result$diagnostics$variance_components
result$diagnostics$sigma
result$diagnostics$nobs
result$diagnostics$aic
result$diagnostics$bic
result$diagnostics$loglik
result$diagnostics$icc
result$diagnostics$singular_fit
result$diagnostics$convergence_messages
```

A compact diagnostics table can be created as follows:

```r
data.frame(
  sigma  = result$diagnostics$sigma,
  nobs   = result$diagnostics$nobs,
  aic    = result$diagnostics$aic,
  bic    = result$diagnostics$bic,
  loglik = result$diagnostics$loglik,
  icc    = result$diagnostics$icc,
  singular_fit = result$diagnostics$singular_fit
)
```

For residual diagnostics, access the fitted `lmerMod` object directly:

```r
fit <- result$fitted_model

plot(fitted(fit), resid(fit),
     xlab = "Fitted values",
     ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, lty = 2)

qqnorm(resid(fit))
qqline(resid(fit))

lme4::ranef(fit)
```

## Multiple domain variables

The `domain` argument accepts a character vector or a one-sided formula.

The following example assumes that both `data_model` and `data_proj` contain the variables `prov` and `kabkot`.

```r
result_multi <- sae_ml_linear(
  formula    = income ~ educ + age + (1 | kabkot),
  data_model = data_model,
  data_proj  = data_proj,
  domain     = c("prov", "kabkot"),
  weight     = "weight"
)

result_multi$estimates
```

The estimates table will include all domain columns:

```r
# Example structure:
#   prov kabkot estimate variance se rse
```

## Survey design with clustering and stratification

The survey design arguments `cluster_ids`, `weight`, and `strata` are used in the aggregation step through `survey::svydesign()`.

The following example assumes that both `data_model` and `data_proj` contain `psu_id`, `weight`, and `stratum` variables.

```r
result_svy <- sae_ml_linear(
  formula     = income ~ educ + age + (1 | kabkot),
  data_model  = data_model,
  data_proj   = data_proj,
  domain      = "kabkot",
  cluster_ids = "psu_id",
  weight      = "weight",
  strata      = "stratum"
)
```

## Output object structure

`sae_ml_linear()` returns an S3 object of class `"sae_ml_linear"`.

| Component | Description |
|---|---|
| `$call` | The matched call |
| `$formula` | The model formula |
| `$estimator` | Estimator type: `"bias_corrected"` or `"synthetic"` |
| `$fitted_model` | The fitted `lmerMod` object from `lme4::lmer()` |
| `$estimates` | Final domain-level estimates: domain column(s), `estimate`, `variance`, `se`, `rse` |
| `$estimation_details` | Estimation components: synthetic estimate, correction, final estimate, and sample sizes per domain |
| `$diagnostics` | Model diagnostics: fixed effects, variance components, AIC, BIC, logLik, ICC, singular fit, and convergence messages |
| `$notes` | Conditional notes generated from model or estimator conditions, such as singular fit, convergence issues, undefined ICC, or the bias-corrected variance assumption |
| `$unit_projection` | Unit-level `data_proj` with `.y_hat`, only if `keep_unit = TRUE` |
| `$unit_model_residual` | Unit-level `data_model` with `.y_hat_model` and `.resid`, only if `keep_unit = TRUE` |

## S3 methods

| Method | Behaviour |
|---|---|
| `print(result)` | Prints a concise output: formula, estimator, number of domains, and a preview of `$estimates` |
| `summary(result)` | Prints a detailed output: fitted model summary, diagnostics, final estimates, and notes |
| `as.data.frame(result)` | Returns `result$estimates` directly |

## Function interface

```r
sae_ml_linear(
  formula,
  data_model,
  data_proj,
  domain,
  cluster_ids = ~1,
  weight      = NULL,
  strata      = NULL,
  estimator   = c("bias_corrected", "synthetic"),
  keep_unit   = FALSE,
  control     = lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)),
  ...
)
```

| Argument | Description |
|---|---|
| `formula` | `lme4::lmer()`-style formula, for example `y ~ x1 + x2 + (1 \| area)` |
| `data_model` | Model survey data frame containing the response, predictors, grouping variables, domain variable(s), and survey design variables |
| `data_proj` | Projection survey data frame containing predictors, domain variable(s), and survey design variables; response is not required |
| `domain` | Domain variable name(s): character scalar, character vector, or one-sided formula |
| `cluster_ids` | PSU or cluster variable for survey design; use `~1` for no clustering |
| `weight` | Survey weight variable; use `NULL` for equal weights |
| `strata` | Stratification variable; use `NULL` if not applicable |
| `estimator` | `"bias_corrected"` or `"synthetic"` |
| `keep_unit` | If `TRUE`, unit-level predictions and residuals are stored in the output |
| `control` | `lme4::lmerControl()` object passed to `lme4::lmer()` |
| `...` | Additional named arguments passed to `survey::svydesign()` only, for example `nest = TRUE`. These are not forwarded to `lme4::lmer()`; use `control` for lmer-specific tuning. |

## Methodological notes

- The working model is a linear mixed-effects model. Random effects are used during model fitting to account for hierarchical structure in `data_model`, such as observations nested within areas or groups.
- Projected predictions on `data_proj` are generated using the fixed-effect component of the fitted model only, with `re.form = NA`. Group-specific BLUPs or random-effect predictions are not carried over to the projection survey.
- This prediction strategy allows `data_proj` to contain new random-effect levels or areas that were not observed in `data_model`, as long as the fixed-effect predictors are available.
- Survey design arguments (`cluster_ids`, `weight`, and `strata`) are used only in the aggregation step through `survey::svydesign()` and `survey::svyby()`, not in model fitting.
- The reported variance is a plug-in estimate from `survey::svyby()`. It is approximate and does not fully account for uncertainty in the estimated model parameters.
- For the bias-corrected estimator, the final variance is computed as `variance_synthetic + variance_correction`.
- Missing values are not removed automatically. The function stops with an informative error if required variables contain missing values.
- Categorical predictors in `data_proj` must not contain levels that are absent from `data_model`.
- AIC, BIC, and logLik in `diagnostics` are REML-based and should not be used to compare models with different fixed-effect structures. For that purpose, refit models with `REML = FALSE`.
- ICC is computed only for models with a single grouping factor and random intercept only. It returns `NA` for random-slope or multi-grouping structures.

## References

Kim, J.K. & Rao, J.N.K. (2012). Combining data from two independent surveys: a model-assisted approach. *Biometrika*, 99(1), 85–100.

Bates, D., Maechler, M., Bolker, B. & Walker, S. (2015). Fitting linear mixed-effects models using lme4. *Journal of Statistical Software*, 67(1), 1–48.

Lumley, T. (2010). *Complex Surveys: A Guide to Analysis Using R*. Wiley.

Rao, J.N.K. & Molina, I. (2015). *Small Area Estimation* (2nd ed.). Wiley.

## License

MIT © Nazlya Rahma Susanto
