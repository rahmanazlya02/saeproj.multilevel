# saeproj.multilevel

## Author
Nazlya Rahma Susanto

## Maintainer
Nazlya Rahma Susanto <susantonazlya@gmail.com>

## Description
The **saeproj.multilevel** package provides tools for *Small Area Estimation* (SAE) using a projection estimator with a linear mixed-effects working model.
The method is designed for two-survey settings:
- a smaller model survey, `data_model`, that contains the response variable and auxiliary predictors;
- a larger projection survey, `data_proj`, that contains auxiliary predictors but does not necessarily contain the response variable.

The main function is:

```r
sae_ml_linear()
```

The function fits a linear mixed-effects model using `lme4::lmer()`, generates unit-level predictions for the projection dataset, aggregates those predictions by domain using survey design information, and applies a design-based residual bias correction.
The final estimator returned by the function is a bias-corrected projection estimator:

```r
estimate = synthetic projection + residual correction
```

The synthetic projection component and residual correction component are stored in:

```r
result$estimation_details
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
The package imports:
- `lme4` — for fitting linear mixed-effects models;
- `survey` — for survey design and domain-level aggregation;
- `dplyr` — for joining estimation components;
- `cli` — for errors and selected warnings;
- `reformulas` — for parsing mixed-effects formulas.

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
  kabkot = sample(c(paste0("A", 1:n_area), "A7"), n_proj, replace = TRUE),
  educ   = sample(1:3, n_proj, replace = TRUE),
  age    = runif(n_proj, 20, 60),
  weight = runif(n_proj, 0.5, 2)
)
```

### Fit the multilevel projection estimator

```r
result <- sae_ml_linear(
  formula    = income ~ educ + age + (1 | kabkot),
  data_model = data_model,
  data_proj  = data_proj,
  domain     = "kabkot",
  weight     = "weight"
)

result
```

## Primary output: domain-level estimates
The final domain-level estimates are stored in:

```r
result$estimates
```

The output contains:

| Column | Description |
|---|---|
| domain variable(s) | Domain identifier column(s), based on the `domain` argument |
| `estimate` | Final bias-corrected projection estimate |
| `variance` | Approximate variance of the final estimate |
| `se` | Standard error, computed as `sqrt(variance)` |
| `rse` | Relative standard error in percent |

`as.data.frame(result)` returns the same table and can be used for downstream workflows.

```r
as.data.frame(result)
```

## Estimation components
Detailed estimation components are stored in:

```r
result$estimation_details
```

This table contains:

| Column | Description |
|---|---|
| domain variable(s) | Domain identifier column(s) |
| `estimate_synthetic` | Synthetic projection estimate |
| `variance_synthetic` | Variance of the synthetic projection estimate |
| `correction` | Design-based residual correction |
| `variance_correction` | Variance of the residual correction |
| `estimate_final` | Final estimate, computed as `estimate_synthetic + correction` |
| `variance_final` | Final variance, computed as `variance_synthetic + variance_correction` |
| `se_final` | Standard error of the final estimate |
| `rse_final` | Relative standard error of the final estimate |
| `n_model` | Number of observations in the domain in `data_model` |
| `n_proj` | Number of observations in the domain in `data_proj` |

Example:

```r
head(result$estimation_details)
```

## Synthetic projection component
The function returns the bias-corrected estimator by default.
The synthetic projection component is available in:

```r
result$estimation_details[, c(
  "kabkot",
  "estimate_synthetic",
  "variance_synthetic"
)]
```

For multiple domain variables, include all domain columns when selecting from `estimation_details`.

## Direct estimator
Set `return_direct = TRUE` to return direct design-based estimates from `data_model`.

```r
result_direct <- sae_ml_linear(
  formula       = income ~ educ + age + (1 | kabkot),
  data_model    = data_model,
  data_proj     = data_proj,
  domain        = "kabkot",
  weight        = "weight",
  return_direct = TRUE
)

result_direct$direct_estimator
```

The direct estimator is stored separately and is not used to replace the projection estimator.

## Print and summary
A concise output can be displayed with:

```r
print(result)
```

A compact summary can be displayed with:

```r
summary(result)
```

The `summary()` method displays formula, estimator type, number of domains, selected model diagnostics, and a preview of the final estimates.

Full model output can be accessed from the fitted `lmerMod` object:

```r
summary(result$fitted_model)
```

## Retaining unit-level predictions

Set `keep_unit = TRUE` to store unit-level projection data and model residual data.

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

- `result_ku$unit_projection` contains `data_proj` with the unit-level prediction column `.prediction`;
- `result_ku$unit_model_residual` contains `data_model` with `.fitted_model` and `.model_residual`.

## Model diagnostics

Model diagnostics are stored in:

```r
result$diagnostics
```

The diagnostics include:

```r
result$diagnostics$icc
result$diagnostics$icc_note
result$diagnostics$singular_fit
result$diagnostics$convergence
result$diagnostics$sigma
result$diagnostics$residual_variance
result$diagnostics$random_effects
result$diagnostics$random_effect_groups
result$diagnostics$random_effect_dims
result$diagnostics$is_random_intercept_only
result$diagnostics$nobs
result$diagnostics$REML
result$diagnostics$logLik
result$diagnostics$AIC
result$diagnostics$BIC
```

A compact diagnostics table can be created as follows:

```r
data.frame(
  icc = result$diagnostics$icc,
  singular_fit = result$diagnostics$singular_fit,
  convergence = result$diagnostics$convergence,
  sigma = result$diagnostics$sigma,
  residual_variance = result$diagnostics$residual_variance,
  REML = result$diagnostics$REML,
  AIC = result$diagnostics$AIC,
  BIC = result$diagnostics$BIC
)
```

For additional model diagnostics, access the fitted `lmerMod` object directly:

```r
fit <- result$fitted_model
summary(fit)
```

Residual diagnostics can be inspected from the fitted model:

```r
plot(
  fitted(fit),
  resid(fit),
  xlab = "Fitted values",
  ylab = "Residuals",
  main = "Residuals vs Fitted"
)
abline(h = 0, lty = 2)

qqnorm(resid(fit))
qqline(resid(fit))

lme4::ranef(fit)
```

## Model parameters
Estimated model parameters are stored in:

```r
result$model_parameters
```

This object includes:

```r
result$model_parameters$fixed_effects
result$model_parameters$random_effects
result$model_parameters$variance_components
result$model_parameters$residual_sd
result$model_parameters$residual_variance
```

## Notes

Run-specific notes are stored in:

```r
result$notes
```

The notes are intentionally concise and are not printed automatically by `summary()`.

They may include information such as:

- removed zero-variance predictors;
- singular model fit;
- convergence issues;
- random-slope or complex random-effect structure where simple ICC is not computed;
- out-of-sample domains with zero residual correction;
- negative plug-in variance clamped to zero.

Out-of-sample domains are not treated as warnings because they are expected in SAE projection. They are recorded in `result$notes`.

## Multiple domain variables
The `domain` argument accepts a character scalar, a character vector, or a one-sided formula.
The following example assumes that both `data_model` and `data_proj` contain the variables `prov` and `kabkot`.

```r
data_model$prov <- substr(data_model$kabkot, 1, 1)
data_proj$prov  <- substr(data_proj$kabkot, 1, 1)

result_multi <- sae_ml_linear(
  formula    = income ~ educ + age + (1 | kabkot),
  data_model = data_model,
  data_proj  = data_proj,
  domain     = c("prov", "kabkot"),
  weight     = "weight"
)

result_multi$estimates
```

## Survey design with clustering and stratification
The survey design arguments `cluster_ids`, `weight`, and `strata` are used in the aggregation step through `survey::svydesign()`.
The following example assumes that both `data_model` and `data_proj` contain `psu_id`, `weight`, and `stratum` variables.

```r
data_model$psu_id  <- sample(1:30, nrow(data_model), replace = TRUE)
data_proj$psu_id   <- sample(1:30, nrow(data_proj), replace = TRUE)
data_model$stratum <- sample(1:5, nrow(data_model), replace = TRUE)
data_proj$stratum  <- sample(1:5, nrow(data_proj), replace = TRUE)

result_svy <- sae_ml_linear(
  formula     = income ~ educ + age + (1 | kabkot),
  data_model  = data_model,
  data_proj   = data_proj,
  domain      = "kabkot",
  cluster_ids = "psu_id",
  weight      = "weight",
  strata      = "stratum",
  nest        = TRUE
)
```

Use `cluster_ids = ~1` when there is no clustering:

```r
result_no_cluster <- sae_ml_linear(
  formula     = income ~ educ + age + (1 | kabkot),
  data_model  = data_model,
  data_proj   = data_proj,
  domain      = "kabkot",
  cluster_ids = ~1,
  weight      = "weight"
)
```

## Output object structure
`sae_ml_linear()` returns an S3 object of class `"sae_ml_linear"`.
Typical components are:
| Component | Description |
|---|---|
| `$call` | The matched function call |
| `$formula` | The model formula used after preprocessing |
| `$estimator` | Estimator type; currently always `"bias_corrected"` |
| `$fitted_model` | The fitted `lmerMod` object from `lme4::lmer()` |
| `$model_parameters` | Fixed effects, random effects, variance components, residual SD, and residual variance |
| `$estimates` | Final domain-level estimates |
| `$estimation_details` | Synthetic estimate, correction, final estimate, and sample sizes per domain |
| `$diagnostics` | Model diagnostics: ICC when applicable, random-effect structure, singular fit, convergence, sigma, residual variance, REML, logLik, AIC, and BIC |
| `$notes` | Concise run-specific notes |
| `$unit_projection` | Unit-level `data_proj` with `.prediction`, only if `keep_unit = TRUE` |
| `$unit_model_residual` | Unit-level `data_model` with `.fitted_model` and `.model_residual`, only if `keep_unit = TRUE` |
| `$direct_estimator` | Direct design-based estimates, only if `return_direct = TRUE` |

## S3 methods

| Method | Behaviour |
|---|---|
| `print(result)` | Prints formula, estimator, number of domains, and a preview of `$estimates` |
| `summary(result)` | Prints selected diagnostics and a preview of final estimates |
| `as.data.frame(result)` | Returns `result$estimates` |

## Function interface

```r
sae_ml_linear(
  formula,
  data_model,
  data_proj,
  domain,
  cluster_ids = ~1,
  weight = NULL,
  strata = NULL,
  summary_function = "mean",
  keep_unit = FALSE,
  seed = 1,
  control = lme4::lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5)
  ),
  return_direct = FALSE,
  ...
)
```

| Argument | Description |
|---|---|
| `formula` | `lme4::lmer()`-style formula, for example `y ~ x1 + x2 + (1 \| area)` |
| `data_model` | Model survey data frame containing the response, predictors, grouping variables, domain variable(s), and survey design variables |
| `data_proj` | Projection survey data frame containing predictors, grouping variables, domain variable(s), and survey design variables; response is not required |
| `domain` | Domain variable name(s): character scalar, character vector, or one-sided formula |
| `cluster_ids` | PSU or cluster variable for survey design; use `~1` for no clustering |
| `weight` | Survey weight variable; use `NULL` for equal weights |
| `strata` | Stratification variable; use `NULL` if not applicable |
| `summary_function` | Domain-level statistic: `"mean"` or `"total"` |
| `keep_unit` | If `TRUE`, unit-level predictions and residuals are stored in the output |
| `seed` | Integer seed used before model fitting |
| `control` | `lme4::lmerControl()` object passed to `lme4::lmer()` |
| `return_direct` | If `TRUE`, direct design-based estimates from `data_model` are returned |
| `...` | Additional named arguments passed to `survey::svydesign()`, for example `nest = TRUE` |

The `weight` argument identifies the survey weight column used in both `data_model` and `data_proj`. The column name must be the same in both datasets, but the weight values may differ.

In `data_model`, weights are used for residual correction and optional direct estimation. In `data_proj`, weights are used for synthetic projection aggregation.

## Methodological notes

- The working model is a linear mixed-effects model fitted with `lme4::lmer()`.
- The user fully specifies the fixed-effect and random-effect structure through the `formula` argument.
- Prediction uses `re.form = NULL` and `allow.new.levels = TRUE`.
- For grouping levels observed in `data_model`, predictions include the estimated random-effect contribution.
- For grouping levels appearing only in `data_proj`, the random-effect contribution is set to zero, so prediction uses the fixed part of the model.
- ICC is diagnostic only and does not determine the prediction rule.
- Simple ICC is computed only for pure random-intercept structures.
- Fixed-effect categorical predictors in `data_proj` must not contain levels that are absent from `data_model`.
- Survey design arguments (`cluster_ids`, `weight`, and `strata`) are used in the aggregation step through `survey::svydesign()` and `survey::svyby()`.
- The final estimate is computed as:

```r
estimate_final = estimate_synthetic + correction
```

- The final variance is computed as:

```r
variance_final = variance_synthetic + variance_correction
```

- This plug-in variance is approximate and assumes that `data_model` and `data_proj` are independent, or treated as independent for variance approximation.
- The reported variance does not fully account for uncertainty in the estimated mixed-model parameters.
- Missing values are not removed automatically. The function stops with an informative error if required variables contain missing values.

## Summary function

- The argument `summary_function` supports `"mean"` and `"total"` because both are linear domain parameters.
- For `"mean"`, the synthetic component and residual correction are aggregated using `survey::svymean`.
- For `"total"`, both components are aggregated using `survey::svytotal`, so the estimate and variance are returned on the total scale.
- The `"total"` option should only be used when the survey weights are appropriate expansion weights for population totals.

## References

Kim, J. K. and Rao, J. N. K. (2012). Combining data from two independent surveys: a model-assisted approach. *Biometrika*, 99(1), 85–100.
Moura, F. A. S. and Holt, D. (1999). Small area estimation using multilevel models. *Survey Methodology*, 25(1), 73–80.
Bates, D., Maechler, M., Bolker, B. and Walker, S. (2015). Fitting linear mixed-effects models using lme4. *Journal of Statistical Software*, 67(1), 1–48.
Lumley, T. (2010). *Complex Surveys: A Guide to Analysis Using R*. Wiley.
Rao, J. N. K. and Molina, I. (2015). *Small Area Estimation* (2nd ed.). Wiley.

## License

MIT © Nazlya Rahma Susanto
