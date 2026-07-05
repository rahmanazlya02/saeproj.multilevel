# R/sae_ml_linear.R
# Small Area Estimation via Projection Estimator using Linear Multilevel Working Model
#
# ---- Setup ------------------------------------------------------------------

utils::globalVariables(c(
  "estimate", "variance", "se", "rse",
  ".prediction", ".fitted_model", ".model_residual"
))

# ---- Internal helpers --------------------------------------------------------

#' @noRd
# Extract the response variable from a mixed-effects formula.
.get_response_var <- function(formula) {
  all.vars(reformulas::nobars(formula))[1L]
}

#' @noRd
# Extract fixed-effect variables from a mixed-effects formula.
.get_fixed_vars <- function(formula) {
  setdiff(all.vars(reformulas::nobars(formula)), .get_response_var(formula))
}

#' @noRd
# Extract grouping variables from random-effect terms.
.get_group_vars <- function(formula) {
  bars <- reformulas::findbars(formula)
  if (length(bars) == 0L) return(character(0L))
  unique(unlist(lapply(bars, function(x) all.vars(x[[3L]]))))
}

#' @noRd
# Check that required variables exist in a data frame.
.check_required_columns <- function(data, vars, data_name) {
  vars <- unique(vars[!is.na(vars) & nzchar(vars)])
  missing <- setdiff(vars, names(data))

  if (length(missing) > 0L) {
    cli::cli_abort(paste0(
      "Required variable(s) not found in `", data_name, "`: ",
      paste(missing, collapse = ", "), "."
    ))
  }

  invisible(TRUE)
}

#' @noRd
# Check required variables for missing values.
.check_missing_values <- function(data, vars, data_name) {
  vars <- intersect(unique(vars), names(data))
  if (length(vars) == 0L) return(invisible(TRUE))

  na_counts <- vapply(vars, function(v) sum(is.na(data[[v]])), integer(1L))
  bad <- na_counts[na_counts > 0L]

  if (length(bad) > 0L) {
    cli::cli_abort(paste0(
      "Missing values detected in ", data_name, ": ",
      paste(names(bad), collapse = ", "),
      ". Please clean or impute data before modeling."
    ))
  }

  invisible(TRUE)
}

#' @noRd
# Validate a survey design or domain variable and convert character input to formula.
.check_variable <- function(variable, data_model, data_proj, arg_name = "variable") {
  if (is.null(variable)) return(NULL)

  if (inherits(variable, "formula")) {
    if (identical(variable, ~1) || identical(variable, ~0)) return(variable)

    tryCatch(
      stats::model.frame(variable, data_model, na.action = NULL),
      error = function(e) cli::cli_abort(paste0(
        "Variable in `", arg_name, "` was not found in `data_model`: ",
        deparse(variable), "."
      ))
    )

    tryCatch(
      stats::model.frame(variable, data_proj, na.action = NULL),
      error = function(e) cli::cli_abort(paste0(
        "Variable in `", arg_name, "` was not found in `data_proj`: ",
        deparse(variable), "."
      ))
    )

    return(variable)
  }

  if (is.character(variable)) {
    missing_model <- setdiff(variable, names(data_model))
    missing_proj <- setdiff(variable, names(data_proj))

    if (length(missing_model) > 0L) {
      cli::cli_abort(paste0(
        "Variable(s) in `", arg_name, "` not found in `data_model`: ",
        paste(missing_model, collapse = ", "), "."
      ))
    }

    if (length(missing_proj) > 0L) {
      cli::cli_abort(paste0(
        "Variable(s) in `", arg_name, "` not found in `data_proj`: ",
        paste(missing_proj, collapse = ", "), "."
      ))
    }

    return(stats::as.formula(paste0("~", paste(variable, collapse = " + "))))
  }

  cli::cli_abort(paste0(
    "`", arg_name, "` must be a one-sided formula, character vector, or NULL."
  ))
}

#' @noRd
# Harmonize factor levels across datasets.
.harmonize_levels <- function(formula, data_model, data_proj) {
  fixed_vars <- .get_fixed_vars(formula)
  group_vars <- .get_group_vars(formula)

  # Fixed-effect levels in data_proj must already exist in data_model.
  for (v in fixed_vars) {
    if (is.factor(data_model[[v]]) || is.character(data_model[[v]])) {
      train_levels <- unique(as.character(data_model[[v]]))
      new_levels <- setdiff(unique(as.character(data_proj[[v]])), train_levels)
      new_levels <- new_levels[!is.na(new_levels)]

      if (length(new_levels) > 0L) {
        cli::cli_abort(paste0(
          "Variable `", v, "` has level(s) in data_proj not found in data_model: ",
          paste(utils::head(new_levels, 10L), collapse = ", "), "."
        ))
      }

      data_model[[v]] <- factor(as.character(data_model[[v]]), levels = train_levels)
      data_proj[[v]] <- factor(as.character(data_proj[[v]]), levels = train_levels)
    }
  }

  # New grouping levels in data_proj are allowed.
  for (v in group_vars) {
    model_levels <- unique(as.character(data_model[[v]]))
    proj_levels <- unique(as.character(data_proj[[v]]))
    all_levels <- union(model_levels, proj_levels)

    data_model[[v]] <- factor(as.character(data_model[[v]]), levels = model_levels)
    data_proj[[v]] <- factor(as.character(data_proj[[v]]), levels = all_levels)
  }

  list(data_model = data_model, data_proj = data_proj)
}

#' @noRd
# Fit the linear multilevel regression model.
.fit_lmer_model <- function(formula, data, control) {
  lme4::lmer(
    formula = formula,
    data = data,
    REML = TRUE,
    control = control
  )
}

#' @noRd
# Extract model diagnostics from a fitted lmerMod object.
.get_lmer_diagnostics <- function(fit) {
  vc <- lme4::VarCorr(fit)
  vc_df <- as.data.frame(vc)

  var_e <- stats::sigma(fit)^2
  singular_fit <- lme4::isSingular(fit)

  vcov_list <- lapply(vc, as.matrix)
  random_effect_groups <- names(vc)
  random_effect_dims <- vapply(vcov_list, ncol, integer(1L))
  is_random_intercept_only <- all(random_effect_dims == 1L)

  if (is_random_intercept_only) {
    var_u <- sum(vapply(vcov_list, function(x) x[1L, 1L], numeric(1L)))
    icc <- var_u / (var_u + var_e)
    icc_note <- "ICC computed for random-intercept structure."
  } else {
    icc <- NA_real_
    icc_note <- "Simple ICC is not computed for random-slope or complex random-effect structure."
  }

  convergence <- fit@optinfo$conv$lme4$messages
  convergence <- if (is.null(convergence) || length(convergence) == 0L) {
    "OK"
  } else {
    paste(convergence, collapse = "; ")
  }

  list(
    sigma = stats::sigma(fit),
    residual_variance = var_e,
    random_effects = vc_df,
    random_effect_groups = random_effect_groups,
    random_effect_dims = random_effect_dims,
    is_random_intercept_only = is_random_intercept_only,
    icc = icc,
    icc_note = icc_note,
    singular_fit = singular_fit,
    convergence = convergence,
    nobs = stats::nobs(fit),
    REML = lme4::isREML(fit),
    logLik = as.numeric(stats::logLik(fit)),
    AIC = stats::AIC(fit),
    BIC = stats::BIC(fit)
  )
}

#' @noRd
# Extract fixed effects, random effects, and variance components.
.get_model_parameters <- function(fit) {
  list(
    fixed_effects = lme4::fixef(fit),
    random_effects = lme4::ranef(fit),
    variance_components = as.data.frame(lme4::VarCorr(fit)),
    residual_sd = stats::sigma(fit),
    residual_variance = stats::sigma(fit)^2
  )
}

#' @noRd
# Build a survey design object for domain aggregation.
.make_survey_design <- function(data, ids, weights, strata, ...) {
  if (is.null(ids)) ids <- ~1

  survey::svydesign(
    ids = ids,
    weights = weights,
    strata = strata,
    data = data,
    ...
  )
}

#' @noRd
# Rename survey::svyby output columns consistently.
.rename_svyby <- function(x, domain_vars, value_col, estimate_name, variance_name) {
  if (value_col %in% names(x)) {
    names(x)[names(x) == value_col] <- estimate_name
  }

  non_domain_cols <- setdiff(names(x), c(domain_vars, estimate_name))
  var_col <- grep("^(var|se2)(\\.|$)", non_domain_cols, value = TRUE)[1L]

  if (is.na(var_col) && length(non_domain_cols) == 1L) {
    var_col <- non_domain_cols
  }

  if (!is.na(var_col)) {
    names(x)[names(x) == var_col] <- variance_name
  }

  x
}

#' @noRd
# Compute optional direct design-based estimates.
.get_direct_estimator <- function(response_var, domain_formula, domain_chr, design, FUN) {

  direct <- survey::svyby(
    formula = stats::as.formula(paste0("~", response_var)),
    by = domain_formula, design = design, FUN = FUN,
    vartype = "var", na.rm = TRUE
  )

  direct <- .rename_svyby(
    direct, domain_vars = domain_chr, value_col = response_var,
    estimate_name = "estimate", variance_name = "variance"
  )

  direct$variance <- pmax(direct$variance, 0)
  direct$se <- sqrt(direct$variance)
  direct$rse <- ifelse(
    direct$estimate == 0, NA_real_,
    100 * direct$se / abs(direct$estimate)
  )

  direct[, c(domain_chr, "estimate", "variance", "se", "rse"), drop = FALSE]
}

#' @noRd
# Compute unweighted row counts by domain.
.domain_counts <- function(data, domain_chr, col_name) {
  out <- stats::aggregate(rep(1L, nrow(data)), data[domain_chr], length)
  names(out)[ncol(out)] <- col_name
  out
}

# ---- Main estimator function -------------------------------------------------

#' Small Area Estimation Using a Projection Estimator with a Linear Multilevel Regression Model
#'
#' @description
#' \code{sae_ml_linear()} implements a projection estimator for small area
#' estimation using a linear multilevel regression working model fitted with \code{\link[lme4]{lmer}}.
#'
#' The function is designed for situations where information is available from
#' two related surveys:
#'
#' \itemize{
#'   \item \strong{Model survey}: a smaller survey containing the response
#'   variable and auxiliary variables.
#'   \item \strong{Projection survey}: a larger survey containing auxiliary
#'   variables and survey design information, but not the response variable.
#' }
#'
#' The function fits a linear multilevel model using the model survey, predicts
#' the response variable for all units in the projection survey, and produces
#' domain-level projection estimates.
#'
#' @details
#' The model formula must follow \code{lme4::lmer()} syntax and must include at
#' least one random-effect term.
#'
#' A random-intercept model can be specified as:
#'
#' \preformatted{
#' Y ~ X1 + X2 + Z1 + Z2 + (1 | kab_kota)
#' }
#'
#' where \code{Y} is the response variable, \code{X1} and \code{X2} are
#' unit-level auxiliary variables, \code{Z1} and \code{Z2} are auxiliary
#' variables, and \code{kab_kota} identifies the random-intercept grouping
#' level.
#'
#' The estimation procedure consists of three main steps:
#'
#' \enumerate{
#'   \item A linear multilevel model is fitted to \code{data_model} using
#'   \code{lme4::lmer()} with restricted maximum likelihood estimation.
#'
#'   \item Unit-level predictions are generated for all observations in
#'   \code{data_proj}. Predictions use \code{re.form = NULL} and
#'   \code{allow.new.levels = TRUE}.
#'
#'   \item Predicted values are aggregated by domain to obtain synthetic
#'   estimates. A design-based residual correction calculated from
#'   \code{data_model} is added to obtain the final projection estimate.
#' }
#'
#' For domain \eqn{d}, the final projection estimator is:
#'
#' \deqn{
#' \hat{Y}^{PR}*{d} =
#' \hat{Y}^{SYN}*{d} +
#' \hat{B}*{d}
#' }
#'
#' where \eqn{\hat{Y}^{SYN}*{d}} is the synthetic estimate obtained from the
#' projection survey and \eqn{\hat{B}*{d}} is the design-based residual
#' correction obtained from the model survey.
#'
#' The plug-in variance estimator is:
#'
#' \deqn{
#' \widehat{Var}(\hat{Y}^{PR}*{d}) =
#' \widehat{Var}(\hat{Y}^{SYN}*{d}) +
#' \widehat{Var}(\hat{B}*{d})
#' }
#'
#' The plug-in variance does not account for uncertainty in the estimated
#' multilevel model parameters.
#'
#' Fixed-effect predictors in \code{data_proj} must have levels that already
#' exist in \code{data_model}. New factor levels for fixed-effect predictors
#' produce an error. In contrast, new grouping levels for random effects are
#' allowed and receive a random-effect contribution of zero.
#'
#' Fixed-effect predictors with zero variance in \code{data_model} are removed
#' automatically before model fitting. A note identifying removed predictors is
#' stored in the returned object.
#'
#' If a domain appears in \code{data_proj} but has no observations in
#' \code{data_model}, the residual correction and its variance are set to zero.
#' The final estimate for that domain is therefore equal to its synthetic
#' estimate.
#'
#' Survey design variables, including cluster identifiers, strata, and sampling
#' weights, are used for domain-level aggregation and residual correction
#' through the \pkg{survey} package. They are not used as fitting weights in
#' \code{lme4::lmer()}.
#'
#' @param formula A two-sided linear multilevel model formula written using
#' \code{lme4::lmer()} syntax. The formula must include at least one
#' random-effect term. For example,
#' \code{Y ~ X1 + X2 + (1 | kab_kota)}.
#'
#' @param data_model A data frame containing model-survey data. It must contain
#' the response variable, all fixed-effect variables, random-effect grouping
#' variables, domain variables, and survey design variables used in estimation.
#'
#' @param data_proj A data frame containing projection-survey data. It must
#' contain all fixed-effect variables, random-effect grouping variables, domain
#' variables, and survey design variables. The response variable is not required
#' in \code{data_proj}.
#'
#' @param domain A character vector or one-sided formula identifying the domain
#' variable or variables used for domain-level aggregation. For example,
#' \code{"kab_kota"}, \code{c("prov", "kab_kota")}, or
#' \code{~prov + kab_kota}.
#'
#' @param cluster_ids A character vector or one-sided formula identifying
#' cluster or primary sampling unit variables used in the survey design. Use
#' \code{~1} when clustering is not used. The default is \code{~1}.
#'
#' @param weight A character string or one-sided formula identifying the survey
#' weight variable. The specification must identify exactly one variable. Use
#' \code{NULL} when sampling weights are not available. The default is
#' \code{NULL}.
#'
#' @param strata A character string or one-sided formula identifying the
#' stratification variable used in the survey design. Use \code{NULL} when
#' stratification is not used. The default is \code{NULL}.
#'
#' @param summary_function Character string specifying the domain-level
#' estimand. Available options are \code{"mean"} for domain means and
#' \code{"total"} for domain totals. The default is \code{"mean"}.
#'
#' @param keep_unit Logical. If \code{TRUE}, the returned object includes
#' unit-level predictions for \code{data_proj}, and unit-level fitted values and
#' model residuals for \code{data_model}. The default is \code{FALSE}.
#'
#' @param seed Integer value used to set the random seed before model fitting.
#' The default is \code{1}.
#'
#' @param control A control object created by \code{lme4::lmerControl()} and
#' passed to \code{lme4::lmer()}. The default uses the \code{"bobyqa"}
#' optimizer with \code{maxfun = 2e5}.
#'
#' @param return_direct Logical. If \code{TRUE}, direct design-based estimates
#' are calculated from \code{data_model} and included in the returned object.
#' The default is \code{FALSE}.
#'
#' @param ... Additional arguments passed to \code{survey::svydesign()}. For
#' example, \code{nest = TRUE} can be supplied for nested cluster designs.
#'
#' @return
#' An object of class \code{"sae_ml_linear"} containing:
#'
#' \describe{
#'   \item{call}{
#'   The matched function call.
#'   }
#'
#'   \item{formula}{
#'   The final model formula used for estimation after removal of any
#'   zero-variance fixed-effect predictors.
#'   }
#'
#'   \item{estimator}{
#'   A character string identifying the estimator as \code{"bias_corrected"}.
#'   }
#'
#'   \item{fitted_model}{
#'   The fitted \code{lmerMod} object returned by \code{lme4::lmer()}.
#'   }
#'
#'   \item{model_parameters}{
#'   A list containing fixed effects, random effects, variance components,
#'   residual standard deviation, and residual variance.
#'   }
#'
#'   \item{estimates}{
#'   A data frame containing one row for each domain. It includes domain
#'   identifiers, \code{estimate}, \code{variance}, \code{se}, and
#'   \code{rse}.
#'   }
#'
#'   \item{estimation_details}{
#'   A data frame containing domain identifiers,
#'   \code{estimate_synthetic}, \code{variance_synthetic},
#'   \code{correction}, \code{variance_correction},
#'   \code{estimate_final}, \code{variance_final}, \code{se_final},
#'   \code{rse_final}, \code{n_model}, and \code{n_proj}.
#'   }
#'
#'   \item{diagnostics}{
#'   A list containing model diagnostics, including residual standard
#'   deviation, residual variance, random-effect variance components,
#'   intraclass correlation coefficient where applicable, singularity status,
#'   convergence information, number of observations, REML status,
#'   log-likelihood, AIC, and BIC.
#'   }
#'
#'   \item{notes}{
#'   A character vector containing notes generated during estimation, such as
#'   removed zero-variance predictors, singular fits, convergence issues,
#'   negative variances clamped to zero, or domains without residual correction.
#'   }
#'
#'   \item{unit_projection}{
#'   Returned only when \code{keep_unit = TRUE}. A data frame containing
#'   \code{data_proj} with an additional \code{.prediction} column.
#'   }
#'
#'   \item{unit_model_residual}{
#'   Returned only when \code{keep_unit = TRUE}. A data frame containing
#'   \code{data_model} with additional \code{.fitted_model} and
#'   \code{.model_residual} columns.
#'   }
#'
#'   \item{direct_estimator}{
#'   Returned only when \code{return_direct = TRUE}. A data frame containing
#'   direct design-based estimates, variances, and relative standard errors for
#'   each domain.
#'   }
#' }
#'
#' @references
#' Kim, J. K., & Rao, J. N. K. (2012). Combining data from two independent
#' surveys: A model-assisted approach. \emph{Biometrika}, 99(1), 85--100.
#'
#' Moura, F. A. S., & Holt, D. (1999). Small area estimation using multilevel
#' models. \emph{Survey Methodology}, 25(1), 73--80.
#'
#' Bates, D., Maechler, M., Bolker, B., & Walker, S. (2015). Fitting linear
#' mixed-effects models using lme4. \emph{Journal of Statistical Software},
#' 67(1), 1--48.
#'
#' Food and Agriculture Organization of the United Nations. (2021).
#' \emph{Guidelines on data disaggregation for SDG indicators using survey
#' data} (1st ed.). https://doi.org/10.4060/cb3253en
#'
#' Finch, W. H., Bolin, J. E., & Kelley, K. (2014).
#' \emph{Multilevel Modeling Using R}. CRC Press.
#'
#' Hox, J. J., Moerbeek, M., & van de Schoot, R. (2018).
#' \emph{Multilevel Analysis: Techniques and Applications} (3rd ed.).
#' Routledge.
#'
#' @seealso
#' \code{\link[lme4]{lmer}} for fitting linear multilevel models,
#' \code{\link[survey]{svydesign}} for survey design specification, and
#' \code{\link[lme4]{isSingular}} for diagnosing singular fits.
#'
#' @examples
#' data("saeml_modelsvy", package = "saeproj.multilevel")
#' data("saeml_projsvy", package = "saeproj.multilevel")
#'
#' result <- sae_ml_linear(
#'   formula = Y ~ X1 + X2 + X3 + X4 + Z1 + Z2 + (1 | kab_kota),
#'   data_model = saeml_modelsvy,
#'   data_proj = saeml_projsvy,
#'   domain = "kab_kota",
#'   cluster_ids = ~1,
#'   weight = "WEIND",
#'   strata = "kab_kota",
#'   summary_function = "mean"
#' )
#'
#' summary(result)
#'
#' @importFrom survey svydesign svyby svymean svytotal
#' @importFrom cli cli_abort cli_warn
#' @importFrom dplyr left_join
#' @importFrom lme4 lmer lmerControl VarCorr fixef ranef isSingular isREML
#' @importFrom reformulas findbars nobars
#' @importFrom stats AIC BIC aggregate as.formula logLik model.frame predict sigma update
#' @importFrom utils head globalVariables
#'
#' @export
sae_ml_linear <- function(
    formula,
    data_model,
    data_proj,
    domain,
    cluster_ids = ~1,
    weight = NULL,
    strata = NULL,
    summary_function = "mean",
    keep_unit = FALSE,
    seed = 1L,
    control = lme4::lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 2e5)
    ),
    return_direct = FALSE,
    ...
) {
  mc <- match.call()
  notes <- character(0L)

  # -- 0. Basic type checks ----------------------------------------------------
  if (!is.data.frame(data_model)) {
    cli::cli_abort("`data_model` must be a data.frame.")
  }

  if (!is.data.frame(data_proj)) {
    cli::cli_abort("`data_proj` must be a data.frame.")
  }

  if (!inherits(formula, "formula") || length(formula) != 3L) {
    cli::cli_abort("`formula` must be a two-sided formula, e.g. `y ~ x + (1 | area)`.")
  }

  # -- 1. Input validation -----------------------------------------------------
  cluster_ids <- .check_variable(cluster_ids, data_model, data_proj, "cluster_ids")
  weight <- .check_variable(weight, data_model, data_proj, "weight")
  strata <- .check_variable(strata, data_model, data_proj, "strata")
  domain_formula <- .check_variable(domain, data_model, data_proj, "domain")

  if (!is.null(weight) && length(all.vars(weight)) != 1L) {
    cli::cli_abort("`weight` must identify exactly one survey weight variable.")
  }

  domain_chr <- all.vars(domain_formula)
  response_var <- .get_response_var(formula)
  fixed_vars <- .get_fixed_vars(formula)
  group_vars <- .get_group_vars(formula)

  if (length(domain_chr) == 0L) {
    cli::cli_abort("`domain` must identify at least one domain variable.")
  }

  if (length(group_vars) == 0L) {
    cli::cli_abort("`formula` must contain at least one random effect, e.g. `(1 | area)`.")
  }

  summary_function <- match.arg(summary_function, c("mean", "total"))
  FUN <- switch(summary_function, mean = survey::svymean, total = survey::svytotal)

  # -- 2. Column presence & missing value checks -------------------------------
  req_model <- unique(c(
    response_var, fixed_vars, group_vars, domain_chr,
    all.vars(cluster_ids), all.vars(weight), all.vars(strata)
  ))

  req_proj <- setdiff(req_model, response_var)

  .check_required_columns(data_model, req_model, "data_model")
  .check_required_columns(data_proj, req_proj, "data_proj")
  .check_missing_values(data_model, req_model, "data_model")
  .check_missing_values(data_proj, req_proj, "data_proj")

  # -- 3. Level harmonization --------------------------------------------------
  harmonized <- .harmonize_levels(formula, data_model, data_proj)
  data_model <- harmonized$data_model
  data_proj <- harmonized$data_proj

  y <- data_model[[response_var]]

  # -- 4. Zero-variance predictor removal -------------------------------------
  zv_vars <- Filter(function(v) length(unique(data_model[[v]])) == 1L, fixed_vars)

  if (length(zv_vars) > 0L) {
    cli::cli_warn(paste0(
      "Removing zero-variance predictor(s): ",
      paste(zv_vars, collapse = ", "), "."
    ))

    notes <- c(notes, paste0(
      "Removed zero-variance predictor(s): ",
      paste(zv_vars, collapse = ", "), "."
    ))

    formula <- stats::update(
      formula,
      stats::as.formula(paste("~ . -", paste(zv_vars, collapse = " - ")))
    )

    fixed_vars <- setdiff(fixed_vars, zv_vars)
  }

  # -- 5. Model fitting --------------------------------------------------------
  set.seed(seed)
  fit <- .fit_lmer_model(formula, data_model, control)

  diagnostics <- .get_lmer_diagnostics(fit)
  model_parameters <- .get_model_parameters(fit)

  if (isTRUE(diagnostics$singular_fit)) {
    cli::cli_warn("Singular fit detected.")
    notes <- c(notes, "Singular fit detected.")
  }

  if (!identical(diagnostics$convergence, "OK")) {
    cli::cli_warn(paste0("Convergence issue: ", diagnostics$convergence))
    notes <- c(notes, paste0("Convergence issue: ", diagnostics$convergence))
  }

  if (!diagnostics$is_random_intercept_only) {
    notes <- c(notes, diagnostics$icc_note)
  }

  # -- 6. Prediction -----------------------------------------------------------
  data_proj$.prediction <- stats::predict(
    fit,
    newdata = data_proj,
    re.form = NULL,
    allow.new.levels = TRUE
  )

  data_model$.fitted_model <- stats::predict(
    fit,
    newdata = data_model,
    re.form = NULL,
    allow.new.levels = TRUE
  )

  data_model$.model_residual <- y - data_model$.fitted_model

  # -- 7. Survey designs -------------------------------------------------------
  svy_model <- .make_survey_design(data_model, cluster_ids, weight, strata, ...)
  svy_proj <- .make_survey_design(data_proj, cluster_ids, weight, strata, ...)

  # -- 8. Synthetic estimate & bias correction ---------------------------------
  est_bias <- survey::svyby(
    formula = ~.model_residual,
    by = domain_formula,
    design = svy_model,
    FUN = FUN,
    vartype = "var",
    na.rm = TRUE
  )

  est_bias <- .rename_svyby(
    est_bias,
    domain_vars = domain_chr,
    value_col = ".model_residual",
    estimate_name = "correction",
    variance_name = "variance_correction"
  )

  est_ypr <- survey::svyby(
    formula = ~.prediction,
    by = domain_formula,
    design = svy_proj,
    FUN = FUN,
    vartype = "var",
    na.rm = TRUE
  )

  est_ypr <- .rename_svyby(
    est_ypr,
    domain_vars = domain_chr,
    value_col = ".prediction",
    estimate_name = "estimate_synthetic",
    variance_name = "variance_synthetic"
  )

  # -- 9. Combine & final estimates --------------------------------------------
  df_result <- dplyr::left_join(est_ypr, est_bias, by = domain_chr)

  no_corr_idx <- is.na(df_result$correction)

  if (any(no_corr_idx)) {
    n_no_corr <- sum(no_corr_idx)

    notes <- c(notes, paste0(
      n_no_corr,
      " out-of-sample domain(s): correction set to zero."
    ))

    df_result$correction[no_corr_idx] <- 0
    df_result$variance_correction[no_corr_idx] <- 0
  }

  estimate_final <- df_result$estimate_synthetic + df_result$correction
  variance_raw <- df_result$variance_synthetic + df_result$variance_correction

  neg_var_idx <- which(!is.na(variance_raw) & variance_raw < 0)

  if (length(neg_var_idx) > 0L) {
    cli::cli_warn("Negative plug-in variance detected; clamped to zero.")

    notes <- c(notes, paste0(
      length(neg_var_idx),
      " domain(s): negative variance clamped to zero."
    ))
  }

  variance_final <- pmax(variance_raw, 0)
  se_final <- sqrt(variance_final)

  rse_final <- ifelse(
    estimate_final == 0,
    NA_real_,
    100 * se_final / abs(estimate_final)
  )

  estimates <- df_result[, domain_chr, drop = FALSE]
  estimates$estimate <- estimate_final
  estimates$variance <- variance_final
  estimates$se <- se_final
  estimates$rse <- rse_final

  n_model <- .domain_counts(data_model, domain_chr, "n_model")
  n_proj <- .domain_counts(data_proj, domain_chr, "n_proj")

  estimation_details <- df_result
  estimation_details$estimate_final <- estimate_final
  estimation_details$variance_final <- variance_final
  estimation_details$se_final <- se_final
  estimation_details$rse_final <- rse_final

  estimation_details <- dplyr::left_join(estimation_details, n_model, by = domain_chr)
  estimation_details <- dplyr::left_join(estimation_details, n_proj, by = domain_chr)

  estimation_details$n_model[is.na(estimation_details$n_model)] <- 0L
  estimation_details$n_proj[is.na(estimation_details$n_proj)] <- 0L

  # -- 10. Output assembly -----------------------------------------------------
  out <- list(
    call = mc,
    formula = formula,
    estimator = "bias_corrected",
    fitted_model = fit,
    model_parameters = model_parameters,
    estimates = estimates,
    estimation_details = estimation_details,
    diagnostics = diagnostics,
    notes = unique(notes)
  )

  if (keep_unit) {
    out$unit_projection <- data_proj
    out$unit_model_residual <- data_model
  }

  if (return_direct) {
    out$direct_estimator <- .get_direct_estimator(
      response_var = response_var,
      domain_formula = domain_formula,
      domain_chr = domain_chr,
      design = svy_model,
      FUN = FUN
    )
  }

  structure(out, class = "sae_ml_linear")
}

# ---- S3 methods --------------------------------------------------------------

#' Print method for sae_ml_linear
#'
#' @param x Object of class \code{"sae_ml_linear"}.
#' @param n Number of rows to display.
#' @param ... Further arguments.
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
#' @method print sae_ml_linear
print.sae_ml_linear <- function(x, n = 6L, ...) {
  cat("SAE Projection Estimator using Linear Multilevel Model\n")
  cat("-------------------------------------------------------\n")
  cat("Formula   :", deparse(x$formula), "\n")
  cat("Estimator :", x$estimator, "\n")
  cat("Domains   :", nrow(x$estimates), "\n\n")

  cat("Estimates:\n")
  print(utils::head(x$estimates, n), row.names = FALSE)

  invisible(x)
}

#' Summary method for sae_ml_linear
#'
#' @param object Object of class \code{"sae_ml_linear"}.
#' @param n Number of rows to display.
#' @param ... Further arguments.
#'
#' @return Invisibly returns \code{object}.
#'
#' @export
#' @method summary sae_ml_linear
summary.sae_ml_linear <- function(object, n = 6L, ...) {
  cat("SAE Projection Estimator using Linear Multilevel Model\n")
  cat("-------------------------------------------------------\n")
  cat("Formula   :", deparse(object$formula), "\n")
  cat("Estimator :", object$estimator, "\n")
  cat("Domains   :", nrow(object$estimates), "\n\n")

  cat("Model diagnostics:\n")
  cat("  nobs        :", object$diagnostics$nobs, "\n")
  cat("  sigma       :", round(object$diagnostics$sigma, 4L), "\n")
  cat("  ICC         :",
      if (is.na(object$diagnostics$icc)) "NA" else round(object$diagnostics$icc, 4L),
      "\n")
  cat("  singular    :", object$diagnostics$singular_fit, "\n")
  cat("  convergence :", object$diagnostics$convergence, "\n\n")

  cat("Estimates:\n")
  print(utils::head(object$estimates, n), row.names = FALSE)

  invisible(object)
}

#' Coerce an sae_ml_linear object to a data frame
#'
#' @param x Object of class \code{"sae_ml_linear"}.
#' @param row.names Ignored.
#' @param optional Ignored.
#' @param ... Further arguments.
#'
#' @return A data frame containing domain-level estimates.
#'
#' @export
#' @method as.data.frame sae_ml_linear
as.data.frame.sae_ml_linear <- function(x, row.names = NULL, optional = FALSE, ...) {
  as.data.frame(x$estimates)
}
