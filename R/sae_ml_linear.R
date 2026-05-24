# R/sae_ml_linear.R
# Small Area Estimation via Projection Estimator using Linear Multilevel Working Model
#
# References:
#   Kim & Rao (2012). Biometrika, 99(1), 85-100.
#   Moura & Holt (1999). Survey Methodology, 25(1), 73-80.
#   Bates et al. (2015). Journal of Statistical Software, 67(1), 1-48.
#   Lumley (2010). Complex Surveys. Wiley.

# ---- Setup ------------------------------------------------------------------

# Suppress R CMD check NOTE for internal data.frame column variables.
utils::globalVariables(c(
  "bias", "var_bias", "ypr", "var_ypr", "estimate",
  "variance", "se", "rse", ".y_hat", ".resid", ".y_hat_model"
))

# ---- Internal helpers --------------------------------------------------------

#' @noRd
# Extracts the response variable from a mixed-effects formula.
.get_response_var <- function(formula) {
  all.vars(reformulas::nobars(formula))[1L]
}

#' @noRd
# Extracts fixed-effect variables from a mixed-effects formula.
.get_fixed_vars <- function(formula) {
  setdiff(all.vars(reformulas::nobars(formula)), .get_response_var(formula))
}

#' @noRd
# Extracts grouping variables from random-effect terms.
.get_group_vars <- function(formula) {
  bars <- reformulas::findbars(formula)
  if (length(bars) == 0L) return(character(0L))
  unique(unlist(lapply(bars, function(x) all.vars(x[[3L]]))))
}

#' @noRd
# Checks required variables for missing values.
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
# Checks that required variables exist in a data frame.
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
# Validates a design/domain variable and converts character input to formula.
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
# Harmonizes categorical fixed-effect and grouping levels across datasets.
.harmonize_levels <- function(formula, data_model, data_proj) {
  fixed_vars <- .get_fixed_vars(formula)
  group_vars <- .get_group_vars(formula)

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
# Fits the linear mixed-effects working model.
# Survey weights are used in the design-based aggregation step, not as lmer weights.
.fit_lmer_model <- function(formula, data, control) {
  lme4::lmer(
    formula = formula,
    data = data,
    REML = FALSE,
    control = control
  )
}

#' @noRd
# Extracts model diagnostics and determines prediction mode.
.get_lmer_diagnostics <- function(fit) {
  vc <- lme4::VarCorr(fit)
  var_u <- sum(unlist(lapply(vc, function(x) diag(as.matrix(x)))))
  var_e <- stats::sigma(fit)^2
  icc <- var_u / (var_u + var_e)
  singular_fit <- lme4::isSingular(fit)

  prediction_mode <- if (!is.na(icc) && icc >= 0.05 && !singular_fit) {
    "conditional"
  } else {
    "fixed_only"
  }

  convergence <- fit@optinfo$conv$lme4$messages
  convergence <- if (is.null(convergence) || length(convergence) == 0L) {
    "OK"
  } else {
    paste(convergence, collapse = "; ")
  }

  list(
    icc = icc,
    singular_fit = singular_fit,
    prediction_mode = prediction_mode,
    convergence = convergence,
    sigma = stats::sigma(fit),
    nobs = stats::nobs(fit),
    AIC = stats::AIC(fit),
    BIC = stats::BIC(fit)
  )
}

#' @noRd
# Builds a survey design object for domain aggregation.
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
# Renames survey::svyby output columns consistently.
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
# Computes optional direct design-based estimates.
.get_direct_estimator <- function(response_var, domain_formula, domain_chr, design, FUN) {
  direct <- survey::svyby(
    formula = stats::as.formula(paste0("~", response_var)),
    by = domain_formula,
    design = design,
    FUN = FUN,
    vartype = c("var", "cvpct"),
    na.rm = TRUE
  )

  non_domain_cols <- setdiff(names(direct), domain_chr)
  var_col <- grep("^(var|se2)", non_domain_cols, value = TRUE, ignore.case = TRUE)[1L]
  cv_col <- grep("^cv", non_domain_cols, value = TRUE, ignore.case = TRUE)[1L]
  est_col <- setdiff(non_domain_cols, c(var_col, cv_col))[1L]

  if (!is.na(est_col)) names(direct)[names(direct) == est_col] <- "estimate"
  if (!is.na(var_col)) names(direct)[names(direct) == var_col] <- "variance"
  if (!is.na(cv_col)) names(direct)[names(direct) == cv_col] <- "rse"

  direct
}

#' @noRd
# Computes unweighted row counts by domain.
.domain_counts <- function(data, domain_chr, col_name) {
  out <- stats::aggregate(rep(1L, nrow(data)), data[domain_chr], length)
  names(out)[ncol(out)] <- col_name
  out
}

# ---- Main estimator function -------------------------------------------------

#' Small Area Estimation via Projection Estimator with Linear Multilevel Model
#'
#' @description
#' Implements the projection estimator for Small Area Estimation (SAE) using a
#' linear mixed-effects working model fitted with \code{\link[lme4]{lmer}}.
#'
#' The function fits the working model on \code{data_model}, generates synthetic
#' predictions on \code{data_proj}, and aggregates predictions to domain-level
#' means or totals using survey design information.
#'
#' @details
#' The estimator combines synthetic prediction from a multilevel model,
#' design-based residual bias correction, adaptive prediction using ICC
#' diagnostics, and survey-weighted domain aggregation.
#'
#' Survey design variables, including \code{cluster_ids}, \code{weight}, and
#' \code{strata}, are used in the design-based aggregation and residual
#' correction steps through \code{\link[survey]{svydesign}}. The linear
#' mixed-effects working model is fitted with \code{\link[lme4]{lmer}} and does
#' not treat survey weights as full survey-design weights.
#'
#' Conditional prediction is used when ICC is at least \code{0.05} and the fitted
#' model is not singular. Otherwise, fixed-effect-only prediction is used
#' through \code{re.form = NA}. New grouping levels in \code{data_proj} are
#' allowed during prediction.
#'
#' The plug-in variance for the bias-corrected estimator is computed as
#' \code{Var(final estimate) = Var(synthetic estimate) + Var(residual correction)}.
#' This approximation assumes \code{data_model} and \code{data_proj} are
#' independent, or treated as independent, and does not fully account for
#' mixed-model parameter uncertainty.
#'
#' @param formula An \code{lme4::lmer()}-style formula, for example
#'   \code{y ~ x1 + x2 + (1 | area)}.
#' @param data_model Data frame for the smaller model survey. Must contain the
#'   response, all predictors, grouping variables, domain variable(s), and survey
#'   design variables.
#' @param data_proj Data frame for the larger projection survey. Must contain all
#'   predictors, domain variable(s), and survey design variables. The response
#'   variable is not required.
#' @param domain Domain variable name(s): a character scalar, character vector,
#'   or a one-sided formula.
#' @param cluster_ids Cluster or PSU variable for survey design. Character,
#'   formula, or \code{~1}.
#' @param weight Survey weight variable. Character scalar, one-sided formula, or
#'   \code{NULL}. The variable name must exist in both \code{data_model} and
#'   \code{data_proj}; the values may differ between datasets.
#' @param strata Stratification variable. Character, formula, or \code{NULL}.
#' @param summary_function Aggregation function for the domain level:
#'   \code{"mean"} or \code{"total"}.
#' @param keep_unit Logical. If \code{TRUE}, unit-level predictions and model
#'   residuals are returned.
#' @param seed Integer seed for reproducibility.
#' @param control Control object passed to \code{\link[lme4]{lmerControl}}.
#' @param return_direct Logical. If \code{TRUE}, returns direct survey estimators
#'   based on \code{data_model}.
#' @param ... Additional arguments passed to \code{\link[survey]{svydesign}}.
#'
#' @return An object of class \code{"sae_ml_linear"}, a list with:
#' \describe{
#'   \item{\code{call}}{The matched call.}
#'   \item{\code{formula}}{The model formula.}
#'   \item{\code{estimator}}{The estimator type used; always
#'   \code{"bias_corrected"}.}
#'   \item{\code{fitted_model}}{The fitted \code{lmerMod} object.}
#'   \item{\code{estimates}}{Domain-level estimates with columns for domain
#'   variable(s), \code{estimate}, \code{variance}, \code{se}, and \code{rse}.}
#'   \item{\code{estimation_details}}{Extended domain-level components including
#'   synthetic estimates, bias corrections, and sample counts.}
#'   \item{\code{diagnostics}}{Model diagnostics including ICC, singular-fit
#'   status, prediction mode, convergence message, sigma, AIC, and BIC.}
#'   \item{\code{notes}}{Concise run-specific notes and warning conditions.}
#'   \item{\code{unit_projection}}{Returned when \code{keep_unit = TRUE}.}
#'   \item{\code{unit_model_residual}}{Returned when \code{keep_unit = TRUE}.}
#'   \item{\code{direct_estimator}}{Returned when
#'   \code{return_direct = TRUE}.}
#' }
#'
#' @references
#' Kim, J. K. and Rao, J. N. K. (2012). Combining data from two independent
#' surveys: a model-assisted approach. \emph{Biometrika}, 99(1), 85--100.
#'
#' Moura, F. A. S. and Holt, D. (1999). Small area estimation using multilevel
#' models. \emph{Survey Methodology}, 25(1), 73--80.
#'
#' Bates, D., Maechler, M., Bolker, B. and Walker, S. (2015). Fitting linear
#' mixed-effects models using lme4. \emph{Journal of Statistical Software},
#' 67(1), 1--48.
#'
#' Lumley, T. (2010). \emph{Complex Surveys: A Guide to Analysis Using R}. Wiley.
#'
#' @examples
#' \dontrun{
#' result <- sae_ml_linear(
#'   formula = income ~ educ + age + (1 | kabkot),
#'   data_model = survey_model,
#'   data_proj = survey_proj,
#'   domain = "kabkot",
#'   weight = "w",
#'   summary_function = "mean",
#'   nest = TRUE
#' )
#'
#' print(result)
#' summary(result)
#' as.data.frame(result)
#' result$estimation_details
#' result$notes
#' }
#'
#' @importFrom survey svydesign svyby svymean svytotal
#' @importFrom cli cli_abort cli_warn
#' @importFrom dplyr left_join
#' @importFrom reformulas findbars nobars
#' @importFrom stats AIC BIC aggregate as.formula model.frame predict sd sigma update
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
    seed = 1,
    control = lme4::lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 2e5)
    ),
    return_direct = FALSE,
    ...
) {
  mc <- match.call()
  estimator_type <- "bias_corrected"
  notes <- character(0L)

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
  if (length(domain_chr) == 0L) {
    cli::cli_abort("`domain` must identify at least one domain variable.")
  }
  response_var <- .get_response_var(formula)
  fixed_vars <- .get_fixed_vars(formula)
  group_vars <- .get_group_vars(formula)

  if (length(group_vars) == 0L) {
    cli::cli_abort("`formula` must contain at least one random effect.")
  }

  y <- data_model[[response_var]]

  summary_function <- match.arg(summary_function, c("mean", "total"))
  FUN <- switch(
    summary_function,
    mean = survey::svymean,
    total = survey::svytotal
  )

  # -- 2. Missing values & level harmonization --------------------------------
  req_model <- unique(c(
    response_var, fixed_vars, group_vars, domain_chr,
    all.vars(cluster_ids), all.vars(weight), all.vars(strata)
  ))
  req_proj <- setdiff(req_model, response_var)

  .check_required_columns(data_model, req_model, "data_model")
  .check_required_columns(data_proj, req_proj, "data_proj")

  .check_missing_values(data_model, req_model, "data_model")
  .check_missing_values(data_proj, req_proj, "data_proj")

  harmonized <- .harmonize_levels(formula, data_model, data_proj)
  data_model <- harmonized$data_model
  data_proj <- harmonized$data_proj

  # -- 3. Zero-variance predictors ---------------------------------------------
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

  # -- 4. Model fitting & diagnostics -----------------------------------------
  set.seed(seed)
  fit <- .fit_lmer_model(formula, data_model, control)
  diagnostics <- .get_lmer_diagnostics(fit)

  re_form_use <- if (identical(diagnostics$prediction_mode, "conditional")) NULL else NA

  if (isTRUE(diagnostics$singular_fit)) {
    cli::cli_warn("Singular fit detected. Consider simplifying the random-effect structure.")
    notes <- c(notes, "Singular fit detected.")
  }

  if (!identical(diagnostics$convergence, "OK")) {
    cli::cli_warn(paste0("Convergence issue detected: ", diagnostics$convergence))
    notes <- c(notes, paste0("Convergence issue: ", diagnostics$convergence))
  }

  if (!is.na(diagnostics$icc) && diagnostics$icc < 0.05) {
    cli::cli_warn(paste0(
      "Low ICC detected (ICC = ", round(diagnostics$icc, 4L),
      "). Fixed-effect-only prediction is used."
    ))
    notes <- c(notes, paste0("Low ICC detected: ", round(diagnostics$icc, 4L), "."))
  }

  if (identical(diagnostics$prediction_mode, "fixed_only")) {
    notes <- c(notes, "Fixed-effect-only prediction used.")
  }

  # -- 5. Prediction & survey design ------------------------------------------
  data_proj$.y_hat <- stats::predict(
    fit,
    newdata = data_proj,
    re.form = re_form_use,
    allow.new.levels = TRUE
  )

  data_model$.y_hat_model <- stats::predict(
    fit,
    newdata = data_model,
    re.form = re_form_use,
    allow.new.levels = TRUE
  )

  data_model$.resid <- y - data_model$.y_hat_model

  svy_model <- .make_survey_design(data_model, cluster_ids, weight, strata, ...)
  svy_proj <- .make_survey_design(data_proj, cluster_ids, weight, strata, ...)

  # -- 6. Synthetic estimator & bias correction --------------------------------
  est_bias <- survey::svyby(
    formula = ~.resid,
    by = domain_formula,
    design = svy_model,
    FUN = FUN,
    vartype = "var",
    na.rm = TRUE
  )
  est_bias <- .rename_svyby(
    est_bias,
    domain_vars = domain_chr,
    value_col = ".resid",
    estimate_name = "correction",
    variance_name = "variance_correction"
  )

  est_ypr <- survey::svyby(
    formula = ~.y_hat,
    by = domain_formula,
    design = svy_proj,
    FUN = FUN,
    vartype = "var",
    na.rm = TRUE
  )
  est_ypr <- .rename_svyby(
    est_ypr,
    domain_vars = domain_chr,
    value_col = ".y_hat",
    estimate_name = "estimate_synthetic",
    variance_name = "variance_synthetic"
  )

  # -- 7. Final estimates ------------------------------------------------------
  df_result <- dplyr::left_join(est_ypr, est_bias, by = domain_chr)

  no_corr <- df_result[is.na(df_result$correction), domain_chr, drop = FALSE]

  if (nrow(no_corr) > 0L) {
    cli::cli_warn(paste0(
      nrow(no_corr),
      " domain(s) have no residual correction from data_model; correction is set to zero."
    ))

    notes <- c(notes, paste0(
      nrow(no_corr),
      " domain(s) have no residual correction; correction set to zero."
    ))
  }

  df_result$correction[is.na(df_result$correction)] <- 0
  df_result$variance_correction[is.na(df_result$variance_correction)] <- 0

  estimate_final <- df_result$estimate_synthetic + df_result$correction
  variance_raw <- df_result$variance_synthetic + df_result$variance_correction
  neg_var_idx <- which(!is.na(variance_raw) & variance_raw < 0)

  if (length(neg_var_idx) > 0L) {
    cli::cli_warn(paste0(
      length(neg_var_idx),
      " domain(s) have negative plug-in variance; values are clamped to zero."
    ))

    notes <- c(notes, paste0(
      length(neg_var_idx),
      " domain(s) have negative variance clamped to zero."
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

  # -- 8. Output assembly ------------------------------------------------------
  out <- list(
    call = mc,
    formula = formula,
    estimator = estimator_type,
    fitted_model = fit,
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
  cat("---------------------------------------------------\n")
  cat("Formula   :", deparse(x$formula), "\n")
  cat("Estimator :", x$estimator, "\n")
  cat("Domains   :", nrow(x$estimates), "\n\n")
  cat("Estimates (first", min(n, nrow(x$estimates)), "rows):\n")
  print(utils::head(x$estimates, n), row.names = FALSE)
  invisible(x)
}

#' Summary method for sae_ml_linear
#'
#' @param object Object of class \code{"sae_ml_linear"}.
#' @param ... Further arguments.
#'
#' @return Invisibly returns \code{object}.
#'
#' @export
#' @method summary sae_ml_linear
summary.sae_ml_linear <- function(object, ...) {
  cat("SAE Projection Estimator using Linear Multilevel Model\n\n")
  cat("Call:\n")
  print(object$call)
  cat("\n")
  cat("Formula   :", deparse(object$formula), "\n")
  cat("Estimator :", object$estimator, "\n\n")

  cat("Working model:\n")
  print(summary(object$fitted_model))

  cat("\nDiagnostics:\n")
  cat("  sigma       :", round(object$diagnostics$sigma, 4L), "\n")
  cat("  nobs        :", object$diagnostics$nobs, "\n")
  cat("  ICC         :", if (is.na(object$diagnostics$icc)) "NA" else round(object$diagnostics$icc, 4L), "\n")
  cat("  pred_mode   :", object$diagnostics$prediction_mode, "\n")
  cat("  singular    :", object$diagnostics$singular_fit, "\n")
  cat("  convergence :", object$diagnostics$convergence, "\n\n")

  cat("Final estimates:\n")
  print(utils::head(object$estimates), row.names = FALSE)
  cat("\n")

  if (!is.null(object$notes) && length(object$notes) > 0L) {
    cat("Notes:\n")
    for (note in object$notes) cat(" *", note, "\n")
    cat("\n")
  }

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
