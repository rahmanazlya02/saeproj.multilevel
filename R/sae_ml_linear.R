# R/sae_ml_linear.R
# Small Area Estimation via Projection Estimator — Linear Multilevel Working Model
#
# References:
#   Kim & Rao (2012). Biometrika, 99(1), 85-100.
#   Bates et al. (2015). J. Statistical Software, 67(1), 1-48.
#   Lumley (2010). Complex Surveys. Wiley.

# Suppress R CMD check NOTE for internal data.frame column variables
utils::globalVariables(c(".y_hat", ".resid", ".y_hat_model"))

# ---- Internal helpers --------------------------------------------------------

#' @noRd
.get_response_var <- function(formula) {
  all.vars(reformulas::nobars(formula))[1L]
}

#' @noRd
.get_fixed_vars <- function(formula) {
  setdiff(all.vars(reformulas::nobars(formula)), .get_response_var(formula))
}

#' @noRd
.get_group_vars <- function(formula) {
  bars <- reformulas::findbars(formula)
  if (length(bars) == 0L) return(character(0L))
  unique(unlist(lapply(bars, function(x) all.vars(x[[3L]]))))
}

#' @noRd
.to_formula <- function(x) {
  if (is.null(x)) return(NULL)
  if (inherits(x, "formula")) return(x)
  stats::as.formula(paste0("~", paste(x, collapse = " + ")))
}

#' @noRd
.var_names <- function(x) {
  if (is.null(x)) return(character(0L))
  if (inherits(x, "formula")) return(all.vars(x))
  if (is.character(x)) return(x)
  character(0L)
}

#' @noRd
.check_cols <- function(data, cols, data_name) {
  missing <- setdiff(unique(cols[nzchar(cols)]), names(data))
  if (length(missing) > 0L)
    stop("Column(s) not found in ", data_name, ": ",
         paste(missing, collapse = ", "), call. = FALSE)
  invisible(TRUE)
}

#' Check for missing values in required variables and stop with an informative message.
#' @noRd
.check_missing_values <- function(data, vars, data_name) {
  vars <- intersect(unique(vars), names(data))
  if (length(vars) == 0L) return(invisible(TRUE))

  na_counts <- vapply(vars, function(v) sum(is.na(data[[v]])), integer(1L))
  bad <- na_counts[na_counts > 0L]
  if (length(bad) == 0L) return(invisible(TRUE))

  lines <- paste0("- ", names(bad), ": ", bad, " missing value(s)")
  stop(
    "Missing values were found in required variables.\n\n",
    "Dataset: ", data_name, "\n",
    paste(lines, collapse = "\n"), "\n\n",
    "Please handle missing values before running sae_ml_linear(). ",
    "Rows are not removed automatically because this may change survey weights, ",
    "domain composition, and SAE estimates.",
    call. = FALSE
  )
}

#' Validate and harmonize categorical predictor levels between datasets.
#' @noRd
.validate_predictor_compatibility <- function(formula, data_model, data_proj) {
  fixed_vars <- .get_fixed_vars(formula)
  group_vars <- .get_group_vars(formula)

  dm <- data_model
  dp <- data_proj

  for (v in fixed_vars) {
    xm <- dm[[v]]
    xp <- dp[[v]]

    if (is.character(xm) || is.factor(xm) || is.character(xp) || is.factor(xp)) {
      lm_vals <- unique(as.character(xm[!is.na(xm)]))
      lp_vals <- unique(as.character(xp[!is.na(xp)]))
      new_lvls <- setdiff(lp_vals, lm_vals)

      if (length(new_lvls) > 0L)
        stop(
          "Variable '", v, "' has category level(s) in data_proj not found in data_model: ",
          paste(utils::head(new_lvls, 10L), collapse = ", "), ". ",
          "Please harmonize category labels before running sae_ml_linear().",
          call. = FALSE
        )

      dm[[v]] <- factor(as.character(xm), levels = lm_vals)
      dp[[v]] <- factor(as.character(xp), levels = lm_vals)

    } else if (is.numeric(xm) && !is.numeric(xp)) {
      stop("Variable '", v, "' is numeric in data_model but not in data_proj.", call. = FALSE)
    }
  }

  # Grouping variable: must exist in both; new levels are allowed (re.form = NA)
  for (v in group_vars) {
    if (!v %in% names(dp)) next
    xm <- dm[[v]]
    xp <- dp[[v]]
    if (is.character(xm) || is.factor(xm)) {
      all_lvls <- union(unique(as.character(xm[!is.na(xm)])),
                        unique(as.character(xp[!is.na(xp)])))
      dm[[v]] <- factor(as.character(xm), levels = all_lvls)
      dp[[v]] <- factor(as.character(xp), levels = all_lvls)
    }
  }

  list(data_model = dm, data_proj = dp)
}

#' @noRd
.build_svy_design <- function(data, ids, weight, strata, ...) {
  ids    <- if (is.null(ids))    stats::as.formula("~1") else .to_formula(ids)
  weight <- if (is.null(weight)) stats::as.formula("~1") else .to_formula(weight)

  args <- list(ids = ids, weights = weight, data = data, ...)
  if (!is.null(strata)) args$strata <- .to_formula(strata)
  do.call(survey::svydesign, args)
}

#' @noRd
.domain_counts <- function(data, domain_chr, col_name) {
  out <- stats::aggregate(
    rep(1L, nrow(data)), data[domain_chr], length
  )
  names(out)[ncol(out)] <- col_name
  out
}

#' Compute ICC for a random-intercept lmer model.
#'
#' For models with a single grouping factor and random intercept only, ICC is
#' computed as var_random / (var_random + var_residual).  For more complex
#' random-effect structures (random slopes, multiple grouping factors) the
#' interpretation is ambiguous and NA is returned.
#'
#' @noRd
.compute_icc <- function(fit) {
  vc      <- lme4::VarCorr(fit)
  groups  <- names(vc)

  # Only support single grouping factor with intercept-only random effect
  if (length(groups) != 1L) return(NA_real_)

  vc_grp  <- vc[[groups[1L]]]
  # Random intercept only: 1x1 matrix with one variance component
  if (!identical(dim(vc_grp), c(1L, 1L))) return(NA_real_)

  var_random  <- as.numeric(vc_grp[1L, 1L])
  var_resid   <- stats::sigma(fit)^2
  var_random / (var_random + var_resid)
}

#' Extract model diagnostics from a fitted lmerMod object.
#'
#' All fields are computable from base R and lme4 without additional
#' dependencies. AIC, BIC, and logLik are based on the REML criterion used
#' during fitting; they should not be used to compare models with different
#' fixed-effect structures (refit with REML = FALSE for that purpose).
#'
#' ICC is computed for models with a single grouping factor and random
#' intercept only. Returns NA for random-slope or multi-group structures,
#' where a single ICC value is not well-defined.
#'
#' @noRd
.extract_diagnostics <- function(fit) {
  list(
    fixed_effects        = lme4::fixef(fit),
    variance_components  = as.data.frame(lme4::VarCorr(fit)),
    random_effect_groups = names(lme4::VarCorr(fit)),
    singular_fit         = lme4::isSingular(fit),
    convergence_messages = fit@optinfo$conv$lme4$messages,
    sigma                = stats::sigma(fit),
    nobs                 = stats::nobs(fit),
    aic                  = stats::AIC(fit),
    bic                  = stats::BIC(fit),
    loglik               = as.numeric(stats::logLik(fit)),
    icc                  = .compute_icc(fit)
  )
}

# ---- Main function -----------------------------------------------------------

#' Small Area Estimation via Projection Estimator with Linear Multilevel Model
#'
#' @description
#' Implements the projection estimator for Small Area Estimation (SAE) using a
#' linear mixed-effects working model fitted with \code{\link[lme4]{lmer}}.
#'
#' The function fits the working model on \code{data_model} (a small survey with
#' response \code{y} and predictors \code{X}), generates synthetic predictions on
#' \code{data_proj} (a large survey with \code{X} but no \code{y}), and aggregates
#' predictions to domain-level means using survey design information.
#'
#' Predictions always use fixed effects only (\code{re.form = NA}), so new
#' random-effect levels in \code{data_proj} are permitted. The \code{cluster_ids},
#' \code{weight}, and \code{strata} arguments are used exclusively in
#' \code{\link[survey]{svydesign}} for the aggregation step, not in model fitting.
#'
#' @param formula An \code{lme4::lmer()}-style formula, e.g.
#'   \code{y ~ x1 + x2 + (1 | area)}.
#' @param data_model Data frame for the smaller model survey. Must contain the
#'   response, all predictors, grouping variables, domain variable(s), and any
#'   survey design variables.
#' @param data_proj Data frame for the larger projection survey. Must contain all
#'   predictors, domain variable(s), and any survey design variables.
#'   The response variable is not required.
#' @param domain Domain variable name(s): a character scalar, character vector,
#'   or a one-sided formula (e.g. \code{~prov + kabkot}).
#' @param cluster_ids Cluster/PSU variable for survey design. Character, formula,
#'   or \code{~1} for no clustering.
#' @param weight Survey weight variable. Character scalar, formula, or \code{NULL}
#'   for equal weights.
#' @param strata Stratification variable. Character, formula, or \code{NULL}.
#' @param estimator \code{"bias_corrected"} (default) adds an empirical residual
#'   correction from \code{data_model}; \code{"synthetic"} uses projection only.
#' @param keep_unit Logical. If \code{TRUE}, unit-level predictions and model
#'   residuals are returned in the output object.
#' @param control Control object passed to \code{\link[lme4]{lmer}}.
#' @param ... Additional named arguments passed to both
#'   \code{\link[lme4]{lmer}} and \code{\link[survey]{svydesign}}.
#'   Supply only arguments accepted by both, or rely on defaults for one of them.
#'
#' @return An object of class \code{"sae_ml_linear"}, a list with:
#' \describe{
#'   \item{\code{call}}{The matched call.}
#'   \item{\code{formula}}{The model formula.}
#'   \item{\code{estimator}}{The estimator type used.}
#'   \item{\code{fitted_model}}{The fitted \code{lmerMod} object.}
#'   \item{\code{estimates}}{Data frame of final domain-level estimates with
#'     columns: domain variable(s), \code{estimate}, \code{variance},
#'     \code{se}, \code{rse}.}
#'   \item{\code{estimation_details}}{Data frame of estimation components for
#'     advanced users and debugging: domain variable(s),
#'     \code{estimate_synthetic}, \code{variance_synthetic},
#'     \code{correction}, \code{variance_correction},
#'     \code{estimate_final}, \code{variance_final},
#'     \code{se_final}, \code{rse_final}, \code{n_model}, \code{n_proj}.}
#'   \item{\code{diagnostics}}{Named list of model diagnostics:
#'     \code{fixed_effects} (named numeric vector of fixed-effect coefficients),
#'     \code{variance_components} (data frame from \code{lme4::VarCorr}),
#'     \code{random_effect_groups} (character vector of grouping factor names),
#'     \code{singular_fit} (logical),
#'     \code{convergence_messages} (character vector, empty if none),
#'     \code{sigma} (residual standard deviation),
#'     \code{nobs} (number of observations used in model fitting),
#'     \code{aic}, \code{bic}, \code{loglik} (REML-based fit statistics),
#'     \code{icc} (intraclass correlation coefficient for single random-intercept
#'     models; \code{NA} for random-slope or multi-group structures).
#'     For residual diagnostics and additional model evaluation, use
#'     \code{result$fitted_model} directly.}
#'   \item{\code{notes}}{Character vector of methodological notes.}
#'   \item{\code{unit_projection}}{Unit-level data with \code{.y_hat}
#'     (only if \code{keep_unit = TRUE}).}
#'   \item{\code{unit_model_residual}}{Unit-level data with \code{.y_hat_model}
#'     and \code{.resid} (only if \code{keep_unit = TRUE}).}
#' }
#'
#' @references
#' Kim, J.K. & Rao, J.N.K. (2012). Combining data from two independent surveys:
#' a model-assisted approach. \emph{Biometrika}, 99(1), 85-100.
#'
#' Bates, D., Maechler, M., Bolker, B. & Walker, S. (2015). Fitting linear
#' mixed-effects models using lme4. \emph{Journal of Statistical Software},
#' 67(1), 1-48.
#'
#' Lumley, T. (2010). \emph{Complex Surveys: A Guide to Analysis Using R}. Wiley.
#'
#' @examples
#' \dontrun{
#' result <- sae_ml_linear(
#'   formula    = income ~ educ + age + (1 | kabkot),
#'   data_model = survey_model,
#'   data_proj  = survey_proj,
#'   domain     = "kabkot",
#'   weight     = "w",
#'   estimator  = "bias_corrected"
#' )
#'
#' # Final estimates only
#' result$estimates
#' as.data.frame(result)
#'
#' # Estimation components for debugging
#' result$estimation_details
#'
#' # Model diagnostics
#' result$diagnostics$aic
#' result$diagnostics$icc
#'
#' # Print and summary
#' print(result)
#' summary(result)
#' }
#'
#' @export
sae_ml_linear <- function(
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
) {
  mc        <- match.call()
  estimator <- match.arg(estimator)

  # -- Input validation -------------------------------------------------------
  if (!is.data.frame(data_model)) stop("`data_model` must be a data.frame.", call. = FALSE)
  if (!is.data.frame(data_proj))  stop("`data_proj` must be a data.frame.",  call. = FALSE)
  if (!inherits(formula, "formula")) stop("`formula` must be a formula.", call. = FALSE)

  response_var <- .get_response_var(formula)
  fixed_vars   <- .get_fixed_vars(formula)
  group_vars   <- .get_group_vars(formula)
  domain_chr   <- .var_names(domain)

  if (length(group_vars) == 0L)
    stop("`formula` must include at least one random effect, e.g. (1 | area).", call. = FALSE)
  if (length(domain_chr) == 0L)
    stop("`domain` must identify at least one domain variable.", call. = FALSE)

  svy_vars_model <- .var_names(cluster_ids)
  svy_vars_proj  <- .var_names(cluster_ids)
  if (!is.null(weight)) {
    svy_vars_model <- c(svy_vars_model, .var_names(weight))
    svy_vars_proj  <- c(svy_vars_proj,  .var_names(weight))
  }
  if (!is.null(strata)) {
    svy_vars_model <- c(svy_vars_model, .var_names(strata))
    svy_vars_proj  <- c(svy_vars_proj,  .var_names(strata))
  }

  req_model <- unique(c(response_var, fixed_vars, group_vars, domain_chr, svy_vars_model))
  req_proj  <- unique(c(fixed_vars, domain_chr, svy_vars_proj))
  # group_vars in data_proj is optional (re.form = NA); check only if present
  grp_in_proj <- intersect(group_vars, names(data_proj))
  req_proj    <- unique(c(req_proj, grp_in_proj))

  .check_cols(data_model, req_model, "data_model")
  .check_cols(data_proj,  req_proj,  "data_proj")
  .check_missing_values(data_model, req_model, "data_model")
  .check_missing_values(data_proj,  req_proj,  "data_proj")

  # -- Predictor compatibility ------------------------------------------------
  compat     <- .validate_predictor_compatibility(formula, data_model, data_proj)
  data_model <- compat$data_model
  data_proj  <- compat$data_proj

  # -- Model fitting ----------------------------------------------------------
  fit <- lme4::lmer(
    formula = formula,
    data    = data_model,
    REML    = TRUE,
    control = control,
    ...
  )

  diag  <- .extract_diagnostics(fit)
  notes <- character(0L)

  if (isTRUE(diag$singular_fit)) {
    warning("Singular fit detected. Consider simplifying the random-effect structure.",
            call. = FALSE)
    notes <- c(notes, "Singular fit detected; interpret random effects with caution.")
  }
  if (length(diag$convergence_messages) > 0L) {
    warning("Convergence issue: ",
            paste(diag$convergence_messages, collapse = "; "), call. = FALSE)
  }

  # -- Predictions ------------------------------------------------------------
  y_hat_proj  <- stats::predict(fit, newdata = data_proj,  re.form = NA)
  y_hat_model <- stats::predict(fit, newdata = data_model, re.form = NA)

  data_proj$.y_hat        <- y_hat_proj
  data_model$.y_hat_model <- y_hat_model
  data_model$.resid       <- data_model[[response_var]] - y_hat_model

  # -- Survey design ----------------------------------------------------------
  design_proj  <- .build_svy_design(data_proj,  cluster_ids, weight, strata, ...)
  design_model <- .build_svy_design(data_model, cluster_ids, weight, strata, ...)

  domain_formula <- .to_formula(domain_chr)

  # -- Aggregation: synthetic estimates ---------------------------------------
  est_proj <- as.data.frame(survey::svyby(
    formula    = ~.y_hat,
    by         = domain_formula,
    design     = design_proj,
    FUN        = survey::svymean,
    vartype    = "var",
    na.rm      = TRUE,
    keep.names = FALSE
  ))

  names(est_proj)[names(est_proj) == ".y_hat"]      <- "estimate_synthetic"
  # svyby vartype="var" may append "var" or "var.varname"; handle both
  if (!"variance_synthetic" %in% names(est_proj)) {
    var_col <- grep("^var", names(est_proj), value = TRUE)[1L]
    if (!is.na(var_col))
      names(est_proj)[names(est_proj) == var_col] <- "variance_synthetic"
  }

  correction_col <- rep(0, nrow(est_proj))
  var_correction <- rep(0, nrow(est_proj))
  detail_note    <- "Synthetic estimator; no residual correction."

  # -- Bias correction --------------------------------------------------------
  if (estimator == "bias_corrected") {
    est_resid <- as.data.frame(survey::svyby(
      formula    = ~.resid,
      by         = domain_formula,
      design     = design_model,
      FUN        = survey::svymean,
      vartype    = "var",
      na.rm      = TRUE,
      keep.names = FALSE
    ))

    names(est_resid)[names(est_resid) == ".resid"] <- "correction"
    var_col_r <- grep("^var", names(est_resid), value = TRUE)[1L]
    if (!is.na(var_col_r))
      names(est_resid)[names(est_resid) == var_col_r] <- "variance_correction"

    merged <- merge(
      est_proj,
      est_resid[, c(domain_chr, "correction", "variance_correction")],
      by    = domain_chr,
      all.x = TRUE
    )
    merged$correction[is.na(merged$correction)]               <- 0
    merged$variance_correction[is.na(merged$variance_correction)] <- 0

    est_proj       <- merged
    correction_col <- est_proj$correction
    var_correction <- est_proj$variance_correction
    detail_note    <- "Bias-corrected estimator; residual correction is added."
  }

  # -- Final estimates --------------------------------------------------------
  estimate_final  <- est_proj$estimate_synthetic + correction_col
  variance_final  <- est_proj$variance_synthetic + var_correction
  se_final        <- sqrt(variance_final)
  rse_final       <- ifelse(
    estimate_final == 0 | is.na(estimate_final),
    NA_real_,
    100 * se_final / abs(estimate_final)
  )

  # Ensure correction and variance_correction columns exist for details
  if (!"correction" %in% names(est_proj))          est_proj$correction          <- 0
  if (!"variance_correction" %in% names(est_proj)) est_proj$variance_correction <- 0

  # Domain counts
  n_proj  <- .domain_counts(data_proj,  domain_chr, "n_proj")
  n_model <- .domain_counts(data_model, domain_chr, "n_model")

  # -- Build estimates (user-facing, concise) ---------------------------------
  estimates <- est_proj[, domain_chr, drop = FALSE]
  estimates$estimate <- estimate_final
  estimates$variance <- variance_final
  estimates$se       <- se_final
  estimates$rse      <- rse_final
  row.names(estimates) <- NULL

  # -- Build estimation_details (components for advanced users) ---------------
  estimation_details <- est_proj[, domain_chr, drop = FALSE]
  estimation_details$estimate_synthetic  <- est_proj$estimate_synthetic
  estimation_details$variance_synthetic  <- est_proj$variance_synthetic
  estimation_details$correction          <- est_proj$correction
  estimation_details$variance_correction <- est_proj$variance_correction
  estimation_details$estimate_final      <- estimate_final
  estimation_details$variance_final      <- variance_final
  estimation_details$se_final            <- se_final
  estimation_details$rse_final           <- rse_final

  estimation_details <- merge(estimation_details, n_model, by = domain_chr, all.x = TRUE)
  estimation_details <- merge(estimation_details, n_proj,  by = domain_chr, all.x = TRUE)
  row.names(estimation_details) <- NULL

  # -- Notes ------------------------------------------------------------------
  notes <- c(notes,
             "Predictions use fixed effects only (re.form = NA).",
             "Plug-in variance from survey::svyby() is approximate.",
             if (estimator == "bias_corrected")
               "Bias-corrected: variance = var(synthetic) + var(residual correction)."
             else
               "Synthetic estimator: no empirical residual correction.",
               "For further model evaluation, use summary(result$fitted_model) or lme4 diagnostic tools.",
             if (is.na(diag$icc))
               "ICC is NA: random-slope or multi-group model detected; ICC requires a single random-intercept structure."
  )
  notes <- notes[!vapply(notes, is.null, logical(1L))]

  # -- Assemble output --------------------------------------------------------
  out <- list(
    call               = mc,
    formula            = formula,
    estimator          = estimator,
    fitted_model       = fit,
    estimates          = estimates,
    estimation_details = estimation_details,
    diagnostics        = diag,
    notes              = notes
  )

  if (keep_unit) {
    out$unit_projection <- data_proj
    out$unit_model_residual <- data_model[
      , c(domain_chr, response_var, ".y_hat_model", ".resid"),
      drop = FALSE
    ]
  }

  structure(out, class = "sae_ml_linear")
}

# ---- S3 methods --------------------------------------------------------------

#' Print method for sae_ml_linear
#'
#' Displays a concise summary: title, formula, estimator, number of domains,
#' and a preview of the final estimates table.
#'
#' @param x An object of class \code{"sae_ml_linear"}.
#' @param n Integer. Number of rows to preview from \code{estimates}.
#' @param ... Further arguments (currently unused).
#'
#' @export
#' @method print sae_ml_linear
print.sae_ml_linear <- function(x, n = 6L, ...) {
  cat("SAE Projection Estimator \u2014 Linear Multilevel Model\n")
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
#' Displays a detailed summary including the model call, formula, estimator,
#' fitted model summary, key diagnostics, final estimates, notes, and a
#' pointer to \code{object$estimation_details}.
#'
#' @param object An object of class \code{"sae_ml_linear"}.
#' @param ... Further arguments (currently unused).
#'
#' @export
#' @method summary sae_ml_linear
summary.sae_ml_linear <- function(object, ...) {
  cat("=== SAE Projection Estimator \u2014 Linear Multilevel Model ===\n\n")
  cat("Call:\n"); print(object$call); cat("\n")
  cat("Formula   :", deparse(object$formula), "\n")
  cat("Estimator :", object$estimator, "\n\n")

  cat("--- Fitted Model ---\n")
  print(summary(object$fitted_model))
  cat("\n")

  d <- object$diagnostics
  cat("--- Diagnostics ---\n")
  cat(" sigma      :", round(d$sigma, 4L), "\n")
  cat(" nobs       :", d$nobs, "\n")
  cat(" AIC        :", round(d$aic,    2L), "\n")
  cat(" BIC        :", round(d$bic,    2L), "\n")
  cat(" logLik     :", round(d$loglik, 2L), "\n")
  icc_str <- if (is.na(d$icc)) "NA (random-slope or multi-group model)" else round(d$icc, 4L)
  cat(" ICC        :", icc_str, "\n")
  cat(" singular   :", d$singular_fit, "\n")
  conv_str <- if (length(d$convergence_messages) > 0L)
    paste(d$convergence_messages, collapse = "; ") else "OK"
  cat(" convergence:", conv_str, "\n\n")

  cat("--- Final Estimates ---\n")
  print(object$estimates, row.names = FALSE)
  cat("\n")

  if (length(object$notes) > 0L) {
    cat("--- Notes ---\n")
    for (n in object$notes) cat(" *", n, "\n")
    cat("\n")
  }

  cat("Additional estimation components are available in object$estimation_details.\n")

  invisible(object)
}

#' Coerce an sae_ml_linear object to a data frame
#'
#' Returns \code{x$estimates}: the final domain-level estimates table with
#' columns domain variable(s), \code{estimate}, \code{variance}, \code{se},
#' and \code{rse}.
#'
#' @param x An object of class \code{"sae_ml_linear"}.
#' @param row.names Passed to \code{\link{as.data.frame}} (unused; included for
#'   S3 compatibility).
#' @param optional Passed to \code{\link{as.data.frame}} (unused; included for
#'   S3 compatibility).
#' @param ... Further arguments (currently unused).
#'
#' @export
#' @method as.data.frame sae_ml_linear
as.data.frame.sae_ml_linear <- function(x, row.names = NULL, optional = FALSE, ...) {
  x$estimates
}
