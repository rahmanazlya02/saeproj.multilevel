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
# Accepts: formula object, formula string ("~kabkot"), character vector ("kabkot"),
# "1" / "~1" (intercept-only), or NULL.
.to_formula <- function(x) {
  if (is.null(x))             return(NULL)
  if (inherits(x, "formula")) return(x)
  if (!is.character(x))
    stop("Cannot convert input to a one-sided formula.", call. = FALSE)

  x <- trimws(x)
  if (length(x) == 0L || all(!nzchar(x)) || identical(x, "1"))
    return(stats::as.formula("~1"))
  if (length(x) == 1L && grepl("^\\s*~", x))
    return(stats::as.formula(x))

  stats::as.formula(paste0("~", paste(x[nzchar(x)], collapse = " + ")))
}

#' @noRd
# Extracts variable names from formula, formula string, or character vector.
.var_names <- function(x) {
  if (is.null(x))             return(character(0L))
  if (inherits(x, "formula")) return(all.vars(x))
  if (!is.character(x))       return(character(0L))
  if (length(x) == 1L && grepl("^\\s*~", x))
    return(all.vars(stats::as.formula(x)))
  unique(x[nzchar(x)])
}

#' @noRd
.check_cols <- function(data, cols, data_name) {
  cols    <- unique(cols[!is.na(cols) & nzchar(cols)])
  missing <- setdiff(cols, names(data))
  if (length(missing) > 0L)
    stop("Column(s) not found in ", data_name, ": ",
         paste(missing, collapse = ", "), call. = FALSE)
  invisible(TRUE)
}

#' @noRd
# Stops if any required variable contains missing values.
# NA in predictors, domain, or design variables will silently distort SAE estimates.
.check_missing_values <- function(data, vars, data_name) {
  vars <- intersect(unique(vars), names(data))
  if (length(vars) == 0L) return(invisible(TRUE))

  na_counts <- vapply(vars, function(v) sum(is.na(data[[v]])), integer(1L))
  bad       <- na_counts[na_counts > 0L]
  if (length(bad) == 0L) return(invisible(TRUE))

  stop(
    "Missing values in ", data_name, ": ",
    paste(paste0(names(bad), " (", bad, ")"), collapse = ", "), ". ",
    "Handle missing values before calling sae_ml_linear().",
    call. = FALSE
  )
}

#' @noRd
# Validates weight variable: must be numeric and not all zero.
# Warns if negative or zero values are present (unusual but not always fatal).
.check_weight <- function(data, weight, data_name) {
  wvars <- .var_names(weight)
  if (length(wvars) == 0L) return(invisible(TRUE))
  if (length(wvars) > 1L)
    stop("`weight` must identify exactly one weight variable.", call. = FALSE)

  w <- data[[wvars]]
  if (!is.numeric(w))
    stop("Weight variable `", wvars, "` in ", data_name, " must be numeric.", call. = FALSE)
  if (all(w == 0, na.rm = TRUE))
    stop("Weight variable `", wvars, "` in ", data_name, " contains only zero values.", call. = FALSE)
  if (any(w < 0, na.rm = TRUE))
    warning("Weight variable `", wvars, "` in ", data_name,
            " contains negative values. This is unusual for standard survey weights; ",
            "verify the weight construction before proceeding.", call. = FALSE)
  if (any(w == 0, na.rm = TRUE))
    warning("Weight variable `", wvars, "` in ", data_name,
            " contains zero values; verify that this is intended.", call. = FALSE)

  invisible(TRUE)
}

#' @noRd
# Renames columns returned by survey::svyby() to consistent names.
# svyby() names the variance column inconsistently across versions
# ("var.<col>" in newer, "var" in older). This handles both.
.rename_svyby_output <- function(x, domain_chr, value_col, estimate_name, variance_name) {
  if (value_col %in% names(x))
    names(x)[names(x) == value_col] <- estimate_name

  if (!estimate_name %in% names(x))
    stop("Cannot find estimate column `", value_col, "` returned by survey::svyby().",
         call. = FALSE)

  if (!variance_name %in% names(x)) {
    non_domain_cols <- setdiff(names(x), c(domain_chr, estimate_name))
    var_cols        <- grep("^(var|se2)(\\.|$)", non_domain_cols, value = TRUE)

    if (length(var_cols) != 1L) {
      if (length(non_domain_cols) == 1L) {
        var_cols <- non_domain_cols
      } else {
        stop("Cannot uniquely identify variance column returned by survey::svyby(); ",
             "columns found: ", paste(non_domain_cols, collapse = ", "), call. = FALSE)
      }
    }
    names(x)[names(x) == var_cols] <- variance_name
  }

  x
}

#' @noRd
# Builds a survey design object for use in svyby() aggregation.
# `...` is passed to svydesign() only — not to lmer().
.build_svy_design <- function(data, ids, weight, strata, ...) {
  ids    <- if (is.null(ids))    stats::as.formula("~1") else .to_formula(ids)
  weight <- if (is.null(weight)) stats::as.formula("~1") else .to_formula(weight)

  args <- list(ids = ids, weights = weight, data = data, ...)
  if (!is.null(strata)) args$strata <- .to_formula(strata)
  do.call(survey::svydesign, args)
}

#' @noRd
# Unweighted row count per domain — reflects model data availability,
# not population size. Stored in estimation_details as n_model / n_proj.
.domain_counts <- function(data, domain_chr, col_name) {
  out <- stats::aggregate(rep(1L, nrow(data)), data[domain_chr], length)
  names(out)[ncol(out)] <- col_name
  out
}

#' @noRd
# Builds a readable domain label for warning messages, e.g. "prov / kabkot".
.domain_labels <- function(data, domain_chr) {
  if (nrow(data) == 0L) return(character(0L))
  apply(data[domain_chr], 1L, function(z) paste(z, collapse = " / "))
}

#' @noRd
# Returns unique sorted levels for categorical variables.
.model_levels <- function(x) {
  if (is.factor(x))  return(levels(x))
  if (is.logical(x)) return(c("FALSE", "TRUE"))
  sort(unique(as.character(x[!is.na(x)])))
}

#' @noRd
# Checks type consistency and harmonizes factor levels between data_model and
# data_proj. Stops if types are incompatible or data_proj has unseen levels
# (predict() would silently produce NA for those observations).
.validate_predictor_compatibility <- function(formula, data_model, data_proj) {
  fixed_vars <- .get_fixed_vars(formula)
  group_vars <- .get_group_vars(formula)
  dm <- data_model
  dp <- data_proj

  is_cat <- function(x) is.factor(x) || is.character(x) || is.logical(x)

  for (v in fixed_vars) {
    xm <- dm[[v]]; xp <- dp[[v]]

    if (is_cat(xm) || is_cat(xp)) {
      if (is.numeric(xm) || is.numeric(xp))
        stop("Variable '", v, "' has incompatible types between data_model and data_proj ",
             "(categorical vs numeric). Harmonize types before calling sae_ml_linear().",
             call. = FALSE)

      lm_vals  <- .model_levels(xm)
      new_lvls <- setdiff(unique(as.character(xp[!is.na(xp)])), lm_vals)
      if (length(new_lvls) > 0L)
        stop("Variable '", v, "' has level(s) in data_proj not in data_model: ",
             paste(utils::head(new_lvls, 10L), collapse = ", "), ". ",
             "Harmonize labels before calling sae_ml_linear().",
             call. = FALSE)

      dm[[v]] <- factor(as.character(xm), levels = lm_vals, ordered = is.ordered(xm))
      dp[[v]] <- factor(as.character(xp), levels = lm_vals, ordered = is.ordered(xm))

    } else if (is.numeric(xm) != is.numeric(xp)) {
      stop("Variable '", v, "' has incompatible types between data_model and data_proj.",
           call. = FALSE)
    }
  }

  # Grouping variable is optional in data_proj (re.form = NA); harmonize only
  # if present, so downstream diagnostics remain stable.
  for (v in group_vars) {
    if (!v %in% names(dp)) next
    xm <- dm[[v]]; xp <- dp[[v]]
    if (is_cat(xm) || is_cat(xp)) {
      all_lvls <- union(.model_levels(xm), unique(as.character(xp[!is.na(xp)])))
      dm[[v]]  <- factor(as.character(xm), levels = all_lvls, ordered = is.ordered(xm))
      dp[[v]]  <- factor(as.character(xp), levels = all_lvls, ordered = is.ordered(xm))
    }
  }

  list(data_model = dm, data_proj = dp)
}

#' @noRd
# Computes ICC for a single random-intercept model.
# Returns NA for random-slope or multi-group models where ICC is not well-defined.
.compute_icc <- function(fit) {
  vc     <- lme4::VarCorr(fit)
  groups <- names(vc)
  if (length(groups) != 1L)               return(NA_real_)
  if (!identical(dim(vc[[1L]]), c(1L, 1L))) return(NA_real_)
  var_u <- as.numeric(vc[[1L]][1L, 1L])
  var_e <- stats::sigma(fit)^2
  var_u / (var_u + var_e)
}

#' @noRd
.extract_diagnostics <- function(fit) {
  conv_msg <- fit@optinfo$conv$lme4$messages
  if (is.null(conv_msg)) conv_msg <- character(0L)

  list(
    fixed_effects        = lme4::fixef(fit),
    variance_components  = as.data.frame(lme4::VarCorr(fit)),
    singular_fit         = lme4::isSingular(fit),
    convergence_messages = conv_msg,
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
#' The working model is fitted with the multilevel random-effects structure
#' specified in \code{formula}. For projection, synthetic predictions are computed
#' with \code{re.form = NA}, so random effects/BLUPs are not added to the
#' projected values. The \code{cluster_ids}, \code{weight}, and \code{strata}
#' arguments are used exclusively in \code{\link[survey]{svydesign}} for the
#' aggregation step, not in model fitting.
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
#' @param ... Additional named arguments passed to \code{\link[survey]{svydesign}}
#'   only (e.g. \code{nest = TRUE} when PSU IDs are not unique across strata).
#'   These are \strong{not} forwarded to \code{\link[lme4]{lmer}};
#'   use \code{control} for lmer-specific tuning.
#'
#' @return An object of class \code{"sae_ml_linear"}, a list with:
#' \describe{
#'   \item{\code{call}}{The matched call.}
#'   \item{\code{formula}}{The model formula.}
#'   \item{\code{estimator}}{The estimator type used.}
#'   \item{\code{fitted_model}}{The fitted \code{lmerMod} object.}
#'   \item{\code{estimates}}{Data frame of final domain-level estimates:
#'     domain variable(s), \code{estimate}, \code{variance}, \code{se}, \code{rse}.}
#'   \item{\code{estimation_details}}{Estimation components for debugging:
#'     \code{estimate_synthetic}, \code{variance_synthetic}, \code{correction},
#'     \code{variance_correction}, \code{estimate_final}, \code{variance_final},
#'     \code{se_final}, \code{rse_final}, \code{n_model}, \code{n_proj}.}
#'   \item{\code{diagnostics}}{Model diagnostics: fixed effects, variance components,
#'     singular fit, convergence messages, sigma, nobs, AIC, BIC, logLik, ICC.}
#'   \item{\code{notes}}{Conditional notes generated from model or estimator
#'     conditions, such as singular fit, undefined ICC, or the bias-corrected
#'     variance assumption. Static methodological notes are documented in the
#'     sections below and are not repeated in the output object.}
#'   \item{\code{unit_projection}}{Unit-level predictions (only if \code{keep_unit = TRUE}).}
#'   \item{\code{unit_model_residual}}{Unit-level residuals (only if \code{keep_unit = TRUE}).}
#' }
#'
#' @section Variance assumptions:
#' The plug-in variance for the bias-corrected estimator is computed as
#' \eqn{V(\hat{Y}_d^{BC}) \approx V(\hat{Y}_d^{syn}) + V(\bar{e}_d)},
#' which assumes \code{data_model} and \code{data_proj} are from
#' \strong{independent surveys} (Kim & Rao, 2012, Section 3). If \code{data_model}
#' is a subsample of \code{data_proj}, the covariance term is non-zero and variance
#' will be underestimated. This plug-in variance is approximate and does not
#' fully account for model-parameter uncertainty.
#'
#' @section Prediction strategy:
#' The working model is fitted using the full multilevel structure specified in
#' \code{formula}, including random effects. However, projection predictions are computed with
#' \code{re.form = NA}, so the synthetic values use only the fixed-effect
#' component: \eqn{\hat{y}_{ij} = x_{ij}^T \hat{\beta}}.
#' Random effects/BLUPs are not added to the projected values. This allows prediction for domains or
#' grouping levels that are absent from \code{data_model}.
#' For the bias-corrected estimator, remaining domain-level discrepancies are accounted for through
#' the design-weighted mean residual correction computed from \code{data_model}.
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
#' result$estimates
#' as.data.frame(result)
#' result$estimation_details
#' summary(result)
#'
#' # With nested PSUs (PSU IDs not unique across strata)
#' result2 <- sae_ml_linear(
#'   formula     = income ~ educ + (1 | prov),
#'   data_model  = survey_model,
#'   data_proj   = survey_proj,
#'   domain      = "kabkot",
#'   cluster_ids = "psu_id",
#'   weight      = "w",
#'   strata      = "strata_id",
#'   nest        = TRUE
#' )
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
  if (!is.data.frame(data_model))
    stop("`data_model` must be a data.frame.", call. = FALSE)
  if (!is.data.frame(data_proj))
    stop("`data_proj` must be a data.frame.", call. = FALSE)
  if (!inherits(formula, "formula") || length(formula) != 3L)
    stop("`formula` must be a two-sided formula, e.g. y ~ x + (1 | area).", call. = FALSE)

  response_var <- .get_response_var(formula)
  fixed_vars   <- .get_fixed_vars(formula)
  group_vars   <- .get_group_vars(formula)
  domain_chr   <- .var_names(domain)

  if (length(group_vars) == 0L)
    stop("`formula` must include at least one random effect, e.g. (1 | area).", call. = FALSE)
  if (length(domain_chr) == 0L)
    stop("`domain` must identify at least one domain variable.", call. = FALSE)

  # Warn if cluster/strata declared but no weight — silently assigns equal weights
  cluster_formula <- if (is.null(cluster_ids)) ~1 else .to_formula(cluster_ids)
  if (is.null(weight) && (!is.null(strata) || !identical(cluster_formula, ~1)))
    warning(
      "Cluster or strata specified but `weight = NULL`; equal weights assumed. ",
      "Set `weight` for design-consistent estimates.",
      call. = FALSE
    )

  # -- Column checks ----------------------------------------------------------
  svy_vars    <- unique(c(.var_names(cluster_ids), .var_names(weight), .var_names(strata)))
  req_model   <- unique(c(response_var, fixed_vars, group_vars, domain_chr, svy_vars))
  req_proj    <- unique(c(fixed_vars, domain_chr, svy_vars,
                          intersect(group_vars, names(data_proj))))

  .check_cols(data_model, req_model, "data_model")
  .check_cols(data_proj,  req_proj,  "data_proj")
  .check_missing_values(data_model, req_model, "data_model")
  .check_missing_values(data_proj,  req_proj,  "data_proj")
  .check_weight(data_model, weight, "data_model")
  .check_weight(data_proj,  weight, "data_proj")

  # Guard against reserved column names that would be silently overwritten
  reserved <- c(".y_hat", ".y_hat_model", ".resid")
  clash_proj  <- intersect(reserved, names(data_proj))
  clash_model <- intersect(reserved, names(data_model))
  if (length(clash_proj) > 0L)
    stop("Reserved column(s) already exist in `data_proj`: ",
         paste(clash_proj, collapse = ", "), ". Rename them before calling sae_ml_linear().",
         call. = FALSE)
  if (length(clash_model) > 0L)
    stop("Reserved column(s) already exist in `data_model`: ",
         paste(clash_model, collapse = ", "), ". Rename them before calling sae_ml_linear().",
         call. = FALSE)

  # -- Predictor compatibility ------------------------------------------------
  compat     <- .validate_predictor_compatibility(formula, data_model, data_proj)
  data_model <- compat$data_model
  data_proj  <- compat$data_proj

  # -- Model fitting ----------------------------------------------------------
  # `...` is NOT forwarded to lmer(); use `control` for lmer-specific tuning.
  fit   <- lme4::lmer(formula = formula, data = data_model, REML = TRUE, control = control)
  diag  <- .extract_diagnostics(fit)
  notes <- character(0L)

  if (isTRUE(diag$singular_fit)) {
    warning("Singular fit detected. Consider simplifying the random-effect structure.",
            call. = FALSE)
    notes <- c(notes, "Singular fit detected; interpret random effects with caution.")
  }
  if (length(diag$convergence_messages) > 0L) {
    conv_text <- paste(diag$convergence_messages, collapse = "; ")
    warning("Convergence issue: ", conv_text, call. = FALSE)
    notes <- c(notes, paste0("Convergence issue: ", conv_text))
  }

  # -- Predictions ------------------------------------------------------------
  # The model is fitted with random effects, but projection predictions use
  # re.form = NA; random effects/BLUPs are not added to synthetic values.
  # Residuals are e_i = y_i - x_i' * beta_hat, not BLUP residuals.
  data_proj$.y_hat        <- stats::predict(fit, newdata = data_proj,  re.form = NA)
  data_model$.y_hat_model <- stats::predict(fit, newdata = data_model, re.form = NA)
  data_model$.resid       <- data_model[[response_var]] - data_model$.y_hat_model

  # -- Survey design ----------------------------------------------------------
  design_args <- list(ids = cluster_ids, weight = weight, strata = strata, ...)

  design_proj  <- do.call(.build_svy_design, c(list(data = data_proj),  design_args))
  design_model <- do.call(.build_svy_design, c(list(data = data_model), design_args))

  domain_formula <- .to_formula(domain_chr)

  # -- Synthetic estimates ----------------------------------------------------
  est_proj <- as.data.frame(survey::svyby(
    formula = ~.y_hat, by = domain_formula, design = design_proj,
    FUN = survey::svymean, vartype = "var", na.rm = TRUE, keep.names = FALSE
  ))
  est_proj <- .rename_svyby_output(est_proj, domain_chr,
                                   ".y_hat", "estimate_synthetic", "variance_synthetic")

  correction_col <- rep(0, nrow(est_proj))
  var_correction <- rep(0, nrow(est_proj))

  # -- Bias correction --------------------------------------------------------
  if (estimator == "bias_corrected") {
    est_resid <- as.data.frame(survey::svyby(
      formula = ~.resid, by = domain_formula, design = design_model,
      FUN = survey::svymean, vartype = "var", na.rm = TRUE, keep.names = FALSE
    ))
    est_resid <- .rename_svyby_output(est_resid, domain_chr,
                                      ".resid", "correction", "variance_correction")

    merged <- merge(
      est_proj,
      est_resid[, c(domain_chr, "correction", "variance_correction")],
      by = domain_chr, all.x = TRUE
    )

    # Domains in data_proj with no data_model observations get no correction
    no_corr <- .domain_labels(merged[is.na(merged$correction), domain_chr, drop = FALSE],
                              domain_chr)
    if (length(no_corr) > 0L)
      warning(
        length(no_corr), " domain(s) have no observations in data_model; ",
        "residual correction set to 0 (synthetic estimator used). Domain(s): ",
        paste(utils::head(no_corr, 5L), collapse = ", "),
        if (length(no_corr) > 5L) " [... and more]" else "",
        call. = FALSE
      )

    merged$correction[is.na(merged$correction)]                   <- 0
    merged$variance_correction[is.na(merged$variance_correction)] <- 0

    est_proj       <- merged
    correction_col <- est_proj$correction
    var_correction <- est_proj$variance_correction
  }

  # -- Final estimates --------------------------------------------------------
  # Var(Y_d^BC) = Var(Y_d^syn) + Var(e_bar_d), assuming data_model and data_proj
  # are independent surveys (Kim & Rao 2012, Section 3).
  estimate_final <- est_proj$estimate_synthetic + correction_col
  variance_final <- est_proj$variance_synthetic  + var_correction

  # Plug-in variance can be numerically negative for very small domains;
  # clamp to 0 before sqrt to avoid NaN.
  neg_var_idx <- which(!is.na(variance_final) & variance_final < 0)
  if (length(neg_var_idx) > 0L) {
    neg_domains <- .domain_labels(est_proj[neg_var_idx, domain_chr, drop = FALSE],
                                  domain_chr)
    warning(
      length(neg_var_idx), " domain(s) have negative plug-in variance; clamped to 0. ",
      "Consider the synthetic estimator for these domains. Domain(s): ",
      paste(utils::head(neg_domains, 5L), collapse = ", "),
      if (length(neg_domains) > 5L) " [... and more]" else "",
      call. = FALSE
    )
    variance_final[neg_var_idx] <- 0
  }

  se_final  <- sqrt(variance_final)
  rse_final <- ifelse(
    estimate_final == 0 | is.na(estimate_final), NA_real_,
    100 * se_final / abs(estimate_final)
  )

  if (!"correction" %in% names(est_proj))          est_proj$correction          <- 0
  if (!"variance_correction" %in% names(est_proj)) est_proj$variance_correction <- 0

  # -- Build output tables ----------------------------------------------------
  n_proj  <- .domain_counts(data_proj,  domain_chr, "n_proj")
  n_model <- .domain_counts(data_model, domain_chr, "n_model")

  estimates <- est_proj[, domain_chr, drop = FALSE]
  estimates$estimate <- estimate_final
  estimates$variance <- variance_final
  estimates$se       <- se_final
  estimates$rse      <- rse_final

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
  estimation_details$n_model[is.na(estimation_details$n_model)] <- 0L
  estimation_details$n_proj[is.na(estimation_details$n_proj)]   <- 0L

  # Re-sort both tables by domain columns (merge() does not guarantee row order)
  estimates          <- estimates[do.call(order, estimates[, domain_chr, drop = FALSE]), ,
                                  drop = FALSE]
  estimation_details <- estimation_details[do.call(order,
                                                   estimation_details[, domain_chr, drop = FALSE]), ,
                                           drop = FALSE]
  row.names(estimates)          <- NULL
  row.names(estimation_details) <- NULL

  # -- Conditional notes ------------------------------------------------------
  # Static methodological notes are documented in roxygen sections above.
  # Only run-specific notes are stored in the output object.
  if (estimator == "bias_corrected") {
    notes <- c(notes, paste0(
      "Bias-corrected: Var = Var(synthetic) + Var(correction). ",
      "Assumes data_model and data_proj are independent surveys ",
      "(Kim & Rao 2012, Section 3)."
    ))
  }
  if (is.na(diag$icc))
    notes <- c(notes,
               "ICC is NA: random-slope or multi-group model; ICC requires single random-intercept.")

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
    out$unit_projection     <- data_proj[,
                                         intersect(c(domain_chr, ".y_hat"), names(data_proj)), drop = FALSE]
    out$unit_model_residual <- data_model[,
                                          intersect(c(domain_chr, response_var, ".y_hat_model", ".resid"), names(data_model)),
                                          drop = FALSE]
  }

  structure(out, class = "sae_ml_linear")
}

# ---- S3 methods --------------------------------------------------------------

#' Print method for sae_ml_linear
#'
#' @param x An object of class \code{"sae_ml_linear"}.
#' @param n Integer. Number of rows to preview from \code{estimates}.
#' @param ... Further arguments (currently unused).
#' @export
#' @method print sae_ml_linear
print.sae_ml_linear <- function(x, n = 6L, ...) {
  if (!is.numeric(n) || length(n) != 1L || is.na(n) || n < 0L)
    stop("`n` must be a single non-negative number.", call. = FALSE)
  if (is.infinite(n)) n <- nrow(x$estimates)
  n <- as.integer(n)

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
#' @param object An object of class \code{"sae_ml_linear"}.
#' @param ... Further arguments (currently unused).
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
  cat(" sigma      :", round(d$sigma,  4L), "\n")
  cat(" nobs       :", d$nobs, "\n")
  cat(" ICC        :", if (is.na(d$icc)) "NA (random-slope or multi-group)" else round(d$icc, 4L), "\n")
  cat(" singular   :", d$singular_fit, "\n")
  cat(" convergence:", if (length(d$convergence_messages) > 0L)
    paste(d$convergence_messages, collapse = "; ") else "OK", "\n\n")

  cat("--- Final Estimates ---\n")
  print(object$estimates, row.names = FALSE)
  cat("\n")

  if (!is.null(object$notes) && length(object$notes) > 0L) {
    cat("--- Notes ---\n")
    for (n in object$notes) cat(" *", n, "\n")
    cat("\n")
  }

  cat("See object$estimation_details for estimation components.\n")
  invisible(object)
}

#' Coerce an sae_ml_linear object to a data frame
#'
#' Returns \code{x$estimates}: domain-level estimates with columns
#' domain variable(s), \code{estimate}, \code{variance}, \code{se}, \code{rse}.
#'
#' @param x An object of class \code{"sae_ml_linear"}.
#' @param row.names Passed to \code{\link{as.data.frame}} (unused; for S3 compatibility).
#' @param optional Passed to \code{\link{as.data.frame}} (unused; for S3 compatibility).
#' @param ... Further arguments (currently unused).
#' @export
#' @method as.data.frame sae_ml_linear
as.data.frame.sae_ml_linear <- function(x, row.names = NULL, optional = FALSE, ...) {
  x$estimates
}
