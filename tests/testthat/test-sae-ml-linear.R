# tests/testthat/test-sae-ml-linear.R

make_test_data <- function(seed = 123) {
  set.seed(seed)

  n_area <- 6
  n_model <- 120
  n_proj <- 600

  area_model <- sample(paste0("A", 1:n_area), n_model, replace = TRUE)

  data_model <- data.frame(
    domain = paste0("D", ceiling(as.integer(substr(area_model, 2, 2)) / 2)),
    area = area_model,
    x1 = rnorm(n_model),
    x2 = runif(n_model),
    education = factor(
      sample(c("SD", "SMP", "SMA"), n_model, replace = TRUE),
      levels = c("SD", "SMP", "SMA")
    ),
    weight = runif(n_model, 1, 3)
  )

  area_eff <- setNames(rnorm(n_area), paste0("A", 1:n_area))

  data_model$y <- with(
    data_model,
    5 + 0.5 * x1 + 2 * x2 + area_eff[area] + rnorm(n_model, 0, 0.5)
  )

  area_proj <- sample(c(paste0("A", 1:n_area), "A7"), n_proj, replace = TRUE)

  data_proj <- data.frame(
    domain = ifelse(
      area_proj == "A7",
      "D3",
      paste0("D", ceiling(as.integer(substr(area_proj, 2, 2)) / 2))
    ),
    area = area_proj,
    x1 = rnorm(n_proj),
    x2 = runif(n_proj),
    education = factor(
      sample(c("SD", "SMP", "SMA"), n_proj, replace = TRUE),
      levels = c("SD", "SMP", "SMA")
    ),
    weight = runif(n_proj, 1, 3)
  )

  list(data_model = data_model, data_proj = data_proj)
}

test_that("runs with random intercept", {
  d <- make_test_data()

  res <- sae_ml_linear(
    formula = y ~ x1 + x2 + (1 | area),
    data_model = d$data_model,
    data_proj = d$data_proj,
    domain = "domain",
    weight = "weight"
  )

  expect_s3_class(res, "sae_ml_linear")
  expect_true(is.data.frame(res$estimates))
  expect_true(all(c("estimate", "variance", "se", "rse") %in% names(res$estimates)))
})

test_that("runs with random slope", {
  d <- make_test_data()

  expect_no_error(
    sae_ml_linear(
      formula = y ~ x1 + x2 + (1 + x1 | area),
      data_model = d$data_model,
      data_proj = d$data_proj,
      domain = "domain",
      weight = "weight"
    )
  )
})

test_that("data_proj can have new random effect levels", {
  d <- make_test_data()
  expect_true("A7" %in% d$data_proj$area)

  res <- sae_ml_linear(
    formula = y ~ x1 + x2 + (1 | area),
    data_model = d$data_model,
    data_proj = d$data_proj,
    domain = "domain",
    weight = "weight"
  )

  expect_s3_class(res, "sae_ml_linear")
})

test_that("formula without random effect fails", {
  d <- make_test_data()

  expect_error(
    sae_ml_linear(
      formula = y ~ x1 + x2,
      data_model = d$data_model,
      data_proj = d$data_proj,
      domain = "domain",
      weight = "weight"
    ),
    regexp = "random effect"
  )
})

test_that("missing values fail", {
  d <- make_test_data()
  d$data_model$x1[1] <- NA

  expect_error(
    sae_ml_linear(
      formula = y ~ x1 + x2 + (1 | area),
      data_model = d$data_model,
      data_proj = d$data_proj,
      domain = "domain",
      weight = "weight"
    ),
    regexp = "Missing values"
  )
})

test_that("missing required column fails", {
  d <- make_test_data()
  d$data_proj$x2 <- NULL

  expect_error(
    sae_ml_linear(
      formula = y ~ x1 + x2 + (1 | area),
      data_model = d$data_model,
      data_proj = d$data_proj,
      domain = "domain",
      weight = "weight"
    ),
    regexp = "not found"
  )
})

test_that("new fixed-effect categorical levels fail", {
  d <- make_test_data()

  d$data_proj$education <- as.character(d$data_proj$education)
  d$data_proj$education[1] <- "S1"
  d$data_proj$education <- factor(d$data_proj$education)

  expect_error(
    sae_ml_linear(
      formula = y ~ x1 + education + (1 | area),
      data_model = d$data_model,
      data_proj = d$data_proj,
      domain = "domain",
      weight = "weight"
    ),
    regexp = "level"
  )
})

test_that("zero variance predictor is removed with warning", {
  d <- make_test_data()
  d$data_model$x_const <- 1
  d$data_proj$x_const <- 1

  expect_warning(
    res <- sae_ml_linear(
      formula = y ~ x1 + x_const + (1 | area),
      data_model = d$data_model,
      data_proj = d$data_proj,
      domain = "domain",
      weight = "weight"
    ),
    regexp = "zero-variance"
  )

  expect_s3_class(res, "sae_ml_linear")
  expect_true(any(grepl("zero-variance", res$notes)))
})

test_that("summary_function total works", {
  d <- make_test_data()

  res <- sae_ml_linear(
    formula = y ~ x1 + x2 + (1 | area),
    data_model = d$data_model,
    data_proj = d$data_proj,
    domain = "domain",
    weight = "weight",
    summary_function = "total"
  )

  expect_s3_class(res, "sae_ml_linear")
  expect_true(is.data.frame(res$estimates))
})

test_that("return_direct returns direct estimator", {
  d <- make_test_data()

  res <- sae_ml_linear(
    formula = y ~ x1 + x2 + (1 | area),
    data_model = d$data_model,
    data_proj = d$data_proj,
    domain = "domain",
    weight = "weight",
    return_direct = TRUE
  )

  expect_false(is.null(res$direct_estimator))
  expect_true(is.data.frame(res$direct_estimator))
  expect_true(all(c("estimate", "variance") %in% names(res$direct_estimator)))
})

test_that("print, summary, and as.data.frame methods work", {
  d <- make_test_data()

  res <- sae_ml_linear(
    formula = y ~ x1 + x2 + (1 | area),
    data_model = d$data_model,
    data_proj = d$data_proj,
    domain = "domain",
    weight = "weight"
  )

  expect_output(print(res), "SAE Projection Estimator")
  expect_output(summary(res), "Diagnostics")
  expect_true(is.data.frame(as.data.frame(res)))
  expect_identical(as.data.frame(res), as.data.frame(res$estimates))
})

test_that("returns bias-corrected estimator and synthetic components", {
  d <- make_test_data()

  res <- sae_ml_linear(
    formula = y ~ x1 + x2 + (1 | area),
    data_model = d$data_model,
    data_proj = d$data_proj,
    domain = "domain",
    weight = "weight"
  )

  expect_equal(res$estimator, "bias_corrected")
  expect_true(all(c(
    "estimate_synthetic",
    "variance_synthetic",
    "correction",
    "variance_correction",
    "estimate_final",
    "variance_final",
    "se_final",
    "rse_final",
    "n_model",
    "n_proj"
  ) %in% names(res$estimation_details)))
})

test_that("keep_unit = TRUE returns unit-level data", {
  d <- make_test_data()

  res <- sae_ml_linear(
    formula = y ~ x1 + x2 + (1 | area),
    data_model = d$data_model,
    data_proj = d$data_proj,
    domain = "domain",
    weight = "weight",
    keep_unit = TRUE
  )

  expect_false(is.null(res$unit_projection))
  expect_false(is.null(res$unit_model_residual))
  expect_true(".y_hat" %in% names(res$unit_projection))
  expect_true(".resid" %in% names(res$unit_model_residual))
})

test_that("cluster_ids = NULL is treated as no clustering", {
  d <- make_test_data()

  res <- sae_ml_linear(
    formula = y ~ x1 + x2 + (1 | area),
    data_model = d$data_model,
    data_proj = d$data_proj,
    domain = "domain",
    cluster_ids = NULL,
    weight = "weight"
  )

  expect_s3_class(res, "sae_ml_linear")
})

test_that("notes are condition-dependent and concise", {
  d <- make_test_data()

  res <- sae_ml_linear(
    formula = y ~ x1 + x2 + (1 | area),
    data_model = d$data_model,
    data_proj = d$data_proj,
    domain = "domain",
    weight = "weight"
  )

  expect_true(is.character(res$notes))
  expect_equal(res$estimator, "bias_corrected")
  expect_false(any(grepl("Bias-corrected estimator returned by default", res$notes)))
  expect_false(any(grepl("Plug-in variance from svyby", res$notes)))
})

test_that("diagnostics has expected components", {
  d <- make_test_data()

  res <- sae_ml_linear(
    formula = y ~ x1 + x2 + (1 | area),
    data_model = d$data_model,
    data_proj = d$data_proj,
    domain = "domain",
    weight = "weight"
  )

  expect_true(all(c(
    "icc",
    "singular_fit",
    "prediction_mode",
    "convergence",
    "sigma",
    "AIC",
    "BIC"
  ) %in% names(res$diagnostics)))

  expect_false("fixed_effects" %in% names(res$diagnostics))
  expect_false("convergence_messages" %in% names(res$diagnostics))
})

test_that("warns when projection domain has no model observations", {
  d <- make_test_data()

  dp_bad <- d$data_proj
  dp_bad$domain[dp_bad$area == "A7"] <- "D_new"

  expect_warning(
    res <- sae_ml_linear(
      formula = y ~ x1 + x2 + (1 | area),
      data_model = d$data_model,
      data_proj = dp_bad,
      domain = "domain",
      weight = "weight"
    ),
    regexp = "no residual correction"
  )

  expect_s3_class(res, "sae_ml_linear")
  expect_true(any(grepl("no residual correction", res$notes)))
})
