# tests/testthat/test-sae-ml-linear.R

# Helper: data simulasi kecil untuk test
make_test_data <- function(seed = 123) {
  set.seed(seed)
  n_area  <- 6
  n_model <- 120
  n_proj  <- 600
  area_model <- sample(paste0("A", 1:n_area), n_model, replace = TRUE)
  data_model <- data.frame(
    domain    = paste0("D", ceiling(as.integer(substr(area_model, 2, 2)) / 2)),
    area      = area_model,
    x1        = rnorm(n_model),
    x2        = runif(n_model),
    education = factor(sample(c("SD","SMP","SMA"), n_model, replace = TRUE),
                       levels = c("SD","SMP","SMA")),
    weight    = runif(n_model, 1, 3)
  )
  area_eff       <- setNames(rnorm(n_area), paste0("A", 1:n_area))
  data_model$y   <- with(data_model,
                         5 + 0.5*x1 + 2*x2 + area_eff[area] + rnorm(n_model, 0, 0.5)
  )
  area_proj <- sample(c(paste0("A", 1:n_area), "A7"), n_proj, replace = TRUE)
  data_proj <- data.frame(
    domain    = ifelse(area_proj == "A7", "D4",
                       paste0("D", ceiling(as.integer(substr(area_proj, 2, 2)) / 2))),
    area      = area_proj,
    x1        = rnorm(n_proj),
    x2        = runif(n_proj),
    education = factor(sample(c("SD","SMP","SMA"), n_proj, replace = TRUE),
                       levels = c("SD","SMP","SMA")),
    weight    = runif(n_proj, 1, 3)
  )
  list(data_model = data_model, data_proj = data_proj)
}

# Test 1: Fungsi berjalan dengan random intercept
test_that("runs with random intercept", {
  d   <- make_test_data()
  res <- sae_ml_linear(
    formula    = y ~ x1 + x2 + (1 | area),
    data_model = d$data_model,
    data_proj  = d$data_proj,
    domain     = "domain",
    weight     = "weight"
  )
  expect_s3_class(res, "sae_ml_linear")
  expect_true(is.data.frame(res$estimates))
})

# Test 2: Fungsi berjalan dengan random slope
test_that("runs with random slope", {
  d <- make_test_data()
  expect_no_error(
    sae_ml_linear(
      formula    = y ~ x1 + x2 + (1 + x1 | area),
      data_model = d$data_model,
      data_proj  = d$data_proj,
      domain     = "domain",
      weight     = "weight"
    )
  )
})

# Test 3: data_proj boleh punya area baru
test_that("data_proj can have new random effect levels", {
  d <- make_test_data()
  expect_true("A7" %in% d$data_proj$area)
  res <- sae_ml_linear(
    formula    = y ~ x1 + x2 + (1 | area),
    data_model = d$data_model,
    data_proj  = d$data_proj,
    domain     = "domain",
    weight     = "weight"
  )
  expect_s3_class(res, "sae_ml_linear")
})

# Test 4: Error jika predictor tidak ada di data_proj
test_that("error if predictor missing in data_proj", {
  d      <- make_test_data()
  dp_bad <- d$data_proj
  dp_bad$x1 <- NULL
  expect_error(
    sae_ml_linear(
      formula    = y ~ x1 + x2 + (1 | area),
      data_model = d$data_model,
      data_proj  = dp_bad,
      domain     = "domain",
      weight     = "weight"
    ),
    regexp = "x1"
  )
})

# Test 5: Error jika categorical punya level baru di data_proj
test_that("error if categorical has new level in data_proj", {
  d      <- make_test_data()
  dp_bad <- d$data_proj
  dp_bad$education <- factor(
    sample(c("SD","SMP","SMA","S1"), nrow(dp_bad), replace = TRUE),
    levels = c("SD","SMP","SMA","S1")
  )
  expect_error(
    sae_ml_linear(
      formula    = y ~ x1 + x2 + education + (1 | area),
      data_model = d$data_model,
      data_proj  = dp_bad,
      domain     = "domain",
      weight     = "weight"
    ),
    regexp = "S1"
  )
})

# Test 6: Error jika ada missing value
test_that("error if missing values present", {
  d      <- make_test_data()
  dm_bad <- d$data_model
  dm_bad$x1[1] <- NA
  expect_error(
    sae_ml_linear(
      formula    = y ~ x1 + x2 + (1 | area),
      data_model = dm_bad,
      data_proj  = d$data_proj,
      domain     = "domain",
      weight     = "weight"
    ),
    regexp = "Missing values"
  )
})

# Test 7: Output punya kolom yang benar
test_that("estimates has required columns", {
  d   <- make_test_data()
  res <- sae_ml_linear(
    formula    = y ~ x1 + x2 + (1 | area),
    data_model = d$data_model,
    data_proj  = d$data_proj,
    domain     = "domain",
    weight     = "weight"
  )
  expect_true(all(c("estimate","variance","se","rse") %in% names(res$estimates)))
})

# Test 8: as.data.frame() mengembalikan data.frame
test_that("as.data.frame returns data.frame", {
  d   <- make_test_data()
  res <- sae_ml_linear(
    formula    = y ~ x1 + x2 + (1 | area),
    data_model = d$data_model,
    data_proj  = d$data_proj,
    domain     = "domain",
    weight     = "weight"
  )
  expect_true(is.data.frame(as.data.frame(res)))
})

# Test 9: estimator synthetic berjalan dan correction = 0
test_that("synthetic estimator sets correction to 0", {
  d   <- make_test_data()
  res <- sae_ml_linear(
    formula    = y ~ x1 + x2 + (1 | area),
    data_model = d$data_model,
    data_proj  = d$data_proj,
    domain     = "domain",
    weight     = "weight",
    estimator  = "synthetic"
  )
  expect_equal(res$estimator, "synthetic")
  expect_true(all(res$estimates$correction == 0))
})

# Test 10: keep_unit = TRUE menyimpan unit data
test_that("keep_unit = TRUE returns unit data", {
  d   <- make_test_data()
  res <- sae_ml_linear(
    formula    = y ~ x1 + x2 + (1 | area),
    data_model = d$data_model,
    data_proj  = d$data_proj,
    domain     = "domain",
    weight     = "weight",
    keep_unit  = TRUE
  )
  expect_false(is.null(res$unit_projection))
  expect_false(is.null(res$unit_model_residual))
})
