# data-raw/saeml_data.R

library(dplyr)

# -------------------------------------------------------------------------
# 1. Generate the fixed population
# -------------------------------------------------------------------------

set.seed(999)

# Population size
n_kabupaten <- 50
n_per_kab <- 1000
N <- n_kabupaten * n_per_kab

# Area-level model parameters
a0 <- 100
a1 <- 0.2
a2 <- -0.2

# Unit-level model parameters
b1 <- 25
b2 <- 18
b3 <- -22
b4 <- 20

# Variance parameters
tau0_2 <- 900
sigma2 <- 80

# Area-level auxiliary variables
Z1 <- rnorm(n_kabupaten)
Z2 <- rnorm(n_kabupaten)

# Random intercept for each district/city
u0 <- rnorm(
  n_kabupaten,
  mean = 0,
  sd = sqrt(tau0_2)
)

# Intercept for each district/city
b0 <- a0 + a1 * Z1 + a2 * Z2 + u0

# Unit-level variables
kab_kota <- rep(1:n_kabupaten, each = n_per_kab)
id_individu <- 1:N

X1 <- rnorm(N)
X2 <- rbinom(N, size = 1, prob = 0.5)
X3 <- rnorm(N)
X4 <- runif(N, min = -2, max = 2)

# Unit-level residual
eij <- rnorm(
  N,
  mean = 0,
  sd = sqrt(sigma2)
)

# Target variable
Y <- rep(b0, each = n_per_kab) +
  b1 * X1 +
  b2 * X2 +
  b3 * X3 +
  b4 * X4 +
  eij

# Fixed population
population <- data.frame(
  prov = "35",
  kab_kota = as.character(kab_kota),
  id_individu = id_individu,
  Z1 = rep(Z1, each = n_per_kab),
  Z2 = rep(Z2, each = n_per_kab),
  X1 = X1,
  X2 = X2,
  X3 = X3,
  X4 = X4,
  Y = Y
)

# -------------------------------------------------------------------------
# 2. Draw one fixed sample replication
# -------------------------------------------------------------------------

# Small model survey: 0.5% = 5 units per district/city
# Large projection survey: 30% = 300 units per district/city
sample_small_prop <- 0.005
sample_large_prop <- 0.30

n_model_per_area <- n_per_kab * sample_small_prop
n_proj_per_area <- n_per_kab * sample_large_prop

# Use the first fixed simulation replication
set.seed(3001)

# Randomize the unit order within each district/city
sample_order <- population %>%
  group_by(kab_kota) %>%
  slice_sample(prop = 1) %>%
  mutate(sample_order = row_number()) %>%
  ungroup()

# -------------------------------------------------------------------------
# 3. Create the small model survey
# -------------------------------------------------------------------------

saeml_modelsvy <- sample_order %>%
  filter(sample_order <= n_model_per_area) %>%
  group_by(kab_kota) %>%
  mutate(WEIND = n_per_kab / n()) %>%
  ungroup() %>%
  select(
    prov,
    kab_kota,
    id_individu,
    Z1,
    Z2,
    X1,
    X2,
    X3,
    X4,
    Y,
    WEIND
  ) %>%
  arrange(as.numeric(kab_kota), id_individu)

# -------------------------------------------------------------------------
# 4. Create the large projection survey
# -------------------------------------------------------------------------

saeml_projsvy <- sample_order %>%
  filter(sample_order > n_per_kab - n_proj_per_area) %>%
  group_by(kab_kota) %>%
  mutate(WEIND = n_per_kab / n()) %>%
  ungroup() %>%
  select(
    prov,
    kab_kota,
    id_individu,
    Z1,
    Z2,
    X1,
    X2,
    X3,
    X4,
    WEIND
  ) %>%
  arrange(as.numeric(kab_kota), id_individu)

# -------------------------------------------------------------------------
# 5. Check the generated datasets
# -------------------------------------------------------------------------

n_overlap <- sum(
  saeml_modelsvy$id_individu %in% saeml_projsvy$id_individu
)

stopifnot(
  nrow(saeml_modelsvy) == 250,
  nrow(saeml_projsvy) == 15000,
  length(unique(saeml_modelsvy$kab_kota)) == 50,
  length(unique(saeml_projsvy$kab_kota)) == 50,
  all(table(saeml_modelsvy$kab_kota) == 5),
  all(table(saeml_projsvy$kab_kota) == 300),
  n_overlap == 0,
  "Y" %in% names(saeml_modelsvy),
  !("Y" %in% names(saeml_projsvy)),
  all(complete.cases(saeml_modelsvy)),
  all(complete.cases(saeml_projsvy))
)

# -------------------------------------------------------------------------
# 6. Save the datasets in the package
# -------------------------------------------------------------------------

usethis::use_data(
  saeml_modelsvy,
  saeml_projsvy,
  compress = "xz",
  overwrite = TRUE
)
