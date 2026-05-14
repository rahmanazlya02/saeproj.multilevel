# dev/simulasi.R — jalankan setelah devtools::load_all()
set.seed(42)

n_area  <- 10
n_model <- 200
n_proj  <- 2000

# ---- data_model ----
area_model <- sample(paste0("A", 1:n_area), n_model, replace = TRUE)
data_model <- data.frame(
  domain    = paste0("D", ceiling(as.integer(substr(area_model, 2, 3)) / 2)),
  area      = area_model,
  x1        = rnorm(n_model, 10, 2),
  x2        = runif(n_model),
  education = factor(sample(c("SD","SMP","SMA"), n_model, replace = TRUE,
                            prob = c(0.3,0.4,0.3)),
                     levels = c("SD","SMP","SMA")),
  weight    = runif(n_model, 1, 3)
)
area_eff       <- setNames(rnorm(n_area, 0, 2), paste0("A", 1:n_area))
data_model$y   <- with(data_model,
                       5 + 0.8*x1 + 3*x2 +
                         ifelse(education=="SMP", 1.5, ifelse(education=="SMA", 3.0, 0)) +
                         area_eff[area] + rnorm(n_model, 0, 1)
)

# ---- data_proj (boleh punya area baru: A11, A12) ----
area_proj <- sample(c(paste0("A", 1:n_area), "A11", "A12"), n_proj, replace = TRUE)
data_proj <- data.frame(
  domain    = ifelse(area_proj %in% c("A11","A12"), "D6",
                     paste0("D", ceiling(as.integer(substr(area_proj, 2, 3)) / 2))),
  area      = area_proj,
  x1        = rnorm(n_proj, 10, 2),
  x2        = runif(n_proj),
  education = factor(sample(c("SD","SMP","SMA"), n_proj, replace = TRUE,
                            prob = c(0.3,0.4,0.3)),
                     levels = c("SD","SMP","SMA")),
  weight    = runif(n_proj, 1, 3)
)

# ---- Jalankan fungsi ----
devtools::load_all()

res <- sae_ml_linear(
  formula     = y ~ x1 + x2 + education + (1 | area),
  data_model  = data_model,
  data_proj   = data_proj,
  domain      = "domain",
  cluster_ids = ~1,
  weight      = "weight",
  estimator   = "bias_corrected"
)

#print(res)  # mencetak final estimates ringkas
summary(res) # mencetak fitted model + diagnostics + final estimates + notes

cat("\n--- Estimation details ---\n")
print(res$estimation_details) #mencetak estimasi detail dengan tambahan variance correction

# ---- Numeric diagnostics ----
fit <- res$fitted_model

resid_model <- resid(fit)
fitted_model <- fitted(fit)
scaled_resid <- resid_model / sigma(fit)
ranef_area <- lme4::ranef(fit)$area[, 1]

diagnostic_numeric <- data.frame(
  mean_residual = mean(resid_model),
  sd_residual = sd(resid_model),
  sigma = sigma(fit),
  max_abs_scaled_residual = max(abs(scaled_resid)),
  prop_abs_scaled_resid_gt_2 = mean(abs(scaled_resid) > 2),
  prop_abs_scaled_resid_gt_3 = mean(abs(scaled_resid) > 3),
  cor_abs_resid_fitted = cor(abs(resid_model), fitted_model),
  shapiro_resid_p = shapiro.test(resid_model)$p.value,
  shapiro_ranef_p = shapiro.test(ranef_area)$p.value,
  singular_fit = lme4::isSingular(fit)
)
cat("\n--- Numeric diagnostics ---\n")
print(diagnostic_numeric, row.names = FALSE)

# ---- Diagnostic plots ----
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))

plot(fitted(fit), resid(fit),
     xlab = "Fitted values",
     ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, lty = 2)

qqnorm(resid(fit), main = "QQ Plot Residuals")
qqline(resid(fit))

hist(resid(fit), breaks = 20, main = "Residual Histogram",
     xlab = "Residuals")

ranef_area <- lme4::ranef(fit)$area[, 1]
qqnorm(ranef_area, main = "QQ Plot Random Effects")
qqline(ranef_area)

par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)
