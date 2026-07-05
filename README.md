
# saeproj.multilevel

## Author

Nazlya Rahma Susanto, Azka Ubaidillah

## Maintainer

Nazlya Rahma Susanto <susantonazlya@gmail.com>

## Description

The **saeproj.multilevel** package provides tools for *Small Area
Estimation* (SAE) using a projection estimator with a multilevel
regression model.

The method is designed for two-survey settings:

- a smaller model survey, `data_model`, that contains the response
  variable and auxiliary predictors;
- a larger projection survey, `data_proj`, that contains auxiliary
  predictors and survey design information, but does not contain the
  response variable.

The main function is:

``` r
sae_ml_linear()
```

The function fits a linear multilevel regression model using
`lme4::lmer()`, generates unit-level predictions for the projection
dataset, aggregates those predictions by domain using survey design
information, and applies a design-based residual correction.

The final projection estimator is:

``` r
estimate_final = estimate_synthetic + correction
```

The plug-in variance is calculated as:

``` r
variance_final = variance_synthetic + variance_correction
```

The synthetic projection component and residual correction component are
stored in:

``` r
result$estimation_details
```

## Installation

The development version of `saeproj.multilevel` can be installed from
GitHub with:

``` r
# install.packages("devtools")
devtools::install_github("rahmanazlya02/saeproj.multilevel")
```

To install the package together with the vignette, use:

``` r
# install.packages("devtools")
devtools::install_github(
  "rahmanazlya02/saeproj.multilevel",
  build_vignettes = TRUE,
  dependencies = TRUE
)
```

After installation, the vignette can be opened with:

``` r
browseVignettes("saeproj.multilevel")
```

Or directly:

``` r
vignette(
  "sae_ml_linear",
  package = "saeproj.multilevel"
)
```

## Dependencies

The package imports:

- `lme4` — for fitting linear multilevel regression models;
- `survey` — for survey design and domain-level aggregation;
- `dplyr` — for joining estimation components;
- `cli` — for errors and selected warnings;
- `reformulas` — for parsing multilevel model formulas.

## Package datasets

The package includes two simulated datasets generated from one fixed
replication of the study-simulation design.

``` r
data("saeml_modelsvy")
data("saeml_projsvy")
```

### `saeml_modelsvy`

`saeml_modelsvy` is a small model-survey dataset containing:

- 250 observations;
- 50 domains identified by `kab_kota`;
- 5 sampled units in each domain;
- the target variable `Y`;
- unit-level auxiliary variables `X1`, `X2`, `X3`, and `X4`;
- area-level auxiliary variables `Z1` and `Z2`;
- the survey weight variable `WEIND`;
- no separate PSU or cluster variable.

### `saeml_projsvy`

`saeml_projsvy` is a large projection-survey dataset containing:

- 15,000 observations;
- 50 domains identified by `kab_kota`;
- 300 sampled units in each domain;
- unit-level auxiliary variables `X1`, `X2`, `X3`, and `X4`;
- area-level auxiliary variables `Z1` and `Z2`;
- the survey weight variable `WEIND`;
- no target variable `Y`;
- no separate PSU or cluster variable.

The two datasets are drawn from the same simulated population and do not
contain overlapping sampled units.

``` r
dim(saeml_modelsvy)
#> [1] 250  11

dim(saeml_projsvy)
#> [1] 15000    10
```

## Example

### Fit the multilevel projection estimator

``` r
result <- sae_ml_linear(
  formula = Y ~ X1 + X2 + X3 + X4 + Z1 + Z2 + (1 | kab_kota),
  data_model = saeml_modelsvy,
  data_proj = saeml_projsvy,
  domain = "kab_kota",
  cluster_ids = ~1,
  weight = "WEIND",
  strata = "kab_kota",
  summary_function = "mean"
)

result
#> SAE Projection Estimator using Linear Multilevel Model
#> -------------------------------------------------------
#> Formula   : Y ~ X1 + X2 + X3 + X4 + Z1 + Z2 + (1 | kab_kota) 
#> Estimator : bias_corrected 
#> Domains   : 50 
#> 
#> Estimates:
#>  kab_kota  estimate variance       se       rse
#>         1  63.63811 28.51940 5.340356  8.391758
#>         2 123.57033 22.92302 4.787799  3.874554
#>         3  72.21099 24.03748 4.902803  6.789553
#>         4  89.15406 25.01544 5.001543  5.610001
#>         5 160.68935 12.41104 3.522931  2.192386
#>         6  27.48805 28.16499 5.307070 19.306825
```

The package datasets do not contain a separate PSU or cluster variable.
Therefore, the example uses:

``` r
cluster_ids = ~1
```

This specifies an unclustered survey-design structure. The variable
`id_individu` is only a unique sampled-unit identifier and is not used
as a PSU or cluster identifier.

### Domain-level estimates

The final domain-level estimates are stored in:

``` r
result$estimates
```

The complete results for all 50 domains are shown below.

``` r
result$estimates
#>    kab_kota  estimate  variance       se       rse
#> 1         1  63.63811 28.519405 5.340356  8.391758
#> 2         2 123.57033 22.923017 4.787799  3.874554
#> 3         3  72.21099 24.037478 4.902803  6.789553
#> 4         4  89.15406 25.015436 5.001543  5.610001
#> 5         5 160.68935 12.411044 3.522931  2.192386
#> 6         6  27.48805 28.164991 5.307070 19.306825
#> 7         7 118.01107 21.596657 4.647220  3.937953
#> 8         8 154.32727 20.918597 4.573685  2.963627
#> 9         9  66.40287 19.156260 4.376787  6.591261
#> 10       10  89.89285 26.197321 5.118332  5.693814
#> 11       11  93.40441 18.723591 4.327077  4.632625
#> 12       12  67.82925 10.371817 3.220531  4.747996
#> 13       13  86.94722 26.877987 5.184398  5.962696
#> 14       14  83.26791 24.559496 4.955754  5.951577
#> 15       15 104.31824 28.927284 5.378409  5.155771
#> 16       16  66.99395 14.566786 3.816646  5.697001
#> 17       17 138.20477 44.228946 6.650485  4.812051
#> 18       18 146.17343 12.079225 3.475518  2.377667
#> 19       19  69.20977 30.322787 5.506613  7.956411
#> 20       20 126.71531 21.429139 4.629162  3.653198
#> 21       21  91.89536  8.232509 2.869235  3.122285
#> 22       22 112.26391 35.805903 5.983803  5.330122
#> 23       23  38.37135 13.554550 3.681650  9.594791
#> 24       24  54.02795 28.115349 5.302391  9.814163
#> 25       25 146.38489 41.543110 6.445395  4.403046
#> 26       26 119.26013 19.467277 4.412174  3.699622
#> 27       27 122.18386 30.941331 5.562493  4.552560
#> 28       28 114.81818 30.886520 5.557564  4.840317
#> 29       29 140.21555 27.966702 5.288355  3.771590
#> 30       30 114.57995 22.802072 4.775151  4.167528
#> 31       31  91.52510 43.345164 6.583704  7.193332
#> 32       32 102.80384 26.687399 5.165985  5.025090
#> 33       33 117.15727  9.750002 3.122499  2.665220
#> 34       34  81.85505 20.017747 4.474120  5.465906
#> 35       35  81.37672  9.574114 3.094207  3.802324
#> 36       36 121.88129  7.630933 2.762414  2.266479
#> 37       37  66.20806 43.124893 6.566955  9.918664
#> 38       38  90.88219 28.563450 5.344478  5.880666
#> 39       39  87.61008 17.579341 4.192772  4.785719
#> 40       40 149.73935 36.185367 6.015427  4.017266
#> 41       41 108.77503 12.750292 3.570755  3.282697
#> 42       42 154.02090 22.883156 4.783634  3.105834
#> 43       43 106.91171  9.739403 3.120802  2.919046
#> 44       44 142.66619 18.383221 4.287566  3.005313
#> 45       45 125.47587 52.515914 7.246786  5.775442
#> 46       46  94.25971 17.091366 4.134171  4.385936
#> 47       47 134.25627 36.781517 6.064777  4.517314
#> 48       48 116.44011 37.351694 6.111603  5.248710
#> 49       49  80.50337 24.459513 4.945656  6.143415
#> 50       50 133.66421 17.656474 4.201961  3.143669
```

The output contains:

| Column | Description |
|----|----|
| domain variable(s) | Domain identifier column(s), based on the `domain` argument |
| `estimate` | Final projection estimate with design-based residual correction |
| `variance` | Plug-in variance of the final estimate |
| `se` | Standard error, computed as `sqrt(variance)` |
| `rse` | Relative standard error in percent |

The same result can be extracted for further analysis with:

``` r
as.data.frame(result)
```

## Estimation components

Detailed estimation components are stored in:

``` r
result$estimation_details
```

The complete synthetic estimate, residual correction, final estimate,
variance, and sample-size information for all 50 domains are shown
below.

``` r
result$estimation_details
#>    kab_kota estimate_synthetic variance_synthetic   correction
#> 1         1           64.51902           6.651059 -0.880903021
#> 2         2          123.47261           6.399579  0.097720077
#> 3         3           72.45427           5.914642 -0.243283358
#> 4         4           89.31018           6.307633 -0.156121828
#> 5         5          159.60585           6.707629  1.083499650
#> 6         6           29.12853           5.924376 -1.640481978
#> 7         7          117.95382           6.655086  0.057244068
#> 8         8          153.70964           5.881752  0.617629498
#> 9         9           67.29964           6.336749 -0.896767634
#> 10       10           90.47768           6.723465 -0.584833009
#> 11       11           93.35554           5.523737  0.048871257
#> 12       12           68.38506           5.808774 -0.555804610
#> 13       13           86.95638           5.906556 -0.009161651
#> 14       14           83.46482           6.117049 -0.196910212
#> 15       15          104.05296           5.361198  0.265282500
#> 16       16           68.09623           6.211223 -1.102279167
#> 17       17          137.62282           6.025534  0.581955378
#> 18       18          145.28490           5.861022  0.888529509
#> 19       19           69.77688           6.319282 -0.567113232
#> 20       20          126.46963           5.735366  0.245679453
#> 21       21           92.32333           5.985796 -0.427971400
#> 22       22          111.93106           6.367506  0.332843695
#> 23       23           39.62914           6.675530 -1.257790218
#> 24       24           54.97937           6.246625 -0.951418482
#> 25       25          145.57317           5.810129  0.811718138
#> 26       26          118.73939           5.431806  0.520735814
#> 27       27          121.93888           5.632564  0.244978471
#> 28       28          114.92443           5.067101 -0.106249134
#> 29       29          139.62463           6.424310  0.590921202
#> 30       30          114.38188           5.972599  0.198068985
#> 31       31           91.24777           5.591878  0.277322035
#> 32       32          102.51121           5.524387  0.292626541
#> 33       33          116.89962           6.666377  0.257653811
#> 34       34           82.00056           6.491236 -0.145510938
#> 35       35           82.03037           5.283640 -0.653657212
#> 36       36          121.87937           5.405497  0.001919269
#> 37       37           66.79017           6.181848 -0.582110617
#> 38       38           91.12162           5.891640 -0.239429982
#> 39       39           87.72478           6.084069 -0.114693188
#> 40       40          149.04234           5.802844  0.697014744
#> 41       41          108.43554           6.201553  0.339490042
#> 42       42          152.87854           5.714695  1.142356288
#> 43       43          107.17430           5.709964 -0.262587216
#> 44       44          141.97962           5.618465  0.686567860
#> 45       45          125.08834           7.204382  0.387531328
#> 46       46           94.20918           5.228550  0.050531652
#> 47       47          134.03805           6.232447  0.218213781
#> 48       48          116.21694           5.058848  0.223165144
#> 49       49           80.82393           5.787536 -0.320562289
#> 50       50          132.92864           7.120847  0.735570186
#>    variance_correction estimate_final variance_final se_final rse_final n_model
#> 1            21.868346       63.63811      28.519405 5.340356  8.391758       5
#> 2            16.523438      123.57033      22.923017 4.787799  3.874554       5
#> 3            18.122836       72.21099      24.037478 4.902803  6.789553       5
#> 4            18.707803       89.15406      25.015436 5.001543  5.610001       5
#> 5             5.703415      160.68935      12.411044 3.522931  2.192386       5
#> 6            22.240615       27.48805      28.164991 5.307070 19.306825       5
#> 7            14.941571      118.01107      21.596657 4.647220  3.937953       5
#> 8            15.036845      154.32727      20.918597 4.573685  2.963627       5
#> 9            12.819511       66.40287      19.156260 4.376787  6.591261       5
#> 10           19.473856       89.89285      26.197321 5.118332  5.693814       5
#> 11           13.199854       93.40441      18.723591 4.327077  4.632625       5
#> 12            4.563044       67.82925      10.371817 3.220531  4.747996       5
#> 13           20.971432       86.94722      26.877987 5.184398  5.962696       5
#> 14           18.442446       83.26791      24.559496 4.955754  5.951577       5
#> 15           23.566086      104.31824      28.927284 5.378409  5.155771       5
#> 16            8.355563       66.99395      14.566786 3.816646  5.697001       5
#> 17           38.203412      138.20477      44.228946 6.650485  4.812051       5
#> 18            6.218204      146.17343      12.079225 3.475518  2.377667       5
#> 19           24.003505       69.20977      30.322787 5.506613  7.956411       5
#> 20           15.693773      126.71531      21.429139 4.629162  3.653198       5
#> 21            2.246713       91.89536       8.232509 2.869235  3.122285       5
#> 22           29.438396      112.26391      35.805903 5.983803  5.330122       5
#> 23            6.879021       38.37135      13.554550 3.681650  9.594791       5
#> 24           21.868725       54.02795      28.115349 5.302391  9.814163       5
#> 25           35.732982      146.38489      41.543110 6.445395  4.403046       5
#> 26           14.035471      119.26013      19.467277 4.412174  3.699622       5
#> 27           25.308766      122.18386      30.941331 5.562493  4.552560       5
#> 28           25.819418      114.81818      30.886520 5.557564  4.840317       5
#> 29           21.542393      140.21555      27.966702 5.288355  3.771590       5
#> 30           16.829473      114.57995      22.802072 4.775151  4.167528       5
#> 31           37.753286       91.52510      43.345164 6.583704  7.193332       5
#> 32           21.163012      102.80384      26.687399 5.165985  5.025090       5
#> 33            3.083625      117.15727       9.750002 3.122499  2.665220       5
#> 34           13.526511       81.85505      20.017747 4.474120  5.465906       5
#> 35            4.290474       81.37672       9.574114 3.094207  3.802324       5
#> 36            2.225436      121.88129       7.630933 2.762414  2.266479       5
#> 37           36.943044       66.20806      43.124893 6.566955  9.918664       5
#> 38           22.671810       90.88219      28.563450 5.344478  5.880666       5
#> 39           11.495273       87.61008      17.579341 4.192772  4.785719       5
#> 40           30.382523      149.73935      36.185367 6.015427  4.017266       5
#> 41            6.548739      108.77503      12.750292 3.570755  3.282697       5
#> 42           17.168461      154.02090      22.883156 4.783634  3.105834       5
#> 43            4.029439      106.91171       9.739403 3.120802  2.919046       5
#> 44           12.764755      142.66619      18.383221 4.287566  3.005313       5
#> 45           45.311531      125.47587      52.515914 7.246786  5.775442       5
#> 46           11.862816       94.25971      17.091366 4.134171  4.385936       5
#> 47           30.549069      134.25627      36.781517 6.064777  4.517314       5
#> 48           32.292847      116.44011      37.351694 6.111603  5.248710       5
#> 49           18.671977       80.50337      24.459513 4.945656  6.143415       5
#> 50           10.535627      133.66421      17.656474 4.201961  3.143669       5
#>    n_proj
#> 1     300
#> 2     300
#> 3     300
#> 4     300
#> 5     300
#> 6     300
#> 7     300
#> 8     300
#> 9     300
#> 10    300
#> 11    300
#> 12    300
#> 13    300
#> 14    300
#> 15    300
#> 16    300
#> 17    300
#> 18    300
#> 19    300
#> 20    300
#> 21    300
#> 22    300
#> 23    300
#> 24    300
#> 25    300
#> 26    300
#> 27    300
#> 28    300
#> 29    300
#> 30    300
#> 31    300
#> 32    300
#> 33    300
#> 34    300
#> 35    300
#> 36    300
#> 37    300
#> 38    300
#> 39    300
#> 40    300
#> 41    300
#> 42    300
#> 43    300
#> 44    300
#> 45    300
#> 46    300
#> 47    300
#> 48    300
#> 49    300
#> 50    300
```

This table contains:

| Column | Description |
|----|----|
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

### Synthetic projection component

The function returns the projection estimator with a design-based
residual correction by default.

The complete synthetic component for all 50 domains is available below.

``` r
result$estimation_details[, c(
  "kab_kota",
  "estimate_synthetic",
  "variance_synthetic"
)]
#>    kab_kota estimate_synthetic variance_synthetic
#> 1         1           64.51902           6.651059
#> 2         2          123.47261           6.399579
#> 3         3           72.45427           5.914642
#> 4         4           89.31018           6.307633
#> 5         5          159.60585           6.707629
#> 6         6           29.12853           5.924376
#> 7         7          117.95382           6.655086
#> 8         8          153.70964           5.881752
#> 9         9           67.29964           6.336749
#> 10       10           90.47768           6.723465
#> 11       11           93.35554           5.523737
#> 12       12           68.38506           5.808774
#> 13       13           86.95638           5.906556
#> 14       14           83.46482           6.117049
#> 15       15          104.05296           5.361198
#> 16       16           68.09623           6.211223
#> 17       17          137.62282           6.025534
#> 18       18          145.28490           5.861022
#> 19       19           69.77688           6.319282
#> 20       20          126.46963           5.735366
#> 21       21           92.32333           5.985796
#> 22       22          111.93106           6.367506
#> 23       23           39.62914           6.675530
#> 24       24           54.97937           6.246625
#> 25       25          145.57317           5.810129
#> 26       26          118.73939           5.431806
#> 27       27          121.93888           5.632564
#> 28       28          114.92443           5.067101
#> 29       29          139.62463           6.424310
#> 30       30          114.38188           5.972599
#> 31       31           91.24777           5.591878
#> 32       32          102.51121           5.524387
#> 33       33          116.89962           6.666377
#> 34       34           82.00056           6.491236
#> 35       35           82.03037           5.283640
#> 36       36          121.87937           5.405497
#> 37       37           66.79017           6.181848
#> 38       38           91.12162           5.891640
#> 39       39           87.72478           6.084069
#> 40       40          149.04234           5.802844
#> 41       41          108.43554           6.201553
#> 42       42          152.87854           5.714695
#> 43       43          107.17430           5.709964
#> 44       44          141.97962           5.618465
#> 45       45          125.08834           7.204382
#> 46       46           94.20918           5.228550
#> 47       47          134.03805           6.232447
#> 48       48          116.21694           5.058848
#> 49       49           80.82393           5.787536
#> 50       50          132.92864           7.120847
```

## Direct estimator

Set `return_direct = TRUE` to return direct design-based estimates from
`data_model`.

``` r
result_direct <- sae_ml_linear(
  formula = Y ~ X1 + X2 + X3 + X4 + Z1 + Z2 + (1 | kab_kota),
  data_model = saeml_modelsvy,
  data_proj = saeml_projsvy,
  domain = "kab_kota",
  cluster_ids = ~1,
  weight = "WEIND",
  strata = "kab_kota",
  summary_function = "mean",
  return_direct = TRUE
)

result_direct$direct_estimator
```

The direct estimator is stored separately and does not replace the
projection estimator.

## Print and summary

A concise output can be displayed with:

``` r
print(result)
#> SAE Projection Estimator using Linear Multilevel Model
#> -------------------------------------------------------
#> Formula   : Y ~ X1 + X2 + X3 + X4 + Z1 + Z2 + (1 | kab_kota) 
#> Estimator : bias_corrected 
#> Domains   : 50 
#> 
#> Estimates:
#>  kab_kota  estimate variance       se       rse
#>         1  63.63811 28.51940 5.340356  8.391758
#>         2 123.57033 22.92302 4.787799  3.874554
#>         3  72.21099 24.03748 4.902803  6.789553
#>         4  89.15406 25.01544 5.001543  5.610001
#>         5 160.68935 12.41104 3.522931  2.192386
#>         6  27.48805 28.16499 5.307070 19.306825
```

A compact summary can be displayed with:

``` r
summary(result)
#> SAE Projection Estimator using Linear Multilevel Model
#> -------------------------------------------------------
#> Formula   : Y ~ X1 + X2 + X3 + X4 + Z1 + Z2 + (1 | kab_kota) 
#> Estimator : bias_corrected 
#> Domains   : 50 
#> 
#> Model diagnostics:
#>   nobs        : 250 
#>   sigma       : 9.6444 
#>   ICC         : 0.9048 
#>   singular    : FALSE 
#>   convergence : OK 
#> 
#> Estimates:
#>  kab_kota  estimate variance       se       rse
#>         1  63.63811 28.51940 5.340356  8.391758
#>         2 123.57033 22.92302 4.787799  3.874554
#>         3  72.21099 24.03748 4.902803  6.789553
#>         4  89.15406 25.01544 5.001543  5.610001
#>         5 160.68935 12.41104 3.522931  2.192386
#>         6  27.48805 28.16499 5.307070 19.306825
```

The `summary()` method displays the formula, estimator type, number of
domains, selected model diagnostics, and a preview of the final
estimates.

Full model output can be accessed from the fitted `lmerMod` object:

``` r
fit <- result$fitted_model

summary(fit)
#> Linear mixed model fit by REML ['lmerMod']
#> Formula: Y ~ X1 + X2 + X3 + X4 + Z1 + Z2 + (1 | kab_kota)
#>    Data: data
#> Control: control
#> 
#> REML criterion at convergence: 2009.9
#> 
#> Scaled residuals: 
#>      Min       1Q   Median       3Q      Max 
#> -2.16963 -0.59140  0.04972  0.52613  2.44401 
#> 
#> Random effects:
#>  Groups   Name        Variance Std.Dev.
#>  kab_kota (Intercept) 884.17   29.735  
#>  Residual              93.01    9.644  
#> Number of obs: 250, groups:  kab_kota, 50
#> 
#> Fixed effects:
#>             Estimate Std. Error t value
#> (Intercept)  92.4309     4.4555  20.745
#> X1           25.1077     0.6709  37.422
#> X2           18.8027     1.3788  13.637
#> X3          -22.5309     0.6800 -33.133
#> X4           20.4876     0.5781  35.439
#> Z1           -9.6485     4.3391  -2.224
#> Z2           -6.1241     4.8901  -1.252
#> 
#> Correlation of Fixed Effects:
#>    (Intr) X1     X2     X3     X4     Z1    
#> X1 -0.003                                   
#> X2 -0.155 -0.023                            
#> X3 -0.019 -0.007  0.010                     
#> X4  0.005  0.003  0.051  0.109              
#> Z1  0.253 -0.018 -0.004 -0.007  0.006       
#> Z2 -0.054 -0.017  0.003 -0.012  0.024 -0.024
```

## Retaining unit-level predictions

Set `keep_unit = TRUE` to store unit-level projection data and model
residual data.

``` r
result_ku <- sae_ml_linear(
  formula = Y ~ X1 + X2 + X3 + X4 + Z1 + Z2 + (1 | kab_kota),
  data_model = saeml_modelsvy,
  data_proj = saeml_projsvy,
  domain = "kab_kota",
  cluster_ids = ~1,
  weight = "WEIND",
  strata = "kab_kota",
  summary_function = "mean",
  keep_unit = TRUE
)

head(result_ku$unit_projection)
head(result_ku$unit_model_residual)
```

When `keep_unit = TRUE`:

- `result_ku$unit_projection` contains `data_proj` with the unit-level
  prediction column `.prediction`;
- `result_ku$unit_model_residual` contains `data_model` with
  `.fitted_model` and `.model_residual`.

## Model diagnostics

Model diagnostics are stored in:

``` r
result$diagnostics
```

``` r
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
#>         icc singular_fit convergence    sigma residual_variance REML      AIC
#> 1 0.9048142        FALSE          OK 9.644356           93.0136 TRUE 2027.936
#>        BIC
#> 1 2059.629
```

The estimated random effects for all domain groups can be inspected
directly:

``` r
lme4::ranef(result$fitted_model)$kab_kota
#>     (Intercept)
#> 1  -41.86827839
#> 2    4.64451965
#> 3  -11.56297019
#> 4   -7.42028578
#> 5   51.49745648
#> 6  -77.97016757
#> 7    2.72074285
#> 8   29.35519933
#> 9  -42.62230470
#> 10 -27.79642104
#> 11   2.32279301
#> 12 -26.41673557
#> 13  -0.43544244
#> 14  -9.35890943
#> 15  12.60856337
#> 16 -52.39002472
#> 17  27.65965064
#> 18  42.23075634
#> 19 -26.95422097
#> 20  11.67685373
#> 21 -20.34097431
#> 22  15.81966705
#> 23 -59.78128098
#> 24 -45.21979483
#> 25  38.58000279
#> 26  24.74995717
#> 27  11.64353691
#> 28  -5.04989563
#> 29  28.08578562
#> 30   9.41398453
#> 31  13.18078825
#> 32  13.90819327
#> 33  12.24598079
#> 34  -6.91596274
#> 35 -31.06755394
#> 36   0.09122059
#> 37 -27.66702893
#> 38 -11.37982376
#> 39  -5.45123152
#> 40  33.12828614
#> 41  16.13556004
#> 42  54.29484285
#> 43 -12.48045973
#> 44  32.63175812
#> 45  18.41890554
#> 46   2.40170965
#> 47  10.37144283
#> 48  10.60677527
#> 49 -15.23594635
#> 50  34.96078070
```

Residual diagnostics can be inspected from the fitted model:

``` r
fit <- result$fitted_model

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
```

## Model parameters

Estimated model parameters are stored in:

``` r
result$model_parameters
```

``` r
result$model_parameters$fixed_effects
#> (Intercept)          X1          X2          X3          X4          Z1 
#>   92.430911   25.107730   18.802686  -22.530931   20.487619   -9.648481 
#>          Z2 
#>   -6.124119

result$model_parameters$variance_components
#>        grp        var1 var2     vcov     sdcor
#> 1 kab_kota (Intercept) <NA> 884.1652 29.734916
#> 2 Residual        <NA> <NA>  93.0136  9.644356

result$model_parameters$residual_variance
#> [1] 93.0136
```

## Notes

Run-specific notes are stored in:

``` r
result$notes
#> character(0)
```

The notes are intentionally concise and are not printed automatically by
`summary()`.

They may include information such as:

- removed zero-variance predictors;
- singular model fit;
- convergence issues;
- random-slope or complex random-effect structures where a simple ICC is
  not computed;
- out-of-sample domains with zero residual correction;
- negative plug-in variance clamped to zero.

Out-of-sample domains are not treated as warnings because they are
expected in SAE projection. They are recorded in `result$notes`.

## Multiple domain variables

The `domain` argument accepts a character scalar, a character vector, or
a one-sided formula.

The following example uses both `prov` and `kab_kota` as domain
identifiers.

``` r
result_multi <- sae_ml_linear(
  formula = Y ~ X1 + X2 + X3 + X4 + Z1 + Z2 + (1 | kab_kota),
  data_model = saeml_modelsvy,
  data_proj = saeml_projsvy,
  domain = c("prov", "kab_kota"),
  cluster_ids = ~1,
  weight = "WEIND",
  strata = "kab_kota",
  summary_function = "mean"
)

result_multi$estimates
```

## Survey design specification

The arguments `cluster_ids`, `weight`, and `strata` are used in the
aggregation step through `survey::svydesign()`.

### Simulated package data

The simulated datasets included in the package do not contain a separate
PSU or cluster variable. Therefore, the package examples use:

``` r
cluster_ids = ~1
weight = "WEIND"
strata = "kab_kota"
```

Here, `cluster_ids = ~1` specifies an unclustered survey-design
structure.

### Survey data with PSU clustering

For a real survey with a PSU or cluster variable, provide the actual PSU
identifier in `cluster_ids`.

The following code is illustrative. Replace `psu_id`, `survey_weight`,
and `stratum` with the corresponding variable names in your data.

``` r
result_clustered <- sae_ml_linear(
  formula = Y ~ X1 + X2 + X3 + X4 + Z1 + Z2 + (1 | kab_kota),
  data_model = data_model,
  data_proj = data_proj,
  domain = "kab_kota",
  cluster_ids = "psu_id",
  weight = "survey_weight",
  strata = "stratum",
  summary_function = "mean",
  nest = TRUE
)
```

In this specification:

- `psu_id` identifies the primary sampling unit or cluster;
- `survey_weight` identifies the sampling weight;
- `stratum` identifies the sampling stratum;
- `nest = TRUE` indicates that PSUs are nested within strata.

Use `cluster_ids = ~1` when the survey design does not include a
separate PSU or cluster variable.

## Output object structure

`sae_ml_linear()` returns an S3 object of class `"sae_ml_linear"`.

Typical components are:

| Component | Description |
|----|----|
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
|----|----|
| `print(result)` | Prints formula, estimator, number of domains, and a preview of `$estimates` |
| `summary(result)` | Prints selected diagnostics and a preview of final estimates |
| `as.data.frame(result)` | Returns `result$estimates` |

## Function interface

``` r
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
|----|----|
| `formula` | `lme4::lmer()`-style formula containing at least one random-effect term |
| `data_model` | Model survey data frame containing the response, predictors, grouping variables, domain variable(s), and survey design variables |
| `data_proj` | Projection survey data frame containing predictors, grouping variables, domain variable(s), and survey design variables; the response is not required |
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

The `weight` argument identifies the survey weight column used in both
`data_model` and `data_proj`. The column name must be the same in both
datasets, but the weight values may differ.

In `data_model`, weights are used for residual correction and optional
direct estimation. In `data_proj`, weights are used for synthetic
projection aggregation.

## Methodological notes

- The working model is a linear multilevel regression model fitted with
  `lme4::lmer()` using restricted maximum likelihood estimation.
- The user fully specifies the fixed-effect and random-effect structure
  through the `formula` argument.
- Prediction uses `re.form = NULL` and `allow.new.levels = TRUE`.
- For grouping levels observed in `data_model`, predictions include the
  estimated random-effect contribution.
- For grouping levels appearing only in `data_proj`, the random-effect
  contribution is set to zero, so prediction uses the fixed part of the
  model.
- ICC is a diagnostic measure and does not determine the prediction
  rule.
- A simple ICC is computed only for pure random-intercept structures.
- Fixed-effect categorical predictors in `data_proj` must not contain
  levels that are absent from `data_model`.
- Fixed-effect predictors with zero variance in `data_model` are removed
  automatically before model fitting.
- Survey design arguments (`cluster_ids`, `weight`, and `strata`) are
  used in the aggregation step through `survey::svydesign()` and
  `survey::svyby()`.
- The plug-in variance is approximate and does not fully account for
  uncertainty in estimated multilevel model parameters.
- Missing values are not removed automatically. The function stops with
  an informative error if required variables contain missing values.

## Summary function

- The argument `summary_function` supports `"mean"` and `"total"`
  because both are linear domain parameters.
- For `"mean"`, the synthetic component and residual correction are
  aggregated using `survey::svymean`.
- For `"total"`, both components are aggregated using
  `survey::svytotal`, so the estimate and variance are returned on the
  total scale.
- The `"total"` option should only be used when survey weights are
  appropriate expansion weights for population totals.

## References

Bates, D., Maechler, M., Bolker, B., & Walker, S. (2015). Fitting linear
mixed-effects models using lme4. *Journal of Statistical Software,
67*(1), 1–48.

Finch, W. H., Bolin, J. E., & Kelley, K. (2014). *Multilevel Modeling
Using R*. CRC Press.

Food and Agriculture Organization of the United Nations. (2021).
*Guidelines on Data Disaggregation for SDG Indicators Using Survey Data*
(1st ed.). <https://doi.org/10.4060/cb3253en>

Hox, J. J., Moerbeek, M., & van de Schoot, R. (2018). *Multilevel
Analysis: Techniques and Applications* (3rd ed.). Routledge.

Kim, J. K., & Rao, J. N. K. (2012). Combining data from two independent
surveys: A model-assisted approach. *Biometrika, 99*(1), 85–100.

Moura, F. A. S., & Holt, D. (1999). Small area estimation using
multilevel models. *Survey Methodology, 25*(1), 73–80.

Rao, J. N. K., & Molina, I. (2015). *Small Area Estimation* (2nd ed.).
Wiley.
