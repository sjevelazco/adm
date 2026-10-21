# Fit and validate Random Forests models

Fit and validate Random Forests models

## Usage

``` r
fit_abund_raf(
  data,
  response,
  predictors,
  predictors_f = NULL,
  fit_formula = NULL,
  partition = NULL,
  hold_out_set = NULL,
  predict_part = FALSE,
  mtry = length(c(predictors, predictors_f))/3,
  ntree = 500,
  verbose = TRUE
)
```

## Arguments

- data:

  tibble or data.frame. Database with response, predictors, and
  partition values

- response:

  character. Column name with species abundance.

- predictors:

  character. Vector with the column names of quantitative predictor
  variables (i.e. continuous variables). Usage predictors = c("temp",
  "precipt", "sand")

- predictors_f:

  character. Vector with the column names of qualitative predictor
  variables (i.e. ordinal or nominal variables type). Usage predictors_f
  = c("landform")

- fit_formula:

  formula. A formula object with response and predictor variables (e.g.
  formula(abund ~ temp + precipt + sand + landform)). Note that the
  variables used here must be consistent with those used in response,
  predictors, and predictors_f arguments. Default NULL

- partition:

  character. Column name with training and validation partition groups.
  If partition = NULL, the model will be validated with the same data
  used for fitting.

- hold_out_set:

  tibble or data.frame. Optional hold-out dataset, never used during
  training, with the response and predictor columns. If provided,
  performance metrics are also calculated on it. Default NULL

- predict_part:

  logical. Save predicted abundance for testing data. Default is FALSE.

- mtry:

  numeric. Number of variables randomly sampled as candidates at each
  split. Default (length(c(predictors, predictors_f))/3)

- ntree:

  numeric. Number of trees to grow. This should not be set to too small
  a number, to ensure that every input row gets predicted at least a few
  times. Default 500

- verbose:

  logical. If FALSE, disables all console messages. Default TRUE

## Value

A list object with:

- model: A "randomForest" class object from randomForest package. This
  object can be used for predicting.

- predictors: A tibble with quantitative (c column names) and
  qualitative (f column names) variables use for modeling.

- performance: Averaged performance metrics (see
  [`adm_eval`](https://sjevelazco.github.io/adm/reference/adm_eval.md)).

- performance_part: Performance metrics for each replica and partition.

- predicted_part: Observed and predicted abundance for each test
  partition.

## See also

[`fit_abund_cnn`](https://sjevelazco.github.io/adm/reference/fit_abund_cnn.md),
[`fit_abund_dnn`](https://sjevelazco.github.io/adm/reference/fit_abund_dnn.md),
[`fit_abund_gam`](https://sjevelazco.github.io/adm/reference/fit_abund_gam.md),
[`fit_abund_gbm`](https://sjevelazco.github.io/adm/reference/fit_abund_gbm.md),
[`fit_abund_glm`](https://sjevelazco.github.io/adm/reference/fit_abund_glm.md),
[`fit_abund_net`](https://sjevelazco.github.io/adm/reference/fit_abund_net.md),
[`fit_abund_svm`](https://sjevelazco.github.io/adm/reference/fit_abund_svm.md),
[`fit_abund_xgb`](https://sjevelazco.github.io/adm/reference/fit_abund_xgb.md)

## Examples

``` r
# \donttest{
require(terra)
require(dplyr)

# Database with species abundance and x and y coordinates
data("sppabund")

# Extract data for a single species
some_sp <- sppabund %>%
  dplyr::filter(species == "Species one") %>%
  dplyr::select(-.part2, -.part3)

# Explore reponse variables
some_sp$ind_ha %>% range()
#> [1]   0.0000 285.2941
some_sp$ind_ha %>% hist()


# Here we balance number of absences
some_sp <-
  balance_dataset(some_sp, response = "ind_ha", absence_ratio = 0.2)

# Fit a RAF model
mraf <- fit_abund_raf(
  data = some_sp,
  response = "ind_ha",
  predictors = c("bio12", "elevation", "sand"),
  predictors_f = c("eco"),
  partition = ".part1",
  mtry = 3,
  ntree = 500,
  predict_part = TRUE
)
#> Formula used for model fitting:
#> ind_ha ~ bio12 + elevation + sand + eco
#> Replica number: 1/1
#> -- Partition number 1/3
#> -- Partition number 2/3
#> -- Partition number 3/3

mraf
#> $model
#> 
#> Call:
#>  randomForest(formula = formula1, data = data, mtry = mtry, ntree = ntree,      importance = FALSE) 
#>                Type of random forest: regression
#>                      Number of trees: 500
#> No. of variables tried at each split: 3
#> 
#>           Mean of squared residuals: 1313.612
#>                     % Var explained: 4.09
#> 
#> $predictors
#> # A tibble: 1 × 6
#>   model response c1    c2        c3    f    
#>   <chr> <chr>    <chr> <chr>     <chr> <chr>
#> 1 raf   ind_ha   bio12 elevation sand  eco  
#> 
#> $performance
#> # A tibble: 1 × 13
#>   model mae_mean mae_sd corr_spear_mean corr_spear_sd corr_pear_mean
#>   <chr>    <dbl>  <dbl>           <dbl>         <dbl>          <dbl>
#> 1 raf       23.1   1.06           0.415        0.0614          0.310
#> # ℹ 7 more variables: corr_pear_sd <dbl>, inter_mean <dbl>, inter_sd <dbl>,
#> #   slope_mean <dbl>, slope_sd <dbl>, pdisp_mean <dbl>, pdisp_sd <dbl>
#> 
#> $performance_part
#> # A tibble: 3 × 9
#>   replica partition model   mae corr_spear corr_pear inter slope pdisp
#>   <chr>   <chr>     <chr> <dbl>      <dbl>     <dbl> <dbl> <dbl> <dbl>
#> 1 1       1         raf    22.6      0.395     0.356  9.88 0.570 0.624
#> 2 1       2         raf    22.4      0.366     0.173 18.7  0.294 0.591
#> 3 1       3         raf    24.3      0.484     0.401  5.66 0.879 0.456
#> 
#> $predicted_part
#> # A tibble: 562 × 4
#>    replica partition observed predicted
#>    <chr>   <chr>        <dbl>     <dbl>
#>  1 1       1             13.3      14.5
#>  2 1       1             20        13.8
#>  3 1       1             10        14.2
#>  4 1       1             23.3      12.9
#>  5 1       1             10        37.9
#>  6 1       1             10        41.5
#>  7 1       1             20        13.9
#>  8 1       1             10        10.3
#>  9 1       1             15        11.6
#> 10 1       1             10        12.1
#> # ℹ 552 more rows
#> 
#> $metadata
#> $metadata$plataform
#> [1] "aarch64-apple-darwin23"
#> 
#> $metadata$R_version
#> [1] "R version 4.6.1 (2026-06-24)"
#> 
#> $metadata$run_date
#> [1] "2026-10-21 13:39:17 PDT"
#> 
#> $metadata$formula
#> ind_ha ~ bio12 + elevation + sand + eco
#> <environment: 0x6de738938>
#> 
#> $metadata$importance
#> [1] FALSE
#> 
#> $metadata$hyperparameters
#> $metadata$hyperparameters$mtry
#> [1] 3
#> 
#> $metadata$hyperparameters$ntree
#> [1] 500
#> 
#> 
#> $metadata$source_function
#> [1] "fit_abund_raf"
#> 
#> 
# }
```
