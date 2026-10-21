# Fit and validate Generalized Boosted Regression models

Fit and validate Generalized Boosted Regression models

## Usage

``` r
fit_abund_gbm(
  data,
  response,
  predictors,
  predictors_f = NULL,
  fit_formula = NULL,
  partition,
  hold_out_set = NULL,
  predict_part = FALSE,
  distribution,
  n.trees = 100,
  interaction.depth = 5,
  n.minobsinnode = 5,
  shrinkage = 0.1,
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

- hold_out_set:

  tibble or data.frame. Optional hold-out dataset, never used during
  training, with the response and predictor columns. If provided,
  performance metrics are also calculated on it. Default NULL

- predict_part:

  logical. Save predicted abundance for testing data. Default = FALSE

- distribution:

  character. A string specifying the distribution to be used. See
  gbm::gbm documentation for details.

- n.trees:

  integer. The total number of trees to fit.

- interaction.depth:

  integer. The maximum depth of each tree. Default 5

- n.minobsinnode:

  integer. The minimum number of observations in the terminal nodes of
  the trees. Default 5

- shrinkage:

  numeric. The learning rate of the algorithm. Default 0.1

- verbose:

  logical. If FALSE, disables all console messages. Default TRUE

## Value

A list object with:

- model: A "gbm" class object from gbm package. This object can be used
  for predicting.

- predictors: A tibble with quantitative (c column names) and
  qualitative (f column names) variables use for modeling.

- performance: Averaged performance metrics (see
  [`adm_eval`](https://sjevelazco.github.io/adm/reference/adm_eval.md)).

- performance_part: Performance metrics for each replica and partition.

- predicted_part: Observed and predicted abundance for each test
  partition.

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

# Fit a GBM model
mgbm <- fit_abund_gbm(
  data = some_sp,
  response = "ind_ha",
  predictors = c("bio12", "elevation", "sand"),
  predictors_f = c("eco"),
  partition = ".part1",
  distribution = "gaussian",
  n.trees = 100,
  interaction.depth = 5,
  n.minobsinnode = 5,
  shrinkage = 0.1,
  predict_part = TRUE
)
#> Formula used for model fitting:
#> ind_ha ~ bio12 + elevation + sand + eco
#> Replica number: 1/1
#> -- Partition number 1/3
#> Using 100 trees...
#> -- Partition number 2/3
#> Using 100 trees...
#> -- Partition number 3/3
#> Using 100 trees...

mgbm
#> $model
#> gbm::gbm(formula = formula1, distribution = distribution, data = data, 
#>     n.trees = n.trees, interaction.depth = interaction.depth, 
#>     n.minobsinnode = n.minobsinnode, shrinkage = shrinkage, bag.fraction = 0.9)
#> A gradient boosted model with gaussian loss function.
#> 100 iterations were performed.
#> There were 4 predictors of which 4 had non-zero influence.
#> 
#> $predictors
#> # A tibble: 1 × 6
#>   model response c1    c2        c3    f    
#>   <chr> <chr>    <chr> <chr>     <chr> <chr>
#> 1 gbm   ind_ha   bio12 elevation sand  eco  
#> 
#> $performance
#> # A tibble: 1 × 13
#>   model mae_mean mae_sd corr_spear_mean corr_spear_sd corr_pear_mean
#>   <chr>    <dbl>  <dbl>           <dbl>         <dbl>          <dbl>
#> 1 gbm       22.9   1.77           0.419        0.0965          0.309
#> # ℹ 7 more variables: corr_pear_sd <dbl>, inter_mean <dbl>, inter_sd <dbl>,
#> #   slope_mean <dbl>, slope_sd <dbl>, pdisp_mean <dbl>, pdisp_sd <dbl>
#> 
#> $performance_part
#> # A tibble: 3 × 9
#>   replica partition model   mae corr_spear corr_pear inter slope pdisp
#>   <chr>   <chr>     <chr> <dbl>      <dbl>     <dbl> <dbl> <dbl> <dbl>
#> 1 1       1         gbm    21.5      0.451     0.395  9.72 0.648 0.610
#> 2 1       2         gbm    22.3      0.311     0.189 17.1  0.343 0.550
#> 3 1       3         gbm    24.9      0.496     0.343 10.2  0.726 0.473
#> 
#> $predicted_part
#> # A tibble: 562 × 4
#>    replica partition observed predicted
#>    <chr>   <chr>        <dbl>     <dbl>
#>  1 1       1             13.3      13.4
#>  2 1       1             20        10.7
#>  3 1       1             10        14.0
#>  4 1       1             23.3      11.5
#>  5 1       1             10        28.5
#>  6 1       1             10        28.5
#>  7 1       1             20        22.3
#>  8 1       1             10        11.8
#>  9 1       1             15        12.7
#> 10 1       1             10        11.8
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
#> [1] "2026-10-21 13:39:14 PDT"
#> 
#> $metadata$formula
#> ind_ha ~ bio12 + elevation + sand + eco
#> <environment: 0x6d1d02270>
#> 
#> $metadata$hyperparameters
#> $metadata$hyperparameters$distribution
#> [1] "gaussian"
#> 
#> $metadata$hyperparameters$n.trees
#> [1] 100
#> 
#> $metadata$hyperparameters$interaction.depth
#> [1] 5
#> 
#> $metadata$hyperparameters$n.minobsinnode
#> [1] 5
#> 
#> $metadata$hyperparameters$shrinkage
#> [1] 0.1
#> 
#> $metadata$hyperparameters$bag.fraction
#> [1] 0.9
#> 
#> 
#> $metadata$source_function
#> [1] "fit_abund_gbm"
#> 
#> 
# }
```
