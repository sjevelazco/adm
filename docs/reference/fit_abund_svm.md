# Fit and validate Support Vector Machine models

Fit and validate Support Vector Machine models

## Usage

``` r
fit_abund_svm(
  data,
  response,
  predictors,
  predictors_f = NULL,
  fit_formula = NULL,
  partition,
  hold_out_set = NULL,
  predict_part = FALSE,
  kernel = "rbfdot",
  sigma = "automatic",
  C = 1,
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

  logical. Save predicted abundance for testing data. Default is FALSE.

- kernel:

  character. A string defining the kernel used in the algorithm. Default
  is "rbfdot".

- sigma:

  numeric or character. Either "automatic" (recommended) or the inverse
  kernel width for the Radial Basis kernel function "rbfdot" and the
  Laplacian kernel "laplacedot". Default is "automatic".

- C:

  numeric. Cost of constraints violation. Default is 1.

- verbose:

  logical. If FALSE, disables all console messages. Default TRUE

## Value

A list object with:

- model: A "ksvm" class object from kernlab package. This object can be
  used for predicting.

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

# Fit a SVM model
msvm <- fit_abund_svm(
  data = some_sp,
  response = "ind_ha",
  predictors = c("bio12", "elevation", "sand"),
  predictors_f = c("eco"),
  partition = ".part1",
  kernel = "rbfdot",
  sigma = "automatic",
  C = 1,
  predict_part = TRUE
)
#> Formula used for model fitting:
#> ind_ha ~ bio12 + elevation + sand + eco
#> Replica number: 1/1
#> -- Partition number 1/3
#> -- Partition number 2/3
#> -- Partition number 3/3

msvm
#> $model
#> Support Vector Machine object of class "ksvm" 
#> 
#> SV type: eps-svr  (regression) 
#>  parameter : epsilon = 0.1  cost C = 1 
#> 
#> Gaussian Radial Basis kernel function. 
#>  Hyperparameter : sigma =  1.0095575217254 
#> 
#> Number of Support Vectors : 426 
#> 
#> Objective Function Value : -229.9293 
#> Training error : 0.818788 
#> 
#> $predictors
#> # A tibble: 1 × 6
#>   model response c1    c2        c3    f    
#>   <chr> <chr>    <chr> <chr>     <chr> <chr>
#> 1 svm   ind_ha   bio12 elevation sand  eco  
#> 
#> $performance
#> # A tibble: 1 × 13
#>   model mae_mean mae_sd corr_spear_mean corr_spear_sd corr_pear_mean
#>   <chr>    <dbl>  <dbl>           <dbl>         <dbl>          <dbl>
#> 1 svm       20.2   2.00           0.475        0.0591          0.377
#> # ℹ 7 more variables: corr_pear_sd <dbl>, inter_mean <dbl>, inter_sd <dbl>,
#> #   slope_mean <dbl>, slope_sd <dbl>, pdisp_mean <dbl>, pdisp_sd <dbl>
#> 
#> $performance_part
#> # A tibble: 3 × 9
#>   replica partition model   mae corr_spear corr_pear inter slope pdisp
#>   <chr>   <chr>     <chr> <dbl>      <dbl>     <dbl> <dbl> <dbl> <dbl>
#> 1 1       1         svm    19.9      0.502     0.400  8.61 0.934 0.428
#> 2 1       2         svm    18.5      0.407     0.274 10.6  0.764 0.358
#> 3 1       3         svm    22.4      0.515     0.458 -2.05 1.58  0.290
#> 
#> $predicted_part
#> # A tibble: 562 × 4
#>    replica partition observed predicted
#>    <chr>   <chr>        <dbl>     <dbl>
#>  1 1       1             13.3     16.2 
#>  2 1       1             20       13.8 
#>  3 1       1             10       14.8 
#>  4 1       1             23.3     14.1 
#>  5 1       1             10       23.3 
#>  6 1       1             10       16.5 
#>  7 1       1             20       15.4 
#>  8 1       1             10        7.50
#>  9 1       1             15       10.9 
#> 10 1       1             10       11.4 
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
#> <environment: 0x6dcf66ee8>
#> 
#> $metadata$type
#> [1] "eps-svr"
#> 
#> $metadata$hyperparameters
#> $metadata$hyperparameters$kernel
#> [1] "rbfdot"
#> 
#> $metadata$hyperparameters$sigma
#> [1] "automatic"
#> 
#> $metadata$hyperparameters$C
#> [1] 1
#> 
#> 
#> $metadata$source_function
#> [1] "fit_abund_svm"
#> 
#> 
# }
```
