# Fit and validate Support Vector Machine models with exploration of hyper-parameters that optimize performance

Fit and validate Support Vector Machine models with exploration of
hyper-parameters that optimize performance

## Usage

``` r
tune_abund_svm(
  data,
  response,
  predictors,
  predictors_f = NULL,
  fit_formula = NULL,
  partition,
  predict_part = FALSE,
  grid = NULL,
  metrics = NULL,
  n_cores = 1,
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

- predict_part:

  logical. Save predicted abundance for testing data. Default is FALSE.

- grid:

  tibble or data.frame. A dataframe with "kernel", "sigma", "C" as
  columns and its values combinations as rows. If now grid is provided,
  funcion will create a default grid combining the next hyperparameters:
  C = seq(0.2, 1, by = 0.2), sigma = "automatic", kernel = c("rbfdot",
  "laplacedot"). In case one or more hyperparameters are provided, the
  function will complete the grid with the default values.

- metrics:

  character. Vector with one or more metrics from
  c("corr_spear","corr_pear","mae","pdisp","inter","slope").

- n_cores:

  numeric. Number of cores used in parallel processing.

- verbose:

  logical. If FALSE, disables all console messages. Default TRUE

## Value

A list object with:

- model: A "ksvm" class object from kernlab package. This object can be
  used for predicting.

- predictors: A tibble with quantitative (c column names) and
  qualitative (f column names) variables use for modeling.

- performance: A tibble with selected model's performance metrics
  calculated in adm_eval.

- performance_part: A tibble with performance metrics for each test
  partition.

- predicted_part: A tibble with predicted abundance for each test
  partition.

- optimal_combination: A tibble with the selected hyperparameter
  combination and its performance.

- all_combinations: A tibble with all hyperparameters combinations and
  its performance.

## Examples

``` r
# \donttest{
require(dplyr)

# Database with species abundance and x and y coordinates
data("sppabund")

# Select data for a single species
some_sp <- sppabund %>%
  dplyr::filter(species == "Species one") %>%
  dplyr::select(-.part2, -.part3)

# Explore response variables
some_sp$ind_ha %>% range()
#> [1]   0.0000 285.2941
some_sp$ind_ha %>% hist()


# Here we balance number of absences
some_sp <-
  balance_dataset(some_sp, response = "ind_ha", absence_ratio = 0.2)

# A small subset keeps this example fast
set.seed(1)
some_sp <- some_sp[sort(sample(nrow(some_sp), 200)), ]

# Create a grid
svm_grid <- expand.grid(
  sigma = "automatic",
  C = c(0.5, 1),
  kernel = c("rbfdot", "laplacedot"),
  stringsAsFactors = FALSE
)

# Tune a SVM model
tuned_svm <- tune_abund_svm(
  data = some_sp,
  response = "ind_ha",
  predictors = c("bio12", "elevation", "sand"),
  predictors_f = c("eco"),
  partition = ".part1",
  predict_part = TRUE,
  metrics = c("corr_pear", "mae"),
  grid = svm_grid,
  n_cores = 1
)
#> Using provided grid.
#> Searching for optimal hyperparameters...
#> 
#> Fitting the best model...
#> Formula used for model fitting:
#> ind_ha ~ bio12 + elevation + sand + eco
#> Replica number: 1/1
#> -- Partition number 1/3
#> -- Partition number 2/3
#> -- Partition number 3/3
#> The best model was achieved with: 
#>  sigma = automatic, kernel = rbfdot and C = 1

tuned_svm
#> $model
#> Support Vector Machine object of class "ksvm" 
#> 
#> SV type: eps-svr  (regression) 
#>  parameter : epsilon = 0.1  cost C = 1 
#> 
#> Gaussian Radial Basis kernel function. 
#>  Hyperparameter : sigma =  0.966088059640013 
#> 
#> Number of Support Vectors : 142 
#> 
#> Objective Function Value : -74.7645 
#> Training error : 0.722657 
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
#> 1 svm       21.5   2.15           0.524         0.109          0.379
#> # ℹ 7 more variables: corr_pear_sd <dbl>, inter_mean <dbl>, inter_sd <dbl>,
#> #   slope_mean <dbl>, slope_sd <dbl>, pdisp_mean <dbl>, pdisp_sd <dbl>
#> 
#> $performance_part
#> # A tibble: 3 × 9
#>   replica partition model   mae corr_spear corr_pear inter slope pdisp
#>   <chr>   <chr>     <chr> <dbl>      <dbl>     <dbl> <dbl> <dbl> <dbl>
#> 1 1       1         svm    19.0      0.510     0.408  7.71 0.877 0.465
#> 2 1       2         svm    22.8      0.422     0.258 10.7  0.639 0.404
#> 3 1       3         svm    22.7      0.640     0.472 -4.14 2.05  0.230
#> 
#> $predicted_part
#> # A tibble: 200 × 4
#>    replica partition observed predicted
#>    <chr>   <chr>        <dbl>     <dbl>
#>  1 1       1             10       11.3 
#>  2 1       1             10        9.16
#>  3 1       1             10       10.8 
#>  4 1       1             15       22.5 
#>  5 1       1             13.3     24.9 
#>  6 1       1             40        6.65
#>  7 1       1             10       19.2 
#>  8 1       1             80       30.3 
#>  9 1       1             30       28.5 
#> 10 1       1             30       29.1 
#> # ℹ 190 more rows
#> 
#> $metadata
#> $metadata$plataform
#> [1] "aarch64-apple-darwin23"
#> 
#> $metadata$R_version
#> [1] "R version 4.6.1 (2026-06-24)"
#> 
#> $metadata$run_date
#> [1] "2026-10-21 13:40:38 PDT"
#> 
#> $metadata$formula
#> ind_ha ~ bio12 + elevation + sand + eco
#> <environment: 0x6dfc72bd8>
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
#> [1] "tune_abund_svm"
#> 
#> 
#> $optimal_combination
#>   comb_id     sigma C kernel model mae_mean  mae_sd corr_spear_mean
#> 2  comb_2 automatic 1 rbfdot   svm 20.79083 1.67909        0.570648
#>   corr_spear_sd corr_pear_mean corr_pear_sd inter_mean inter_sd slope_mean
#> 2    0.08967876      0.4301277   0.05942366     1.3867  7.69838   1.433585
#>    slope_sd pdisp_mean  pdisp_sd
#> 2 0.7235752  0.3380812 0.1197406
#> 
#> $all_combinations
#>   comb_id     sigma   C     kernel model mae_mean   mae_sd corr_spear_mean
#> 1  comb_1 automatic 0.5     rbfdot   svm 20.76384 1.772161       0.5790594
#> 2  comb_2 automatic 1.0     rbfdot   svm 20.79083 1.679090       0.5706480
#> 3  comb_3 automatic 0.5 laplacedot   svm 20.92380 2.036223       0.5696354
#> 4  comb_4 automatic 1.0 laplacedot   svm 20.75776 1.691549       0.5577426
#>   corr_spear_sd corr_pear_mean corr_pear_sd inter_mean inter_sd slope_mean
#> 1    0.10900659      0.4275871   0.05838918 -0.1943213 7.881332   1.516354
#> 2    0.08967876      0.4301277   0.05942366  1.3866995 7.698380   1.433585
#> 3    0.11434882      0.4099282   0.09958389 -2.6842485 9.821532   1.638330
#> 4    0.04875718      0.4251820   0.05863243 -1.3603613 7.211581   1.458333
#>    slope_sd pdisp_mean   pdisp_sd
#> 1 0.7558400  0.3167474 0.11366592
#> 2 0.7235752  0.3380812 0.11974058
#> 3 0.8674095  0.2765669 0.09023825
#> 4 0.6163046  0.3177822 0.10653144
#> 
# }
```
