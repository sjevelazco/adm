# Fit and validate Random Forest models with exploration of hyper-parameters that optimize performance

Fit and validate Random Forest models with exploration of
hyper-parameters that optimize performance

## Usage

``` r
tune_abund_raf(
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

  logical. Save predicted abundance for testing data. Default = FALSE

- grid:

  tibble or data.frame. A dataframe with "mtry" and "ntree" as columns
  and its values combinations as rows. If no grid is provided, function
  will create a default grid combining the next hyperparameters: mtry =
  seq(2, length(predictors), by = 1), ntree = seq(500, 1000, by = 100).
  In case one or more hyperparameters are provided, the function will
  complete the grid with the default values.

- metrics:

  character. Vector with one or more metrics from
  c("corr_spear","corr_pear","mae","pdisp","inter","slope").

- n_cores:

  numeric. Number of cores used in parallel processing.

- verbose:

  logical. If FALSE, disables all console messages. Default TRUE

## Value

A list object with:

- model: A "randomForest" class object from randomForest package. This
  object can be used for predicting.

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
  dplyr::filter(species == "Species two") %>%
  dplyr::select(-.part2, -.part3)

# Explore response variables
some_sp$ind_ha %>% range()
#> [1]   0.0000 186.8627
some_sp$ind_ha %>% hist()


# Here we balance number of absences
some_sp <-
  balance_dataset(some_sp, response = "ind_ha", absence_ratio = 0.2)

# A small subset keeps this example fast
set.seed(1)
some_sp <- some_sp[sort(sample(nrow(some_sp), 200)), ]

# Create a grid
raf_grid <- expand.grid(
  mtry = 2:3,
  ntree = c(50, 100),
  stringsAsFactors = FALSE
)

# Tune a RAF model
tuned_raf <- tune_abund_raf(
  data = some_sp,
  response = "ind_ha",
  predictors = c("bio12", "elevation", "sand"),
  predictors_f = c("eco"),
  partition = ".part1",
  predict_part = TRUE,
  metrics = c("corr_pear", "mae"),
  grid = raf_grid,
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
#>  mtry = 2 and ntree = 100

tuned_raf
#> $model
#> 
#> Call:
#>  randomForest(formula = formula1, data = data, mtry = mtry, ntree = ntree,      importance = FALSE) 
#>                Type of random forest: regression
#>                      Number of trees: 100
#> No. of variables tried at each split: 2
#> 
#>           Mean of squared residuals: 728.9428
#>                     % Var explained: 4.35
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
#> 1 raf       20.3  0.558           0.270        0.0770          0.215
#> # ℹ 7 more variables: corr_pear_sd <dbl>, inter_mean <dbl>, inter_sd <dbl>,
#> #   slope_mean <dbl>, slope_sd <dbl>, pdisp_mean <dbl>, pdisp_sd <dbl>
#> 
#> $performance_part
#> # A tibble: 3 × 9
#>   replica partition model   mae corr_spear corr_pear inter slope pdisp
#>   <chr>   <chr>     <chr> <dbl>      <dbl>     <dbl> <dbl> <dbl> <dbl>
#> 1 1       1         raf    19.8      0.329    0.279   8.36 0.399 0.700
#> 2 1       2         raf    20.9      0.183    0.0789 20.7  0.215 0.367
#> 3 1       3         raf    20.4      0.297    0.287   8.12 0.832 0.345
#> 
#> $predicted_part
#> # A tibble: 200 × 4
#>    replica partition observed predicted
#>    <chr>   <chr>        <dbl>     <dbl>
#>  1 1       1               50      35.9
#>  2 1       1               30      33.5
#>  3 1       1               10      12.7
#>  4 1       1               10      19.5
#>  5 1       1               10      35.5
#>  6 1       1               10      30.1
#>  7 1       1               40      32.9
#>  8 1       1               10      26.8
#>  9 1       1               50      48.7
#> 10 1       1               20      32.1
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
#> [1] "2026-10-21 13:40:35 PDT"
#> 
#> $metadata$formula
#> ind_ha ~ bio12 + elevation + sand + eco
#> <environment: 0x6dda9e0e8>
#> 
#> $metadata$importance
#> [1] FALSE
#> 
#> $metadata$hyperparameters
#> $metadata$hyperparameters$mtry
#> [1] 2
#> 
#> $metadata$hyperparameters$ntree
#> [1] 100
#> 
#> 
#> $metadata$source_function
#> [1] "tune_abund_raf"
#> 
#> 
#> $optimal_combination
#>   comb_id mtry ntree model mae_mean   mae_sd corr_spear_mean corr_spear_sd
#> 3  comb_3    2   100   raf 20.26949 0.858931       0.2715684     0.1008452
#>   corr_pear_mean corr_pear_sd inter_mean inter_sd slope_mean  slope_sd
#> 3      0.2271266    0.1423189    11.6307 8.483081  0.5058793 0.3820859
#>   pdisp_mean  pdisp_sd
#> 3  0.4704724 0.1968739
#> 
#> $all_combinations
#>   comb_id mtry ntree model mae_mean    mae_sd corr_spear_mean corr_spear_sd
#> 1  comb_1    2    50   raf 20.42383 0.8788225       0.2596367    0.11386300
#> 2  comb_2    3    50   raf 20.99984 0.7110636       0.2425709    0.07487930
#> 3  comb_3    2   100   raf 20.26949 0.8589310       0.2715684    0.10084522
#> 4  comb_4    3   100   raf 20.63031 0.3296048       0.2442177    0.07739837
#>   corr_pear_mean corr_pear_sd inter_mean  inter_sd slope_mean  slope_sd
#> 1      0.2196949    0.1813738   11.43846 11.869116  0.5390408 0.5849263
#> 2      0.1907682    0.1409322   14.07839  8.515779  0.4057176 0.4035233
#> 3      0.2271266    0.1423189   11.63070  8.483081  0.5058793 0.3820859
#> 4      0.1978030    0.1351581   13.56967  8.135619  0.4319145 0.3789949
#>   pdisp_mean  pdisp_sd
#> 1  0.4692900 0.1973203
#> 2  0.5365805 0.2830610
#> 3  0.4704724 0.1968739
#> 4  0.4975310 0.2123992
#> 
# }
```
