# Fit and validate Extreme Gradient Boosting models with exploration of hyper-parameters that optimize performance

Fit and validate Extreme Gradient Boosting models with exploration of
hyper-parameters that optimize performance

## Usage

``` r
tune_abund_xgb(
  data,
  response,
  predictors,
  predictors_f = NULL,
  partition,
  predict_part = FALSE,
  hold_out_set = NULL,
  grid = NULL,
  objective = "reg:squarederror",
  metrics = NULL,
  early_stopping = list(cv_strategy = 10, fm_strategy = "median"),
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

- partition:

  character. Column name with training and validation partition groups.

- predict_part:

  logical. Save predicted abundance for testing data. Default = FALSE

- hold_out_set:

  tibble or data.frame. A hold-out dataset used for evaluation and early
  stopping. This data is never used during the training phase.

- grid:

  tibble or data.frame. A dataframe with "n.trees", "interaction.depth",
  "n.minobsinnode" and "shrinkage" as columns and its values
  combinations as rows. If no grid is provided, function will create a
  default grid combining the next hyperparameters: nrounds = c(100, 200,
  300), max_depth = c(4, 6, 8), learning_rate = c(0.2, 0.4, 0.5),
  min_split_loss = c(1, 5, 10), colsample_bytree = c(0.5, 1, 2),
  min_child_weight = c(0.5, 1, 2), subsample = c(0.5, 0.75, 1). In case
  one or more hyperparameters are provided, the function will complete
  the grid with the default values.

- objective:

  character. The learning task and the corresponding learning objective.
  Default is "reg:squarederror", regression with squared loss.

- metrics:

  character. Vector with one or more metrics from
  c("corr_spear","corr_pear","mae","pdisp","inter","slope").

- early_stopping:

  NULL or a list containing two elements:

  - `cv_strategy`: Numerical. Specifies the number of rounds without
    improvement before training stops during the cross-validation stage.

  - `fm_strategy`: A vector defining the strategy for the final model:

    - `c("hold_out", n)`: Stops training after `n` rounds without
      improvement, using the `hold_out_set` as the evaluation set.

    - `c("mean")`: Trains the final model using the average number of
      rounds reached across all cross-validation folds.

    - `c("median")`: Uses the median number of rounds from the
      cross-validation stage.

    - `c("max")`: Uses the maximum number of rounds reached in any fold.

    - `c("min")`: Uses the minimum number of rounds reached in any fold.

- n_cores:

  numeric. Number of cores used in parallel processing.

- verbose:

  logical. If FALSE, disables all console messages. Default TRUE

## Value

A list object with:

- model: A "xgb.Booster" object from xgboost package. This object can be
  used to predicting.

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
xgb_grid <- expand.grid(
  nrounds = c(20, 50),
  max_depth = c(4, 6),
  learning_rate = 0.3,
  min_split_loss = 1,
  colsample_bytree = 1,
  min_child_weight = 1,
  subsample = 0.5,
  stringsAsFactors = FALSE
)
# Tune a XGB model
tuned_xgb <- tune_abund_xgb(
  data = some_sp,
  response = "ind_ha",
  predictors = c("bio12", "elevation", "sand"),
  predictors_f = c("eco"),
  partition = ".part1",
  predict_part = TRUE,
  metrics = c("corr_pear", "mae"),
  grid = xgb_grid,
  objective = "reg:squarederror",
  n_cores = 1
)
#> Using provided grid.
#> Searching for optimal hyperparameters...
#> 
#> Fitting the best model...
#> Warning: Categorical variables aren't available for XGB and will be ignored.
#> Replica number: 1/1
#> -- Partition number 1/3
#> -- Partition number 2/3
#> -- Partition number 3/3
#> The best model was achieved with: 
#>  max_depth = 4
#>  learning_rate = 0.3
#>  min_split_loss = 1
#>  colsample_bytree = 1
#>  min_child_weight = 1
#>  subsample = 0.5
#>  nrounds = 20
tuned_xgb
#> $model
#> ##### xgb.Booster
#> call:
#>   xgboost::xgb.train(params = list(objective = objective, max_depth = max_depth, 
#>     learning_rate = learning_rate, min_child_weight = min_child_weight, 
#>     subsample = subsample, colsample_bytree = colsample_bytree, 
#>     min_split_loss = min_split_loss, seed = 13), data = dfull, 
#>     nrounds = ifelse(early_stopping$fm_strategy[[1]] == "hold_out", 
#>         yes = nrounds, no = early_stop_interpreter(early_stopping, 
#>             train_observer, nrounds)), evals = fm_evals, verbose = 0, 
#>     early_stopping_rounds = fm_early_stopping)
#> # of features: 3 
#> # of rounds:  2 
#> 
#> $predictors
#> # A tibble: 1 × 5
#>   model response c1    c2        c3   
#>   <chr> <chr>    <chr> <chr>     <chr>
#> 1 xgb   ind_ha   bio12 elevation sand 
#> 
#> $performance
#> # A tibble: 1 × 13
#>   model mae_mean mae_sd corr_spear_mean corr_spear_sd corr_pear_mean
#>   <chr>    <dbl>  <dbl>           <dbl>         <dbl>          <dbl>
#> 1 xgb       18.8   1.81           0.213         0.248          0.138
#> # ℹ 7 more variables: corr_pear_sd <dbl>, inter_mean <dbl>, inter_sd <dbl>,
#> #   slope_mean <dbl>, slope_sd <dbl>, pdisp_mean <dbl>, pdisp_sd <dbl>
#> 
#> $performance_part
#> # A tibble: 3 × 9
#>   replica partition model   mae corr_spear corr_pear inter  slope pdisp
#>   <chr>   <chr>     <chr> <dbl>      <dbl>     <dbl> <dbl>  <dbl> <dbl>
#> 1 1       1         xgb    16.9     0.287     0.150  10.1   0.437 0.343
#> 2 1       2         xgb    20.6    -0.0633   -0.0514 36.9  -0.373 0.138
#> 3 1       3         xgb    19.0     0.416     0.315  -8.81  1.45  0.217
#> 
#> $predicted_part
#> # A tibble: 200 × 4
#>    replica partition observed predicted
#>    <chr>   <chr>        <dbl>     <dbl>
#>  1 1       1               50      23.6
#>  2 1       1               30      23.6
#>  3 1       1               10      19.2
#>  4 1       1               10      19.2
#>  5 1       1               10      23.2
#>  6 1       1               10      25.1
#>  7 1       1               40      25.1
#>  8 1       1               10      25.1
#>  9 1       1               50      25.1
#> 10 1       1               20      25.8
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
#> [1] "2026-10-21 13:40:41 PDT"
#> 
#> $metadata$evals
#> NULL
#> 
#> $metadata$early_stopping_rounds
#> NULL
#> 
#> $metadata$boosted_rounds
#> [1] 2
#> 
#> $metadata$hyperparameters
#> $metadata$hyperparameters$objective
#> [1] "reg:squarederror"
#> 
#> $metadata$hyperparameters$max_depth
#> [1] 4
#> 
#> $metadata$hyperparameters$learning_rate
#> [1] 0.3
#> 
#> $metadata$hyperparameters$min_child_weight
#> [1] 1
#> 
#> $metadata$hyperparameters$subsample
#> [1] 0.5
#> 
#> $metadata$hyperparameters$colsample_bytree
#> [1] 1
#> 
#> $metadata$hyperparameters$min_split_loss
#> [1] 1
#> 
#> 
#> $metadata$source_function
#> [1] "tune_abund_xgb"
#> 
#> 
#> $optimal_combination
#>   comb_id nrounds max_depth learning_rate min_split_loss colsample_bytree
#> 1  comb_1      20         4           0.3              1                1
#>   min_child_weight subsample model mae_mean   mae_sd corr_spear_mean
#> 1                1       0.5   xgb 18.83134 1.809452       0.2131009
#>   corr_spear_sd corr_pear_mean corr_pear_sd inter_mean inter_sd slope_mean
#> 1     0.2478548      0.1379899     0.183625   12.72398 22.96598   0.505014
#>    slope_sd pdisp_mean  pdisp_sd
#> 1 0.9136537  0.2328386 0.1036295
#> 
#> $all_combinations
#>   comb_id nrounds max_depth learning_rate min_split_loss colsample_bytree
#> 1  comb_1      20         4           0.3              1                1
#> 2  comb_2      50         4           0.3              1                1
#> 3  comb_3      20         6           0.3              1                1
#> 4  comb_4      50         6           0.3              1                1
#>   min_child_weight subsample model mae_mean   mae_sd corr_spear_mean
#> 1                1       0.5   xgb 18.83134 1.809452       0.2131009
#> 2                1       0.5   xgb 18.83134 1.809452       0.2131009
#> 3                1       0.5   xgb 19.16952 1.937502       0.2365540
#> 4                1       0.5   xgb 19.16952 1.937502       0.2365540
#>   corr_spear_sd corr_pear_mean corr_pear_sd inter_mean inter_sd slope_mean
#> 1     0.2478548      0.1379899    0.1836250   12.72398 22.96598   0.505014
#> 2     0.2478548      0.1379899    0.1836250   12.72398 22.96598   0.505014
#> 3     0.1802596      0.1347968    0.1659993   14.55690 18.75889   0.426665
#> 4     0.1802596      0.1347968    0.1659993   14.55690 18.75889   0.426665
#>    slope_sd pdisp_mean  pdisp_sd
#> 1 0.9136537  0.2328386 0.1036295
#> 2 0.9136537  0.2328386 0.1036295
#> 3 0.7085001  0.2666522 0.1153665
#> 4 0.7085001  0.2666522 0.1153665
#> 
# }
```
