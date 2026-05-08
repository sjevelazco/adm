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

- hold_out_evaluation:

  logical. If `TRUE`, performance metrics will also be calculated for
  the `hold_out_set`.

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
if (FALSE) { # \dontrun{
require(dplyr)
# Database with species abundance and x and y coordinates
data("sppabund")
# Select data for a single species
some_sp <- sppabund %>%
  dplyr::filter(species == "Species two") %>%
  dplyr::select(-.part2, -.part3)
# Explore response variables
some_sp$ind_ha %>% range()
some_sp$ind_ha %>% hist()
# Here we balance number of absences
some_sp <-
  balance_dataset(some_sp, response = "ind_ha", absence_ratio = 0.2)
# Create a grid
xgb_grid <- expand.grid(
  nrounds = c(100, 300),
  max_depth = c(4, 6, 8),
  learning_rate = c(0.2, 0.5),
  min_split_loss = c(1, 5, 10),
  colsample_bytree = c(0.5, 1),
  min_child_weight = c(0.5, 1, 2),
  subsample = c(0.5, 1),
  stringsAsFactors = FALSE
)
# Tune a XGB model
tuned_xgb <- tune_abund_xgb(
  data = some_sp,
  response = "ind_ha",
  predictors = c("bio12", "elevation", "sand"),
  predictors_f = c("eco"),
  partition = ".part",
  predict_part = TRUE,
  metrics = c("corr_pear", "mae"),
  grid = xgb_grid,
  objective = "reg:squarederror",
  n_cores = 3
)
tuned_xgb
} # }
```
