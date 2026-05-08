# Fit and validate Extreme Gradient Boosting models

Fit and validate Extreme Gradient Boosting models

## Usage

``` r
fit_abund_xgb(
  data,
  response,
  predictors,
  predictors_f = NULL,
  partition,
  hold_out_set = NULL,
  predict_part = FALSE,
  nrounds = 1000,
  max_depth = 5,
  learning_rate = 0.1,
  min_split_loss = 1,
  colsample_bytree = 1,
  min_child_weight = 1,
  subsample = 0.5,
  objective = "reg:squarederror",
  early_stopping = list(cv_strategy = 10, fm_strategy = "median"),
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

- hold_out_set:

  tibble or data.frame. A hold-out dataset used for evaluation and early
  stopping. This data is never used during the training phase.

- predict_part:

  logical. Save predicted abundance for testing data. Default = FALSE.

- nrounds:

  integer. Max number of boosting iterations. Default is 100.

- max_depth:

  integer. The maximum depth of each tree. Default 5

- learning_rate:

  numeric. The learning rate of the algorithm. Default 0.1

- min_split_loss:

  numeric. Minimum loss reduction required to make a further partition
  on a leaf node of the tree. Default is 1.

- colsample_bytree:

  numeric. Subsample ratio of columns when constructing each tree.
  Default is 1.

- min_child_weight:

  numeric. Minimum sum of instance weight needed in a child. Default is
  1.

- subsample:

  numeric. Subsample ratio of the training instance. Default is 0.5.

- objective:

  character. The learning task and the corresponding learning objective.
  Default is "reg:squarederror", regression with squared loss.

- early_stopping:

  A list containing two elements:

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

- verbose:

  logical. If FALSE, disables all console messages. Default TRUE.

- hold_out_evaluation:

  logical. If `TRUE`, performance metrics will also be calculated for
  the `hold_out_set`.

## Value

A list object with:

- model: A "xgb.Booster" class object from xgboost package. This object
  can be used for predicting.

- predictors: A tibble with quantitative (c column names) and
  qualitative (f column names) variables use for modeling.

- performance: Averaged performance metrics (see
  [`adm_eval`](https://sjevelazco.github.io/adm/reference/adm_eval.md)).

- performance_part: Performance metrics for each replica and partition.

- predicted_part: Observed and predicted abundance for each test
  partition.

## Examples

``` r
if (FALSE) { # \dontrun{
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
some_sp$ind_ha %>% hist()

# Here we balance number of absences
some_sp <-
  balance_dataset(some_sp, response = "ind_ha", absence_ratio = 0.2)

# Fit a XGB model
mxgb <- fit_abund_xgb(
  data = some_sp,
  response = "ind_ha",
  predictors = c("bio12", "elevation", "sand"),
  predictors_f = NULL,
  partition = ".part",
  nrounds = 200,
  max_depth = 5,
  learning_rate = 0.1,
  min_split_loss = 1,
  colsample_bytree = 0.7,
  min_child_weight = 2,
  subsample = 0.3,
  objective = "reg:squarederror",
  predict_part = TRUE
)

mxgb
} # }
```
