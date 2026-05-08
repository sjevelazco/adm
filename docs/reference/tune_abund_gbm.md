# Fit and validate Generalized Boosted Regression models with exploration of hyper-parameters that optimize performance

Fit and validate Generalized Boosted Regression models with exploration
of hyper-parameters that optimize performance

## Usage

``` r
tune_abund_gbm(
  data,
  response,
  predictors,
  predictors_f = NULL,
  fit_formula = NULL,
  partition,
  predict_part = FALSE,
  grid = NULL,
  distribution,
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

  tibble or data.frame. A dataframe with "n.trees", "interaction.depth",
  "n.minobsinnode" and "shrinkage" as columns and its values as rows. If
  no grid is provided, function will create a default grid combining the
  next hyperparameters: n.trees = c(100, 200, 300), interaction.depth =
  c(1, 2, 3), n.minobsinnode = c(5, 10, 15), shrinkage = seq(0.001, 0.1,
  by = 0.05). In case one or more hyperparameters are provided, the
  function will complete the grid with the default values.

- distribution:

  character. A string specifying the distribution to be used. See
  gbm::gbm documentation for details.

- metrics:

  character. Vector with one or more metrics from
  c("corr_spear","corr_pear","mae","pdisp","inter","slope").

- n_cores:

  numeric. Number of cores used in parallel processing.

- verbose:

  logical. If FALSE, disables all console messages. Default TRUE

## Value

A list object with:

- model: A "gbm" object from gbm package. This object can be used to
  predicting.

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

- selected_arch: A numeric vector describing the selected architecture
  layers.

## Examples

``` r
if (FALSE) { # \dontrun{
require(dplyr)

# Database with species abundance and x and y coordinates
data("sppabund")

# Select data for a single species
some_sp <- sppabund %>%
  dplyr::filter(species == "Species one") %>%
  dplyr::select(-.part2, -.part3)

# Explore response variables
some_sp$ind_ha %>% range()
some_sp$ind_ha %>% hist()

# Here we balance number of absences
some_sp <-
  balance_dataset(some_sp, response = "ind_ha", absence_ratio = 0.2)

# Create a grid
gbm_grid <- expand.grid(
  interaction.depth = c(2, 4, 8, 16),
  n.trees = c(100, 500, 1000),
  n.minobsinnode = c(2, 5, 8),
  shrinkage = c(0.1, 0.5, 0.7),
  stringsAsFactors = FALSE
)

tuned_gbm <- tune_abund_gbm(
  data = some_sp,
  response = "ind_ha",
  predictors = c("bio12", "elevation", "sand"),
  predictors_f = c("eco"),
  partition = ".part",
  predict_part = TRUE,
  metrics = c("corr_pear", "mae"),
  grid = gbm_grid,
  distribution = "gaussian",
  n_cores = 3
)

tuned_gbm
} # }
```
