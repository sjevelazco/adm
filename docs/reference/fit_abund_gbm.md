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

# Fit a GBM model
mgbm <- fit_abund_gbm(
  data = some_sp,
  response = "ind_ha",
  predictors = c("bio12", "elevation", "sand"),
  predictors_f = c("eco"),
  partition = ".part",
  distribution = "gaussian",
  n.trees = 100,
  interaction.depth = 5,
  n.minobsinnode = 5,
  shrinkage = 0.1,
  predict_part = TRUE
)

mgbm
} # }
```
