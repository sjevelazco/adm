# Fit and validate Quantile Regression Random Forests models

Fit and validate Quantile Regression Random Forests models

## Usage

``` r
fit_abund_qrf(
  data,
  response,
  predictors,
  predictors_f = NULL,
  fit_formula = NULL,
  partition,
  predict_part = FALSE,
  framework = "quantregForest",
  train_quantiles = c(0.5),
  eval_quantile = 0.5,
  mtry = length(c(predictors, predictors_f))/3,
  ntree = 2000,
  nodesize = 5,
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

- framework:

  character. Specifies the quantile regression framework to use. Either
  "quantregForest" (default) for traditional quantile regression forest
  from quantregForest package, or "grf" for generalized random forests
  from grf package.

- train_quantiles:

  numeric vector. Quantiles to be estimated during model training.
  Default c(0.5) for median.

- eval_quantile:

  numeric. Specific quantile to use for model evaluation metrics. Must
  be one of the values in train_quantiles. Default 0.5.

- mtry:

  numeric. Number of variables randomly sampled as candidates at each
  split. Default (length(c(predictors, predictors_f))/3)

- ntree:

  numeric. Number of trees to grow. This should not be set to too small
  a number, to ensure that every input row gets predicted at least a few
  times. Default 500

- nodesize:

  numeric. Minimum size of terminal nodes. Controls tree depth - larger
  values produce smaller trees. Default 5 for quantregForest framework.

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
} # }
```
