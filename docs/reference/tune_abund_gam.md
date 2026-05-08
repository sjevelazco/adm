# Fit and validate Generalized Additive Models with exploration of hyper-parameters that optimize performance

Fit and validate Generalized Additive Models with exploration of
hyper-parameters that optimize performance

## Usage

``` r
tune_abund_gam(
  data,
  response,
  predictors,
  predictors_f = NULL,
  fit_formula = NULL,
  sigma_formula = ~1,
  nu_formula = ~1,
  tau_formula = ~1,
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

- sigma_formula:

  formula. formula for fitting a model to the nu parameter. Usage
  sigma_formula = ~ precipt + temp

- nu_formula:

  formula. formula for fitting a model to the nu parameter. Usage
  nu_formula = ~ precipt + temp

- tau_formula:

  formula. formula for fitting a model to the tau parameter. Usage
  tau_formula = ~ precipt + temp

- partition:

  character. Column name with training and validation partition groups.

- predict_part:

  logical. Save predicted abundance for testing data. Default = FALSE

- grid:

  tibble or data.frame. A dataframe with 'distribution', 'inter' as
  columns and its values combinations as rows. If no grid is provided,
  function will create a default grid combining the next
  hyperparameters: distribution = all families selected by
  family_selector, inter = "automatic". In case one or more
  hyperparameters are provided, the function will complete the grid with
  the default values.

- metrics:

  character. Vector with one or more metrics from
  c("corr_spear","corr_pear","mae","pdisp","inter","slope")

- n_cores:

  numeric. Number of cores used in parallel processing.

- verbose:

  logical. If FALSE, disables all console messages. Default TRUE

## Value

A list object with:

- model: A "gamlss" object from gamlss package. This object can be used
  to predicting.

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
require(gamlss)

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
# Explore different family distributions
suitable_distributions <- family_selector(data = some_sp, response = "ind_ha")
suitable_distributions
# Create a grid
gam_grid <- expand.grid(
  distribution = suitable_distributions$family_call,
  stringsAsFactors = FALSE
)
# Tune a GAM model
tuned_gam <- tune_abund_gam(
  data = some_sp,
  response = "ind_ha",
  predictors = c("bio12", "elevation", "sand"),
  fit_formula = formula("ind_ha ~ bio12 + elevation + sand + eco"),
  sigma_formula = formula("ind_ha ~ bio12 + elevation"),
  nu_formula = formula("ind_ha ~ bio12 + elevation"),
  predictors_f = c("eco"),
  partition = ".part",
  predict_part = TRUE,
  metrics = c("corr_pear", "mae"),
  grid = gam_grid,
  n_cores = 3
)

tuned_gam
} # }
```
